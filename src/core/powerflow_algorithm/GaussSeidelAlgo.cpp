// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GaussSeidelAlgo.hpp"

namespace ls2g {

bool GaussSeidelAlgo::compute_pf(
    const EigenRefConstCplxSpMat     & Ybus,
    const Eigen::Ref<const CplxVect> & V,
    const Eigen::Ref<const CplxVect> & Sbus,
    const Eigen::Ref<const IntVect>  & slack_ids,
    const Eigen::Ref<const RealVect> & /*slack_weights*/,  // currently unused
    const Eigen::Ref<const IntVect>  & pv,
    const Eigen::Ref<const IntVect>  & pq,
    int                              max_iter,
    real_type                        tol
)
{
    /**
    pv: id of the pv buses
    pq: id of the pq buses

    **/
    // TODO check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
    // TODO Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
    reset_timer();
    err_ = ErrorType::NoError;
    auto timer = CustTimer();

    // TODO SLACK (for now i put all slacks as PV, except the first one)
    Eigen::VectorXi my_pv = retrieve_pv_with_slack(slack_ids, pv);

    // Row-major view + diagonal of Ybus, rebuilt on EVERY call: an O(nnz) transpose
    // plus a scan, even when Ybus has not moved since the last solve.
    //
    // Gating it on the change flags -- a reset, a dimension change, a new sparsity
    // pattern, new coefficients, the way LSGrid gates the rebuild of Ybus itself --
    // is NOT safe, and was tried. The batch sweeps mutate Ybus in place between
    // solves (YbusPolicy::Contingency::remove_from_Ybus does
    // `Ybus.coeffRef(i, j) -= value`) and then hand the algorithm a control saying
    // `tell_none_changed` from the second step onwards (BaseBatchSweep,
    // `needs_solver_init`). Against that pattern a gated cache answers from a stale
    // copy: reproduced directly -- solve, edit one coefficient, solve again with a
    // "nothing changed" control -- for 0.38 pu of voltage error on case118, silently.
    // nnz does not catch it either: an in-place `-=` leaves the sparsity pattern
    // alone.
    //
    // The DC family avoids all of this by holding its Ybus internally and being
    // handed each edit through update_internal_Ybus, which is DC-only
    // (AlgorithmSelector::update_internal_Ybus throws otherwise). Gating this would
    // need the same hook on the AC side. Whether a Gauss-Seidel can usefully reach
    // those sweeps at all is a separate question -- it is not a solver anyone runs a
    // contingency analysis with -- but "the caller cannot reach it" is not something
    // this function can check, so it rebuilds.
    refresh_ybus_cache(Ybus);

    V_ = V;
    Vm_ = V_.array().abs();  // update Vm and Va again in case
    Va_ = V_.array().arg();  // we wrapped around with a negative Vm

    // first check, if the problem is already solved, i stop there
    RealVect F = _evaluate_Fx(Ybus, V, Sbus, my_pv, pq);
    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    CplxVect tmp_Sbus = Sbus;
    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;

        // ###########################
        // the Gauss Seidel Algorithm
        // ###########################
        // https://www.sanfoundry.com/cpp-program-implement-gauss-seidel-method/
        // https://fr.mathworks.com/matlabcentral/fileexchange/14089-gauss-seidel-load-flow-analysis
        // https://github.com/rwl/PYPOWER/blob/master/pypower/gausspf.py

        auto timer2 = CustTimer();
        one_iter(tmp_Sbus, Ybus, my_pv, pq);
        timer_solve_ += timer2.duration();

        // #####################
        // stopping criteria
        // #####################
        F = _evaluate_Fx(Ybus, V_, tmp_Sbus, my_pv, pq);
        bool tmp = F.allFinite();
        if(!tmp){
            err_ = ErrorType::InifiniteValue;
            break; // divergence due to Nans
        }
        converged = _check_for_convergence(F, tol);
    }
    if(!converged){
        if (err_ == ErrorType::NoError) err_ = ErrorType::TooManyIterations;
        res = false;
    }
    Vm_ = V_.array().abs();  // update Vm and Va again in case
    Va_ = V_.array().arg();  // we wrapped around with a negative Vm

    // The per-bus mismatch LSGrid::compute_results reads back to work out what the
    // generators actually produced (BaseAlgo::mis_bus_, see fills_bus_mismatch()).
    // Written from the ORIGINAL Sbus, not from `tmp_Sbus`: the caller's convention
    // is `V .* conj(Ybus . V) - Sbus` against the injections IT handed us, and
    // tmp_Sbus is this algorithm's own working copy. This family carries no
    // extension state (no distributed-slack unknown, no voltage-control group), so
    // there is nothing to fold in and the plain formula IS the answer -- which is
    // also why the reference bus keeps its full imbalance here, exactly where a
    // single-slack solve wants it.
    //
    // Filled even when the solve did not converge: process_results discards the
    // results of a diverged solve anyway, and leaving a stale buffer from a
    // previous solve behind would be worse than an honest one.
    ybus_v_.noalias() = Ybus * V_;
    mis_bus_ = V_.array() * ybus_v_.array().conjugate() - Sbus.array();

    timer_total_nr_ += timer.duration();
    return res;
}

void GaussSeidelAlgo::refresh_ybus_cache(const EigenRefConstCplxSpMat & Ybus)
{
    ybus_rowmajor_ = Ybus;  // column-major -> row-major (one O(nnz) transpose)
    const Eigen::Index n = Ybus.rows();
    if(ybus_diag_.size() != n) ybus_diag_.resize(n);
    ybus_diag_.setZero();
    // Ybus is column-major and compressed, so a column's row indices ascend: the
    // moment one passes `col` there is no diagonal coefficient left to find in it.
    // Scanning on to the end of every column read the whole lower triangle for
    // nothing.
    for(int col = 0; col < Ybus.outerSize(); ++col){
        for(EigenRefConstCplxSpMat::InnerIterator it(Ybus, col); it; ++it){
            if(it.row() > col) break;   // past the diagonal
            if(it.row() == col){ ybus_diag_.coeffRef(col) = it.value(); break; }
        }
    }
}

void GaussSeidelAlgo::one_iter(
    Eigen::Ref<CplxVect> tmp_Sbus,
    const EigenRefConstCplxSpMat & /*Ybus*/,  // read through the caches below
    const Eigen::Ref<const IntVect> & pv,
    const Eigen::Ref<const IntVect> & pq)
{
    // do an update with the standard GS algorithm
    cplx_type tmp;

    const int n_pv = static_cast<int>(pv.size());
    const int n_pq = static_cast<int>(pq.size());

    // update PQ buses
    for(int k_tmp=0; k_tmp < n_pq; ++k_tmp)
    {
        int k = pq.coeff(k_tmp);
        tmp = tmp_Sbus.coeff(k) / V_.coeff(k);
        tmp = std::conj(tmp);
        tmp -= ybus_row_times_V(k);
        tmp /= ybus_diag_.coeff(k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = ybus_row_times_V(k);  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = tmp_Sbus.coeff(k) / V_.coeff(k);
        tmp = std::conj(tmp);
        tmp -= ybus_row_times_V(k);
        tmp /= ybus_diag_.coeff(k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}

} // namespace ls2g
