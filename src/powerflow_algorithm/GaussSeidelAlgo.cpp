// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GaussSeidelAlgo.h"

bool GaussSeidelAlgo::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                   CplxVect & V,
                                   const CplxVect & Sbus,
                                   const Eigen::VectorXi & slack_ids,
                                   const RealVect & slack_weights,  // currently unused
                                   const Eigen::VectorXi & pv,
                                   const Eigen::VectorXi & pq,
                                   int max_iter,
                                   real_type tol
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
    timer_total_nr_ += timer.duration();
    return res;
}

void GaussSeidelAlgo::one_iter(CplxVect & tmp_Sbus,
                                 const Eigen::SparseMatrix<cplx_type> & Ybus,
                                 const Eigen::VectorXi & pv,
                                 const Eigen::VectorXi & pq)
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
        tmp -= static_cast<cplx_type>(Ybus.row(k) * V_);
        tmp /= Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = static_cast<cplx_type>(Ybus.row(k) * V_);  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = tmp_Sbus.coeff(k) / V_.coeff(k);
        tmp = std::conj(tmp);
        tmp -= static_cast<cplx_type>(Ybus.row(k) * V_);
        tmp /= Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}
