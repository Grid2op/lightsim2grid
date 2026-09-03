// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// #include "DCSolver.h"
#include <cassert>
#include <limits>  // for nans
#include <cmath>  // for nans

template<class LinearSolver>
bool BaseDCAlgo<LinearSolver>::compute_pf_dc(
    const EigenRefConstRealSpMat & Bbus,
    const Eigen::Ref<const CplxVect> & V,
    const Eigen::Ref<const RealVect> & Pbus,
    const Eigen::Ref<const IntVect> & slack_ids,
    const Eigen::Ref<const RealVect> & slack_weights,
    const Eigen::Ref<const IntVect> & pv,
    const Eigen::Ref<const IntVect> & pq
)
{
    // V is used the following way: at pq buses it's completely ignored. For pv bus only the magnitude is used,
    //   and for the slack bus both the magnitude and the angle are used.
    if(!is_linear_solver_valid()) {
        return false;
    }
    reset_timer();
    // std::cout << "\n\n____________________________________________________\n";
    // std::cout << "need_factorize " <<  need_factorize_ <<"\n"; 
    // std::cout << "need_refactor " <<  need_refactor_ <<"\n"; 

    auto timer = CustTimer();
    if(need_factorize_ ||
       _solver_control.need_reset_solver() ||
       _solver_control.has_dimension_changed() ||
       _solver_control.has_slack_participate_changed() ||  // the full "Bbus without slack" has changed, everything needs to be recomputed
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero()
       ){
       reset();
    //    std::cout << "need_factorize_ 1\n";
       // at this stage need_factorize_ is set also to true
    }

    sizeYbus_with_slack_ = static_cast<int>(Bbus.rows());

    // "handle disconnected grid" mode: when masked_buses_ is non-empty the rhs and the
    // factorized matrix are built in a masked-aware way (see below). We never touch the
    // persistent dcYbus_noslack_ / dcSbus_noslack_ members so the next (possibly
    // un-masked) contingency keeps reusing them as today.
    const bool has_mask = !masked_buses_.empty();

    // hvdc angle-droop: refresh the droop data; a change (eg a `status_droop`
    // flip between two solves) modifies the dc matrix values, so it forces a
    // rebuild + refactorization
    const bool droop_changed = update_hvdc_droop_data();

    #ifdef __COUT_TIMES
        auto timer_preproc = CustTimer();
    #endif // __COUT_TIMES

    if(need_factorize_ ||
       _solver_control.has_pv_changed() || _solver_control.has_pq_changed()) {

        // The first slack bus is the angle reference (removed from the linear system); the other
        // slack buses are kept as PV. The active power imbalance is distributed across all the
        // participating slack buses (by slack_weights) further down (see "distributed slack").
        auto timer_pre = CustTimer();
        my_pv_ = retrieve_pv_with_slack(slack_ids, pv);

        // find the slack buses
        slack_buses_ids_solver_ = extract_slack_bus_id(my_pv_, pq, sizeYbus_with_slack_);
        sizeYbus_without_slack_ = sizeYbus_with_slack_ - slack_buses_ids_solver_.size();

        // corresp bus -> solverbus
        fill_mat_bus_id(sizeYbus_with_slack_);
        timer_pre_proc_ += timer_pre.duration();
    }

    // remove the slack bus from Ybus
    if(need_factorize_ ||
       droop_changed ||
       _solver_control.need_recompute_ybus() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero()) {
        auto timer_pre = CustTimer();
        fill_dcYbus_noslack(sizeYbus_with_slack_, Bbus);
        timer_pre_proc_ += timer_pre.duration();
        need_factorize_ = true;  // force a call to "factor" the linear solver as the lhs (Bbus) changed
        // std::cout << "need_factorize_ 2\n"; 
        // no need to refactor if ybus did not change
    }
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc: preproc: " << 1000. * timer_preproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES
    
    // initialize the solver (only if needed)
    #ifdef __COUT_TIMES
        auto timer_solve = CustTimer();
    #endif // __COUT_TIMES

    // remove the slack bus from Pbus
    // (skipped in mask mode: a masked-aware rhs is built into a local vector further
    //  down, leaving the persistent dcSbus_noslack_ untouched for un-masked solves)
    if(!has_mask && (need_factorize_ || _solver_control.need_recompute_sbus() || _solver_control.has_slack_weight_changed())){
        // std::cout << "\t\t\tneed to dcSbus_noslack_\n";
        auto timer_pre = CustTimer();
        dcSbus_noslack_ = RealVect::Constant(sizeYbus_without_slack_, my_zero_);
        for (int k=0; k < sizeYbus_with_slack_; ++k){
            if(mat_bus_id_(k) == -1) continue;  // I don't add anything to the slack bus
            const int col_res = mat_bus_id_(k);
            dcSbus_noslack_(col_res) = Pbus(k);
        }

        // distributed slack: spread the global active power imbalance T = -sum(Pbus) over the
        // participating slack buses, proportionally to slack_weights (normalized, sum = 1).
        // The non-reference slack buses (still in the system) get their share added to the reduced
        // injection so the angles reflect the distribution; the reference bus (removed from the
        // system) absorbs its share implicitly. With a single slack (weight 1 at the reference)
        // this is a no-op and reduces to the historical "reference absorbs everything" behaviour.
        if(slack_weights.size() == sizeYbus_with_slack_){
            const real_type imbalance = -Pbus.sum();
            for (int k=0; k < sizeYbus_with_slack_; ++k){
                if(slack_weights(k) <= my_zero_) continue;  // not a participating slack bus
                const int col_res = mat_bus_id_(k);
                if(col_res == -1) continue;  // reference bus: removed from the system, share is implicit
                dcSbus_noslack_(col_res) += slack_weights(k) * imbalance;
            }
        }
        add_droop_to_dcSbus();
        timer_pre_proc_ += timer_pre.duration();
    }

    // "handle disconnected grid" mode: build the masked working matrix (identity rows
    // for the masked buses) on a copy, so the persistent dcYbus_noslack_ is untouched.
    // Identity rows only overwrite *existing* structural entries (the diagonal and the
    // within-island off-diagonals), so the sparsity pattern -- and thus the symbolic
    // factorization -- is unchanged (a numeric refactorize is enough, never analyze).
    std::vector<char> is_masked;  // size sizeYbus_with_slack_ when has_mask
    Eigen::SparseMatrix<real_type> masked_mat;
    if(has_mask){
        is_masked.assign(sizeYbus_with_slack_, 0);
        for(int b : masked_buses_) if(b >= 0 && b < sizeYbus_with_slack_) is_masked[b] = 1;
        // which *reduced* rows are masked (the reference slack maps to -1 and is skipped)
        std::vector<char> masked_row(sizeYbus_without_slack_, 0);
        for(int b : masked_buses_){
            if(b < 0 || b >= sizeYbus_with_slack_) continue;
            const int mr = mat_bus_id_(b);
            if(mr != -1) masked_row[mr] = 1;
        }
        masked_mat = dcYbus_noslack_;  // copy, leaving the persistent matrix intact
        // single in-place pass over the (column-major) copy: every masked row becomes
        // the identity row e_mr (off-diagonals -> 0, diagonal -> 1). Only *existing*
        // structural entries are written (via valueRef), so the sparsity pattern -- and
        // therefore the symbolic factorization -- is left unchanged (refactorize only).
        for(int col = 0; col < masked_mat.outerSize(); ++col){
            for(typename Eigen::SparseMatrix<real_type>::InnerIterator it(masked_mat, col); it; ++it){
                const int row = static_cast<int>(it.row());
                if(!masked_row[row]) continue;
                it.valueRef() = (row == col) ? 1. : 0.;
            }
        }
    }
    // the system matrix actually handed to the linear solver
    Eigen::SparseMatrix<real_type> & sys_mat = has_mask ? masked_mat : dcYbus_noslack_;

    // analyze (structure) + factorize (values) if topology changed
    bool factorized_now = false;
    if(need_factorize_){
        // std::cout << "\t\t\tneed to factorize\n";
        ErrorType status_init = _linear_solver.analyze(sys_mat);
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            timer_total_nr_ += timer.duration();
            return false;
        }

        status_init = _linear_solver.factorize(sys_mat);
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            timer_total_nr_ += timer.duration();
            return false;
        }
        need_factorize_ = false;  // done just above
        need_refactor_ = false;  // no need to refactor as a factor as been called just now
        factorized_now = true;
        // std::cout << "need_factorize_ 3\n";
    }

    // solve for theta: Sbus = dcY . theta (make a copy to keep dcSbus_noslack_)
    auto timer_pre = CustTimer();
    RealVect Va_dc_without_slack;
    if(has_mask){
        // masked-aware rhs: the masked island is not simulated, so its buses inject
        // nothing and are excluded from the slack imbalance; their (identity) rows get a
        // 0 rhs => theta = 0. Built locally so dcSbus_noslack_ stays valid for un-masked
        // contingencies.
        Va_dc_without_slack = RealVect::Constant(sizeYbus_without_slack_, my_zero_);
        for(int k = 0; k < sizeYbus_with_slack_; ++k){
            if(is_masked[k]) continue;
            const int col_res = mat_bus_id_(k);
            if(col_res == -1) continue;  // reference slack: removed from the system
            Va_dc_without_slack(col_res) = Pbus(k);
        }
        if(slack_weights.size() == sizeYbus_with_slack_){
            real_type imbalance = my_zero_;
            for(int k = 0; k < sizeYbus_with_slack_; ++k){
                if(!is_masked[k]) imbalance -= Pbus(k);  // -sum(Pbus) over live buses only
            }
            for(int k = 0; k < sizeYbus_with_slack_; ++k){
                if(is_masked[k] || slack_weights(k) <= my_zero_) continue;
                const int col_res = mat_bus_id_(k);
                if(col_res == -1) continue;  // reference slack: share is implicit
                Va_dc_without_slack(col_res) += slack_weights(k) * imbalance;
            }
        }
        // angle-droop hvdc: only lines fully inside the live component are stamped (a
        // droop line crossing into a masked island is out of scope for v1, see header).
        const int nb_droop_m = hvdc_droop_data_.size();
        for(int kk = 0; kk < nb_droop_m; ++kk){
            if(hvdc_droop_data_.status(kk) != 0) continue;  // saturated: fixed injection (in Pbus)
            const int b1 = hvdc_droop_data_.bus1(kk);
            const int b2 = hvdc_droop_data_.bus2(kk);
            if(is_masked[b1] || is_masked[b2]) continue;
            const int m1 = mat_bus_id_(b1);
            const int m2 = mat_bus_id_(b2);
            if(m1 != -1) Va_dc_without_slack(m1) -= hvdc_droop_data_.p0(kk);
            if(m2 != -1) Va_dc_without_slack(m2) += hvdc_droop_data_.p0(kk);
        }
    } else {
        Va_dc_without_slack = dcSbus_noslack_;
    }
    timer_pre_proc_ += timer_pre.duration();

    // refactorize (numeric only) when Ybus changed (n-1) or in mask mode (the working
    // matrix differs from the one currently factorized). Skipped right after a factorize.
    if(!factorized_now && (need_refactor_ || has_mask)){
        // we should end-up here only in case of n-1 simulation (handled in contingency analysis)
        // set to true in update_internal_Ybus
        // std::cout << "\t\t\tneed to refactorize\n";
        ErrorType error = _linear_solver.refactorize(sys_mat);
        if(error != ErrorType::NoError){
            err_ = error;
            timer_total_nr_ += timer.duration();
            return false;
        }
    }
    {
        ErrorType error = _linear_solver.solve(Va_dc_without_slack);
        if(error != ErrorType::NoError){
            err_ = error;
            timer_total_nr_ += timer.duration();
            return false;
        }
    }
    
    auto timer_mismatch = CustTimer();
    // single-pass validity check: maxCoeff<PropagateNaN> propagates NaN (the default maxCoeff
    // skips it, because `NaN > x` is false), and +/-Inf yields a non-finite max, so both fail
    // the isfinite test below. This is what catches a singular system (e.g. an islanded grid).
    const real_type maxAbsVal = Va_dc_without_slack.cwiseAbs().maxCoeff<Eigen::PropagateNaN>();
    if(!std::isfinite(maxAbsVal) || maxAbsVal >= _max_dc_voltage_angle){
        // for convergence, all values should be finite
        // and it's not realistic if some Va are too high
        err_ = ErrorType::SolverSolve;
        V_ = CplxVect();
        Vm_ = RealVect();
        Va_ = RealVect();
        timer_total_nr_ += timer.duration();
        return false;
    }

    #ifdef __COUT_TIMES
        std::cout << "\t dc solve: " << 1000. * timer_solve.duration() << "ms" << std::endl;
        auto timer_postproc = CustTimer();
    #endif // __COUT_TIMES

    // retrieve back the results in the proper shape (add back the slack bus)
    // write directly into Va_: slack positions stay 0, non-slack scattered via precomputed indices
    if(Va_.size() != sizeYbus_with_slack_) Va_.resize(sizeYbus_with_slack_);
    Va_.setZero();
    Va_(nonslack_ybus_ids_) = Va_dc_without_slack;
    Va_.array() += std::arg(V(slack_buses_ids_solver_(0)));

    // batch DC theta-only fast path (see BaseAlgo::set_lazy_v): the caller only
    // wants Va_ (theta) from this call, and will reconstruct Vm_/V_ itself, later,
    // only if/when it actually needs them. Vm_/V_ are left untouched (possibly
    // stale) below -- the caller must not read get_Vm()/get_V() while lazy.
    if(!_lazy_v_){
        // DC solves for angles only: the magnitudes it reports are the ones it was
        // GIVEN, in `V` -- the caller's Vinit, with the regulated buses snapped to
        // their setpoint by pre_process_solver. So they must be taken from `V`
        // every time.
        // This used to be guarded by `_solver_control.has_v_changed() || ...`,
        // which asks whether the GRID's voltage setpoints changed -- a different
        // question, and one that says nothing about the vector just passed in. Two
        // dc_pf() calls with different Vinit magnitudes therefore returned the
        // first call's magnitudes for every bus not pinned by a controller,
        // whenever the solver control said "nothing changed" (ie after every
        // `unset_changes()` -- LightSimBackend calls it after every step -- and now
        // after every DC powerflow, since cache reuse is automatic).
        // What the guard saved was one `abs()` over nb_bus doubles, next to a
        // sparse triangular solve of the same dimension: nothing measurable.
        Vm_ = V.array().abs();

        // compute complex voltages with std::polar: uses hardware sincos, no temporaries, fills V and V_ in one pass
        V_.resize(sizeYbus_with_slack_);
        for(int i = 0; i < sizeYbus_with_slack_; ++i){
            V_[i] = std::polar(Vm_[i], Va_[i]);
        }
    }
    nr_iter_ = 1;
    need_refactor_ = false;  // no need to redo it in general cases
    timer_mismatch_ = timer_mismatch.duration();
    // std::cout << "need_refactor " <<  need_refactor_ <<"\n"; 
    // std::cout << "end powerflow\n";
    // std::cout << "---------------------------------------\n";
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc postproc: " << 1000. * timer_postproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES
    timer_total_nr_ += timer.duration();
    return true;
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::fill_mat_bus_id(int nb_bus_solver){
    mat_bus_id_ = Eigen::VectorXi::Constant(nb_bus_solver, -1);
    nonslack_ybus_ids_.resize(sizeYbus_without_slack_);
    int solver_id = 0;
    for (int ybus_id=0; ybus_id < nb_bus_solver; ++ybus_id){
        if(isin(ybus_id, slack_buses_ids_solver_)) continue;  // I don't add anything to the slack bus
        mat_bus_id_(ybus_id) = solver_id;
        nonslack_ybus_ids_(solver_id) = ybus_id;
        ++solver_id;
    }
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::fill_dcYbus_noslack(int nb_bus_solver, const Eigen::Ref<const Eigen::SparseMatrix<real_type>> & ref_mat){
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    remove_slack_buses(nb_bus_solver, ref_mat, dcYbus_noslack_);
    add_droop_to_dcYbus();
}

template<class LinearSolver>
bool BaseDCAlgo<LinearSolver>::update_hvdc_droop_data(){
    HvdcDroopSolverData new_data;
    fill_hvdc_droop_data_from_grid(BaseAlgo::lsgrid_ptr_, new_data, false);
    bool changed = (new_data.size() != hvdc_droop_data_.size());
    if(!changed && (new_data.size() > 0)){
        changed = (new_data.bus1.array() != hvdc_droop_data_.bus1.array()).any() ||
                  (new_data.bus2.array() != hvdc_droop_data_.bus2.array()).any() ||
                  (new_data.status.array() != hvdc_droop_data_.status.array()).any() ||
                  (new_data.k.array() != hvdc_droop_data_.k.array()).any() ||
                  (new_data.p0.array() != hvdc_droop_data_.p0.array()).any();
    }
    hvdc_droop_data_ = new_data;
    return changed;
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::add_droop_to_dcYbus(){
    // the k * (theta1 - theta2) term of a linear-mode droop line behaves like
    // a branch of susceptance k between the two buses
    const int nb_droop = hvdc_droop_data_.size();
    if(nb_droop == 0) return;
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(4 * nb_droop);
    for(int k = 0; k < nb_droop; ++k){
        if(hvdc_droop_data_.status(k) != 0) continue;  // saturated: fixed injection, in Sbus
        const int m1 = mat_bus_id_(hvdc_droop_data_.bus1(k));
        const int m2 = mat_bus_id_(hvdc_droop_data_.bus2(k));
        const real_type k_droop = hvdc_droop_data_.k(k);
        if(m1 != -1) tripletList.push_back({m1, m1, k_droop});
        if(m2 != -1) tripletList.push_back({m2, m2, k_droop});
        if((m1 != -1) && (m2 != -1)){
            tripletList.push_back({m1, m2, -k_droop});
            tripletList.push_back({m2, m1, -k_droop});
        }
    }
    if(tripletList.empty()) return;
    Eigen::SparseMatrix<real_type> droop_mat(sizeYbus_without_slack_, sizeYbus_without_slack_);
    droop_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    dcYbus_noslack_ += droop_mat;
    dcYbus_noslack_.makeCompressed();
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::add_droop_to_dcSbus(){
    // the p0 constant of a linear-mode droop line leaves bus 1 and enters bus 2
    const int nb_droop = hvdc_droop_data_.size();
    for(int k = 0; k < nb_droop; ++k){
        if(hvdc_droop_data_.status(k) != 0) continue;  // saturated: fixed injection, in Sbus
        const int m1 = mat_bus_id_(hvdc_droop_data_.bus1(k));
        const int m2 = mat_bus_id_(hvdc_droop_data_.bus2(k));
        if(m1 != -1) dcSbus_noslack_(m1) -= hvdc_droop_data_.p0(k);
        if(m2 != -1) dcSbus_noslack_(m2) += hvdc_droop_data_.p0(k);
    }
}

template<class LinearSolver>
template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
void BaseDCAlgo<LinearSolver>::remove_slack_buses(int nb_bus_solver, const Eigen::Ref<const Eigen::SparseMatrix<ref_mat_type>> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat){
    // Deleting the slack rows and columns of an already compressed matrix needs
    // no triplet list. mat_bus_id_ is a monotone compaction -- fill_mat_bus_id
    // hands out consecutive ids in bus order, skipping the slack buses -- and
    // the inner iterator walks each column in increasing row order, so the
    // surviving coefficients come out in exactly the order a compressed
    // column-major matrix stores them: sorted, and one per coefficient. There
    // is nothing to sort and nothing to merge, so res_mat is written in place,
    // in two passes over ref_mat -- one to size each column, one to fill it.
    //
    // This used to collect the survivors into a vector of Eigen::Triplet and
    // hand that to setFromTriplets, which re-derived the order it was already
    // in the hard way (bucket by row, then transpose into columns, both through
    // a temporary SparseMatrix) for 5% of a DC solve on a 9241-bus grid.
    const int dim = sizeYbus_without_slack_;  // TODO dist slack: -1 or -mat_bus_id_.size() here ????
    typedef typename Eigen::Ref<const Eigen::SparseMatrix<ref_mat_type>>::InnerIterator RefInnerIt;

    // resize() drops res_mat into compressed mode with a zeroed outer array and
    // keeps what the previous build allocated. The counts are accumulated
    // straight into that outer array, which a prefix sum then turns into the
    // column starts -- the shape it has to end up in anyway.
    res_mat.resize(dim, dim);
    int * outer = res_mat.outerIndexPtr();
    for (int k=0; k < nb_bus_solver; ++k){
        const int col_res = mat_bus_id_(k);
        if(col_res == -1) continue;  // slack column, dropped whole
        assert(col_res < dim && "remove_slack_buses: mat_bus_id_ labels more columns than the matrix has");
        int nb_kept = 0;
        for (RefInnerIt it(ref_mat, k); it; ++it)
            if(mat_bus_id_(static_cast<int>(it.row())) != -1) ++nb_kept;
        outer[col_res + 1] = nb_kept;
    }
    for (int c = 0; c < dim; ++c) outer[c + 1] += outer[c];

    res_mat.resizeNonZeros(outer[dim]);
    outer = res_mat.outerIndexPtr();
    int * inner = res_mat.innerIndexPtr();
    real_type * values = res_mat.valuePtr();
    for (int k=0; k < nb_bus_solver; ++k){
        const int col_res = mat_bus_id_(k);
        if(col_res == -1) continue;
        int pos = outer[col_res];
        for (RefInnerIt it(ref_mat, k); it; ++it){
            const int row_res = mat_bus_id_(static_cast<int>(it.row()));
            if(row_res == -1) continue;  // slack row, dropped
            assert(row_res < dim && "remove_slack_buses: mat_bus_id_ labels more rows than the matrix has");
            inner[pos] = row_res;
            values[pos] = std::real(it.value());
            ++pos;
        }
        assert(pos == outer[col_res + 1] &&
               "remove_slack_buses: the two passes disagree on a column's size");
    }

#ifndef NDEBUG
    // What the linear solvers are entitled to assume, now that the arrays are
    // written here rather than by Eigen. The strictly-increasing check is also
    // what would catch mat_bus_id_ losing its monotonicity, which is the
    // property this function is built on.
    assert(res_mat.isCompressed() && "remove_slack_buses must leave a compressed matrix");
    for (int c = 0; c < dim; ++c)
        for (int p = outer[c] + 1; p < outer[c + 1]; ++p)
            assert(inner[p - 1] < inner[p] &&
                   "remove_slack_buses: row indices must be strictly increasing");
#endif
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::reset(){
    BaseAlgo::reset();
    _linear_solver.reset();
    need_factorize_ = true;
    need_refactor_ = true;
    sizeYbus_with_slack_ = 0;
    sizeYbus_without_slack_ = 0;
    dcSbus_noslack_ = RealVect();
    dcYbus_noslack_ = Eigen::SparseMatrix<real_type>();
    my_pv_ = Eigen::VectorXi();
    slack_buses_ids_solver_ = Eigen::VectorXi();
    mat_bus_id_ = Eigen::VectorXi();
    nonslack_ybus_ids_ = Eigen::VectorXi();
    hvdc_droop_data_.clear();
    masked_buses_.clear();
}

template<class LinearSolver>
RealMat BaseDCAlgo<LinearSolver>::get_ptdf(){
    auto timer = CustTimer();
    Eigen::SparseMatrix<real_type> Bf_T_with_slack;
    RealMat PTDF;
    RealVect rhs = RealVect::Zero(sizeYbus_without_slack_);  // TODO dist slack: -1 or -mat_bus_id_.size() here ????
    // TODO PTDF: sparse matrix ?
    // TODO PTDF: distributed slack
    // TODO PTDF: check that the solver has converged


    //extract the Bf matrix
    BaseAlgo::get_Bf_transpose(Bf_T_with_slack);  // Bf_T_with_slack : [bus_id, line_or_trafo_id]
    const int nb_bus = Bf_T_with_slack.rows();
    const int nb_pow_tr = Bf_T_with_slack.cols();  // cols and not rows because Bf_T_with_slack is transposed
    
    // get the index of buses without slacks
    std::vector<int> ind_no_slack_;
    ind_no_slack_.reserve(nb_bus);
    for(int bus_id = 0; bus_id < nb_bus; ++bus_id){
        if(mat_bus_id_(bus_id) == BaseConstants::_deactivated_bus_id) continue;
        ind_no_slack_.push_back(bus_id);
    }
    const Eigen::VectorXi ind_no_slack = Eigen::VectorXi::Map(&ind_no_slack_[0], ind_no_slack_.size());

    // solve iteratively the linear systems (one per powerline)
    PTDF = RealMat::Zero(Bf_T_with_slack.cols(), Bf_T_with_slack.rows());  // rows and cols are "inverted" because the matrix Bf is transposed
    for (int row_id=0; row_id < nb_pow_tr; ++row_id){
        // build the rhs vector
        for (typename Eigen::SparseMatrix<real_type>::InnerIterator it(Bf_T_with_slack, row_id); it; ++it)
        {
            const auto bus_id = it.row();
            if(mat_bus_id_(bus_id) == BaseConstants::_deactivated_bus_id) continue;  // I don't add anything if it's the slack
            const auto col_res = mat_bus_id_(bus_id);
            rhs[col_res] = it.value();
        }
        // solve the linear system
        _linear_solver.solve(rhs);

        // assign results to the PTDF matrix
        PTDF(row_id, ind_no_slack) = rhs;

        // reset the rhs vector to 0.
        rhs.array() = 0.;
        // rhs = RealVect::Zero(sizeYbus_without_slack_);
    }
    timer_ptdf_ = timer.duration();
    // TODO PTDF: if the solver can solve the MAT directly, do that instead
    return PTDF;
}

template<class LinearSolver>
RealMat BaseDCAlgo<LinearSolver>::get_lodf(const Eigen::Ref<const IntVect> & from_bus,
                                           const Eigen::Ref<const IntVect> & to_bus){
    auto timer = CustTimer();
    const RealMat PTDF = get_ptdf();  // size n_line x n_bus
    RealMat LODF = RealMat::Zero(from_bus.size(), from_bus.rows());  // nb_line, nb_line
    const real_type tol_equal_float = _tol_equal_float;
    for(Eigen::Index line_id=0; line_id < from_bus.size(); ++line_id){
        auto f_bus = from_bus(line_id);
        auto t_bus = to_bus(line_id);
        if ((f_bus == BaseConstants::_deactivated_bus_id) || (t_bus == BaseConstants::_deactivated_bus_id)){
            // element carries no DC flow (disconnected, or half-open -- see
            // TwoSidesContainer_rxh_A::fillBdc's "disco on one side == disco on
            // both sides" convention): "outaging" it has no effect anywhere in
            // the grid, i.e. an identity column. Its own row is already all-0
            // too, since its PTDF row is 0 (fillBf_for_PTDF excludes it from Bf).
            LODF(line_id, line_id) = 1.;
            continue;
        }
        LODF.col(line_id).array() = PTDF.col(f_bus).array() - PTDF.col(t_bus).array();
        const real_type diag_coeff = LODF(line_id, line_id);
        if (abs(diag_coeff - 1.) > tol_equal_float){
            LODF.col(line_id).array() /= (1. - diag_coeff);
            LODF(line_id, line_id) = -1.;
        }else{
            LODF.col(line_id).array() = std::numeric_limits<real_type>::quiet_NaN();
        }
    }
    timer_lodf_ = timer.duration();
    return LODF;
}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCAlgo<LinearSolver>::get_bsdf(){
    // TODO
    return dcYbus_noslack_;

}
