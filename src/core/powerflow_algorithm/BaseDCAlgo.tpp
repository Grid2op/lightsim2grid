// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// #include "DCSolver.h"
#include <limits>  // for nans
#include <cmath>  // for nans

template<class LinearSolver>
bool BaseDCAlgo<LinearSolver>::compute_pf_dc(const Eigen::SparseMatrix<real_type> & Bbus,
                                             CplxVect & V,
                                             const RealVect & Pbus,
                                             Eigen::Ref<const IntVect> slack_ids,
                                             const RealVect & slack_weights,
                                             Eigen::Ref<const IntVect> pv,
                                             Eigen::Ref<const IntVect> pq
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

    // remove the slack bus from Bbus
    if(need_factorize_ ||
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
    if(need_factorize_ || _solver_control.need_recompute_sbus() || _solver_control.has_slack_weight_changed()){
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
        timer_pre_proc_ += timer_pre.duration();
    }

    // analyze (structure) + factorize (values) if topology changed
    if(need_factorize_){
        // std::cout << "\t\t\tneed to factorize\n";
        auto timer_an = CustTimer();
        ErrorType status_init = _linear_solver.analyze(dcYbus_noslack_);
        const double dur_an = timer_an.duration();
        timer_initialize_ += dur_an;
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            timer_total_nr_ += timer.duration();
            return false;
        }

        auto timer_fac = CustTimer();
        status_init = _linear_solver.factorize(dcYbus_noslack_);
        const double dur_fact = timer_fac.duration();
        timer_factor_ += dur_fact;
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            timer_total_nr_ += timer.duration();
            return false;
        }
        need_factorize_ = false;  // done just above
        need_refactor_ = false;  // no need to refactor as a factor as been called just now
        // std::cout << "need_factorize_ 3\n"; 
    }

    // solve for theta: Sbus = dcY . theta (make a copy to keep dcSbus_noslack_)
    auto timer_pre = CustTimer();
    RealVect Va_dc_without_slack = dcSbus_noslack_;      
    timer_pre_proc_ += timer_pre.duration(); 
    
    // std::cout << "\t\tBaseDCAlgo.tpp: dcYbus_noslack_ (max): " << dcYbus_noslack_.coeffs().maxCoeff() << std::endl;  // TODO DEBUG WINDOWS
    // std::cout << "\t\tBaseDCAlgo.tpp: dcYbus_noslack_ (sum): " << dcYbus_noslack_.coeffs().abs().sum() << std::endl;  // TODO DEBUG WINDOWS
    // std::cout << "\t\tBaseDCAlgo.tpp: Va_dc_without_slack (inf norm): " << Va_dc_without_slack.lpNorm<Eigen::Infinity>() << std::endl;  // TODO DEBUG WINDOWS
    // std::cout << "\t\tBaseDCAlgo.tpp: Va_dc_without_slack (l1 norm): " << Va_dc_without_slack.lpNorm<1>() << std::endl;  // TODO DEBUG WINDOWS
    // std::cout << "\t\tBaseDCAlgo.tpp:  V (l1 norm): " <<  V.lpNorm<1>() << std::endl;  // TODO DEBUG WINDOWS
    // std::cout << "\t\tBaseDCAlgo.tpp:  Sbus (l1 norm): " <<  Sbus.lpNorm<1>() << std::endl;  // TODO DEBUG WINDOWS
    if(need_refactor_){
        // we should end-up here only in case of n-1 simulation (handled in contingency analysis)
        // set to true in update_internal_Ybus
        // std::cout << "\t\t\tneed to refactorize\n";
        auto timer_s = CustTimer();
        ErrorType error = _linear_solver.refactorize(dcYbus_noslack_);
        const double dur_refacto = timer_s.duration();
        timer_refactor_ += dur_refacto;
        if(error != ErrorType::NoError){
            err_ = error;
            timer_total_nr_ += timer.duration();
            return false;
        }
    }
    {
        auto timer_s = CustTimer();
        ErrorType error = _linear_solver.solve(Va_dc_without_slack);
        const double dur_solve = timer_s.duration();
        timer_solve_ += dur_solve;
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
        V = CplxVect();
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

    // add the Voltage setpoints of the generator (only recompute when magnitudes may have changed)
    // size check guards against Vm_ being cleared by a previous failed contingency
    if(need_factorize_ ||
        Vm_.size() != sizeYbus_with_slack_ ||  // in contingency analysis, Vm is cleared if a cont diverges, so I need to recompute it
       _solver_control.has_v_changed() ||
       _solver_control.has_dimension_changed() ||
       _solver_control.has_pv_changed() || _solver_control.has_pq_changed()){
        Vm_ = V.array().abs();
    }

    // compute complex voltages with std::polar: uses hardware sincos, no temporaries, fills V and V_ in one pass
    V_.resize(sizeYbus_with_slack_);
    V.resize(sizeYbus_with_slack_);
    for(int i = 0; i < sizeYbus_with_slack_; ++i){
        V_[i] = V[i] = std::polar(Vm_[i], Va_[i]);
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
void BaseDCAlgo<LinearSolver>::fill_dcYbus_noslack(int nb_bus_solver, const Eigen::SparseMatrix<real_type> & ref_mat){
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    remove_slack_buses(nb_bus_solver, ref_mat, dcYbus_noslack_);
}

template<class LinearSolver>
template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
void BaseDCAlgo<LinearSolver>::remove_slack_buses(int nb_bus_solver, const Eigen::SparseMatrix<ref_mat_type> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat){
    res_mat = Eigen::SparseMatrix<real_type>(sizeYbus_without_slack_, sizeYbus_without_slack_);  // TODO dist slack: -1 or -mat_bus_id_.size() here ????
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(ref_mat.nonZeros());
    for (size_t k=0; k < nb_bus_solver; ++k){
        if(mat_bus_id_(k) == -1) continue;  // I don't add anything to the slack bus
        for (typename Eigen::SparseMatrix<ref_mat_type>::InnerIterator it(ref_mat, k); it; ++it)
        {
            size_t row_res = static_cast<size_t>(it.row());  // TODO Eigen::Index here ?
            row_res = mat_bus_id_(row_res);
            size_t col_res = static_cast<size_t>(it.col());  // should be k   // TODO Eigen::Index here ?
            col_res = mat_bus_id_(col_res);
            if(row_res == -1) continue;
            if(col_res == -1) continue;
            tripletList.push_back(Eigen::Triplet<real_type> (row_res, col_res, std::real(it.value())));
        }
    }
    res_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    res_mat.makeCompressed();
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
RealMat BaseDCAlgo<LinearSolver>::get_lodf(const IntVect & from_bus,
                                           const IntVect & to_bus){
    auto timer = CustTimer();
    const RealMat PTDF = get_ptdf();  // size n_line x n_bus
    RealMat LODF = RealMat::Zero(from_bus.size(), from_bus.rows());  // nb_line, nb_line
    const real_type tol_equal_float = _tol_equal_float;
    for(size_t line_id=0; line_id < from_bus.size(); ++line_id){
        auto f_bus = from_bus(line_id);
        auto t_bus = to_bus(line_id);
        if ((f_bus == BaseConstants::_deactivated_bus_id) || (t_bus == BaseConstants::_deactivated_bus_id)){
            // element is disconnected
            LODF.col(line_id).array() = std::numeric_limits<real_type>::quiet_NaN();
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
