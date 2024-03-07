// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// #include "DCSolver.h"

// TODO SLACK !!!
template<class LinearSolver>
bool BaseDCAlgo<LinearSolver>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                          CplxVect & V,
                                          const CplxVect & Sbus,
                                          const Eigen::VectorXi & slack_ids,
                                          const RealVect & slack_weights,
                                          const Eigen::VectorXi & pv,
                                          const Eigen::VectorXi & pq,
                                          int max_iter,
                                          real_type tol
                                          )
{
    // max_iter is ignored
    // tol is ignored
    // V is used the following way: at pq buses it's completely ignored. For pv bus only the magnitude is used,
    //   and for the slack bus both the magnitude and the angle are used.

    if(!is_linear_solver_valid()) {
        // std::cout << "!is_linear_solver_valid()\n";
        return false;
    }
    BaseAlgo::reset_timer();
    bool doesnt_need_refactor = true;

    auto timer = CustTimer();
    if(need_factorize_ ||
       _solver_control.need_reset_solver() || 
       _solver_control.has_dimension_changed() ||
       _solver_control.has_slack_participate_changed() ||  // the full "ybus without slack" has changed, everything needs to be recomputed_solver_control.ybus_change_sparsity_pattern()
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero()
       ){
       reset();
    }
    
    sizeYbus_with_slack_ = static_cast<int>(Ybus.rows());

    #ifdef __COUT_TIMES
        auto timer_preproc = CustTimer();
    #endif // __COUT_TIMES

    if(need_factorize_ || 
       _solver_control.has_pv_changed() || 
       _solver_control.has_pq_changed()) {    

        // TODO SLACK (for now i put all slacks as PV, except the first one)
        // this should be handled in Sbus, because we know the amount of power absorbed by the slack
        // so we can compute it correctly !
        // std::cout << "\tneed to retrieve slack\n";
        my_pv_ = retrieve_pv_with_slack(slack_ids, pv);

        // find the slack buses
        slack_buses_ids_solver_ = extract_slack_bus_id(my_pv_, pq, sizeYbus_with_slack_);
        sizeYbus_without_slack_ = sizeYbus_with_slack_ - slack_buses_ids_solver_.size();

        // corresp bus -> solverbus
        fill_mat_bus_id(sizeYbus_with_slack_);
    }

    // remove the slack bus from Ybus
    if(need_factorize_ || 
       _solver_control.need_recompute_ybus() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero()) {
        // std::cout << "\tneed to change Ybus\n";
        fill_dcYbus_noslack(sizeYbus_with_slack_, Ybus);
        doesnt_need_refactor = false;  // force a call to "factor" the linear solver as the lhs (ybus) changed
        // no need to refactor if ybus did not change
    }
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc: preproc: " << 1000. * timer_preproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES
    
    // initialize the solver (only if needed)
    #ifdef __COUT_TIMES
        auto timer_solve = CustTimer();
    #endif // __COUT_TIMES

    // remove the slack bus from Sbus
    if(need_factorize_ || _solver_control.need_recompute_sbus()){
        // std::cout << "\tneed to compute Sbus\n";
        dcSbus_noslack_ = RealVect::Constant(sizeYbus_without_slack_, my_zero_);
        for (int k=0; k < sizeYbus_with_slack_; ++k){
            if(mat_bus_id_(k) == -1) continue;  // I don't add anything to the slack bus
            const int col_res = mat_bus_id_(k);
            dcSbus_noslack_(col_res) = std::real(Sbus(k));
        }
    }

    // initialize the solver if needed
    if(need_factorize_){
        // std::cout << "\tneed to factorize\n";
        ErrorType status_init = _linear_solver.initialize(dcYbus_noslack_);
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            // std::cout << "_linear_solver.initialize\n";
            return false;
        }
        need_factorize_ = false;
        doesnt_need_refactor = true;
    }

    // solve for theta: Sbus = dcY . theta (make a copy to keep dcSbus_noslack_)
    RealVect Va_dc_without_slack = dcSbus_noslack_;
    ErrorType error = _linear_solver.solve(dcYbus_noslack_, Va_dc_without_slack, doesnt_need_refactor);
    if(error != ErrorType::NoError){
        err_ = error;
        timer_total_nr_ += timer.duration();
            // std::cout << "_linear_solver.solve\n";
        return false;
    }

    if(!Va_dc_without_slack.array().allFinite() || (Va_dc_without_slack.lpNorm<Eigen::Infinity>() >= 1e6)){
        // for convergence, all values should be finite
        // and it's not realistic if some Va are too high
        err_ = ErrorType::SolverSolve;
        V = CplxVect();
        V_ = CplxVect();
        Vm_ = RealVect();
        Va_ = RealVect();
        timer_total_nr_ += timer.duration();
        // std::cout << "_linear_solver.allFinite" << Va_dc_without_slack.array().allFinite() <<", " << Va_dc_without_slack.lpNorm<Eigen::Infinity>() <<"\n";
        return false;
    }
    // std::cout << "\t " << Va_dc_without_slack.lpNorm<Eigen::Infinity>() << "\n";

    #ifdef __COUT_TIMES
        std::cout << "\t dc solve: " << 1000. * timer_solve.duration() << "ms" << std::endl;
        auto timer_postproc = CustTimer();
    #endif // __COUT_TIMES

    // retrieve back the results in the proper shape (add back the slack bus)
    // TODO have a better way for this, for example using `.segment(0,npv)`
    // see the BaseAlgo.cpp: _evaluate_Fx
    RealVect Va_dc = RealVect::Constant(sizeYbus_with_slack_, my_zero_);
    // fill Va from dc approx
    for (int ybus_id=0; ybus_id < sizeYbus_with_slack_; ++ybus_id){
        if(mat_bus_id_(ybus_id) == -1) continue;  // slack bus is handled elsewhere
        const int bus_me = mat_bus_id_(ybus_id);
        Va_dc(ybus_id) = Va_dc_without_slack(bus_me);
    }
    Va_dc.array() += std::arg(V(slack_buses_ids_solver_(0)));

    // save the results
    Va_ = Va_dc;

    // add the Voltage setpoints of the generator
    Vm_ = V.array().abs();
    Vm_(slack_buses_ids_solver_) = V(slack_buses_ids_solver_).array().abs();

    // now compute the resulting complex voltage
    V_ = (Va_.array().cos().template cast<cplx_type>() + my_i * Va_.array().sin().template cast<cplx_type>());

    V_.array() *= Vm_.array();
    nr_iter_ = 1;
    V = V_;
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc postproc: " << 1000. * timer_postproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES
    _solver_control.tell_none_changed();
    timer_total_nr_ += timer.duration();
    return true;
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::fill_mat_bus_id(int nb_bus_solver){
    mat_bus_id_ = Eigen::VectorXi::Constant(nb_bus_solver, -1);
    // Eigen::VectorXi me_to_ybus = Eigen::VectorXi::Constant(nb_bus_solver - slack_bus_ids_solver.size(), -1);
    int solver_id = 0;
    for (int ybus_id=0; ybus_id < nb_bus_solver; ++ybus_id){
        if(isin(ybus_id, slack_buses_ids_solver_)) continue;  // I don't add anything to the slack bus
        mat_bus_id_(ybus_id) = solver_id;
        // me_to_ybus(solver_id) = ybus_id;
        ++solver_id;
    }
}

template<class LinearSolver>
void BaseDCAlgo<LinearSolver>::fill_dcYbus_noslack(int nb_bus_solver, const Eigen::SparseMatrix<cplx_type> & ref_mat){
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    remove_slack_buses(nb_bus_solver, ref_mat, dcYbus_noslack_);
}

template<class LinearSolver>
template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
void BaseDCAlgo<LinearSolver>::remove_slack_buses(int nb_bus_solver, const Eigen::SparseMatrix<ref_mat_type> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat){
    res_mat = Eigen::SparseMatrix<real_type>(sizeYbus_without_slack_, sizeYbus_without_slack_);  // TODO dist slack: -1 or -mat_bus_id_.size() here ????
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(ref_mat.nonZeros());
    for (int k=0; k < nb_bus_solver; ++k){
        if(mat_bus_id_(k) == -1) continue;  // I don't add anything to the slack bus
        for (typename Eigen::SparseMatrix<ref_mat_type>::InnerIterator it(ref_mat, k); it; ++it)
        {
            int row_res = static_cast<int>(it.row());  // TODO Eigen::Index here ?
            if(mat_bus_id_(row_res) == -1) continue;
            row_res = mat_bus_id_(row_res);
            int col_res = static_cast<int>(it.col());  // should be k   // TODO Eigen::Index here ?
            col_res = mat_bus_id_(col_res);
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
    sizeYbus_with_slack_ = 0;
    sizeYbus_without_slack_ = 0;
    dcSbus_noslack_ = RealVect();
    dcYbus_noslack_ = Eigen::SparseMatrix<real_type>();
    my_pv_ = Eigen::VectorXi();
    slack_buses_ids_solver_ = Eigen::VectorXi();
    mat_bus_id_ = Eigen::VectorXi();
}

template<class LinearSolver>
RealMat BaseDCAlgo<LinearSolver>::get_ptdf(const Eigen::SparseMatrix<cplx_type> & dcYbus){
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
        if(mat_bus_id_(bus_id) == -1) continue;
        ind_no_slack_.push_back(bus_id);
    }
    const Eigen::VectorXi ind_no_slack = Eigen::VectorXi::Map(&ind_no_slack_[0], ind_no_slack_.size());

    // solve iteratively the linear systems (one per powerline)
    PTDF = RealMat::Zero(Bf_T_with_slack.cols(), Bf_T_with_slack.rows());  // rows and cols are "inverted" because the matrix Bf is transposed
    for (int line_id=0; line_id < nb_pow_tr; ++line_id){
        // build the rhs vector
        for (typename Eigen::SparseMatrix<real_type>::InnerIterator it(Bf_T_with_slack, line_id); it; ++it)
        {
            const auto bus_id = it.row();
            if(mat_bus_id_(bus_id) == -1) continue;  // I don't add anything if it's the slack
            const auto col_res = mat_bus_id_(bus_id);
            rhs[col_res] = it.value();
        }
        // solve the linear system
        _linear_solver.solve(dcYbus_noslack_, rhs, true);  // I don't need to refactorize the matrix (hence the `true`)

        // assign results to the PTDF matrix
        PTDF(line_id, ind_no_slack) = rhs;

        // reset the rhs vector to 0.
        rhs.array() = 0.;
        // rhs = RealVect::Zero(sizeYbus_without_slack_);
    }
    // TODO PTDF: if the solver can solve the MAT directly, do that instead
    return PTDF;
}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCAlgo<LinearSolver>::get_lodf(){
    // TODO
    return dcYbus_noslack_;

}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCAlgo<LinearSolver>::get_bsdf(){
    // TODO
    return dcYbus_noslack_;

}
