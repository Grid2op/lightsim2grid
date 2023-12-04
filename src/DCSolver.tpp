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
bool BaseDCSolver<LinearSolver>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
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

    auto timer = CustTimer();
    BaseSolver::reset_timer();
    const int nb_bus_solver = static_cast<int>(Ybus.rows());

    #ifdef __COUT_TIMES
        auto timer_preproc = CustTimer();
    #endif // __COUT_TIMES

    const CplxVect & Sbus_tmp = Sbus;

    // TODO SLACK (for now i put all slacks as PV, except the first one)
    // this should be handled in Sbus, because we know the amount of power absorbed by the slack
    // so we can compute it correctly !
    my_pv_ = retrieve_pv_with_slack(slack_ids, pv);
    // const Eigen::VectorXi & my_pv = pv;

    // find the slack buses
    slack_buses_ids_solver_ = extract_slack_bus_id(my_pv_, pq, nb_bus_solver);

    // corresp bus -> solverbus
    fill_mat_bus_id(nb_bus_solver);

    // remove the slack bus from Ybus
    // and extract only real part
    fill_dcYbus_noslack(nb_bus_solver, Ybus);
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc: preproc: " << 1000. * timer_preproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES

    // initialize the solver (only if needed)
    #ifdef __COUT_TIMES
        auto timer_solve = CustTimer();
    #endif // __COUT_TIMES
    bool just_factorize = false;
    if(need_factorize_){
        ErrorType status_init = _linear_solver.initialize(dcYbus_noslack_);
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            return false;
        }
        need_factorize_ = false;
        just_factorize = true;
    }

    // remove the slack bus from Sbus
    dcSbus_noslack_ = RealVect::Constant(nb_bus_solver - slack_buses_ids_solver_.size(), my_zero_);
    for (int k=0; k < nb_bus_solver; ++k){
        if(mat_bus_id_(k) == -1) continue;  // I don't add anything to the slack bus
        const int col_res = mat_bus_id_(k);
        dcSbus_noslack_(col_res) = std::real(Sbus_tmp(k));
    }

    // solve for theta: Sbus = dcY . theta (make a copy to keep dcSbus_noslack_)
    RealVect Va_dc_without_slack = dcSbus_noslack_;
    ErrorType error = _linear_solver.solve(dcYbus_noslack_, Va_dc_without_slack, just_factorize);
    if(error != ErrorType::NoError){
        err_ = error;
        timer_total_nr_ += timer.duration();
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
        return false;
    }

    #ifdef __COUT_TIMES
        std::cout << "\t dc solve: " << 1000. * timer_solve.duration() << "ms" << std::endl;
        auto timer_postproc = CustTimer();
    #endif // __COUT_TIMES

    // retrieve back the results in the proper shape (add back the slack bus)
    // TODO have a better way for this, for example using `.segment(0,npv)`
    // see the BaseSolver.cpp: _evaluate_Fx
    RealVect Va_dc = RealVect::Constant(nb_bus_solver, my_zero_);
    // fill Va from dc approx
    for (int ybus_id=0; ybus_id < nb_bus_solver; ++ybus_id){
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

    timer_total_nr_ += timer.duration();
    return true;
}

template<class LinearSolver>
void BaseDCSolver<LinearSolver>::fill_mat_bus_id(int nb_bus_solver){
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
void BaseDCSolver<LinearSolver>::fill_dcYbus_noslack(int nb_bus_solver, const Eigen::SparseMatrix<cplx_type> & ref_mat){
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    remove_slack_buses(nb_bus_solver, ref_mat, dcYbus_noslack_);
}

template<class LinearSolver>
template<typename ref_mat_type>  // ref_mat_type should be `real_type` or `cplx_type`
void BaseDCSolver<LinearSolver>::remove_slack_buses(int nb_bus_solver, const Eigen::SparseMatrix<ref_mat_type> & ref_mat, Eigen::SparseMatrix<real_type> & res_mat){
    res_mat = Eigen::SparseMatrix<real_type>(nb_bus_solver - 1, nb_bus_solver - 1);  // TODO dist slack: -1 or -mat_bus_id_.size() here ????
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
void BaseDCSolver<LinearSolver>::reset(){
    BaseSolver::reset();
    _linear_solver.reset();
    need_factorize_ = true;
    dcSbus_noslack_ = RealVect();
    dcYbus_noslack_ = Eigen::SparseMatrix<real_type>();
    my_pv_ = Eigen::VectorXi();
    slack_buses_ids_solver_ = Eigen::VectorXi();
    mat_bus_id_ = Eigen::VectorXi();
}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCSolver<LinearSolver>::get_ptdf(){
    // TODO
    // Bf (nb_branch, nb_bus) : en dc un truc du genre 1 / x / tap for (1..nb_branch, from_bus)
    // and -1. / x / tap for (1..nb_branch, to_bus) 
    return dcYbus_noslack_;
}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCSolver<LinearSolver>::get_lodf(){
    // TODO
    return dcYbus_noslack_;

}

template<class LinearSolver>
Eigen::SparseMatrix<real_type> BaseDCSolver<LinearSolver>::get_bsdf(){
    // TODO
    return dcYbus_noslack_;

}
