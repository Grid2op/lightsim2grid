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
    Eigen::SparseMatrix<real_type> dcYbus = Eigen::SparseMatrix<real_type>(nb_bus_solver - 1, nb_bus_solver - 1);

    const CplxVect & Sbus_tmp = Sbus;

    // TODO SLACK (for now i put all slacks as PV, except the first one)
    // this should be handled in Sbus, because we know the amount of power absorbed by the slack
    // so we can compute it correctly !
    Eigen::VectorXi my_pv = retrieve_pv_with_slack(slack_ids, pv);

    // find the slack bus
    int slack_bus_id_solver = extract_slack_bus_id(my_pv, pq, nb_bus_solver);

    // remove the slack bus from Ybus
    // and extract only real part
    // TODO see if "prune" might work here https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title29
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(Ybus.nonZeros());
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        for (Eigen::SparseMatrix<cplx_type>::InnerIterator it(Ybus, k); it; ++it)
        {
            int row_res = static_cast<int>(it.row());
            if(row_res == slack_bus_id_solver) continue;
            row_res = row_res > slack_bus_id_solver ? row_res - 1 : row_res;
            int col_res = static_cast<int>(it.col());
            col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
            tripletList.push_back(Eigen::Triplet<real_type> (row_res, col_res, std::real(it.value())));
        }
    }
    dcYbus.setFromTriplets(tripletList.begin(), tripletList.end());
    dcYbus.makeCompressed();
    
    #ifdef __COUT_TIMES
        std::cout << "\t dc: preproc: " << 1000. * timer_preproc.duration() << "ms" << std::endl;
    #endif // __COUT_TIMES

    // initialize the solver (only if needed)
    #ifdef __COUT_TIMES
        auto timer_solve = CustTimer();
    #endif // __COUT_TIMES
    bool just_factorize = false;
    if(need_factorize_){
        ErrorType status_init = _linear_solver.initialize(dcYbus);
        if(status_init != ErrorType::NoError){
            err_ = status_init;
            return false;
        }
        need_factorize_ = false;
        just_factorize = true;
    }

    // remove the slack bus from Sbus
    RealVect dcSbus = RealVect::Constant(nb_bus_solver - 1, my_zero_);
    for (int k=0; k < nb_bus_solver; ++k){
        if(k == slack_bus_id_solver) continue;  // I don't add anything to the slack bus
        int col_res = k;
        col_res = col_res > slack_bus_id_solver ? col_res - 1 : col_res;
        dcSbus(col_res) = std::real(Sbus_tmp(k));
    }

    // solve for theta: Sbus = dcY . theta
    RealVect Va_dc_without_slack = dcSbus;

    ErrorType error = _linear_solver.solve(dcYbus, Va_dc_without_slack, just_factorize);
    if(error != ErrorType::NoError){
        err_ = error;
        timer_total_nr_ += timer.duration();
        return false;
    }

    if(!Va_dc_without_slack.array().allFinite()){
        // for convergence, all values should be finite
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
    for (int bus_id_with_slack=0; bus_id_with_slack < nb_bus_solver; ++bus_id_with_slack){
        if(bus_id_with_slack == slack_bus_id_solver) continue;  // slack bus is handled elsewhere
        int bus_id_without_slack = bus_id_with_slack > slack_bus_id_solver ? bus_id_with_slack - 1 : bus_id_with_slack;
        Va_dc(bus_id_with_slack) = Va_dc_without_slack(bus_id_without_slack);
    }
    Va_dc.array() += std::arg(V(slack_bus_id_solver));

    // save the results
    Va_ = Va_dc;

    // add the Voltage setpoints of the generator
    Vm_ = V.array().abs();
    Vm_(slack_bus_id_solver) = std::abs(V(slack_bus_id_solver));

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
void BaseDCSolver<LinearSolver>::reset(){
    BaseSolver::reset();
    _linear_solver.reset();
    need_factorize_ = true;
}
