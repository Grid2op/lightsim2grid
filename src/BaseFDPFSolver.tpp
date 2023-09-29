// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// inspired from pypower https://github.com/rwl/PYPOWER/blob/master/pypower/fdpf.py

template<class LinearSolver, FDPFMethod XB_BX>
bool BaseFDPFSolver<LinearSolver, XB_BX>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
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
    /**
    This method uses the newton raphson algorithm to compute voltage angles and magnitudes at each bus
    of the system.
    If the Ybus matrix changed, please uses the appropriate method to recomptue it!
    **/
    // TODO DEBUG MODE check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
    // TODO DEBUG MODE Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
    // TODO DEBUG MODE: check that all nodes is either pv, pq or slack
    // TODO DEBUG MODE: check that the slack_weights sum to 1.0
    if(Sbus.size() != Ybus.rows() || Sbus.size() != Ybus.cols() ){
        // TODO DEBUG MODE
        std::ostringstream exc_;
        exc_ << "BaseFDPFSolver::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.rows() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(V.size() != Ybus.rows() || V.size() != Ybus.cols() ){
        // TODO DEBUG MODE
        std::ostringstream exc_;
        exc_ << "BaseFDPFSolver::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows()<< ", " << Ybus.rows() << ").";
        throw std::runtime_error(exc_.str());
    }
    reset_timer();
    auto timer = CustTimer();
    if(!is_linear_solver_valid()) return false;

    err_ = ErrorType::NoError;  // reset the error if previous error happened

    Eigen::VectorXi my_pv = retrieve_pv_with_slack(slack_ids, pv);  // retrieve_pv_with_slack (not all), add_slack_to_pv (all)
    real_type slack_absorbed = std::real(Sbus.sum());  // initial guess for slack_absorbed
    const auto slack_bus_id = slack_ids(0);
    
    // initialize once and for all the "inverse" of these vectors
    const auto n_pv = my_pv.size();
    const auto n_pq = pq.size();
    Eigen::VectorXi pvpq(n_pv + n_pq);
    pvpq << my_pv, pq;  // pvpq = r_[pv, pq]

    // TODO FDPF inheritance (or specialization for FDXB or FDBX)
    Eigen::SparseMatrix<real_type> Bp;
    Eigen::SparseMatrix<real_type> Bpp;
    fillBp(Bp);
    fillBpp(Bpp);

    // TODO FDPF 
    // Bp = Bp[array([pvpq]).T, pvpq].tocsc() # splu requires a CSC matrix
    // Bpp = Bpp[array([pq]).T, pq].tocsc()
    // then init the solvers !  TODO FDPF 
    
    // some clever tricks are used in the making of the Jacobian to handle the slack bus 
    // (in case there is a distributed slack bus)
    const auto n_pvpq = pvpq.size();
    std::vector<int> pvpq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
    std::vector<int> pq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

    V_ = V;  // V = V0
    Vm_ = V_.array().abs();   // Vm = abs(V)
    Va_ = V_.array().arg();  // Va = angle(V)

    // first check, if the problem is already solved, i stop there
    // compute a first time the mismatch to initialize the slack bus
    CplxVect mis = evaluate_mismatch(Ybus, V, Sbus, slack_bus_id, slack_absorbed, slack_weights);
    p_ = mis.real()(pvpq);  // P = mis[pvpq].real
    q_ = mis.imag()(pq);  // Q = mis[pq].imag
    
    CplxVect tmp_va;

    bool converged = _check_for_convergence(p_, q_, tol);
    nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    bool ls_bp_has_just_been_initialized = false; // TODO !
    bool ls_bpp_has_just_been_initialized = false; // TODO !
    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;

        // do the P iteration (for Va)
        solve(_linear_solver_Bp, Bp_, p_, ls_bp_has_just_been_initialized);  //  dVa = -Bp_solver.solve(P)        
        if(err_ != ErrorType::NoError){
            // I got an error during the solving of the linear system, i need to stop here
            res = false;
            break;
        }
        Va_(pvpq) -= p_;  // Va[pvpq] = Va[pvpq] + dVa
        tmp_va.array() = (Va_.array().cos().template cast<cplx_type>() + my_i * Va_.array().sin().template cast<cplx_type>() );  // reused for Q iteration

        if(has_converged(tmp_va, Ybus, Sbus, slack_bus_id, slack_absorbed, slack_weights, pvpq, pq, tol)) break;

        // do the Q iterations (for Vm)
        solve(_linear_solver_Bpp, Bpp_, q_, ls_bpp_has_just_been_initialized);  //  dVm = -Bpp_solver.solve(Q)   
        if(err_ != ErrorType::NoError){
            // I got an error during the solving of the linear system, i need to stop here
            res = false;
            break;
        }
        Vm_(pq) -= q_;
        converged = has_converged(tmp_va, Ybus, Sbus, slack_bus_id, slack_absorbed, slack_weights, pvpq, pq, tol);
    }
    if(!converged){
        if (err_ == ErrorType::NoError) err_ = ErrorType::TooManyIterations;
        res = false;
    }
    timer_total_nr_ += timer.duration();
    #ifdef __COUT_TIMES
        std::cout << "Computation time: " << "\n\t timer_initialize_: " << timer_initialize_
                  << "\n\t timer_dSbus_ (called in _fillJ_): " << timer_dSbus_
                  << "\n\t timer_fillJ_: " << timer_fillJ_
                  << "\n\t timer_Fx_: " << timer_Fx_
                  << "\n\t timer_check_: " << timer_check_
                  << "\n\t timer_solve_: " << timer_solve_
                  << "\n\t timer_total_nr_: " << timer_total_nr_
                  << "\n\n";
    #endif // __COUT_TIMES
    return res;
}
