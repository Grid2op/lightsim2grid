// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// inspired from pypower https://github.com/rwl/PYPOWER/blob/master/pypower/fdpf.py

template<class LinearSolver, FDPFMethod XB_BX>
bool BaseFDPFAlgo<LinearSolver, XB_BX>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
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
        exc_ << "BaseFDPFAlgo::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(V.size() != Ybus.rows() || V.size() != Ybus.cols() ){
        // TODO DEBUG MODE
        std::ostringstream exc_;
        exc_ << "BaseFDPFAlgo::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows()<< ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    reset_timer();
    if(!is_linear_solver_valid()) return false;
    if(_solver_control.need_reset_solver() || 
       _solver_control.has_dimension_changed() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero() ||
       _solver_control.has_slack_participate_changed() ||
       _solver_control.has_pv_changed() ||
       _solver_control.has_pq_changed()){
       reset();
    }

    err_ = ErrorType::NoError;  // reset the error if previous error happened
    auto timer = CustTimer();

    Eigen::VectorXi my_pv = retrieve_pv_with_slack(slack_ids, pv);  // retrieve_pv_with_slack (not all), add_slack_to_pv (all)
    real_type slack_absorbed = std::real(Sbus.sum());  // initial guess for slack_absorbed
    const auto slack_bus_id = slack_ids(0);
    
    // initialize once and for all the "inverse" of these vectors
    const auto n_pv = my_pv.size();
    const auto n_pq = pq.size();
    Eigen::VectorXi pvpq(n_pv + n_pq);
    pvpq << my_pv, pq;  // pvpq = r_[pv, pq]
    const auto n_pvpq = pvpq.size();

    // fill the sparse matrix Bp and Bpp (depends on the method: BX or XB)
    if(need_factorize_ ||
       _solver_control.need_reset_solver() || 
       _solver_control.has_dimension_changed() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero() ||
       _solver_control.need_recompute_ybus()
       ){
        // need to extract Bp and Bpp for the whole grid
        grid_Bp_ = Eigen::SparseMatrix<real_type> ();
        grid_Bpp_ = Eigen::SparseMatrix<real_type>(); 
        fillBp_Bpp(grid_Bp_, grid_Bpp_);
    }

    // init "my" matrices
    // fill the solver matrices Bp_ and Bpp_
    // Bp_ = Bp[array([pvpq]).T, pvpq].tocsc()
    // Bpp_ = Bpp[array([pq]).T, pq].tocsc()
    if(need_factorize_ ||
       _solver_control.need_reset_solver() || 
       _solver_control.has_dimension_changed() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero() ||
       _solver_control.need_recompute_ybus() ||
       _solver_control.has_slack_participate_changed() ||
       _solver_control.has_pv_changed() ||
       _solver_control.has_pq_changed()
       ){            
        std::vector<int> pvpq_inv(V.size(), -1);
        for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
        std::vector<int> pq_inv(V.size(), -1);
        for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
        fill_sparse_matrices(grid_Bp_, grid_Bpp_, pvpq_inv, pq_inv, n_pvpq, n_pq);
    }

    V_ = V;  // V = V0
    Vm_ = V_.array().abs();   // Vm = abs(V)
    Va_ = V_.array().arg();  // Va = angle(V)

    // first check, if the problem is already solved, i stop there
    // compute a first time the mismatch to initialize the slack bus
    CplxVect mis = evaluate_mismatch(Ybus, V, Sbus, slack_bus_id, slack_absorbed, slack_weights);
    mis.array() /= Vm_.array();  // mis = (V * conj(Ybus * V) - Sbus) / Vm
    p_ = mis(pvpq).real();  // P = mis[pvpq].real
    q_ = mis(pq).imag();  // Q = mis[pq].imag
    
    CplxVect tmp_va;
    nr_iter_ = 0; //current step
    bool converged = _check_for_convergence(p_, q_, tol);
    bool res = true;  // have i converged or not
    bool ls_bp_has_just_been_initialized = false;
    bool ls_bpp_has_just_been_initialized = false;

    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;

        if(need_factorize_){
            initialize();
            if(err_ != ErrorType::NoError){
                // I got an error during the initialization of the linear system, i need to stop here
                res = false;
                break;
            }
            ls_bp_has_just_been_initialized = true;
            ls_bpp_has_just_been_initialized = true;
        }

        // do the P iteration (for Va)
        RealVect x = p_;
        solve(_linear_solver_Bp, Bp_, x, ls_bp_has_just_been_initialized);  //  dVa = -Bp_solver.solve(P)     
        if(err_ != ErrorType::NoError){
            // I got an error during the solving of the linear system, i need to stop here
            res = false;
            break;
        }
        Va_(pvpq) -= x;  // Va[pvpq] = Va[pvpq] + dVa
        tmp_va.array() = (Va_.array().cos().template cast<cplx_type>() + my_i * Va_.array().sin().template cast<cplx_type>() );  // reused for Q iteration
        if(has_converged(tmp_va, Ybus, Sbus, slack_bus_id, slack_absorbed, slack_weights, pvpq, pq, tol)){
            converged = true;
            break;
        }
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

template<class LinearSolver, FDPFMethod XB_BX>
void BaseFDPFAlgo<LinearSolver, XB_BX>::fill_sparse_matrices(const Eigen::SparseMatrix<real_type> & grid_Bp,
                                                             const Eigen::SparseMatrix<real_type> & grid_Bpp,
                                                             const std::vector<int> & pvpq_inv,
                                                             const std::vector<int> & pq_inv,
                                                             Eigen::Index n_pvpq,
                                                             Eigen::Index n_pq)
{
  /**
   Init Bp_ and Bpp_ such that:
    // Bp_ = grid_Bp[array([pvpq]).T, pvpq].tocsc() # splu requires a CSC matrix
    // Bpp_ = grid_Bpp[array([pq]).T, pq].tocsc()
   **/
  aux_fill_sparse_matrices(grid_Bp, pvpq_inv, n_pvpq, Bp_);
  aux_fill_sparse_matrices(grid_Bpp, pq_inv, n_pq, Bpp_);
}

template<class LinearSolver, FDPFMethod XB_BX>
void BaseFDPFAlgo<LinearSolver, XB_BX>::aux_fill_sparse_matrices(const Eigen::SparseMatrix<real_type> & grid_Bp_Bpp,
                                                                   const std::vector<int> & ind_inv,
                                                                   Eigen::Index mat_dim,
                                                                   Eigen::SparseMatrix<real_type> & res)
{
    // clear previous matrix
    res = Eigen::SparseMatrix<real_type>(mat_dim, mat_dim);
    // compute the coefficients
    std::vector<Eigen::Triplet<real_type> > tripletList;
    tripletList.reserve(grid_Bp_Bpp.nonZeros());
    for (int k = 0; k < grid_Bp_Bpp.outerSize(); ++k){
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(grid_Bp_Bpp, k); it; ++it){
            if ((ind_inv[it.row()] >= 0) && (ind_inv[it.col()] >= 0)){
                auto tmp = Eigen::Triplet<real_type>(ind_inv[it.row()], ind_inv[it.col()], it.value());
                tripletList.push_back(tmp);
            }
        }
    }
    // put them in the matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    res.makeCompressed();
}
