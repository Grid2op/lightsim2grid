// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// #include "BaseNRSolverSingleSlack.h"
// #include "BaseNRAlgo.h"

template<class LinearSolver>
bool BaseNRSingleSlackAlgo<LinearSolver>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                                     CplxVect & V,
                                                     const CplxVect & Sbus,
                                                     const Eigen::VectorXi & slack_ids,
                                                     const RealVect & slack_weights,  // unused here
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
    // TODO check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
    // TODO Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
    if(Sbus.size() != Ybus.rows() || Sbus.size() != Ybus.cols() ){
        std::ostringstream exc_;
        exc_ << "BaseNRSingleSlackAlgo::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(V.size() != Ybus.rows() || V.size() != Ybus.cols() ){
        std::ostringstream exc_;
        exc_ << "BaseNRSingleSlackAlgo::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows()<<", "<<Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(!BaseNRAlgo<LinearSolver>::is_linear_solver_valid()){
        return false;
    }
    BaseNRAlgo<LinearSolver>::reset_timer();
    BaseNRAlgo<LinearSolver>::reset_if_needed();
    BaseNRAlgo<LinearSolver>::err_ = ErrorType::NoError;  // reset the error if previous error happened
    
    auto timer = CustTimer();
    auto timer_pre_proc = CustTimer();
    // initialize once and for all the "inverse" of these vectors
    // Eigen::VectorXi my_pv = BaseNRAlgo<LinearSolver>::retrieve_pv_with_slack(slack_ids, pv);
    Eigen::VectorXi my_pv = pv;
    // Eigen::VectorXi my_pv = pv; // BaseNRAlgo<LinearSolver>::retrieve_pv_with_slack(slack_ids, pv);

    const int n_pv = static_cast<int>(my_pv.size());
    const int n_pq = static_cast<int>(pq.size());
    Eigen::VectorXi pvpq(n_pv + n_pq);
    pvpq << my_pv, pq; 
    const int n_pvpq = static_cast<int>(pvpq.size());
    std::vector<int> pvpq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
    std::vector<int> pq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

    BaseNRAlgo<LinearSolver>::V_ = V;
    BaseNRAlgo<LinearSolver>::Vm_ = BaseNRAlgo<LinearSolver>::V_.array().abs();  // update Vm and Va again in case
    BaseNRAlgo<LinearSolver>::Va_ = BaseNRAlgo<LinearSolver>::V_.array().arg();  // we wrapped around with a negative Vm
    BaseNRAlgo<LinearSolver>::timer_pre_proc_ += timer_pre_proc.duration();

    // first check, if the problem is already solved, i stop there
    RealVect F = BaseNRAlgo<LinearSolver>::_evaluate_Fx(Ybus, V, Sbus, my_pv, pq);
    bool converged = BaseNRAlgo<LinearSolver>::_check_for_convergence(F, tol);
    BaseNRAlgo<LinearSolver>::nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    bool has_just_been_initialized = false;  // to avoid a call to klu_refactor follow a call to klu_factor in the same loop

    const cplx_type m_i = BaseNRAlgo<LinearSolver>::my_i;  // otherwise it does not compile
    if(BaseNRAlgo<LinearSolver>::need_factorize_ ||
       BaseNRAlgo<LinearSolver>::_solver_control.need_reset_solver() || 
       BaseNRAlgo<LinearSolver>::_solver_control.has_dimension_changed() ||
       BaseNRAlgo<LinearSolver>::_solver_control.has_slack_participate_changed() ||  // the full "ybus without slack" has changed, everything needs to be recomputed_solver_control.ybus_change_sparsity_pattern()
       BaseNRAlgo<LinearSolver>::_solver_control.ybus_change_sparsity_pattern() ||
       BaseNRAlgo<LinearSolver>::_solver_control.has_ybus_some_coeffs_zero() ||
       BaseNRAlgo<LinearSolver>::_solver_control.need_recompute_ybus() ||
    // BaseNRAlgo<LinearSolver>::   _solver_control.has_slack_participate_changed() ||
       BaseNRAlgo<LinearSolver>::_solver_control.has_pv_changed() ||
       BaseNRAlgo<LinearSolver>::_solver_control.has_pq_changed()
       )
       {
        BaseNRAlgo<LinearSolver>::value_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::col_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::row_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        BaseNRAlgo<LinearSolver>::dS_dVm_.resize(0,0);  // TODO smarter solver: only needed if ybus has changed
        BaseNRAlgo<LinearSolver>::dS_dVa_.resize(0,0);  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::dS_dVm_.setZero();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::dS_dVa_.setZero();  // TODO smarter solver: only needed if ybus has changed

       }
    while ((!converged) & (BaseNRAlgo<LinearSolver>::nr_iter_ < max_iter)){
        BaseNRAlgo<LinearSolver>::nr_iter_++;
        // std::cout << "\tnr_iter_ " << BaseNRAlgo<LinearSolver>::nr_iter_ << std::endl;
        fill_jacobian_matrix(Ybus, BaseNRAlgo<LinearSolver>::V_, pq, pvpq, pq_inv, pvpq_inv);
        if(BaseNRAlgo<LinearSolver>::need_factorize_){
            BaseNRAlgo<LinearSolver>::initialize();
            if(BaseNRAlgo<LinearSolver>::err_ != ErrorType::NoError){
                // I got an error during the initialization of the linear system, i need to stop here
                // std::cout << BaseNRAlgo<LinearSolver>::err_ << std::endl;
                res = false;
                break;
            }
            has_just_been_initialized = true;
            // std::cout << "I just factorized" << std::endl;
        }else{
            // std::cout << "no need to factorize" << std::endl;
        }

        BaseNRAlgo<LinearSolver>::solve(F, has_just_been_initialized);

        has_just_been_initialized = false;
        if(BaseNRAlgo<LinearSolver>::err_ != ErrorType::NoError){
            // I got an error during the solving of the linear system, i need to stop here
            // std::cout << BaseNRAlgo<LinearSolver>::err_ << std::endl;
            res = false;
            break;
        }
        // auto dx = -F;
        auto timer_va_vm = CustTimer();
        BaseNRAlgo<LinearSolver>::Vm_ = BaseNRAlgo<LinearSolver>::V_.array().abs();  // update Vm and Va again in case
        BaseNRAlgo<LinearSolver>::Va_ = BaseNRAlgo<LinearSolver>::V_.array().arg();  // we wrapped around with a negative Vm

        // update voltage (this should be done consistently with "klu_solver._evaluate_Fx")
        if (n_pv > 0) BaseNRAlgo<LinearSolver>::Va_(my_pv) -= F.segment(0, n_pv);
        if (n_pq > 0){
            BaseNRAlgo<LinearSolver>::Va_(pq) -= F.segment(n_pv,n_pq);
            BaseNRAlgo<LinearSolver>::Vm_(pq) -= F.segment(n_pv+n_pq, n_pq);
        }

        // TODO change here for not having to cast all the time ... maybe
        const RealVect & Vm = BaseNRAlgo<LinearSolver>::Vm_;  // I am forced to redefine the type for it to compile properly
        const RealVect & Va = BaseNRAlgo<LinearSolver>::Va_;
        BaseNRAlgo<LinearSolver>::V_ = Vm.array() * (Va.array().cos().cast<cplx_type>() + m_i * Va.array().sin().cast<cplx_type>() );
        if(BaseNRAlgo<LinearSolver>::Vm_.minCoeff() < 0.)
        {
            // update Vm and Va again in case
            // we wrapped around with a negative Vm TODO more efficient way maybe ?
            BaseNRAlgo<LinearSolver>::Vm_ = BaseNRAlgo<LinearSolver>::V_.array().abs();  
            BaseNRAlgo<LinearSolver>::Va_ = BaseNRAlgo<LinearSolver>::V_.array().arg();  
        }
        BaseNRAlgo<LinearSolver>::timer_Va_Vm_ += timer_va_vm.duration();

        F = BaseNRAlgo<LinearSolver>::_evaluate_Fx(Ybus, BaseNRAlgo<LinearSolver>::V_, Sbus, my_pv, pq);
        bool tmp = F.allFinite();
        if(!tmp){
            BaseNRAlgo<LinearSolver>::err_ = ErrorType::InifiniteValue;
            // std::cout << BaseNRAlgo<LinearSolver>::err_ << std::endl;
            break; // divergence due to Nans
        }
        converged = BaseNRAlgo<LinearSolver>::_check_for_convergence(F, tol);
    }
    if(!converged){
        if (BaseNRAlgo<LinearSolver>::err_ == ErrorType::NoError) BaseNRAlgo<LinearSolver>::err_ = ErrorType::TooManyIterations;
        res = false;
    }
    BaseNRAlgo<LinearSolver>::timer_total_nr_ += timer.duration();
    #ifdef __COUT_TIMES
        std::cout << "Computation time: " << "\n\t timer_initialize_: " << BaseNRAlgo<LinearSolver>::timer_initialize_
                  << "\n\t timer_dSbus_ (called in _fillJ_): " << BaseNRAlgo<LinearSolver>::timer_dSbus_
                  << "\n\t timer_fillJ_: " << BaseNRAlgo<LinearSolver>::timer_fillJ_
                  << "\n\t timer_Fx_: " << BaseNRAlgo<LinearSolver>::timer_Fx_
                  << "\n\t timer_check_: " << BaseNRAlgo<LinearSolver>::timer_check_
                  << "\n\t timer_solve_: " << BaseNRAlgo<LinearSolver>::timer_solve_
                  << "\n\t timer_total_nr_: " << BaseNRAlgo<LinearSolver>::timer_total_nr_
                  << "\n\n";
    #endif // __COUT_TIMES
    // update Vm and Va again in case
    // we wrapped around with a negative Vm TODO more efficient way maybe ?
    BaseNRAlgo<LinearSolver>::Vm_ = BaseNRAlgo<LinearSolver>::V_.array().abs();  
    BaseNRAlgo<LinearSolver>::Va_ = BaseNRAlgo<LinearSolver>::V_.array().arg();  
    BaseNRAlgo<LinearSolver>::_solver_control.tell_none_changed();
    return res;
}

template<class LinearSolver>
void BaseNRSingleSlackAlgo<LinearSolver>::fill_jacobian_matrix(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                                                 const CplxVect & V,
                                                                 const Eigen::VectorXi & pq,
                                                                 const Eigen::VectorXi & pvpq,
                                                                 const std::vector<int> & pq_inv,
                                                                 const std::vector<int> & pvpq_inv
                                                                 )
{
    /**
    J has the shape
    | J11 | J12 |               | (pvpq, pvpq) | (pvpq, pq) |
    | --------- | = dimensions: | ------------------------- |
    | J21 | J22 |               |  (pq, pvpq)  | (pq, pq) |
    python implementation:
    J11 = dS_dVa[array([pvpq]).T, pvpq].real
    J12 = dS_dVm[array([pvpq]).T, pq].real
    J21 = dS_dVa[array([pq]).T, pvpq].imag
    J22 = dS_dVm[array([pq]).T, pq].imag
    **/

    BaseNRAlgo<LinearSolver>::_dSbus_dV(Ybus, V);

    auto timer = CustTimer();
    const int n_pvpq = static_cast<int>(pvpq.size());
    const int n_pq = static_cast<int>(pq.size());
    const int size_j = n_pvpq + n_pq;
    // TODO to gain a bit more time below, try to compute directly, in _dSbus_dV(Ybus, V);
    // TODO the `dS_dVa_[pvpq, pvpq]`
    // TODO so that it's easier to retrieve in the next few lines !
    if(BaseNRAlgo<LinearSolver>::J_.cols() != size_j)
    // if(true)
    {
        #ifdef __COUT_TIMES
            auto timer2 = CustTimer();
        #endif  // __COUT_TIMES
        // first time i initialized the matrix, so i need to compute its sparsity pattern
        fill_jacobian_matrix_unkown_sparsity_pattern(Ybus, V, pq, pvpq, pq_inv, pvpq_inv);
        fill_value_map(pq, pvpq, false);
        // std::cout << "\t\tfill_jacobian_matrix_unkown_sparsity_pattern" << std::endl;
        #ifdef __COUT_TIMES
            std::cout << "\t\t fill_jacobian_matrix_unkown_sparsity_pattern : " << timer2.duration() << std::endl;
        #endif  // __COUT_TIMES
    }else{
        // the sparsity pattern of J_ is already known, i can reuse it to fill it
        // properly and faster (or not...)
        #ifdef __COUT_TIMES
            auto timer3 = CustTimer();
        #endif  // __COUT_TIMES
        if (BaseNRAlgo<LinearSolver>::value_map_.size() == 0){
            // std::cout << "\t\tfill_value_map called" << std::endl;
            fill_value_map(pq, pvpq, true);
        }
        fill_jacobian_matrix_kown_sparsity_pattern(pq, pvpq);
        // std::cout << "\t\tfill_jacobian_matrix_kown_sparsity_pattern" << std::endl;
        #ifdef __COUT_TIMES
            std::cout << "\t\t fill_jacobian_matrix_kown_sparsity_pattern : " << timer3.duration() << std::endl;
        #endif  // __COUT_TIMES
    }
    BaseNRAlgo<LinearSolver>::timer_fillJ_ += timer.duration();
}

template<class LinearSolver>
void BaseNRSingleSlackAlgo<LinearSolver>::fill_jacobian_matrix_unkown_sparsity_pattern(
        const Eigen::SparseMatrix<cplx_type> & Ybus,
        const CplxVect & V,
        const Eigen::VectorXi & pq,
        const Eigen::VectorXi & pvpq,
        const std::vector<int> & pq_inv,
        const std::vector<int> & pvpq_inv
    )
{
    /**
    This functions fills the jacobian matrix when its sparsity pattern is not know in advance (typically
    the first iteration of the Newton Raphson)
    For that we need to perform relatively expensive computation from dS_dV* in order to retrieve it.
    This function is NOT optimized for speed...
    Remember:
    J has the shape
    | J11 | J12 |               | (pvpq, pvpq) | (pvpq, pq) |
    | --------- | = dimensions: | ------------------------- |
    | J21 | J22 |               |  (pq, pvpq)  | (pq, pq) |
    python implementation:
    J11 = dS_dVa[array([pvpq]).T, pvpq].real
    J12 = dS_dVm[array([pvpq]).T, pq].real
    J21 = dS_dVa[array([pq]).T, pvpq].imag
    J22 = dS_dVm[array([pq]).T, pq].imag
    **/
    bool need_insert = false;  // i optimization: i don't need to insert the coefficient in the matrix
    const int n_pvpq = static_cast<int>(pvpq.size());
    const int n_pq = static_cast<int>(pq.size());
    const int size_j = n_pvpq + n_pq;

    const Eigen::SparseMatrix<real_type> dS_dVa_r = BaseNRAlgo<LinearSolver>::dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVa_i = BaseNRAlgo<LinearSolver>::dS_dVa_.imag();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = BaseNRAlgo<LinearSolver>::dS_dVm_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_i = BaseNRAlgo<LinearSolver>::dS_dVm_.imag();

    // Method (1) seems to be faster than the others

    // optim : if the matrix was already computed, i don't initialize it, i instead reuse as much as i can
    // i can do that because the matrix will ALWAYS have the same non zero coefficients.
    // in this if, i allocate it in a "large enough" place to avoid copy when first filling it
    if(BaseNRAlgo<LinearSolver>::J_.cols() != size_j)
    {
        need_insert = true;
        BaseNRAlgo<LinearSolver>::J_ = Eigen::SparseMatrix<real_type>(size_j, size_j);
        // pre allocate a large enough matrix
        BaseNRAlgo<LinearSolver>::J_.reserve(2*(BaseNRAlgo<LinearSolver>::dS_dVa_.nonZeros() + BaseNRAlgo<LinearSolver>::dS_dVm_.nonZeros()));
        // from an experiment, outerIndexPtr is initialized, with the number of columns
        // innerIndexPtr and valuePtr are not.
    }

    // std::vector<Eigen::Triplet<double> >coeffs;  // HERE FOR PERF OPTIM (3)
    // coeffs.reserve(2*(dS_dVa_.nonZeros()+dS_dVm_.nonZeros()));  // HERE FOR PERF OPTIM (3)

    // i fill the buffer columns per columns
    int nb_obj_this_col = 0;
    std::vector<Eigen::Index> inner_index;
    std::vector<real_type> values;

    // TODO use the loop provided above (in dS) if J is already initialized
    // fill n_pvpq leftmost columns
    for(int col_id=0; col_id < n_pvpq; ++col_id){
        // reset from the previous column
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
        // and check the row order !
        BaseNRAlgo<LinearSolver>::_get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVa_r,
                      pvpq_inv, pvpq,
                      col_id,
                      0, 
                      0);
        // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
        BaseNRAlgo<LinearSolver>::_get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVa_i,
                      pq_inv, pvpq,
                      col_id,
                      n_pvpq,
                      0);

        // "efficient" insert of the element in the matrix
        for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            int row_id = inner_index[in_ind];
            if(need_insert) BaseNRAlgo<LinearSolver>::J_.insert(row_id, col_id) = values[in_ind];  // HERE FOR PERF OPTIM (1)
            else BaseNRAlgo<LinearSolver>::J_.coeffRef(row_id, col_id) = values[in_ind];  // HERE FOR PERF OPTIM (1)
            // J_.insert(row_id, col_id) = values[in_ind];  // HERE FOR PERF OPTIM (2)
            // coeffs.push_back(Eigen::Triplet<double>(row_id, col_id, values[in_ind]));   // HERE FOR PERF OPTIM (3)
        }
    }

    //TODO make same for the second part (have a funciton for previous loop)
    // fill the remaining n_pq columns
    for(int col_id=0; col_id < n_pq; ++col_id){
        // reset from the previous column
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

          // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
        // and check the row order !
        BaseNRAlgo<LinearSolver>::_get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVm_r,
                      pvpq_inv, pq,
                      col_id,
                      0,
                      0);

        // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
        BaseNRAlgo<LinearSolver>::_get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVm_i,
                      pq_inv, pq,
                      col_id,
                      n_pvpq,
                      0);

        // "efficient" insert of the element in the matrix
        for(int in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            int row_id = inner_index[in_ind];
            if(need_insert) BaseNRAlgo<LinearSolver>::J_.insert(row_id, col_id + n_pvpq) = values[in_ind];  // HERE FOR PERF OPTIM (1)
            else BaseNRAlgo<LinearSolver>::J_.coeffRef(row_id, col_id + n_pvpq) = values[in_ind];  // HERE FOR PERF OPTIM (1)
            // J_.insert(row_id, col_id + n_pvpq) = values[in_ind];  // HERE FOR PERF OPTIM (2)
            // coeffs.push_back(Eigen::Triplet<double>(row_id, col_id + n_pvpq, values[in_ind]));   // HERE FOR PERF OPTIM (3)
        }
    }
    // J_.setFromTriplets(coeffs.begin(), coeffs.end());  // HERE FOR PERF OPTIM (3)
    BaseNRAlgo<LinearSolver>::J_.makeCompressed();
}

/**
fill the value of the `value_map_` that stores pointers to the elements of
dS_dVa_ and dS_dVm_ to be used to fill J_
it requires that J_ is initialized, in compressed mode.
**/
template<class LinearSolver>
void BaseNRSingleSlackAlgo<LinearSolver>::fill_value_map(
        const Eigen::VectorXi & pq,
        const Eigen::VectorXi & pvpq,
        bool reset_J
        )
{
    const int n_pvpq = static_cast<int>(pvpq.size());
    BaseNRAlgo<LinearSolver>::value_map_.clear();
    // std::cout << "BaseNRAlgo<LinearSolver>::J_.nonZeros(): " << BaseNRAlgo<LinearSolver>::J_.nonZeros() << std::endl;
    BaseNRAlgo<LinearSolver>::value_map_.reserve(BaseNRAlgo<LinearSolver>::J_.nonZeros());

    const int n_col = static_cast<int>(BaseNRAlgo<LinearSolver>::J_.cols());
    unsigned int pos_el = 0;
    for (int col_=0; col_ < n_col; ++col_){
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(BaseNRAlgo<LinearSolver>::J_, col_); it; ++it)
        {
            const int row_id = static_cast<int>(it.row());
            const int col_id = static_cast<int>(it.col());  // it's equal to "col_"
            if(reset_J) it.valueRef() = 0.; // "forget" previous J value in this setting
            // real_type & this_el = J_x_ptr[pos_el];
            if((col_id < n_pvpq) && (row_id < n_pvpq)){
                // this is the J11 part (dS_dVa_r)
                const int row_id_dS_dVa_r = pvpq[row_id];
                const int col_id_dS_dVa_r = pvpq[col_id];
                // this_el = dS_dVa_r.coeff(row_id_dS_dVa_r, col_id_dS_dVa_r);
                BaseNRAlgo<LinearSolver>::value_map_.push_back(&BaseNRAlgo<LinearSolver>::dS_dVa_.coeffRef(row_id_dS_dVa_r, col_id_dS_dVa_r));

                // I don't need to perform these checks: if they failed, the element would not be in J_ in the first place
                // const int is_row_non_null = pq_inv[row_id_dS_dVa_r];
                // const int is_col_non_null = pq_inv[col_id_dS_dVa_r];
                // if(is_row_non_null >= 0 && is_col_non_null >= 0)
                //     this_el = dS_dVa_r.coeff(row_id_dS_dVa_r, col_id_dS_dVa_r);
                // else
                //     std::cout << "dS_dVa_r: missed" << std::endl;

            }else if((col_id < n_pvpq) && (row_id >= n_pvpq)){
                // this is the J21 part (dS_dVa_i)
                const int row_id_dS_dVa_i = pq[row_id - n_pvpq];
                const int col_id_dS_dVa_i = pvpq[col_id];
                // this_el = dS_dVa_i.coeff(row_id_dS_dVa_i, col_id_dS_dVa_i);
                BaseNRAlgo<LinearSolver>::value_map_.push_back(&BaseNRAlgo<LinearSolver>::dS_dVa_.coeffRef(row_id_dS_dVa_i, col_id_dS_dVa_i));
            }else if((col_id >= n_pvpq) && (row_id < n_pvpq)){
                // this is the J12 part (dS_dVm_r)
                const int row_id_dS_dVm_r = pvpq[row_id];
                const int col_id_dS_dVm_r = pq[col_id - n_pvpq];
                // this_el = dS_dVm_r.coeff(row_id_dS_dVm_r, col_id_dS_dVm_r);
                BaseNRAlgo<LinearSolver>::value_map_.push_back(&BaseNRAlgo<LinearSolver>::dS_dVm_.coeffRef(row_id_dS_dVm_r, col_id_dS_dVm_r));
            }else if((col_id >= n_pvpq) && (row_id >= n_pvpq)){
                // this is the J22 part (dS_dVm_i)
                const int row_id_dS_dVm_i = pq[row_id - n_pvpq];
                const int col_id_dS_dVm_i = pq[col_id - n_pvpq];
                // this_el = dS_dVm_i.coeff(row_id_dS_dVm_i, col_id_dS_dVm_i);
                BaseNRAlgo<LinearSolver>::value_map_.push_back(&BaseNRAlgo<LinearSolver>::dS_dVm_.coeffRef(row_id_dS_dVm_i, col_id_dS_dVm_i));
            }

            // go to the next element
            ++pos_el;
        }
    }
    // BaseNRAlgo<LinearSolver>::dS_dVa_.makeCompressed();
    // BaseNRAlgo<LinearSolver>::dS_dVm_.makeCompressed();
}

template<class LinearSolver>
void BaseNRSingleSlackAlgo<LinearSolver>::fill_jacobian_matrix_kown_sparsity_pattern(
        const Eigen::VectorXi & pq,
        const Eigen::VectorXi & pvpq
    )
{
    /**
    This functions fills the jacobian matrix when its sparsity pattern is KNOWN in advance (typically
    the second and next iterations of the Newton Raphson)
    We don't need to perform heavy computation but just to read the right data from the right matrix!
    This function is optimized for speed. In particular, it does not "uncompress" the J_ matrix and only
    change the value pointer (not the inner nor outer pointer)
    It is optimized only if J_ is in default Eigen format (column)
    Remember:
    J has the shape
    | J11 | J12 |               | (pvpq, pvpq) | (pvpq, pq) |
    | --------- | = dimensions: | ------------------------- |
    | J21 | J22 |               |  (pq, pvpq)  | (pq, pq) |
    python implementation:
    J11 = dS_dVa[array([pvpq]).T, pvpq].real
    J12 = dS_dVm[array([pvpq]).T, pq].real
    J21 = dS_dVa[array([pq]).T, pvpq].imag
    J22 = dS_dVm[array([pq]).T, pq].imag
    **/

    const int n_pvpq = static_cast<int>(pvpq.size());

    // real_type * J_x_ptr = J_.valuePtr();
    const int n_cols = static_cast<int>(BaseNRAlgo<LinearSolver>::J_.cols());  // equal to nrow
    unsigned int pos_el = 0;
    for (int col_id=0; col_id < n_cols; ++col_id){
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(BaseNRAlgo<LinearSolver>::J_, col_id); it; ++it)
        {
            const auto row_id = it.row();
            // only one if is necessary (magic !)
            // top rows are "real" part and bottom rows are imaginary part (you can check)
            it.valueRef() = row_id < n_pvpq ? std::real(*BaseNRAlgo<LinearSolver>::value_map_[pos_el]) : std::imag(*BaseNRAlgo<LinearSolver>::value_map_[pos_el]);
            // go to the next element
            ++pos_el;
        }
    }
}
