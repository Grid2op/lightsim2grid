// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// #include "BaseNRAlgo.h"  // now a template class, so this file will be included instead !

// TODO get rid of the pvpq, pv, pq etc and put the jacobian "in the right order"
// to ease and make way faster the filling of the sparse matrix J

template<class LinearSolver>
bool BaseNRAlgo<LinearSolver>::compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
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
        exc_ << "BaseNRAlgo::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(V.size() != Ybus.rows() || V.size() != Ybus.cols() ){
        // TODO DEBUG MODE
        std::ostringstream exc_;
        exc_ << "BaseNRAlgo::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows()<< ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(!is_linear_solver_valid()) {
        // err_ = ErrorType::NotInitError;
        return false;
    }
    reset_timer();
    reset_if_needed();
    err_ = ErrorType::NoError;  // reset the error if previous error happened
    auto timer = CustTimer();
    auto timer_pre_proc = CustTimer();

    Eigen::VectorXi my_pv = retrieve_pv_with_slack(slack_ids, pv);  // retrieve_pv_with_slack (not all), add_slack_to_pv (all)
    real_type slack_absorbed = std::real(Sbus.sum());  // initial guess for slack_absorbed
    const auto slack_bus_id = slack_ids(0);
    
    // initialize once and for all the "inverse" of these vectors
    const auto n_pv = my_pv.size();
    const auto n_pq = pq.size();
    Eigen::VectorXi pvpq(n_pv + n_pq);
    pvpq << my_pv, pq; 

    // some clever tricks are used in the making of the Jacobian to handle the slack bus 
    // (in case there is a distributed slack bus)
    const auto n_pvpq = pvpq.size();
    std::vector<int> pvpq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
    std::vector<int> pq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
    // const auto last_index = n_pvpq + n_pq;  // unused

    V_ = V;
    Vm_ = V_.array().abs();  // update Vm and Va again in case
    Va_ = V_.array().arg();  // we "wrapped around" with a negative Vm
    timer_pre_proc_ += timer_pre_proc.duration();

    // first check, if the problem is already solved, i stop there
    // compute a first time the mismatch to initialize the slack bus
    RealVect F = _evaluate_Fx(Ybus, V, Sbus, slack_bus_id, slack_absorbed, slack_weights, my_pv, pq);

    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    bool has_just_been_initialized = false;  // to avoid a call to klu_refactor follow a call to klu_factor in the same loop
    // std::cout << "iter " << nr_iter_ << " dx(0): " << -F(0) << " dx(1): " << -F(1) << std::endl;
    // std::cout << "slack_absorbed " << slack_absorbed << std::endl;
    if(need_factorize_ ||
       _solver_control.need_reset_solver() || 
       _solver_control.has_dimension_changed() ||
       _solver_control.has_slack_participate_changed() ||  // the full "ybus without slack" has changed, everything needs to be recomputed_solver_control.ybus_change_sparsity_pattern()
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero() ||
       _solver_control.need_recompute_ybus() ||
       _solver_control.has_slack_participate_changed() ||
       _solver_control.has_pv_changed() ||
       _solver_control.has_pq_changed()
       )
       {
        value_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::col_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::row_map_.clear();  // TODO smarter solver: only needed if ybus has changed
        dS_dVm_.resize(0,0);  // TODO smarter solver: only needed if ybus has changed
        dS_dVa_.resize(0,0);  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::dS_dVm_.setZero();  // TODO smarter solver: only needed if ybus has changed
        // BaseNRAlgo<LinearSolver>::dS_dVa_.setZero();  // TODO smarter solver: only needed if ybus has changed

       }
    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;
        fill_jacobian_matrix(Ybus, V_, slack_bus_id, slack_weights, pq, pvpq, pq_inv, pvpq_inv);

        if(need_factorize_){
            initialize();
            if(err_ != ErrorType::NoError){
                // I got an error during the initialization of the linear system, i need to stop here
                res = false;
                break;
            }
            has_just_been_initialized = true;
        }
        solve(F, has_just_been_initialized);

        has_just_been_initialized = false;
        if(err_ != ErrorType::NoError){
            // I got an error during the solving of the linear system, i need to stop here
            res = false;
            break;
        }
        // const auto dx = -F;  // removed for speed optimization (-= used below)

        auto timer_va_vm = CustTimer();
        slack_absorbed -= F(0); // by convention in fill_jacobian_matrix the slack bus is the first component
        // update voltage (this should be done consistently with "_evaluate_Fx")
        if (n_pv > 0) Va_(my_pv) -= F.segment(1, n_pv);
        if (n_pq > 0){
            Va_(pq) -= F.segment(n_pv + 1, n_pq);
            Vm_(pq) -= F.segment(n_pv + n_pq + 1, n_pq);
        }
        

        // std::cout << "iter " << nr_iter_ << " dx(0): " << -F(0) << " dx(1): " << -F(1) << std::endl;
        // std::cout << "slack_absorbed " << slack_absorbed << std::endl;
        // TODO change here for not having to cast all the time ... maybe
        V_ = Vm_.array() * (Va_.array().cos().template cast<cplx_type>() + my_i * Va_.array().sin().template cast<cplx_type>() );
        // V_ = Vm_.array() *  (my_i * Va_.array().template cast<cplx_type>()).exp() ;
        if(Vm_.minCoeff() < 0.)
        {
            // update Vm and Va again in case
            // we wrapped around with a negative Vm 
            // TODO more efficient way maybe ?
            Vm_ = V_.array().abs();  
            Va_ = V_.array().arg();  
            // TODO do I need to change V here ???
        }
        timer_Va_Vm_ += timer_va_vm.duration();

        F = _evaluate_Fx(Ybus, V_, Sbus, slack_bus_id, slack_absorbed, slack_weights, my_pv, pq);
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
    timer_total_nr_ += timer.duration();
    #ifdef __COUT_TIMES
        std::cout << "Computation time: " << "\n\t timer_initialize_: " << timer_initialize_
                  << "\n\t timer_dSbus_: " << timer_dSbus_
                  << "\n\t timer_fillJ_: " << timer_fillJ_
                  << "\n\t timer_Fx_: " << timer_Fx_
                  << "\n\t timer_check_: " << timer_check_
                  << "\n\t timer_solve_: " << timer_solve_
                  << "\n\t timer_total_nr_: " << timer_total_nr_
                  << "\n\n";
    #endif // __COUT_TIMES
    Vm_ = V_.array().abs();  
    Va_ = V_.array().arg();  
    _solver_control.tell_none_changed();
    return res;
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::reset(){
    BaseAlgo::reset();
    // reset specific attributes
    J_ = Eigen::SparseMatrix<real_type>();  // the jacobian matrix
    dS_dVm_ = Eigen::SparseMatrix<cplx_type>();
    dS_dVa_ = Eigen::SparseMatrix<cplx_type>();
    need_factorize_ = true;
    n_ = -1;
    // reset linear solver
    ErrorType reset_status = _linear_solver.reset();
    if(reset_status != ErrorType::NoError) err_ = reset_status;
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::_dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cplx_type> > & Ybus,
                                           const Eigen::Ref<const CplxVect > & V){
    // std::cout << "Ybus.nonZeros(): " << Ybus.nonZeros() << std::endl;
    auto timer = CustTimer();
    const auto size_dS = V.size();
    const CplxVect Vnorm = V.array() / V.array().abs();
    const CplxVect Ibus = Ybus * V;
    const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

    if (dS_dVm_.cols() != Ybus.cols())
    {
        // initiliaze the matrices once (especially for the sparsity pattern)
        // as I will use pointers to reference these elements later one, once initialize they need to be
        // fixed
        dS_dVm_ = Ybus;
        dS_dVa_ = Ybus;
    }

    cplx_type * ds_dvm_x_ptr = dS_dVm_.valuePtr();
    cplx_type * ds_dva_x_ptr = dS_dVa_.valuePtr();

    // there can be more value in dS_dVm_ or dS_dVa_ than in Ybus, if a line is disconnected for example
    // TODO smarter solver
    // const Eigen::Index * col_ptr = dS_dVm_.outerIndexPtr();
    // const Eigen::Index * row_ptr = dS_dVm_.innerIndexPtr();

    unsigned int pos_el = 0;
    for (int col_id=0; col_id < size_dS; ++col_id){
        for (Eigen::Ref<const Eigen::SparseMatrix<cplx_type> >::InnerIterator it(Ybus, col_id); it; ++it)
        {
            const int row_id = static_cast<int>(it.row());
            const cplx_type el_ybus = it.value();
            cplx_type & ds_dvm_el = ds_dvm_x_ptr[pos_el];
            cplx_type & ds_dva_el = ds_dva_x_ptr[pos_el];

            // assign the right value (the one in ybus)
            ds_dvm_el = el_ybus;
            ds_dva_el = el_ybus;

            // now compute the derivatives properly
            ds_dvm_el *= Vnorm(col_id);  // dS_dVm[k] *= Vnorm[Yj[k]]
            ds_dvm_el = std::conj(ds_dvm_el) * V(row_id);  // dS_dVm[k] = conj(dS_dVm[k]) * V[r]

            ds_dva_el *= V(col_id);  // dS_dVa[k] *= V[Yj[k]]

            if(col_id == row_id)
            {
                ds_dvm_el += conjIbus_Vnorm(row_id); // std::conj(Ibus(row_id)) * Vnorm(row_id); // dS_dVm[k] += conj(Ibus) * Vnorm
                ds_dva_x_ptr[pos_el] -= Ibus(row_id);  // dS_dVa[k] = dS_dVa[k] - Ibus[r]
            }
            cplx_type tmp = my_i * V(row_id);
            ds_dva_el = std::conj(-ds_dva_el) * tmp;  // dS_dVa[k] = conj(-dS_dVa[k]) * (1j * V[r])

            // go to next element
            ++pos_el;
        }
    }
    timer_dSbus_ += timer.duration();
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::_get_values_J(int & nb_obj_this_col,
                                 std::vector<Eigen::Index> & inner_index,
                                 std::vector<real_type> & values,
                                 const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & mat,  // ex. dS_dVa_r
                                 const std::vector<int> & index_row_inv, // ex. pvpq_inv
                                 const Eigen::VectorXi & index_col, // ex. pvpq
                                 Eigen::Index col_id,
                                 Eigen::Index row_lag,  // 0 for J11 for example, n_pvpq for J12
                                 Eigen::Index col_lag  // to remove the ref slack bus from this
                                 )
{
    /**
    This function will fill the "inner_index" and "values" with the non zero values
    present in the matrix "mat" for the column of the J matrix with id "col_id"
    which corresponds to the column "index_col(col_id)" of the matrix mat.
    The rows need to be converted using another vector too. For example, row "j" of J
    need to be filled with element k of matrix "mat" with k such that "index_row[j] = k"
    Hence, we pass as the argument of this function the "inverse" of index_row, which is such
    that : "j = index_row_inv[k]" is easily computable given k.
    **/
    int col_id_mat = index_col(col_id + col_lag);

    _get_values_J(nb_obj_this_col, inner_index, values, mat, index_row_inv, col_id_mat,
                  row_lag, col_lag);
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::_get_values_J(int & nb_obj_this_col,
                                 std::vector<Eigen::Index> & inner_index,
                                 std::vector<real_type> & values,
                                 const Eigen::Ref<const Eigen::SparseMatrix<real_type> > & mat,  // ex. dS_dVa_r
                                 const std::vector<int> & index_row_inv, // ex. pvpq_inv
                                 Eigen::Index col_id_mat, // ex. pvpq
                                 Eigen::Index row_lag,  // 0 for J11 for example, n_pvpq for J12
                                 Eigen::Index col_lag  // to remove the ref slack bus from this
                                 )
{
    /**
    This function will fill the "inner_index" and "values" with the non zero values
    present in the matrix "mat" for the column of the J matrix with id "col_id"
    which corresponds to the column col_id_mat  of the matrix mat.
    The rows need to be converted using another vector too. For example, row "j" of J
    need to be filled with element k of matrix "mat" with k such that "index_row[j] = k"
    Hence, we pass as the argument of this function the "inverse" of index_row, which is such
    that : "j = index_row_inv[k]" is easily computable given k.
    **/
    const Eigen::Index start_id = mat.outerIndexPtr()[col_id_mat];
    const Eigen::Index end_id = mat.outerIndexPtr()[col_id_mat+1];
    const real_type * val_prt = mat.valuePtr();
    for(Eigen::Index obj_id = start_id; obj_id < end_id; ++obj_id)
    {
        const Eigen::Index row_id_dS_dVa = mat.innerIndexPtr()[obj_id];
        // I add the value only if the rows was selected in the indexes
        const Eigen::Index row_id = index_row_inv[row_id_dS_dVa];
        if(row_id >= 0)
        {
            inner_index.push_back(row_id+row_lag);
            values.push_back(val_prt[obj_id]);
            nb_obj_this_col++;
        }
    }
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::fill_jacobian_matrix(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                                      const CplxVect & V,
                                                      Eigen::Index slack_bus_id,
                                                      const RealVect & slack_weights,
                                                      const Eigen::VectorXi & pq,
                                                      const Eigen::VectorXi & pvpq,
                                                      const std::vector<int> & pq_inv,
                                                      const std::vector<int> & pvpq_inv
                                                      )
{
    /**
    Remember, J has the shape:
    
    | s | slack_bus |               | (pvpq+1,1) |   (1, pvpq)  |  (1, pq)   |
    | l |  -------  |               |            | ------------------------- |
    | a | J11 | J12 | = dimensions: |            | (pvpq, pvpq) | (pvpq, pq) |
    | c | --------- |               |   ------   | ------------------------- |
    | k | J21 | J22 |               |  (pq, 1)   |  (pq, pvpq)  | (pq, pq)   |

    python implementation:
    `J11` = dS_dVa[array([pvpq]).T, pvpq].real
    `J12` = dS_dVm[array([pvpq]).T, pq].real
    `J21` = dS_dVa[array([pq]).T, pvpq].imag
    `J22` = dS_dVm[array([pq]).T, pq].imag

    `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real
    and dS_dVm[slack_bus_id, pq].real

    `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
    the remaining pq components are all 0.
    **/
    _dSbus_dV(Ybus, V);

    auto timer = CustTimer();
    const auto n_pvpq = pvpq.size();
    const auto n_pq = pq.size();
    const auto size_j = n_pvpq + n_pq + 1;   // +1 because i add the slack bus

    // TODO to gain a bit more time below, try to compute directly, in _dSbus_dV(Ybus, V);
    // TODO the `dS_dVa_[pvpq, pvpq]`
    // TODO so that it's easier to retrieve in the next few lines !
    if(J_.cols() != size_j)
    {
        #ifdef __COUT_TIMES
            auto timer2 = CustTimer();
        #endif  // __COUT_TIMES
        // first time i initialized the matrix, so i need to compute its sparsity pattern
        fill_jacobian_matrix_unkown_sparsity_pattern(Ybus, V, slack_bus_id, slack_weights, pq, pvpq, pq_inv, pvpq_inv);
        fill_value_map(slack_bus_id, pq, pvpq, false);
        #ifdef __COUT_TIMES
            std::cout << "\t\t fill_jacobian_matrix_unkown_sparsity_pattern : " << timer2.duration() << std::endl;
        #endif  // __COUT_TIMES
    }else{
        // the sparsity pattern of J_ is already known, i can reuse it to fill it
        // properly and faster (approx 3 times faster than the previous one)
        #ifdef __COUT_TIMES
            auto timer3 = CustTimer();
        #endif  // 
        if (BaseNRAlgo<LinearSolver>::value_map_.size() == 0) fill_value_map(slack_bus_id, pq, pvpq, true);
        fill_jacobian_matrix_kown_sparsity_pattern(slack_bus_id,
                                                   pq, pvpq
                                                   );
        #ifdef __COUT_TIMES
            std::cout << "\t\t fill_jacobian_matrix_kown_sparsity_pattern : " << timer3.duration() << std::endl;
        #endif  // __COUT_TIMES
    }
    timer_fillJ_ += timer.duration();;
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::fill_jacobian_matrix_unkown_sparsity_pattern(
        const Eigen::SparseMatrix<cplx_type> & Ybus,
        const CplxVect & V,
        Eigen::Index slack_bus_id,
        const RealVect & slack_weights,
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

    Remember, J has the shape:
    
    | s | slack_bus |               | (pvpq+1,1) |   (1, pvpq)  |  (1, pq)   |
    | l |  -------  |               |            | ------------------------- |
    | a | J11 | J12 | = dimensions: |            | (pvpq, pvpq) | (pvpq, pq) |
    | c | --------- |               |   ------   | ------------------------- |
    | k | J21 | J22 |               |  (pq, 1)   |  (pq, pvpq)  | (pq, pq)   |

    python implementation:
    `J11` = dS_dVa[array([pvpq]).T, pvpq].real
    `J12` = dS_dVm[array([pvpq]).T, pq].real
    `J21` = dS_dVa[array([pq]).T, pvpq].imag
    `J22` = dS_dVm[array([pq]).T, pq].imag

    `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real
    and dS_dVm[slack_bus_id, pq].real

    `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
    the remaining pq components are all 0.
    **/
   typedef Eigen::SparseMatrix<cplx_type>::StorageIndex StorageIndex;

    const Eigen::Index n_pvpq = pvpq.size();
    const Eigen::Index n_pq = pq.size();
    const auto size_j = n_pvpq + n_pq + 1;  // the +1 here to represent the equation for slack bus

    const Eigen::SparseMatrix<real_type> dS_dVa_r = dS_dVa_.real();
    const Eigen::SparseMatrix<real_type> dS_dVa_i = dS_dVa_.imag();
    const Eigen::SparseMatrix<real_type> dS_dVm_r = dS_dVm_.real();
    const Eigen::SparseMatrix<real_type> dS_dVm_i = dS_dVm_.imag();

    // Method (1) seems to be faster than the others

    // optim : if the matrix was already computed, i don't initialize it, i instead reuse as much as i can
    // i can do that because the matrix will ALWAYS have the same non zero coefficients.
    // in this if, i allocate it in a "large enough" place to avoid copy when first filling it
    if(J_.cols() != size_j) J_ = Eigen::SparseMatrix<real_type>(size_j, size_j);

    std::vector<Eigen::Triplet<double> > coeffs;  // HERE FOR PERF OPTIM (3)
    coeffs.reserve(2*(dS_dVa_.nonZeros()+dS_dVm_.nonZeros())  + slack_weights.size());  // HERE FOR PERF OPTIM (3)

    // i fill the buffer columns per columns
    int nb_obj_this_col = 0;
    std::vector<Eigen::Index> inner_index;
    std::vector<real_type> values;

    // fill n_pvpq leftmost columns
    for(Eigen::Index col_id=0; col_id < n_pvpq; ++col_id){ 
        // reset from the previous column
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
        // and check the row order !
        _get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVa_r,
                      pvpq_inv, pvpq,
                      col_id,
                      1,  // added because the first row is for the slack bus
                      0);
        // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
        _get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVa_i,
                      pq_inv, pvpq,
                      col_id,
                      n_pvpq + 1, // + 1 added because the first row is for the slack bus
                      0);

        // "efficient" insert of the element in the matrix
        for(Eigen::Index in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            StorageIndex row_id = static_cast<StorageIndex>(inner_index[in_ind]);
            coeffs.push_back(Eigen::Triplet<double>(row_id, static_cast<StorageIndex>(col_id) + 1, values[in_ind]));   // HERE FOR PERF OPTIM (3)
        }
    }

    //TODO make same for the second part (have a funciton for previous loop)
    // fill the remaining n_pq columns
    for(Eigen::Index col_id=0; col_id < n_pq; ++col_id){
        // reset from the previous column
        nb_obj_this_col = 0;
        inner_index.clear();
        values.clear();

        // fill with the first column with the column of dS_dVa[:,pvpq[col_id]]
        // and check the row order !
        _get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVm_r,
                      pvpq_inv, pq,
                      col_id,
                      1,  // 1 added because the first row is for the slack bus
                      0);

        // fill the rest of the rows with the first column of dS_dVa_imag[:,pq[col_id]]
        _get_values_J(nb_obj_this_col, inner_index, values,
                      dS_dVm_i,
                      pq_inv, pq,
                      col_id,
                      n_pvpq + 1,   // + 1 added because the first row is for the slack bus
                      0);

        // "efficient" insert of the element in the matrix
        for(Eigen::Index in_ind=0; in_ind < nb_obj_this_col; ++in_ind){
            auto row_id = static_cast<StorageIndex>(inner_index[in_ind]);
            coeffs.push_back(Eigen::Triplet<double>(row_id, static_cast<StorageIndex>(col_id + n_pvpq) + 1, values[in_ind]));   // HERE FOR PERF OPTIM (3)
        }
    }

    // add the slack bus coefficients for the first bus (the "real" slack bus)
    // first row (which corresponds to slack_bus_id)
    const int nb_bus = dS_dVa_r.cols();
    for (Eigen::Index col_id=0; col_id < nb_bus; ++col_id){
        const auto J_col = pvpq_inv[col_id];
        if(J_col < 0) continue;
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVa_r, col_id); it; ++it)
        {
            if(it.row() != slack_bus_id) continue;   // don't add it if it's not the ref slack bus          
            coeffs.push_back(Eigen::Triplet<double>(0, J_col + 1, it.value()));   // HERE FOR PERF OPTIM (3)
        }
    }
    for (Eigen::Index col_id=0; col_id < nb_bus; ++col_id){
        const auto J_col = pq_inv[col_id];
        if(J_col < 0) continue;
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(dS_dVm_r, col_id); it; ++it)
        {
            if(it.row() != slack_bus_id) continue;   // don't add it if it's not the ref slack bus
            coeffs.push_back(Eigen::Triplet<double>(0, static_cast<StorageIndex>(J_col + n_pvpq) + 1, it.value()));   // HERE FOR PERF OPTIM (3)
        }
    }

    // add later on the last column which corresponds to the slack bus equation
    const StorageIndex last_col = 0; // static_cast<StorageIndex>(size_j) - 1;
    // add the ref slack bus coeff
    coeffs.push_back(Eigen::Triplet<double>(0, last_col, slack_weights[slack_bus_id]));   // HERE FOR PERF OPTIM (3)
    // add the other coeffs (for other buses)
    auto row_j = 1;
    for(auto ind: pvpq){
        auto sl_w  = slack_weights(ind);
        if(sl_w != 0.) coeffs.push_back(Eigen::Triplet<double>(row_j, last_col, sl_w));   // HERE FOR PERF OPTIM (3)
        ++row_j;
    }
    J_.setFromTriplets(coeffs.begin(), coeffs.end());  // HERE FOR PERF OPTIM (3)
    // std::cout << "end fill jacobian unknown " << std::endl;
    J_.makeCompressed();
    // std::cout << "end fill_value_map" << std::endl;
}

/**
fill the value of the `value_map_` that stores pointers to the elements of
dS_dVa_ and dS_dVm_ to be used to fill J_
it requires that J_ is initialized, in compressed mode.
**/
template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::fill_value_map(
        Eigen::Index slack_bus_id,
        const Eigen::VectorXi & pq,
        const Eigen::VectorXi & pvpq,
        bool reset_J
        )
{
    const int n_pvpq = static_cast<int>(pvpq.size());
    value_map_ = std::vector<cplx_type*> ();
    value_map_.reserve(BaseNRAlgo<LinearSolver>::J_.nonZeros());
    // col_map_ = std::vector<int> (J_.nonZeros());
    // row_map_ = std::vector<int> (J_.nonZeros());

    const auto n_row = J_.cols();
    unsigned int pos_el = 0;
    for (Eigen::Index col_=1; col_ < n_row; ++col_){  // last column is never updated (slack equation)
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_); it; ++it)
        {
            auto row_id = it.row();
            const auto col_id = it.col() - 1;  // it's equal to "col_"
            if(reset_J) it.valueRef() = 0.; // "forget" previous J value in this setting

            if(row_id==0){
                // this is the row of the slack bus
                const Eigen::Index row_id_dS_dVx_r = slack_bus_id;  // same for both matrices
                if(col_id < n_pvpq){
                    const int col_id_dS_dVa_r = pvpq[col_id];
                    value_map_.push_back(&dS_dVa_.coeffRef(row_id_dS_dVx_r, col_id_dS_dVa_r));
                }
                else{
                    const int col_id_dS_dVm_r = pq[col_id - n_pvpq];
                    value_map_.push_back(&dS_dVm_.coeffRef(row_id_dS_dVx_r, col_id_dS_dVm_r));
                }
            }else{
                row_id -= 1;  // "do not consider" the row for slack bus (handled above)
                // this is consistent with the "+1" added in fill_jacobian_matrix_unkown_sparsity_pattern
                if((col_id < n_pvpq) && (row_id < n_pvpq)){
                    // this is the J11 part (dS_dVa_r)
                    const int row_id_dS_dVa_r = pvpq[row_id];
                    const int col_id_dS_dVa_r = pvpq[col_id];
                    // this_el = dS_dVa_r.coeff(row_id_dS_dVa_r, col_id_dS_dVa_r);
                    value_map_.push_back(&dS_dVa_.coeffRef(row_id_dS_dVa_r, col_id_dS_dVa_r));

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
                    value_map_.push_back(&dS_dVa_.coeffRef(row_id_dS_dVa_i, col_id_dS_dVa_i));
                }else if((col_id >= n_pvpq) && (row_id < n_pvpq)){
                    // this is the J12 part (dS_dVm_r)
                    const int row_id_dS_dVm_r = pvpq[row_id];
                    const int col_id_dS_dVm_r = pq[col_id - n_pvpq];
                    // this_el = dS_dVm_r.coeff(row_id_dS_dVm_r, col_id_dS_dVm_r);
                    value_map_.push_back(&dS_dVm_.coeffRef(row_id_dS_dVm_r, col_id_dS_dVm_r));
                }else if((col_id >= n_pvpq) && (row_id >= n_pvpq)){
                    // this is the J22 part (dS_dVm_i)
                    const int row_id_dS_dVm_i = pq[row_id - n_pvpq];
                    const int col_id_dS_dVm_i = pq[col_id - n_pvpq];
                    // this_el = dS_dVm_i.coeff(row_id_dS_dVm_i, col_id_dS_dVm_i);
                    value_map_.push_back(&dS_dVm_.coeffRef(row_id_dS_dVm_i, col_id_dS_dVm_i));
                }
            }
            // go to the next element
            ++pos_el;
        }
    }
    dS_dVa_.makeCompressed();
    dS_dVm_.makeCompressed();
}

template<class LinearSolver>
void BaseNRAlgo<LinearSolver>::fill_jacobian_matrix_kown_sparsity_pattern(
        Eigen::Index slack_bus_id,
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

    Remember, J has the shape:
    
    | s | slack_bus |               | (pvpq+1,1) |   (1, pvpq)  |  (1, pq)   |
    | l |  -------  |               |            | ------------------------- |
    | a | J11 | J12 | = dimensions: |            | (pvpq, pvpq) | (pvpq, pq) |
    | c | --------- |               |   ------   | ------------------------- |
    | k | J21 | J22 |               |  (pq, 1)   |  (pq, pvpq)  | (pq, pq)   |

    python implementation:
    `J11` = dS_dVa[array([pvpq]).T, pvpq].real
    `J12` = dS_dVm[array([pvpq]).T, pq].real
    `J21` = dS_dVa[array([pq]).T, pvpq].imag
    `J22` = dS_dVm[array([pq]).T, pq].imag

    `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real
    and dS_dVm[slack_bus_id, pq].real

    `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
    the remaining pq components are all 0.

    **/

    const auto n_pvpq_1 = pvpq.size() + 1;
    const auto n_cols = J_.cols();  // equal to nrow
    unsigned int pos_el = 0;
    for (Eigen::Index col_id=1; col_id < n_cols; ++col_id){  // last column is not updated (slack equation)
        for (Eigen::SparseMatrix<real_type>::InnerIterator it(J_, col_id); it; ++it)
        {
            const auto row_id = it.row();
            // only one if is necessary (magic !)
            // top rows are "real" part and bottom rows are imaginary part (you can check)
            it.valueRef() = row_id < n_pvpq_1 ? std::real(*value_map_[pos_el]) : std::imag(*value_map_[pos_el]);
            // go to the next element
            ++pos_el;
        }
    }
}
