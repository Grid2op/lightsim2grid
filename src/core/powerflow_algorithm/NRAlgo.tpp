// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

template<class LinearSolver, class SlackPolicy, class ScalingPolicy>
bool NRAlgo<LinearSolver, SlackPolicy, ScalingPolicy>::compute_pf(
        const Eigen::SparseMatrix<cplx_type>& Ybus,
        CplxVect& V,
        const CplxVect& Sbus,
        Eigen::Ref<const IntVect> slack_ids,
        const RealVect& slack_weights,
        Eigen::Ref<const IntVect> pv,
        Eigen::Ref<const IntVect> pq,
        int max_iter,
        real_type tol)
{
    if(Sbus.size() != Ybus.rows() || Sbus.size() != Ybus.cols()){
        std::ostringstream exc_;
        exc_ << "NRAlgo::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(V.size() != Ybus.rows() || V.size() != Ybus.cols()){
        std::ostringstream exc_;
        exc_ << "NRAlgo::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if(!is_linear_solver_valid()){
        return false;
    }
    reset_timer();
    reset_if_needed();
    if(!((err_ == ErrorType::NotInitError) || (err_ == ErrorType::NoError))){
        // in case the reset fails
        return false;
    }

    err_ = ErrorType::NoError;
    auto timer = CustTimer();
    auto timer_pre_proc = CustTimer();

    // policy-dependent setup: pv buses and lag
    Eigen::VectorXi my_pv = SlackPolicy::get_my_pv(*this, slack_ids, pv);
    real_type slack_absorbed = SlackPolicy::initial_slack_absorbed(Sbus);
    const auto slack_bus_id = slack_ids(0);

    const int n_pv = static_cast<int>(my_pv.size());
    const int n_pq = static_cast<int>(pq.size());
    Eigen::VectorXi pvpq(n_pv + n_pq);
    pvpq << my_pv, pq;
    const int n_pvpq = static_cast<int>(pvpq.size());
    std::vector<int> pvpq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
    std::vector<int> pq_inv(V.size(), -1);
    for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;

    _layout = NRLayout(n_pv, n_pq, SlackPolicy::lag);

    V_ = V;
    Vm_ = V_.array().abs();
    Va_ = V_.array().arg();
    timer_pre_proc_ += timer_pre_proc.duration();

    // initial residual (also initialises slack_absorbed for multi-slack)
    RealVect F = SlackPolicy::evaluate_Fx(*this, Ybus, V, Sbus,
                                           slack_bus_id, slack_absorbed, slack_weights,
                                           my_pv, pq);
    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0;
    bool res = true;

    bool do_refactorize = false;
    if(need_factorize_ ||
       _solver_control.need_reset_solver() ||
       _solver_control.has_dimension_changed() ||
       _solver_control.has_slack_participate_changed() ||
       _solver_control.ybus_change_sparsity_pattern() ||
       _solver_control.has_ybus_some_coeffs_zero() ||
       _solver_control.need_recompute_ybus() ||
       _solver_control.has_pv_changed() ||
       _solver_control.has_pq_changed())
    {
        value_map_.clear();
        dS_dVm_.resize(0, 0);
        dS_dVa_.resize(0, 0);
        do_refactorize = true;
    }

    while((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;

        // corresponds to actual policy (always refactor)
        // but should be updated with the integration of RefactorPolicies
        do_refactorize = true;

        if(do_refactorize){
            SlackPolicy::fill_jacobian_matrix(*this, Ybus, V_, slack_bus_id, slack_weights,
                                              pq, pvpq, pq_inv, pvpq_inv);
            if(need_factorize_){
                auto timer_i = CustTimer();
                n_ = static_cast<int>(J_.cols());
                err_ = _linear_solver.initialize(J_);
                need_factorize_ = false;
                timer_initialize_ += timer_i.duration();
            } else {
                auto timer_r = CustTimer();
                err_ = _linear_solver.refactor(J_);
                timer_refactor_ += timer_r.duration();
            }
            if(err_ != ErrorType::NoError){ res = false; break; }
        }

        {
            auto timer_s = CustTimer();
            err_ = _linear_solver.solve(F);
            timer_solve_ += timer_s.duration();
        }
        if(err_ != ErrorType::NoError){ res = false; break; }

        auto timer_va_vm = CustTimer();

        // scale the Newton step (alpha=1.0 for NoScalingPolicy, <1.0 for limiting policies)
        const real_type alpha = _scaling_policy_.compute_scale(*this, F);
        F *= alpha;

        // policy-dependent slack update (no-op for single-slack)
        SlackPolicy::update_slack_absorbed(*this, F, slack_absorbed);

        // update voltage angles and magnitudes
        if(_layout.nb_pv() > 0) Va_(my_pv) -= _layout.theta(F).segment(0, _layout.nb_pv());
        if(_layout.nb_pq() > 0){
            Va_(pq) -= _layout.theta(F).segment(_layout.nb_pv(), _layout.nb_pq());
            Vm_(pq) -= _layout.vm(F);
        }

        const cplx_type m_i = BaseConstants::my_i;
        V_ = Vm_.array() * (Va_.array().cos().template cast<cplx_type>()
                           + m_i * Va_.array().sin().template cast<cplx_type>());
        if(Vm_.minCoeff() < 0.){
            // wrapped around with a negative Vm
            Vm_ = V_.array().abs();
            Va_ = V_.array().arg();
        }
        timer_Va_Vm_ += timer_va_vm.duration();

        F = SlackPolicy::evaluate_Fx(*this, Ybus, V_, Sbus,
                                      slack_bus_id, slack_absorbed, slack_weights,
                                      my_pv, pq);
        bool tmp = F.allFinite();
        if(!tmp){
            err_ = ErrorType::InifiniteValue;
            break;
        }
        converged = _check_for_convergence(F, tol);
    }

    if(!converged){
        if(err_ == ErrorType::NoError) err_ = ErrorType::TooManyIterations;
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
    return res;
}

template<class LinearSolver, class SlackPolicy, class ScalingPolicy>
void NRAlgo<LinearSolver, SlackPolicy, ScalingPolicy>::reset()
{
    BaseAlgo::reset();
    J_ = Eigen::SparseMatrix<real_type>();
    dS_dVm_ = Eigen::SparseMatrix<cplx_type>();
    dS_dVa_ = Eigen::SparseMatrix<cplx_type>();
    need_factorize_ = true;
    n_ = -1;
    ErrorType reset_status = _linear_solver.reset();
    if(reset_status != ErrorType::NoError) err_ = reset_status;
}

template<class LinearSolver, class SlackPolicy, class ScalingPolicy>
void NRAlgo<LinearSolver, SlackPolicy, ScalingPolicy>::_dSbus_dV(
        const Eigen::Ref<const Eigen::SparseMatrix<cplx_type>>& Ybus,
        const Eigen::Ref<const CplxVect>& V)
{
    auto timer = CustTimer();
    const auto size_dS = V.size();
    const CplxVect Vnorm = V.array() / V.array().abs();
    const CplxVect Ibus = Ybus * V;
    const CplxVect conjIbus_Vnorm = Ibus.array().conjugate() * Vnorm.array();

    if(dS_dVm_.cols() != Ybus.cols()){
        dS_dVm_ = Ybus;
        dS_dVa_ = Ybus;
    }

    cplx_type* ds_dvm_x_ptr = dS_dVm_.valuePtr();
    cplx_type* ds_dva_x_ptr = dS_dVa_.valuePtr();

    unsigned int pos_el = 0;
    for(int col_id=0; col_id < size_dS; ++col_id){
        for(Eigen::Ref<const Eigen::SparseMatrix<cplx_type>>::InnerIterator it(Ybus, col_id); it; ++it)
        {
            const int row_id = static_cast<int>(it.row());
            const cplx_type el_ybus = it.value();
            cplx_type& ds_dvm_el = ds_dvm_x_ptr[pos_el];
            cplx_type& ds_dva_el = ds_dva_x_ptr[pos_el];

            ds_dvm_el = el_ybus;
            ds_dva_el = el_ybus;

            ds_dvm_el *= Vnorm(col_id);
            ds_dvm_el = std::conj(ds_dvm_el) * V(row_id);

            ds_dva_el *= V(col_id);

            if(col_id == row_id){
                ds_dvm_el += conjIbus_Vnorm(row_id);
                ds_dva_x_ptr[pos_el] -= Ibus(row_id);
            }
            cplx_type tmp = BaseConstants::my_i * V(row_id);
            ds_dva_el = std::conj(-ds_dva_el) * tmp;

            ++pos_el;
        }
    }
    timer_dSbus_ += timer.duration();
}

template<class LinearSolver, class SlackPolicy, class ScalingPolicy>
void NRAlgo<LinearSolver, SlackPolicy, ScalingPolicy>::_get_values_J(
        int& nb_obj_this_col,
        std::vector<Eigen::Index>& inner_index,
        std::vector<real_type>& values,
        const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
        const std::vector<int>& index_row_inv,
        const Eigen::VectorXi& index_col,
        size_t col_id,
        size_t row_lag,
        size_t col_lag)
{
    int col_id_mat = index_col(col_id + col_lag);
    _get_values_J(nb_obj_this_col, inner_index, values, mat, index_row_inv,
                  col_id_mat, row_lag, col_lag);
}

template<class LinearSolver, class SlackPolicy, class ScalingPolicy>
void NRAlgo<LinearSolver, SlackPolicy, ScalingPolicy>::_get_values_J(
        int& nb_obj_this_col,
        std::vector<Eigen::Index>& inner_index,
        std::vector<real_type>& values,
        const Eigen::Ref<const Eigen::SparseMatrix<real_type>>& mat,
        const std::vector<int>& index_row_inv,
        size_t col_id_mat,
        size_t row_lag,
        size_t col_lag)
{
    const int start_id = mat.outerIndexPtr()[col_id_mat];
    const int end_id = mat.outerIndexPtr()[col_id_mat + 1];
    const real_type* val_prt = mat.valuePtr();
    for(size_t obj_id = start_id; obj_id < (size_t)end_id; ++obj_id)
    {
        const int row_id_dS = mat.innerIndexPtr()[obj_id];
        const int row_id = index_row_inv[row_id_dS];
        if(row_id >= 0){
            inner_index.push_back(row_id + row_lag);
            values.push_back(val_prt[obj_id]);
            nb_obj_this_col++;
        }
    }
}
