// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

template<class LinearSolver, class NRSystem>
bool NRAlgo<LinearSolver, NRSystem>::compute_pf(
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
    if (Sbus.size() != Ybus.rows() || Sbus.size() != Ybus.cols()) {
        std::ostringstream exc_;
        exc_ << "NRAlgo::compute_pf: Size of the Sbus should be the same as the size of Ybus. Currently: ";
        exc_ << "Sbus  (" << Sbus.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if (V.size() != Ybus.rows() || V.size() != Ybus.cols()) {
        std::ostringstream exc_;
        exc_ << "NRAlgo::compute_pf: Size of V (init voltages) should be the same as the size of Ybus. Currently: ";
        exc_ << "V  (" << V.size() << ") and Ybus (" << Ybus.rows() << ", " << Ybus.cols() << ").";
        throw std::runtime_error(exc_.str());
    }
    if (!is_linear_solver_valid()) return false;

    reset_timer();
    err_ = ErrorType::NoError;
    auto timer     = CustTimer();
    auto timer_pre = CustTimer();

    // Determine whether topology rebuild is required.
    const bool need_rebuild = (
        need_factorize_ ||
        _solver_control.need_reset_solver() ||
        _solver_control.has_dimension_changed() ||
        _solver_control.has_slack_participate_changed() ||
        _solver_control.ybus_change_sparsity_pattern() ||
        _solver_control.has_ybus_some_coeffs_zero() ||
        _solver_control.need_recompute_ybus() ||
        _solver_control.has_pv_changed() ||
        _solver_control.has_pq_changed()
    );

    if (need_rebuild) {
        // Reset linear solver state before re-initialization.
        ErrorType rs = _linear_solver.reset();
        if (rs != ErrorType::NoError) err_ = rs;
    }

    // linear solver failed to reset or I had a previous error
    // without a explicit call to "reset"
    if (!((err_ == ErrorType::NotInitError) || (err_ == ErrorType::NoError))) {
        timer_pre_proc_ += timer_pre.duration();
        timer_total_nr_ += timer.duration();
        return false;
    }

    // Phase 1: rebuild pvpq maps, lag, etc. (skipped when topology is unchanged).
    // std::cout << "need_factorize_ " << need_factorize_ << "\n";
    // std::cout << "_solver_control.need_reset_solver() " << _solver_control.need_reset_solver() << "\n";
    // std::cout << "_solver_control.has_dimension_changed() " << _solver_control.has_dimension_changed() << "\n";
    // std::cout << "_solver_control.has_slack_participate_changed() " << _solver_control.has_slack_participate_changed() << "\n";
    // std::cout << "_solver_control.ybus_change_sparsity_pattern() " << _solver_control.ybus_change_sparsity_pattern() << "\n";
    // std::cout << "_solver_control.has_ybus_some_coeffs_zero() " << _solver_control.has_ybus_some_coeffs_zero() << "\n";
    // std::cout << "_solver_control.need_recompute_ybus() " << _solver_control.need_recompute_ybus() << "\n";
    // std::cout << "_solver_control.has_pv_changed() " << _solver_control.has_pv_changed() << "\n";
    // std::cout << "_solver_control.has_pq_changed() " << _solver_control.has_pq_changed() << "\n";
    // std::cout << "need_rebuild " << need_rebuild << "\n";
    if (need_rebuild)
        _system.init_topology(Ybus, Sbus, slack_ids, slack_weights, pv, pq);

    // Phase 1.5: update V/Sbus pointers and initial voltage state (cheap, always).
    _system.update_state(Ybus, V, Sbus);

    // Phase 2: rebuild J sparsity + value_map, then initialize linear solver.
    bool need_init = need_rebuild;
    if (need_rebuild) {
        _system.build_J_sparsity();
        n_ = static_cast<int>(_system.J().cols());
    }

    timer_pre_proc_ += timer_pre.duration();

    // Initial mismatch (negated: Sbus - Scomp)
    RealVect F = _system.mismatch();
    bool converged = _check_for_convergence(F, tol);
    nr_iter_       = 0;
    bool res       = true;
    bool need_factorize = true;   // factorize on first iteration

    while ((!converged) & (nr_iter_ < max_iter)) {
        nr_iter_++;

        need_factorize = (need_factorize || should_refactor_policy(nr_iter_));
        if (need_factorize) {
            // Phase 3: fill J numerically with current V.
            _system.fill_J();

            if (need_init) {
                // New sparsity pattern: (re-)initialize factorization.
                auto timer_i = CustTimer();
                err_ = _linear_solver.initialize(_system.J());
                need_init = false;
                timer_initialize_ += timer_i.duration();
            } else {
                auto timer_r = CustTimer();
                err_ = _linear_solver.refactor(_system.J());
                timer_refactor_ += timer_r.duration();
            }
            if (err_ != ErrorType::NoError) { res = false; break; }
            need_factorize = false;
        }

        // Solve J * dx = F  (F = mismatch, negated convention; F overwritten with dx)
        {
            auto timer_s = CustTimer();
            err_ = _linear_solver.solve(F);
            timer_solve_ += timer_s.duration();
        }
        if (err_ != ErrorType::NoError) { res = false; break; }

        // Apply scaling policy (runtime dispatch)
        real_type coeff = scaling_policy_->scale(_system, F);
        auto timer_va_vm = CustTimer();
        _system.apply_step(coeff * F);
        timer_Va_Vm_ += timer_va_vm.duration();

        // New mismatch
        F = _system.mismatch();
        if (!F.allFinite()) { err_ = ErrorType::InifiniteValue; break; }
        converged = _check_for_convergence(F, tol);
        need_factorize = false;
    }

    if (!converged) {
        if (err_ == ErrorType::NoError) err_ = ErrorType::TooManyIterations;
        res = false;
    }
    timer_total_nr_ += timer.duration();

    // Synchronise BaseAlgo's voltage state
    V_  = _system.V();
    Vm_ = _system.Vm();
    Va_ = _system.Va();
    V   = V_;

    // Propagate NRSystem timers to NRAlgo
    timer_dSbus_ += _system.timer_dSbus();
    timer_fillJ_ += _system.timer_fillJ();
    _system.reset_timers();

    #ifdef __COUT_TIMES
        std::cout << "Computation time: "
                  << "\n\t timer_initialize_: " << timer_initialize_
                  << "\n\t timer_dSbus_: " << timer_dSbus_
                  << "\n\t timer_fillJ_: " << timer_fillJ_
                  << "\n\t timer_Fx_: " << timer_Fx_
                  << "\n\t timer_check_: " << timer_check_
                  << "\n\t timer_solve_: " << timer_solve_
                  << "\n\t timer_total_nr_: " << timer_total_nr_
                  << "\n\n";
    #endif
    return res;
}

template<class LinearSolver, class NRSystem>
void NRAlgo<LinearSolver, NRSystem>::reset()
{
    BaseAlgo::reset();
    _system.clear_jacobian();
    need_factorize_ = true;
    n_ = -1;
    ErrorType reset_status = _linear_solver.reset();
    if (reset_status != ErrorType::NoError) err_ = reset_status;
}
