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
    reset_if_needed();
    if (!((err_ == ErrorType::NotInitError) || (err_ == ErrorType::NoError))) return false;

    err_ = ErrorType::NoError;
    auto timer     = CustTimer();
    auto timer_pre = CustTimer();

    // Build pvpq maps, initialise voltage state
    _system.setup(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq);

    // Pre-loop: clear J caches if structural changes detected
    if (need_factorize_ ||
        _solver_control.need_reset_solver() ||
        _solver_control.has_dimension_changed() ||
        _solver_control.has_slack_participate_changed() ||
        _solver_control.ybus_change_sparsity_pattern() ||
        _solver_control.has_ybus_some_coeffs_zero() ||
        _solver_control.need_recompute_ybus() ||
        _solver_control.has_pv_changed() ||
        _solver_control.has_pq_changed())
    {
        _system.clear_jacobian();
    }

    timer_pre_proc_ += timer_pre.duration();

    // Initial mismatch (negated: Sbus - Scomp)
    RealVect F = _system.mismatch();
    bool converged = _check_for_convergence(F, tol);
    nr_iter_       = 0;
    bool res       = true;

    while ((!converged) & (nr_iter_ < max_iter)) {
        nr_iter_++;

        if (should_refactor(nr_iter_)) {
            // Build / update Jacobian
            _system.assemble_jacobian();

            // Factorize: initialize (new pattern) or refactor (same pattern)
            if (need_factorize_ || _system.pattern_changed()) {
                auto timer_i = CustTimer();
                n_    = static_cast<int>(_system.J().cols());
                err_  = _linear_solver.initialize(_system.J());
                need_factorize_ = false;
                _system.clear_pattern_changed();
                timer_initialize_ += timer_i.duration();
            } else {
                auto timer_r = CustTimer();
                err_ = _linear_solver.refactor(_system.J());
                timer_refactor_ += timer_r.duration();
            }
            if (err_ != ErrorType::NoError) { res = false; break; }
        }

        // Solve J * dx = F  (F = mismatch, negated convention; F overwritten with dx)
        {
            auto timer_s = CustTimer();
            err_ = _linear_solver.solve(F);
            timer_solve_ += timer_s.duration();
        }
        if (err_ != ErrorType::NoError) { res = false; break; }

        auto timer_va_vm = CustTimer();

        // Apply scaling policy (runtime dispatch)
        switch (scaling_policy_) {

        case ScalingPolicyType::NoScaling:
            _system.apply_step(F);
            break;

        case ScalingPolicyType::MaxVoltageChange: {
            real_type alpha = static_cast<real_type>(1.0);
            const NRLayout& lyt = _system.layout();
            if (lyt.theta_size() > 0) {
                const real_type max_abs_dtheta = lyt.theta(F).array().abs().maxCoeff();
                if (max_abs_dtheta > max_dVa_)
                    alpha = std::min(alpha, max_dVa_ / max_abs_dtheta);
            }
            if (lyt.vm_size() > 0) {
                const real_type max_abs_dvm = lyt.vm(F).array().abs().maxCoeff();
                if (max_abs_dvm > max_dVm_)
                    alpha = std::min(alpha, max_dVm_ / max_abs_dvm);
            }
            F *= alpha;
            _system.apply_step(F);
            break;
        }

        case ScalingPolicyType::LineSearch: {
            // Save current merit (||mismatch||^2 before step)
            const real_type F_norm_sq_0 = F.squaredNorm();
            real_type alpha = static_cast<real_type>(1.0);
            for (int k = 0; k < ls_max_iter_; ++k) {
                const real_type threshold =
                    (static_cast<real_type>(1.0)
                     - static_cast<real_type>(2.0) * ls_c_ * alpha) * F_norm_sq_0;
                if (_system.mismatch_sq_norm_at(alpha * F) <= threshold) break;
                alpha *= ls_rho_;
            }
            F *= alpha;
            _system.apply_step(F);
            break;
        }

        case ScalingPolicyType::Iwamoto: {
            const real_type g0 = F.squaredNorm();
            const real_type g1 = _system.mismatch_sq_norm_at(F);
            real_type mu = (g0 + g1 > static_cast<real_type>(0.))
                           ? g0 / (g0 + g1)
                           : static_cast<real_type>(1.0);
            mu = std::max(iw_mu_min_, std::min(iw_mu_max_, mu));
            F *= mu;
            _system.apply_step(F);
            break;
        }

        default:
            _system.apply_step(F);
            break;
        }

        timer_Va_Vm_ += timer_va_vm.duration();

        // New mismatch
        F = _system.mismatch();
        if (!F.allFinite()) { err_ = ErrorType::InifiniteValue; break; }
        converged = _check_for_convergence(F, tol);
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
