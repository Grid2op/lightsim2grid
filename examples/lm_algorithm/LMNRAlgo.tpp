// Copyright (c) 2026
// Experimental plugin, not part of the lightsim2grid core distribution.

template<class LinearSolver, class NRSystemT>
bool LMNRAlgo<LinearSolver, NRSystemT>::compute_pf(
        const EigenRefConstCplxSpMat     & Ybus,
        const Eigen::Ref<const CplxVect> & V,
        const Eigen::Ref<const CplxVect> & Sbus,
        const Eigen::Ref<const IntVect>  & slack_ids,
        const Eigen::Ref<const RealVect> & slack_weights,
        const Eigen::Ref<const IntVect>  & pv,
        const Eigen::Ref<const IntVect>  & pq,
        int                              max_iter,
        real_type                        tol)
{
    if (!is_linear_solver_valid()) return false;

    reset_timer();
    err_ = ErrorType::NoError;
    auto timer = CustTimer();

    // Determine whether topology rebuild is required (same logic as NRAlgo).
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
        ErrorType rs = _linear_solver.reset();
        if (rs != ErrorType::NoError) err_ = rs;
    }
    if (!((err_ == ErrorType::NotInitError) || (err_ == ErrorType::NoError))) {
        timer_total_nr_ += timer.duration();
        return false;
    }

    _system.update_state(BaseAlgo::lsgrid_ptr_, Ybus, V, Sbus, slack_weights, _solver_control);
    if (need_rebuild) _system.init_topology(slack_ids, slack_weights, pv, pq);

    bool need_init = need_rebuild;
    if (need_rebuild) {
        _system.build_J_sparsity();
        n_ = static_cast<int>(_system.J().cols());
    }

    // Same contract as NRAlgo::compute_pf: update_state runs every solve but
    // init_topology / register_in only on a rebuild, so an extension whose
    // element count changed in between would index its stale row / column
    // caches out of bounds. Cheap, and must come before the first mismatch().
    _system.validate_registration();

    RealVect F = _system.mismatch();
    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0;
    bool res = true;
    lambda_ = lambda_init_;

    while ((!converged) && (nr_iter_ < max_iter)) {
        ++nr_iter_;

        // Baseline (undamped) Jacobian for this outer iteration, once -- LM
        // always rebuilds it fresh (lambda_ changes every outer iteration,
        // so there is no stable-across-iterations J to reuse the way
        // NRAlgo's RefactorPolicy does).
        _system.fill_internal_variables();
        _system.fill_J();

        if (need_init) {
            err_ = _linear_solver.analyze(_system.J());
            if (err_ == ErrorType::NoError) err_ = _linear_solver.factorize(_system.J());
            need_init = false;
            need_factorize_ = false;
        } else {
            err_ = _linear_solver.refactorize(_system.J());
        }
        if (err_ != ErrorType::NoError) { res = false; break; }

        const real_type F_sq = F.squaredNorm();
        real_type lambda_applied = static_cast<real_type>(0.);  // matches the just-built undamped baseline
        bool accepted = false;
        int trial = 0;
        RealVect dx;

        while (!accepted && trial < max_trials_) {
            ++trial;

            const real_type delta = lambda_ - lambda_applied;
            if (delta != static_cast<real_type>(0.)) {
                // Fast path: only the n diagonal positions move, not the
                // whole Jacobian -- see NRSystem::update_trailing_feature_values.
                _system.update_trailing_feature_values(n_, RealVect::Constant(n_, delta));
                lambda_applied = lambda_;
                err_ = _linear_solver.refactorize(_system.J());
                if (err_ != ErrorType::NoError) break;
            }

            dx = F;  // solve overwrites in place (mismatch, negated convention)
            err_ = _linear_solver.solve(dx);
            if (err_ != ErrorType::NoError) break;

            // Evaluate the TRUE nonlinear residual at the trial step without
            // mutating _system's V/Va/Vm -- only commit (apply_step) once accepted.
            const real_type trial_sq = _system.mismatch_sq_norm_at(dx);

            if (trial_sq < F_sq) {
                _system.apply_step(dx);
                accepted = true;
                lambda_ = std::max(lambda_ / nu_, lambda_min_);
            } else {
                lambda_ = std::min(lambda_ * nu_ + lambda_eps_, lambda_max_);
            }
        }
        if (err_ != ErrorType::NoError) { res = false; break; }
        if (!accepted) {
            // Exhausted max_trials_ without finding a step that reduces the
            // mismatch, even at lambda_max_: treat like NRAlgo's "linear
            // solver failed" bailout.
            err_ = ErrorType::SingularMatrix;
            res = false;
            break;
        }

        F = _system.mismatch();
        if (!F.allFinite()) { err_ = ErrorType::InifiniteValue; res = false; break; }
        converged = _check_for_convergence(F, tol);
    }

    if (!converged && res) {
        if (err_ == ErrorType::NoError) err_ = ErrorType::TooManyIterations;
        res = false;
    }

    V_  = _system.V();
    Vm_ = _system.Vm();
    Va_ = _system.Va();

    timer_dSbus_ += _system.timer_dSbus();
    timer_fillJ_ += _system.timer_fillJ();
    _system.reset_timers();
    timer_total_nr_ += timer.duration();

    return res;
}
