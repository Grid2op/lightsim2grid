// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

template<class Inner>
bool OuterLoopAlgo<Inner>::compute_pf(
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
    // what changed since the last solve is told to this object (LSGrid, the batches)
    inner_.algo().tell_solver_control(_solver_control);
    if (!inner_.begin()) {
        err_ = inner_.algo().get_error();
        return false;
    }

    reset_timer();
    err_ = ErrorType::NoError;
    auto timer = CustTimer();
    stats_.clear();
    total_nr_iterations_ = 0;
    unrealistic_ = false;

    // the loops edit these, never the grid's
    Sbus_init_ = Sbus;
    Sbus_target_ = Sbus;
    Sbus_ = Sbus;

    // ... and Ybus, whose values the phase shifters patch (see BranchControl)
    Ybus_ = Ybus;

    bool need_init = false;
    if (!inner_.setup(Ybus_, V, Sbus_, slack_ids, slack_weights, pv, pq, need_init,
                      [this](){ return _declare(); })) {
        err_ = inner_.algo().get_error();
        timer_total_nr_ += timer.duration();
        return false;
    }

    slack_bus_ = slack_ids.size() > 0 ? slack_ids(0) : -1;
    state_.Sbus = &Sbus_;
    state_.Sbus_init = &Sbus_init_;
    state_.Sbus_target = &Sbus_target_;
    state_.gen_target_p.clear();
    state_.storage_target_p.clear();
    state_.shunt_control.clear();
    state_.shunt_sections.clear();
    inner_.controls().reset_states();
    OuterState & state = state_;

    // OpenLoadFlow's isNeeded filter, then initialize, both before the first solve
    std::vector<BaseOuterLoop *> active;
    {
        OuterContext ctx = _context(&state);
        ctx.V = nullptr;  // nothing solved yet
        ctx.Va = nullptr;
        ctx.Vm = nullptr;
        for (auto & loop : loops_) {
            if (loop->is_needed(ctx)) active.push_back(loop.get());
        }
        for (BaseOuterLoop * loop : active) loop->initialize(ctx);
    }
    for (BaseOuterLoop * loop : active) stats_.loop_iterations.emplace_back(loop->name(), 0);

    // robust mode: the unrealistic-voltage check waits for the last loop able to fix it
    int last_fixing = -1;
    if (params_.voltage_remote_control_robust_mode) {
        for (std::size_t i = 0; i < active.size(); ++i) {
            if (active[i]->can_fix_unrealistic_state()) last_fixing = static_cast<int>(i);
        }
    }
    bool check_unrealistic = last_fixing < 0;

    bool solver_ok = _solve(max_iter, tol, need_init, check_unrealistic);

    OuterLoopStatus last_status = OuterLoopStatus::STABLE;
    std::string last_name;
    const BaseOuterLoop * last_unstable = nullptr;
    if (solver_ok) {
        int old_total;
        do {
            check_unrealistic = last_fixing < 0;
            old_total = total_nr_iterations_;
            ++stats_.nb_passes;
            for (std::size_t i = 0; i < active.size(); ++i) {
                BaseOuterLoop * loop = active[i];
                if (loop == last_unstable || !solver_ok || last_status == OuterLoopStatus::FAILED ||
                    stats_.nb_outer_iterations >= params_.max_outer_iterations) {
                    break;
                }
                // runOuterLoop: re-solve until this loop is stable
                OuterLoopStatus status;
                do {
                    OuterContext ctx = _context(&state);
                    ctx.iteration = stats_.loop_iterations[i].second;
                    const std::size_t first_decision = stats_.decisions.size();
                    status = loop->check(ctx);
                    for (std::size_t d = first_decision; d < stats_.decisions.size(); ++d) {
                        stats_.decisions[d].loop = loop->name();
                        stats_.decisions[d].outer_iteration = stats_.nb_outer_iterations;
                    }
                    last_status = status;
                    last_name = loop->name();
                    if (status == OuterLoopStatus::UNSTABLE) {
                        last_unstable = loop;
                        solver_ok = _solve(max_iter, tol, need_init, check_unrealistic);
                        ++stats_.nb_outer_iterations;
                        ++stats_.loop_iterations[i].second;
                    }
                } while (status == OuterLoopStatus::UNSTABLE && solver_ok &&
                         stats_.nb_outer_iterations < params_.max_outer_iterations);

                if (!check_unrealistic && static_cast<int>(i) == last_fixing && solver_ok &&
                    lsgrid_ptr_ != nullptr) {
                    // the deferred check, once the last loop able to fix it is done
                    if (is_state_unrealistic(*lsgrid_ptr_, V_, inner_.vm_unknown(), params_)) {
                        unrealistic_ = true;
                        solver_ok = false;
                    }
                }
                if (static_cast<int>(i) == last_fixing) check_unrealistic = true;
            }
        } while (total_nr_iterations_ > old_total && solver_ok &&
                 last_status != OuterLoopStatus::FAILED &&
                 stats_.nb_outer_iterations < params_.max_outer_iterations);
    }

    {
        OuterContext ctx = _context(&state);
        for (auto it = active.rbegin(); it != active.rend(); ++it) (*it)->cleanup(ctx);
    }

    if (last_status == OuterLoopStatus::FAILED) {
        stats_.status = OuterLoopStatus::FAILED;
        stats_.failed_loop = last_name;
    } else {
        stats_.status = stats_.nb_outer_iterations < params_.max_outer_iterations
                        ? OuterLoopStatus::STABLE : OuterLoopStatus::UNSTABLE;
    }
    stats_.unrealistic_state = unrealistic_;

    bool res = solver_ok && stats_.status == OuterLoopStatus::STABLE;
    if (solver_ok) {
        if (stats_.status == OuterLoopStatus::FAILED) err_ = ErrorType::OuterLoopFailed;
        else if (stats_.status == OuterLoopStatus::UNSTABLE) err_ = ErrorType::TooManyIterations;
    } else if (unrealistic_) {
        err_ = ErrorType::UnrealisticState;
    }
    // the inner algorithm's timers, accumulated over every solve
    std::tie(timer_Fx_, timer_solve_, timer_check_, std::ignore) = inner_.algo().get_timers();
    timer_total_nr_ += timer.duration();
    return res;
}

template<class Inner>
bool OuterLoopAlgo<Inner>::_solve(int max_iter, real_type tol, bool & need_init, bool check_unrealistic)
{
    inner_.apply_state(state_);
    bool converged = inner_.newton(max_iter, tol, need_init);
    // the first iteration analyzed (or tried to): every later solve only refactorizes. A
    // solve that converged in zero iterations factorized nothing, so it does not count.
    nr_iter_ = inner_.algo().get_nb_iter();
    if (nr_iter_ > 0) need_init = false;
    inner_.finalize();
    // the results, here: what LSGrid and the loops read
    V_ = inner_.algo().get_V();
    Va_ = inner_.algo().get_Va();
    Vm_ = inner_.algo().get_Vm();
    mis_bus_ = inner_.algo().get_bus_mismatch();
    err_ = inner_.algo().get_error();
    // the residual is published against the units' targets: the reactive power a loop froze
    // into the algorithm's (a bus held at a limit) is its units' output, which the results
    // and the physical checks read off this residual. Overwritten by the next solve's copy.
    // The active part is not touched: the targets a loop moved are published as such
    // (get_outer_target_p).
    if (mis_bus_.size() == Sbus_.size() && Sbus_target_.size() == Sbus_.size()) {
        for (Eigen::Index b = 0; b < Sbus_.size(); ++b) {
            const real_type dq = std::imag(Sbus_(b)) - std::imag(Sbus_target_(b));
            if (dq != 0.) mis_bus_(b) += cplx_type(0., dq);
        }
    }
    stats_.nr_iterations.push_back(nr_iter_);
    // OpenLoadFlow's Newton tests convergence only after a step, so a solve is at least one
    // iteration there, and "a pass added Newton iterations" (whether another pass runs)
    // means "a pass re-solved". Counted the same way here, where a solve already within
    // tolerance takes no step at all.
    total_nr_iterations_ += std::max(nr_iter_, 1);
    if (converged && check_unrealistic && lsgrid_ptr_ != nullptr &&
        is_state_unrealistic(*lsgrid_ptr_, V_, inner_.vm_unknown(), params_)) {
        unrealistic_ = true;
        converged = false;
    }
    return converged;
}

template<class Inner>
OuterDeclaration OuterLoopAlgo<Inner>::_declare()
{
    state_.Sbus = &Sbus_;
    state_.Sbus_init = &Sbus_init_;
    state_.Sbus_target = &Sbus_target_;
    OuterContext ctx = _context(&state_);
    ctx.V = nullptr;
    ctx.Va = nullptr;
    ctx.Vm = nullptr;
    OuterDeclaration decl;
    for (const auto & loop : loops_) loop->declare(ctx, decl);
    return decl;
}

template<class Inner>
OuterContext OuterLoopAlgo<Inner>::_context(OuterState * state)
{
    OuterContext ctx;
    ctx.grid = lsgrid_ptr_;
    ctx.V = &V_;
    ctx.Va = &Va_;
    ctx.Vm = &Vm_;
    ctx.bus_mismatch = &mis_bus_;
    ctx.slack_bus = slack_bus_;
    ctx.state = state;
    ctx.controls = &inner_.controls();
    ctx.trace = &stats_.decisions;
    inner_.fill_context(ctx, controller_q_);
    return ctx;
}
