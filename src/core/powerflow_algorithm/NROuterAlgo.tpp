// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

template<class LinearSolver>
bool NROuterAlgo<LinearSolver>::compute_pf(
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
    if (!this->is_linear_solver_valid()) return false;

    this->reset_timer();
    this->err_ = ErrorType::NoError;
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
    this->_system.set_branch_mutable_ybus(&Ybus_);

    bool need_init = false;
    if (!this->_setup(Ybus_, V, Sbus_, slack_ids, slack_weights, pv, pq, need_init)) {
        this->timer_total_nr_ += timer.duration();
        return false;
    }

    slack_bus_ = slack_ids.size() > 0 ? slack_ids(0) : -1;
    state_.Sbus = &Sbus_;
    state_.Sbus_init = &Sbus_init_;
    state_.Sbus_target = &Sbus_target_;
    state_.gen_target_p.clear();
    state_.storage_target_p.clear();
    state_.hvdc_status.clear();
    state_.svc_target_vm.clear();
    state_.pq_buses.clear();
    state_.vm_set.clear();
    state_.controller_hold_q.clear();
    state_.phase_tap.clear();
    state_.phase_control.clear();
    state_.ratio_tap.clear();
    state_.ratio_control.clear();
    state_.suspended_buses.clear();
    OuterState & state = state_;

    // OpenLoadFlow's isNeeded filter, then initialize, both before the first solve
    std::vector<BaseOuterLoop *> active;
    {
        OuterContext ctx = _context(&state);
        ctx.V = nullptr;  // nothing solved yet
        ctx.Va = nullptr;
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
                    status = loop->check(ctx);
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
                    this->lsgrid_ptr_ != nullptr) {
                    // the deferred check, once the last loop able to fix it is done
                    if (is_state_unrealistic(*this->lsgrid_ptr_, this->V_, _vm_unknown(), params_)) {
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
        if (stats_.status == OuterLoopStatus::FAILED) this->err_ = ErrorType::OuterLoopFailed;
        else if (stats_.status == OuterLoopStatus::UNSTABLE) this->err_ = ErrorType::TooManyIterations;
    } else if (unrealistic_) {
        this->err_ = ErrorType::UnrealisticState;
    }
    this->timer_total_nr_ += timer.duration();
    return res;
}

template<class LinearSolver>
bool NROuterAlgo<LinearSolver>::_solve(int max_iter, real_type tol, bool & need_init, bool check_unrealistic)
{
    // what the loops changed outside the injection: the hvdc lines' regimes
    if (!state_.hvdc_status.empty()) {
        this->_system.set_hvdc_status_override(state_.hvdc_status, OuterState::HVDC_KEEP);
    }
    // ... and the idle standby SVCs switched on (a release is never undone in a solve)
    if (!state_.svc_target_vm.empty()) {
        this->_system.release_held_svcs(state_.svc_target_vm);
    }
    // the voltage controllers a loop holds (none: the plan's own state)
    this->_system.set_held_voltage_controllers(state_.controller_hold_q);
    // the phase taps a loop moved, then the shifts it lets the Newton solve for
    BranchControl * phase = this->_system.branch_control();
    if (phase != nullptr) {
        for (std::size_t t = 0; t < state_.phase_tap.size(); ++t) {
            if (state_.phase_tap[t] != OuterState::TAP_KEEP) phase->set_tap(static_cast<int>(t), state_.phase_tap[t]);
        }
        for (std::size_t t = 0; t < state_.phase_control.size(); ++t) {
            if (state_.phase_control[t] >= 0) phase->set_control_on(static_cast<int>(t), state_.phase_control[t] == 1);
        }
        // the ratio taps (a move once), then the voltage controls
        for (std::size_t t = 0; t < state_.ratio_tap.size(); ++t) {
            if (state_.ratio_tap[t] == OuterState::TAP_KEEP) continue;
            phase->set_ratio_tap(static_cast<int>(t), state_.ratio_tap[t]);
            state_.ratio_tap[t] = OuterState::TAP_KEEP;
        }
        for (std::size_t t = 0; t < state_.ratio_control.size(); ++t) {
            if (state_.ratio_control[t] >= 0) phase->set_ratio_control_on(static_cast<int>(t), state_.ratio_control[t] == 1);
        }
    }
    // the switchable buses: PV (pinned) unless a loop made them PQ, a caller's on top
    pinned_ = caller_pinned_;
    for (int bus : declared_switchable_) {
        if (!state_.pq_buses.count(bus)) pinned_.push_back(bus);
    }
    this->_system.set_pv_pinned_buses(pinned_);
    // the magnitudes a loop reset (a bus back to PV at its set-point, robust mode)
    if (!state_.vm_set.empty()) {
        std::vector<int> buses;
        std::vector<real_type> vm;
        for (const auto & bv : state_.vm_set) { buses.push_back(bv.first); vm.push_back(bv.second); }
        this->_system.set_vm_at(buses, vm);
        state_.vm_set.clear();
    }
    bool converged = this->_newton(max_iter, tol, need_init);
    // the first iteration analyzed (or tried to): every later solve only refactorizes. A
    // solve that converged in zero iterations factorized nothing, so it does not count.
    if (this->nr_iter_ > 0) need_init = false;
    this->_finalize();
    // the residual is published against the units' targets: the reactive power a loop froze
    // into the algorithm's (a bus held at a limit) is its units' output, which the results
    // and the physical checks read off this residual. Overwritten by the next
    // Newton's first mismatch. The active part is not touched: the targets a loop moved are
    // published as such (get_outer_target_p).
    if (this->mis_bus_.size() == Sbus_.size() && Sbus_target_.size() == Sbus_.size()) {
        for (Eigen::Index b = 0; b < Sbus_.size(); ++b) {
            const real_type dq = std::imag(Sbus_(b)) - std::imag(Sbus_target_(b));
            if (dq != 0.) this->mis_bus_(b) += cplx_type(0., dq);
        }
    }
    stats_.nr_iterations.push_back(this->nr_iter_);
    // OpenLoadFlow's Newton tests convergence only after a step, so a solve is at least one
    // iteration there, and "a pass added Newton iterations" (whether another pass runs)
    // means "a pass re-solved". Counted the same way here, where a solve already within
    // tolerance takes no step at all.
    total_nr_iterations_ += std::max(this->nr_iter_, 1);
    if (converged && check_unrealistic && this->lsgrid_ptr_ != nullptr &&
        is_state_unrealistic(*this->lsgrid_ptr_, this->V_, _vm_unknown(), params_)) {
        unrealistic_ = true;
        converged = false;
    }
    return converged;
}

template<class LinearSolver>
void NROuterAlgo<LinearSolver>::_before_init_topology()
{
    state_.Sbus = &Sbus_;
    state_.Sbus_init = &Sbus_init_;
    state_.Sbus_target = &Sbus_target_;
    OuterContext ctx = _context(&state_);
    ctx.V = nullptr;
    ctx.Va = nullptr;
    OuterDeclaration decl;
    for (const auto & loop : loops_) loop->declare(ctx, decl);
    std::vector<int> switchable = caller_switchable_;
    declared_switchable_ = decl.switchable_vm_buses();
    std::sort(declared_switchable_.begin(), declared_switchable_.end());
    declared_switchable_.erase(std::unique(declared_switchable_.begin(), declared_switchable_.end()),
                               declared_switchable_.end());
    switchable.insert(switchable.end(), declared_switchable_.begin(), declared_switchable_.end());
    this->_system.set_switchable_vm_buses(switchable);  // a set there: duplicates are fine
    this->_system.set_may_hold_voltage_controllers(decl.holds_voltage_controllers());
    this->_system.set_phase_controllers(decl.phase_shifters(), decl.phase_shifter_column());
    std::vector<BranchControl::RatioGroupDecl> groups;
    for (const auto & g : decl.ratio_groups()) {
        BranchControl::RatioGroupDecl d;
        d.bus_solver = g.bus_solver;
        d.target_vm = g.target_vm;
        d.trafos = g.trafos;
        d.solved = g.solved;
        groups.push_back(d);
    }
    this->_system.set_ratio_groups(groups);
}

template<class LinearSolver>
std::vector<bool> NROuterAlgo<LinearSolver>::_vm_unknown() const
{
    const std::vector<int> & vm_col = this->_system.vm_to_J_col();
    std::vector<bool> res(vm_col.size(), false);
    for (std::size_t bus = 0; bus < vm_col.size(); ++bus) res[bus] = vm_col[bus] >= 0;
    // a pinned Q row makes the bus PV: its magnitude is not an unknown
    for (int bus : pinned_) {
        if (bus >= 0 && static_cast<std::size_t>(bus) < res.size()) res[static_cast<std::size_t>(bus)] = false;
    }
    return res;
}
