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
    Sbus_ = Sbus;

    bool need_init = false;
    if (!this->_setup(Ybus, V, Sbus_, slack_ids, slack_weights, pv, pq, need_init)) {
        this->timer_total_nr_ += timer.duration();
        return false;
    }

    slack_bus_ = slack_ids.size() > 0 ? slack_ids(0) : -1;
    state_.Sbus = &Sbus_;
    state_.Sbus_init = &Sbus_init_;
    state_.gen_target_p.clear();
    state_.storage_target_p.clear();
    state_.hvdc_status.clear();
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
    bool converged = this->_newton(max_iter, tol, need_init);
    // the first iteration analyzed (or tried to): every later solve only refactorizes. A
    // solve that converged in zero iterations factorized nothing, so it does not count.
    if (this->nr_iter_ > 0) need_init = false;
    this->_finalize();
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
    OuterContext ctx = _context(&state_);
    ctx.V = nullptr;
    ctx.Va = nullptr;
    OuterDeclaration decl;
    for (const auto & loop : loops_) loop->declare(ctx, decl);
    std::vector<int> switchable = caller_switchable_;
    switchable.insert(switchable.end(), decl.switchable_vm_buses().begin(), decl.switchable_vm_buses().end());
    this->_system.set_switchable_vm_buses(switchable);  // a set there: duplicates are fine
}

template<class LinearSolver>
std::vector<bool> NROuterAlgo<LinearSolver>::_vm_unknown() const
{
    const std::vector<int> & vm_col = this->_system.vm_to_J_col();
    std::vector<bool> res(vm_col.size(), false);
    for (std::size_t bus = 0; bus < vm_col.size(); ++bus) res[bus] = vm_col[bus] >= 0;
    return res;
}
