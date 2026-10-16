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
    solve_cache_.reset();  // a new solve: a new view (the row, its grid)

    // the loops edit these, never the grid's
    Sbus_init_ = Sbus;
    Sbus_target_ = Sbus;
    Sbus_ = Sbus;

    // ... and Ybus, whose values the phase shifters patch (see BranchControl)
    Ybus_ = Ybus;

    // a row resuming from another solve (a batch's "n" one) starts from that solve's loops
    const OuterResumePoint * from = (row_ != nullptr && row_->resume) ? row_->resume_from : nullptr;
    if (from != nullptr) _take_loops(*from);

    // OpenLoadFlow's Newton has one slack bus, the reference (slack_ids(0)): sharing the
    // imbalance is the DistributedSlack loop's job. Any other slack bus of the grid keeps
    // its P equation, as a PV bus -- Base gives it a Vm unknown and a Q equation as well
    // when no local regulator pins it (free_vm_slack_buses). Left as slack buses, they
    // would take whatever P closes the balance, and a target the loop moves onto one of
    // their units would change nothing.
    slack_ref_ = slack_ids.head(std::min<Eigen::Index>(slack_ids.size(), 1));
    const Eigen::Index nb_other_slack = slack_ids.size() - slack_ref_.size();
    pv_with_slack_.resize(pv.size() + nb_other_slack);
    pv_with_slack_.head(pv.size()) = pv;
    pv_with_slack_.tail(nb_other_slack) = slack_ids.tail(nb_other_slack);

    bool need_init = false;
    if (!inner_.setup(Ybus_, V, Sbus_, slack_ref_, slack_weights, pv_with_slack_, pq, need_init,
                      [this](){ _declare(); })) {
        err_ = inner_.algo().get_error();
        timer_total_nr_ += timer.duration();
        return false;
    }

    slack_bus_ = slack_ids.size() > 0 ? slack_ids(0) : -1;
    injections_.Sbus = &Sbus_;
    injections_.Sbus_init = &Sbus_init_;
    injections_.Sbus_target = &Sbus_target_;
    injections_.gen_target_p.clear();
    injections_.storage_target_p.clear();
    inner_.controls().reset_states();
    // a controlled transformer the solve does not have (a batch row's contingency, or one the
    // grid opened since the controls were built) holds nothing in Ybus: nothing of it may be
    // patched in this solve, a tap carried from another solve included
    if (lsgrid_ptr_ != nullptr) {
        const OuterGridView view(*lsgrid_ptr_, row_);
        inner_.drop_out_of_row([&view](int trafo, int b1, int b2) {
            int bus1 = -1, bus2 = -1;
            return view.trafo_closed(trafo, bus1, bus2) && bus1 == b1 && bus2 == b2;
        });
    }
    // a resumed solve (a time series' next step, a batch row from its "n" solve) starts from
    // where that solve left the loops; any other from the grid
    const bool resume = from != nullptr ||
                        (has_carry_ && (row_ != nullptr ? row_->resume : params_.resume_outer_state));
    if (from != nullptr) _take_targets(*from);
    if (resume) _apply_carried_positions();

    // OpenLoadFlow's isNeeded filter, then initialize, both before the first solve
    std::vector<BaseOuterLoop *> active;
    {
        OuterContext ctx = _context();
        ctx.resume = resume;
        if (!resume) {
            ctx.V = nullptr;  // nothing solved yet (resuming: the previous solve's)
            ctx.Va = nullptr;
            ctx.Vm = nullptr;
        }
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
                    OuterContext ctx = _context();
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
        OuterContext ctx = _context();
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
    // what the next solve may resume from
    has_carry_ = res;
    if (res) {
        inner_.phase_tap(carried_phase_tap_);
        inner_.ratio_tap(carried_ratio_tap_);
        inner_.shunt_sections(carried_sections_);
    }
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
    inner_.apply_controls();
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
void OuterLoopAlgo<Inner>::_declare()
{
    injections_.Sbus = &Sbus_;
    injections_.Sbus_init = &Sbus_init_;
    injections_.Sbus_target = &Sbus_target_;
    OuterContext ctx = _context();
    ctx.V = nullptr;
    ctx.Va = nullptr;
    ctx.Vm = nullptr;
    for (const auto & loop : loops_) loop->declare(ctx);
}

template<class Inner>
void OuterLoopAlgo<Inner>::_apply_carried_positions()
{
    OuterControls & controls = inner_.controls();
    const int none = std::numeric_limits<int>::min();
    auto carried = [none](const std::vector<int> & v, int id) {
        return (id >= 0 && static_cast<std::size_t>(id) < v.size()) ? v[static_cast<std::size_t>(id)] : none;
    };
    for (const auto & ts : controls.phase_shifters()) {
        const int pos = carried(carried_phase_tap_, ts.first);
        if (pos != none) ts.second->move_tap(pos);
    }
    for (const auto & tr : controls.ratio_taps()) {
        const int pos = carried(carried_ratio_tap_, tr.first);
        if (pos != none) tr.second->move_tap(pos);
    }
    for (const auto & g : controls.shunt_groups()) {
        for (std::size_t c = 0; c < g.controller_buses.size() && c < g.shunts.size(); ++c) {
            ShuntSectionControl * shunt = controls.shunt_controller(g.controller_buses[c]);
            if (shunt == nullptr) continue;
            std::vector<int> counts;
            for (int s : g.shunts[c]) counts.push_back(carried(carried_sections_, s));
            if (std::find(counts.begin(), counts.end(), none) == counts.end()) shunt->set_sections(counts);
        }
    }
}

template<class Inner>
std::shared_ptr<const OuterResumePoint> OuterLoopAlgo<Inner>::get_outer_resume_point() const
{
    if (!has_carry_) return nullptr;
    auto point = std::make_shared<OuterResumePoint>();
    point->loops.reserve(loops_.size());
    for (const auto & loop : loops_) point->loops.push_back(std::shared_ptr<const BaseOuterLoop>(loop->clone()));
    point->phase_tap = carried_phase_tap_;
    point->ratio_tap = carried_ratio_tap_;
    point->sections = carried_sections_;
    point->gen_target_p = injections_.gen_target_p;
    point->storage_target_p = injections_.storage_target_p;
    point->V = V_;
    point->Va = Va_;
    point->Vm = Vm_;
    return point;
}

template<class Inner>
void OuterLoopAlgo<Inner>::_take_loops(const OuterResumePoint & from)
{
    // the same loops, in the same order: the point was taken on an algorithm given the same
    // list (the batch's, see set_outer_loops)
    bool same = from.loops.size() == loops_.size();
    for (std::size_t k = 0; same && k < loops_.size(); ++k) same = from.loops[k]->name() == loops_[k]->name();
    if (!same) {
        throw std::runtime_error("OuterLoopAlgo: the solve to resume from ran other outer loops than this one.");
    }
    // a clone: a loop's state is per solve, and the point is shared by every worker. A loop
    // re-binds what it holds of the controls in initialize(), so a clone taken on another
    // algorithm is fine
    for (std::size_t k = 0; k < loops_.size(); ++k) loops_[k] = from.loops[k]->clone();
    carried_phase_tap_ = from.phase_tap;
    carried_ratio_tap_ = from.ratio_tap;
    carried_sections_ = from.sections;
    V_ = from.V;
    Va_ = from.Va;
    Vm_ = from.Vm;
}

template<class Inner>
void OuterLoopAlgo<Inner>::_take_targets(const OuterResumePoint & from)
{
    if (lsgrid_ptr_ == nullptr) return;
    // the units start at the targets the loops moved them to there -- but where the row sets
    // its own: that is the row's input
    const OuterGridView view(*lsgrid_ptr_, row_);
    for (std::size_t g = 0; g < from.gen_target_p.size(); ++g) {
        const int gen_id = static_cast<int>(g);
        if (std::isnan(from.gen_target_p[g]) || view.row_sets_gen_target_p(gen_id)) continue;
        const int bus = view.gen_bus_solver(gen_id);
        if (bus >= 0) injections_.move_gen_target_p(view, gen_id, bus, from.gen_target_p[g]);
    }
    for (std::size_t s = 0; s < from.storage_target_p.size(); ++s) {
        const int storage_id = static_cast<int>(s);
        if (std::isnan(from.storage_target_p[s]) || view.row_sets_storage_target_p(storage_id)) continue;
        const int bus = view.storage_bus_solver(storage_id);
        if (bus >= 0) injections_.move_storage_target_p(view, storage_id, bus, from.storage_target_p[s]);
    }
}

template<class Inner>
OuterContext OuterLoopAlgo<Inner>::_context()
{
    OuterContext ctx;
    ctx.grid = lsgrid_ptr_;
    ctx.V = &V_;
    ctx.Va = &Va_;
    ctx.Vm = &Vm_;
    ctx.bus_mismatch = &mis_bus_;
    ctx.slack_bus = slack_bus_;
    ctx.injections = &injections_;
    ctx.controls = &inner_.controls();
    ctx.trace = &stats_.decisions;
    ctx.row = row_;
    ctx.solve_cache = &solve_cache_;
    inner_.fill_context(ctx, controller_q_);
    return ctx;
}
