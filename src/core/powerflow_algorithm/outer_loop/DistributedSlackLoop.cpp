// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DistributedSlackLoop.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "LSGrid.hpp"

namespace ls2g {

namespace {

// the units taking part in the slack, in the row the view shows: connected, in the solved
// grid, with a "can participate" weight. A unit on a bus the row strands (masked) is not in
// the solved grid: what it would take never reaches the slack bus. Storage targets are in the
// load convention, their P limits in the generator one.
void collect_units(const OuterGridView & view,
                   std::vector<slack_redistribution::Participant> & units,
                   std::vector<int> & solver_bus)
{
    const GeneratorContainer & gens = view.grid().get_generators();
    const StorageContainer & storages = view.grid().get_storages();
    const GlobalBusIdVect & to_me = view.id_solver_to_me();
    auto add = [&](slack_redistribution::UnitKind kind, int el_id, int bus_solver, real_type weight,
                   real_type injection_mw, real_type min_p, real_type max_p){
        if(!(weight > 0.) || bus_solver < 0 || view.is_masked(bus_solver)) return;
        slack_redistribution::Participant unit;
        unit.kind = kind;
        unit.el_id = el_id;
        unit.bus = to_me[bus_solver].cast_int();
        unit.injection_mw = injection_mw;
        unit.weight = weight;
        unit.min_p_mw = min_p;
        unit.max_p_mw = max_p;
        unit.in_slack = false;
        units.push_back(unit);
        solver_bus.push_back(bus_solver);
    };
    for(int el_id = 0; el_id < gens.nb(); ++el_id){
        add(slack_redistribution::UnitKind::GENERATOR, el_id, view.gen_bus_solver(el_id),
            gens.get_can_participate_slack_weight(el_id), view.gen_target_p(el_id),
            gens.get_min_p(el_id), gens.get_max_p(el_id));
    }
    for(int el_id = 0; el_id < storages.nb(); ++el_id){
        add(slack_redistribution::UnitKind::STORAGE, el_id, view.storage_bus_solver(el_id),
            storages.get_can_participate_slack_weight(el_id), -view.storage_target_p(el_id),
            storages.get_min_p(el_id), storages.get_max_p(el_id));
    }
}

// the solver bus of grid bus `bus_me` in the solve's labelling, -1 when it is not solved
int solver_bus_of(const OuterGridView & view, int bus_me)
{
    const SolverBusIdVect & to_solver = view.id_me_to_solver();
    if (bus_me < 0 || bus_me >= static_cast<int>(to_solver.size())) return -1;
    const int bus = to_solver[bus_me].cast_int();
    return bus == BaseConstants::_deactivated_bus_id ? -1 : bus;
}

}  // namespace

DistributedSlackLoop::DistributedSlackLoop() : DistributedSlackLoop(Params()) {}

DistributedSlackLoop::DistributedSlackLoop(const Params & params) : params_(params)
{
    if(!(params.slack_bus_p_max_mismatch_mw >= 0.)){
        throw std::runtime_error("DistributedSlackLoop: slack_bus_p_max_mismatch_mw must be >= 0.");
    }
    if(!(params.p_residue_eps_mw >= 0.)){
        throw std::runtime_error("DistributedSlackLoop: p_residue_eps_mw must be >= 0.");
    }
    if(!(params.moved_fraction > 0.)){
        throw std::runtime_error("DistributedSlackLoop: moved_fraction must be > 0.");
    }
}

bool DistributedSlackLoop::_in_newton(const OuterContext & ctx)
{
    return ctx.controls != nullptr && ctx.controls->slack_share() != nullptr;
}

void DistributedSlackLoop::_declare(const OuterContext & ctx) const
{
    if (ctx.controls == nullptr || !ctx.controls->supports_slack_share()) return;
    // the bus every unit flagged "can participate" may share the slack from: its own in the
    // grid, connected or not in the solve that declares, and the one this solve places it on.
    // A batch declares once for all its rows.
    const OuterGridView view = ctx.view();
    const GeneratorContainer & gens = view.grid().get_generators();
    const StorageContainer & storages = view.grid().get_storages();
    std::vector<int> buses;
    const auto add = [&buses](int bus){ if (bus >= 0) buses.push_back(bus); };
    for (int el_id = 0; el_id < gens.nb(); ++el_id) {
        if (!(gens.get_can_participate_slack_weight(el_id) > 0.)) continue;
        add(solver_bus_of(view, gens.get_bus_id()(el_id).cast_int()));
        add(view.gen_bus_solver(el_id));
    }
    for (int el_id = 0; el_id < storages.nb(); ++el_id) {
        if (!(storages.get_can_participate_slack_weight(el_id) > 0.)) continue;
        add(solver_bus_of(view, storages.get_bus_id()(el_id).cast_int()));
        add(view.storage_bus_solver(el_id));
    }
    ctx.controls->reserve_slack_share(buses);
}

bool DistributedSlackLoop::_is_needed(const OuterContext & ctx) const
{
    // a unit off the slack bus to take the mismatch: the slack bus' own units would only take
    // what it already absorbs (see the class comment)
    return _has_participant(ctx.view(), ctx.slack_bus);
}

void DistributedSlackLoop::_initialize(OuterContext & ctx)
{
    // the units the solve this one resumes from left out of the Newton's share (this object is
    // that solve's loop, or a clone of it): they stay out, rather than being pushed past their
    // bound again by the first solve. Matched by unit: a batch row may solve other units.
    std::vector<std::pair<std::pair<slack_redistribution::UnitKind, int>, char> > was_parked;
    if (ctx.resume) {
        for (std::size_t k = 0; k < units_.size() && k < parked_.size() && k < parked_high_.size(); ++k) {
            if (parked_[k]) was_parked.emplace_back(std::make_pair(units_[k].kind, units_[k].el_id), parked_high_[k]);
        }
    }
    units_.clear();
    unit_solver_bus_.clear();
    collect_units(ctx.view(), units_, unit_solver_bus_);
    sn_mva_ = ctx.grid->get_sn_mva();
    // a unit starts from its target in force: the one a solve resumed from left it at (a batch
    // row from its "n" solve, as OpenLoadFlow's contingencies start from the targetP of the
    // pre-contingency solve), the row's or the grid's otherwise
    if(ctx.injections != nullptr){
        const OuterGridView view = ctx.view();
        for(auto & unit : units_){
            unit.injection_mw = unit.kind == slack_redistribution::UnitKind::GENERATOR
                ? ctx.injections->gen_target_p_in_force(view, unit.el_id)
                : -ctx.injections->storage_target_p_in_force(view, unit.el_id);
        }
    }
    current_mw_.resize(units_.size());
    for(std::size_t k = 0; k < units_.size(); ++k) current_mw_[k] = units_[k].injection_mw;
    lost_p_saturated_.clear();
    _share_lost_p(ctx);
    // every unit takes its share in the Newton to start with, but those a resumed solve kept out,
    // and those the row's lost power put at a bound: where the losses grow with the contingency,
    // the first solve would push them past it
    parked_.assign(units_.size(), 0);
    parked_high_.assign(units_.size(), 0);
    if (_in_newton(ctx)) {
        _keep_parked(ctx, was_parked);
        for (std::size_t k = 0; k < units_.size() && k < lost_p_saturated_.size(); ++k) {
            if (!lost_p_saturated_[k]) continue;
            parked_[k] = 1;
            parked_high_[k] = lost_p_high_ ? 1 : 0;
        }
        _set_newton_weights(ctx);
    }
}

real_type DistributedSlackLoop::_bound(std::size_t k, bool high) const
{
    const slack_redistribution::Participant & unit = units_[k];
    if (high) {
        real_type hi = unit.max_p_mw;
        if (unit.injection_mw < 0.) hi = std::isnan(hi) ? 0. : std::min(hi, static_cast<real_type>(0.));
        return hi;
    }
    real_type lo = unit.min_p_mw;
    if (!(unit.injection_mw < 0.)) lo = std::isnan(lo) ? 0. : std::max(lo, static_cast<real_type>(0.));
    return lo;
}

void DistributedSlackLoop::_keep_parked(
        OuterContext & ctx,
        const std::vector<std::pair<std::pair<slack_redistribution::UnitKind, int>, char> > & was_parked)
{
    // At the bound it was parked at: a time series' step restarts its units from its own targets.
    // Safe either way: each check shares from the initial targets, so a unit kept out that should
    // take a share gets its weight back, at the cost of a re-solve.
    if (was_parked.empty() || ctx.injections == nullptr) return;
    OuterInjections & state = *ctx.injections;
    const OuterGridView view = ctx.view();
    for (std::size_t k = 0; k < units_.size(); ++k) {
        const std::pair<slack_redistribution::UnitKind, int> unit(units_[k].kind, units_[k].el_id);
        for (const auto & parked : was_parked) {
            if (parked.first != unit) continue;
            const real_type bound = _bound(k, parked.second != 0);
            if (std::isnan(bound)) break;
            if (units_[k].kind == slack_redistribution::UnitKind::GENERATOR) {
                state.move_gen_target_p(view, units_[k].el_id, unit_solver_bus_[k], bound);
            } else {
                state.move_storage_target_p(view, units_[k].el_id, unit_solver_bus_[k], -bound);
            }
            current_mw_[k] = bound;
            parked_[k] = 1;
            parked_high_[k] = parked.second;
            break;
        }
    }
}

void DistributedSlackLoop::_targets_in_force(const OuterContext & ctx, std::vector<real_type> & out) const
{
    const OuterGridView view = ctx.view();
    out.resize(units_.size());
    for (std::size_t k = 0; k < units_.size(); ++k) {
        out[k] = units_[k].kind == slack_redistribution::UnitKind::GENERATOR
            ? ctx.injections->gen_target_p_in_force(view, units_[k].el_id)
            : -ctx.injections->storage_target_p_in_force(view, units_[k].el_id);
    }
}

void DistributedSlackLoop::_set_newton_weights(OuterContext & ctx) const
{
    std::vector<SlackShareControl::Unit> units;
    units.reserve(units_.size());
    for (std::size_t k = 0; k < units_.size(); ++k) {
        SlackShareControl::Unit unit;
        unit.kind = units_[k].kind;
        unit.el_id = units_[k].el_id;
        unit.bus_solver = unit_solver_bus_[k];
        unit.weight = parked_[k] ? 0. : units_[k].weight;
        units.push_back(unit);
    }
    ctx.controls->slack_share()->set_units(units);
}

real_type DistributedSlackLoop::lost_p_to_share_mw(const OuterGridView & view,
                                                   const std::vector<real_type> & resumed_gen_p,
                                                   const std::vector<real_type> & resumed_sto_p)
{
    const OuterRowView * row = view.row();
    if(row == nullptr || std::isnan(row->lost_p_mw)) return std::numeric_limits<real_type>::quiet_NaN();
    // a unit the row loses: off in the row, or on a bus it strands
    const auto lost = [&view](bool on, int bus_solver){ return !on || view.is_masked(bus_solver); };
    real_type res = row->lost_p_mw;
    // only a unit a solve moved has a resumed target: one in the solved grid
    for(std::size_t g = 0; g < resumed_gen_p.size(); ++g){
        const int gen_id = static_cast<int>(g);
        if(std::isnan(resumed_gen_p[g]) || view.row_sets_gen_target_p(gen_id)) continue;
        if(!lost(view.gen_on(gen_id), view.gen_bus_solver(gen_id))) continue;
        res += resumed_gen_p[g] - view.gen_target_p(gen_id);
    }
    for(std::size_t s = 0; s < resumed_sto_p.size(); ++s){
        const int storage_id = static_cast<int>(s);
        if(std::isnan(resumed_sto_p[s]) || view.row_sets_storage_target_p(storage_id)) continue;
        if(!lost(view.storage_on(storage_id), view.storage_bus_solver(storage_id))) continue;
        // load convention: what it consumed more is generation lost less
        res -= resumed_sto_p[s] - view.storage_target_p(storage_id);
    }
    return res;
}

void DistributedSlackLoop::_share_lost_p(OuterContext & ctx)
{
    if(ctx.injections == nullptr || units_.empty()) return;
    const OuterGridView view = ctx.view();
    const OuterResumePoint * from = (ctx.row != nullptr && ctx.row->resume) ? ctx.row->resume_from : nullptr;
    static const std::vector<real_type> none;
    const real_type lost = lost_p_to_share_mw(view, from != nullptr ? from->gen_target_p : none,
                                              from != nullptr ? from->storage_target_p : none);
    if(!(std::abs(lost) > params_.p_residue_eps_mw)) return;

    std::vector<real_type> new_mw;
    std::vector<char> saturated;
    // what cannot be shared (every unit at a bound) stays on the slack bus: the checks then
    // find it in the mismatch, as OpenLoadFlow's DistributedSlack loop does
    const slack_redistribution::Report report =
        slack_redistribution::distribute(units_, lost, params_.p_residue_eps_mw, new_mw, saturated);
    // every unit at a bound: distribute clears the flags, all are
    lost_p_saturated_ = report.all_saturated ? std::vector<char>(units_.size(), 1) : saturated;
    lost_p_high_ = lost > 0.;
    OuterInjections & state = *ctx.injections;
    for(std::size_t k = 0; k < units_.size(); ++k){
        if(units_[k].kind == slack_redistribution::UnitKind::GENERATOR){
            state.move_gen_target_p(view, units_[k].el_id, unit_solver_bus_[k], new_mw[k]);
        } else {
            state.move_storage_target_p(view, units_[k].el_id, unit_solver_bus_[k], -new_mw[k]);
        }
        units_[k].injection_mw = new_mw[k];
        current_mw_[k] = new_mw[k];
    }
    ctx.record("SHARE_LOST_P", true, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
               lost, params_.p_residue_eps_mw);
}

real_type DistributedSlackLoop::_mismatch_mw(const OuterContext & ctx) const
{
    if(ctx.bus_mismatch == nullptr || ctx.slack_bus < 0 || ctx.slack_bus >= ctx.bus_mismatch->size()){
        return std::numeric_limits<real_type>::quiet_NaN();
    }
    // what the slack bus injects beyond its target: what the units must inject more. The
    // same quantity LSGrid::compute_results books on the slack generator: the bus' residual,
    // less what an in-Newton distributed slack carried in its own unknown (one slack bus,
    // so all of it)
    return (std::real((*ctx.bus_mismatch)(ctx.slack_bus)) - ctx.slack_absorbed) * ctx.grid->get_sn_mva();
}

bool DistributedSlackLoop::_triggered(const OuterContext & ctx, real_type & mismatch_mw) const
{
    mismatch_mw = _mismatch_mw(ctx);
    return std::abs(mismatch_mw) > params_.slack_bus_p_max_mismatch_mw && std::abs(mismatch_mw) > params_.p_residue_eps_mw;
}

void DistributedSlackLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    real_type mismatch = 0.;
    if(!_triggered(ctx, mismatch)) return;
    // the same rule as is_needed; the outer-loop mode never gets here otherwise
    if(ctx.is_detection() && !_has_participant(ctx.view(), ctx.slack_bus)) return;
    out.push_back(LimitViolation{ViolationElementType::GRID, -1, 0, LimitViolationType::SLACK_MISMATCH,
                                 mismatch, params_.slack_bus_p_max_mismatch_mw, std::string()});
}

bool DistributedSlackLoop::_has_participant(const OuterGridView & view, int skip_solver_bus)
{
    std::vector<slack_redistribution::Participant> units;
    std::vector<int> solver_bus;
    collect_units(view, units, solver_bus);
    return std::any_of(solver_bus.begin(), solver_bus.end(),
                       [skip_solver_bus](int bus){ return bus != skip_solver_bus; });
}

OuterLoopStatus DistributedSlackLoop::_check(OuterContext & ctx)
{
    if (_in_newton(ctx)) return _check_in_newton(ctx);
    real_type mismatch = 0.;
    const bool triggered = _triggered(ctx, mismatch);
    ctx.record("DISTRIBUTE", triggered, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
               mismatch, std::max(params_.slack_bus_p_max_mismatch_mw, params_.p_residue_eps_mw));
    if(!triggered) return OuterLoopStatus::STABLE;

    // OpenLoadFlow shares the cumulative mismatch from the initial targets every time: what
    // the units already took is given back first
    real_type remaining = mismatch;
    for(std::size_t k = 0; k < units_.size(); ++k) remaining += current_mw_[k] - units_[k].injection_mw;

    std::vector<real_type> new_mw;
    std::vector<char> saturated;
    const slack_redistribution::Report report = slack_redistribution::distribute(
        units_, remaining, params_.p_residue_eps_mw, new_mw, saturated);
    // with no unit at all, nothing was shared and everything is left
    const real_type residue = units_.empty() ? remaining : report.not_distributed_mw;
    if(params_.fail_on_residue && std::abs(residue) > params_.p_residue_eps_mw) {
        ctx.record("FAIL_RESIDUE", true, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
                   residue, params_.p_residue_eps_mw);
        return OuterLoopStatus::FAILED;
    }

    OuterInjections & state = *ctx.injections;
    const OuterGridView view = ctx.view();
    real_type moved = 0.;
    for(std::size_t k = 0; k < units_.size(); ++k){
        moved += std::abs(units_[k].kind == slack_redistribution::UnitKind::GENERATOR
            ? state.move_gen_target_p(view, units_[k].el_id, unit_solver_bus_[k], new_mw[k])
            : state.move_storage_target_p(view, units_[k].el_id, unit_solver_bus_[k], -new_mw[k]));
        current_mw_[k] = new_mw[k];
    }
    // OpenLoadFlow's PreviousStateInfo.moved, with its margin against rounding
    const bool unstable = moved > params_.moved_fraction * params_.p_residue_eps_mw;
    ctx.record("UNITS_MOVED", unstable, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
               moved, params_.moved_fraction * params_.p_residue_eps_mw);
    return unstable ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

OuterLoopStatus DistributedSlackLoop::_check_in_newton(OuterContext & ctx)
{
    // where the solve left the units: the driver moved their targets by their share
    _targets_in_force(ctx, current_mw_);
    // the same cumulative sharing from the initial targets as around a single-slack Newton; what
    // the reference bus kept is not nothing only when every unit is parked
    real_type mismatch = _mismatch_mw(ctx);
    if (std::isnan(mismatch)) mismatch = 0.;
    real_type remaining = mismatch;
    for (std::size_t k = 0; k < units_.size(); ++k) remaining += current_mw_[k] - units_[k].injection_mw;

    std::vector<real_type> new_mw;
    std::vector<char> saturated;
    const slack_redistribution::Report report = slack_redistribution::distribute(
        units_, remaining, params_.p_residue_eps_mw, new_mw, saturated);
    // every unit at a bound: none of them shares any more, the reference bus takes the rest
    if (report.all_saturated) saturated.assign(units_.size(), 1);
    const real_type residue = units_.empty() ? remaining : report.not_distributed_mw;
    if (params_.fail_on_residue && std::abs(residue) > params_.p_residue_eps_mw) {
        ctx.record("FAIL_RESIDUE", true, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
                   residue, params_.p_residue_eps_mw);
        return OuterLoopStatus::FAILED;
    }

    OuterInjections & state = *ctx.injections;
    const OuterGridView view = ctx.view();
    real_type moved = 0.;
    for (std::size_t k = 0; k < units_.size(); ++k) {
        const bool is_gen = units_[k].kind == slack_redistribution::UnitKind::GENERATOR;
        if (saturated[k] != parked_[k]) {
            // a unit at a bound leaves the Newton's share; one back inside its bounds rejoins it
            const real_type bound = new_mw[k];
            ctx.record(saturated[k] ? "PARK" : "UNPARK", true,
                       is_gen ? ViolationElementType::GENERATOR : ViolationElementType::STORAGE, units_[k].el_id,
                       remaining > 0. ? LimitViolationType::HIGH_P : LimitViolationType::LOW_P, current_mw_[k], bound);
        }
        moved += std::abs(is_gen
            ? state.move_gen_target_p(view, units_[k].el_id, unit_solver_bus_[k], new_mw[k])
            : state.move_storage_target_p(view, units_[k].el_id, unit_solver_bus_[k], -new_mw[k]));
        current_mw_[k] = new_mw[k];
    }
    parked_.assign(saturated.begin(), saturated.end());
    parked_high_.assign(units_.size(), remaining > 0. ? 1 : 0);
    _set_newton_weights(ctx);
    // OpenLoadFlow's PreviousStateInfo.moved, with its margin against rounding: the units moved
    // only when some unit was pushed past a bound, or came back
    const real_type threshold = params_.moved_fraction * params_.p_residue_eps_mw;
    bool unstable = moved > threshold;
    ctx.record("UNITS_MOVED", unstable, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
               moved, threshold);
    // ... or when the Newton's share, moved into the targets after the solve, took a unit's
    // target Q with it: OpenLoadFlow re-solves after every move of the targets, and the solve
    // has not seen that one
    const real_type unsolved_q = state.unsolved_q_mvar;
    state.unsolved_q_mvar = 0.;
    if (unsolved_q > threshold) {
        ctx.record("TARGET_Q_FOLLOWED", true, ViolationElementType::GRID, -1, LimitViolationType::SLACK_MISMATCH,
                   unsolved_q, threshold);
        unstable = true;
    }
    return unstable ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

AlgoConfig DistributedSlackLoop::_get_params() const
{
    AlgoConfig cfg;
    cfg.int_params = {params_.fail_on_residue ? 1 : 0};
    cfg.real_params = {static_cast<double>(params_.slack_bus_p_max_mismatch_mw),
                       static_cast<double>(params_.p_residue_eps_mw),
                       static_cast<double>(params_.moved_fraction)};
    return cfg;
}

void DistributedSlackLoop::_adjoint_terms(OuterAdjointTerms & out) const
{
    // a unit parked at a bound sits there whatever its set-point (the Newton gives it no
    // share, and the sharing from the initial targets puts it back at the bound) -- unless
    // that "bound" is its own set-point: a unit already beyond a limit is never pulled back
    // to it (slack_redistribution::distribute), it stays where its target puts it
    if (sn_mva_ > 0.) out.inv_sn = 1. / sn_mva_;
    for (std::size_t k = 0; k < units_.size() && k < parked_.size() && k < unit_solver_bus_.size(); ++k) {
        if (!parked_[k] || units_[k].kind != slack_redistribution::UnitKind::GENERATOR) continue;
        if (k < current_mw_.size() && current_mw_[k] == units_[k].injection_mw) continue;
        out.fixed_p.emplace_back(units_[k].el_id, unit_solver_bus_[k]);
    }
}

}  // namespace ls2g
