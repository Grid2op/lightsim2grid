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

bool DistributedSlackLoop::_is_needed(const OuterContext & ctx) const
{
    // a unit off the slack bus to take the mismatch: the slack bus' own units would only take
    // what it already absorbs (see the class comment)
    return _has_participant(ctx.view(), ctx.slack_bus);
}

void DistributedSlackLoop::_initialize(OuterContext & ctx)
{
    units_.clear();
    unit_solver_bus_.clear();
    collect_units(ctx.view(), units_, unit_solver_bus_);
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
    _share_lost_p(ctx);
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
    slack_redistribution::distribute(units_, lost, params_.p_residue_eps_mw, new_mw, saturated);
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

AlgoConfig DistributedSlackLoop::_get_params() const
{
    AlgoConfig cfg;
    cfg.int_params = {params_.fail_on_residue ? 1 : 0};
    cfg.real_params = {static_cast<double>(params_.slack_bus_p_max_mismatch_mw),
                       static_cast<double>(params_.p_residue_eps_mw),
                       static_cast<double>(params_.moved_fraction)};
    return cfg;
}

}  // namespace ls2g
