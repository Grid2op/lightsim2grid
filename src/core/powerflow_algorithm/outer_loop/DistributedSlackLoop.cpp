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

// the units of one container taking part in the slack: connected, in the solved grid, with a
// "can participate" weight; `sign` turns the container's target into the generator convention
template<class Container>
void collect_units(const Container & container,
                   slack_redistribution::UnitKind kind,
                   real_type sign,
                   const SolverBusIdVect & id_me_to_solver,
                   std::vector<slack_redistribution::Participant> & units,
                   std::vector<int> & solver_bus)
{
    const std::vector<bool> & status = container.get_status();
    const GlobalBusIdVect & bus_id = container.get_bus_id();
    Eigen::Ref<const RealVect> target_p = container.get_target_p();
    for(int el_id = 0; el_id < container.nb(); ++el_id){
        if(!status[el_id]) continue;
        const real_type weight = container.get_can_participate_slack_weight(el_id);
        if(!(weight > 0.)) continue;
        const int bus_me = bus_id(el_id).cast_int();
        if(bus_me == BaseConstants::_deactivated_bus_id) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver == BaseConstants::_deactivated_bus_id) continue;
        slack_redistribution::Participant unit;
        unit.kind = kind;
        unit.el_id = el_id;
        unit.bus = bus_me;
        unit.injection_mw = sign * target_p(el_id);
        unit.weight = weight;
        unit.min_p_mw = container.get_min_p(el_id);
        unit.max_p_mw = container.get_max_p(el_id);
        unit.in_slack = false;
        units.push_back(unit);
        solver_bus.push_back(bus_solver);
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
    return _has_participant(*ctx.grid, ctx.grid->id_me_to_ac_solver());
}

void DistributedSlackLoop::_initialize(OuterContext & ctx)
{
    units_.clear();
    unit_solver_bus_.clear();
    const LSGrid & grid = *ctx.grid;
    const SolverBusIdVect & id_me_to_solver = grid.id_me_to_ac_solver();
    collect_units(grid.get_generators(), slack_redistribution::UnitKind::GENERATOR, 1.,
                  id_me_to_solver, units_, unit_solver_bus_);
    // a storage unit's target is in the load convention, its P limits in the generator one
    collect_units(grid.get_storages(), slack_redistribution::UnitKind::STORAGE, -1.,
                  id_me_to_solver, units_, unit_solver_bus_);
    current_mw_.resize(units_.size());
    for(std::size_t k = 0; k < units_.size(); ++k) current_mw_[k] = units_[k].injection_mw;
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
    // only on a grid that says who would share the slack (see is_needed); the outer-loop
    // mode never gets here otherwise
    if(ctx.is_detection() &&
       !_has_participant(*ctx.grid, ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver
                                                                   : ctx.grid->id_me_to_ac_solver())) return;
    out.push_back(LimitViolation{ViolationElementType::GRID, -1, 0, LimitViolationType::SLACK_MISMATCH,
                                 mismatch, params_.slack_bus_p_max_mismatch_mw, std::string()});
}

bool DistributedSlackLoop::_has_participant(const LSGrid & grid, const SolverBusIdVect & id_me_to_solver)
{
    std::vector<slack_redistribution::Participant> units;
    std::vector<int> solver_bus;
    collect_units(grid.get_generators(), slack_redistribution::UnitKind::GENERATOR, 1.,
                  id_me_to_solver, units, solver_bus);
    if(!units.empty()) return true;
    collect_units(grid.get_storages(), slack_redistribution::UnitKind::STORAGE, -1.,
                  id_me_to_solver, units, solver_bus);
    return !units.empty();
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
    const LSGrid & grid = *ctx.grid;
    const real_type sn_mva = grid.get_sn_mva();
    if(state.gen_target_p.empty() && grid.get_generators().nb() > 0){
        state.gen_target_p.assign(grid.get_generators().nb(), std::numeric_limits<real_type>::quiet_NaN());
    }
    if(state.storage_target_p.empty() && grid.get_storages().nb() > 0){
        state.storage_target_p.assign(grid.get_storages().nb(), std::numeric_limits<real_type>::quiet_NaN());
    }
    real_type moved = 0.;
    for(std::size_t k = 0; k < units_.size(); ++k){
        const real_type delta = new_mw[k] - current_mw_[k];
        moved += std::abs(delta);
        cplx_type ds = {delta / sn_mva, 0.};
        current_mw_[k] = new_mw[k];
        if(units_[k].kind == slack_redistribution::UnitKind::GENERATOR){
            const int gen_id = units_[k].el_id;
            // a unit that does not regulate: its target Q follows its limits at the new P
            // (forceTargetQInReactiveLimits, LfGeneratorImpl.getTargetQ)
            const bool pq = !grid.get_generators().get_voltage_regulator_on(gen_id);
            const real_type q_before = pq ? grid.gen_target_q_at_outer_target(gen_id, state.gen_target_p) : 0.;
            state.gen_target_p[static_cast<std::size_t>(gen_id)] = new_mw[k];
            if(pq) ds += cplx_type(0., (grid.gen_target_q_at_outer_target(gen_id, state.gen_target_p) - q_before) / sn_mva);
        } else {
            state.storage_target_p[units_[k].el_id] = -new_mw[k];
        }
        if(ds != cplx_type(0., 0.)){
            (*state.Sbus)(unit_solver_bus_[k]) += ds;
            if(state.Sbus_target != nullptr) (*state.Sbus_target)(unit_solver_bus_[k]) += ds;
        }
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
