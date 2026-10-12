// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "VoltageMonitoringLoop.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "LSGrid.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver : ctx.grid->id_me_to_ac_solver();
}

bool is_released(const OuterContext & ctx, int svc_id)
{
    const StandbySvcControl * svc = ctx.controls != nullptr ? ctx.controls->standby_svc(svc_id) : nullptr;
    return svc != nullptr && svc->is_on();
}

}  // namespace

std::vector<int> VoltageMonitoringLoop::_held_svcs(const LSGrid & grid)
{
    std::vector<int> res;
    const VoltageControlSolverData & data = grid.get_ac_voltage_control_plan().controllers();
    for(int j = 0; j < data.n_controllers(); ++j){
        if(data.kind(j) == VoltageControlSolverData::SVC && data.is_held(j)) res.push_back(data.elem_id(j));
    }
    std::sort(res.begin(), res.end());
    return res;
}

void VoltageMonitoringLoop::_declare(const OuterContext & ctx) const
{
    // every idle standby SVC may be switched on
    for(int svc_id : _held_svcs(*ctx.grid)) ctx.controls->reserve_standby_svc(svc_id);
}

bool VoltageMonitoringLoop::_is_needed(const OuterContext & ctx) const
{
    const SvcContainer & svcs = ctx.grid->get_svcs();
    for(int svc_id : _held_svcs(*ctx.grid)){
        if(svcs.get_regulated_bus_id(svc_id) == svcs.get_bus_id()(svc_id).cast_int()) return true;
    }
    return false;
}

void VoltageMonitoringLoop::_initialize(OuterContext & ctx)
{
    svc_standby_check::SvcStandbyPlan plan;
    svc_standby_check::build_svc_standby_plan(*ctx.grid, solver_map(ctx), plan);
    const std::vector<int> held = _held_svcs(*ctx.grid);
    monitors_.clear();
    for(const auto & entry : plan.svcs){
        if(std::binary_search(held.begin(), held.end(), entry.svc_id)) monitors_.svcs.push_back(entry);
    }
}

void VoltageMonitoringLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    if(!ctx.vm_checks || ctx.V == nullptr) return;
    if(ctx.is_detection()){
        // every idle SVC the grid flags standby, in the labelling of the solve that was checked
        svc_standby_check::SvcStandbyPlan plan;
        svc_standby_check::build_svc_standby_plan(*ctx.grid, solver_map(ctx), plan);
        svc_standby_check::check_svc_standby_violations(plan, *ctx.V, ctx.tol_vm_pu, ctx.masked, out);
        return;
    }
    // the monitors still idle, compared strictly
    svc_standby_check::SvcStandbyPlan idle;
    for(const auto & entry : monitors_.svcs){
        if(!is_released(ctx, entry.svc_id)) idle.svcs.push_back(entry);
    }
    svc_standby_check::check_svc_standby_violations(idle, *ctx.V, ctx.tol_vm_pu, ctx.masked, out);
}

OuterLoopStatus VoltageMonitoringLoop::_check(OuterContext & ctx)
{
    std::vector<LimitViolation> trigger;
    _detect(ctx, trigger);
    if(trigger.empty()) return OuterLoopStatus::STABLE;

    const SvcContainer & svcs = ctx.grid->get_svcs();
    bool changed = false;
    for(const LimitViolation & v : trigger){
        if(v.element_type != ViolationElementType::SVC) continue;
        real_type target;
        if(v.violation_type == LimitViolationType::LOW_VOLTAGE_SVC_STANDBY) target = svcs.get_standby_low_target_vm_pu(v.element_id);
        else if(v.violation_type == LimitViolationType::HIGH_VOLTAGE_SVC_STANDBY) target = svcs.get_standby_high_target_vm_pu(v.element_id);
        else continue;
        // without a set-point, or without the control a reservation gives, it stays idle
        StandbySvcControl * svc = std::isfinite(target) ? ctx.controls->standby_svc(v.element_id) : nullptr;
        ctx.record("SWITCH_ON", svc != nullptr, ViolationElementType::SVC, v.element_id,
                   v.violation_type, v.value, v.limit);
        if(svc == nullptr) continue;
        svc->switch_on(target);
        changed = true;
    }
    return changed ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

}  // namespace ls2g
