// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "HvdcAcEmulationLimitsLoop.hpp"

#include "LSGrid.hpp"

namespace ls2g {

namespace {

RealVect angles(const CplxVect & V)
{
    return V.array().arg().matrix();
}

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver : ctx.grid->id_me_to_ac_solver();
}

}  // namespace

hvdc_p_check::HvdcPPlan HvdcAcEmulationLimitsLoop::_ac_emulation_lines(const OuterContext & ctx)
{
    hvdc_p_check::HvdcPPlan plan;
    hvdc_p_check::build_hvdc_p_plan(*ctx.grid, solver_map(ctx), plan);
    // a line the caller froze at a limit (droop off) is not in AC emulation
    hvdc_p_check::HvdcPPlan res;
    for(const auto & line : plan.lines) if(line.frozen_dir == 0) res.lines.push_back(line);
    return res;
}

void HvdcAcEmulationLimitsLoop::_declare(const OuterContext & ctx) const
{
    // every line in AC emulation may saturate -- and so may a line the caller holds at a
    // limit now: putting it back in AC emulation (set_status_droop) changes values only, and
    // does not reach this reservation again
    const HvdcLineContainer & hvdcs = ctx.grid->get_dclines();
    for(int hvdc_id = 0; hvdc_id < hvdcs.nb(); ++hvdc_id){
        if(hvdcs.is_droop_active(hvdc_id)) ctx.controls->reserve_hvdc_regime(hvdc_id);
    }
}

bool HvdcAcEmulationLimitsLoop::_is_needed(const OuterContext & ctx) const
{
    return !_ac_emulation_lines(ctx).empty();
}

void HvdcAcEmulationLimitsLoop::_initialize(OuterContext & ctx)
{
    // every line starts in its linear regime (frozen_dir 0), as in OpenLoadFlow
    plan_ = _ac_emulation_lines(ctx);
}

void HvdcAcEmulationLimitsLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    if(ctx.Va == nullptr && ctx.V == nullptr) return;
    const RealVect Va = ctx.Va != nullptr ? *ctx.Va : angles(*ctx.V);
    if(ctx.is_detection()){
        // the grid's own lines and regimes, in the labelling of the solve that was checked
        hvdc_p_check::HvdcPPlan plan;
        hvdc_p_check::build_hvdc_p_plan(*ctx.grid, solver_map(ctx), plan);
        hvdc_p_check::check_hvdc_p_violations(plan, Va, ctx.tol_mw, ctx.masked, out);
        return;
    }
    hvdc_p_check::check_hvdc_p_violations(plan_, Va, ctx.tol_mw, ctx.masked, out);
}

OuterLoopStatus HvdcAcEmulationLimitsLoop::_check(OuterContext & ctx)
{
    std::vector<LimitViolation> trigger;
    _detect(ctx, trigger);
    if(trigger.empty()) return OuterLoopStatus::STABLE;

    bool changed = false;
    for(const LimitViolation & v : trigger){
        if(v.element_type != ViolationElementType::HVDC) continue;
        int regime;
        if(v.violation_type == LimitViolationType::HIGH_P) regime = v.side == 1 ? 1 : -1;
        else if(v.violation_type == LimitViolationType::HVDC_AC_EMULATION_RELEASE) regime = 0;
        else continue;
        bool acted = false;
        for(auto & line : plan_.lines){
            if(line.hvdc_id != v.element_id || line.frozen_dir == regime) continue;
            HvdcRegimeControl * control = ctx.controls->hvdc_regime(line.hvdc_id);
            if(control == nullptr) continue;
            line.frozen_dir = regime;
            control->set(regime);
            changed = true;
            acted = true;
        }
        ctx.record(regime == 0 ? "RELEASE" : "SATURATE", acted, ViolationElementType::HVDC, v.element_id,
                   v.violation_type, v.value, v.limit);
    }
    return changed ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

}  // namespace ls2g
