// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "PhaseControlLoop.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "LSGrid.hpp"
#include "powerflow_algorithm/NRSystem.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver : ctx.grid->id_me_to_ac_solver();
}

// the solver buses of a branch connected at both ends, false otherwise
template<class Container>
bool closed_buses(const Container & branches, int el, const SolverBusIdVect & to_solver, int & b1, int & b2)
{
    if (!branches.get_status_global()[el] || !branches.get_status_side_1()[el] || !branches.get_status_side_2()[el]) {
        return false;
    }
    b1 = to_solver[branches.get_bus_side_1(el).cast_int()].cast_int();
    b2 = to_solver[branches.get_bus_side_2(el).cast_int()].cast_int();
    return b1 >= 0 && b2 >= 0;
}

int find_root(std::vector<int> & parent, int x)
{
    while (parent[static_cast<std::size_t>(x)] != x) {
        parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        x = parent[static_cast<std::size_t>(x)];
    }
    return x;
}

// whether buses b1 and b2 stay connected through the lines and transformers once
// transformer `skip` is removed (OpenLoadFlow's connectivity: its AC branches)
bool connected_without(const LSGrid & grid, const SolverBusIdVect & to_solver, int skip, int b1, int b2)
{
    std::vector<int> parent(to_solver.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto join = [&](int a, int b) {
        if (a < 0 || b < 0 || static_cast<std::size_t>(a) >= parent.size() || static_cast<std::size_t>(b) >= parent.size()) return;
        parent[static_cast<std::size_t>(find_root(parent, a))] = find_root(parent, b);
    };
    int x1, x2;
    const LineContainer & lines = grid.get_lines();
    for (int el = 0; el < lines.nb(); ++el) {
        if (closed_buses(lines, el, to_solver, x1, x2)) join(x1, x2);
    }
    const TrafoContainer & trafos = grid.get_trafos();
    for (int el = 0; el < trafos.nb(); ++el) {
        if (el != skip && closed_buses(trafos, el, to_solver, x1, x2)) join(x1, x2);
    }
    return find_root(parent, b1) == find_root(parent, b2);
}

// a current in A from pu of the base of a bus of nominal voltage `vn_kv`
real_type to_amps(real_type i_pu, real_type sn_mva, real_type vn_kv)
{
    return i_pu * sn_mva * 1000. / (std::sqrt(3.) * vn_kv);
}

}  // namespace

std::vector<PhaseControlLoop::Shifter> PhaseControlLoop::shifters(const LSGrid & grid)
{
    std::vector<Shifter> res;
    const TrafoContainer & trafos = grid.get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    const SolverBusIdVect & to_solver = grid.id_me_to_ac_solver();
    if (ptc.nb() != trafos.nb() || to_solver.size() == 0) return res;
    for (int t = 0; t < trafos.nb(); ++t) {
        if (!ptc.has(t) || !ptc.regulating(t)) continue;
        const RegulationMode mode = ptc.mode(t);
        if (mode != RegulationMode::ACTIVE_POWER && mode != RegulationMode::CURRENT_LIMITER) continue;
        const int side = ptc.regulated(t);
        if (side != 1 && side != 2) continue;
        int b1, b2;
        if (!closed_buses(trafos, t, to_solver, b1, b2) || b1 == b2) continue;
        const bool regulates = mode != RegulationMode::ACTIVE_POWER || connected_without(grid, to_solver, t, b1, b2);
        res.push_back(Shifter{t, mode, side, regulates});
    }
    return res;
}

void PhaseControlLoop::_declare(const OuterContext & ctx, OuterDeclaration & /*decl*/) const
{
    for (const Shifter & s : shifters(*ctx.grid)) {
        ctx.controls->reserve_phase_shifter(s.trafo, s.mode == RegulationMode::ACTIVE_POWER && s.regulates);
    }
}

bool PhaseControlLoop::_is_needed(const OuterContext & /*ctx*/) const
{
    return true;  // OpenLoadFlow's: the loop always runs, stable at once without a phase shifter
}

void PhaseControlLoop::_initialize(OuterContext & ctx)
{
    shifters_ = shifters(*ctx.grid);
    // the active power controllers regulate from the first solve on
    for (const Shifter & s : shifters_) {
        if (s.mode != RegulationMode::ACTIVE_POWER || !s.regulates) continue;
        PhaseShifterControl * control = ctx.controls->phase_shifter(s.trafo);
        if (control != nullptr) control->set_control(true);
    }
}

OuterLoopStatus PhaseControlLoop::_check(OuterContext & ctx)
{
    // nothing to do without a phase shifter the inner algorithm can act on
    bool any = false;
    for (const Shifter & s : shifters_) any = any || ctx.controls->phase_shifter(s.trafo) != nullptr;
    if (!any) return OuterLoopStatus::STABLE;
    // the control of a shifter the solve handles, null otherwise
    auto handled = [&](const Shifter & s) -> PhaseShifterControl * {
        PhaseShifterControl * control = ctx.controls->phase_shifter(s.trafo);
        return control != nullptr && control->handled() ? control : nullptr;
    };
    const TapChangers & ptc = ctx.grid->get_trafos().get_tap_changers(true);
    if (ctx.iteration == 0) {
        // the active power controllers are switched off, their shift rounded to the
        // closest tap (the current one unless another is strictly closer)
        for (const Shifter & s : shifters_) {
            if (s.mode != RegulationMode::ACTIVE_POWER || !s.regulates) continue;
            PhaseShifterControl * phase = handled(s);
            if (phase == nullptr) continue;
            phase->set_control(false);
            const real_type a = phase->shift();
            int best = phase->position();
            real_type best_distance = std::abs(a - ptc.alpha_at(s.trafo, best));
            for (int pos = ptc.low_tap(s.trafo); pos <= ptc.high_tap(s.trafo); ++pos) {
                const real_type distance = std::abs(a - ptc.alpha_at(s.trafo, pos));
                if (distance < best_distance) {
                    best = pos;
                    best_distance = distance;
                }
            }
            phase->move_tap(best);
            // the shift (rad) and the one of the tap it is rounded to
            ctx.record("ROUND_TAP", true, ViolationElementType::TRAFO, s.trafo, LimitViolationType::PHASE_CONTROL_P,
                       a, ptc.alpha_at(s.trafo, best));
        }
        // OpenLoadFlow re-solves whenever there is a phase shifter, limiters included
        return OuterLoopStatus::UNSTABLE;
    }
    // the current limiters above their limit move one tap, the way that lowers the current
    bool moved = false;
    const TrafoContainer & trafos = ctx.grid->get_trafos();
    const Eigen::Ref<const RealVect> vn_kv = ctx.grid->get_bus_vn_kv();
    for (const Shifter & s : shifters_) {
        if (s.mode != RegulationMode::CURRENT_LIMITER) continue;
        PhaseShifterControl * phase = handled(s);
        if (phase == nullptr) continue;
        real_type i_pu, di_da;
        phase->current(s.side, i_pu, di_da);
        const int bus = (s.side == 1 ? trafos.get_bus_side_1(s.trafo) : trafos.get_bus_side_2(s.trafo)).cast_int();
        const real_type i_a = to_amps(i_pu, ctx.grid->get_sn_mva(), vn_kv(bus));
        if (!(ptc.target(s.trafo) < i_a)) {
            ctx.record("KEPT_TAP", false, ViolationElementType::TRAFO, s.trafo, LimitViolationType::PHASE_LIMITER_CURRENT,
                       i_a, ptc.target(s.trafo));
            continue;
        }
        // PiModelArray.shiftOneTapPositionToChangeA1: increase (decrease) the shift when a
        // larger (smaller) one lowers the current
        const bool increase = !(di_da > 0.);
        const real_type a = phase->shift();
        int pos = phase->position();
        const int old_pos = pos;
        if (pos < ptc.high_tap(s.trafo)) {
            const real_type next = ptc.alpha_at(s.trafo, pos + 1);
            if ((increase && next > a) || (!increase && next < a)) ++pos;
        }
        if (pos > ptc.low_tap(s.trafo)) {
            const real_type previous = ptc.alpha_at(s.trafo, pos - 1);
            if ((increase && previous > a) || (!increase && previous < a)) --pos;
        }
        ctx.record("MOVE_TAP", pos != old_pos, ViolationElementType::TRAFO, s.trafo,
                   LimitViolationType::PHASE_LIMITER_CURRENT, i_a, ptc.target(s.trafo));
        if (pos != old_pos) {
            phase->move_tap(pos);
            moved = true;
        }
    }
    return moved ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

void PhaseControlLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    if (ctx.V == nullptr) return;
    const LSGrid & grid = *ctx.grid;
    const TrafoContainer & trafos = grid.get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    const SolverBusIdVect & to_solver = solver_map(ctx);
    const Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    const real_type sn = grid.get_sn_mva();
    for (const Shifter & s : shifters(grid)) {
        if (!s.regulates) continue;
        int b1, b2;
        if (!closed_buses(trafos, s.trafo, to_solver, b1, b2) || b1 >= ctx.V->size() || b2 >= ctx.V->size()) continue;
        const TrafoInfo info = trafos[s.trafo];
        const cplx_type V1 = (*ctx.V)(b1);
        const cplx_type V2 = (*ctx.V)(b2);
        const cplx_type I = s.side == 1 ? info.yac_eff_11 * V1 + info.yac_eff_12 * V2
                                        : info.yac_eff_21 * V1 + info.yac_eff_22 * V2;
        const cplx_type V = s.side == 1 ? V1 : V2;
        if (s.mode == RegulationMode::ACTIVE_POWER) {
            const real_type p_mw = std::real(V * std::conj(I)) * sn;
            const real_type target = ptc.target(s.trafo);
            const real_type threshold = std::max(ptc.deadband(s.trafo), ctx.tol_mw);
            if (std::abs(p_mw - target) > threshold) {
                out.push_back(LimitViolation{ViolationElementType::TRAFO, s.trafo, s.side,
                                             LimitViolationType::PHASE_CONTROL_P, p_mw, target, info.name});
            }
        } else {
            const int bus = (s.side == 1 ? trafos.get_bus_side_1(s.trafo) : trafos.get_bus_side_2(s.trafo)).cast_int();
            const real_type i_a = to_amps(std::abs(I), sn, vn_kv(bus));
            if (i_a > ptc.target(s.trafo)) {
                out.push_back(LimitViolation{ViolationElementType::TRAFO, s.trafo, s.side,
                                             LimitViolationType::PHASE_LIMITER_CURRENT, i_a, ptc.target(s.trafo),
                                             info.name});
            }
        }
    }
}

}  // namespace ls2g
