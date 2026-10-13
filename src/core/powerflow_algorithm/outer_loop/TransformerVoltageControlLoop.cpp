// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "TransformerVoltageControlLoop.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

#include "LSGrid.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.view().id_me_to_solver();
}

int find_root(std::vector<int> & parent, int x)
{
    while (parent[static_cast<std::size_t>(x)] != x) {
        parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        x = parent[static_cast<std::size_t>(x)];
    }
    return x;
}

// the phase tap position of transformer `t` in the solve: a phase shifter's (PhaseControl may
// move it), the grid's otherwise
int phase_position(const OuterContext & ctx, int t)
{
    const PhaseShifterControl * phase = ctx.controls->phase_shifter(t);
    if (phase != nullptr && phase->handled()) return phase->position();
    const TapChangers & ptc = ctx.grid->get_trafos().get_tap_changers(true);
    return ptc.has(t) ? ptc.position(t) : 0;
}

// the ratio control of transformer `t` when the solve handles it, null otherwise
RatioTapControl * handled_ratio(const OuterContext & ctx, int t)
{
    RatioTapControl * tap = ctx.controls->ratio_tap(t);
    return tap != nullptr && tap->handled() ? tap : nullptr;
}

}  // namespace

TransformerVoltageControlLoop::TransformerVoltageControlLoop() : TransformerVoltageControlLoop(Params()) {}

TransformerVoltageControlLoop::TransformerVoltageControlLoop(const Params & params) : params_(params)
{
    if (!(params.min_target_deadband_kv >= 0.)) {
        throw std::runtime_error("TransformerVoltageControlLoop: min_target_deadband_kv must be >= 0.");
    }
}

std::vector<TransformerVoltageControlLoop::Group> TransformerVoltageControlLoop::groups(const OuterGridView & view) const
{
    std::vector<Group> res;
    const LSGrid & grid = view.grid();
    const TrafoContainer & trafos = grid.get_trafos();
    const TapChangers & rtc = trafos.get_tap_changers(false);
    const SolverBusLayout & layout = view.layout();
    const SolverBusIdVect & to_solver = layout.id_me_to_solver;
    if (rtc.nb() != trafos.nb() || to_solver.size() == 0) return res;
    const Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    // the buses a generator, an SVC or a VSC station regulates (OpenLoadFlow's GENERATOR and
    // VOLTAGE_SOURCE_CONVERTER controls, of a higher priority than a transformer's)
    std::set<int> generator_controlled;
    for (const SolverBusId & b : layout.bus_pv) generator_controlled.insert(b.cast_int());
    for (const SolverBusId & b : layout.slack_bus_id_solver) generator_controlled.insert(b.cast_int());
    const VoltageControlSolverData & ctrl = layout.voltage_control.controllers();
    for (int g = 0; g < ctrl.n_groups(); ++g) generator_controlled.insert(ctrl.reg_bus(g));

    std::map<int, std::size_t> group_of_bus;
    for (int t = 0; t < trafos.nb(); ++t) {
        if (!rtc.has(t) || !rtc.regulating(t) || rtc.mode(t) != RegulationMode::VOLTAGE) continue;
        const int reg = rtc.regulated(t);
        int b1, b2;
        if (reg < 0 || reg >= static_cast<int>(to_solver.size()) || !view.trafo_closed(t, b1, b2)) continue;
        const int bus = to_solver[reg].cast_int();
        if (bus < 0 || b1 == b2) continue;
        auto it = group_of_bus.find(bus);
        if (it == group_of_bus.end()) {
            Group g;
            g.bus_solver = bus;
            g.bus_grid = reg;
            g.target = rtc.target(t);  // the first one's holds
            g.half_deadband = std::numeric_limits<real_type>::infinity();
            g.hidden = generator_controlled.count(bus) > 0;
            it = group_of_bus.emplace(bus, res.size()).first;
            res.push_back(g);
        }
        Group & g = res[it->second];
        g.trafos.push_back(t);
        if (rtc.deadband(t) > 0.) g.half_deadband = std::min(g.half_deadband, rtc.deadband(t) / 2.);
    }
    for (Group & g : res) {
        if (!std::isfinite(g.half_deadband)) g.half_deadband = params_.min_target_deadband_kv / vn_kv(g.bus_grid) / 2.;
    }
    return res;
}

real_type TransformerVoltageControlLoop::_limit(const OuterGridView & view) const
{
    if (params_.max_controlled_nominal_voltage >= 0.) return params_.max_controlled_nominal_voltage;
    // GeneratorVoltageControlManager.computeDefaultMinNominalVoltageLimit
    real_type res = std::numeric_limits<real_type>::min();
    const Eigen::Ref<const RealVect> vn_kv = view.grid().get_bus_vn_kv();
    const GlobalBusIdVect & to_me = view.id_solver_to_me();
    for (const Group & g : groups_) {
        bool valid = true;
        for (int t : g.trafos) {
            int b1, b2;
            if (!view.trafo_closed(t, b1, b2) ||
                vn_kv(to_me[b1].cast_int()) == vn_kv(to_me[b2].cast_int())) {
                valid = false;
                break;
            }
        }
        if (valid) res = std::max(res, vn_kv(g.bus_grid));
    }
    return res;
}

void TransformerVoltageControlLoop::_declare(const OuterContext & ctx) const
{
    const std::vector<Group> gs = groups(ctx.view());
    if (gs.empty()) return;
    // a hidden group never acts (OpenLoadFlow's getControllerElements keeps the visible
    // controls only): nothing to reserve for it; nor for a group the inner algorithm cannot
    // solve the ratios of (no BranchControl)
    bool any = false;
    for (const Group & g : gs) {
        if (g.hidden) continue;
        if (ctx.controls->reserve_ratio_group(g.bus_solver, g.target, g.trafos, true)) any = true;
    }
    if (!any) return;
    // the generators the INITIAL step may freeze: every one (declaring the lot costs a few
    // pinned rows)
    bus_q_check::BusQPlan plan;
    bus_q_check::build_bus_q_plan(ctx.view(), plan);
    for (const auto & entry : plan.buses) {
        if (entry.ctrl_pos.empty() && entry.svc_ids.empty()) ctx.controls->reserve_bus_voltage(entry.bus_solver);
        for (int c : entry.ctrl_pos) ctx.controls->reserve_controller_hold(c);
    }
}

bool TransformerVoltageControlLoop::_is_needed(const OuterContext & ctx) const
{
    for (const Group & g : groups(ctx.view())) {
        if (!g.hidden) return true;
    }
    return false;
}

void TransformerVoltageControlLoop::_initialize(OuterContext & ctx)
{
    groups_ = groups(ctx.view());
    plan_ = bus_q_check::BusQPlan();
    bus_q_check::build_bus_q_plan(ctx.view(), plan_);
    step_ = Step::INITIAL;
    frozen_.clear();
    const std::size_t nb = static_cast<std::size_t>(ctx.grid->get_trafos().nb());
    enabled_.assign(nb, 0);
    controller_.assign(nb, 0);
    ratios_.assign(nb, Ratio());
    for (const Group & g : groups_) {
        for (int t : g.trafos) controller_[static_cast<std::size_t>(t)] = 1;
    }
    // every transformer voltage control off for the first solve (no control to act through
    // in detection mode)
    if (ctx.controls == nullptr) return;
    for (std::size_t t = 0; t < nb; ++t) {
        if (!controller_[t]) continue;
        RatioTapControl * tap = ctx.controls->ratio_tap(static_cast<int>(t));
        if (tap != nullptr) tap->set_control(false);
    }
}

void TransformerVoltageControlLoop::_set_on(OuterContext & ctx, int t, bool on)
{
    enabled_[static_cast<std::size_t>(t)] = on ? 1 : 0;
    RatioTapControl * tap = ctx.controls->ratio_tap(t);
    if (tap != nullptr) tap->set_control(on);
}

int TransformerVoltageControlLoop::_closest_tap(const OuterContext & ctx, int t, real_type value) const
{
    // PiModelArray.roundR1ToClosestTap: the current position unless another is strictly closer
    const TrafoContainer & trafos = ctx.grid->get_trafos();
    const TapChangers & rtc = trafos.get_tap_changers(false);
    const int ppos = phase_position(ctx, t);
    int best = ctx.controls->ratio_tap(t)->position();
    real_type best_distance = std::abs(value - trafos.ratio_at(t, best, ppos));
    for (int pos = rtc.low_tap(t); pos <= rtc.high_tap(t); ++pos) {
        const real_type distance = std::abs(value - trafos.ratio_at(t, pos, ppos));
        if (distance < best_distance) {
            best = pos;
            best_distance = distance;
        }
    }
    return best;
}

std::vector<TransformerVoltageControlLoop::ClosedBranch> TransformerVoltageControlLoop::_closed_branches(const OuterGridView & view)
{
    std::vector<ClosedBranch> res;
    const LSGrid & grid = view.grid();
    int b1, b2;
    const int nb_lines = grid.get_lines().nb();
    for (int el = 0; el < nb_lines; ++el) {
        if (view.line_closed(el, b1, b2)) res.push_back(ClosedBranch{b1, b2, -1});
    }
    const int nb_trafos = grid.get_trafos().nb();
    for (int el = 0; el < nb_trafos; ++el) {
        if (view.trafo_closed(el, b1, b2)) res.push_back(ClosedBranch{b1, b2, el});
    }
    return res;
}

bool TransformerVoltageControlLoop::_step_up(const bus_q_check::BusQEntry & entry, real_type limit,
                                             const Eigen::Ref<const RealVect> & vn_kv,
                                             const std::vector<char> & on_controller,
                                             const std::vector<real_type> & min_connected) const
{
    // GeneratorVoltageControlManager.hasStepUpTransformers
    if (!entry.station_ids.empty()) return true;  // a VSC converter station
    if (on_controller[static_cast<std::size_t>(entry.bus_grid)]) return false;
    // NaN (no closed branch at the bus) compares false: not a step-up unit
    const real_type min_v = min_connected[static_cast<std::size_t>(entry.bus_solver)];
    return min_v > vn_kv(entry.bus_grid) && min_v > limit;
}

namespace {

// a controller bus whose voltage control a loop took away: switched PQ, or suspended
bool voltage_control_off(const OuterControls & controls, int bus)
{
    const BusVoltageControl * voltage = controls.bus_voltage(bus);
    return (voltage != nullptr && voltage->is_pq()) || controls.suspended(bus);
}

}  // namespace

void TransformerVoltageControlLoop::_freeze_generators(OuterContext & ctx, real_type limit,
                                                       const std::vector<ClosedBranch> & closed)
{
    // disableGeneratorVoltageControlsUnderMaxControlledNominalVoltage: a controller bus frozen at
    // its bus' Q equation, the injection, which OpenLoadFlow takes as the generation (its load
    // is then counted twice)
    OuterInjections & st = *ctx.injections;
    const LSGrid & grid = *ctx.grid;
    const real_type sn = grid.get_sn_mva();
    const Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    const OuterGridView view = ctx.view();
    const VoltageControlSolverData & ctrl = view.layout().voltage_control.controllers();
    const GlobalBusIdVect & id_solver_to_me = view.id_solver_to_me();
    // what hasStepUpTransformers reads (see _step_up), once for every bus: the grid buses a
    // controller transformer ends at, and per solver bus the lowest of max(vn side 1, vn side 2)
    // over the closed branches at it, NaN when none
    std::vector<char> on_controller(static_cast<std::size_t>(vn_kv.size()), 0);
    {
        const TrafoContainer & trafos = grid.get_trafos();
        const int nb_grid_bus = static_cast<int>(vn_kv.size());
        for (int t = 0; t < trafos.nb(); ++t) {
            if (!controller_[static_cast<std::size_t>(t)]) continue;
            for (const int b : {trafos.get_bus_side_1(t).cast_int(), trafos.get_bus_side_2(t).cast_int()}) {
                if (b >= 0 && b < nb_grid_bus) on_controller[static_cast<std::size_t>(b)] = 1;
            }
        }
    }
    std::vector<real_type> min_connected(id_solver_to_me.size(), std::numeric_limits<real_type>::quiet_NaN());
    for (const ClosedBranch & br : closed) {
        const real_type v = std::max(vn_kv(id_solver_to_me[br.bus1].cast_int()), vn_kv(id_solver_to_me[br.bus2].cast_int()));
        for (const int b : {br.bus1, br.bus2}) {
            real_type & m = min_connected[static_cast<std::size_t>(b)];
            m = std::isnan(m) ? v : std::min(m, v);
        }
    }
    // the load of every bus, pu
    std::vector<real_type> load_q(static_cast<std::size_t>(ctx.V->size()), 0.);
    {
        const int nb_loads = grid.get_loads().nb();
        for (int l = 0; l < nb_loads; ++l) {
            const int b = view.load_bus_solver(l);  // -1 when off
            if (b >= 0 && b < static_cast<int>(load_q.size())) load_q[static_cast<std::size_t>(b)] += view.load_target_q(l) / sn;
        }
    }
    for (std::size_t k = 0; k < plan_.buses.size(); ++k) {
        const auto & entry = plan_.buses[k];
        const bool local = entry.ctrl_pos.empty() && entry.svc_ids.empty();
        int controlled_grid = entry.bus_grid;
        if (!local) {
            const int reg = ctrl.reg_bus(ctrl.group(entry.ctrl_pos.front()));
            controlled_grid = id_solver_to_me[reg].cast_int();
        }
        if (controlled_grid < 0 || vn_kv(controlled_grid) > limit) continue;
        const int b = entry.bus_solver;
        // enabled: not switched PQ by another loop
        if (voltage_control_off(*ctx.controls, b)) continue;
        if (!local) {
            bool held = false;
            for (int c : entry.ctrl_pos) {
                if (ctx.controls->is_held(c)) held = true;
            }
            if (held) continue;
        }
        if (_step_up(entry, limit, vn_kv, on_controller, min_connected)) continue;
        real_type target_vm = ctx.vm(b);
        if (!entry.gen_ids.empty()) target_vm = view.gen_target_vm(entry.gen_ids.front());
        else if (!entry.storage_ids.empty()) target_vm = view.storage_target_vm(entry.storage_ids.front());
        const real_type load = load_q[static_cast<std::size_t>(b)];
        if (local) {
            // the regulating units' output now (they are out of Sbus): the bus' residual
            const real_type q_gen = std::imag((*ctx.bus_mismatch)(b));
            CplxVect & Sbus = *st.Sbus;
            Sbus(b) = cplx_type(std::real(Sbus(b)), std::imag((*st.Sbus_target)(b)) + q_gen - load);
            BusVoltageControl * voltage = ctx.controls->bus_voltage(b);
            if (voltage != nullptr) voltage->set_pq();
        } else {
            bool first = true;
            for (int c : entry.ctrl_pos) {
                real_type q = ctx.controller_q != nullptr && c < ctx.controller_q->size() ? (*ctx.controller_q)(c) : 0.;
                if (first) q -= load;
                first = false;
                VoltageControllerHold * hold = ctx.controls->controller_hold(c);
                if (hold != nullptr) hold->hold(q);
            }
        }
        ctx.controls->set_suspended(b, true);
        frozen_.push_back(Frozen{static_cast<int>(k), b, local, target_vm});
    }
}

void TransformerVoltageControlLoop::_release_generators(OuterContext & ctx)
{
    OuterInjections & st = *ctx.injections;
    for (const Frozen & f : frozen_) {
        if (f.local) {
            CplxVect & Sbus = *st.Sbus;
            Sbus(f.bus_solver) = cplx_type(std::real(Sbus(f.bus_solver)), std::imag((*st.Sbus_target)(f.bus_solver)));
            BusVoltageControl * voltage = ctx.controls->bus_voltage(f.bus_solver);
            if (voltage != nullptr) voltage->set_pv();
            ctx.controls->reset_vm(f.bus_solver, f.target_vm);
        } else {
            for (int c : plan_.buses[static_cast<std::size_t>(f.entry)].ctrl_pos) {
                VoltageControllerHold * hold = ctx.controls->controller_hold(c);
                if (hold != nullptr) hold->release();
            }
        }
        ctx.controls->set_suspended(f.bus_solver, false);
    }
    frozen_.clear();
}

void TransformerVoltageControlLoop::_fix_controls(OuterContext & ctx, const std::vector<ClosedBranch> & closed)
{
    // LfNetwork.fixTransformerVoltageControls: a transformer whose other side, once every
    // transformer switched on is taken out, keeps no PV bus is switched off
    const OuterGridView view = ctx.view();
    std::vector<int> parent(static_cast<std::size_t>(ctx.V->size()));
    std::iota(parent.begin(), parent.end(), 0);
    auto join = [&](int a, int b) {
        if (a < 0 || b < 0 || a >= static_cast<int>(parent.size()) || b >= static_cast<int>(parent.size())) return;
        parent[static_cast<std::size_t>(find_root(parent, a))] = find_root(parent, b);
    };
    for (const ClosedBranch & br : closed) {
        if (br.trafo >= 0 && enabled_[static_cast<std::size_t>(br.trafo)]) continue;
        join(br.bus1, br.bus2);
    }
    // the components holding a PV bus (a generator regulating its own bus, enabled)
    std::set<int> with_pv;
    for (const auto & entry : plan_.buses) {
        if (!(entry.ctrl_pos.empty() && entry.svc_ids.empty())) continue;
        const int b = entry.bus_solver;
        if (voltage_control_off(*ctx.controls, b)) continue;
        with_pv.insert(find_root(parent, b));
    }
    for (const Group & g : groups_) {
        for (int t : g.trafos) {
            if (!enabled_[static_cast<std::size_t>(t)]) continue;
            int b1, b2;
            if (!view.trafo_closed(t, b1, b2)) continue;
            int other;
            if (g.bus_solver == b1) other = b2;
            else if (g.bus_solver == b2) other = b1;
            else continue;
            if (!with_pv.count(find_root(parent, other))) _set_on(ctx, t, false);
        }
    }
}

OuterLoopStatus TransformerVoltageControlLoop::_check(OuterContext & ctx)
{
    if (step_ == Step::COMPLETE || groups_.empty()) return OuterLoopStatus::STABLE;
    // nothing to do without a transformer the inner algorithm can act on
    bool any = false;
    for (const Group & g : groups_) {
        for (int t : g.trafos) any = any || ctx.controls->ratio_tap(t) != nullptr;
    }
    if (!any) return OuterLoopStatus::STABLE;
    const TrafoContainer & trafos = ctx.grid->get_trafos();
    const TapChangers & ptc = trafos.get_tap_changers(true);
    const TapChangers & rtc = trafos.get_tap_changers(false);

    if (step_ == Step::INITIAL) {
        bool need_run = false;
        for (const Group & g : groups_) {
            if (g.hidden) continue;
            const real_type v = ctx.vm(g.bus_solver);
            // the distance to the target against the half deadband, pu
            const bool outside = !(std::abs(g.target - v) <= g.half_deadband);
            ctx.record_bus("SWITCH_ON", outside, g.bus_solver, LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND,
                           std::abs(g.target - v), g.half_deadband);
            if (!outside) continue;
            for (int t : g.trafos) {
                if (handled_ratio(ctx, t) == nullptr) continue;
                _set_on(ctx, t, true);
                need_run = true;
            }
        }
        // TransformerRatioManager: the ratios of the transformers switched on
        for (const Group & g : groups_) {
            std::vector<int> on;
            for (int t : g.trafos) {
                if (enabled_[static_cast<std::size_t>(t)]) on.push_back(t);
            }
            if (on.empty()) continue;
            real_type a = 0., b = 0., sum_min = 0., sum_max = 0.;
            for (int t : on) {
                Ratio & r = ratios_[static_cast<std::size_t>(t)];
                const int ppos = ptc.has(t) ? phase_position(ctx, t) : 0;
                r.initial = ctx.controls->ratio_tap(t)->ratio();
                r.min = std::numeric_limits<real_type>::infinity();
                r.max = -std::numeric_limits<real_type>::infinity();
                for (int pos = rtc.low_tap(t); pos <= rtc.high_tap(t); ++pos) {
                    const real_type value = trafos.ratio_at(t, pos, ppos);
                    r.min = std::min(r.min, value);
                    r.max = std::max(r.max, value);
                }
                if (r.max != r.min) {
                    a += (r.max - r.initial) / (r.max - r.min);
                    b += (r.initial - r.min) / (r.max - r.min);
                } else {
                    a += 1.;
                    b += 1.;
                }
                sum_min += r.min;
                sum_max += r.max;
            }
            const real_type n = static_cast<real_type>(on.size());
            for (int t : on) {
                Ratio & r = ratios_[static_cast<std::size_t>(t)];
                if (params_.use_initial_tap_position) {
                    r.shared_min = sum_min / n;
                    r.shared_max = sum_max / n;
                    r.shared_initial = (r.shared_min * a + r.shared_max * b) / n;
                } else {
                    r.shared_min = r.min;
                    r.shared_max = r.max;
                    r.shared_initial = r.initial;
                }
            }
        }
        if (!need_run) {
            step_ = Step::COMPLETE;
            return OuterLoopStatus::STABLE;
        }
        const std::vector<ClosedBranch> closed = _closed_branches(ctx.view());
        _freeze_generators(ctx, _limit(ctx.view()), closed);
        _fix_controls(ctx, closed);
        step_ = Step::CONTROL;
        return OuterLoopStatus::UNSTABLE;
    }

    // CONTROL: the ratios out of their range rounded to their extreme tap
    bool out_of_range = false;
    for (const Group & g : groups_) {
        for (int t : g.trafos) {
            if (!enabled_[static_cast<std::size_t>(t)]) continue;
            const Ratio & r = ratios_[static_cast<std::size_t>(t)];
            RatioTapControl * tap = ctx.controls->ratio_tap(t);
            const real_type value = tap->ratio();
            if (value < r.shared_min || value > r.shared_max) {
                // the ratio, and the end of the range it left
                ctx.record("ROUND_TO_RANGE", true, ViolationElementType::TRAFO, t, LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND,
                           value, value > r.shared_max ? r.shared_max : r.shared_min);
                tap->move_tap(_closest_tap(ctx, t, value > r.shared_max ? r.shared_max : r.shared_min));
                _set_on(ctx, t, false);
                out_of_range = true;
            }
        }
    }
    if (!out_of_range) {
        // updateContinuousRatio, then every transformer rounded and switched off
        for (const Group & g : groups_) {
            if (g.hidden) continue;
            for (int t : g.trafos) {
                RatioTapControl * tap = handled_ratio(ctx, t);
                if (tap == nullptr) continue;
                real_type value = tap->ratio();
                if (enabled_[static_cast<std::size_t>(t)] && params_.use_initial_tap_position) {
                    const Ratio & r = ratios_[static_cast<std::size_t>(t)];
                    value = value >= r.shared_initial
                        ? r.initial + (value - r.shared_initial) * (r.max - r.initial) / (r.shared_max - r.shared_initial)
                        : r.initial - (r.shared_initial - value) * (r.initial - r.min) / (r.shared_initial - r.shared_min);
                }
                const int position = _closest_tap(ctx, t, value);
                tap->move_tap(position);
                // the continuous ratio, and the tap position it is rounded to
                ctx.record("ROUND_TAP", true, ViolationElementType::TRAFO, t, LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND,
                           value, static_cast<real_type>(position));
                _set_on(ctx, t, false);
            }
        }
        _release_generators(ctx);
        step_ = Step::COMPLETE;
    }
    // in any case, the loop must run again
    return OuterLoopStatus::UNSTABLE;
}

void TransformerVoltageControlLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    if (ctx.V == nullptr || !ctx.vm_checks) return;
    const Eigen::Ref<const RealVect> vn_kv = ctx.grid->get_bus_vn_kv();
    const SolverBusIdVect & to_solver = solver_map(ctx);
    for (const Group & g : groups(ctx.view())) {
        if (g.hidden) continue;
        const int b = to_solver[g.bus_grid].cast_int();
        if (b < 0 || b >= ctx.V->size()) continue;
        const real_type v = ctx.vm(b);
        if (std::abs(g.target - v) > std::max(g.half_deadband, ctx.tol_vm_pu)) {
            const real_type vn = vn_kv(g.bus_grid);
            out.push_back(LimitViolation{ViolationElementType::BUS, g.bus_grid, 0,
                                         LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND, v * vn, g.target * vn, ""});
        }
    }
}

AlgoConfig TransformerVoltageControlLoop::_get_params() const
{
    AlgoConfig cfg;
    cfg.int_params = {params_.use_initial_tap_position ? 1 : 0};
    cfg.real_params = {static_cast<double>(params_.max_controlled_nominal_voltage),
                       static_cast<double>(params_.min_target_deadband_kv)};
    return cfg;
}

}  // namespace ls2g
