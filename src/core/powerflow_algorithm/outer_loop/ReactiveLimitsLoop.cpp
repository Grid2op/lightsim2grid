// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ReactiveLimitsLoop.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "LSGrid.hpp"
#include "OuterLoopDriver.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.view().id_me_to_solver();
}

}  // namespace

ReactiveLimitsLoop::ReactiveLimitsLoop() : ReactiveLimitsLoop(Params()) {}

ReactiveLimitsLoop::ReactiveLimitsLoop(const Params & params) : params_(params)
{
    if (params.max_pq_pv_switch < 0) throw std::runtime_error("ReactiveLimitsLoop: max_pq_pv_switch must be >= 0.");
    if (!(params.max_reactive_power_mismatch >= 0.)) {
        throw std::runtime_error("ReactiveLimitsLoop: max_reactive_power_mismatch must be >= 0.");
    }
    if (!(params.mismatch_base_mva > 0.)) throw std::runtime_error("ReactiveLimitsLoop: mismatch_base_mva must be > 0.");
    if (params.robust_mode < -1 || params.robust_mode > 1) {
        throw std::runtime_error("ReactiveLimitsLoop: robust_mode must be 1 (on), 0 (off) or -1 (the driver's).");
    }
    // NaN: the driver's
    if (std::isinf(params.min_realistic_voltage) || std::isinf(params.max_realistic_voltage)) {
        throw std::runtime_error("ReactiveLimitsLoop: min_realistic_voltage / max_realistic_voltage must be finite "
                                 "(or NaN for the driver's).");
    }
    if (!std::isnan(params.min_realistic_voltage) && !std::isnan(params.max_realistic_voltage) &&
        !(params.min_realistic_voltage < params.max_realistic_voltage)) {
        throw std::runtime_error("ReactiveLimitsLoop: min_realistic_voltage must be lower than max_realistic_voltage.");
    }
    if (!(params.realistic_voltage_margin > 0.)) {
        throw std::runtime_error("ReactiveLimitsLoop: realistic_voltage_margin must be > 0.");
    }
    if (!(params.robust_restart_vm_pu > 0.)) {
        throw std::runtime_error("ReactiveLimitsLoop: robust_restart_vm_pu must be > 0.");
    }
}

void ReactiveLimitsLoop::_declare(const OuterContext & ctx) const
{
    // every bus holding its own voltage through the PV path may become PQ
    OuterControls & controls = *ctx.controls;
    bus_q_check::BusQPlan own;
    const bus_q_check::BusQPlan & plan = bus_q_check::solve_bus_q_plan(ctx, own);
    for (const auto & entry : plan.buses) {
        if (entry.ctrl_pos.empty() && entry.svc_ids.empty()) controls.reserve_bus_voltage(entry.bus_solver);
        // ... and the controllers of a group may be held at a limit
        for (int c : entry.ctrl_pos) controls.reserve_controller_hold(c);
    }
    // so may a voltage monitor, once switched on (see _initialize)
    const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
    for (int c = 0; c < ctrl.n_controllers(); ++c) {
        if (ctrl.is_held(c) && ctrl.kind(c) == VoltageControlSolverData::SVC) controls.reserve_controller_hold(c);
    }
}

bool ReactiveLimitsLoop::_has_monitor(const OuterContext & ctx)
{
    const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
    for (int c = 0; c < ctrl.n_controllers(); ++c) {
        if (ctrl.is_held(c) && ctrl.kind(c) == VoltageControlSolverData::SVC) return true;
    }
    return false;
}

bool ReactiveLimitsLoop::_is_needed(const OuterContext & ctx) const
{
    bus_q_check::BusQPlan own;
    return !bus_q_check::solve_bus_q_plan(ctx, own).empty() || _has_monitor(ctx);
}

void ReactiveLimitsLoop::_initialize(OuterContext & ctx)
{
    {
        const bus_q_check::BusQPlan & plan = bus_q_check::solve_bus_q_plan(ctx, plan_);
        if (&plan != &plan_) plan_ = plan;
    }
    // a voltage monitor is a held controller of the plan, which bus_q_check leaves out: once
    // switched on it holds its bus like any SVC, so it gets a bus of its own here
    std::vector<int> monitor_of_entry(plan_.buses.size(), -1);
    {
        const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
        const GlobalBusIdVect & id_solver_to_me = ctx.view().id_solver_to_me();
        for (int c = 0; c < ctrl.n_controllers(); ++c) {
            if (!ctrl.is_held(c) || ctrl.kind(c) != VoltageControlSolverData::SVC) continue;
            bus_q_check::BusQEntry entry;
            entry.bus_solver = ctrl.bus(c);
            entry.bus_grid = id_solver_to_me[ctrl.bus(c)].cast_int();
            entry.svc_ids.push_back(ctrl.elem_id(c));
            entry.ctrl_pos.push_back(c);
            plan_.buses.push_back(entry);
            monitor_of_entry.push_back(ctrl.elem_id(c));
        }
    }
    // resuming: the buses the previous solve left frozen, by grid bus (0: none; of two at one
    // bus, the first), frozen again below for this solve. Not by solver bus: a topology change
    // since that solve gives the buses other solver ids. A voltage monitor by its SVC: it
    // shares its grid bus with the entry of that bus' own units, whose state is not its own
    std::vector<int> carried_state;
    std::vector<char> carried_realistic;
    std::vector<int> carried_monitor_state;
    std::vector<char> carried_monitor_realistic;
    if (ctx.resume) {
        const std::size_t nb_grid_bus = static_cast<std::size_t>(ctx.grid->total_bus());
        const std::size_t nb_svc = static_cast<std::size_t>(ctx.grid->get_svcs().nb());
        for (const ControllerBus & bus : buses_) {
            if (bus.state == 0) continue;
            if (bus.monitor_svc >= 0) {
                const std::size_t s = static_cast<std::size_t>(bus.monitor_svc);
                if (s >= nb_svc) continue;
                if (carried_monitor_state.empty()) {
                    carried_monitor_state.assign(nb_svc, 0);
                    carried_monitor_realistic.assign(nb_svc, 0);
                }
                carried_monitor_state[s] = bus.state;
                carried_monitor_realistic[s] = bus.realistic ? 1 : 0;
                continue;
            }
            if (bus.bus_grid < 0 || static_cast<std::size_t>(bus.bus_grid) >= nb_grid_bus) continue;
            if (carried_state.empty()) {
                carried_state.assign(nb_grid_bus, 0);
                carried_realistic.assign(nb_grid_bus, 0);
            }
            const std::size_t b = static_cast<std::size_t>(bus.bus_grid);
            if (carried_state[b] != 0) continue;
            carried_state[b] = bus.state;
            carried_realistic[b] = bus.realistic ? 1 : 0;
        }
    }
    buses_.clear();
    const LSGrid & grid = *ctx.grid;
    const OuterGridView view = ctx.view();
    const GeneratorContainer & gens = grid.get_generators();
    const StorageContainer & storages = grid.get_storages();
    const HvdcLineContainer & hvdcs = grid.get_dclines();
    const Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    // the robust mode's: the target Q of every generator, summed per solver bus once (each bus
    // in generator order)
    std::vector<real_type> target_q_of_bus(view.id_solver_to_me().size(), 0.);
    for (int gen_id = 0; gen_id < gens.nb(); ++gen_id) {
        const int b = view.gen_bus_solver(gen_id);  // -1 when off
        if (b >= 0 && static_cast<std::size_t>(b) < target_q_of_bus.size()) {
            target_q_of_bus[static_cast<std::size_t>(b)] += gens.get_target_q()(gen_id);
        }
    }
    for (std::size_t k = 0; k < plan_.buses.size(); ++k) {
        const auto & entry = plan_.buses[k];
        ControllerBus bus;
        bus.entry = static_cast<int>(k);
        bus.bus_solver = entry.bus_solver;
        bus.bus_grid = entry.bus_grid;
        bus.local = entry.ctrl_pos.empty() && entry.svc_ids.empty();
        if (bus.local && ctx.controls != nullptr) bus.voltage = ctx.controls->bus_voltage(entry.bus_solver);
        bus.monitor_svc = monitor_of_entry[k];
        bus.reg_bus_solver = entry.bus_solver;
        bus.nominal_v = vn_kv(entry.bus_grid);
        bool has_target = false;
        if (!bus.local) {
            // the group's: the bus it regulates, its set-point
            const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
            const int g = ctrl.group(entry.ctrl_pos.front());
            bus.group = g;
            bus.reg_bus_solver = ctrl.reg_bus(g);
            bus.target_vm = ctrl.v_set(g);
            has_target = true;
            const int reg_grid = ctx.view().id_solver_to_me()[bus.reg_bus_solver].cast_int();
            if (reg_grid >= 0) bus.nominal_v = vn_kv(reg_grid);
        }
        bus.target_q = target_q_of_bus[static_cast<std::size_t>(entry.bus_solver)];
        // (a generator's limits follow its target P: read at each check, see _limits). The
        // target in force: one a resumed solve's DistributedSlack moved (see
        // OuterInjections::gen_target_p_in_force) ranks the bus (KEPT_PV_STRONGEST)
        for (int gen_id : entry.gen_ids) {
            bus.target_p += ctx.injections != nullptr ? ctx.injections->gen_target_p_in_force(view, gen_id)
                                                      : view.gen_target_p(gen_id);
            if (!has_target) { bus.target_vm = view.gen_target_vm(gen_id); has_target = true; }
        }
        for (int storage_id : entry.storage_ids) {
            bus.min_q += storages.get_min_q(storage_id);
            bus.max_q += storages.get_max_q(storage_id);
            bus.target_p -= ctx.injections != nullptr ? ctx.injections->storage_target_p_in_force(view, storage_id)
                                                      : view.storage_target_p(storage_id);  // load convention
            if (!has_target) { bus.target_vm = view.storage_target_vm(storage_id); has_target = true; }
        }
        for (const auto & station : entry.station_ids) {
            bus.min_q += hvdcs.get_station_min_q_mvar(station.first, station.second);
            bus.max_q += hvdcs.get_station_max_q_mvar(station.first, station.second);
            if (!has_target) {
                bus.target_vm = hvdcs.get_station_target_vm_pu(station.first, station.second);
                has_target = true;
            }
        }
        buses_.push_back(bus);
    }
    // ... each at the limit it holds in this solve (its target P may have moved), or at its
    // target Q for the robust mode's; the release test lets it go when this solve calls for it.
    // A monitor only once switched on (VoltageMonitoring, ahead in the list, switches back on
    // the ones the previous solve did): an idle one is held at Q = 0, and nothing would release it
    if ((!carried_state.empty() || !carried_monitor_state.empty()) && ctx.injections != nullptr &&
        ctx.controls != nullptr) {
        for (ControllerBus & bus : buses_) {
            int state = 0;
            bool realistic = false;
            if (bus.monitor_svc >= 0) {
                const std::size_t s = static_cast<std::size_t>(bus.monitor_svc);
                if (s >= carried_monitor_state.size()) continue;
                const StandbySvcControl * standby = ctx.controls->standby_svc(bus.monitor_svc);
                if (standby == nullptr || !standby->is_on()) continue;
                state = carried_monitor_state[s];
                realistic = carried_monitor_realistic[s] != 0;
            } else {
                if (bus.bus_grid < 0 || static_cast<std::size_t>(bus.bus_grid) >= carried_state.size()) continue;
                state = carried_state[static_cast<std::size_t>(bus.bus_grid)];
                realistic = carried_realistic[static_cast<std::size_t>(bus.bus_grid)] != 0;
            }
            if (state == 0) continue;
            bus.realistic = realistic;
            real_type q_min, q_max;
            _limits(ctx, bus, q_min, q_max);
            _freeze(ctx, bus, bus.realistic ? bus.target_q : (state < 0 ? q_min : q_max), state);
        }
    }
}

real_type ReactiveLimitsLoop::_bus_q(const OuterContext & ctx, const ControllerBus & bus) const
{
    // bus_q_check's: the bus' residual, plus what its controllers solved for
    const auto & entry = plan_.buses[static_cast<std::size_t>(bus.entry)];
    real_type q = std::imag((*ctx.bus_mismatch)(bus.bus_solver));
    if (ctx.controller_q != nullptr) {
        for (int pos : entry.ctrl_pos) {
            if (pos < ctx.controller_q->size()) q += (*ctx.controller_q)(pos);
        }
    }
    return q * ctx.grid->get_sn_mva();
}

void ReactiveLimitsLoop::_limits(const OuterContext & ctx, const ControllerBus & bus,
                                 real_type & q_min, real_type & q_max) const
{
    q_min = bus.min_q;
    q_max = bus.max_q;
    const auto & entry = plan_.buses[static_cast<std::size_t>(bus.entry)];
    // a generator: its curve at the target P an outer loop gave it, its fixed limits otherwise
    for (int gen_id : entry.gen_ids) {
        real_type lo, hi;
        ctx.grid->gen_limits_at_outer_target(gen_id, _gen_target_p(ctx, gen_id), lo, hi);
        q_min += lo;
        q_max += hi;
    }
    if (entry.svc_ids.empty() || ctx.V == nullptr) return;
    // an SVC: its susceptance range at the bus' voltage (generator convention)
    const SvcContainer & svcs = ctx.grid->get_svcs();
    const real_type v2_sn = std::norm((*ctx.V)(bus.bus_solver)) * ctx.grid->get_sn_mva();
    for (int svc_id : entry.svc_ids) {
        q_min += svcs.get_b_min(svc_id) * v2_sn;
        q_max += svcs.get_b_max(svc_id) * v2_sn;
    }
}

real_type ReactiveLimitsLoop::_gen_target_p(const OuterContext & ctx, int gen_id)
{
    static const std::vector<real_type> none;
    return ctx.view().gen_target_p_override(ctx.injections != nullptr ? ctx.injections->gen_target_p : none, gen_id);
}

real_type ReactiveLimitsLoop::_gen_limit_slope(const OuterContext & ctx, int gen_id, bool max)
{
    real_type dlo, dhi;
    ctx.grid->gen_limit_slopes_at_outer_target(gen_id, _gen_target_p(ctx, gen_id), dlo, dhi);
    return max ? dhi : dlo;
}

void ReactiveLimitsLoop::_adjoint_terms(OuterAdjointTerms & out) const
{
    for (const ControllerBus & bus : buses_) {
        if (bus.state != 0) out.q_slopes.insert(out.q_slopes.end(), bus.slopes.begin(), bus.slopes.end());
    }
}

real_type ReactiveLimitsLoop::_controller_limit(const OuterContext & ctx, int ctrl_pos, bool max) const
{
    const LSGrid & grid = *ctx.grid;
    const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
    const int el = ctrl.elem_id(ctrl_pos);
    switch (ctrl.kind(ctrl_pos)) {
        case VoltageControlSolverData::GEN: {
            real_type lo, hi;
            grid.gen_limits_at_outer_target(el, _gen_target_p(ctx, el), lo, hi);
            return max ? hi : lo;
        }
        case VoltageControlSolverData::SVC: {
            const real_type v2_sn = std::norm((*ctx.V)(ctrl.bus(ctrl_pos))) * grid.get_sn_mva();
            return (max ? grid.get_svcs().get_b_max(el) : grid.get_svcs().get_b_min(el)) * v2_sn;
        }
        case VoltageControlSolverData::HVDC_SIDE_1:
        case VoltageControlSolverData::HVDC_SIDE_2: {
            const int side = ctrl.kind(ctrl_pos) == VoltageControlSolverData::HVDC_SIDE_1 ? 1 : 2;
            return max ? grid.get_dclines().get_station_max_q_mvar(el, side)
                       : grid.get_dclines().get_station_min_q_mvar(el, side);
        }
        default:
            return 0.;
    }
}

std::vector<char> ReactiveLimitsLoop::_groups_holding(const OuterContext & ctx) const
{
    const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
    std::vector<char> holding(static_cast<std::size_t>(ctrl.n_groups()), 0);
    std::vector<char> sloped(holding.size(), 0);
    for (int c = 0; c < ctrl.n_controllers(); ++c) {
        const std::size_t g = static_cast<std::size_t>(ctrl.group(c));
        if (ctrl.slope(c) != 0.) sloped[g] = 1;
        // held at a limit by this loop
        if (ctx.controls != nullptr && ctx.controls->is_held(c)) continue;
        // held by the plan: a frozen regulator, or a voltage monitor not switched on
        if (ctrl.is_held(c)) {
            const StandbySvcControl * standby = ctx.controls != nullptr ? ctx.controls->standby_svc(ctrl.elem_id(c)) : nullptr;
            const bool switched_on = ctrl.kind(c) == VoltageControlSolverData::SVC &&
                                     standby != nullptr && standby->is_on();
            if (!switched_on) continue;
        }
        holding[g] = 1;
    }
    // with a slope, the regulated voltage is off its set-point by design
    for (std::size_t g = 0; g < holding.size(); ++g) holding[g] = holding[g] && !sloped[g];
    return holding;
}

void ReactiveLimitsLoop::_evaluate(const OuterContext & ctx, std::vector<Switch> & to_pq,
                                   std::vector<Switch> & to_pv, std::vector<int> & moved,
                                   int & remaining_pv, std::vector<Switch> * kept) const
{
    // OpenLoadFlow's tolerance, plus a detection's own (0 in the outer-loop mode)
    const real_type eps = params_.max_reactive_power_mismatch * params_.mismatch_base_mva + ctx.tol_mw;
    real_type v_low, v_high;
    _realistic_band(ctx, v_low, v_high);
    v_low -= ctx.tol_vm_pu;
    v_high += ctx.tol_vm_pu;
    const bool robust = _robust(ctx);
    remaining_pv = 0;
    const std::vector<char> groups_holding = _groups_holding(ctx);
    for (std::size_t k = 0; k < buses_.size(); ++k) {
        const ControllerBus & bus = buses_[k];
        const int ki = static_cast<int>(k);
        // a bus whose voltage control another loop suspended (TransformerVoltageControl): neither
        // a PV bus to check nor one this loop switched, as OpenLoadFlow's (a frozen bus with no
        // reactive limit type)
        if (ctx.controls != nullptr && ctx.controls->suspended(bus.bus_solver)) continue;
        real_type target_vm = bus.target_vm;
        if (bus.monitor_svc >= 0) {
            const StandbySvcControl * standby = ctx.controls != nullptr ? ctx.controls->standby_svc(bus.monitor_svc) : nullptr;
            if (standby == nullptr || !standby->is_on()) continue;  // still idle
            target_vm = standby->target_vm();  // the set-point it was switched on at
        }
        real_type q_min, q_max;
        _limits(ctx, bus, q_min, q_max);
        if (bus.state == 0) {
            // PV: checkControllerBus
            const real_type q = _bus_q(ctx, bus);
            if (q < q_min - eps) {
                to_pq.push_back(Switch{ki, LimitViolationType::LOW_Q, q, q_min});
            } else if (q > q_max + eps) {
                to_pq.push_back(Switch{ki, LimitViolationType::HIGH_Q, q, q_max});
            } else if (robust && bus.reg_bus_solver != bus.bus_solver) {
                // a remote controller within its limits, but its own bus' voltage unrealistic
                const real_type v = ctx.vm(bus.bus_solver);
                if (v < v_low) {
                    to_pq.push_back(Switch{ki, LimitViolationType::LOW_Q, q, bus.target_q, true});
                } else if (v > v_high) {
                    to_pq.push_back(Switch{ki, LimitViolationType::HIGH_Q, q, bus.target_q, true});
                } else {
                    ++remaining_pv;
                }
            } else {
                ++remaining_pv;
            }
            continue;
        }
        // PQ: checkPqBus, the voltage it regulates against its set-point. A bus another
        // controller of its group still holds is AT the set-point: the Newton's voltage row
        // makes it so, to within a rounding that differs from one CPU to the next, and a
        // strict comparison of the two would release the bus on the sign of that rounding.
        // Taken as equal, the bus stays frozen, which is what OpenLoadFlow's strict test
        // says whenever the two are equal.
        const bool held_by_group = bus.group >= 0 && groups_holding[static_cast<std::size_t>(bus.group)];
        const real_type vm_read = ctx.vm(bus.reg_bus_solver);
        const real_type vm = held_by_group ? target_vm : vm_read;
        const real_type vn = bus.nominal_v;
        const LimitViolationType release = bus.state < 0 ? LimitViolationType::LOW_VOLTAGE_AT_MIN_Q
                                                         : LimitViolationType::HIGH_VOLTAGE_AT_MAX_Q;
        if (bus.state < 0 ? vm < target_vm : vm > target_vm) {
            to_pv.push_back(Switch{ki, release, vm * vn, target_vm * vn});
            continue;
        }
        // not released: the test as it read (the magnitude the solve gave, even where the
        // group's regulation settled it), for the decision trace
        if (kept != nullptr) {
            Switch sw{ki, release, vm_read * vn, target_vm * vn};
            sw.group_holds = held_by_group;
            kept->push_back(sw);
        }
        const real_type q_lim = bus.state < 0 ? q_min : q_max;
        if (!bus.realistic && std::abs(q_lim - bus.frozen_q) > eps) moved.push_back(ki);
    }
}

void ReactiveLimitsLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    // the switches check() would make, from the state the loop is in: in detection mode the
    // one initialize() set, every controller bus PV -- so the PV -> PQ tests of OpenLoadFlow's
    // first check
    if (!ctx.vm_checks || ctx.V == nullptr || ctx.bus_mismatch == nullptr) return;
    std::vector<Switch> to_pq, to_pv;
    std::vector<int> moved;
    int remaining_pv = 0;
    _evaluate(ctx, to_pq, to_pv, moved, remaining_pv);
    const OuterGridView view = ctx.view();
    for (const std::vector<Switch> * list : {&to_pq, &to_pv}) {
        for (const Switch & sw : *list) {
            const ControllerBus & bus = buses_[static_cast<std::size_t>(sw.k)];
            const auto & entry = plan_.buses[static_cast<std::size_t>(bus.entry)];
            if (view.is_masked(bus.bus_solver)) continue;
            if (!sw.realistic) {
                out.push_back(LimitViolation{ViolationElementType::BUS, entry.bus_grid, 0, sw.type,
                                             sw.value, sw.limit, entry.sub_name});
                continue;
            }
            // the robust mode's: a remote controller whose own bus is at an unrealistic
            // voltage, on each generator holding that bus (kV of the bus)
            const bool low = sw.type == LimitViolationType::LOW_Q;
            const real_type vn = ctx.grid->get_bus_vn_kv()(entry.bus_grid);
            real_type v_low, v_high;
            _realistic_band(ctx, v_low, v_high);
            const real_type bound = low ? v_low : v_high;
            const std::vector<std::string> & names = ctx.grid->get_generators().get_names();  // empty if never set
            for (int gen_id : entry.gen_ids) {
                const std::string name = static_cast<std::size_t>(gen_id) < names.size() ? names[static_cast<std::size_t>(gen_id)] : std::string();
                out.push_back(LimitViolation{ViolationElementType::GENERATOR, gen_id, 0,
                                             low ? LimitViolationType::LOW_VOLTAGE_REMOTE_CONTROL
                                                 : LimitViolationType::HIGH_VOLTAGE_REMOTE_CONTROL,
                                             ctx.vm(bus.bus_solver) * vn, bound * vn, name});
            }
        }
    }
}

void ReactiveLimitsLoop::_freeze(OuterContext & ctx, ControllerBus & bus, real_type q_mvar, int state) const
{
    // its units now inject their limit, on top of what the bus' Sbus held already (their
    // own Q is not in it while they regulate)
    OuterInjections & st = *ctx.injections;
    const real_type sn = ctx.grid->get_sn_mva();
    const auto & entry = plan_.buses[static_cast<std::size_t>(bus.entry)];
    // a generator frozen at a limit of its curve follows it: what the adjoint adds to its
    // gen_p gradient (the robust mode's target Q does not move with P)
    bus.slopes.clear();
    if (!bus.local) {
        // each of its controllers held at its own limit: the bus' at that limit
        const VoltageControlSolverData & ctrl = ctx.view().layout().voltage_control.controllers();
        for (int c : entry.ctrl_pos) {
            real_type q = _controller_limit(ctx, c, state > 0);
            if (ctrl.kind(c) == VoltageControlSolverData::GEN) {
                if (bus.realistic) {
                    q = ctx.grid->get_generators().get_target_q()(ctrl.elem_id(c));  // the robust mode's
                } else {
                    const real_type slope = _gen_limit_slope(ctx, ctrl.elem_id(c), state > 0);
                    if (slope != 0.) bus.slopes.push_back(OuterQSlope{ctrl.elem_id(c), bus.bus_solver, c, slope / sn});
                }
            }
            VoltageControllerHold * hold = ctx.controls->controller_hold(c);
            if (hold != nullptr) hold->hold(q / sn);
        }
        bus.state = state;
        bus.frozen_q = q_mvar;
        return;
    }
    if (!bus.realistic) {
        for (int gen_id : entry.gen_ids) {
            const real_type slope = _gen_limit_slope(ctx, gen_id, state > 0);
            if (slope != 0.) bus.slopes.push_back(OuterQSlope{gen_id, bus.bus_solver, -1, slope / sn});
        }
    }
    CplxVect & Sbus = *st.Sbus;
    const real_type q_init = std::imag((*st.Sbus_target)(bus.bus_solver));
    Sbus(bus.bus_solver) = cplx_type(std::real(Sbus(bus.bus_solver)), q_init + q_mvar / sn);
    if (bus.voltage != nullptr) bus.voltage->set_pq();
    bus.state = state;
    bus.frozen_q = q_mvar;
}

void ReactiveLimitsLoop::_release(OuterContext & ctx, ControllerBus & bus) const
{
    OuterInjections & st = *ctx.injections;
    if (!bus.local) {
        const auto & entry = plan_.buses[static_cast<std::size_t>(bus.entry)];
        for (int c : entry.ctrl_pos) {
            VoltageControllerHold * hold = ctx.controls->controller_hold(c);
            if (hold != nullptr) hold->release();
        }
        bus.state = 0;
        bus.realistic = false;
        bus.frozen_q = 0.;
        bus.slopes.clear();
        return;
    }
    CplxVect & Sbus = *st.Sbus;
    Sbus(bus.bus_solver) = cplx_type(std::real(Sbus(bus.bus_solver)), std::imag((*st.Sbus_target)(bus.bus_solver)));
    if (bus.voltage != nullptr) bus.voltage->set_pv();
    // a pinned row keeps the magnitude the bus has: back at its set-point
    ctx.controls->reset_vm(bus.bus_solver, bus.target_vm);
    bus.state = 0;
    bus.realistic = false;
    bus.frozen_q = 0.;
    bus.slopes.clear();
}

OuterLoopStatus ReactiveLimitsLoop::_check(OuterContext & ctx)
{
    std::vector<Switch> to_pq, to_pv, kept;
    std::vector<int> moved;
    int remaining_pv = 0;
    _evaluate(ctx, to_pq, to_pv, moved, remaining_pv, &kept);
    bool changed = false;
    // the decision trace: a controller bus by its own bus
    auto record = [&](const char * action, bool taken, const Switch & sw) {
        ctx.record_bus(action, taken, buses_[static_cast<std::size_t>(sw.k)].bus_solver, sw.type, sw.value, sw.limit);
    };
    for (const Switch & sw : kept) record(sw.group_holds ? "KEPT_PQ_GROUP_HOLDS" : "KEPT_PQ", false, sw);

    // PV -> PQ, keeping the strongest bus PV when every one of them would switch
    if (!to_pq.empty() && remaining_pv == 0) {
        auto stronger = [&](const Switch & a, const Switch & b) {
            const ControllerBus & ba = buses_[static_cast<std::size_t>(a.k)];
            const ControllerBus & bb = buses_[static_cast<std::size_t>(b.k)];
            if (ba.nominal_v != bb.nominal_v) return ba.nominal_v > bb.nominal_v;
            if (ba.target_p != bb.target_p) return ba.target_p > bb.target_p;
            return ba.bus_solver < bb.bus_solver;
        };
        const auto strongest = std::min_element(to_pq.begin(), to_pq.end(), stronger);
        record("KEPT_PV_STRONGEST", false, *strongest);
        to_pq.erase(strongest);
    }
    for (const Switch & sw : to_pq) {
        record(sw.realistic ? "PV_TO_PQ_UNREALISTIC" : "PV_TO_PQ", true, sw);
        ControllerBus & bus = buses_[static_cast<std::size_t>(sw.k)];
        bus.realistic = sw.realistic;
        _freeze(ctx, bus, sw.limit, sw.type == LimitViolationType::LOW_Q ? -1 : +1);
        ++bus.nb_pv_pq;
        // the robust mode: a remote controller with an unrealistic voltage of its own
        // restarts from 1 pu
        if (_robust(ctx) && bus.reg_bus_solver != bus.bus_solver) {
            const real_type v = ctx.vm(bus.bus_solver);
            real_type v_low, v_high;
            _realistic_band(ctx, v_low, v_high);
            if (sw.realistic || v < v_low || v > v_high) {
                ctx.controls->reset_vm(bus.bus_solver, params_.robust_restart_vm_pu);
            }
        }
        changed = true;
    }
    // PQ -> PV, but not past max_pq_pv_switch
    for (const Switch & sw : to_pv) {
        ControllerBus & bus = buses_[static_cast<std::size_t>(sw.k)];
        if (bus.nb_pv_pq >= params_.max_pq_pv_switch) {
            record("KEPT_PQ_MAX_SWITCH", false, sw);
            continue;
        }
        record("PQ_TO_PV", true, sw);
        _release(ctx, bus);
        changed = true;
    }
    // a frozen bus whose limit moved: frozen at the new one
    for (int k : moved) {
        ControllerBus & bus = buses_[static_cast<std::size_t>(k)];
        real_type q_min, q_max;
        _limits(ctx, bus, q_min, q_max);
        ctx.record_bus("LIMIT_MOVED", true, bus.bus_solver, bus.state < 0 ? LimitViolationType::LOW_Q : LimitViolationType::HIGH_Q,
                       bus.state < 0 ? q_min : q_max, bus.frozen_q);
        _freeze(ctx, bus, bus.state < 0 ? q_min : q_max, bus.state);
        changed = true;
    }
    return changed ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

bool ReactiveLimitsLoop::_robust(const OuterContext & ctx) const
{
    return params_.robust_mode < 0 ? ctx.driver_params().voltage_remote_control_robust_mode : params_.robust_mode != 0;
}

void ReactiveLimitsLoop::_realistic_band(const OuterContext & ctx, real_type & v_low, real_type & v_high) const
{
    const OuterLoopDriverParams & driver = ctx.driver_params();
    const real_type v_min = std::isnan(params_.min_realistic_voltage) ? driver.min_realistic_voltage : params_.min_realistic_voltage;
    const real_type v_max = std::isnan(params_.max_realistic_voltage) ? driver.max_realistic_voltage : params_.max_realistic_voltage;
    v_low = v_min * params_.realistic_voltage_margin;
    v_high = v_max / params_.realistic_voltage_margin;
}

AlgoConfig ReactiveLimitsLoop::_get_params() const
{
    AlgoConfig res;
    res.int_params = {params_.max_pq_pv_switch, params_.robust_mode};
    res.real_params = {params_.max_reactive_power_mismatch, params_.mismatch_base_mva,
                       params_.min_realistic_voltage, params_.max_realistic_voltage,
                       params_.realistic_voltage_margin, params_.robust_restart_vm_pu};
    return res;
}

}  // namespace ls2g
