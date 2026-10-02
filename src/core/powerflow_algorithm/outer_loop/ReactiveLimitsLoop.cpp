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
#include <sstream>
#include <stdexcept>

#include "LSGrid.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver : ctx.grid->id_me_to_ac_solver();
}

}  // namespace

bus_q_check::BusQPlan ReactiveLimitsLoop::_plan(const OuterContext & ctx)
{
    bus_q_check::BusQPlan plan;
    bus_q_check::build_bus_q_plan(*ctx.grid, solver_map(ctx),
                                  ctx.grid->get_ac_voltage_control_plan().controllers(), plan);
    return plan;
}

void ReactiveLimitsLoop::_declare(const OuterContext & ctx, OuterDeclaration & decl) const
{
    // every bus holding its own voltage through the PV path may become PQ
    const bus_q_check::BusQPlan plan = _plan(ctx);
    for (const auto & entry : plan.buses) {
        if (entry.ctrl_pos.empty() && entry.svc_ids.empty()) decl.add_switchable_vm_bus(entry.bus_solver);
    }
}

bool ReactiveLimitsLoop::_is_needed(const OuterContext & ctx) const
{
    return !_plan(ctx).empty();
}

void ReactiveLimitsLoop::_initialize(OuterContext & ctx)
{
    plan_ = _plan(ctx);
    buses_.clear();
    const LSGrid & grid = *ctx.grid;
    const GeneratorContainer & gens = grid.get_generators();
    const StorageContainer & storages = grid.get_storages();
    const HvdcLineContainer & hvdcs = grid.get_dclines();
    const Eigen::Ref<const RealVect> vn_kv = grid.get_bus_vn_kv();
    for (std::size_t k = 0; k < plan_.buses.size(); ++k) {
        const auto & entry = plan_.buses[k];
        ControllerBus bus;
        bus.entry = static_cast<int>(k);
        bus.bus_solver = entry.bus_solver;
        bus.local = entry.ctrl_pos.empty() && entry.svc_ids.empty();
        bus.nominal_v = vn_kv(entry.bus_grid);
        bool has_target = false;
        for (int gen_id : entry.gen_ids) {
            bus.min_q += gens.get_min_q(gen_id);
            bus.max_q += gens.get_max_q(gen_id);
            bus.target_p += gens.get_target_p()(gen_id);
            if (!has_target) { bus.target_vm = gens.get_target_vm_pu(gen_id); has_target = true; }
        }
        for (int storage_id : entry.storage_ids) {
            bus.min_q += storages.get_min_q(storage_id);
            bus.max_q += storages.get_max_q(storage_id);
            bus.target_p -= storages.get_target_p()(storage_id);  // load convention
            if (!has_target) { bus.target_vm = storages.get_target_vm_pu(storage_id); has_target = true; }
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

void ReactiveLimitsLoop::_evaluate(const OuterContext & ctx, std::vector<Switch> & to_pq,
                                   std::vector<Switch> & to_pv, std::vector<int> & moved,
                                   int & remaining_pv) const
{
    const real_type eps = max_reactive_power_mismatch * OLF_SB_MVA;
    remaining_pv = 0;
    for (std::size_t k = 0; k < buses_.size(); ++k) {
        const ControllerBus & bus = buses_[k];
        const int ki = static_cast<int>(k);
        if (bus.state == 0) {
            // PV: checkControllerBus (only the local ones switch for now)
            if (!bus.local) { ++remaining_pv; continue; }
            const real_type q = _bus_q(ctx, bus);
            if (q < bus.min_q - eps) {
                to_pq.push_back(Switch{ki, LimitViolationType::LOW_Q, q, bus.min_q});
            } else if (q > bus.max_q + eps) {
                to_pq.push_back(Switch{ki, LimitViolationType::HIGH_Q, q, bus.max_q});
            } else {
                ++remaining_pv;
            }
            continue;
        }
        // PQ: checkPqBus, the voltage it regulates against its set-point
        const real_type vm = std::abs((*ctx.V)(bus.bus_solver));
        const real_type vn = bus.nominal_v;
        if (bus.state < 0) {
            if (vm < bus.target_vm) {
                to_pv.push_back(Switch{ki, LimitViolationType::LOW_VOLTAGE_AT_MIN_Q, vm * vn, bus.target_vm * vn});
            } else if (std::abs(bus.min_q - bus.frozen_q) > eps) {
                moved.push_back(ki);
            }
        } else {
            if (vm > bus.target_vm) {
                to_pv.push_back(Switch{ki, LimitViolationType::HIGH_VOLTAGE_AT_MAX_Q, vm * vn, bus.target_vm * vn});
            } else if (std::abs(bus.max_q - bus.frozen_q) > eps) {
                moved.push_back(ki);
            }
        }
    }
}

void ReactiveLimitsLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    // detection mode: the physical checks (bus_q_check, gen_pv_release_check) already
    // report these on the grid's own solve
    if (ctx.is_detection() || ctx.V == nullptr || ctx.bus_mismatch == nullptr) return;
    std::vector<Switch> to_pq, to_pv;
    std::vector<int> moved;
    int remaining_pv = 0;
    _evaluate(ctx, to_pq, to_pv, moved, remaining_pv);
    for (const std::vector<Switch> * list : {&to_pq, &to_pv}) {
        for (const Switch & sw : *list) {
            const auto & entry = plan_.buses[static_cast<std::size_t>(buses_[static_cast<std::size_t>(sw.k)].entry)];
            out.push_back(LimitViolation{ViolationElementType::BUS, entry.bus_grid, 0, sw.type,
                                         sw.value, sw.limit, entry.sub_name});
        }
    }
}

void ReactiveLimitsLoop::_freeze(OuterContext & ctx, ControllerBus & bus, real_type q_mvar, int state) const
{
    // its units now inject their limit, on top of what the bus' Sbus held already (their
    // own Q is not in it while they regulate)
    OuterState & st = *ctx.state;
    const real_type sn = ctx.grid->get_sn_mva();
    CplxVect & Sbus = *st.Sbus;
    const real_type q_init = std::imag((*st.Sbus_init)(bus.bus_solver));
    Sbus(bus.bus_solver) = cplx_type(std::real(Sbus(bus.bus_solver)), q_init + q_mvar / sn);
    st.pq_buses.insert(bus.bus_solver);
    bus.state = state;
    bus.frozen_q = q_mvar;
}

void ReactiveLimitsLoop::_release(OuterContext & ctx, ControllerBus & bus) const
{
    OuterState & st = *ctx.state;
    CplxVect & Sbus = *st.Sbus;
    Sbus(bus.bus_solver) = cplx_type(std::real(Sbus(bus.bus_solver)), std::imag((*st.Sbus_init)(bus.bus_solver)));
    st.pq_buses.erase(bus.bus_solver);
    // a pinned row keeps the magnitude the bus has: back at its set-point
    st.vm_set.push_back(std::make_pair(bus.bus_solver, bus.target_vm));
    bus.state = 0;
    bus.frozen_q = 0.;
}

OuterLoopStatus ReactiveLimitsLoop::_check(OuterContext & ctx)
{
    std::vector<Switch> to_pq, to_pv;
    std::vector<int> moved;
    int remaining_pv = 0;
    _evaluate(ctx, to_pq, to_pv, moved, remaining_pv);
    bool changed = false;

    // PV -> PQ, keeping the strongest bus PV when every one of them would switch
    if (!to_pq.empty() && remaining_pv == 0) {
        auto stronger = [&](const Switch & a, const Switch & b) {
            const ControllerBus & ba = buses_[static_cast<std::size_t>(a.k)];
            const ControllerBus & bb = buses_[static_cast<std::size_t>(b.k)];
            if (ba.nominal_v != bb.nominal_v) return ba.nominal_v > bb.nominal_v;
            if (ba.target_p != bb.target_p) return ba.target_p > bb.target_p;
            return ba.bus_solver < bb.bus_solver;
        };
        to_pq.erase(std::min_element(to_pq.begin(), to_pq.end(), stronger));
    }
    for (const Switch & sw : to_pq) {
        ControllerBus & bus = buses_[static_cast<std::size_t>(sw.k)];
        _freeze(ctx, bus, sw.limit, sw.type == LimitViolationType::LOW_Q ? -1 : +1);
        ++bus.nb_pv_pq;
        changed = true;
    }
    // PQ -> PV, but not past max_pq_pv_switch
    for (const Switch & sw : to_pv) {
        ControllerBus & bus = buses_[static_cast<std::size_t>(sw.k)];
        if (bus.nb_pv_pq >= max_pq_pv_switch) continue;
        _release(ctx, bus);
        changed = true;
    }
    // a frozen bus whose limit moved: frozen at the new one
    for (int k : moved) {
        ControllerBus & bus = buses_[static_cast<std::size_t>(k)];
        _freeze(ctx, bus, bus.state < 0 ? bus.min_q : bus.max_q, bus.state);
        changed = true;
    }
    return changed ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

AlgoConfig ReactiveLimitsLoop::_get_params() const
{
    AlgoConfig res;
    res.int_params = {max_pq_pv_switch};
    res.real_params = {max_reactive_power_mismatch};
    return res;
}

void ReactiveLimitsLoop::_set_params(const AlgoConfig & params)
{
    if (params.int_params.size() != 1 || params.real_params.size() != 1) {
        throw std::runtime_error("ReactiveLimitsLoop::set_params: expects 1 int and 1 real parameter.");
    }
    if (params.int_params[0] < 0) throw std::runtime_error("ReactiveLimitsLoop: max_pq_pv_switch must be >= 0.");
    if (!(params.real_params[0] >= 0.)) throw std::runtime_error("ReactiveLimitsLoop: max_reactive_power_mismatch must be >= 0.");
    max_pq_pv_switch = params.int_params[0];
    max_reactive_power_mismatch = params.real_params[0];
}

}  // namespace ls2g
