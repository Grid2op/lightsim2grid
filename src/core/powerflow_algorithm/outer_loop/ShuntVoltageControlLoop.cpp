// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ShuntVoltageControlLoop.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

#include "LSGrid.hpp"
#include "powerflow_algorithm/NRSystem.hpp"

namespace ls2g {

namespace {

const SolverBusIdVect & solver_map(const OuterContext & ctx)
{
    return ctx.id_me_to_solver != nullptr ? *ctx.id_me_to_solver : ctx.grid->id_me_to_ac_solver();
}

}  // namespace

std::vector<ShuntVoltageControlLoop::Group> ShuntVoltageControlLoop::groups(const LSGrid & grid,
                                                                             const std::set<int> & transformer_buses)
{
    std::vector<Group> res;
    const ShuntContainer & shunts = grid.get_shunts();
    const SolverBusIdVect & to_solver = grid.id_me_to_ac_solver();
    if (to_solver.size() == 0) return res;
    // the buses a control of a higher priority holds: a generator, an SVC, a VSC station
    std::set<int> higher(transformer_buses);
    for (const SolverBusId & b : grid.get_ac_pv_solver()) higher.insert(b.cast_int());
    for (const SolverBusId & b : grid.get_slack_ids_solver()) higher.insert(b.cast_int());
    const VoltageControlSolverData & ctrl = grid.get_ac_voltage_control_plan().controllers();
    for (int g = 0; g < ctrl.n_groups(); ++g) higher.insert(ctrl.reg_bus(g));

    // the controllers: the regulating shunts of each bus, in shunt order
    std::vector<int> controller_order;
    std::map<int, std::vector<int> > shunts_of_bus;
    std::map<int, int> regulated_of_bus;  // the bus its first shunt regulates (grid id)
    for (int s = 0; s < shunts.nb(); ++s) {
        if (!shunts.has_sections(s) || !shunts.get_section_regulating(s) || !shunts.get_status(s)) continue;
        const int reg = shunts.get_section_regulated_bus(s);
        if (reg < 0 || reg >= static_cast<int>(to_solver.size()) || to_solver[reg].cast_int() < 0) continue;
        const int bus = to_solver[shunts.get_bus_id()(s).cast_int()].cast_int();
        if (bus < 0) continue;
        if (!shunts_of_bus.count(bus)) {
            controller_order.push_back(bus);
            regulated_of_bus[bus] = reg;
        }
        shunts_of_bus[bus].push_back(s);
    }
    // the groups, by the bus regulated (a controller in one group only: its first shunt's)
    std::map<int, std::size_t> group_of_bus;
    for (int bus : controller_order) {
        const int reg = regulated_of_bus[bus];
        const int reg_solver = to_solver[reg].cast_int();
        auto it = group_of_bus.find(reg_solver);
        const int first = shunts_of_bus[bus].front();
        if (it == group_of_bus.end()) {
            Group g;
            g.bus_solver = reg_solver;
            g.bus_grid = reg;
            g.target = shunts.get_section_target_vm_pu(first);
            g.half_deadband = 0.;
            g.hidden = higher.count(reg_solver) > 0;
            it = group_of_bus.emplace(reg_solver, res.size()).first;
            res.push_back(g);
        }
        Group & g = res[it->second];
        g.controller_buses.push_back(bus);
        g.shunts.push_back(shunts_of_bus[bus]);
        // the smallest deadband
        for (int s : shunts_of_bus[bus]) {
            const real_type deadband = shunts.get_section_deadband_pu(s);
            if (deadband > 0.) g.half_deadband = g.half_deadband > 0. ? std::min(g.half_deadband, deadband / 2.) : deadband / 2.;
        }
    }
    return res;
}

void ShuntVoltageControlLoop::_declare(const OuterContext & ctx, OuterDeclaration & decl) const
{
    // the buses the transformers regulate: a transformer voltage control declared before this
    // loop (OpenLoadFlow's order) hides a shunt one
    transformer_buses_.clear();
    for (const auto & g : decl.ratio_groups()) transformer_buses_.insert(g.bus_solver);
    // a hidden group never acts (OpenLoadFlow's getControllerElements keeps the visible
    // controls only): nothing to reserve for it
    for (const Group & g : groups(*ctx.grid, transformer_buses_)) {
        if (!g.hidden) decl.add_shunt_group(g.bus_solver, g.target, g.controller_buses, g.shunts, true);
    }
}

bool ShuntVoltageControlLoop::_is_needed(const OuterContext & ctx) const
{
    for (const Group & g : groups(*ctx.grid, transformer_buses_)) {
        if (!g.hidden) return true;
    }
    return false;
}

void ShuntVoltageControlLoop::_initialize(OuterContext & ctx)
{
    groups_ = groups(*ctx.grid, transformer_buses_);
    OuterState & st = *ctx.state;
    st.shunt_control.assign(static_cast<std::size_t>(ctx.grid->id_ac_solver_to_me().size()), -1);
    // every controller on for the first solve
    for (const Group & g : groups_) {
        if (g.hidden) continue;
        for (int bus : g.controller_buses) st.shunt_control[static_cast<std::size_t>(bus)] = 1;
    }
}

std::vector<int> ShuntVoltageControlLoop::_dispatch(const LSGrid & grid, const std::vector<int> & ids, real_type b) const
{
    const ShuntContainer & shunts = grid.get_shunts();
    const real_type sn = grid.get_sn_mva();
    auto b_at = [&](int s, int count) { return -shunts.section_q(s, count) / sn; };
    // the largest first (the B between its fewest sections and its most), a stable order
    std::vector<std::size_t> order(ids.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<real_type> magnitude(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const int s = ids[i];
        magnitude[i] = std::abs(b_at(s, shunts.get_max_section_count(s)) - b_at(s, shunts.get_min_section_count(s)));
    }
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) { return magnitude[a] > magnitude[c]; });
    std::vector<int> counts(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i) counts[i] = shunts.get_section_count(ids[i]);
    real_type residue = b;
    std::size_t remaining = ids.size();
    for (std::size_t i : order) {
        const int s = ids[i];
        const real_type share = residue / static_cast<real_type>(remaining--);
        // roundBToClosestSection: the current count unless another is strictly closer
        int best = counts[i];
        real_type best_distance = std::abs(share - b_at(s, best));
        for (int count = shunts.get_min_section_count(s); count <= shunts.get_max_section_count(s); ++count) {
            const real_type distance = std::abs(share - b_at(s, count));
            if (distance < best_distance) {
                best = count;
                best_distance = distance;
            }
        }
        counts[i] = best;
        residue -= b_at(s, best);
    }
    return counts;
}

OuterLoopStatus ShuntVoltageControlLoop::_check(OuterContext & ctx)
{
    if (ctx.iteration != 0 || groups_.empty() || ctx.shunt_control == nullptr) return OuterLoopStatus::STABLE;
    OuterState & st = *ctx.state;
    const ShuntControl & shunt = *ctx.shunt_control;
    bool any = false;
    for (const Group & g : groups_) {
        if (g.hidden) continue;
        for (std::size_t c = 0; c < g.controller_buses.size(); ++c) {
            const int bus = g.controller_buses[c];
            if (!shunt.handles(bus)) continue;
            st.shunt_control[static_cast<std::size_t>(bus)] = 0;
            st.shunt_sections.emplace_back(bus, _dispatch(*ctx.grid, g.shunts[c], shunt.b(bus)));
            any = true;
        }
    }
    return any ? OuterLoopStatus::UNSTABLE : OuterLoopStatus::STABLE;
}

void ShuntVoltageControlLoop::_detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const
{
    if (ctx.V == nullptr || !ctx.vm_checks) return;
    const Eigen::Ref<const RealVect> vn_kv = ctx.grid->get_bus_vn_kv();
    const SolverBusIdVect & to_solver = solver_map(ctx);
    for (const Group & g : groups(*ctx.grid, transformer_buses_)) {
        if (g.hidden) continue;
        const int b = to_solver[g.bus_grid].cast_int();
        if (b < 0 || b >= ctx.V->size()) continue;
        const real_type v = std::abs((*ctx.V)(b));
        if (std::abs(g.target - v) > std::max(g.half_deadband, ctx.tol_vm_pu)) {
            const real_type vn = vn_kv(g.bus_grid);
            out.push_back(LimitViolation{ViolationElementType::BUS, g.bus_grid, 0,
                                         LimitViolationType::SHUNT_VOLTAGE_CONTROL, v * vn, g.target * vn, ""});
        }
    }
}

}  // namespace ls2g
