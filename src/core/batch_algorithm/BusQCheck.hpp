// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BUSQCHECK_H
#define BUSQCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"
#include "powerflow_algorithm/outer_loop/OuterRowView.hpp"
#include "powerflow_algorithm/outer_loop/BaseOuterLoop.hpp"
#include "BusQPlan.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace ls2g {

/**
 * The reactive capability of a bus, PER BUS: which machines hold each bus' voltage and where
 * their reactive output is read from. The plan `ReactiveLimitsLoop` (and
 * `TransformerVoltageControlLoop`) builds its controller buses from, on the grid -- or the batch
 * row -- an OuterGridView shows; that loop's trigger is what `compute_physical_violations`
 * reports (see outer_detect).
 *
 * WHAT THIS ANSWERS. A voltage-regulating machine has no reactive setpoint: its Q is
 * whatever it takes to hold the set-point, and lightsim2grid never clamps it. So a row
 * may converge holding a bus at 1.05 pu with its machines having to produce three times
 * the reactive power they own. That is not a limit somebody chose to exceed -- it is a
 * solution the grid cannot reach at all, because a machine cannot produce reactive power
 * it does not have. See ViolationCategory::PHYSICAL. It is the condition PowSyBl
 * OpenLoadFlow's `ReactiveLimits` outer loop acts on (it switches the bus PV -> PQ and
 * re-solves).
 *
 * WHAT HOLDS A BUS' VOLTAGE, AND WHAT IT CAN DO. Exactly four families have a reactive
 * output the solver computes rather than reads, and all four are checked -- the sum of
 * their capability is the bus' capability:
 *
 *   - a voltage-regulating GENERATOR: [min_q_mvar, max_q_mvar], fixed;
 *   - a voltage-regulating STORAGE UNIT: the same [min_q_mvar, max_q_mvar] in generator
 *     convention (what it can inject, see StorageContainer::init_full), fixed. It only
 *     ever pins its own bus through the PV path, never as a controller of the plan;
 *   - a voltage-regulating hvdc CONVERTER STATION: its own [min_q, max_q], also in MVAr
 *     and also fixed (`get_station_min_q_mvar` / `get_station_max_q_mvar`);
 *   - a voltage-mode SVC: a SUSCEPTANCE range `[b_min, b_max]` (pu, base sn_mva), so its
 *     capability depends on the solved voltage -- `q = b . |V|^2` in generator convention
 *     (a shunt admittance `jb` consumes `-j.b.|V|^2`, ie injects `+b.|V|^2`), hence
 *     `[b_min, b_max] . |V|^2 . sn_mva` MVAr, re-evaluated every row.
 *
 * Anything else standing on the bus has a reactive injection that is INPUT data, part of
 * `Sbus` (a load, a shunt, a non-regulating storage unit or generator, a REACTIVE_POWER
 * mode SVC, a fixed-Q station): it is not what holds the voltage, and a violation of its
 * own limits would be an input error rather than something a solve produced. Those are
 * not in the picture at all.
 *
 * WHY PER BUS AND NOT PER MACHINE. The reactive power a bus needs is a fact about the
 * solution; how it is divided between several machines standing on that bus is not -- the
 * solver never decides it, `LSGrid::_split_q_residual_per_bus` does, by a sharing
 * convention (proportional to each machine's reactive range). Checking machine by machine
 * would therefore report the convention rather than the grid: two 20 MVAr machines
 * covering 30 MVAr together is fine, and a per-machine check on the same solution says
 * both are 5 MVAr over. So the question asked here is the one that has an answer
 * independent of the convention: did this bus need more reactive power than the sum of
 * what the machines holding it can produce?
 *
 * WHY IT LIVES HERE AND NOT ON LSGrid. A single solve publishes every element's reactive
 * output (`GenInfo::res_q_mvar` and its SVC / station counterparts, see
 * LSGrid::compute_results), so the check is a few lines in Python. A batch publishes no
 * such thing: it keeps the voltages and drops everything else, precisely because
 * reconstructing the per-element results of ten thousand rows is what the batch classes
 * exist not to do. So for a batch the bus' reactive power has to be re-derived, row by
 * row, which is what this does.
 *
 * HOW THE BUS' REACTIVE POWER IS RE-DERIVED. From the algorithm's own per-bus mismatch
 * (`BaseAlgo::get_bus_mismatch`, see LSGrid::_fill_bus_mismatch_ac for why that is the
 * right thing to read and not merely the cheap one) plus the reactive injections it
 * solved for explicitly:
 *
 *     q_bus = imag(mis_bus(b)) + sum of the Q of the controllers standing on b
 *
 * `mis_bus` is the raw residual `V .* conj(Ybus . V) - Sbus` with each component's own
 * injection added back -- including `- i * Q_c` at every `VoltageControl` controller's
 * bus. Adding those `Q_c` back (`BaseAlgo::get_controller_q`) therefore restores the raw
 * reactive residual at `b`, and the raw residual IS what the machines holding that bus
 * produced, since a regulating machine's Q is never part of `Sbus`. One expression covers
 * an ordinary PV bus (no controller, so just the mismatch), a bus a control group holds
 * (every machine on it is a controller, so just their Q), and the mixture.
 *
 * WHICH MACHINES COUNT, exactly. One rule, and it is the complement of the above: a
 * machine's capability enters the sum iff its reactive output is part of `q_bus`, ie iff
 * it is a controller of the plan (its Q comes back through `get_controller_q`) or it pins
 * its own bus through the classical PV path (its Q is inside the residual). A machine that
 * is neither -- a remote-regulating generator absent from the controller list, which is
 * what an algorithm that cannot do remote voltage control leaves behind -- is left out
 * rather than counted against a bus that never saw its reactive power.
 */
namespace bus_q_check {

/**
 * Work out which buses can be checked at all, which machines hold each one, and where its
 * reactive power will have to be read from -- in the grid as `view` sees it: a batch row's
 * own labelling (its layout) and its own machines (a generator it disconnects or places
 * elsewhere, see OuterGridView).
 *
 * The controller list is the layout's voltage-control plan, whose order is the order of
 * `BaseAlgo::get_controller_q()`. Entries come in grid bus order.
 */
inline void build_bus_q_plan(const OuterGridView & view, BusQPlan & out)
{
    out.clear();
    const LSGrid & grid_model = view.grid();
    const SolverBusIdVect & id_me_to_solver = view.id_me_to_solver();
    const VoltageControlSolverData & ctrl = view.layout().voltage_control.controllers();
    const int nb_bus_grid = static_cast<int>(grid_model.total_bus());
    if(nb_bus_grid == 0) return;
    const int nb_bus_solver = static_cast<int>(view.id_solver_to_me().size());

    // Whether each element is a controller of the plan: the only way its reactive output
    // reaches `q_bus` other than through the residual. Keyed the way the controller list
    // is (kind, element id); for an hvdc line the kind says which end.
    const GeneratorContainer & generators = grid_model.get_generators();
    const SvcContainer & svcs = grid_model.get_svcs();
    const HvdcLineContainer & hvdcs = grid_model.get_dclines();
    const int nb_gen = generators.nb();
    const int nb_svc = svcs.nb();
    const int nb_hvdc = hvdcs.nb();
    const int nb_ctrl = ctrl.n_controllers();

    std::vector<char> gen_is_ctrl(static_cast<std::size_t>(nb_gen > 0 ? nb_gen : 0), 0);
    std::vector<char> svc_is_ctrl(static_cast<std::size_t>(nb_svc > 0 ? nb_svc : 0), 0);
    std::vector<char> st1_is_ctrl(static_cast<std::size_t>(nb_hvdc > 0 ? nb_hvdc : 0), 0);
    std::vector<char> st2_is_ctrl(static_cast<std::size_t>(nb_hvdc > 0 ? nb_hvdc : 0), 0);
    for(int c = 0; c < nb_ctrl; ++c){
        // a held controller (LSGrid::set_hold_frozen_regulators) holds no bus: its
        // frozen output is its bus' Sbus, exactly as without the option
        if(ctrl.is_held(c)) continue;
        const int el_id = ctrl.elem_id(c);
        switch(ctrl.kind(c)){
            case VoltageControlSolverData::GEN:
                if(el_id >= 0 && el_id < nb_gen) gen_is_ctrl[static_cast<std::size_t>(el_id)] = 1;
                break;
            case VoltageControlSolverData::SVC:
                if(el_id >= 0 && el_id < nb_svc) svc_is_ctrl[static_cast<std::size_t>(el_id)] = 1;
                break;
            case VoltageControlSolverData::HVDC_SIDE_1:
                if(el_id >= 0 && el_id < nb_hvdc) st1_is_ctrl[static_cast<std::size_t>(el_id)] = 1;
                break;
            default:  // HVDC_SIDE_2
                if(el_id >= 0 && el_id < nb_hvdc) st2_is_ctrl[static_cast<std::size_t>(el_id)] = 1;
                break;
        }
    }

    // per SOLVER bus, who holds it -- see "WHICH MACHINES COUNT" above: a controller of the
    // plan, or an element pinning its own bus through the classical PV path
    auto in_range = [nb_bus_solver](int b){ return b >= 0 && b < nb_bus_solver; };
    auto solver_of = [&id_me_to_solver, nb_bus_grid](int bus_me){
        if(bus_me < 0 || bus_me >= nb_bus_grid) return -1;
        const int b = id_me_to_solver[bus_me].cast_int();
        return b == BaseConstants::_deactivated_bus_id ? -1 : b;
    };
    const std::size_t nb = static_cast<std::size_t>(nb_bus_solver > 0 ? nb_bus_solver : 0);
    std::vector<std::vector<int> > gens_of_bus(nb);
    std::vector<std::vector<std::pair<int, int> > > stations_of_bus(nb);
    std::vector<std::vector<int> > svcs_of_bus(nb);
    std::vector<std::vector<int> > storages_of_bus(nb);

    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        const bool is_ctrl = gen_is_ctrl[static_cast<std::size_t>(gen_id)] != 0;
        // a local voltage controller: connected in this row, regulating, not treated as off,
        // AND pinning its own bus -- the residual path
        if(is_ctrl ? !view.gen_on(gen_id) : !view.gen_is_local_voltage_controller(gen_id)) continue;
        const int b = view.gen_bus_solver(gen_id);
        if(!in_range(b)) continue;
        gens_of_bus[static_cast<std::size_t>(b)].push_back(gen_id);
    }
    // a storage unit only regulates its own bus (LSGrid::check_grid refuses anything else),
    // so it is never a controller of the plan: the PV path is the whole rule for it
    const int nb_storage = grid_model.get_storages().nb();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!view.storage_is_local_voltage_controller(storage_id)) continue;
        const int b = view.storage_bus_solver(storage_id);
        if(!in_range(b)) continue;
        storages_of_bus[static_cast<std::size_t>(b)].push_back(storage_id);
    }
    // a station always regulates the bus it stands on, so it is either a group controller
    // or an ordinary PV pin -- `station_is_voltage_controller` covers both
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            const bool is_ctrl = (side == 1 ? st1_is_ctrl : st2_is_ctrl)[static_cast<std::size_t>(hvdc_id)] != 0;
            if(!is_ctrl && !hvdcs.station_is_voltage_controller(hvdc_id, side)) continue;
            const int b = solver_of(hvdcs.get_station_bus(hvdc_id, side).cast_int());
            if(!in_range(b)) continue;
            stations_of_bus[static_cast<std::size_t>(b)].push_back(std::make_pair(hvdc_id, side));
        }
    }
    // an SVC never takes the PV path (see SvcContainer): it holds a bus only as a
    // controller, so the controller list is the whole rule for it
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(svc_is_ctrl[static_cast<std::size_t>(svc_id)] == 0) continue;
        const int b = solver_of(svcs.get_bus_id()(svc_id).cast_int());
        if(!in_range(b)) continue;
        svcs_of_bus[static_cast<std::size_t>(b)].push_back(svc_id);
    }

    const SubstationContainer & subs = grid_model.get_substations();
    const std::vector<std::string> & sub_names = subs.get_sub_names();  // empty if never set

    for(int bus_me = 0; bus_me < nb_bus_grid; ++bus_me){
        const int bus_solver = solver_of(bus_me);
        if(!in_range(bus_solver)) continue;  // not in the solved system
        // stranded by the row: its mismatch means nothing, and nothing holds its voltage
        if(view.is_masked(bus_solver)) continue;
        const std::size_t b = static_cast<std::size_t>(bus_solver);
        if(gens_of_bus[b].empty() && storages_of_bus[b].empty() && stations_of_bus[b].empty() &&
           svcs_of_bus[b].empty()) continue;

        BusQEntry entry;
        entry.bus_solver = bus_solver;
        entry.bus_grid = bus_me;
        entry.gen_ids = gens_of_bus[b];
        entry.storage_ids = storages_of_bus[b];
        entry.station_ids = stations_of_bus[b];
        entry.svc_ids = svcs_of_bus[b];
        for(int c = 0; c < nb_ctrl; ++c){
            if(ctrl.bus(c) != bus_solver || ctrl.is_held(c)) continue;
            entry.ctrl_pos.push_back(c);
        }
        if(!entry.ctrl_pos.empty()) out.needs_controller_q = true;
        const int sub_id = subs.sub_id_of_bus(bus_me);
        if(sub_id >= 0 && static_cast<std::size_t>(sub_id) < sub_names.size()){
            entry.sub_name = sub_names[static_cast<std::size_t>(sub_id)];
        }
        out.buses.push_back(entry);
    }
}

/// build_bus_q_plan of ctx.view(): the solve's own when the driver keeps a cache (built by
/// the first loop asking, then shared by the others), built into `own` otherwise
inline const BusQPlan & solve_bus_q_plan(const OuterContext & ctx, BusQPlan & own)
{
    if(ctx.solve_cache == nullptr){
        build_bus_q_plan(ctx.view(), own);
        return own;
    }
    OuterSolveCache & cache = *ctx.solve_cache;
    if(!cache.has_bus_q_plan){
        build_bus_q_plan(ctx.view(), cache.bus_q_plan);
        cache.has_bus_q_plan = true;
    }
    return cache.bus_q_plan;
}

}  // namespace bus_q_check
}  // namespace ls2g

#endif  // BUSQCHECK_H
