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

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace ls2g {

/**
 * Post-solve reactive-capability check, PER BUS -- one half of
 * `compute_physical_violations` (the other is HvdcPCheck.hpp).
 *
 * WHAT THIS ANSWERS. A voltage-regulating machine has no reactive setpoint: its Q is
 * whatever it takes to hold the set-point, and lightsim2grid never clamps it. So a row
 * may converge holding a bus at 1.05 pu with its machines having to produce three times
 * the reactive power they own. That is not a limit somebody chose to exceed -- it is a
 * solution the grid cannot reach at all, because a machine cannot produce reactive power
 * it does not have. See ViolationCategory::PHYSICAL. It is the condition PowSyBl
 * OpenLoadFlow's `ReactiveLimits` outer loop acts on (it would switch the bus PV -> PQ
 * and re-solve); nothing here re-solves anything, it only reports.
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

/// One bus whose voltage is held by machines, and everything a row needs to check it.
/// Built once per compute() by `build_bus_q_plan`; nothing here varies from row to row.
struct BusQEntry
{
    int bus_solver = -1;          ///< the bus, solver numbering
    int bus_grid = -1;            ///< the same bus, grid numbering (what is reported)
    /// the generators holding this bus: a row that disconnects one takes its limits out of
    /// the sum (`min_q_mvar` / `max_q_mvar`, fixed)
    std::vector<int> gen_ids;
    /// the voltage-regulating storage units holding it (always through the PV path, so
    /// never part of `ctrl_pos`): a fixed [min_q, max_q] in MVAr, generator convention
    std::vector<int> storage_ids;
    /// the hvdc converter stations holding it, as (hvdc line id, side): also a fixed
    /// [min_q, max_q] in MVAr
    std::vector<std::pair<int, int> > station_ids;
    /// the voltage-mode SVCs holding it: a susceptance range, so their MVAr capability is
    /// re-evaluated at every row's own voltage
    std::vector<int> svc_ids;
    /// positions in the algorithm's controller list of the controllers standing on this
    /// bus (`BaseAlgo::get_controller_q`), empty for an ordinary PV bus
    std::vector<int> ctrl_pos;
    std::string sub_name;         ///< substation of the bus, empty if the grid has no names
};

struct BusQPlan
{
    std::vector<BusQEntry> buses;
    /// whether any bus of the plan needs the controller list at all: false on an ordinary
    /// grid of local PV machines, and then a row never asks for it
    bool needs_controller_q = false;

    bool empty() const { return buses.empty(); }
    void clear() { buses.clear(); needs_controller_q = false; }
};

/**
 * Work out, once, which buses can be checked at all, which machines hold each one, and
 * where its reactive power will have to be read from.
 *
 * `id_me_to_solver` and `ctrl` must both describe the labelling the batch solves in
 * (`active_layout()`); `ctrl` is the plan's controller list, whose order is the order of
 * `BaseAlgo::get_controller_q()`.
 */
inline void build_bus_q_plan(const LSGrid & grid_model,
                             const SolverBusIdVect & id_me_to_solver,
                             const VoltageControlSolverData & ctrl,
                             BusQPlan & out)
{
    out.clear();
    const int nb_bus_grid = static_cast<int>(grid_model.total_bus());
    if(nb_bus_grid == 0) return;

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

    // per GRID bus, who holds it -- see "WHICH MACHINES COUNT" above: a controller of the
    // plan, or an element pinning its own bus through the classical PV path
    std::vector<std::vector<int> > gens_of_bus(static_cast<std::size_t>(nb_bus_grid));
    std::vector<std::vector<std::pair<int, int> > > stations_of_bus(static_cast<std::size_t>(nb_bus_grid));
    std::vector<std::vector<int> > svcs_of_bus(static_cast<std::size_t>(nb_bus_grid));

    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        const bool is_ctrl = gen_is_ctrl[static_cast<std::size_t>(gen_id)] != 0;
        // is_local_voltage_controller: active (connected, regulating, not treated as off)
        // AND pinning its own bus -- the residual path
        if(!is_ctrl && !generators.is_local_voltage_controller(gen_id)) continue;
        const int bus_me = gen_buses(gen_id).cast_int();
        if(bus_me < 0 || bus_me >= nb_bus_grid) continue;
        gens_of_bus[static_cast<std::size_t>(bus_me)].push_back(gen_id);
    }
    // a storage unit only regulates its own bus (LSGrid::check_grid refuses anything else),
    // so it is never a controller of the plan: the PV path is the whole rule for it
    const StorageContainer & storages = grid_model.get_storages();
    const int nb_storage = storages.nb();
    std::vector<std::vector<int> > storages_of_bus(static_cast<std::size_t>(nb_bus_grid));
    const GlobalBusIdVect & storage_buses = storages.get_bus_id();
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        if(!storages.is_local_voltage_controller(storage_id)) continue;
        const int bus_me = storage_buses(storage_id).cast_int();
        if(bus_me < 0 || bus_me >= nb_bus_grid) continue;
        storages_of_bus[static_cast<std::size_t>(bus_me)].push_back(storage_id);
    }
    // a station always regulates the bus it stands on, so it is either a group controller
    // or an ordinary PV pin -- `station_is_voltage_controller` covers both
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            const bool is_ctrl = (side == 1 ? st1_is_ctrl : st2_is_ctrl)[static_cast<std::size_t>(hvdc_id)] != 0;
            if(!is_ctrl && !hvdcs.station_is_voltage_controller(hvdc_id, side)) continue;
            const int bus_me = hvdcs.get_station_bus(hvdc_id, side).cast_int();
            if(bus_me < 0 || bus_me >= nb_bus_grid) continue;
            stations_of_bus[static_cast<std::size_t>(bus_me)].push_back(std::make_pair(hvdc_id, side));
        }
    }
    // an SVC never takes the PV path (see SvcContainer): it holds a bus only as a
    // controller, so the controller list is the whole rule for it
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(svc_is_ctrl[static_cast<std::size_t>(svc_id)] == 0) continue;
        const int bus_me = svcs.get_bus_id()(svc_id).cast_int();
        if(bus_me < 0 || bus_me >= nb_bus_grid) continue;
        svcs_of_bus[static_cast<std::size_t>(bus_me)].push_back(svc_id);
    }

    const SubstationContainer & subs = grid_model.get_substations();
    const std::vector<std::string> & sub_names = subs.get_sub_names();  // empty if never set

    for(int bus_me = 0; bus_me < nb_bus_grid; ++bus_me){
        const std::size_t b = static_cast<std::size_t>(bus_me);
        if(gens_of_bus[b].empty() && storages_of_bus[b].empty() && stations_of_bus[b].empty() &&
           svcs_of_bus[b].empty()) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        BusQEntry entry;
        entry.bus_solver = bus_solver;
        entry.bus_grid = bus_me;
        entry.gen_ids = gens_of_bus[b];
        entry.storage_ids = storages_of_bus[b];
        entry.station_ids = stations_of_bus[b];
        entry.svc_ids = svcs_of_bus[b];
        for(int c = 0; c < nb_ctrl; ++c){
            if(ctrl.bus(c) != bus_solver) continue;
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

/**
 * Append to `out` one LimitViolation per bus whose machines had to produce more (or less)
 * reactive power than they own, for ONE converged row.
 *
 * `bus_mismatch` is the algorithm's own per-bus mismatch (solver numbering, pu, see
 * BaseAlgo::get_bus_mismatch), `controller_q` its converged per-controller reactive
 * injection (pu), and `V` the row's converged complex voltage (solver numbering, pu --
 * needed for the SVCs, whose capability is a susceptance). `masked_solver_ids` is this
 * row's masked (stranded) solver buses -- sorted, may be nullptr -- whose voltage is
 * reported as 0 and whose mismatch means nothing. `is_gen_off(gen_id)` tells whether this
 * row disconnected that generator (a generator contingency): it produces nothing, so it
 * leaves the summed capability, and a bus whose every machine is gone is not checked at
 * all (it is an ordinary PQ bus in that row, and the residual there is nobody's reactive
 * output).
 */
template<class IsGenOff>
inline void check_bus_q_violations(const BusQPlan & plan,
                                   const LSGrid & grid_model,
                                   const Eigen::Ref<const CplxVect> & bus_mismatch,
                                   const Eigen::Ref<const CplxVect> & V,
                                   const RealVect & controller_q,
                                   real_type sn_mva,
                                   real_type tol_mvar,
                                   const std::vector<int> * masked_solver_ids,
                                   IsGenOff is_gen_off,
                                   std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;
    if(bus_mismatch.size() == 0) return;

    const GeneratorContainer & generators = grid_model.get_generators();
    const StorageContainer & storages = grid_model.get_storages();
    const SvcContainer & svcs = grid_model.get_svcs();
    const HvdcLineContainer & hvdcs = grid_model.get_dclines();
    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();

    for(std::size_t k = 0; k < plan.buses.size(); ++k){
        const BusQEntry & entry = plan.buses[k];
        if(entry.bus_solver >= static_cast<int>(bus_mismatch.size())) continue;
        // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
        if(has_masked && std::binary_search(masked_solver_ids->begin(), masked_solver_ids->end(),
                                            entry.bus_solver)) continue;

        // ---- what this bus' machines can produce, this row ------------------------
        real_type q_min = 0.;
        real_type q_max = 0.;
        int nb_live = 0;
        for(std::size_t g = 0; g < entry.gen_ids.size(); ++g){
            const int gen_id = entry.gen_ids[g];
            if(is_gen_off(gen_id)) continue;  // disconnected by this row: produces nothing
            ++nb_live;
            q_min += generators.get_min_q(gen_id);
            q_max += generators.get_max_q(gen_id);
        }
        for(std::size_t s = 0; s < entry.storage_ids.size(); ++s){
            // no row disconnects a storage unit: always live
            const int storage_id = entry.storage_ids[s];
            ++nb_live;
            q_min += storages.get_min_q(storage_id);
            q_max += storages.get_max_q(storage_id);
        }
        for(std::size_t s = 0; s < entry.station_ids.size(); ++s){
            const int hvdc_id = entry.station_ids[s].first;
            const int side = entry.station_ids[s].second;
            ++nb_live;
            q_min += hvdcs.get_station_min_q_mvar(hvdc_id, side);
            q_max += hvdcs.get_station_max_q_mvar(hvdc_id, side);
        }
        if(!entry.svc_ids.empty() && entry.bus_solver < static_cast<int>(V.size())){
            // an SVC's capability is a susceptance range, so it is worth what the voltage
            // makes it worth: q = b . |V|^2 (generator convention), pu -> MVAr
            const real_type v2_sn = std::norm(V(entry.bus_solver)) * sn_mva;
            for(std::size_t s = 0; s < entry.svc_ids.size(); ++s){
                const int svc_id = entry.svc_ids[s];
                ++nb_live;
                q_min += svcs.get_b_min(svc_id) * v2_sn;
                q_max += svcs.get_b_max(svc_id) * v2_sn;
            }
        }
        if(nb_live == 0) continue;  // nothing holds this bus in this row

        // ---- ... and what it had to produce ---------------------------------------
        // (see the file comment: the raw reactive residual, which is the mismatch with
        // every controller's own injection added back)
        real_type q_bus = std::imag(bus_mismatch(entry.bus_solver));
        for(std::size_t c = 0; c < entry.ctrl_pos.size(); ++c){
            const int pos = entry.ctrl_pos[c];
            if(pos >= static_cast<int>(controller_q.size())) continue;  // no controller data
            q_bus += controller_q(pos);
        }
        q_bus *= sn_mva;
        if(!std::isfinite(q_bus)) continue;

        if(std::isfinite(q_min) && (q_bus < q_min - tol_mvar)){
            out.push_back(LimitViolation{ViolationElementType::BUS, entry.bus_grid, 0,
                                         LimitViolationType::LOW_Q, q_bus, q_min,
                                         entry.sub_name});
        } else if(std::isfinite(q_max) && (q_bus > q_max + tol_mvar)){
            out.push_back(LimitViolation{ViolationElementType::BUS, entry.bus_grid, 0,
                                         LimitViolationType::HIGH_Q, q_bus, q_max,
                                         entry.sub_name});
        }
    }
}

}  // namespace bus_q_check
}  // namespace ls2g

#endif  // BUSQCHECK_H
