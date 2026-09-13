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
#include <vector>

namespace ls2g {

/**
 * Post-solve reactive-capability check, PER BUS, for the buses whose voltage is held by
 * generators.
 *
 * WHAT THIS ANSWERS. A voltage-regulating generator has no reactive setpoint: its Q is
 * whatever it takes to hold the set-point, and lightsim2grid never clamps it. So a row
 * may converge holding a bus at 1.05 pu with its machines having to produce three times
 * the reactive power they own. That is not a limit somebody chose to exceed -- it is a
 * solution the grid cannot reach at all, because a machine cannot produce reactive power
 * it does not have. See ViolationCategory::PHYSICAL. It is the condition PowSyBl
 * OpenLoadFlow's `ReactiveLimits` outer loop acts on (it would switch the bus PV -> PQ
 * and re-solve); nothing here re-solves anything, it only reports.
 *
 * WHY PER BUS AND NOT PER MACHINE. The reactive power a bus needs is a fact about the
 * solution; how it is divided between several machines standing on that bus is not -- the
 * solver never decides it, `LSGrid::_split_q_residual_per_bus` does, by a sharing
 * convention (proportional to each machine's reactive range). Checking machine by machine
 * would therefore report the convention rather than the grid: two 20 MVAr machines
 * covering 30 MVAr together is fine, and a per-machine check on the same solution says
 * both are 5 MVAr over. So the question asked here is the one that has an answer
 * independent of the convention: did this bus need more reactive power than the sum of
 * what its machines can produce?
 *
 * WHY IT LIVES HERE AND NOT ON LSGrid. A single solve publishes every generator's
 * reactive output (`GenInfo::res_q_mvar`, see LSGrid::compute_results), so the check is a
 * few lines in Python. A batch publishes no such thing: it keeps the voltages and drops
 * everything else, precisely because reconstructing the per-element results of ten
 * thousand rows is what the batch classes exist not to do. So for a batch the bus'
 * reactive power has to be re-derived, row by row, which is what this does.
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
 * reactive residual at `b`, and the raw residual IS what the machines pinning that bus
 * produced, since a regulating machine's Q is never part of `Sbus`. One expression covers
 * both cases -- an ordinary PV bus (no controller, so just the mismatch) and a bus a
 * control group holds (every machine on it is a controller, so just their Q).
 *
 * WHAT IS DELIBERATELY NOT CHECKED. A bus whose reactive power is not produced by
 * generators alone: one carrying a voltage-mode SVC or a voltage-regulating hvdc
 * converter station is skipped entirely, because `q_bus` would then include that
 * element's contribution while the summed limits would not (an SVC's capability is a
 * susceptance range that is never enforced anywhere -- see SvcContainer -- and a station's
 * is its own [min_q, max_q]). Skipping is the honest answer; extending the sum to those
 * two families is the natural next step. Generators that do not regulate voltage are not
 * in the picture at all: their Q is their own setpoint, part of `Sbus`, and a violation
 * there would be an input error rather than something the solve produced.
 */
namespace bus_q_check {

/// One bus whose voltage is held by generators, and everything a row needs to check it.
/// Built once per compute() by `build_bus_q_plan`; nothing here varies from row to row.
struct BusQEntry
{
    int bus_solver = -1;          ///< the bus, solver numbering
    int bus_grid = -1;            ///< the same bus, grid numbering (what is reported)
    /// the generators whose reactive output this bus' Q is made of: a row that
    /// disconnects one takes its limits out of the sum
    std::vector<int> gen_ids;
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
 * Work out, once, which buses can be checked at all and what each one's reactive power
 * will have to be read from.
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
    const GeneratorContainer & generators = grid_model.get_generators();
    const int nb_gen = generators.nb();
    if(nb_gen == 0) return;

    const int nb_bus_grid = static_cast<int>(grid_model.total_bus());
    // per GRID bus: the regulating generators standing on it, and whether anything else
    // standing on it also has a free reactive output (which makes the bus unanswerable
    // with generator limits alone -- see the note on what is not checked)
    std::vector<std::vector<int> > gens_of_bus(static_cast<std::size_t>(nb_bus_grid));
    std::vector<char> has_other_source(static_cast<std::size_t>(nb_bus_grid), 0);

    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        // an ACTIVE voltage regulator: connected, regulating, and not treated as off
        if(!generators.is_voltage_controller(gen_id)) continue;
        const int bus_me = gen_buses(gen_id).cast_int();
        if(bus_me < 0 || bus_me >= nb_bus_grid) continue;
        gens_of_bus[static_cast<std::size_t>(bus_me)].push_back(gen_id);
    }

    // voltage-mode SVCs and voltage-regulating hvdc converter stations: not checked, and
    // they disqualify the bus they stand on
    const SvcContainer & svcs = grid_model.get_svcs();
    const int nb_svc = svcs.nb();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svcs.is_voltage_controller(svc_id)) continue;
        const int bus_me = svcs.get_bus_id()(svc_id).cast_int();
        if(bus_me >= 0 && bus_me < nb_bus_grid) has_other_source[static_cast<std::size_t>(bus_me)] = 1;
    }
    const HvdcLineContainer & hvdcs = grid_model.get_dclines();
    const int nb_hvdc = hvdcs.nb();
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            if(!hvdcs.station_is_voltage_controller(hvdc_id, side)) continue;
            const int bus_me = hvdcs.get_station_bus(hvdc_id, side).cast_int();
            if(bus_me >= 0 && bus_me < nb_bus_grid) has_other_source[static_cast<std::size_t>(bus_me)] = 1;
        }
    }

    // controller list, by the bus the controller STANDS on (`ctrl.bus`, which is where the
    // extension injects its -i.Q_c), so a bus' entry can add its controllers' Q back
    const int nb_ctrl = ctrl.n_controllers();
    std::vector<std::vector<int> > ctrl_of_bus_solver;

    const SubstationContainer & subs = grid_model.get_substations();
    const std::vector<std::string> & sub_names = subs.get_sub_names();  // empty if never set

    for(int bus_me = 0; bus_me < nb_bus_grid; ++bus_me){
        const std::vector<int> & gens_here = gens_of_bus[static_cast<std::size_t>(bus_me)];
        if(gens_here.empty()) continue;
        if(has_other_source[static_cast<std::size_t>(bus_me)]) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        BusQEntry entry;
        entry.bus_solver = bus_solver;
        entry.bus_grid = bus_me;
        entry.gen_ids = gens_here;
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
 * Append to `out` one LimitViolation per bus whose voltage-regulating generators had to
 * produce more (or less) reactive power than they own, for ONE converged row.
 *
 * `bus_mismatch` is the algorithm's own per-bus mismatch (solver numbering, pu, see
 * BaseAlgo::get_bus_mismatch) and `controller_q` its converged per-controller reactive
 * injection (pu). `masked_solver_ids` is this row's masked (stranded) solver buses --
 * sorted, may be nullptr -- whose voltage is reported as 0 and whose mismatch means
 * nothing. `is_gen_off(gen_id)` tells whether this row disconnected that generator (a
 * generator contingency): it produces nothing, so it leaves the summed capability, and a
 * bus whose every generator is off is not checked at all (it is an ordinary PQ bus in
 * that row, and the residual there is nobody's reactive output).
 */
template<class IsGenOff>
inline void check_bus_q_violations(const BusQPlan & plan,
                                   const GeneratorContainer & generators,
                                   const Eigen::Ref<const CplxVect> & bus_mismatch,
                                   const RealVect & controller_q,
                                   real_type sn_mva,
                                   real_type tol_mvar,
                                   const std::vector<int> * masked_solver_ids,
                                   IsGenOff is_gen_off,
                                   std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;
    if(bus_mismatch.size() == 0) return;

    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();

    for(std::size_t k = 0; k < plan.buses.size(); ++k){
        const BusQEntry & entry = plan.buses[k];
        if(entry.bus_solver >= static_cast<int>(bus_mismatch.size())) continue;
        // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
        if(has_masked && std::binary_search(masked_solver_ids->begin(), masked_solver_ids->end(),
                                            entry.bus_solver)) continue;

        // what this bus' machines can produce, this row
        real_type q_min = 0.;
        real_type q_max = 0.;
        int nb_live = 0;
        for(std::size_t g = 0; g < entry.gen_ids.size(); ++g){
            const int gen_id = entry.gen_ids[g];
            if(is_gen_off(gen_id)) continue;
            ++nb_live;
            q_min += generators.get_min_q(gen_id);
            q_max += generators.get_max_q(gen_id);
        }
        if(nb_live == 0) continue;  // nothing regulates this bus in this row

        // ... and what it had to produce (see the file comment: the raw reactive residual,
        // which is the mismatch with every controller's own injection added back)
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
