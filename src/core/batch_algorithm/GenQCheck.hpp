// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENQCHECK_H
#define GENQCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace ls2g {

/**
 * Post-solve reactive-limit check for the VOLTAGE-REGULATING generators of a batch.
 *
 * WHAT THIS ANSWERS. A voltage-regulating generator has no reactive setpoint: its Q is
 * whatever it takes to hold its set-point, and lightsim2grid never clamps it. So a
 * converged row may well hold a bus at 1.05 pu with a generator producing three times
 * the reactive power its Qmax allows -- a perfectly converged, physically impossible
 * answer. This reports those, per row, exactly as PowSyBl OpenLoadFlow's
 * `ReactiveLimits` outer loop would DETECT them (it would then switch the bus PV -> PQ
 * and re-solve; nothing here re-solves anything, it only reports).
 *
 * WHY IT LIVES HERE AND NOT ON LSGrid. A single solve already publishes every
 * generator's reactive output (`GenInfo::res_q_mvar`, see LSGrid::compute_results), so
 * the check is a two-line comparison in Python. A batch publishes no such thing: it
 * keeps the voltages and drops everything else, precisely because reconstructing the
 * per-element results of ten thousand rows is what the batch classes exist not to do.
 * So for a batch the reactive output has to be re-derived, row by row, which is what
 * this does -- from exactly the two things LSGrid::compute_results itself uses.
 *
 * THE TWO SOURCES OF A GENERATOR'S Q, and they are the same two as in
 * LSGrid::compute_results (which is the point -- a row of a TimeSeries whose injection
 * is the grid's own state must report the number a plain `ac_pf` publishes):
 *
 *   1. the algorithm solved for it. A generator that is a controller of a
 *      `VoltageControl` group (it regulates a remote bus, or it stands on a bus a group
 *      regulates) has its own reactive injection as a Jacobian unknown. Its converged
 *      value is read straight off `BaseAlgo::get_controller_q()` (pu -> MVAr), the same
 *      write-back LSGrid::_write_back_controller_q does.
 *   2. it takes a share of its bus' reactive residual. The ordinary PV generator: the
 *      bus' Q equation is absent from the system, so `mis_bus_.imag()` at that bus IS
 *      what the machines standing on it produced (see LSGrid::_fill_bus_mismatch_ac for
 *      why the algorithm's own mismatch is the right, and not merely the cheap, thing to
 *      read). Several machines on one bus split it proportionally to their reactive
 *      RANGE, the rule mirrored from LSGrid::_split_q_residual_per_bus -- including the
 *      hvdc converter stations, which take their own share of that same residual and so
 *      have to be in the denominator even though no limit is checked on them here.
 *
 * WHAT IS DELIBERATELY NOT CHECKED: a generator that does not regulate voltage (its Q
 * is its own setpoint -- a violation there is an input error, not a control action), an
 * SVC (its `b_min`/`b_max` are stored but never enforced, see SvcContainer) and an hvdc
 * converter station. Generators only, which is what `ReactiveLimits` is about.
 */
namespace gen_q_check {

/// One voltage-regulating generator, and where its reactive output comes from. Built
/// once per compute() by `build_gen_q_plan`; nothing here varies from row to row.
struct GenQEntry
{
    int gen_id = -1;
    int bus_solver = -1;      ///< the generator's OWN bus, solver numbering
    real_type span = 0.;      ///< qmax - qmin (MVAr): the residual-sharing key
    real_type min_q_mvar = 0.;
    real_type max_q_mvar = 0.;
    /// index in the algorithm's controller list (`get_controller_q()`), or -1 for a
    /// generator that takes a share of its bus' reactive residual instead
    int ctrl_pos = -1;
};

/// The generators of one bus, as one contiguous range of `GenQPlan::gens` plus what
/// else on that bus takes a share of the same reactive residual.
struct GenQBusGroup
{
    int bus_solver = -1;
    int first = 0;            ///< range [first, last) into GenQPlan::gens
    int last = 0;
    real_type other_span = 0.;  ///< summed reactive range of the non-generator sharers
    int nb_other = 0;           ///< how many they are (hvdc converter stations)
    bool other_all_finite = true;
};

struct GenQPlan
{
    std::vector<GenQEntry> gens;      ///< sorted by bus_solver, then by generator id
    std::vector<GenQBusGroup> groups; ///< one per distinct bus of `gens`, same order
    /// whether any generator of the plan reads its reactive output off the algorithm's
    /// controller list: false on an ordinary grid of local PV machines, and then a row
    /// never has to ask for that list at all
    bool needs_controller_q = false;

    bool empty() const { return gens.empty(); }
    void clear() { gens.clear(); groups.clear(); needs_controller_q = false; }
};

/**
 * Work out, once, which generators can violate a reactive limit at all and where each
 * one's reactive output will have to be read from.
 *
 * `id_me_to_solver` and `ctrl` must both describe the labelling the batch solves in
 * (`active_layout()`); `ctrl` is the plan's controller list, whose order is the order of
 * `BaseAlgo::get_controller_q()`.
 */
inline void build_gen_q_plan(const LSGrid & grid_model,
                             const SolverBusIdVect & id_me_to_solver,
                             const VoltageControlSolverData & ctrl,
                             GenQPlan & out)
{
    out.clear();
    const GeneratorContainer & generators = grid_model.get_generators();
    const int nb_gen = generators.nb();
    if(nb_gen == 0) return;

    // gen id -> its position in the controller list (-1: not a controller)
    std::vector<int> ctrl_pos_of_gen(static_cast<std::size_t>(nb_gen), -1);
    const int nb_ctrl = ctrl.n_controllers();
    for(int c = 0; c < nb_ctrl; ++c){
        if(ctrl.kind(c) != VoltageControlSolverData::GEN) continue;
        const int gen_id = ctrl.elem_id(c);
        if(gen_id >= 0 && gen_id < nb_gen) ctrl_pos_of_gen[static_cast<std::size_t>(gen_id)] = c;
    }

    const GlobalBusIdVect & gen_buses = generators.get_bus_id();
    out.gens.reserve(static_cast<std::size_t>(nb_gen));
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        // an ACTIVE voltage regulator: connected, regulating, not treated as off. A
        // generator that does not regulate injects its own reactive setpoint and is not
        // this check's business.
        if(!generators.is_voltage_controller(gen_id)) continue;

        const int bus_me = gen_buses(gen_id).cast_int();
        if(bus_me == BaseConstants::_deactivated_bus_id) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        const int ctrl_pos = ctrl_pos_of_gen[static_cast<std::size_t>(gen_id)];
        // A remote controller MUST be in the controller list (LSGrid::ac_pf refuses to
        // solve otherwise, naming the offending element). If one is not -- an algorithm
        // that cannot do remote voltage control, whose plan stays empty -- there is no
        // way to know what it produced, so it is left out rather than reported wrong.
        if(ctrl_pos < 0 && !generators.is_local_voltage_controller(gen_id)) continue;

        GenQEntry entry;
        entry.gen_id = gen_id;
        entry.bus_solver = bus_solver;
        entry.min_q_mvar = generators.get_min_q(gen_id);
        entry.max_q_mvar = generators.get_max_q(gen_id);
        entry.span = entry.max_q_mvar - entry.min_q_mvar;
        entry.ctrl_pos = ctrl_pos;
        if(ctrl_pos >= 0) out.needs_controller_q = true;
        out.gens.push_back(entry);
    }
    if(out.gens.empty()) return;

    // grouped by bus, exactly like LSGrid::_split_q_residual_per_bus does it: one sort
    // over a packed (bus, position) key, so the order within a bus stays the insertion
    // (generator id) order.
    std::vector<std::uint64_t> by_bus;
    by_bus.reserve(out.gens.size());
    for(std::size_t k = 0; k < out.gens.size(); ++k){
        by_bus.push_back((static_cast<std::uint64_t>(out.gens[k].bus_solver) << 32) |
                         static_cast<std::uint64_t>(k));
    }
    std::sort(by_bus.begin(), by_bus.end());
    std::vector<GenQEntry> sorted;
    sorted.reserve(out.gens.size());
    for(std::size_t k = 0; k < by_bus.size(); ++k){
        sorted.push_back(out.gens[static_cast<std::size_t>(by_bus[k] & 0xffffffffu)]);
    }
    out.gens.swap(sorted);

    // the other machines that share a bus' reactive residual: an hvdc converter station
    // regulating the bus it stands on WITHOUT being a group controller (a station that
    // joined a group has its own Q unknown instead, like a remote generator). Voltage
    // mode SVCs are always group controllers, so they are never in this list.
    std::vector<std::pair<int, real_type> > others;
    const HvdcLineContainer & hvdcs = grid_model.get_dclines();
    const int nb_hvdc = hvdcs.nb();
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        for(int side = 1; side <= 2; ++side){
            if(!hvdcs.station_is_voltage_controller(hvdc_id, side)) continue;
            bool is_ctrl = false;
            const int want_kind = (side == 1) ? VoltageControlSolverData::HVDC_SIDE_1
                                              : VoltageControlSolverData::HVDC_SIDE_2;
            for(int c = 0; c < nb_ctrl; ++c){
                if(ctrl.kind(c) == want_kind && ctrl.elem_id(c) == hvdc_id){ is_ctrl = true; break; }
            }
            if(is_ctrl) continue;
            const int bus_me = hvdcs.get_station_bus(hvdc_id, side).cast_int();
            if(bus_me == BaseConstants::_deactivated_bus_id) continue;
            const int bus_solver = id_me_to_solver[bus_me].cast_int();
            if(bus_solver == BaseConstants::_deactivated_bus_id) continue;
            others.push_back(std::make_pair(bus_solver, hvdcs.get_station_q_range_mvar(hvdc_id, side)));
        }
    }
    std::sort(others.begin(), others.end());

    for(std::size_t first = 0; first < out.gens.size(); ){
        const int bus_solver = out.gens[first].bus_solver;
        std::size_t last = first + 1;
        while((last < out.gens.size()) && (out.gens[last].bus_solver == bus_solver)) ++last;

        GenQBusGroup group;
        group.bus_solver = bus_solver;
        group.first = static_cast<int>(first);
        group.last = static_cast<int>(last);
        for(std::size_t k = 0; k < others.size(); ++k){
            if(others[k].first != bus_solver) continue;
            ++group.nb_other;
            if(!std::isfinite(others[k].second)) group.other_all_finite = false;
            else group.other_span += others[k].second;
        }
        out.groups.push_back(group);
        first = last;
    }
}

/**
 * Append to `out` one LimitViolation per voltage-regulating generator whose reactive
 * output leaves [min_q, max_q] by more than `tol_mvar`, for ONE converged row.
 *
 * `bus_mismatch` is the algorithm's own per-bus mismatch (solver numbering, pu, see
 * BaseAlgo::get_bus_mismatch) and `controller_q` its converged per-controller reactive
 * injection (pu). `masked_solver_ids` is this row's masked (stranded) solver buses --
 * sorted, may be nullptr -- whose voltage is reported as 0 and whose mismatch means
 * nothing. `is_gen_off(gen_id)` tells whether this row disconnected that generator (a
 * generator contingency); such a generator produces nothing and leaves the sharing
 * denominator with its span.
 */
template<class IsGenOff>
inline void check_gen_q_violations(const GenQPlan & plan,
                                   const Eigen::Ref<const CplxVect> & bus_mismatch,
                                   const RealVect & controller_q,
                                   real_type sn_mva,
                                   real_type tol_mvar,
                                   const std::vector<int> * masked_solver_ids,
                                   const std::vector<std::string> & gen_names,
                                   IsGenOff is_gen_off,
                                   std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;

    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();
    // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
    auto is_masked = [&](int bus){
        return has_masked && std::binary_search(masked_solver_ids->begin(), masked_solver_ids->end(), bus);
    };
    auto name_of = [&](int gen_id){
        return static_cast<std::size_t>(gen_id) < gen_names.size() ? gen_names[static_cast<std::size_t>(gen_id)]
                                                                  : std::string();
    };
    auto emit = [&](const GenQEntry & entry, real_type q_mvar){
        if(!std::isfinite(q_mvar)) return;
        if(std::isfinite(entry.min_q_mvar) && (q_mvar < entry.min_q_mvar - tol_mvar)){
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::LOW_Q, q_mvar, entry.min_q_mvar,
                                         name_of(entry.gen_id)});
        } else if(std::isfinite(entry.max_q_mvar) && (q_mvar > entry.max_q_mvar + tol_mvar)){
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::HIGH_Q, q_mvar, entry.max_q_mvar,
                                         name_of(entry.gen_id)});
        }
    };

    // ---- 1. the generators the algorithm solved for ------------------------------
    for(std::size_t k = 0; k < plan.gens.size(); ++k){
        const GenQEntry & entry = plan.gens[k];
        if(entry.ctrl_pos < 0) continue;
        if(entry.ctrl_pos >= static_cast<int>(controller_q.size())) continue;  // no controller data
        if(is_gen_off(entry.gen_id)) continue;
        if(is_masked(entry.bus_solver)) continue;
        emit(entry, controller_q(entry.ctrl_pos) * sn_mva);
    }

    // ---- 2. the ordinary PV generators, sharing their bus' reactive residual ------
    if(bus_mismatch.size() == 0) return;
    const real_type eps_q = 1e-8;
    for(std::size_t g = 0; g < plan.groups.size(); ++g){
        const GenQBusGroup & group = plan.groups[g];
        if(group.bus_solver >= static_cast<int>(bus_mismatch.size())) continue;
        if(is_masked(group.bus_solver)) continue;

        // pass 1: who is live here, and what the sharing denominator is. Mirrors
        // LSGrid::_split_q_residual_per_bus: the two degenerate cases (every range
        // zero, any range not finite) are one and the same -- nothing to weigh the
        // participants against -- and collapse to an equal split.
        int nb_here = group.nb_other;
        real_type total_span = group.other_span;
        bool all_finite = group.other_all_finite;
        for(int k = group.first; k < group.last; ++k){
            const GenQEntry & entry = plan.gens[static_cast<std::size_t>(k)];
            if(entry.ctrl_pos >= 0) continue;  // its Q is a Jacobian unknown, not a share
            if(is_gen_off(entry.gen_id)) continue;
            ++nb_here;
            if(!std::isfinite(entry.span)) all_finite = false;
            else total_span += entry.span;
        }
        if(nb_here == 0) continue;  // nothing regulating this bus in this row

        const real_type q_to_absorb = std::imag(bus_mismatch(group.bus_solver)) * sn_mva;
        const real_type nb_here_r = static_cast<real_type>(nb_here);
        for(int k = group.first; k < group.last; ++k){
            const GenQEntry & entry = plan.gens[static_cast<std::size_t>(k)];
            if(entry.ctrl_pos >= 0) continue;
            if(is_gen_off(entry.gen_id)) continue;
            real_type q;
            if(nb_here == 1) q = q_to_absorb;
            else if(!all_finite) q = q_to_absorb / nb_here_r;
            else q = q_to_absorb * (entry.span + eps_q) / (total_span + nb_here_r * eps_q);
            emit(entry, q);
        }
    }
}

}  // namespace gen_q_check
}  // namespace ls2g

#endif  // GENQCHECK_H
