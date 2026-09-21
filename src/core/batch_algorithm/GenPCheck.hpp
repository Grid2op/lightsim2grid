// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENPCHECK_H
#define GENPCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace ls2g {

/**
 * Post-solve active-power check of the machines that carry the DISTRIBUTED SLACK -- one
 * third of `compute_physical_violations` (the others are BusQCheck.hpp and HvdcPCheck.hpp).
 *
 * BOTH FAMILIES, IN ONE PLACE. Generators and storage units take a share of the slack
 * under the same rule (SlackParticipation, which LSGrid drives for both), so they are
 * checked together and CANNOT be checked separately: the share a machine takes is a
 * fraction of the raw participation of its whole bus, families included, so leaving the
 * storage units out of the sum would overstate every generator's share (see GenPPlan).
 * A storage unit's setpoints are in the LOAD convention in this library; everything here
 * works in the GENERATOR one, so its target is negated on the way in and the value
 * reported is what the unit injects -- the convention its limits are given in.
 *
 * WHAT THIS ANSWERS. lightsim2grid does not distribute the slack between solves: it solves
 * it, inside the Newton system (`MultiSlack` -- the absorbed power is an unknown of the
 * Jacobian, shared out by fixed participation factors). Nothing in that formulation knows
 * what a machine can deliver, so a participating generator's converged active power
 *
 *     p = target_p + its share of the imbalance
 *
 * can land anywhere, including beyond `max_p_mw` or below `min_p_mw`. The row then
 * describes a machine producing power it does not have -- the same kind of statement as a
 * bus needing reactive power nobody can give it (ViolationCategory::PHYSICAL), and the
 * condition OpenLoadFlow's `DistributedSlack` outer loop acts on: it would take the
 * saturated unit out of the distribution, re-share what is left, and re-solve. Nothing
 * here re-solves anything and nothing is clamped: a row is reported, not fixed.
 *
 * WHY PER MACHINE, where the reactive check is per bus. The two look symmetric and are
 * not. A bus' reactive power is split between its machines by a CONVENTION
 * (`LSGrid::_split_q_residual_per_bus`, proportional to each one's reactive range), so a
 * per-machine reactive check would report that convention. The active split is not a
 * convention: each generator's participation factor is input data the caller chose
 * (`add_gen_slackbus(gen_id, weight)`), and the solver distributes by it -- so
 * `target_p + share` is that machine's own power, and its own limits are the right thing
 * to compare it against. Mirrors `GeneratorContainer::set_p_slack`, which is what
 * publishes the very same number after a single solve.
 *
 * WHICH MACHINES ARE LOOKED AT. Those that actually take part in the distribution:
 * connected, flagged slack, with a nonzero weight, and given at least one finite limit. A
 * machine that does not participate keeps its target power exactly -- so a violation
 * there would be an input error (the caller asked for a power the machine does not have),
 * not something the solve produced, and it is left to the caller to notice. Non-finite
 * limits (the default, where the grid was never given any) are skipped per machine and per
 * side, exactly like a branch with no thermal rating.
 */
namespace gen_p_check {

/// One machine taking part in the distributed slack, and everything a row needs to check
/// it. Built once per compute() by `build_gen_p_plan`.
struct GenPEntry
{
    /// GENERATOR or STORAGE -- which container `el_id` indexes, and what the violation
    /// is reported on
    ViolationElementType el_type = ViolationElementType::GENERATOR;
    int el_id = -1;
    int bus_solver = -1;          ///< the bus it stands on, solver numbering
    real_type slack_weight = 0.;  ///< its own participation factor
    real_type min_p_mw = 0.;      ///< NaN where no limit was given, generator convention
    real_type max_p_mw = 0.;
    /// STORAGE only: its active setpoint (MW, generator convention -- ie the negated
    /// load-convention one the container stores). A generator's is the ROW's, which only
    /// the caller knows, and comes through `target_p_of`; no batch varies a storage unit's
    /// injection or disconnects one, so the grid's own is every row's.
    real_type target_p_mw = 0.;
    /// LSGrid::set_gen_names / set_storage_names, empty if never set
    std::string name;
};

/// One machine taking part in the distribution, limits or not. Needed because the share
/// each machine gets is a fraction of the RAW participation factors (see `GenPPlan`), and
/// a machine with no limit still takes part in that total.
struct GenPParticipant
{
    ViolationElementType el_type = ViolationElementType::GENERATOR;
    int el_id = -1;
    int bus_solver = -1;
    real_type slack_weight = 0.;
};

/**
 * WHY TWO LISTS, and the one subtlety of this whole file: there are two per-bus weight
 * vectors and they are not the same one.
 *
 *   * the vector the SOLVER is given (`SolverBusLayout::slack_weights`, or a row's own
 *     re-derived one) is NORMALIZED -- it sums to 1 over the participating buses, which is
 *     what makes `slack_absorbed * w(bus)` the share that bus absorbed;
 *   * the one `GeneratorContainer::set_p_slack` divides by is the RAW per-bus sum of the
 *     participation factors the caller gave (`bus_slack_weight_`), so that
 *     `raw_w(gen) / raw_total(bus)` is a fraction of ONE BUS.
 *
 * Only the normalized one reaches a batch row, so the raw total of a bus is recovered as
 * `w_norm(bus) * (the row's total raw weight)` -- and that total is what `participants`
 * is for: the sum, over every machine the row leaves participating, of its own factor. A
 * machine with no limit can never be reported and is still in it, because it still takes
 * its share -- and so is a machine of the OTHER family, for the same reason: both
 * families are in `bus_slack_weight_`, so both must be in the total it is recovered with.
 */
struct GenPPlan
{
    std::vector<GenPEntry> gens;               ///< those that CAN be reported
    std::vector<GenPParticipant> participants; ///< ... and everything that takes a share

    bool empty() const { return gens.empty(); }
    void clear() { gens.clear(); participants.clear(); }
};

/**
 * One family's contribution to the plan: every participating element goes into
 * `participants`, and the ones with a limit into `gens` as well. `Container` is a
 * GeneratorContainer or a StorageContainer -- the two answer the same questions, and the
 * only thing that differs is `target_sign` (+1 for a container in generator convention,
 * -1 for the load-convention storage units, whose target is negated into the generator
 * one this whole file works in).
 */
template<class Container>
inline void add_family_to_gen_p_plan(const Container & container,
                                     ViolationElementType el_type,
                                     real_type target_sign,
                                     const SolverBusIdVect & id_me_to_solver,
                                     GenPPlan & out)
{
    const int nb_el = container.nb();
    if(nb_el == 0) return;
    // nothing to compare THIS family against: the grid was never given limits for it. Its
    // elements still take their share, so they still go into `participants`.
    const bool has_limits = (container.get_p_min_mw().size() != 0) ||
                            (container.get_p_max_mw().size() != 0);

    const std::vector<bool> & status = container.get_status();
    const GlobalBusIdVect & buses = container.get_bus_id();
    const std::vector<std::string> & names = container.get_names();  // empty if never set
    const Eigen::Ref<const RealVect> targets = container.get_target_p();

    for(int el_id = 0; el_id < nb_el; ++el_id){
        if(!status[el_id]) continue;
        if(!container.is_slack(el_id)) continue;              // takes no part in the distribution
        const real_type weight = container.get_slack_weight(el_id);
        if(std::abs(weight) < 1e-12) continue;                // ... nor does a zero weight

        const int bus_me = buses(el_id).cast_int();
        if(bus_me == BaseConstants::_deactivated_bus_id) continue;
        const int bus_solver = id_me_to_solver[bus_me].cast_int();
        if(bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        // every participant counts towards the row's total, limits or not, family or not
        GenPParticipant part;
        part.el_type = el_type;
        part.el_id = el_id;
        part.bus_solver = bus_solver;
        part.slack_weight = weight;
        out.participants.push_back(part);

        if(!has_limits) continue;
        const real_type min_p = container.get_min_p(el_id);
        const real_type max_p = container.get_max_p(el_id);
        if(!std::isfinite(min_p) && !std::isfinite(max_p)) continue;  // no limit on this one

        GenPEntry entry;
        entry.el_type = el_type;
        entry.el_id = el_id;
        entry.bus_solver = bus_solver;
        entry.slack_weight = weight;
        entry.min_p_mw = min_p;
        entry.max_p_mw = max_p;
        entry.target_p_mw = target_sign * targets.coeff(el_id);
        if(static_cast<std::size_t>(el_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(el_id)];
        }
        out.gens.push_back(entry);
    }
}

/**
 * Work out, once, which machines can be reported at all.
 *
 * `id_me_to_solver` must describe the labelling the batch solves in (`active_layout()`) --
 * the same one the per-bus vectors handed to `check_gen_p_violations` are indexed in.
 */
inline void build_gen_p_plan(const LSGrid & grid_model,
                             const SolverBusIdVect & id_me_to_solver,
                             GenPPlan & out)
{
    out.clear();
    // the generators first, so that a report's order follows the container order of each
    // family and stays stable
    add_family_to_gen_p_plan(grid_model.get_generators(), ViolationElementType::GENERATOR,
                             1., id_me_to_solver, out);
    add_family_to_gen_p_plan(grid_model.get_storages(), ViolationElementType::STORAGE,
                             -1., id_me_to_solver, out);
    // no machine can be reported: the participants alone are of no use to anyone
    if(out.gens.empty()) out.clear();
}

/**
 * What a row leaves behind about its distributed slack, as the three solver families
 * express it -- enough to recompute any bus' share without materializing a vector (this is
 * read once per generator, from several threads at a time).
 */
struct SlackShareInputs
{
    /// The two per-bus vectors are REFERENCES into the row's own state, so this is built
    /// where they are alive and never outlives them. `bus_mismatch` may legitimately be
    /// empty (a DC solve leaves none); `ac` says which of the two branches below applies.
    SlackShareInputs(const Eigen::Ref<const CplxVect> & bus_mismatch_,
                     const Eigen::Ref<const RealVect> & bus_slack_weight_,
                     real_type sn_mva_,
                     bool ac_)
        : bus_mismatch(bus_mismatch_),
          bus_slack_weight(bus_slack_weight_),
          sn_mva(sn_mva_),
          ac(ac_)
    {}

    /// AC: the algorithm's per-bus mismatch (pu, solver numbering) and the converged value
    /// of the `MultiSlack` unknown. `LSGrid::_fill_bus_mismatch_ac` subtracts the second
    /// from the first to recover the raw active residual, which is what the slack machines
    /// of that bus produced on top of their targets.
    Eigen::Ref<const CplxVect> bus_mismatch;
    real_type slack_absorbed = 0.;
    /// DC: there is no mismatch to read -- the whole imbalance `-sum(Sbus)` is shared out by
    /// the (normalized) per-bus weights, see `LSGrid::_fill_bus_mismatch_dc`. In MW.
    real_type dc_imbalance_mw = 0.;
    /// this row's per-bus participation, NORMALIZED -- the vector the solver was given, the
    /// row's own rather than the layout's where a generator contingency re-derived it
    /// (BaseBatchSweep::_row_slack_weights / _masked_slack_weights). See GenPPlan for why
    /// its per-bus entry is not the denominator a machine's share is taken over.
    Eigen::Ref<const RealVect> bus_slack_weight;
    real_type sn_mva = 1.;
    bool ac = true;

    /// the active power the slack machines of `bus` produced on top of their targets (MW)
    real_type node_mismatch_mw(int bus) const {
        const real_type w = (bus < static_cast<int>(bus_slack_weight.size())) ? bus_slack_weight(bus) : 0.;
        if(!ac) return w * dc_imbalance_mw;
        if(bus >= static_cast<int>(bus_mismatch.size())) return 0.;
        return (std::real(bus_mismatch(bus)) - slack_absorbed * w) * sn_mva;
    }
};

/**
 * Append to `out` one LimitViolation per participating machine whose converged active
 * power left its limits, for ONE converged row.
 *
 * `slack` is what that row's solve left about the distribution (see SlackShareInputs);
 * `target_p_of(gen_id)` is this row's own active set-point for a GENERATOR, and
 * `is_gen_off(gen_id)` whether the row disconnected it. Both are asked about generators
 * only: no batch varies a storage unit's injection or disconnects one, so its target is
 * the grid's own, read into the entry when the plan was built. `masked_solver_ids` is this
 * row's masked (stranded) solver buses -- sorted, may be nullptr.
 */
template<class TargetPOf, class IsGenOff>
inline void check_gen_p_violations(const GenPPlan & plan,
                                   const SlackShareInputs & slack,
                                   real_type tol_mw,
                                   const std::vector<int> * masked_solver_ids,
                                   TargetPOf target_p_of,
                                   IsGenOff is_gen_off,
                                   std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;

    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();
    // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
    auto is_masked = [&](int bus){
        return has_masked && std::binary_search(masked_solver_ids->begin(),
                                                masked_solver_ids->end(), bus);
    };

    // a row takes a GENERATOR out of the distribution; a storage unit is never varied nor
    // disconnected by one, so it always keeps its share
    auto is_off = [&](ViolationElementType el_type, int el_id){
        return (el_type == ViolationElementType::GENERATOR) && is_gen_off(el_id);
    };

    // the row's total raw participation, over the machines it actually leaves participating
    // -- what turns the normalized per-bus weight back into a raw one, see GenPPlan. BOTH
    // families are in it: the per-bus weights the solver was given sum them together.
    real_type total_raw_w = 0.;
    for(std::size_t k = 0; k < plan.participants.size(); ++k){
        const GenPParticipant & part = plan.participants[k];
        if(is_masked(part.bus_solver)) continue;
        if(is_off(part.el_type, part.el_id)) continue;
        total_raw_w += part.slack_weight;
    }
    if(!(std::abs(total_raw_w) > 1e-12)) return;  // nothing left distributing anything

    for(std::size_t k = 0; k < plan.gens.size(); ++k){
        const GenPEntry & entry = plan.gens[k];
        if(is_masked(entry.bus_solver)) continue;
        if(is_off(entry.el_type, entry.el_id)) continue;  // disconnected by this row: produces nothing

        // its target, plus its share of what its bus had to make up -- exactly as
        // GeneratorContainer::set_p_slack computes it after a single solve, the raw per-bus
        // total written as `w_norm(bus) * total_raw_w`. Generator convention throughout, so
        // a storage unit's share is added to its target the same way round (see the file's
        // header: the target was negated when the plan was built).
        real_type p_mw = (entry.el_type == ViolationElementType::GENERATOR)
                         ? target_p_of(entry.el_id) : entry.target_p_mw;
        if(entry.bus_solver < static_cast<int>(slack.bus_slack_weight.size())){
            const real_type bus_raw_w = slack.bus_slack_weight(entry.bus_solver) * total_raw_w;
            if(std::abs(bus_raw_w) > 1e-12){
                p_mw += slack.node_mismatch_mw(entry.bus_solver) * entry.slack_weight / bus_raw_w;
            }
        }
        if(!std::isfinite(p_mw)) continue;

        if(std::isfinite(entry.min_p_mw) && (p_mw < entry.min_p_mw - tol_mw)){
            out.push_back(LimitViolation{entry.el_type, entry.el_id, 0,
                                         LimitViolationType::LOW_P, p_mw, entry.min_p_mw,
                                         entry.name});
        } else if(std::isfinite(entry.max_p_mw) && (p_mw > entry.max_p_mw + tol_mw)){
            out.push_back(LimitViolation{entry.el_type, entry.el_id, 0,
                                         LimitViolationType::HIGH_P, p_mw, entry.max_p_mw,
                                         entry.name});
        }
    }
}

}  // namespace gen_p_check
}  // namespace ls2g

#endif  // GENPCHECK_H
