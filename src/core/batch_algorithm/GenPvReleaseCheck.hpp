// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENPVRELEASECHECK_H
#define GENPVRELEASECHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace ls2g {

/**
 * Post-solve check of the PQ generators an outer loop pinned at a reactive limit: would
 * that loop switch them back to PV? One of the physical checks of
 * `compute_physical_violations` (with BusQCheck.hpp, HvdcPCheck.hpp and GenPCheck.hpp).
 *
 * WHAT THIS ANSWERS. OpenLoadFlow's `ReactiveLimits` outer loop has two directions. The
 * PV -> PQ one -- a regulating machine asked for more reactive power than it has -- is what
 * BusQCheck.hpp detects (LOW_Q / HIGH_Q). The PQ -> PV one is its mirror image: a machine
 * the loop pinned at its minimum (resp. maximum) reactive power, whose regulated voltage
 * then sits BELOW (resp. ABOVE) its target, is a machine that absorbs (resp. produces) too
 * much for that target and would regulate again if the grid let it. The converged solution
 * then assumes a control that the loop would not leave in place: the same kind of statement
 * as LOW_Q / HIGH_Q (ViolationCategory::PHYSICAL), and nothing here re-solves or switches
 * anything -- a row is reported, not fixed.
 *
 * ONLY THE MACHINES A CALLER FLAGGED. lightsim2grid never pins a machine itself, so it
 * cannot tell one that an outer loop pinned from one that was PQ to begin with -- a
 * fixed-Q unit is an input, and a voltage it does not regulate is nobody's business. The
 * flag is `GeneratorContainer::can_be_pv` (LSGrid::set_gen_can_be_pv), which the pypowsybl
 * converter fills from what `bake_outer_loops` froze; everything else is skipped.
 *
 * WHAT IS CHECKED. A flagged PQ generator whose reactive setpoint sits at one of its limits
 * (within `tol_mva`, the side being decided once per plan), with a reactive range of at
 * least MIN_REACTIVE_RANGE_MVAR (OpenLoadFlow's own plausibility floor: a narrower range
 * never regulates) and a meaningful target voltage. Per row, the magnitude of its regulated
 * bus is compared with that target: below it at the minimum, above it at the maximum, by
 * more than `tol_vm_pu`, is reported -- on the GENERATOR, `value` the regulated voltage and
 * `limit` the target, both in kV of the regulated bus. AC only: a DC solve has no voltage
 * magnitude to compare.
 */
namespace gen_pv_release_check {

/// A voltage controller with a narrower reactive range never regulates in OpenLoadFlow
/// (`PlausibleValues.MIN_REACTIVE_RANGE`), so a flagged machine below it is not a candidate
constexpr real_type MIN_REACTIVE_RANGE_MVAR = 1.;

/// One flagged PQ machine pinned at a reactive limit, and everything a row needs to check
/// it. Built once per compute() by `build_gen_pv_release_plan`.
struct GenPvReleaseEntry
{
    int gen_id = -1;
    int reg_bus_grid = -1;        ///< the bus it regulates, grid numbering
    int reg_bus_solver = -1;      ///< ... solver numbering (what a row's V is indexed in)
    bool at_min = true;           ///< pinned at min_q (else at max_q)
    real_type target_vm_pu = 0.;  ///< the GRID's target; a row may hand its own (see check)
    real_type vn_kv = 0.;         ///< nominal voltage of the regulated bus: value / limit in kV
    /// LSGrid::set_gen_names, empty if never set
    std::string name;
};

struct GenPvReleasePlan
{
    std::vector<GenPvReleaseEntry> gens;
    bool empty() const { return gens.empty(); }
    void clear() { gens.clear(); }
};

/**
 * Work out, once, which machines can be reported at all, and at which limit each sits.
 *
 * `id_me_to_solver` must describe the labelling the batch solves in (`active_layout()`) --
 * the one a row's `V` is indexed in. `tol_mva` decides "sits at a limit": the reactive
 * setpoint is within it of `min_q` (checked first) or of `max_q`.
 */
inline void build_gen_pv_release_plan(const LSGrid & grid_model,
                                      const SolverBusIdVect & id_me_to_solver,
                                      real_type tol_mva,
                                      GenPvReleasePlan & out)
{
    out.clear();
    const GeneratorContainer & generators = grid_model.get_generators();
    const int nb_gen = generators.nb();
    const std::vector<bool> & status = generators.get_status();
    const std::vector<bool> & can_be_pv = generators.get_can_be_pv();
    const std::vector<std::string> & names = generators.get_names();  // empty if never set
    const Eigen::Ref<const RealVect> vn_kv = grid_model.get_bus_vn_kv();

    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!status[gen_id]) continue;
        if(static_cast<std::size_t>(gen_id) >= can_be_pv.size() || !can_be_pv[gen_id]) continue;
        if(generators.get_voltage_regulator_on(gen_id)) continue;  // regulating: BusQCheck's

        const real_type min_q = generators.get_min_q(gen_id);
        const real_type max_q = generators.get_max_q(gen_id);
        if(!std::isfinite(min_q) || !std::isfinite(max_q)) continue;
        if(max_q - min_q < MIN_REACTIVE_RANGE_MVAR) continue;  // never a voltage controller
        const real_type target_q = generators.get_target_q_mvar(gen_id);
        bool at_min;
        if(std::abs(target_q - min_q) <= tol_mva) at_min = true;
        else if(std::abs(target_q - max_q) <= tol_mva) at_min = false;
        else continue;  // inside its range: not pinned, nothing to release

        const real_type target_vm_pu = generators.get_target_vm_pu(gen_id);
        if(!std::isfinite(target_vm_pu) || target_vm_pu <= 0.) continue;  // no target to hold

        const int reg_bus_grid = generators.get_regulated_bus_id(gen_id);
        if(reg_bus_grid < 0 || reg_bus_grid >= vn_kv.size()) continue;
        const int reg_bus_solver = id_me_to_solver[reg_bus_grid].cast_int();
        if(reg_bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        GenPvReleaseEntry entry;
        entry.gen_id = gen_id;
        entry.reg_bus_grid = reg_bus_grid;
        entry.reg_bus_solver = reg_bus_solver;
        entry.at_min = at_min;
        entry.target_vm_pu = target_vm_pu;
        entry.vn_kv = vn_kv(reg_bus_grid);
        if(static_cast<std::size_t>(gen_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(gen_id)];
        }
        out.gens.push_back(entry);
    }
}

/**
 * Append to `out` one LimitViolation per flagged machine whose regulated voltage sits on
 * the release side of its target, for ONE converged row.
 *
 * `V` is the row's converged complex voltage (solver numbering, pu). `target_vm_of(gen_id)`
 * is the target that machine would hold in THIS row (a sweep may vary it, see
 * BaseBatchSweep::modify_gen_v; the grid's own otherwise) and `is_gen_off(gen_id)` whether
 * the row disconnected it (a generator contingency: nothing to release). `masked_solver_ids`
 * is this row's masked (stranded) solver buses -- sorted, may be nullptr -- whose voltage
 * means nothing.
 */
template<class TargetVmOf, class IsGenOff>
inline void check_gen_pv_release_violations(const GenPvReleasePlan & plan,
                                            const Eigen::Ref<const CplxVect> & V,
                                            real_type tol_vm_pu,
                                            const std::vector<int> * masked_solver_ids,
                                            TargetVmOf target_vm_of,
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

    for(std::size_t k = 0; k < plan.gens.size(); ++k){
        const GenPvReleaseEntry & entry = plan.gens[k];
        if(entry.reg_bus_solver < 0 || entry.reg_bus_solver >= V.size()) continue;
        if(is_masked(entry.reg_bus_solver)) continue;
        if(is_gen_off(entry.gen_id)) continue;
        const real_type target = target_vm_of(entry.gen_id);
        if(!std::isfinite(target) || target <= 0.) continue;
        const real_type vm = std::abs(V(entry.reg_bus_solver));
        if(!std::isfinite(vm)) continue;

        if(entry.at_min && (vm < target - tol_vm_pu)){
            // absorbing as much as it can, and the voltage is still below the target: it
            // absorbs too much for that target, the loop would let it regulate again
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::LOW_VOLTAGE_AT_MIN_Q,
                                         vm * entry.vn_kv, target * entry.vn_kv, entry.name});
        } else if(!entry.at_min && (vm > target + tol_vm_pu)){
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::HIGH_VOLTAGE_AT_MAX_Q,
                                         vm * entry.vn_kv, target * entry.vn_kv, entry.name});
        }
    }
}

}  // namespace gen_pv_release_check
}  // namespace ls2g

#endif  // GENPVRELEASECHECK_H
