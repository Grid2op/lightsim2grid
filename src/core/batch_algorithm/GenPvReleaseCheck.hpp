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
 * SVCS TOO. An SVC an outer loop froze at a reactive limit (a voltage-mode SVC turned fixed-Q
 * at the edge of its susceptance range -- flagged with LSGrid::set_svc_can_be_pv, which the
 * converter fills from what `bake_outer_loops` froze) is the same statement: absorbing all it
 * can while the bus it regulates is still below its target (or producing all it can while
 * that bus is above it), OpenLoadFlow would let it regulate again. Reported on the SVC with
 * the same two types. Its range is a susceptance, worth `b * V^2` MVAr: which limit it sits at
 * is the nearer end of that range at its target voltage (a frozen SVC keeps the output it had
 * while holding its target, a hair inside the limit at most). No batch axis varies an SVC's
 * target, and no contingency disconnects one: a row checks it against the grid's target.
 *
 * WHAT IS CHECKED. A flagged PQ generator, pinned at the NEARER of its two reactive limits
 * (decided once per plan): the flag already says an outer loop froze it at a limit, and a
 * bake freezes it at the output it had, which may sit a hair inside that limit -- so no
 * tolerance on the setpoint gates it. With a reactive range of at
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
    /// GENERATOR, or SVC for a flagged frozen SVC (then `gen_id` is the svc id)
    ViolationElementType el_type = ViolationElementType::GENERATOR;
    int gen_id = -1;
    int reg_bus_grid = -1;        ///< the bus it regulates, grid numbering
    int reg_bus_solver = -1;      ///< ... solver numbering (what a row's V is indexed in)
    int gen_bus_solver = -1;      ///< the machine's own bus, solver numbering
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
 * the one a row's `V` is indexed in. A flagged machine sits at the nearer of its limits
 * (`min_q` on a tie): the flag says it was frozen at one.
 */
inline void build_gen_pv_release_plan(const LSGrid & grid_model,
                                      const SolverBusIdVect & id_me_to_solver,
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
        // the nearer limit: the flag says an outer loop froze it at one, and a bake keeps the
        // output it had there, which may sit a hair inside it
        const bool at_min = std::abs(target_q - min_q) <= std::abs(target_q - max_q);

        const real_type target_vm_pu = generators.get_target_vm_pu(gen_id);
        if(!std::isfinite(target_vm_pu) || target_vm_pu <= 0.) continue;  // no target to hold

        const int reg_bus_grid = generators.get_regulated_bus_id(gen_id);
        if(reg_bus_grid < 0 || reg_bus_grid >= vn_kv.size()) continue;
        const int reg_bus_solver = id_me_to_solver[reg_bus_grid].cast_int();
        if(reg_bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system
        const int gen_bus_grid = generators.get_bus_id()(gen_id).cast_int();
        if(gen_bus_grid < 0 || gen_bus_grid >= vn_kv.size()) continue;
        const int gen_bus_solver = id_me_to_solver[gen_bus_grid].cast_int();
        if(gen_bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        GenPvReleaseEntry entry;
        entry.gen_id = gen_id;
        entry.reg_bus_grid = reg_bus_grid;
        entry.reg_bus_solver = reg_bus_solver;
        entry.gen_bus_solver = gen_bus_solver;
        entry.at_min = at_min;
        entry.target_vm_pu = target_vm_pu;
        entry.vn_kv = vn_kv(reg_bus_grid);
        if(static_cast<std::size_t>(gen_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(gen_id)];
        }
        out.gens.push_back(entry);
    }

    // the flagged SVCs an outer loop froze at a limit (LSGrid::set_svc_can_be_pv)
    const SvcContainer & svcs = grid_model.get_svcs();
    const int nb_svc = svcs.nb();
    const std::vector<bool> & svc_status = svcs.get_status();
    const std::vector<bool> & svc_can_be_pv = svcs.get_can_be_pv();
    const std::vector<std::string> & svc_names = svcs.get_names();
    const real_type sn_mva = grid_model.get_sn_mva();
    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!svc_status[svc_id]) continue;
        if(static_cast<std::size_t>(svc_id) >= svc_can_be_pv.size() || !svc_can_be_pv[svc_id]) continue;
        // regulating (already a controller) or OFF: nothing to release
        if(svcs.get_regulation_mode(svc_id) != SvcContainer::RegulationMode::REACTIVE_POWER) continue;

        const real_type target_vm_pu = svcs.get_target_vm_pu(svc_id);
        if(!std::isfinite(target_vm_pu) || target_vm_pu <= 0.) continue;  // no target to hold
        // its reactive range at its target voltage, MVAr (b in pu of sn_mva)
        const real_type v2 = target_vm_pu * target_vm_pu * sn_mva;
        const real_type min_q = svcs.get_b_min(svc_id) * v2;
        const real_type max_q = svcs.get_b_max(svc_id) * v2;
        if(!std::isfinite(min_q) || !std::isfinite(max_q)) continue;
        if(max_q - min_q < MIN_REACTIVE_RANGE_MVAR) continue;  // never a voltage controller
        const real_type target_q = svcs.get_target_q()(svc_id);  // generator convention, MVAr
        const bool at_min = std::abs(target_q - min_q) <= std::abs(target_q - max_q);

        const int reg_bus_grid = svcs.get_regulated_bus_id(svc_id);
        if(reg_bus_grid < 0 || reg_bus_grid >= vn_kv.size()) continue;
        const int reg_bus_solver = id_me_to_solver[reg_bus_grid].cast_int();
        if(reg_bus_solver == BaseConstants::_deactivated_bus_id) continue;
        const int svc_bus_grid = svcs.get_bus_id()(svc_id).cast_int();
        if(svc_bus_grid < 0 || svc_bus_grid >= vn_kv.size()) continue;
        const int svc_bus_solver = id_me_to_solver[svc_bus_grid].cast_int();
        if(svc_bus_solver == BaseConstants::_deactivated_bus_id) continue;

        GenPvReleaseEntry entry;
        entry.el_type = ViolationElementType::SVC;
        entry.gen_id = svc_id;
        entry.reg_bus_grid = reg_bus_grid;
        entry.reg_bus_solver = reg_bus_solver;
        entry.gen_bus_solver = svc_bus_solver;
        entry.at_min = at_min;
        entry.target_vm_pu = target_vm_pu;
        entry.vn_kv = vn_kv(reg_bus_grid);
        if(static_cast<std::size_t>(svc_id) < svc_names.size()){
            entry.name = svc_names[static_cast<std::size_t>(svc_id)];
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
 * the row disconnected it (a generator contingency: nothing to release). Both are asked
 * about the GENERATOR entries only: an SVC entry keeps the grid's target and is never
 * disconnected by a row. `masked_solver_ids`
 * is this row's masked (stranded) solver buses -- sorted, may be nullptr -- whose voltage
 * means nothing: a machine is skipped when its regulated bus OR its own bus is masked.
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
        // a machine stranded outside the main component (its own bus masked) is out of
        // the solve, as when the contingency disconnects it: it releases nothing, even
        // when the bus it regulates stays in the main component
        if(is_masked(entry.gen_bus_solver)) continue;
        const bool is_gen = (entry.el_type == ViolationElementType::GENERATOR);
        if(is_gen && is_gen_off(entry.gen_id)) continue;
        const real_type target = is_gen ? target_vm_of(entry.gen_id) : entry.target_vm_pu;
        if(!std::isfinite(target) || target <= 0.) continue;
        const real_type vm = std::abs(V(entry.reg_bus_solver));
        if(!std::isfinite(vm)) continue;

        if(entry.at_min && (vm < target - tol_vm_pu)){
            // absorbing as much as it can, and the voltage is still below the target: it
            // absorbs too much for that target, the loop would let it regulate again
            out.push_back(LimitViolation{entry.el_type, entry.gen_id, 0,
                                         LimitViolationType::LOW_VOLTAGE_AT_MIN_Q,
                                         vm * entry.vn_kv, target * entry.vn_kv, entry.name});
        } else if(!entry.at_min && (vm > target + tol_vm_pu)){
            out.push_back(LimitViolation{entry.el_type, entry.gen_id, 0,
                                         LimitViolationType::HIGH_VOLTAGE_AT_MAX_Q,
                                         vm * entry.vn_kv, target * entry.vn_kv, entry.name});
        }
    }
}

}  // namespace gen_pv_release_check
}  // namespace ls2g

#endif  // GENPVRELEASECHECK_H
