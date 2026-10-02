// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef HVDCPCHECK_H
#define HVDCPCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// HvdcPEntry::frozen_dir and the HVDC_AC_EMULATION_RELEASE report exist
#define LS2G_HAS_HVDC_AC_EMULATION_FROZEN 1

namespace ls2g {

/**
 * Post-solve active-power check of the angle-droop ("AC emulation") hvdc lines -- one half
 * of `compute_physical_violations` (the other is BusQCheck.hpp).
 *
 * WHAT THIS ANSWERS. An angle-droop hvdc line transmits what the angle difference asks
 * of it,
 *
 *     p = p0 + k . (theta_1 - theta_2)          (MW, positive 1 -> 2)
 *
 * and lightsim2grid solves that equation as written: the `Hvdc` extension of the
 * Newton-Raphson system (and the DC algorithm) uses `status_droop` as an INPUT -- 0
 * linear, +1 saturated 1->2, -1 saturated 2->1 -- constant across one solve. So in the
 * linear regime nothing stops `p` from coming out beyond `pmax_1to2_mw` /
 * `pmax_2to1_mw`, and the converged row then describes a converter transmitting power it
 * cannot transmit. `LSGrid::set_status_droop_hvdc`'s own documentation says where the
 * saturation belongs -- *"meant to be run between two solves"* -- which is exactly
 * PowSyBl OpenLoadFlow's `HvdcAcEmulationLimits` outer loop, and exactly what this
 * detects. Nothing here re-solves anything, and nothing is clamped: a row is reported,
 * not fixed.
 *
 * WHY IT IS A PHYSICAL LIMIT (ViolationCategory::PHYSICAL, like a machine's reactive
 * capability and unlike a branch's current rating). A line above its thermal rating is a
 * state the grid reaches and should not sit in; a converter beyond its maximum power is a
 * state it does not reach at all -- its own control saturates first, which is what
 * `status_droop = +/-1` models. So the converged solution is not the grid's answer to the
 * question asked, it is the answer to a question with a constraint left out.
 *
 * WHAT IS COMPARED. The SENDING end, in each direction, against that direction's own
 * maximum -- `pmax_1to2_mw` bounds `p` and `pmax_2to1_mw` bounds `-p`, so the whole
 * condition is `-pmax_2to1 <= p <= pmax_1to2`. That is the same quantity the solver's own
 * saturation pins (`HvdcDroopSolverData::flows_pu` sets the controller-side flow to
 * `pmax12` / `pmax21` when `status_droop` says saturated), so a row reported here is a row
 * whose flow the saturated regime would have moved. The receiving end carries less (the
 * converter and dc-line losses) and has no separate limit of its own.
 *
 * WHICH LINES ARE LOOKED AT. Those whose droop is active (`is_droop_active`: enabled and
 * the line in service) AND still in the LINEAR regime (`status_droop == 0`). A line the
 * caller has already saturated is skipped: its flow is pinned at the very limit this
 * would compare against, so there is nothing left to detect -- the outer loop, in effect,
 * already fired. A line the CALLER says an outer loop froze at its limit
 * (LSGrid::set_hvdc_ac_emulation_frozen: droop off, set-point at the limit, droop parameters
 * kept) is checked for the other half of that loop instead: the flow its droop would ask for
 * from this row's angles, below that limit by more than the tolerance, says the loop would
 * have left it in AC emulation -- reported as HVDC_AC_EMULATION_RELEASE, `side` the direction
 * it is frozen in. (The same test as the PQ -> PV release of GenPvReleaseCheck.hpp: the
 * frozen state's own angles, not the ones the released line would lead to.)
 */
namespace hvdc_p_check {

/// One angle-droop hvdc line in the linear regime, and everything a row needs to check
/// it. Built once per compute() by `build_hvdc_p_plan`; nothing here varies from row to
/// row (a contingency edits `Ybus`, never the droop configuration).
struct HvdcPEntry
{
    int hvdc_id = -1;
    int bus1_solver = -1;
    int bus2_solver = -1;
    real_type p0_mw = 0.;            ///< the droop's offset
    real_type k_mw_per_rad = 0.;     ///< ... and its slope
    real_type pmax_1to2_mw = 0.;
    real_type pmax_2to1_mw = 0.;
    /// 0: a line in its linear regime (the check above); +1 / -1: a line an outer loop froze
    /// at its limit in the direction 1 -> 2 / 2 -> 1 (LSGrid::set_hvdc_ac_emulation_frozen),
    /// checked for its release instead
    int frozen_dir = 0;
    std::string name;                ///< LSGrid::set_dcline_names, empty if never set
};

struct HvdcPPlan
{
    std::vector<HvdcPEntry> lines;

    bool empty() const { return lines.empty(); }
    void clear() { lines.clear(); }
};

/**
 * Work out, once, which angle-droop hvdc lines can be reported at all.
 *
 * `id_me_to_solver` must describe the labelling the batch solves in (`active_layout()`),
 * the same one `BaseAlgo::get_Va()` is indexed in -- reading the angles of one labelling
 * through another is a wrong answer, not an error, so the two come from the same place.
 */
inline void build_hvdc_p_plan(const LSGrid & grid_model,
                              const SolverBusIdVect & id_me_to_solver,
                              HvdcPPlan & out)
{
    out.clear();
    const HvdcLineContainer & hvdcs = grid_model.get_dclines();
    const int nb_hvdc = hvdcs.nb();
    if(nb_hvdc == 0) return;
    const std::vector<std::string> & names = hvdcs.get_names();  // empty if never set

    const std::vector<bool> & status_global = hvdcs.get_status_global();
    const std::vector<bool> & frozen = hvdcs.get_ac_emulation_frozen();
    for(int hvdc_id = 0; hvdc_id < nb_hvdc; ++hvdc_id){
        // a line an outer loop froze at its limit (droop off, set-point there): its release
        int frozen_dir = 0;
        if(static_cast<std::size_t>(hvdc_id) < frozen.size() && frozen[hvdc_id] &&
           !hvdcs.is_droop_active(hvdc_id) && status_global[hvdc_id]){
            frozen_dir = (hvdcs.get_converters_mode(hvdc_id) == HvdcLineContainer::ConvertersMode::SIDE_1_RECTIFIER)
                         ? 1 : -1;
        } else {
            if(!hvdcs.is_droop_active(hvdc_id)) continue;        // no droop, or out of service
            if(hvdcs.get_status_droop(hvdc_id) != 0) continue;   // already saturated by the caller
        }
        // A droop-active line always has both converters connected (LSGrid opens one only
        // through deactivate_dcline_sideX / disconnect_if_not_in_main_component, both of
        // which disable_droop, and fill_hvdc_droop_solver_data throws otherwise). Checked
        // anyway rather than angle-deriving a flow across an open converter.
        if(!hvdcs.get_connected_side_1(hvdc_id) || !hvdcs.get_connected_side_2(hvdc_id)) continue;

        const int bus1_me = hvdcs.get_bus_side_1(hvdc_id).cast_int();
        const int bus2_me = hvdcs.get_bus_side_2(hvdc_id).cast_int();
        if(bus1_me == BaseConstants::_deactivated_bus_id) continue;
        if(bus2_me == BaseConstants::_deactivated_bus_id) continue;
        const int bus1_solver = id_me_to_solver[bus1_me].cast_int();
        const int bus2_solver = id_me_to_solver[bus2_me].cast_int();
        if(bus1_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system
        if(bus2_solver == BaseConstants::_deactivated_bus_id) continue;

        HvdcPEntry entry;
        entry.hvdc_id = hvdc_id;
        entry.bus1_solver = bus1_solver;
        entry.bus2_solver = bus2_solver;
        // MW and MW/rad, exactly as LSGrid::fill_hvdc_droop_solver_data reads them before
        // dividing by sn_mva -- so `p` comes out in MW and the two pmax are already there
        entry.p0_mw = hvdcs.get_droop_p0_mw(hvdc_id);
        entry.k_mw_per_rad = hvdcs.get_droop_k_mw_per_rad(hvdc_id);
        entry.pmax_1to2_mw = hvdcs.get_pmax_1to2_mw(hvdc_id);
        entry.pmax_2to1_mw = hvdcs.get_pmax_2to1_mw(hvdc_id);
        entry.frozen_dir = frozen_dir;
        if(static_cast<std::size_t>(hvdc_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(hvdc_id)];
        }
        out.lines.push_back(entry);
    }
}

/**
 * Append to `out` one LimitViolation per angle-droop hvdc line whose active power left
 * what its converters can transmit, for ONE converged row.
 *
 * `Va` is the row's converged bus angles (radians, solver numbering -- `BaseAlgo::get_Va`,
 * which every algorithm fills, AC and DC alike: this check needs nothing else from the
 * solve). `masked_solver_ids` is this row's masked (stranded) solver buses -- sorted, may
 * be nullptr -- whose angle means nothing.
 *
 * `side` on the reported violation says which direction, and therefore which limit: 1 for
 * a flow leaving side 1 (1 -> 2, against `pmax_1to2_mw`), 2 for the other way. `value` is
 * that flow, always positive, in MW.
 */
inline void check_hvdc_p_violations(const HvdcPPlan & plan,
                                    const Eigen::Ref<const RealVect> & Va,
                                    real_type tol_mw,
                                    const std::vector<int> * masked_solver_ids,
                                    std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;
    if(Va.size() == 0) return;

    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();
    auto is_masked = [&](int bus){
        // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
        return has_masked && std::binary_search(masked_solver_ids->begin(),
                                                masked_solver_ids->end(), bus);
    };

    for(std::size_t k = 0; k < plan.lines.size(); ++k){
        const HvdcPEntry & entry = plan.lines[k];
        if(entry.bus1_solver >= static_cast<int>(Va.size())) continue;
        if(entry.bus2_solver >= static_cast<int>(Va.size())) continue;
        if(is_masked(entry.bus1_solver) || is_masked(entry.bus2_solver)) continue;

        const real_type p_mw = entry.p0_mw + entry.k_mw_per_rad *
                               (Va(entry.bus1_solver) - Va(entry.bus2_solver));
        if(!std::isfinite(p_mw)) continue;

        if(entry.frozen_dir != 0){
            // saturated (frozen) at its limit in that direction: OpenLoadFlow's
            // AbstractHvdcAcEmulationLimitsOuterLoop.checkSaturationMode. Its droop flow,
            // strictly inside the limit of the direction it now flows in, releases it; a flow
            // reversed beyond the other direction's limit saturates it on that side instead.
            const bool inside = (p_mw > 0. && p_mw < entry.pmax_1to2_mw - tol_mw) ||
                                (p_mw < 0. && -p_mw < entry.pmax_2to1_mw - tol_mw);
            if(inside){
                const int side = entry.frozen_dir > 0 ? 1 : 2;
                const real_type flow = entry.frozen_dir > 0 ? p_mw : -p_mw;
                const real_type limit = entry.frozen_dir > 0 ? entry.pmax_1to2_mw : entry.pmax_2to1_mw;
                out.push_back(LimitViolation{ViolationElementType::HVDC, entry.hvdc_id, side,
                                             LimitViolationType::HVDC_AC_EMULATION_RELEASE, flow,
                                             limit, entry.name});
            } else if(entry.frozen_dir > 0 && p_mw < 0.){
                out.push_back(LimitViolation{ViolationElementType::HVDC, entry.hvdc_id, 2,
                                             LimitViolationType::HIGH_P, -p_mw,
                                             entry.pmax_2to1_mw, entry.name});
            } else if(entry.frozen_dir < 0 && p_mw > 0.){
                out.push_back(LimitViolation{ViolationElementType::HVDC, entry.hvdc_id, 1,
                                             LimitViolationType::HIGH_P, p_mw,
                                             entry.pmax_1to2_mw, entry.name});
            }
            continue;
        }
        if(std::isfinite(entry.pmax_1to2_mw) && (p_mw > entry.pmax_1to2_mw + tol_mw)){
            out.push_back(LimitViolation{ViolationElementType::HVDC, entry.hvdc_id, 1,
                                         LimitViolationType::HIGH_P, p_mw,
                                         entry.pmax_1to2_mw, entry.name});
        } else if(std::isfinite(entry.pmax_2to1_mw) && (-p_mw > entry.pmax_2to1_mw + tol_mw)){
            out.push_back(LimitViolation{ViolationElementType::HVDC, entry.hvdc_id, 2,
                                         LimitViolationType::HIGH_P, -p_mw,
                                         entry.pmax_2to1_mw, entry.name});
        }
    }
}

}  // namespace hvdc_p_check
}  // namespace ls2g

#endif  // HVDCPCHECK_H
