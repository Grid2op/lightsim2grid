// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SVCSTANDBYCHECK_H
#define SVCSTANDBYCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace ls2g {

/**
 * Post-solve check of the idle SVCs carrying a standby automaton: would that automaton
 * switch them to voltage control? One of the physical checks of
 * `compute_physical_violations` (with BusQCheck.hpp, GenPvReleaseCheck.hpp, HvdcPCheck.hpp
 * and GenPCheck.hpp).
 *
 * WHAT THIS ANSWERS. A voltage-regulating SVC whose powsybl `standbyAutomaton` extension
 * says `standby` does not regulate at first: OpenLoadFlow's MonitoringVoltageOuterLoop
 * keeps it idle (a fixed `b0` shunt) while the voltage of the bus it controls stays inside
 * [low_voltage_threshold, high_voltage_threshold], and switches it to voltage control as
 * soon as that voltage leaves them. An outer-loop-free model of such an SVC (what
 * `bake_outer_loops` produces: a fixed-Q SVC) cannot do that switch, so a solution whose
 * regulated voltage is outside the thresholds assumes an SVC the loop would not leave idle
 * -- the same kind of statement as the PQ -> PV release of GenPvReleaseCheck.hpp
 * (ViolationCategory::CONTROL: an automaton not triggered), and nothing here re-solves or switches anything: a row
 * is reported, not fixed.
 *
 * ONLY THE SVCS A CALLER FLAGGED. After a bake, an idle standby SVC is a plain fixed-Q
 * SVC, and OpenLoadFlow ignores the automaton of an SVC that does not regulate voltage in
 * the first place: lightsim2grid cannot tell the two apart. The flag and the thresholds
 * are `SvcContainer::set_standby` (LSGrid::set_svc_standby), which the pypowsybl converter
 * fills for the SVCs `bake_outer_loops` left idle; everything else is skipped, and so is a
 * flagged SVC that regulates voltage (it is already switched on).
 *
 * WHAT IS CHECKED. Per row, the magnitude of the bus a flagged, connected, non-regulating
 * SVC regulates (the one OpenLoadFlow monitors) is compared with the thresholds: below the
 * low one or above the high one by more than `tol_vm_pu` is reported -- on the SVC,
 * `value` the regulated voltage and `limit` the threshold, both in kV of the regulated
 * bus. AC only: a DC solve has no voltage magnitude to compare.
 */
namespace svc_standby_check {

/// One flagged idle SVC and everything a row needs to check it. Built once per compute()
/// by `build_svc_standby_plan`.
struct SvcStandbyEntry
{
    int svc_id = -1;
    int reg_bus_grid = -1;        ///< the bus it regulates, grid numbering
    int reg_bus_solver = -1;      ///< ... solver numbering (what a row's V is indexed in)
    int svc_bus_solver = -1;      ///< the SVC's own bus, solver numbering
    real_type low_vm_pu = 0.;     ///< the automaton's low threshold, pu of the regulated bus
    real_type high_vm_pu = 0.;    ///< the automaton's high threshold, pu of the regulated bus
    real_type vn_kv = 0.;         ///< nominal voltage of the regulated bus: value / limit in kV
    /// LSGrid::set_svc_names, empty if never set
    std::string name;
};

struct SvcStandbyPlan
{
    std::vector<SvcStandbyEntry> svcs;
    bool empty() const { return svcs.empty(); }
    void clear() { svcs.clear(); }
};

/**
 * Work out, once, which SVCs can be reported at all.
 *
 * `id_me_to_solver` must describe the labelling the batch solves in (`active_layout()`) --
 * the one a row's `V` is indexed in.
 */
inline void build_svc_standby_plan(const LSGrid & grid_model,
                                   const SolverBusIdVect & id_me_to_solver,
                                   SvcStandbyPlan & out)
{
    out.clear();
    const SvcContainer & svcs = grid_model.get_svcs();
    const int nb_svc = svcs.nb();
    const std::vector<bool> & status = svcs.get_status();
    const std::vector<bool> & standby = svcs.get_standby();
    const std::vector<std::string> & names = svcs.get_names();  // empty if never set
    const Eigen::Ref<const RealVect> vn_kv = grid_model.get_bus_vn_kv();

    for(int svc_id = 0; svc_id < nb_svc; ++svc_id){
        if(!status[svc_id]) continue;
        if(static_cast<std::size_t>(svc_id) >= standby.size() || !standby[svc_id]) continue;
        // regulating: already switched on, nothing left to switch
        if(svcs.get_regulation_mode(svc_id) == SvcContainer::RegulationMode::VOLTAGE) continue;

        const real_type low_vm_pu = svcs.get_standby_low_vm_pu(svc_id);
        const real_type high_vm_pu = svcs.get_standby_high_vm_pu(svc_id);
        if(!std::isfinite(low_vm_pu) || !std::isfinite(high_vm_pu)) continue;

        const int reg_bus_grid = svcs.get_regulated_bus_id(svc_id);
        if(reg_bus_grid < 0 || reg_bus_grid >= vn_kv.size()) continue;
        const int reg_bus_solver = id_me_to_solver[reg_bus_grid].cast_int();
        if(reg_bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system
        const int svc_bus_grid = svcs.get_bus_id()(svc_id).cast_int();
        if(svc_bus_grid < 0 || svc_bus_grid >= vn_kv.size()) continue;
        const int svc_bus_solver = id_me_to_solver[svc_bus_grid].cast_int();
        if(svc_bus_solver == BaseConstants::_deactivated_bus_id) continue;  // not in the solved system

        SvcStandbyEntry entry;
        entry.svc_id = svc_id;
        entry.reg_bus_grid = reg_bus_grid;
        entry.reg_bus_solver = reg_bus_solver;
        entry.svc_bus_solver = svc_bus_solver;
        entry.low_vm_pu = low_vm_pu;
        entry.high_vm_pu = high_vm_pu;
        entry.vn_kv = vn_kv(reg_bus_grid);
        if(static_cast<std::size_t>(svc_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(svc_id)];
        }
        out.svcs.push_back(entry);
    }
}

/**
 * Append to `out` one LimitViolation per flagged SVC whose regulated voltage is outside
 * its automaton's thresholds, for ONE converged row.
 *
 * `V` is the row's converged complex voltage (solver numbering, pu). `masked_solver_ids`
 * is this row's masked (stranded) solver buses -- sorted, may be nullptr -- whose voltage
 * means nothing: an SVC is skipped when its regulated bus OR its own bus is masked (an SVC
 * stranded outside the main component is out of the solve, as when it is disconnected).
 */
inline void check_svc_standby_violations(const SvcStandbyPlan & plan,
                                         const Eigen::Ref<const CplxVect> & V,
                                         real_type tol_vm_pu,
                                         const std::vector<int> * masked_solver_ids,
                                         std::vector<LimitViolation> & out)
{
    if(plan.empty()) return;
    const bool has_masked = (masked_solver_ids != nullptr) && !masked_solver_ids->empty();
    // every _li_masked entry is sorted (see BaseBatchSweep::_prepare_connectivity)
    auto is_masked = [&](int bus){
        return has_masked && std::binary_search(masked_solver_ids->begin(),
                                                masked_solver_ids->end(), bus);
    };

    for(std::size_t k = 0; k < plan.svcs.size(); ++k){
        const SvcStandbyEntry & entry = plan.svcs[k];
        if(entry.reg_bus_solver < 0 || entry.reg_bus_solver >= V.size()) continue;
        if(is_masked(entry.reg_bus_solver)) continue;
        if(is_masked(entry.svc_bus_solver)) continue;
        const real_type vm = std::abs(V(entry.reg_bus_solver));
        if(!std::isfinite(vm)) continue;

        if(vm < entry.low_vm_pu - tol_vm_pu){
            out.push_back(LimitViolation{ViolationElementType::SVC, entry.svc_id, 0,
                                         LimitViolationType::LOW_VOLTAGE_SVC_STANDBY,
                                         vm * entry.vn_kv, entry.low_vm_pu * entry.vn_kv, entry.name});
        } else if(vm > entry.high_vm_pu + tol_vm_pu){
            out.push_back(LimitViolation{ViolationElementType::SVC, entry.svc_id, 0,
                                         LimitViolationType::HIGH_VOLTAGE_SVC_STANDBY,
                                         vm * entry.vn_kv, entry.high_vm_pu * entry.vn_kv, entry.name});
        }
    }
}

}  // namespace svc_standby_check
}  // namespace ls2g

#endif  // SVCSTANDBYCHECK_H
