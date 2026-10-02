// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef REMOTEVOLTAGECONTROLCHECK_H
#define REMOTEVOLTAGECONTROLCHECK_H

#include "LSGrid.hpp"
#include "LimitViolation.hpp"
#include "VoltageControlData.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace ls2g {

/**
 * Post-solve check of the generators regulating a REMOTE bus: does holding that remote
 * target take their own bus to a voltage nobody would let a controller sit at? One of the
 * physical checks of `compute_physical_violations` (with BusQCheck.hpp,
 * GenPvReleaseCheck.hpp, SvcStandbyCheck.hpp, HvdcPCheck.hpp and GenPCheck.hpp).
 *
 * WHAT THIS ANSWERS. A generator can hold the voltage of a bus other than its own (through
 * its step-up transformer, typically). Holding a remote target may take a lot of reactive
 * power, and with it the generator's own bus far from anything realistic, while its
 * reactive output stays inside its limits -- so BusQCheck has nothing to say. OpenLoadFlow's
 * ReactiveLimits loop, in its "robust" remote voltage control mode (the default,
 * `voltageRemoteControlRobustMode`), switches such a controller bus to PQ at its target
 * reactive power as soon as its voltage leaves
 * [minRealisticVoltage * 1.02, maxRealisticVoltage / 1.02] (pu of its nominal voltage).
 * An outer-loop-free solve cannot do that switch, so a solution whose remote controller sits
 * outside that range assumes a control the loop would not leave in place, at a voltage a
 * controller's own bus is not meant to reach (ViolationCategory::PHYSICAL), and nothing
 * here re-solves or switches anything: a row is reported, not fixed.
 *
 * WHICH GENERATORS. The controllers of the voltage-control plan the solve used whose own bus
 * is not the bus their group regulates: the generators that actually hold a remote bus in
 * that solve. A held controller (LSGrid::set_hold_frozen_regulators) does not regulate and is
 * skipped, and so is any other kind of controller (OpenLoadFlow's robust mode is about
 * generator voltage control). Nothing is checked unless the caller set the range
 * (LSGrid::set_remote_voltage_control_vm_range).
 *
 * WHAT IS CHECKED. Per row, the magnitude of the generator's own bus: below the low bound or
 * above the high one by more than `tol_vm_pu` is reported -- on the GENERATOR, `value` that
 * voltage and `limit` the bound, both in kV of its own bus. AC only: a DC solve has no
 * voltage magnitude to compare.
 */
namespace remote_voltage_control_check {

/// One generator regulating a remote bus, and everything a row needs to check it. Built once
/// per compute() by `build_remote_voltage_control_plan`.
struct RemoteVoltageControlEntry
{
    int gen_id = -1;
    int gen_bus_solver = -1;      ///< the generator's own bus, solver numbering (what is checked)
    int reg_bus_solver = -1;      ///< the remote bus it regulates, solver numbering
    real_type vn_kv = 0.;         ///< nominal voltage of its own bus: value / limit in kV
    /// LSGrid::set_gen_names, empty if never set
    std::string name;
};

struct RemoteVoltageControlPlan
{
    std::vector<RemoteVoltageControlEntry> gens;
    real_type min_vm_pu = 0.;  ///< LSGrid::get_remote_voltage_control_min_vm_pu
    real_type max_vm_pu = 0.;  ///< LSGrid::get_remote_voltage_control_max_vm_pu
    bool empty() const { return gens.empty(); }
    void clear() { gens.clear(); }
};

/**
 * Work out, once, which generators regulate a remote bus in the solve.
 *
 * `id_me_to_solver` and `controllers` must describe the labelling and the voltage-control
 * plan the batch solves with (`active_layout()`) -- the ones a row's `V` is indexed in.
 */
inline void build_remote_voltage_control_plan(const LSGrid & grid_model,
                                              const SolverBusIdVect & id_me_to_solver,
                                              const VoltageControlSolverData & controllers,
                                              RemoteVoltageControlPlan & out)
{
    out.clear();
    out.min_vm_pu = grid_model.get_remote_voltage_control_min_vm_pu();
    out.max_vm_pu = grid_model.get_remote_voltage_control_max_vm_pu();
    // not set: nothing to check
    if(!std::isfinite(out.min_vm_pu) && !std::isfinite(out.max_vm_pu)) return;

    const GeneratorContainer & generators = grid_model.get_generators();
    const std::vector<std::string> & names = generators.get_names();  // empty if never set
    const Eigen::Ref<const RealVect> vn_kv = grid_model.get_bus_vn_kv();
    const int nb_gen = generators.nb();
    for(int c = 0; c < controllers.n_controllers(); ++c){
        if(controllers.kind(c) != VoltageControlSolverData::GEN) continue;
        if(controllers.is_held(c)) continue;  // does not regulate
        const int gen_bus_solver = controllers.bus(c);
        const int reg_bus_solver = controllers.reg_bus(controllers.group(c));
        if(gen_bus_solver == reg_bus_solver) continue;  // a local controller
        const int gen_id = controllers.elem_id(c);
        if(gen_id < 0 || gen_id >= nb_gen) continue;
        const int gen_bus_grid = generators.get_bus_id()(gen_id).cast_int();
        if(gen_bus_grid < 0 || gen_bus_grid >= vn_kv.size()) continue;
        if(id_me_to_solver[gen_bus_grid].cast_int() != gen_bus_solver) continue;  // not this labelling

        RemoteVoltageControlEntry entry;
        entry.gen_id = gen_id;
        entry.gen_bus_solver = gen_bus_solver;
        entry.reg_bus_solver = reg_bus_solver;
        entry.vn_kv = vn_kv(gen_bus_grid);
        if(static_cast<std::size_t>(gen_id) < names.size()){
            entry.name = names[static_cast<std::size_t>(gen_id)];
        }
        out.gens.push_back(entry);
    }
}

/**
 * Append to `out` one LimitViolation per remote controller whose own bus is outside the
 * realistic range, for ONE converged row.
 *
 * `V` is the row's converged complex voltage (solver numbering, pu). `is_gen_off(gen_id)`
 * says whether the row disconnected it (a generator contingency: it regulates nothing).
 * `masked_solver_ids` is this row's masked (stranded) solver buses -- sorted, may be
 * nullptr -- whose voltage means nothing: a generator is skipped when its own bus OR the
 * bus it regulates is masked (it holds nothing then).
 */
template<class IsGenOff>
inline void check_remote_voltage_control_violations(const RemoteVoltageControlPlan & plan,
                                                    const Eigen::Ref<const CplxVect> & V,
                                                    real_type tol_vm_pu,
                                                    const std::vector<int> * masked_solver_ids,
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
        const RemoteVoltageControlEntry & entry = plan.gens[k];
        if(entry.gen_bus_solver < 0 || entry.gen_bus_solver >= V.size()) continue;
        if(is_masked(entry.gen_bus_solver)) continue;
        if(is_masked(entry.reg_bus_solver)) continue;
        if(is_gen_off(entry.gen_id)) continue;
        const real_type vm = std::abs(V(entry.gen_bus_solver));
        if(!std::isfinite(vm)) continue;

        if(std::isfinite(plan.min_vm_pu) && (vm < plan.min_vm_pu - tol_vm_pu)){
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::LOW_VOLTAGE_REMOTE_CONTROL,
                                         vm * entry.vn_kv, plan.min_vm_pu * entry.vn_kv, entry.name});
        } else if(std::isfinite(plan.max_vm_pu) && (vm > plan.max_vm_pu + tol_vm_pu)){
            out.push_back(LimitViolation{ViolationElementType::GENERATOR, entry.gen_id, 0,
                                         LimitViolationType::HIGH_VOLTAGE_REMOTE_CONTROL,
                                         vm * entry.vn_kv, plan.max_vm_pu * entry.vn_kv, entry.name});
        }
    }
}

}  // namespace remote_voltage_control_check
}  // namespace ls2g

#endif  // REMOTEVOLTAGECONTROLCHECK_H
