// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "PowerFlowAdjoint.hpp"

#include <algorithm>
#include <map>

#include "LSGrid.hpp"

namespace ls2g {
namespace power_flow_adjoint {

namespace {

bool contains(const std::vector<int> & sorted, int value)
{
    return std::binary_search(sorted.begin(), sorted.end(), value);
}

}  // namespace

GenVLayout gen_v_layout(const LSGrid & grid, const SolverBusLayout & layout, const VoltageControlPlan & plan)
{
    const Eigen::Index nb_gen = static_cast<Eigen::Index>(grid.get_generators_as_data().nb());
    GenVLayout res;
    res.target_bus = IntVect::Constant(nb_gen, -1);
    res.vc_group = IntVect::Constant(nb_gen, -1);
    res.share = RealVect::Zero(nb_gen);

    std::map<int, std::vector<int> > gens_of_bus;   // solver bus -> the generators writing its |V|
    grid.get_generators().vm_writers_by_bus(layout.id_me_to_solver, gens_of_bus);
    std::map<int, real_type> pinned;                // solver buses an SVC / a station holds
    grid.get_svcs().vm_targets_by_bus(layout.id_me_to_solver, pinned);
    grid.get_dclines().vm_targets_by_bus(layout.id_me_to_solver, pinned);

    // the buses a voltage-control group holds: a generator regulating one fixes no |V|, its
    // set-point is the group's v_set -- unless the group holds an SVC or a station too, whose
    // set-point no input moves (-1: no derivative for anyone there)
    const VoltageControlSolverData & ctrl = plan.controllers();
    std::map<int, int> group_of_bus;
    for(int grp = 0; grp < ctrl.n_groups(); ++grp){
        bool gens_only = true;
        for(int j = ctrl.grp_start(grp); j < ctrl.grp_start(grp) + ctrl.grp_count(grp); ++j){
            if(ctrl.kind(j) != VoltageControlSolverData::GEN) gens_only = false;
        }
        group_of_bus[ctrl.reg_bus(grp)] = gens_only ? grp : -1;
    }

    const Eigen::Index nb_bus_solver = static_cast<Eigen::Index>(layout.id_solver_to_me.size());
    std::map<int, int> count_bus, count_grp;
    for(const auto & bus_and_gens : gens_of_bus){
        const int b = bus_and_gens.first;
        if(b < 0 || b >= nb_bus_solver) continue;
        const auto grp = group_of_bus.find(b);
        if(grp != group_of_bus.end()){
            if(grp->second < 0) continue;
            for(int g : bus_and_gens.second) res.vc_group[g] = grp->second;
            count_grp[grp->second] += static_cast<int>(bus_and_gens.second.size());
            continue;
        }
        if(pinned.find(b) != pinned.end()) continue;   // held by something else
        for(int g : bus_and_gens.second) res.target_bus[g] = b;
        count_bus[b] += static_cast<int>(bus_and_gens.second.size());
        res.vm_buses.push_back(b);                     // in increasing order (a std::map)
    }
    for(Eigen::Index g = 0; g < nb_gen; ++g){
        if(res.target_bus[g] >= 0) res.share[g] = 1. / static_cast<real_type>(count_bus[res.target_bus[g]]);
        else if(res.vc_group[g] >= 0) res.share[g] = 1. / static_cast<real_type>(count_grp[res.vc_group[g]]);
    }
    return res;
}

void identity_rows(const AdjointSolveInfo & info, const std::vector<int> * masked,
                   const IntVect & p_row, const IntVect & q_row, std::vector<int> & out)
{
    const auto push = [&out](const IntVect & map, int bus){
        if(bus >= 0 && bus < map.size() && map[bus] >= 0) out.push_back(map[bus]);
    };
    if(masked != nullptr){
        for(int bus : *masked){
            push(p_row, bus);
            push(q_row, bus);
        }
    }
    for(int bus : info.pinned_buses) push(q_row, bus);
}

real_type share_in_solve(const GenVLayout & layout, const AdjointSolveInfo & info,
                         const IntVect & vm_col, Eigen::Index g)
{
    if(g < 0 || g >= layout.share.size()) return 0.;
    const int b = layout.target_bus[g];
    if(b >= 0){
        // a magnitude unknown: fixed by the set-point only while the bus is pinned PV
        const bool unknown = b < vm_col.size() && vm_col[b] >= 0;
        return unknown && !contains(info.pinned_buses, b) ? 0. : layout.share[g];
    }
    const int grp = layout.vc_group[g];
    if(grp >= 0) return contains(info.vset_lost_groups, grp) ? 0. : layout.share[g];
    return 0.;
}

void gen_v_indirect(const GenVLayout & layout, const AdjointSolveInfo & info,
                    const IntVect & vm_col, const IntVect & group_v_row,
                    const real_type * lambda, real_type * out)
{
    const Eigen::Index nb_gen = layout.share.size();
    // -lambda^T dF/d|V_b| for each bus whose column the solve gave, once per bus
    std::vector<real_type> acc(info.vm_col_bus.size(), 0.);
    for(std::size_t k = 0; k < info.vm_col_bus.size(); ++k){
        real_type a = 0.;
        for(int p = info.vm_col_start[k]; p < info.vm_col_start[k + 1]; ++p){
            a -= lambda[info.vm_col_row[static_cast<std::size_t>(p)]] * info.vm_col_val[static_cast<std::size_t>(p)];
        }
        acc[k] = a;
    }
    for(Eigen::Index g = 0; g < nb_gen; ++g){
        out[g] = 0.;
        const real_type s = share_in_solve(layout, info, vm_col, g);
        if(s == 0.) continue;
        const int b = layout.target_bus[g];
        if(b >= 0){
            // the set-point fixes |V_b|: through the starting voltage (no unknown), or through
            // the pinned row |V_b| = v (the unknown's own column, whose identity entry meets a
            // dropped multiplier)
            const auto it = std::lower_bound(info.vm_col_bus.begin(), info.vm_col_bus.end(), b);
            if(it == info.vm_col_bus.end() || *it != b) continue;
            out[g] = acc[static_cast<std::size_t>(it - info.vm_col_bus.begin())] * s;
            continue;
        }
        // a group's v_set: dF_v/dv_set = -1, so lambda at the group's voltage row
        const int grp = layout.vc_group[g];
        if(grp < group_v_row.size() && group_v_row[grp] >= 0) out[g] = lambda[group_v_row[grp]] * s;
    }
}

void gen_p_indirect(const AdjointSolveInfo & info, const real_type * lambda, real_type * out, Eigen::Index nb_gen)
{
    for(const AdjointSolveInfo::GenPTerm & term : info.gen_p_terms){
        if(term.gen_id < 0 || term.gen_id >= nb_gen) continue;
        out[term.gen_id] += lambda[term.row] * term.coeff;
    }
}

IntVect bus_map_to_grid(const SolverBusLayout & layout, const IntVect & solver_keyed)
{
    const auto me_to_solver = layout.id_me_to_solver.as_eigen();
    IntVect res = IntVect::Constant(me_to_solver.size(), -1);
    for(Eigen::Index bus_me = 0; bus_me < me_to_solver.size(); ++bus_me){
        const int bus_solver = me_to_solver[bus_me];
        if(bus_solver < 0) continue;   // bus not in the solver (deactivated, isolated)
        if(bus_solver >= solver_keyed.size()) continue;
        res[bus_me] = solver_keyed[bus_solver];
    }
    return res;
}

}  // namespace power_flow_adjoint
}  // namespace ls2g
