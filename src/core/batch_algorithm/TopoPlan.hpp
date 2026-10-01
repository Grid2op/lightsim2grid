// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TOPOPLAN_H
#define TOPOPLAN_H

#include "LSGrid.hpp"
#include "light_env/topo_action.hpp"

#include <map>
#include <vector>

namespace ls2g {

/**
 * What a `TopoAction` changes on a grid, as a diff against the grid's base state and
 * WITHOUT mutating it: the batch sweeps (ScenarioSweep::set_topo_actions) replay a row's
 * action as value edits on a fixed solver layout, so they need to know, element by
 * element, where the row puts it -- not a mutated grid.
 *
 * Buses are gridmodel bus ids (`sub_id + (busbar - 1) * n_sub`); -1 means "off". Only
 * real changes are listed: an element the action leaves where it is, or disconnects
 * while it already is, does not appear.
 */
struct TopoElPlacement {
    int el_id;         // id inside its container
    int base_bus_me;   // the grid's own placement, -1 if the element is off in the base grid
    int row_bus_me;    // the row's, -1 if the row disconnects it
};
struct TopoBranchPlacement {
    int branch_id;     // gridmodel numbering: powerlines then transformers
    int base_bus1_me;  // the buses the branch stands on (kept while off: a reconnection
    int base_bus2_me;  // without set_bus goes back there)
    bool base_on;
    int row_bus1_me;
    int row_bus2_me;
    bool row_on;
};
struct RowTopoPlan {
    std::vector<TopoElPlacement> loads;
    std::vector<TopoElPlacement> gens;
    std::vector<TopoElPlacement> storages;
    std::vector<TopoBranchPlacement> branches;
    bool empty() const {
        return loads.empty() && gens.empty() && storages.empty() && branches.empty();
    }
};

namespace topo_plan_detail {

template<class ContainerType>
inline void resolve_one_side(const std::vector<TopoAction::SetBusEntry> & entries,
                             TopoAction::Kind kind,
                             const ContainerType & container,
                             std::vector<TopoElPlacement> & out)
{
    // the action may name an element twice (last wins in TopoAction::add_element,
    // but the resolved list is one entry per element anyway); a map keeps one
    std::map<int, TopoElPlacement> placements;
    const std::vector<bool> & status = container.get_status();
    const GlobalBusIdVect & bus_id = container.get_bus_id();
    for(const auto & sb : entries){
        if(sb.kind != kind) continue;
        TopoElPlacement & p = placements[sb.el_id];
        p.el_id = sb.el_id;
        p.base_bus_me = status[sb.el_id] ? bus_id(sb.el_id).cast_int() : BaseConstants::_deactivated_bus_id;
        p.row_bus_me = sb.local_bus > 0 ? sb.global_bus.cast_int() : BaseConstants::_deactivated_bus_id;
    }
    for(const auto & el : placements){
        if(el.second.row_bus_me != el.second.base_bus_me) out.push_back(el.second);
    }
}

}  // namespace topo_plan_detail

/**
 * Replays `action` (checked against `grid`, see TopoAction::check_validity) on the base
 * state of `grid` and lists what it changes. Same semantics as
 * TopoAction::apply_to_gridmodel: set_line_status first, then set_bus; -1 on either end
 * of a branch disconnects the whole branch; a busbar on one end of a disconnected branch
 * reconnects it, the other end going back to the bus it was last on.
 */
inline void resolve_row_topo(const TopoAction & action, const LSGrid & grid, RowTopoPlan & out)
{
    out = RowTopoPlan();
    const int nb_line = static_cast<int>(grid.nb_powerline());
    const LineContainer & lines = grid.get_powerlines_as_data();
    const TrafoContainer & trafos = grid.get_trafos_as_data();

    std::map<int, TopoBranchPlacement> branches;
    auto touch_branch = [&](int branch_id) -> TopoBranchPlacement & {
        auto it = branches.find(branch_id);
        if(it != branches.end()) return it->second;
        const bool is_trafo = branch_id >= nb_line;
        const int internal_id = is_trafo ? branch_id - nb_line : branch_id;
        const BranchContainer & container = is_trafo ? static_cast<const BranchContainer &>(trafos)
                                                     : static_cast<const BranchContainer &>(lines);
        TopoBranchPlacement p;
        p.branch_id = branch_id;
        p.base_bus1_me = container.get_bus_id_side_1()(internal_id).cast_int();
        p.base_bus2_me = container.get_bus_id_side_2()(internal_id).cast_int();
        p.base_on = container.get_status_global()[internal_id];
        p.row_bus1_me = p.base_bus1_me;
        p.row_bus2_me = p.base_bus2_me;
        p.row_on = p.base_on;
        return branches.emplace(branch_id, p).first->second;
    };

    for(const auto & ls : action.line_status_entries()){
        TopoBranchPlacement & p = touch_branch(ls.g2op_line_id);
        p.row_on = ls.status == 1;
    }
    for(const auto & sb : action.set_bus_entries()){
        const bool side1 = sb.kind == TopoAction::Kind::line_side1 || sb.kind == TopoAction::Kind::trafo_side1;
        const bool side2 = sb.kind == TopoAction::Kind::line_side2 || sb.kind == TopoAction::Kind::trafo_side2;
        if(!side1 && !side2) continue;
        TopoBranchPlacement & p = touch_branch(sb.g2op_line_id);
        if(sb.local_bus == -1){
            p.row_on = false;
        }else{
            p.row_on = true;
            (side1 ? p.row_bus1_me : p.row_bus2_me) = sb.global_bus.cast_int();
        }
    }
    for(const auto & br : branches){
        const TopoBranchPlacement & p = br.second;
        const bool same = (p.row_on == p.base_on) &&
                          (!p.row_on || (p.row_bus1_me == p.base_bus1_me && p.row_bus2_me == p.base_bus2_me));
        if(!same) out.branches.push_back(p);
    }

    topo_plan_detail::resolve_one_side(action.set_bus_entries(), TopoAction::Kind::load, grid.get_loads_as_data(), out.loads);
    topo_plan_detail::resolve_one_side(action.set_bus_entries(), TopoAction::Kind::gen, grid.get_generators_as_data(), out.gens);
    topo_plan_detail::resolve_one_side(action.set_bus_entries(), TopoAction::Kind::storage, grid.get_storages(), out.storages);
}

} // namespace ls2g

#endif  // TOPOPLAN_H
