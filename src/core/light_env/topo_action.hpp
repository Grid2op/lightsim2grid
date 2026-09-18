// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TOPO_ACTION_H
#define TOPO_ACTION_H

#include "LSGrid.hpp"

#include "light_env_utils.hpp"
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ls2g {

/**
 * A grid2op-like topological action for the light environment.
 *
 * It is filled with grid2op semantics:
 *
 * - `add_element(el_type, el_id, local_bus)` is grid2op's `set_bus`: `local_bus` is -1
 *   (disconnect the element), 0 (leave it unchanged) or 1..n_busbar_per_sub (connect it to
 *   that busbar of the substation the element belongs to). For `line_or` / `line_ex`,
 *   `el_id` is the grid2op line id: powerlines first, then transformers. `trafo_hv` /
 *   `trafo_lv` take the transformer id directly.
 * - `set_line_status(line_id, status)` is grid2op's `set_line_status`: -1 disconnects the
 *   line, +1 reconnects it on the buses it was last connected to, 0 leaves it unchanged.
 *   `line_id` is the grid2op line id (powerlines then transformers).
 *
 * Nothing is checked at fill time. `check_validity(grid)` validates every entry against the
 * grid (the element exists, the busbar exists, no contradictory entries, supported element
 * type) and resolves the local busbars into gridmodel bus ids. It throws
 * `std::invalid_argument` (bad busbar, contradiction, unsupported type) or
 * `std::out_of_range` (bad element id) on an invalid action. `compute_impact` and
 * `apply_to_gridmodel` require a checked action.
 */
class TopoAction
{
    protected:
        typedef std::map<int, int> ElBusMapping;

        enum class Kind {load, gen, storage, line_side1, line_side2, trafo_side1, trafo_side2};

        struct SetBusEntry {
            Kind kind;
            int el_id;         // id inside its LSGrid container
            int g2op_line_id;  // for line / trafo sides: the grid2op line id, else -1
            int sub_id;
            int local_bus;     // -1 or 1..n_busbar_per_sub (0 entries are dropped)
            GridModelBusId global_bus;  // meaningful only when local_bus > 0
        };
        struct LineStatusEntry {
            int g2op_line_id;
            bool is_trafo;
            int internal_id;
            int status;  // -1 or +1 (0 entries are dropped)
        };

    public:
        TopoAction(): checked_(false) {}

        /**
         * grid2op `set_bus` on one element. See the class documentation for the meaning of
         * `local_bus_id`. Calling it twice on the same element keeps the last value.
         */
        void add_element(ElementType el_type, int el_id, int local_bus_id){
            ElBusMapping & buses_tab = choose_bus_table(el_type);
            buses_tab[el_id] = local_bus_id;
            checked_ = false;
        }

        /**
         * grid2op `set_line_status` on one line (grid2op numbering: powerlines then
         * transformers): -1 disconnect, +1 reconnect, 0 nothing.
         */
        void set_line_status(int line_id, int status){
            lines_status_[line_id] = status;
            checked_ = false;
        }

        // true if applying this action would not touch anything (before or after the check)
        bool is_do_nothing() const {
            for(const ElBusMapping * tab : {&loads_bus_, &gens_bus_, &storages_bus_, &lines_or_bus_,
                                            &lines_ex_bus_, &trafo_hv_bus_, &trafo_lv_bus_, &lines_status_}){
                for(const auto & el : *tab) if(el.second != 0) return false;
            }
            return true;
        }

        bool has_been_checked() const {return checked_;}
        int nb_set_bus() const {return static_cast<int>(set_bus_.size());}
        int nb_set_line_status() const {return static_cast<int>(line_status_.size());}

        /**
         * Check every entry against `grid` and resolve the busbars. Throws on an invalid
         * action (see the class documentation). Must be called again if the action is
         * modified afterwards.
         */
        void check_validity(const LSGrid & grid) {
            checked_ = false;
            set_bus_.clear();
            line_status_.clear();

            const int n_busbar = grid.get_max_nb_bus_per_sub();
            if(n_busbar <= 0){
                throw std::runtime_error("TopoAction::check_validity: the grid does not know its number "
                                         "of busbars per substation (LSGrid.set_max_nb_bus_per_sub was never called).");
            }
            const int nb_line = static_cast<int>(grid.nb_powerline());
            const int nb_trafo = static_cast<int>(grid.nb_trafo());
            const int n_line_g2op = nb_line + nb_trafo;

            aux_resolve_one_side(loads_bus_, Kind::load, "load", grid.get_loads().nb(),
                                 grid.get_loads().get_subid(), n_busbar, grid);
            aux_resolve_one_side(gens_bus_, Kind::gen, "generator", grid.get_generators().nb(),
                                 grid.get_generators().get_subid(), n_busbar, grid);
            aux_resolve_one_side(storages_bus_, Kind::storage, "storage", grid.get_storages().nb(),
                                 grid.get_storages().get_subid(), n_busbar, grid);

            // line ends, grid2op numbering (powerlines then transformers)
            for(const auto & el : lines_or_bus_) aux_resolve_line_side(el.first, el.second, true, "line_or", nb_line, n_line_g2op, n_busbar, grid);
            for(const auto & el : lines_ex_bus_) aux_resolve_line_side(el.first, el.second, false, "line_ex", nb_line, n_line_g2op, n_busbar, grid);
            // transformer ends, transformer numbering
            for(const auto & el : trafo_hv_bus_){
                aux_check_el_id(el.first, nb_trafo, "trafo_hv");
                aux_resolve_line_side(el.first + nb_line, el.second, true, "trafo_hv", nb_line, n_line_g2op, n_busbar, grid);
            }
            for(const auto & el : trafo_lv_bus_){
                aux_check_el_id(el.first, nb_trafo, "trafo_lv");
                aux_resolve_line_side(el.first + nb_line, el.second, false, "trafo_lv", nb_line, n_line_g2op, n_busbar, grid);
            }

            // line status
            for(const auto & el : lines_status_){
                const int line_id = el.first;
                const int status = el.second;
                if(status == 0) continue;
                aux_check_el_id(line_id, n_line_g2op, "line (set_line_status)");
                if(status != -1 && status != 1){
                    std::ostringstream exc_;
                    exc_ << "TopoAction::check_validity: set_line_status for line " << line_id
                         << " should be -1 (disconnect), 0 (nothing) or +1 (reconnect), you provided " << status << ".";
                    throw std::invalid_argument(exc_.str());
                }
                LineStatusEntry e;
                e.g2op_line_id = line_id;
                e.is_trafo = line_id >= nb_line;
                e.internal_id = e.is_trafo ? line_id - nb_line : line_id;
                e.status = status;
                line_status_.push_back(e);
            }

            // contradictions between set_line_status and set_bus on the line ends
            for(const auto & ls : line_status_){
                for(const auto & sb : set_bus_){
                    if(sb.g2op_line_id != ls.g2op_line_id) continue;
                    if((ls.status == -1 && sb.local_bus > 0) || (ls.status == 1 && sb.local_bus == -1)){
                        std::ostringstream exc_;
                        exc_ << "TopoAction::check_validity: line " << ls.g2op_line_id << " is "
                             << (ls.status == -1 ? "disconnected" : "reconnected")
                             << " with set_line_status and one of its ends is at the same time "
                             << (sb.local_bus > 0 ? "connected to a busbar" : "disconnected")
                             << " with set_bus. This is ambiguous.";
                        throw std::invalid_argument(exc_.str());
                    }
                }
            }
            checked_ = true;
        }

        /**
         * Which substations and which lines (grid2op numbering) this action impacts, given the
         * current state of the grid. Follows grid2op's `get_topological_impact`: a line
         * whose status is changed (by set_line_status, or by set_bus on one of its ends that
         * disconnects a connected line or reconnects a disconnected one) impacts the line and
         * not its substations; every other set_bus impacts the substation of the element.
         */
        void compute_impact(const LSGrid & grid,
                            std::vector<bool> & subs_impacted,
                            std::vector<bool> & lines_impacted) const {
            aux_require_checked("compute_impact");
            const int nb_line = static_cast<int>(grid.nb_powerline());
            subs_impacted.assign(grid.get_n_sub(), false);
            lines_impacted.assign(nb_line + static_cast<int>(grid.nb_trafo()), false);
            for(const auto & ls : line_status_) lines_impacted[ls.g2op_line_id] = true;
            const std::vector<bool> & line_status = grid.get_lines_status();
            const std::vector<bool> & trafo_status = grid.get_trafo_status();
            for(const auto & sb : set_bus_){
                switch (sb.kind)
                {
                case Kind::load:
                case Kind::gen:
                case Kind::storage:
                    subs_impacted[sb.sub_id] = true;
                    break;
                default:
                    {
                        const bool is_trafo = (sb.kind == Kind::trafo_side1) || (sb.kind == Kind::trafo_side2);
                        const bool connected = is_trafo ? trafo_status[sb.el_id] : line_status[sb.el_id];
                        if(sb.local_bus == -1){
                            // disconnecting a line: a line status change (nothing if already disconnected)
                            if(connected) lines_impacted[sb.g2op_line_id] = true;
                        }else if(connected){
                            // moving the end of a connected line: a topological change of its substation
                            subs_impacted[sb.sub_id] = true;
                        }else{
                            // reconnecting a line: a line status change
                            lines_impacted[sb.g2op_line_id] = true;
                        }
                    }
                }
            }
        }

        void apply_to_gridmodel(LSGrid & grid) const {
            aux_require_checked("apply_to_gridmodel");
            for(const auto & ls : line_status_){
                if(ls.status == -1){
                    if(ls.is_trafo) grid.deactivate_trafo(ls.internal_id);
                    else grid.deactivate_powerline(ls.internal_id);
                }else{
                    if(ls.is_trafo) grid.reactivate_trafo(ls.internal_id);
                    else grid.reactivate_powerline(ls.internal_id);
                }
            }
            for(const auto & sb : set_bus_){
                switch (sb.kind)
                {
                case Kind::load:
                    aux_apply_one_side(grid, sb, grid.get_loads(), &LSGrid::deactivate_load,
                                       &LSGrid::reactivate_load, &LSGrid::change_bus_load);
                    break;
                case Kind::gen:
                    aux_apply_one_side(grid, sb, grid.get_generators(), &LSGrid::deactivate_gen,
                                       &LSGrid::reactivate_gen, &LSGrid::change_bus_gen);
                    break;
                case Kind::storage:
                    aux_apply_one_side(grid, sb, grid.get_storages(), &LSGrid::deactivate_storage,
                                       &LSGrid::reactivate_storage, &LSGrid::change_bus_storage);
                    break;
                case Kind::line_side1:
                    aux_apply_branch(grid, sb, grid.get_lines_status(), &LSGrid::deactivate_powerline,
                                     &LSGrid::reactivate_powerline, &LSGrid::change_bus1_powerline);
                    break;
                case Kind::line_side2:
                    aux_apply_branch(grid, sb, grid.get_lines_status(), &LSGrid::deactivate_powerline,
                                     &LSGrid::reactivate_powerline, &LSGrid::change_bus2_powerline);
                    break;
                case Kind::trafo_side1:
                    aux_apply_branch(grid, sb, grid.get_trafo_status(), &LSGrid::deactivate_trafo,
                                     &LSGrid::reactivate_trafo, &LSGrid::change_bus1_trafo);
                    break;
                case Kind::trafo_side2:
                    aux_apply_branch(grid, sb, grid.get_trafo_status(), &LSGrid::deactivate_trafo,
                                     &LSGrid::reactivate_trafo, &LSGrid::change_bus2_trafo);
                    break;
                }
            }
        }

    protected:
        void aux_require_checked(const std::string & fun_name) const {
            if(!checked_){
                throw std::runtime_error("TopoAction::" + fun_name + ": the action has not been checked "
                                         "against a grid (call check_validity first).");
            }
        }

        static void aux_check_el_id(int el_id, int nb_el, const std::string & name){
            if(el_id < 0 || el_id >= nb_el){
                std::ostringstream exc_;
                exc_ << "TopoAction::check_validity: " << name << " with id " << el_id
                     << " does not exist, the grid has " << nb_el << " of them (valid ids: 0 to " << nb_el - 1 << ").";
                throw std::out_of_range(exc_.str());
            }
        }

        static void aux_check_local_bus(int local_bus, int n_busbar, const std::string & name, int el_id){
            if(local_bus < -1 || local_bus > n_busbar){
                std::ostringstream exc_;
                exc_ << "TopoAction::check_validity: " << name << " " << el_id << " is assigned to busbar "
                     << local_bus << ". A busbar should be -1 (disconnected), 0 (unchanged) or between 1 and "
                     << n_busbar << " (number of busbars per substation on this grid).";
                throw std::invalid_argument(exc_.str());
            }
        }

        static void aux_check_subid(const IntVect & subids, int nb_el, const std::string & name){
            if(subids.size() != nb_el){
                std::ostringstream exc_;
                exc_ << "TopoAction::check_validity: the grid does not know the substation of its " << name
                     << "s (set_..._to_subid was never called), a set_bus action cannot be resolved.";
                throw std::runtime_error(exc_.str());
            }
        }

        void aux_push_entry(Kind kind, int el_id, int g2op_line_id, int sub_id, int local_bus, const LSGrid & grid){
            SetBusEntry e;
            e.kind = kind;
            e.el_id = el_id;
            e.g2op_line_id = g2op_line_id;
            e.sub_id = sub_id;
            e.local_bus = local_bus;
            e.global_bus = (local_bus > 0) ?
                           grid.get_substations().local_to_gridmodel(sub_id, LocalBusId(local_bus)) :
                           GridModelBusId(BaseConstants::_deactivated_bus_id);
            set_bus_.push_back(e);
        }

        void aux_resolve_one_side(const ElBusMapping & tab, Kind kind, const std::string & name,
                                  int nb_el, const IntVect & subids, int n_busbar, const LSGrid & grid){
            for(const auto & el : tab){
                const int el_id = el.first;
                const int local_bus = el.second;
                if(local_bus == 0) continue;
                aux_check_el_id(el_id, nb_el, name);
                aux_check_local_bus(local_bus, n_busbar, name, el_id);
                aux_check_subid(subids, nb_el, name);
                aux_push_entry(kind, el_id, -1, subids(el_id), local_bus, grid);
            }
        }

        void aux_resolve_line_side(int g2op_line_id, int local_bus, bool side1, const std::string & name,
                                   int nb_line, int n_line_g2op, int n_busbar, const LSGrid & grid){
            if(local_bus == 0) return;
            aux_check_el_id(g2op_line_id, n_line_g2op, name);
            aux_check_local_bus(local_bus, n_busbar, name, g2op_line_id);
            const bool is_trafo = g2op_line_id >= nb_line;
            const int internal_id = is_trafo ? g2op_line_id - nb_line : g2op_line_id;
            const IntVect & subids = is_trafo ?
                (side1 ? grid.get_trafos().get_subid_side_1() : grid.get_trafos().get_subid_side_2()) :
                (side1 ? grid.get_lines().get_subid_side_1() : grid.get_lines().get_subid_side_2());
            aux_check_subid(subids, is_trafo ? static_cast<int>(grid.nb_trafo()) : nb_line, name);
            Kind kind = is_trafo ? (side1 ? Kind::trafo_side1 : Kind::trafo_side2)
                                 : (side1 ? Kind::line_side1 : Kind::line_side2);
            aux_push_entry(kind, internal_id, g2op_line_id, subids(internal_id), local_bus, grid);
        }

        template<class ContainerType, class FunDeact, class FunReact, class FunChange>
        static void aux_apply_one_side(LSGrid & grid, const SetBusEntry & sb, const ContainerType & container,
                                       FunDeact fun_deact, FunReact fun_react, FunChange fun_change){
            if(sb.local_bus == -1){
                (grid.*fun_deact)(sb.el_id);
            }else{
                if(!container.get_status(sb.el_id)) (grid.*fun_react)(sb.el_id);
                (grid.*fun_change)(sb.el_id, sb.global_bus);
            }
        }

        template<class FunDeact, class FunReact, class FunChange>
        static void aux_apply_branch(LSGrid & grid, const SetBusEntry & sb, const std::vector<bool> & status,
                                     FunDeact fun_deact, FunReact fun_react, FunChange fun_change){
            if(sb.local_bus == -1){
                // grid2op: disconnecting one end disconnects the whole line
                (grid.*fun_deact)(sb.el_id);
            }else{
                // grid2op: connecting one end of a disconnected line reconnects it
                // (the other end goes back to the bus it was last connected to)
                if(!status[sb.el_id]) (grid.*fun_react)(sb.el_id);
                (grid.*fun_change)(sb.el_id, sb.global_bus);
            }
        }

        ElBusMapping & choose_bus_table(ElementType el_type){
            switch (el_type)
            {
            case ElementType::load:
                return loads_bus_;
            case ElementType::gen:
                return gens_bus_;
            case ElementType::line_or:
                return lines_or_bus_;
            case ElementType::line_ex:
                return lines_ex_bus_;
            case ElementType::storage:
                return storages_bus_;
            case ElementType::trafo_hv:
                return trafo_hv_bus_;
            case ElementType::trafo_lv:
                return trafo_lv_bus_;
            default:
                std::ostringstream exc_;
                exc_ << "TopoAction: modification of element type ";
                exc_ << el_type << " is not supported at the moment.";
                throw std::invalid_argument(exc_.str());
            }
        }

    protected:
        // raw entries, grid2op semantics (local busbar ids)
        ElBusMapping loads_bus_;
        ElBusMapping gens_bus_;
        ElBusMapping storages_bus_;
        ElBusMapping lines_or_bus_;  // grid2op line numbering
        ElBusMapping lines_ex_bus_;  // grid2op line numbering
        ElBusMapping trafo_hv_bus_;  // transformer numbering
        ElBusMapping trafo_lv_bus_;  // transformer numbering
        ElBusMapping lines_status_;  // grid2op line numbering

        // resolved entries, filled by check_validity
        bool checked_;
        std::vector<SetBusEntry> set_bus_;
        std::vector<LineStatusEntry> line_status_;
};

} // namespace ls2g

#endif // TOPO_ACTION_H
