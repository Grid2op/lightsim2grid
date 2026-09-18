// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SUBSTATION_TOPOLOGY_H
#define SUBSTATION_TOPOLOGY_H

#include <string>
#include <tuple>
#include <vector>

#include "Eigen/Core"

#include "Utils.hpp"
#include "ls2g_api.hpp"

namespace ls2g {

/**
 * One terminal of an element inside a substation: which element, which end of
 * it, and the (substation-local) node it stands on.
 *
 * Derived from the element containers (each terminal carries its own node id,
 * see OneSideContainer::set_node_id) and rebuilt from them, never serialized:
 * it is the per-substation index the projection of the switches onto the
 * elements walks, nothing more.
 */
struct Terminal
{
    TerminalKind kind;
    int el_id;
    int node;
};

/**
 * The detailed topology of ONE substation (a pypowsybl voltage level), in the
 * node-breaker view: `nb_nodes_` connectivity nodes numbered 0..nb_nodes_-1 like
 * pypowsybl's, busbar sections standing on some of them, and switches joining
 * pairs of them. Everything here is local to the substation -- a node id, a
 * switch id or a busbar-section id only means something together with the
 * substation it belongs to. SubstationContainer owns one of these per
 * substation and is what turns local ids into the grid-wide ones the python
 * surface speaks.
 *
 * This is NOT an element container: it stamps nothing into Ybus or Sbus and
 * takes part in no per-element loop. Its only job is to hold the switch graph
 * and, later, to label its connected components as electrical buses.
 *
 * What is held here is the DECLARATION of the topology. Whether the elements'
 * buses currently agree with the switch positions is a separate question
 * (the projection), answered elsewhere.
 */
class LS2G_API SubstationTopology final
{
    public:
        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)
        using StateRes = std::tuple<
            int,                       // nb_nodes_
            std::vector<int>,          // bbs_node_
            std::vector<std::string>,  // bbs_names_ (optional: empty, or one per busbar section)
            std::vector<int>,          // sw_node1_
            std::vector<int>,          // sw_node2_
            std::vector<int>,          // sw_kind_ (SwitchKind, as its integer value)
            std::vector<bool>,         // sw_open_
            std::vector<bool>,         // sw_retained_
            std::vector<std::string>   // sw_names_ (optional: empty, or one per switch)
            >;
        enum StateResIdx {
            NB_NODES = 0,
            BBS_NODE,
            BBS_NAMES,
            SW_NODE1,
            SW_NODE2,
            SW_KIND,
            SW_OPEN,
            SW_RETAINED,
            SW_NAMES,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "SubstationTopology::StateRes and StateResIdx do not match");

        SubstationTopology() noexcept = default;

        /**
         * Declare the topology: `nb_nodes` nodes, one busbar section per entry of
         * `bbs_node`, one switch per entry of the `sw_*` vectors. Everything is
         * validated first (see check_valid), so a rejected call leaves the object
         * untouched.
         */
        void init(int nb_nodes,
                  const IntVect & bbs_node,
                  const IntVect & sw_node1,
                  const IntVect & sw_node2,
                  const std::vector<SwitchKind> & sw_kind,
                  const std::vector<bool> & sw_open,
                  const std::vector<bool> & sw_retained);

        int nb_nodes() const { return nb_nodes_; }
        int nb_switches() const { return static_cast<int>(sw_node1_.size()); }
        int nb_busbar_sections() const { return static_cast<int>(bbs_node_.size()); }

        // busbar sections (ids local to this substation, range-checked)
        int bbs_node(int bbs_id) const { return bbs_node_(_checked_bbs_id(bbs_id, "bbs_node")); }
        std::string bbs_name(int bbs_id) const;  // "" when no names were set
        const std::vector<std::string> & bbs_names() const { return bbs_names_; }
        void set_bbs_names(const std::vector<std::string> & names);  // one per busbar section

        // switches (ids local to this substation, range-checked)
        int sw_node1(int sw_id) const { return sw_node1_(_checked_sw_id(sw_id, "sw_node1")); }
        int sw_node2(int sw_id) const { return sw_node2_(_checked_sw_id(sw_id, "sw_node2")); }
        SwitchKind sw_kind(int sw_id) const { return sw_kind_[_checked_sw_id(sw_id, "sw_kind")]; }
        bool is_open(int sw_id) const { return sw_open_[_checked_sw_id(sw_id, "is_open")]; }
        bool is_retained(int sw_id) const { return sw_retained_[_checked_sw_id(sw_id, "is_retained")]; }
        const std::vector<bool> & sw_open() const { return sw_open_; }
        std::string sw_name(int sw_id) const;  // "" when no names were set
        const std::vector<std::string> & sw_names() const { return sw_names_; }
        void set_sw_names(const std::vector<std::string> & names);  // one per switch
        /**
         * Open or close a switch. Returns whether anything changed. An internal
         * connection cannot be opened (it is not a switch anybody operates).
         *
         * Changes the DECLARED position only: the labels below are stale until
         * label() runs again, and whoever owns the elements is the one to push the
         * new labels onto them (see LSGrid).
         */
        bool set_open(int sw_id, bool open);

        // ---- labels: the electrical buses the closed switches make ---------------
        /**
         * Label the connected components of the closed-switch graph as electrical
         * buses, pypowsybl's way (powsybl-core `Networks.isBusValid`): a component
         * is a bus iff
         *
         *     (n_busbar_sections >= 1 && n_feeders >= 1) || (n_branches >= 1 && n_feeders >= 2)
         *
         * where every terminal is a feeder, and a branch is a line end, a
         * transformer end or an HVDC converter station. An isolated busbar section,
         * a lone line end behind an open breaker, a load with a generator and
         * nothing else: none of these is a bus, and their elements are disconnected.
         *
         * Numbering, 1-based like LocalBusId: the components holding a busbar
         * section first, in busbar-section order (what grid2op's
         * `from_switches_position` does too), then the remaining valid ones by
         * lowest node. Deterministic: a function of the switch positions only.
         *
         * `nmax_busbar_per_sub` is the substation's capacity in the grid's bus
         * layout: more buses than that cannot be numbered, and this throws BEFORE
         * touching the current labels (the loader sizes every substation for its
         * true maximum, so that is a hand-built grid's error). `sub_id` only names
         * the substation in that message.
         */
        void label(int nmax_busbar_per_sub, int sub_id = -1);
        /// have the labels been computed for the current switch positions?
        bool labels_ready() const { return labels_ready_; }
        /// how many electrical buses the current labels count
        int nb_buses() const { return nb_buses_; }
        /// the local bus (1-based) of every node, -1 for a node in no valid component
        const IntVect & node_bus() const { return node_bus_; }
        int node_bus(int node) const;  // range-checked
        /// the local bus of a busbar section (-1: isolated, or no feeder reaches it)
        int bbs_bus(int bbs_id) const { return node_bus(bbs_node(bbs_id)); }

        // ---- results: the voltage of each busbar section after a powerflow --------
        // Filled by LSGrid::compute_results (a section reads its bus' voltage),
        // cleared by reset_results; derived, never serialized.
        bool has_bbs_results() const { return has_bbs_res_; }
        real_type bbs_res_v_kv(int bbs_id) const { return bbs_res_v_kv_(_checked_bbs_id(bbs_id, "bbs_res_v_kv")); }
        real_type bbs_res_theta_deg(int bbs_id) const { return bbs_res_theta_deg_(_checked_bbs_id(bbs_id, "bbs_res_theta_deg")); }
        void set_bbs_results(const RealVect & v_kv, const RealVect & theta_deg);  // one per section
        void reset_bbs_results();

        // ---- terminals: which element ends stand on which node -----------------
        // Derived, rebuilt by LSGrid from the containers (see
        // LSGrid::_rebuild_terminal_lists); not part of StateRes.
        void clear_terminals() { terminals_.clear(); _invalidate_labels(); }
        void add_terminal(TerminalKind kind, int el_id, int node);  // node is range-checked
        const std::vector<Terminal> & terminals() const { return terminals_; }

        // ---- state -------------------------------------------------------------
        StateRes get_state() const;
        // validates the whole state before assigning anything (a rejected state
        // leaves the object untouched)
        void set_state(StateRes & my_state);
        /**
         * Consistency of the declaration: every node referenced exists, no switch
         * joins a node to itself, no two busbar sections share a node, an internal
         * connection is closed, optional names are absent or complete. `sub_id`
         * only names the substation in error messages. Throws std::runtime_error /
         * std::out_of_range.
         */
        void check_valid(int sub_id = -1) const;

    private:
        static void _check(int nb_nodes,
                           const IntVect & bbs_node,
                           const std::vector<std::string> & bbs_names,
                           const IntVect & sw_node1,
                           const IntVect & sw_node2,
                           const std::vector<SwitchKind> & sw_kind,
                           const std::vector<bool> & sw_open,
                           const std::vector<bool> & sw_retained,
                           const std::vector<std::string> & sw_names,
                           int sub_id);
        int _checked_bbs_id(int bbs_id, const char * fun_name) const;
        int _checked_sw_id(int sw_id, const char * fun_name) const;
        /// the labels no longer describe the switch positions (a switch moved, a
        /// terminal was added, the declaration changed)
        void _invalidate_labels() { labels_ready_ = false; }
        // union-find over the nodes (scratch for label())
        int _find(int node);

    private:
        int nb_nodes_ = 0;

        // busbar sections
        IntVect bbs_node_;
        std::vector<std::string> bbs_names_;  // optional

        // switches
        IntVect sw_node1_;
        IntVect sw_node2_;
        std::vector<SwitchKind> sw_kind_;
        std::vector<bool> sw_open_;
        std::vector<bool> sw_retained_;
        std::vector<std::string> sw_names_;  // optional

        // derived, never serialized
        std::vector<Terminal> terminals_;
        // labels (derived, never serialized): see label()
        bool labels_ready_ = false;
        int nb_buses_ = 0;
        IntVect node_bus_;
        // union-find scratch, sized nb_nodes_ by label()
        std::vector<int> uf_parent_;
        // busbar-section results (derived, never serialized): see set_bbs_results
        bool has_bbs_res_ = false;
        RealVect bbs_res_v_kv_;
        RealVect bbs_res_theta_deg_;
};

}  // namespace ls2g

#endif  // SUBSTATION_TOPOLOGY_H
