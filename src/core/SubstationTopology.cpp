// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SubstationTopology.hpp"

#include <sstream>
#include <stdexcept>

namespace ls2g {

namespace {

std::string sub_label(int sub_id)
{
    if(sub_id < 0) return "SubstationTopology";
    std::ostringstream res;
    res << "SubstationTopology (substation " << sub_id << ")";
    return res.str();
}

bool is_known_kind(SwitchKind kind)
{
    switch(kind){
        case SwitchKind::BREAKER:
        case SwitchKind::DISCONNECTOR:
        case SwitchKind::LOAD_BREAK_SWITCH:
        case SwitchKind::INTERNAL_CONNECTION:
            return true;
    }
    return false;
}

}  // anonymous namespace

void SubstationTopology::_check(int nb_nodes,
                                const IntVect & bbs_node,
                                const std::vector<std::string> & bbs_names,
                                const IntVect & sw_node1,
                                const IntVect & sw_node2,
                                const std::vector<SwitchKind> & sw_kind,
                                const std::vector<bool> & sw_open,
                                const std::vector<bool> & sw_retained,
                                const std::vector<std::string> & sw_names,
                                int sub_id)
{
    const std::string who = sub_label(sub_id);
    if(nb_nodes < 0){
        std::ostringstream exc_;
        exc_ << who << ": the number of nodes must be >= 0, got " << nb_nodes << ".";
        throw std::runtime_error(exc_.str());
    }
    const auto check_node = [&](int node, const char * what, Eigen::Index id){
        if((node < 0) || (node >= nb_nodes)){
            std::ostringstream exc_;
            exc_ << who << ": " << what << " " << id << " stands on node " << node
                 << ", out of range [0, " << nb_nodes << ").";
            throw std::out_of_range(exc_.str());
        }
    };

    // busbar sections: on an existing node, one per node at most (IIDM allows
    // one connectable per node, and two sections on one node would be the same
    // section twice for every purpose here)
    std::vector<char> node_has_bbs(static_cast<std::size_t>(nb_nodes), 0);
    for(Eigen::Index i = 0; i < bbs_node.size(); ++i){
        check_node(bbs_node(i), "busbar section", i);
        if(node_has_bbs[static_cast<std::size_t>(bbs_node(i))]){
            std::ostringstream exc_;
            exc_ << who << ": two busbar sections stand on node " << bbs_node(i)
                 << " (busbar section " << i << " and an earlier one); a node holds at most one.";
            throw std::runtime_error(exc_.str());
        }
        node_has_bbs[static_cast<std::size_t>(bbs_node(i))] = 1;
    }
    if(!bbs_names.empty() && (static_cast<Eigen::Index>(bbs_names.size()) != bbs_node.size())){
        std::ostringstream exc_;
        exc_ << who << ": " << bbs_names.size() << " busbar section names for "
             << bbs_node.size() << " busbar sections (names must be either absent or complete).";
        throw std::runtime_error(exc_.str());
    }

    // switches: the four parallel vectors have one entry per switch
    const Eigen::Index nb_sw = sw_node1.size();
    const auto check_sw_len = [&](std::size_t actual, const char * name){
        if(static_cast<Eigen::Index>(actual) != nb_sw){
            std::ostringstream exc_;
            exc_ << who << ": '" << name << "' has " << actual << " entries for " << nb_sw
                 << " switches.";
            throw std::runtime_error(exc_.str());
        }
    };
    check_sw_len(static_cast<std::size_t>(sw_node2.size()), "sw_node2");
    check_sw_len(sw_kind.size(), "sw_kind");
    check_sw_len(sw_open.size(), "sw_open");
    check_sw_len(sw_retained.size(), "sw_retained");
    if(!sw_names.empty()) check_sw_len(sw_names.size(), "sw_names");
    for(Eigen::Index i = 0; i < nb_sw; ++i){
        check_node(sw_node1(i), "switch (side 1)", i);
        check_node(sw_node2(i), "switch (side 2)", i);
        if(sw_node1(i) == sw_node2(i)){
            std::ostringstream exc_;
            exc_ << who << ": switch " << i << " joins node " << sw_node1(i)
                 << " to itself.";
            throw std::runtime_error(exc_.str());
        }
        if(!is_known_kind(sw_kind[static_cast<std::size_t>(i)])){
            std::ostringstream exc_;
            exc_ << who << ": switch " << i << " has an unknown kind ("
                 << static_cast<int>(sw_kind[static_cast<std::size_t>(i)]) << ").";
            throw std::runtime_error(exc_.str());
        }
        // an internal connection is a link nobody can open: "open" is not a
        // state it has, and a file claiming otherwise is inconsistent
        if((sw_kind[static_cast<std::size_t>(i)] == SwitchKind::INTERNAL_CONNECTION) &&
           sw_open[static_cast<std::size_t>(i)]){
            std::ostringstream exc_;
            exc_ << who << ": switch " << i << " is an internal connection and is open; "
                 << "an internal connection is always closed.";
            throw std::runtime_error(exc_.str());
        }
    }
}

void SubstationTopology::init(int nb_nodes,
                              const IntVect & bbs_node,
                              const IntVect & sw_node1,
                              const IntVect & sw_node2,
                              const std::vector<SwitchKind> & sw_kind,
                              const std::vector<bool> & sw_open,
                              const std::vector<bool> & sw_retained)
{
    const std::vector<std::string> no_names;
    _check(nb_nodes, bbs_node, no_names, sw_node1, sw_node2, sw_kind, sw_open, sw_retained, no_names, -1);
    nb_nodes_ = nb_nodes;
    bbs_node_ = bbs_node;
    bbs_names_.clear();
    sw_node1_ = sw_node1;
    sw_node2_ = sw_node2;
    sw_kind_ = sw_kind;
    sw_open_ = sw_open;
    sw_retained_ = sw_retained;
    sw_names_.clear();
    terminals_.clear();
    node_bus_ = IntVect::Constant(nb_nodes_, -1);
    nb_buses_ = 0;
    _invalidate_labels();
}

bool SubstationTopology::set_open(int sw_id, bool open)
{
    _checked_sw_id(sw_id, "set_open");
    const std::size_t idx = static_cast<std::size_t>(sw_id);
    if(sw_kind_[idx] == SwitchKind::INTERNAL_CONNECTION){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::set_open: switch " << sw_id;
        if(!sw_names_.empty()) exc_ << " ('" << sw_names_[idx] << "')";
        exc_ << " is an internal connection: it is always closed and cannot be operated.";
        throw std::runtime_error(exc_.str());
    }
    if(sw_open_[idx] == open) return false;
    sw_open_[idx] = open;
    _invalidate_labels();
    return true;
}

int SubstationTopology::_find(int node)
{
    // path halving: every other node on the way up is re-pointed at its grandparent
    while(uf_parent_[static_cast<std::size_t>(node)] != node){
        const int parent = uf_parent_[static_cast<std::size_t>(node)];
        const int grand_parent = uf_parent_[static_cast<std::size_t>(parent)];
        uf_parent_[static_cast<std::size_t>(node)] = grand_parent;
        node = grand_parent;
    }
    return node;
}

void SubstationTopology::label(int nmax_busbar_per_sub, int sub_id)
{
    const std::size_t n = static_cast<std::size_t>(nb_nodes_);

    // 1. the components of the closed-switch graph
    uf_parent_.resize(n);
    for(std::size_t i = 0; i < n; ++i) uf_parent_[i] = static_cast<int>(i);
    const int nb_sw = nb_switches();
    for(int sw = 0; sw < nb_sw; ++sw){
        if(sw_open_[static_cast<std::size_t>(sw)]) continue;  // an internal connection is never open
        const int root1 = _find(sw_node1_(sw));
        const int root2 = _find(sw_node2_(sw));
        if(root1 != root2) uf_parent_[static_cast<std::size_t>(root1)] = root2;
    }

    // 2. what each component holds, indexed by root
    std::vector<int> n_bbs(n, 0), n_branch(n, 0), n_feeder(n, 0);
    for(Eigen::Index i = 0; i < bbs_node_.size(); ++i) ++n_bbs[static_cast<std::size_t>(_find(bbs_node_(i)))];
    for(const Terminal & term : terminals_){
        const std::size_t root = static_cast<std::size_t>(_find(term.node));
        ++n_feeder[root];
        switch(term.kind){
            case TerminalKind::LINE_1:
            case TerminalKind::LINE_2:
            case TerminalKind::TRAFO_1:
            case TerminalKind::TRAFO_2:
            case TerminalKind::HVDC_1:
            case TerminalKind::HVDC_2:
                ++n_branch[root];
                break;
            case TerminalKind::LOAD:
            case TerminalKind::GEN:
            case TerminalKind::SGEN:
            case TerminalKind::STORAGE:
            case TerminalKind::SHUNT:
            case TerminalKind::SVC:
                break;
        }
    }

    // 3. number the valid components: busbar-section holders first, in busbar-
    //    section order, then the rest by lowest node. Into scratch, so that a
    //    component count the layout cannot hold leaves the current labels alone.
    std::vector<int> root_bus(n, -1);
    int nb_buses = 0;
    const auto number_root = [&](int root){
        const std::size_t r = static_cast<std::size_t>(root);
        if(root_bus[r] != -1) return;
        const bool valid = ((n_bbs[r] >= 1) && (n_feeder[r] >= 1)) ||
                           ((n_branch[r] >= 1) && (n_feeder[r] >= 2));
        if(!valid) return;
        root_bus[r] = ++nb_buses;  // LocalBusId is 1-based
    };
    for(Eigen::Index i = 0; i < bbs_node_.size(); ++i) number_root(_find(bbs_node_(i)));
    for(std::size_t node = 0; node < n; ++node) number_root(_find(static_cast<int>(node)));

    if(nb_buses > nmax_busbar_per_sub){
        std::ostringstream exc_;
        exc_ << sub_label(sub_id) << "::label: the switch positions make " << nb_buses
             << " electrical buses, but the grid's bus layout holds at most " << nmax_busbar_per_sub
             << " per substation (nmax_busbar_per_sub). Declare the substations with a larger "
             << "capacity (init_bus / n_busbar_per_sub).";
        throw std::runtime_error(exc_.str());
    }

    // 4. commit
    node_bus_.resize(nb_nodes_);
    for(std::size_t node = 0; node < n; ++node){
        node_bus_(static_cast<Eigen::Index>(node)) = root_bus[static_cast<std::size_t>(_find(static_cast<int>(node)))];
    }
    nb_buses_ = nb_buses;
    labels_ready_ = true;
}

int SubstationTopology::node_bus(int node) const
{
    if((node < 0) || (node >= nb_nodes_)){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::node_bus: node " << node << " is out of range [0, "
             << nb_nodes_ << ").";
        throw std::out_of_range(exc_.str());
    }
    if(!labels_ready_){
        throw std::runtime_error("SubstationTopology::node_bus: the labels have not been computed for "
                                 "the current switch positions (see label()).");
    }
    return node_bus_(node);
}

int SubstationTopology::_checked_bbs_id(int bbs_id, const char * fun_name) const
{
    if((bbs_id < 0) || (bbs_id >= nb_busbar_sections())){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::" << fun_name << ": busbar section id " << bbs_id
             << " is out of range [0, " << nb_busbar_sections() << ").";
        throw std::out_of_range(exc_.str());
    }
    return bbs_id;
}

int SubstationTopology::_checked_sw_id(int sw_id, const char * fun_name) const
{
    if((sw_id < 0) || (sw_id >= nb_switches())){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::" << fun_name << ": switch id " << sw_id
             << " is out of range [0, " << nb_switches() << ").";
        throw std::out_of_range(exc_.str());
    }
    return sw_id;
}

std::string SubstationTopology::bbs_name(int bbs_id) const
{
    _checked_bbs_id(bbs_id, "bbs_name");
    if(bbs_names_.empty()) return "";
    return bbs_names_[static_cast<std::size_t>(bbs_id)];
}

std::string SubstationTopology::sw_name(int sw_id) const
{
    _checked_sw_id(sw_id, "sw_name");
    if(sw_names_.empty()) return "";
    return sw_names_[static_cast<std::size_t>(sw_id)];
}

void SubstationTopology::set_bbs_names(const std::vector<std::string> & names)
{
    if(static_cast<int>(names.size()) != nb_busbar_sections()){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::set_bbs_names: " << names.size() << " names for "
             << nb_busbar_sections() << " busbar sections.";
        throw std::runtime_error(exc_.str());
    }
    bbs_names_ = names;
}

void SubstationTopology::set_sw_names(const std::vector<std::string> & names)
{
    if(static_cast<int>(names.size()) != nb_switches()){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::set_sw_names: " << names.size() << " names for "
             << nb_switches() << " switches.";
        throw std::runtime_error(exc_.str());
    }
    sw_names_ = names;
}

void SubstationTopology::add_terminal(TerminalKind kind, int el_id, int node)
{
    if((node < 0) || (node >= nb_nodes_)){
        std::ostringstream exc_;
        exc_ << "SubstationTopology::add_terminal: element " << el_id << " (terminal kind "
             << static_cast<int>(kind) << ") stands on node " << node << ", out of range [0, "
             << nb_nodes_ << ").";
        throw std::out_of_range(exc_.str());
    }
    terminals_.push_back(Terminal{kind, el_id, node});
    _invalidate_labels();
}

SubstationTopology::StateRes SubstationTopology::get_state() const
{
    std::vector<int> bbs_node(bbs_node_.begin(), bbs_node_.end());
    std::vector<int> sw_node1(sw_node1_.begin(), sw_node1_.end());
    std::vector<int> sw_node2(sw_node2_.begin(), sw_node2_.end());
    std::vector<int> sw_kind;
    sw_kind.reserve(sw_kind_.size());
    for(const SwitchKind kind : sw_kind_) sw_kind.push_back(static_cast<int>(kind));
    return StateRes(nb_nodes_, bbs_node, bbs_names_, sw_node1, sw_node2, sw_kind,
                    sw_open_, sw_retained_, sw_names_);
}

void SubstationTopology::set_state(StateRes & my_state)
{
    // every field comes from a pickle or a binary file: validate everything
    // before assigning anything
    const int nb_nodes = std::get<NB_NODES>(my_state);
    const std::vector<int> & bbs_node_v = std::get<BBS_NODE>(my_state);
    const std::vector<std::string> & bbs_names = std::get<BBS_NAMES>(my_state);
    const std::vector<int> & sw_node1_v = std::get<SW_NODE1>(my_state);
    const std::vector<int> & sw_node2_v = std::get<SW_NODE2>(my_state);
    const std::vector<int> & sw_kind_v = std::get<SW_KIND>(my_state);
    const std::vector<bool> & sw_open = std::get<SW_OPEN>(my_state);
    const std::vector<bool> & sw_retained = std::get<SW_RETAINED>(my_state);
    const std::vector<std::string> & sw_names = std::get<SW_NAMES>(my_state);

    const IntVect bbs_node = bbs_node_v.empty() ? IntVect() : IntVect::Map(bbs_node_v.data(), bbs_node_v.size());
    const IntVect sw_node1 = sw_node1_v.empty() ? IntVect() : IntVect::Map(sw_node1_v.data(), sw_node1_v.size());
    const IntVect sw_node2 = sw_node2_v.empty() ? IntVect() : IntVect::Map(sw_node2_v.data(), sw_node2_v.size());
    std::vector<SwitchKind> sw_kind;
    sw_kind.reserve(sw_kind_v.size());
    // the cast is what _check validates (an unknown integer becomes an unknown
    // kind, which it rejects), so nothing is trusted here either
    for(const int kind : sw_kind_v) sw_kind.push_back(static_cast<SwitchKind>(kind));

    _check(nb_nodes, bbs_node, bbs_names, sw_node1, sw_node2, sw_kind, sw_open, sw_retained, sw_names, -1);

    nb_nodes_ = nb_nodes;
    bbs_node_ = bbs_node;
    bbs_names_ = bbs_names;
    sw_node1_ = sw_node1;
    sw_node2_ = sw_node2;
    sw_kind_ = sw_kind;
    sw_open_ = sw_open;
    sw_retained_ = sw_retained;
    sw_names_ = sw_names;
    terminals_.clear();  // derived: whoever restores the elements rebuilds it
    node_bus_ = IntVect::Constant(nb_nodes_, -1);
    nb_buses_ = 0;
    _invalidate_labels();
}

void SubstationTopology::check_valid(int sub_id) const
{
    _check(nb_nodes_, bbs_node_, bbs_names_, sw_node1_, sw_node2_, sw_kind_, sw_open_, sw_retained_, sw_names_, sub_id);
    // the terminals are derived, but they are indexed against nb_nodes_ by
    // whoever consumes them: an entry out of range would be a bug in the
    // rebuild, and cheap to refuse here
    for(const Terminal & term : terminals_){
        if((term.node < 0) || (term.node >= nb_nodes_)){
            std::ostringstream exc_;
            exc_ << sub_label(sub_id) << ": a terminal (element " << term.el_id << ", kind "
                 << static_cast<int>(term.kind) << ") stands on node " << term.node
                 << ", out of range [0, " << nb_nodes_ << ").";
            throw std::out_of_range(exc_.str());
        }
    }
}

}  // namespace ls2g
