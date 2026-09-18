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
