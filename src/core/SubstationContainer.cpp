// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SubstationContainer.hpp"
#include "BinaryArchive.hpp"

#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>

namespace ls2g {

SubstationContainer::StateRes SubstationContainer::get_state() const
{
     std::vector<real_type> sub_vn_kv(sub_vn_kv_.begin(), sub_vn_kv_.end());
     std::vector<real_type> bus_vn_kv(bus_vn_kv_.begin(), bus_vn_kv_.end());
     std::vector<real_type> bus_vmin_kv(bus_vmin_kv_.begin(), bus_vmin_kv_.end());
     std::vector<real_type> bus_vmax_kv(bus_vmax_kv_.begin(), bus_vmax_kv_.end());
     std::vector<SubstationTopology::StateRes> topologies;
     topologies.reserve(topologies_.size());
     for(const SubstationTopology & topo : topologies_) topologies.push_back(topo.get_state());
     SubstationContainer::StateRes res(
        n_sub_,
        nmax_busbar_per_sub_,
        sub_vn_kv,
        bus_vn_kv,
        sub_names_,
        bus_vmin_kv,
        bus_vmax_kv,
        topologies);
     return res;
}

void SubstationContainer::set_state(SubstationContainer::StateRes & my_state)
{
    // NB every field below comes straight from a pickle or a binary file, ie from
    // outside the C++ core. This container is the *root* of the grid's index space:
    // `nb_bus()` (= bus_vn_kv_.size()) is the upper bound LSGrid::check_grid()
    // validates every element's bus id against, it is what sizes
    // `nb_elements_per_bus_` (which `is_bus_connected()` indexes with those same
    // ids), and `n_sub_` is both the substation-id bound and a modulo divisor
    // (sub_id_of_bus). Nothing downstream re-derives these from each other, so if
    // they are allowed to disagree here, a *validated* bus id still reads and writes
    // past the end of the count vector -- release wheels are built -O3 -DNDEBUG, so
    // neither Eigen nor the standard library catches it. Validate them all now,
    // before anything is assigned.
    const int n_sub = std::get<0>(my_state);
    const int nmax_busbar_per_sub = std::get<1>(my_state);
    std::vector<real_type> & sub_vn_kv = std::get<2>(my_state);
    std::vector<real_type> & bus_vn_kv = std::get<3>(my_state);
    const std::vector<std::string> & sub_names = std::get<4>(my_state);
    std::vector<real_type> & bus_vmin_kv = std::get<5>(my_state);
    std::vector<real_type> & bus_vmax_kv = std::get<6>(my_state);
    std::vector<SubstationTopology::StateRes> & topo_states = std::get<7>(my_state);

    // a default-constructed container has n_sub_ == nmax_busbar_per_sub_ == -1 and
    // no bus at all: that state must still round-trip, so it is the one case where
    // negative counts are legal.
    const bool is_empty = bus_vn_kv.empty();
    // positivity + no int overflow on n_sub * nmax_busbar_per_sub, the same guarantee
    // init_bus() / init_sub() get from the same helper.
    const std::int64_t n_bus_max = is_empty
        ? std::int64_t{0}
        : static_cast<std::int64_t>(checked_nb_bus(n_sub, nmax_busbar_per_sub,
                                                   "SubstationContainer::set_state"));

    // check sizes
    // bus_vn_kv_ defines nb_bus(), the bound every other index is checked against.
    const auto check_len = [](std::size_t actual, std::int64_t expected, const char * name){
        if(static_cast<std::int64_t>(actual) != expected){
            std::ostringstream exc_;
            exc_ << "SubstationContainer::set_state: '" << name << "' has " << actual
                 << " elements but " << expected << " were expected. The grid state is inconsistent "
                 << "(this field is indexed with ids validated against another one, so a mismatch "
                 << "would cause an out-of-bounds access).";
            throw std::runtime_error(exc_.str());
        }
    };
    if(!is_empty) check_len(bus_vn_kv.size(), n_bus_max, "bus_vn_kv");
    // Optional fields: either absent, or fully sized. sub_vn_kv_ and sub_names_ are
    // only filled by init_sub() / init_sub_names(), which the usual init_bus() path
    // does not call -- so "empty" is the normal state for most grids, NOT a defect.
    if(!sub_vn_kv.empty()) check_len(sub_vn_kv.size(), n_sub, "sub_vn_kv");
    if(!sub_names.empty()) check_len(sub_names.size(), n_sub, "sub_names");
    if(!bus_vmin_kv.empty()) check_len(bus_vmin_kv.size(), static_cast<std::int64_t>(bus_vn_kv.size()), "bus_vmin_kv");
    if(!bus_vmax_kv.empty()) check_len(bus_vmax_kv.size(), static_cast<std::int64_t>(bus_vn_kv.size()), "bus_vmax_kv");
    // The detailed topology: absent, or one entry per substation. Each entry
    // validates itself (SubstationTopology::set_state) before anything is
    // assigned here, so a poisoned switch table is refused with the rest.
    if(!topo_states.empty()) check_len(topo_states.size(), n_sub, "detailed topology");
    std::vector<SubstationTopology> topologies(topo_states.size());
    for(std::size_t sub_id = 0; sub_id < topo_states.size(); ++sub_id){
        topologies[sub_id].set_state(topo_states[sub_id]);
    }

    // assign data (nothing above has modified `this`, so a rejected state leaves
    // the container untouched)
    n_sub_ = n_sub;
    nmax_busbar_per_sub_ = nmax_busbar_per_sub;
    n_bus_max_ = static_cast<int>(n_bus_max);
    sub_vn_kv_ = RealVect::Map(sub_vn_kv.data(), sub_vn_kv.size());
    bus_vn_kv_ = RealVect::Map(bus_vn_kv.data(), bus_vn_kv.size());
    // which buses are in the solved system is not restored, it is recounted: the
    // elements' own status IS in the file, and they are what holds a bus. Sized to
    // the grid we just restored (so it must come after bus_vn_kv_) and left zeroed
    // and disarmed; LSGrid counts from the elements before anything reads it.
    reset_bus_element_counts();
    sub_names_ = sub_names;
    bus_vmin_kv_ = bus_vmin_kv.empty() ? RealVect() : RealVect::Map(bus_vmin_kv.data(), bus_vmin_kv.size());
    bus_vmax_kv_ = bus_vmax_kv.empty() ? RealVect() : RealVect::Map(bus_vmax_kv.data(), bus_vmax_kv.size());
    topologies_ = std::move(topologies);
    _rebuild_topology_offsets();
}

// ---- detailed topology ------------------------------------------------------

void SubstationContainer::init_detailed_topology(std::vector<SubstationTopology> topologies)
{
    if(topologies.empty()){
        clear_detailed_topology();
        return;
    }
    if(n_sub_ <= 0){
        throw std::runtime_error("SubstationContainer::init_detailed_topology: the substations must be "
                                 "declared first (see init_bus).");
    }
    if(static_cast<int>(topologies.size()) != n_sub_){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::init_detailed_topology: " << topologies.size()
             << " substation topologies for " << n_sub_ << " substations (one per substation, "
             << "possibly with zero nodes, is required).";
        throw std::runtime_error(exc_.str());
    }
    for(std::size_t sub_id = 0; sub_id < topologies.size(); ++sub_id){
        topologies[sub_id].check_valid(static_cast<int>(sub_id));
    }
    topologies_ = std::move(topologies);
    _rebuild_topology_offsets();
}

void SubstationContainer::clear_detailed_topology()
{
    topologies_.clear();
    _rebuild_topology_offsets();
}

void SubstationContainer::_rebuild_topology_offsets()
{
    node_offset_.clear();
    switch_offset_.clear();
    bbs_offset_.clear();
    switch_sub_.clear();
    bbs_sub_.clear();
    if(topologies_.empty()) return;
    const std::size_t n_sub = topologies_.size();
    node_offset_.assign(n_sub + 1, 0);
    switch_offset_.assign(n_sub + 1, 0);
    bbs_offset_.assign(n_sub + 1, 0);
    for(std::size_t sub_id = 0; sub_id < n_sub; ++sub_id){
        const SubstationTopology & topo = topologies_[sub_id];
        node_offset_[sub_id + 1] = node_offset_[sub_id] + topo.nb_nodes();
        switch_offset_[sub_id + 1] = switch_offset_[sub_id] + topo.nb_switches();
        bbs_offset_[sub_id + 1] = bbs_offset_[sub_id] + topo.nb_busbar_sections();
    }
    switch_sub_.reserve(static_cast<std::size_t>(switch_offset_.back()));
    bbs_sub_.reserve(static_cast<std::size_t>(bbs_offset_.back()));
    for(std::size_t sub_id = 0; sub_id < n_sub; ++sub_id){
        const SubstationTopology & topo = topologies_[sub_id];
        switch_sub_.insert(switch_sub_.end(), static_cast<std::size_t>(topo.nb_switches()), static_cast<int>(sub_id));
        bbs_sub_.insert(bbs_sub_.end(), static_cast<std::size_t>(topo.nb_busbar_sections()), static_cast<int>(sub_id));
    }
}

int SubstationContainer::_checked_sub_id(int sub_id, const char * fun_name) const
{
    if(topologies_.empty()){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::" << fun_name << ": this grid has no detailed topology "
             << "(see init_detailed_topology).";
        throw std::runtime_error(exc_.str());
    }
    if((sub_id < 0) || (sub_id >= static_cast<int>(topologies_.size()))){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::" << fun_name << ": substation id " << sub_id
             << " is out of range [0, " << topologies_.size() << ").";
        throw std::out_of_range(exc_.str());
    }
    return sub_id;
}

int SubstationContainer::_checked_switch_id(int switch_id, const char * fun_name) const
{
    if((switch_id < 0) || (switch_id >= nb_switches())){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::" << fun_name << ": switch id " << switch_id
             << " is out of range [0, " << nb_switches() << ").";
        throw std::out_of_range(exc_.str());
    }
    return switch_id;
}

int SubstationContainer::_checked_bbs_id(int bbs_id, const char * fun_name) const
{
    if((bbs_id < 0) || (bbs_id >= nb_busbar_sections())){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::" << fun_name << ": busbar section id " << bbs_id
             << " is out of range [0, " << nb_busbar_sections() << ").";
        throw std::out_of_range(exc_.str());
    }
    return bbs_id;
}

const SubstationTopology & SubstationContainer::topology(int sub_id) const
{
    return topologies_[static_cast<std::size_t>(_checked_sub_id(sub_id, "topology"))];
}

SubstationTopology & SubstationContainer::topology(int sub_id)
{
    return topologies_[static_cast<std::size_t>(_checked_sub_id(sub_id, "topology"))];
}

void SubstationContainer::set_switch_names(const std::vector<std::string> & names)
{
    if(static_cast<int>(names.size()) != nb_switches()){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::set_switch_names: " << names.size() << " names for "
             << nb_switches() << " switches.";
        throw std::runtime_error(exc_.str());
    }
    for(std::size_t sub_id = 0; sub_id < topologies_.size(); ++sub_id){
        const auto first = names.begin() + switch_offset_[sub_id];
        const auto last = names.begin() + switch_offset_[sub_id + 1];
        topologies_[sub_id].set_sw_names(std::vector<std::string>(first, last));
    }
}

void SubstationContainer::label_all()
{
    for(std::size_t sub_id = 0; sub_id < topologies_.size(); ++sub_id){
        topologies_[sub_id].label(nmax_busbar_per_sub_, static_cast<int>(sub_id));
    }
}

IntVect SubstationContainer::get_node_bus() const
{
    IntVect res = IntVect::Constant(nb_nodes(), BaseConstants::_deactivated_bus_id);
    for(std::size_t sub_id = 0; sub_id < topologies_.size(); ++sub_id){
        const SubstationTopology & topo = topologies_[sub_id];
        if(!topo.labels_ready()){
            std::ostringstream exc_;
            exc_ << "SubstationContainer::get_node_bus: the labels of substation " << sub_id
                 << " have not been computed for its current switch positions (see label_all).";
            throw std::runtime_error(exc_.str());
        }
        const IntVect & local = topo.node_bus();
        const int first = node_offset_[sub_id];
        for(Eigen::Index node = 0; node < local.size(); ++node){
            if(local(node) == BaseConstants::_deactivated_bus_id) continue;
            res(first + node) = local_to_gridmodel(static_cast<int>(sub_id), LocalBusId(local(node))).cast_int();
        }
    }
    return res;
}

void SubstationContainer::set_busbar_section_names(const std::vector<std::string> & names)
{
    if(static_cast<int>(names.size()) != nb_busbar_sections()){
        std::ostringstream exc_;
        exc_ << "SubstationContainer::set_busbar_section_names: " << names.size() << " names for "
             << nb_busbar_sections() << " busbar sections.";
        throw std::runtime_error(exc_.str());
    }
    for(std::size_t sub_id = 0; sub_id < topologies_.size(); ++sub_id){
        const auto first = names.begin() + bbs_offset_[sub_id];
        const auto last = names.begin() + bbs_offset_[sub_id + 1];
        topologies_[sub_id].set_bbs_names(std::vector<std::string>(first, last));
    }
}

void SubstationContainer::check_valid() const
{
    // a default-constructed / never-initialized container is consistent by
    // definition: it has no bus at all, so nothing can index into it.
    if(bus_vn_kv_.size() == 0) return;

    if(n_sub_ <= 0){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the grid declares " << n_sub_ << " substation(s) but has "
             << bus_vn_kv_.size() << " buses. The number of substations must be strictly positive "
             << "(it is also used as a modulo divisor in sub_id_of_bus).";
        throw std::runtime_error(exc_.str());
    }
    if(nmax_busbar_per_sub_ <= 0){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the grid declares " << nmax_busbar_per_sub_
             << " busbar(s) per substation; it must be strictly positive.";
        throw std::runtime_error(exc_.str());
    }
    // nb_elements_per_bus_ is indexed with the very bus ids check_grid() validates
    // against nb_bus() (= bus_vn_kv_.size()): the two must have exactly the same
    // length. It is derived, not loaded, so this cannot be violated by a crafted
    // file -- it is here to catch a container mutated into an inconsistent state
    // some other way.
    if(nb_elements_per_bus_.size() != static_cast<std::size_t>(bus_vn_kv_.size())){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the per-bus element count vector has "
             << nb_elements_per_bus_.size()
             << " entries while the grid has " << bus_vn_kv_.size() << " buses. Both are indexed "
             << "by bus id and must have the same length.";
        throw std::runtime_error(exc_.str());
    }
    if(static_cast<std::int64_t>(bus_vn_kv_.size()) !=
       static_cast<std::int64_t>(n_sub_) * static_cast<std::int64_t>(nmax_busbar_per_sub_)){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the grid has " << bus_vn_kv_.size() << " buses but declares "
             << n_sub_ << " substations of at most " << nmax_busbar_per_sub_ << " busbars ("
             << static_cast<std::int64_t>(n_sub_) * static_cast<std::int64_t>(nmax_busbar_per_sub_)
             << " buses). Bus ids are laid out as `sub_id + (busbar - 1) * n_sub`, so both must match.";
        throw std::runtime_error(exc_.str());
    }
    // optional fields: either absent, or one entry per substation / per bus.
    // sub_vn_kv_ / sub_names_ are only filled by init_sub() / init_sub_names(),
    // which the usual init_bus() path does not call: empty is the normal state.
    if((sub_vn_kv_.size() != 0) && (sub_vn_kv_.size() != n_sub_)){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the substation nominal-voltage vector has " << sub_vn_kv_.size()
             << " entries for " << n_sub_ << " substations (it must be either empty or complete).";
        throw std::runtime_error(exc_.str());
    }
    if(!sub_names_.empty() && (sub_names_.size() != static_cast<std::size_t>(n_sub_))){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the substation names vector has " << sub_names_.size()
             << " entries for " << n_sub_ << " substations (it must be either empty or complete).";
        throw std::runtime_error(exc_.str());
    }
    if((bus_vmin_kv_.size() != 0) && (bus_vmin_kv_.size() != bus_vn_kv_.size())){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the per-bus min voltage vector has " << bus_vmin_kv_.size()
             << " entries for " << bus_vn_kv_.size() << " buses (it must be either empty or complete).";
        throw std::runtime_error(exc_.str());
    }
    if((bus_vmax_kv_.size() != 0) && (bus_vmax_kv_.size() != bus_vn_kv_.size())){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the per-bus max voltage vector has " << bus_vmax_kv_.size()
             << " entries for " << bus_vn_kv_.size() << " buses (it must be either empty or complete).";
        throw std::runtime_error(exc_.str());
    }
    // vmin_ and vmax_ are consumed together, one indexed by the other's length
    // (ContingencyAnalysis::check_bus_voltage_violations loops to vmin.size() then
    // reads vmax(grid_id)). Each is individually allowed to be empty-or-complete
    // above, but one present while the other is empty makes that loop read past the
    // end of the empty one. They must have the same presence.
    if((bus_vmin_kv_.size() == 0) != (bus_vmax_kv_.size() == 0)){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the per-bus min and max voltage vectors must both be set or "
             << "both be empty (min has " << bus_vmin_kv_.size() << " entries, max has "
             << bus_vmax_kv_.size() << "). They are consumed together, one indexed by the other's "
             << "length, so a mismatch would cause an out-of-bounds read.";
        throw std::runtime_error(exc_.str());
    }

    // the detailed topology: absent, or one entry per substation, each consistent
    // on its own, with the derived grid-wide numbering sized for it
    if(!topologies_.empty()){
        if(topologies_.size() != static_cast<std::size_t>(n_sub_)){
            std::ostringstream exc_;
            exc_ << "LSGrid::check_grid: the detailed topology describes " << topologies_.size()
                 << " substations while the grid has " << n_sub_ << " (it must describe all of them "
                 << "or none).";
            throw std::runtime_error(exc_.str());
        }
        for(std::size_t sub_id = 0; sub_id < topologies_.size(); ++sub_id){
            topologies_[sub_id].check_valid(static_cast<int>(sub_id));
        }
        const std::size_t expected = topologies_.size() + 1;
        if((node_offset_.size() != expected) || (switch_offset_.size() != expected) ||
           (bbs_offset_.size() != expected) ||
           (switch_sub_.size() != static_cast<std::size_t>(switch_offset_.back())) ||
           (bbs_sub_.size() != static_cast<std::size_t>(bbs_offset_.back()))){
            throw std::runtime_error("LSGrid::check_grid: the grid-wide numbering of the detailed "
                                     "topology is out of step with the substations (a bug: it is "
                                     "derived, see SubstationContainer::_rebuild_topology_offsets).");
        }
    }

    // every nominal voltage must be a finite, strictly positive number
    check_vn_kv_positive();

    // All busbars of a substation share its nominal voltage. init_bus() is supposed
    // to enforce this at build time, but set_state() bypasses init_bus() entirely,
    // so a pickle / binary file could always carry a grid violating it. Re-check
    // here, which is the point of check_grid().
    check_bus_vn_kv_uniform_per_sub();

    // sub_vn_kv_ is optional (see the note on StateRes), but when it IS present it
    // is by definition the common nominal voltage of each substation's busbars,
    // i.e. bus_vn_kv_ restricted to the first n_sub_ bus ids (busbar 1 of each
    // substation). A stored value disagreeing with that is an inconsistent state.
    if(sub_vn_kv_.size() != 0){
        for(int sub_id = 0; sub_id < n_sub_; ++sub_id){
            if(std::abs(sub_vn_kv_(sub_id) - bus_vn_kv_(sub_id)) > BaseConstants::_tol_equal_float){
                std::ostringstream exc_;
                exc_ << "LSGrid::check_grid: substation " << sub_id << " declares a nominal voltage of "
                     << sub_vn_kv_(sub_id) << " kV but its buses carry " << bus_vn_kv_(sub_id)
                     << " kV. Both must agree.";
                throw std::runtime_error(exc_.str());
            }
        }
    }
}

void SubstationContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

SubstationContainer SubstationContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<SubstationContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
