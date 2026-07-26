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
     SubstationContainer::StateRes res(
        n_sub_,
        nmax_busbar_per_sub_,
        sub_vn_kv,
        bus_status_,
        bus_vn_kv,
        sub_names_,
        bus_vmin_kv,
        bus_vmax_kv);
     return res;
}

void SubstationContainer::set_state(SubstationContainer::StateRes & my_state)
{
    // NB every field below comes straight from a pickle or a binary file, ie from
    // outside the C++ core. This container is the *root* of the grid's index space:
    // `nb_bus()` (= bus_vn_kv_.size()) is the upper bound LSGrid::check_grid()
    // validates every element's bus id against, while `bus_status_` is what
    // `is_bus_connected()` / `disconnect_all_buses()` actually index with those same
    // ids, and `n_sub_` is both the substation-id bound and a modulo divisor
    // (sub_id_of_bus). Nothing downstream re-derives these from each other, so if
    // they are allowed to disagree here, a *validated* bus id still reads and writes
    // past the end of bus_status_ -- release wheels are built -O3 -DNDEBUG, so
    // neither Eigen nor the standard library catches it. Validate them all now,
    // before anything is assigned.
    const int n_sub = std::get<0>(my_state);
    const int nmax_busbar_per_sub = std::get<1>(my_state);
    std::vector<real_type> & sub_vn_kv = std::get<2>(my_state);
    std::vector<bool> & bus_status = std::get<3>(my_state);
    std::vector<real_type> & bus_vn_kv = std::get<4>(my_state);
    const std::vector<std::string> & sub_names = std::get<5>(my_state);
    std::vector<real_type> & bus_vmin_kv = std::get<6>(my_state);
    std::vector<real_type> & bus_vmax_kv = std::get<7>(my_state);

    // a default-constructed container has n_sub_ == nmax_busbar_per_sub_ == -1 and
    // no bus at all: that state must still round-trip, so it is the one case where
    // negative counts are legal.
    const bool is_empty = bus_vn_kv.empty() && bus_status.empty();
    if(!(is_empty && (n_sub <= 0) && (nmax_busbar_per_sub <= 0)))
    {
        if(n_sub <= 0){
            std::ostringstream exc_;
            exc_ << "SubstationContainer::set_state: the number of substations must be strictly "
                 << "positive on a grid that has buses, got " << n_sub << " (it is also used as a "
                 << "modulo divisor in sub_id_of_bus, so 0 would be a division by zero).";
            throw std::runtime_error(exc_.str());
        }
        if(nmax_busbar_per_sub <= 0){
            std::ostringstream exc_;
            exc_ << "SubstationContainer::set_state: the maximum number of busbars per substation "
                 << "must be strictly positive, got " << nmax_busbar_per_sub << ".";
            throw std::runtime_error(exc_.str());
        }
        // n_sub_ * nmax_busbar_per_sub_ is stored in an int: reject a product that
        // does not fit rather than letting it overflow (UB, and a negative bus count).
        if(nmax_busbar_per_sub > std::numeric_limits<int>::max() / n_sub){
            std::ostringstream exc_;
            exc_ << "SubstationContainer::set_state: n_sub (" << n_sub << ") * nmax_busbar_per_sub ("
                 << nmax_busbar_per_sub << ") overflows the maximum number of buses representable.";
            throw std::runtime_error(exc_.str());
        }
    }
    const std::int64_t n_bus_max = static_cast<std::int64_t>(n_sub) * static_cast<std::int64_t>(nmax_busbar_per_sub);

    // bus_vn_kv_ defines nb_bus(), the bound every other index is checked against;
    // bus_status_ is indexed with those same ids and MUST match it exactly.
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
    check_len(bus_status.size(), static_cast<std::int64_t>(bus_vn_kv.size()), "bus_status");
    // Optional fields: either absent, or fully sized. sub_vn_kv_ and sub_names_ are
    // only filled by init_sub() / init_sub_names(), which the usual init_bus() path
    // does not call -- so "empty" is the normal state for most grids, NOT a defect.
    if(!sub_vn_kv.empty()) check_len(sub_vn_kv.size(), n_sub, "sub_vn_kv");
    if(!sub_names.empty()) check_len(sub_names.size(), n_sub, "sub_names");
    if(!bus_vmin_kv.empty()) check_len(bus_vmin_kv.size(), static_cast<std::int64_t>(bus_vn_kv.size()), "bus_vmin_kv");
    if(!bus_vmax_kv.empty()) check_len(bus_vmax_kv.size(), static_cast<std::int64_t>(bus_vn_kv.size()), "bus_vmax_kv");

    // assign data (nothing above has modified `this`, so a rejected state leaves
    // the container untouched)
    n_sub_ = n_sub;
    nmax_busbar_per_sub_ = nmax_busbar_per_sub;
    n_bus_max_ = static_cast<int>(n_bus_max);
    sub_vn_kv_ = RealVect::Map(sub_vn_kv.data(), sub_vn_kv.size());
    bus_status_ = bus_status;
    bus_vn_kv_ = RealVect::Map(bus_vn_kv.data(), bus_vn_kv.size());
    sub_names_ = sub_names;
    bus_vmin_kv_ = bus_vmin_kv.empty() ? RealVect() : RealVect::Map(bus_vmin_kv.data(), bus_vmin_kv.size());
    bus_vmax_kv_ = bus_vmax_kv.empty() ? RealVect() : RealVect::Map(bus_vmax_kv.data(), bus_vmax_kv.size());
}

void SubstationContainer::check_valid() const
{
    // a default-constructed / never-initialized container is consistent by
    // definition: it has no bus at all, so nothing can index into it.
    if((bus_vn_kv_.size() == 0) && bus_status_.empty()) return;

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
    // bus_status_ is indexed with the very bus ids check_grid() validates against
    // nb_bus() (= bus_vn_kv_.size()), and disconnect_all_buses() writes nb_bus()
    // entries into it: the two must have exactly the same length.
    if(bus_status_.size() != static_cast<std::size_t>(bus_vn_kv_.size())){
        std::ostringstream exc_;
        exc_ << "LSGrid::check_grid: the bus status vector has " << bus_status_.size()
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
