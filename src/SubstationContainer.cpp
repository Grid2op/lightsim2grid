// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SubstationContainer.hpp"

#include <iostream>
#include <sstream>

SubstationInfo::SubstationInfo(const SubstationContainer & r_data, int my_id):
    id(my_id),
    name(""),
    nb_max_busbars(-1),
    vn_kv(-1.)
{
    if(my_id < 0) return;
    if(my_id > r_data.nb()) return;
    name = r_data.sub_names_[my_id];
    nb_max_busbars = r_data.nmax_busbar_per_sub_;
    vn_kv = r_data.bus_vn_kv_[my_id];
}

SubstationContainer::StateRes SubstationContainer::get_state() const
{
     std::vector<real_type> sub_vn_kv(sub_vn_kv_.begin(), sub_vn_kv_.end());
     std::vector<real_type> bus_vn_kv(bus_vn_kv_.begin(), bus_vn_kv_.end());
     SubstationContainer::StateRes res(
        n_sub_,
        nmax_busbar_per_sub_,
        sub_vn_kv,
        bus_status_,
        bus_vn_kv,
        sub_names_);
     return res;
}

void SubstationContainer::set_state(SubstationContainer::StateRes & my_state)
{
    n_sub_ = std::get<0>(my_state);
    nmax_busbar_per_sub_ = std::get<1>(my_state);
    n_bus_max_ = n_sub_ * nmax_busbar_per_sub_;

    // the generators themelves
    std::vector<real_type> & sub_vn_kv = std::get<2>(my_state);
    std::vector<bool> & bus_status = std::get<3>(my_state);
    std::vector<real_type> & bus_vn_kv = std::get<4>(my_state);

    // check sizes
    // TODO dev switches

    // assign data
    sub_vn_kv_ = RealVect::Map(&sub_vn_kv[0], sub_vn_kv.size());
    bus_status_ = bus_status;
    bus_vn_kv_ = RealVect::Map(&bus_vn_kv[0], bus_vn_kv.size());
    sub_names_ = std::get<5>(my_state);
}
