// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LoadContainer.hpp"

#include <sstream>
#include <iostream>

LoadInfo::LoadInfo(const LoadContainer & r_data_load, int my_id): 
        OneSidePQInfo(r_data_load, my_id) {}

LoadContainer::StateRes LoadContainer::get_state() const
{
    const auto tmp = get_osc_pq_state();  // osc : one side container
    LoadContainer::StateRes res(tmp);
    return res;
}

void LoadContainer::set_state(LoadContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<0>(my_state));  // osc : one side container
    reset_results();
}

void LoadContainer::fillSbus(CplxVect & Sbus,
                             const std::vector<int> & id_grid_to_solver,
                             bool ac) const
{
    int nb_load = nb();
    int bus_id_me, bus_id_solver;
    cplx_type tmp;
    for(int load_id = 0; load_id < nb_load; ++load_id){
        //  i don't do anything if the load is disconnected
        if(!status_[load_id]) continue;

        bus_id_me = bus_id_(load_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "LoadContainer::fillSbus: the load with id ";
            exc_ << load_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        tmp = {target_p_mw_(load_id), target_q_mvar_(load_id)};
        Sbus.coeffRef(bus_id_solver) -= tmp;
    }
}
