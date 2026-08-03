// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "StorageContainer.hpp"
#include "BinaryArchive.hpp"

#include <sstream>
#include <iostream>

namespace ls2g {

StorageContainer::StateRes StorageContainer::get_state() const
{
    const auto tmp = get_osc_pq_state();  // osc : one side container
    StorageContainer::StateRes res(tmp);
    return res;
}

void StorageContainer::set_state(StorageContainer::StateRes & my_state)
{
    set_osc_pq_state(std::get<0>(my_state));  // osc : one side container
    reset_results();
}

void StorageContainer::fillSbus(Eigen::Ref<CplxVect> Sbus,
                                const SolverBusIdVect & id_grid_to_solver,
                                bool /*ac*/) const
{
    int nb_storage = nb();
    GlobalBusId bus_id_me;
    SolverBusId bus_id_solver;
    cplx_type tmp;
    for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
        //  i don't do anything if the storage is disconnected
        if(!status_[storage_id]) continue;

        bus_id_me = bus_id_(storage_id);
        if(bus_id_me.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "StorageContainer::fillSbus: the storage with id ";
            exc_ << storage_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        bus_id_solver = id_grid_to_solver[bus_id_me.cast_int()];
        if(bus_id_solver.cast_int() == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "StorageContainer::fillSbus: the storage with id ";
            exc_ << storage_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        // load convention: positive target_p means power drawn from the grid
        tmp = {target_p_mw_(storage_id), target_q_mvar_(storage_id)};
        Sbus.coeffRef(bus_id_solver.cast_int()) -= tmp;
    }
}

void StorageContainer::save_binary(const std::string & path, bool atomic) const {
    ls2g::save_binary_generic(*this, path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR, atomic);
}

StorageContainer StorageContainer::load_binary(const std::string & path) {
    return ls2g::load_binary_generic<StorageContainer>(path, VERSION_MAJOR, VERSION_MEDIUM, VERSION_MINOR);
}

} // namespace ls2g
