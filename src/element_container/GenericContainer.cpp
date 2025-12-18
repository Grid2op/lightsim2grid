// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GenericContainer.hpp"

#include <iostream>
#include <sstream>

const int GenericContainer::_deactivated_bus_id = BaseConstants::_deactivated_bus_id;

// TODO all functions bellow are generic ! Make a base class for that
void GenericContainer::_get_amps(RealVect & a, const RealVect & p, const RealVect & q, const RealVect & v) const {
    RealVect p2q2 = p.array() * p.array() + q.array() * q.array();
    p2q2 = p2q2.array().cwiseSqrt();

    // modification in case of disconnected powerlines
    // because i don't want to divide by 0. below
    RealVect v_tmp = v;
    for(auto & el: v_tmp){
        if(abs(el) < _tol_equal_float) el = 1.0;
    }
    a = p2q2.array() * _1_sqrt_3 / v_tmp.array();
}

void GenericContainer::_generic_reactivate(int global_bus_id, SubstationContainer & substation){
    _check_in_range(static_cast<std::vector<bool>::size_type>(global_bus_id),
                    substation.get_bus_status(),
                    "_generic_reactivate");
    substation.reconnect_bus(global_bus_id);
    // status[el_id] = true;  //TODO why it's needed to do that again
}

void GenericContainer::_generic_deactivate(int global_bus_id, SubstationContainer & substation){
    _check_in_range(static_cast<std::vector<bool>::size_type>(global_bus_id),
                    substation.get_bus_status(),
                    "_generic_deactivate");
    substation.disconnect_bus(global_bus_id);
    // status[el_id] = false;
}

void GenericContainer::_generic_reactivate(int el_id, std::vector<bool> & eltype_status){
    _check_in_range(static_cast<std::vector<bool>::size_type>(el_id),
                    eltype_status,
                    "_generic_reactivate");
    eltype_status[el_id] = true;  //TODO why it's needed to do that again
}

void GenericContainer::_generic_deactivate(int el_id, std::vector<bool> & eltype_status){
    _check_in_range(static_cast<std::vector<bool>::size_type>(el_id),
                    eltype_status,
                    "_generic_deactivate");
    eltype_status[el_id] = false;   //TODO why it's needed to do that again
}

void GenericContainer::_generic_change_bus(
    int el_id,
    GridModelBusId new_gridmodel_bus_id,
    Eigen::Ref<GlobalBusIdVect> el_bus_ids,
    SolverControl & solver_control,
    int nb_max_bus) const {
    // bus id here "me_id" and NOT "solver_id"

    // throw error: object id does not exist
    _check_in_range(static_cast<Eigen::Index>(el_id),
                    el_bus_ids,
                    "_change_bus");

    int bus_id_int = new_gridmodel_bus_id;
    // throw error: bus id does not exist
    if(bus_id_int >= nb_max_bus)
    {
        // TODO DEBUG MODE: only check in debug mode
        std::ostringstream exc_;
        exc_ << "GenericContainer::_change_bus: Cannot change an element to bus ";
        exc_ << bus_id_int;
        exc_ << " There are only ";
        exc_ << nb_max_bus;
        exc_ << " distinct buses on this grid.";
        throw std::out_of_range(exc_.str());
    }
    if(bus_id_int < 0)
    {
        // TODO DEBUG MODE: only check in debug mode
        std::ostringstream exc_;
        exc_ << "GenericContainer::_change_bus: new bus id should be >=0 and not ";
        exc_ << bus_id_int;
        throw std::out_of_range(exc_.str());
    }
    GlobalBusId & bus_me_id = el_bus_ids(el_id);
    bus_me_id = new_gridmodel_bus_id;
}

GridModelBusId GenericContainer::_get_bus(int el_id, const std::vector<bool> & status_, const GlobalBusIdVect & bus_id_) const
{
    _check_in_range(static_cast<std::vector<bool>::size_type>(el_id),
                    status_,
                    "_get_bus");
    int res;
    bool val = status_[el_id];  // also check if the el_id is out of bound
    if(!val) res = _deactivated_bus_id;
    else{
        res = bus_id_(el_id);
    }
    return res;
}

void GenericContainer::v_kv_from_vpu(const Eigen::Ref<const RealVect> & Va,
                                     const Eigen::Ref<const RealVect> & Vm,
                                     const std::vector<bool> & status,
                                     int nb_element,
                                     const GlobalBusIdVect & bus_me_id,
                                     const std::vector<SolverBusId> & id_grid_to_solver,
                                     const RealVect & bus_vn_kv,
                                     RealVect & v) const
{
    for(int el_id = 0; el_id < nb_element; ++el_id){
        // if the element is disconnected, i leave it like that
        if(!status[el_id]) {
            v(el_id) = v_disco_el_;
            continue;
        }
        int el_bus_me_id = bus_me_id(el_id);
        int bus_solver_id = id_grid_to_solver[el_bus_me_id];
        if(bus_solver_id == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GenericContainer::v_kv_from_vpu: The element of id ";
            exc_ << bus_solver_id;
            exc_ << " is connected to a disconnected bus";
            throw std::runtime_error(exc_.str());
        }
        real_type bus_vn_kv_me = bus_vn_kv(el_bus_me_id);
        v(el_id) = Vm(bus_solver_id) * bus_vn_kv_me;
    }
}

void GenericContainer::v_deg_from_va(const Eigen::Ref<const RealVect> & Va,
                                     const Eigen::Ref<const RealVect> & Vm,
                                     const std::vector<bool> & status,
                                     int nb_element,
                                     const GlobalBusIdVect & bus_me_id,
                                     const std::vector<SolverBusId> & id_grid_to_solver,
                                     const RealVect & bus_vn_kv,
                                     RealVect & theta) const
{
    for(int el_id = 0; el_id < nb_element; ++el_id){
        // if the element is disconnected, i leave it like that
        if(!status[el_id]) {
            theta(el_id) = theta_disco_el_;
            continue;
        }
        int el_bus_me_id = bus_me_id(el_id);
        int bus_solver_id = id_grid_to_solver[el_bus_me_id];
        if(bus_solver_id == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GenericContainer::v_deg_from_va: The element of id ";
            exc_ << bus_solver_id;
            exc_ << " is connected to a disconnected bus";
            throw std::runtime_error(exc_.str());
        }
        theta(el_id) = Va(bus_solver_id) * my_180_pi_;
    }
}
