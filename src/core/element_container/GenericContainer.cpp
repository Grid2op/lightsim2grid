// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GenericContainer.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

namespace ls2g {

const int GenericContainer::_deactivated_bus_id = BaseConstants::_deactivated_bus_id;

// TODO all functions bellow are generic ! Make a base class for that
void GenericContainer::_get_amps(Eigen::Ref<RealVect> a,
                                  const Eigen::Ref<const RealVect> & p,
                                  const Eigen::Ref<const RealVect> & q,
                                  const Eigen::Ref<const RealVect> & v) const {
    // One pass, no temporaries. This was four: sum of squares into a vector,
    // square root over it, a copy of v, then a scan of that copy to replace the
    // zeros -- two full-length allocations per call, and this is called four
    // times per solve (both ends of the powerlines and of the transformers).
    // Same arithmetic in the same order, so the values are bit-identical; the
    // guard on v is what stops a disconnected element (v = 0) dividing by zero.
    const Eigen::Index nb_el = a.size();
    for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
        const real_type v_el = v(el_id);
        const real_type v_div = (std::abs(v_el) < _tol_equal_float) ? 1.0 : v_el;
        const real_type p_el = p(el_id);
        const real_type q_el = q(el_id);
        a(el_id) = std::sqrt(p_el * p_el + q_el * q_el) * _1_sqrt_3 / v_div;
    }
}

void GenericContainer::_generic_reactivate(const GlobalBusId & global_bus_id, SubstationContainer & substation){
    _check_in_range(static_cast<std::vector<bool>::size_type>(global_bus_id.cast_int()),
                    substation.get_bus_status(),
                    "_generic_reactivate");
    substation.reconnect_bus(global_bus_id);
    // status[el_id] = true;  //TODO why it's needed to do that again
}

void GenericContainer::_generic_deactivate(const GlobalBusId & global_bus_id, SubstationContainer & substation){
    _check_in_range(static_cast<std::vector<bool>::size_type>(global_bus_id.cast_int()),
                    substation.get_bus_status(),
                    "_generic_deactivate");
    substation.disconnect_bus(global_bus_id);
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
    const GridModelBusId & new_gridmodel_bus_id,
    GlobalBusIdVect & el_bus_ids,
    DualAlgoControl & /*solver_control*/,
    int nb_max_bus) const {
    // bus id here "me_id" and NOT "solver_id"

    // throw error: object id does not exist
    _check_in_range(static_cast<Eigen::Index>(el_id),
                    el_bus_ids,
                    "_change_bus");

    // throw error: bus id does not exist
    if(new_gridmodel_bus_id.cast_int() >= nb_max_bus)
    {
        // TODO DEBUG MODE: only check in debug mode
        std::ostringstream exc_;
        exc_ << "GenericContainer::_change_bus: Cannot change an element to bus ";
        exc_ << new_gridmodel_bus_id.cast_int();
        exc_ << " There are only ";
        exc_ << nb_max_bus;
        exc_ << " distinct buses on this grid.";
        throw std::out_of_range(exc_.str());
    }
    if(new_gridmodel_bus_id.cast_int() < 0)
    {
        // TODO DEBUG MODE: only check in debug mode
        std::ostringstream exc_;
        exc_ << "GenericContainer::_change_bus: new bus id should be >=0 and not ";
        exc_ << new_gridmodel_bus_id.cast_int();
        throw std::out_of_range(exc_.str());
    }
    auto bus_me_id = el_bus_ids(el_id);
    bus_me_id = new_gridmodel_bus_id;
}

GridModelBusId GenericContainer::_get_bus(int el_id, const std::vector<bool> & status_, const GlobalBusIdVect & bus_id_) const
{
    _check_in_range(static_cast<std::vector<bool>::size_type>(el_id),
                    status_,
                    "_get_bus");
    GridModelBusId res;
    bool val = status_[el_id];  // also check if the el_id is out of bound
    if(!val) res = GridModelBusId(_deactivated_bus_id);
    else{
        res = bus_id_(el_id);
    }
    return res;
}

void GenericContainer::v_kv_theta_from_vpu(const Eigen::Ref<const RealVect> & Va,
                                           const Eigen::Ref<const RealVect> & Vm,
                                           const std::vector<bool> & status,
                                           int nb_element,
                                           const GlobalBusIdVect & bus_me_id,
                                           const SolverBusIdVect & id_grid_to_solver,
                                           const Eigen::Ref<const RealVect> & bus_vn_kv,
                                           Eigen::Ref<RealVect> v,
                                           Eigen::Ref<RealVect> theta) const
{
    for(int el_id = 0; el_id < nb_element; ++el_id){
        // if the element is disconnected, i leave it like that
        if(!status[el_id]) {
            v(el_id) = v_disco_el_;
            theta(el_id) = theta_disco_el_;
            continue;
        }
        const GlobalBusId el_bus_me_id = bus_me_id(el_id);
        if(el_bus_me_id.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GenericContainer::v_kv_theta_from_vpu: element with id ";
            exc_ << el_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        const SolverBusId bus_solver_id = id_grid_to_solver[el_bus_me_id.cast_int()];
        if(bus_solver_id.cast_int() == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "GenericContainer::v_kv_theta_from_vpu: The element of id ";
            exc_ << el_id;
            exc_ << " is connected to a disconnected bus";
            throw std::runtime_error(exc_.str());
        }
        v(el_id) = Vm(bus_solver_id.cast_int()) * bus_vn_kv(el_bus_me_id.cast_int());
        theta(el_id) = Va(bus_solver_id.cast_int()) * my_180_pi_;
    }
}

} // namespace ls2g
