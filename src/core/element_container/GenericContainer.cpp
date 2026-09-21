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

void GenericContainer::_throw_on_disconnected_bus(const char * fun_name, int el_id, bool grid_side)
{
    std::ostringstream exc_;
    exc_ << fun_name << ": the element with id " << el_id
         << " is connected to a disconnected bus while being connected to the grid ("
         << (grid_side ? "GlobalBusId" : "SolverBusId") << " is the deactivated bus id).";
    throw std::runtime_error(exc_.str());
}

void GenericContainer::_get_amps(Eigen::Ref<RealVect> a,
                                 const Eigen::Ref<const RealVect> & p,
                                 const Eigen::Ref<const RealVect> & q,
                                 const Eigen::Ref<const RealVect> & v) {
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
                                           Eigen::Ref<RealVect> theta)
{
    for(int el_id = 0; el_id < nb_element; ++el_id){
        // if the element is disconnected, i leave it like that
        if(!status[el_id]) {
            v(el_id) = v_disco_el_;
            theta(el_id) = theta_disco_el_;
            continue;
        }
        const GlobalBusId el_bus_me_id = bus_me_id(el_id);
        const SolverBusId bus_solver_id = _solver_bus(el_id, el_bus_me_id, id_grid_to_solver,
                                                     "GenericContainer::v_kv_theta_from_vpu");
        v(el_id) = Vm(bus_solver_id.cast_int()) * bus_vn_kv(el_bus_me_id.cast_int());
        theta(el_id) = Va(bus_solver_id.cast_int()) * my_180_pi_;
    }
}

} // namespace ls2g
