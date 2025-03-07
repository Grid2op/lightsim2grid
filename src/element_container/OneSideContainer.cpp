// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "OneSideContainer.h"
#include <sstream>
#include <iostream>

// void OneSideContainer::init_base(const RealVect & els_p,
//                                  const RealVect & els_q,
//                                  const Eigen::VectorXi & els_bus_id,
//                                  const std::string & name_el
//                                  )
// {
//     int size = static_cast<int>(els_p.size());
//     check_size(els_p, size, name_el + "_p");
//     check_size(els_q, size, name_el + "_q");
//     check_size(els_bus_id, size, name_el + "_bus_id");

//     p_mw_ = els_p;
//     q_mvar_ = els_q;
//     bus_id_ = els_bus_id;
//     status_ = std::vector<bool>(els_p.size(), true);
// }

// OneSideContainer::StateRes OneSideContainer::get_base_state() const
// {
//      std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
//      std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
//      std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
//      std::vector<bool> status = status_;
//      OneSideContainer::StateRes res(names_, p_mw, q_mvar, bus_id, status);
//      return res;
// }

// void OneSideContainer::set_base_state(OneSideContainer::StateRes & my_state)
// {
//     // read data
//     names_ = std::get<0>(my_state);
//     std::vector<real_type> & p_mw = std::get<1>(my_state);
//     std::vector<real_type> & q_mvar = std::get<2>(my_state);
//     std::vector<int> & bus_id = std::get<3>(my_state);
//     std::vector<bool> & status = std::get<4>(my_state);

//     // check sizes
//     const auto size = p_mw.size();
//     check_size(names_, size, "names");
//     check_size(p_mw, size, "p_mw");
//     check_size(q_mvar, size, "q_mvar");
//     check_size(bus_id, size, "bus_id");
//     check_size(status, size, "status");

//     // input data
//     p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
//     q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
//     bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
//     status_ = status;
// }

void OneSideContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                       const Eigen::Ref<const RealVect> & Vm,
                                       const Eigen::Ref<const CplxVect> & V,
                                       const std::vector<int> & id_grid_to_solver,
                                       const RealVect & bus_vn_kv,
                                       real_type sn_mva,
                                       bool ac)
{
    const int nb_els = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_els, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    this->_compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
}

void OneSideContainer::reset_osc_results(){
    // std::cout << "Loads reset_results \n";
    res_p_ = RealVect(nb());  // in MW
    res_q_ =  RealVect(nb());  // in MVar
    res_v_ = RealVect(nb());  // in kV
    res_theta_ = RealVect(nb());  // in deg
    this->_reset_results();
}

// void OneSideContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const {
//     const int nb_els = nb();
//     for(int el_id = 0; el_id < nb_els; ++el_id)
//     {
//         if(!status_[el_id]) continue;
//         const auto my_bus = bus_id_(el_id);
//         if(my_bus == _deactivated_bus_id){
//             // TODO DEBUG MODE only this in debug mode
//             std::ostringstream exc_;
//             exc_ << "OneSideContainer::reconnect_connected_buses: element with id ";
//             exc_ << el_id;
//             exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_xxx(...)` ?.";
//             throw std::runtime_error(exc_.str());
//         }
//         bus_status[my_bus] = true;  // this bus is connected
//     }
// }

// void OneSideContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
//     const int nb_el = nb();
//     SolverControl unused_solver_control;
//     for(int el_id = 0; el_id < nb_el; ++el_id)
//     {
//         if(!status_[el_id]) continue;
//         const auto my_bus = bus_id_(el_id);
//         if(!busbar_in_main_component[my_bus]){
//             deactivate(el_id, unused_solver_control);
//         }
//     }    
// }
