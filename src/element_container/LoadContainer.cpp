// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "LoadContainer.h"
#include <sstream>
#include <iostream>

void LoadContainer::init(const RealVect & loads_p,
                         const RealVect & loads_q,
                         const Eigen::VectorXi & loads_bus_id)
{
    int size = static_cast<int>(loads_p.size());
    GenericContainer::check_size(loads_p, size, "loads_p");
    GenericContainer::check_size(loads_q, size, "loads_q");
    GenericContainer::check_size(loads_bus_id, size, "loads_bus_id");

    p_mw_ = loads_p;
    q_mvar_ = loads_q;
    bus_id_ = loads_bus_id;
    status_ = std::vector<bool>(loads_p.size(), true);
    reset_results();
}

LoadContainer::StateRes LoadContainer::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     LoadContainer::StateRes res(names_, p_mw, q_mvar, bus_id, status);
     return res;
}

void LoadContainer::set_state(LoadContainer::StateRes & my_state )
{
    names_ = std::get<0>(my_state);
    std::vector<real_type> & p_mw = std::get<1>(my_state);
    std::vector<real_type> & q_mvar = std::get<2>(my_state);
    std::vector<int> & bus_id = std::get<3>(my_state);
    std::vector<bool> & status = std::get<4>(my_state);
    // TODO check sizes

    // input data
    p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
    q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
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
        tmp = static_cast<cplx_type>(p_mw_(load_id));
        tmp += my_i * q_mvar_(load_id);
        Sbus.coeffRef(bus_id_solver) -= tmp;
    }
}

void LoadContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                    const Eigen::Ref<const RealVect> & Vm,
                                    const Eigen::Ref<const CplxVect> & V,
                                    const std::vector<int> & id_grid_to_solver,
                                    const RealVect & bus_vn_kv,
                                    real_type sn_mva,
                                    bool ac)
{
    const int nb_load = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_load, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_load, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    res_p_ = p_mw_;
    if(ac) res_q_ = q_mvar_;
    else{
        // no q in DC mode
        for(int load_id = 0; load_id < nb_load; ++load_id) res_q_(load_id) = 0.;
    }
}

void LoadContainer::reset_results(){
    // std::cout << "Loads reset_results \n";
    res_p_ = RealVect(nb());  // in MW
    res_q_ =  RealVect(nb());  // in MVar
    res_v_ = RealVect(nb());  // in kV
    res_theta_ = RealVect(nb());  // in deg
}

void LoadContainer::change_p(int load_id, real_type new_p, SolverControl & solver_control)
{
    bool my_status = status_.at(load_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "LoadContainer::change_p: Impossible to change the active value of a disconnected load (check load id ";
        exc_ << load_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    if (p_mw_(load_id) != new_p) {
        solver_control.tell_recompute_sbus();
        p_mw_(load_id) = new_p;
    }
}

void LoadContainer::change_q(int load_id, real_type new_q, SolverControl & solver_control)
{
    bool my_status = status_.at(load_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "LoadContainer::change_q: Impossible to change the reactive value of a disconnected load (check load id ";
        exc_ << load_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    if (q_mvar_(load_id) != new_q) {
        solver_control.tell_recompute_sbus();
        q_mvar_(load_id) = new_q;
    }
}

void LoadContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const {
    const int nb_load = nb();
    for(int load_id = 0; load_id < nb_load; ++load_id)
    {
        if(!status_[load_id]) continue;
        const auto my_bus = bus_id_(load_id);
        if(my_bus == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "LoadContainer::reconnect_connected_buses: Load with id ";
            exc_ << load_id;
            exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_load(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[my_bus] = true;  // this bus is connected
    }
}

void LoadContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
    const int nb_el = nb();
    SolverControl unused_solver_control;
    for(int el_id = 0; el_id < nb_el; ++el_id)
    {
        if(!status_[el_id]) continue;
        const auto my_bus = bus_id_(el_id);
        if(!busbar_in_main_component[my_bus]){
            deactivate(el_id, unused_solver_control);
        }
    }    
}
