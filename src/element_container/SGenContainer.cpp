// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SGenContainer.h"
#include <iostream>

void SGenContainer::init(const RealVect & sgen_p,
                         const RealVect & sgen_q,
                         const RealVect & sgen_pmin,
                         const RealVect & sgen_pmax,
                         const RealVect & sgen_qmin,
                         const RealVect & sgen_qmax,
                         const Eigen::VectorXi & sgen_bus_id)
{
    int size = static_cast<int>(sgen_p.size());
    GenericContainer::check_size(sgen_p, size, "sgen_p");
    GenericContainer::check_size(sgen_q, size, "sgen_q");
    GenericContainer::check_size(sgen_pmin, size, "sgen_pmin");
    GenericContainer::check_size(sgen_pmax, size, "sgen_pmax");
    GenericContainer::check_size(sgen_qmin, size, "sgen_qmin");
    GenericContainer::check_size(sgen_qmax, size, "sgen_qmax");
    GenericContainer::check_size(sgen_bus_id, size, "sgen_bus_id");

    p_mw_ = sgen_p;
    q_mvar_ = sgen_q;
    p_min_mw_ = sgen_pmin;
    p_max_mw_ = sgen_pmax;
    q_min_mvar_ = sgen_qmin;
    q_max_mvar_ = sgen_qmax;
    bus_id_ = sgen_bus_id;
    status_ = std::vector<bool>(sgen_p.size(), true);
    reset_results();
}

SGenContainer::StateRes SGenContainer::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<real_type> p_min(p_min_mw_.begin(), p_min_mw_.end());
     std::vector<real_type> p_max(p_max_mw_.begin(), p_max_mw_.end());
     std::vector<real_type> q_min(q_min_mvar_.begin(), q_min_mvar_.end());
     std::vector<real_type> q_max(q_max_mvar_.begin(), q_max_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     SGenContainer::StateRes res(names_, p_mw, q_mvar, p_min, p_max, q_min, q_max, bus_id, status);
     return res;
}

void SGenContainer::set_state(SGenContainer::StateRes & my_state )
{    
    names_ = std::get<0>(my_state);
    std::vector<real_type> & p_mw = std::get<1>(my_state);
    std::vector<real_type> & q_mvar = std::get<2>(my_state);
    std::vector<real_type> & p_min = std::get<3>(my_state);
    std::vector<real_type> & p_max = std::get<4>(my_state);
    std::vector<real_type> & q_min = std::get<5>(my_state);
    std::vector<real_type> & q_max = std::get<6>(my_state);
    std::vector<int> & bus_id = std::get<7>(my_state);
    std::vector<bool> & status = std::get<8>(my_state);
    auto size = p_mw.size();
    GenericContainer::check_size(p_mw, size, "p_mw");
    GenericContainer::check_size(q_mvar, size, "q_mvar");
    GenericContainer::check_size(p_min, size, "p_min");
    GenericContainer::check_size(p_max, size, "p_max");
    GenericContainer::check_size(q_min, size, "q_min");
    GenericContainer::check_size(q_max, size, "q_max");
    GenericContainer::check_size(bus_id, size, "bus_id");
    GenericContainer::check_size(status, size, "status");

    p_mw_ = RealVect::Map(&p_mw[0], size);
    q_mvar_ = RealVect::Map(&q_mvar[0], size);
    q_mvar_ = RealVect::Map(&q_mvar[0], size);
    p_min_mw_ = RealVect::Map(&p_min[0], size);
    p_max_mw_ = RealVect::Map(&p_max[0], size);
    q_min_mvar_ = RealVect::Map(&q_min[0], size);
    q_max_mvar_ = RealVect::Map(&q_max[0], size);
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
    reset_results();
}

void SGenContainer::fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const {
    const int nb_sgen = nb();
    int bus_id_me, bus_id_solver;
    cplx_type tmp;
    for(int sgen_id = 0; sgen_id < nb_sgen; ++sgen_id){
        //  i don't do anything if the static generator is disconnected
        if(!status_[sgen_id]) continue;

        bus_id_me = bus_id_(sgen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "SGenContainer::fillSbus: Static Generator with id ";
            exc_ << sgen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = {p_mw_(sgen_id), q_mvar_(sgen_id)};
        Sbus.coeffRef(bus_id_solver) += tmp;
    }
}

void SGenContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
                                    const Eigen::Ref<const RealVect> & Vm,
                                    const Eigen::Ref<const CplxVect> & V,
                                    const std::vector<int> & id_grid_to_solver,
                                    const RealVect & bus_vn_kv,
                                    real_type sn_mva,
                                    bool ac)
{
    const int nb_sgen = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_sgen, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_sgen, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    res_p_ = p_mw_;
    if(ac) res_q_ = q_mvar_;
    else{
        // no q in DC mode
        for(int sgen_id = 0; sgen_id < nb_sgen; ++sgen_id) res_q_(sgen_id) = 0.;
    }
}

void SGenContainer::reset_results(){
    res_p_ = RealVect(nb());  // in MW
    res_q_ =  RealVect(nb());  // in MVar
    res_v_ = RealVect(nb());  // in kV
    res_theta_ = RealVect(nb());  // in deg
}

void SGenContainer::change_p(int sgen_id, real_type new_p, SolverControl & solver_control)
{
    bool my_status = status_.at(sgen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "SGenContainer::change_p: Impossible to change the active value of a disconnected static generator (check sgen id ";
        exc_ << sgen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    if (p_mw_(sgen_id) != new_p){
        solver_control.tell_recompute_sbus();
        p_mw_(sgen_id) = new_p;
    }
}

void SGenContainer::change_q(int sgen_id, real_type new_q, SolverControl & solver_control)
{
    bool my_status = status_.at(sgen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "SGenContainer::change_q: Impossible to change the reactive value of a disconnected static generator (check sgen id ";
        exc_ << sgen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    if (q_mvar_(sgen_id) != new_q){
        solver_control.tell_recompute_sbus();
        q_mvar_(sgen_id) = new_q;
    }
}

void SGenContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const {
    const int nb_sgen = nb();
    for(int sgen_id = 0; sgen_id < nb_sgen; ++sgen_id)
    {
        if(!status_[sgen_id]) continue;
        const auto my_bus = bus_id_(sgen_id);
        if(my_bus == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "SGenContainer::reconnect_connected_buses: Static Generator with id ";
            exc_ << sgen_id;
            exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_sgen(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[my_bus] = true;  // this bus is connected
    }
}

void SGenContainer::gen_p_per_bus(std::vector<real_type> & res) const
{
    const int nb_gen = nb();
    for(int sgen_id = 0; sgen_id < nb_gen; ++sgen_id)
    {
        if(!status_[sgen_id]) continue;
        const auto my_bus = bus_id_(sgen_id);
        res[my_bus] += p_mw_(sgen_id);
    }
}

void SGenContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
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
