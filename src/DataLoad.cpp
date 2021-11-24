// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataLoad.h"
#include <sstream>

void DataLoad::init(const RealVect & loads_p,
                    const RealVect & loads_q,
                    const Eigen::VectorXi & loads_bus_id)
{
    p_mw_ = loads_p;
    q_mvar_ = loads_q;
    bus_id_ = loads_bus_id;
    status_ = std::vector<bool>(loads_p.size(), true);
}


DataLoad::StateRes DataLoad::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     DataLoad::StateRes res(p_mw, q_mvar, bus_id, status);
     return res;
}
void DataLoad::set_state(DataLoad::StateRes & my_state )
{
    reset_results();

    std::vector<real_type> & p_mw = std::get<0>(my_state);
    std::vector<real_type> & q_mvar = std::get<1>(my_state);
    std::vector<int> & bus_id = std::get<2>(my_state);
    std::vector<bool> & status = std::get<3>(my_state);
    // TODO check sizes

    // input data
    p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
    q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
}


void DataLoad::fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){
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
            exc_ << "DataLoad::fillSbus: the load with id ";
            exc_ << load_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        tmp = static_cast<cplx_type>(p_mw_(load_id));
        if(ac) tmp += my_i * q_mvar_(load_id);
        Sbus.coeffRef(bus_id_solver) -= tmp;
    }
}

void DataLoad::compute_results(const Eigen::Ref<const RealVect> & Va,
                               const Eigen::Ref<const RealVect> & Vm,
                               const Eigen::Ref<const CplxVect> & V,
                               const std::vector<int> & id_grid_to_solver,
                               const RealVect & bus_vn_kv,
                               real_type sn_mva,
                               bool ac)
{
    int nb_load = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_load, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_load, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    res_p_ = p_mw_;
    res_q_ = q_mvar_;
}

void DataLoad::reset_results(){
    res_p_ = RealVect();  // in MW
    res_q_ =  RealVect();  // in MVar
    res_v_ = RealVect();  // in kV
}

void DataLoad::change_p(int load_id, real_type new_p, bool & need_reset)
{
    bool my_status = status_.at(load_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "DataLoad::change_p: Impossible to change the active value of a disconnected load (check load id ";
        exc_ << load_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    p_mw_(load_id) = new_p;
}

void DataLoad::change_q(int load_id, real_type new_q, bool & need_reset)
{
    bool my_status = status_.at(load_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        std::ostringstream exc_;
        exc_ << "DataLoad::change_q: Impossible to change the reactive value of a disconnected load (check load id ";
        exc_ << load_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    q_mvar_(load_id) = new_q;
}
