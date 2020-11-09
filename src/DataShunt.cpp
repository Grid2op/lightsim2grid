// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataShunt.h"
#include <iostream>

void DataShunt::init(const Eigen::VectorXd & shunt_p_mw,
                     const Eigen::VectorXd & shunt_q_mvar,
                     const Eigen::VectorXi & shunt_bus_id)
{
    p_mw_ = shunt_p_mw;
    q_mvar_ = shunt_q_mvar;
    bus_id_ = shunt_bus_id;
    status_ = std::vector<bool>(p_mw_.size(), true); // by default everything is connected
}


DataShunt::StateRes DataShunt::get_state() const
{
     std::vector<double> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<double> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     DataShunt::StateRes res(p_mw, q_mvar, bus_id, status);
     return res;
}
void DataShunt::set_state(DataShunt::StateRes & my_state )
{
    reset_results();

    std::vector<double> & p_mw = std::get<0>(my_state);
    std::vector<double> & q_mvar = std::get<1>(my_state);
    std::vector<int> & bus_id = std::get<2>(my_state);
    std::vector<bool> & status = std::get<3>(my_state);
    // TODO check sizes

    // input data
    p_mw_ = Eigen::VectorXd::Map(&p_mw[0], p_mw.size());
    q_mvar_ = Eigen::VectorXd::Map(&q_mvar[0], q_mvar.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
}

void DataShunt::fillYbus(std::vector<Eigen::Triplet<cdouble> > & res, bool ac, const std::vector<int> & id_grid_to_solver){
    int nb_shunt = q_mvar_.size();
    cdouble tmp;
    int bus_id_me, bus_id_solver;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        // assign diagonal coefficient
        tmp = p_mw_(shunt_id) + my_i * q_mvar_(shunt_id);
        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            throw std::runtime_error("GridModel::fillYbusShunt: A shunt is connected to a disconnected bus.");
        }
        res.push_back(Eigen::Triplet<cdouble> (bus_id_solver, bus_id_solver, -tmp));
    }
}
void DataShunt::fillYbus_spmat(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int> & id_grid_to_solver){
    int nb_shunt = q_mvar_.size();
    cdouble tmp;
    int bus_id_me, bus_id_solver;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        // assign diagonal coefficient
        tmp = p_mw_(shunt_id) + my_i * q_mvar_(shunt_id);
        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            throw std::runtime_error("GridModel::fillYbusShunt: A shunt is connected to a disconnected bus.");
        }
        res.coeffRef(bus_id_solver, bus_id_solver) -= tmp;
    }
}

void DataShunt::compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                                const Eigen::Ref<Eigen::VectorXd> & Vm,
                                const Eigen::Ref<Eigen::VectorXcd> & V,
                                const std::vector<int> & id_grid_to_solver,
                                const Eigen::VectorXd & bus_vn_kv)
{
    int nb_shunt = p_mw_.size();
    v_kv_from_vpu(Va, Vm, status_, nb_shunt, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    res_p_ = Eigen::VectorXd::Constant(nb_shunt, 0.);
    res_q_ = Eigen::VectorXd::Constant(nb_shunt, 0.);
    for(int shunt_id = 0; shunt_id < nb_shunt; ++shunt_id){
        if(!status_[shunt_id]) continue;
        int bus_id_me = bus_id_(shunt_id);
        int bus_solver_id = id_grid_to_solver[bus_id_me];
        if(bus_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataShunt::compute_results: A shunt is connected to a disconnected bus.");
        }
        cdouble E = V(bus_solver_id);
        cdouble y = -1.0 * (p_mw_(shunt_id) + my_i * q_mvar_(shunt_id));
        cdouble I = y * E;
        I = std::conj(I);
        cdouble s = E * I;
        res_p_(shunt_id) = std::real(s);
        res_q_(shunt_id) = std::imag(s);
    }
}

void DataShunt::reset_results(){
    res_p_ = Eigen::VectorXd();  // in MW
    res_q_ = Eigen::VectorXd();  // in MVar
    res_v_ = Eigen::VectorXd();  // in kV
}

void DataShunt::change_p(int shunt_id, double new_p, bool & need_reset)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the active value of a disconnected shunt");
    if(p_mw_(shunt_id) != new_p) need_reset = true;
    p_mw_(shunt_id) = new_p;

}

void DataShunt::change_q(int shunt_id, double new_q, bool & need_reset)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the reactive value of a disconnected shunt");
    if(q_mvar_(shunt_id) != new_q) need_reset = true;
    q_mvar_(shunt_id) = new_q;
}

double DataShunt::get_p_slack(int slack_bus_id)
{
    int nb_element = nb();
    double res = 0.;
    for(int shunt_id = 0; shunt_id < nb_element; ++shunt_id)
    {
        if(!status_[shunt_id]) continue;
        if(bus_id_(shunt_id) == slack_bus_id) res += res_p_(shunt_id);  //TODO plus here or minus ???
    }
    return res;
}

void DataShunt::get_q(std::vector<double>& q_by_bus)
{
    int nb_element = nb();
    for(int shunt_id = 0; shunt_id < nb_element; ++shunt_id)
    {
        if(!status_[shunt_id]) continue;
        int bus_id = bus_id_[shunt_id];
        q_by_bus[bus_id] += res_q_(shunt_id);  //TODO plus here or minus ???
    }
}


