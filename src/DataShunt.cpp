// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataShunt.h"
#include <iostream>

void DataShunt::init(const RealVect & shunt_p_mw,
                     const RealVect & shunt_q_mvar,
                     const Eigen::VectorXi & shunt_bus_id)
{
    p_mw_ = shunt_p_mw;
    q_mvar_ = shunt_q_mvar;
    bus_id_ = shunt_bus_id;
    status_ = std::vector<bool>(p_mw_.size(), true); // by default everything is connected
}


DataShunt::StateRes DataShunt::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     DataShunt::StateRes res(p_mw, q_mvar, bus_id, status);
     return res;
}
void DataShunt::set_state(DataShunt::StateRes & my_state )
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

void DataShunt::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                         bool ac,
                         const std::vector<int> & id_grid_to_solver,
                         real_type sn_mva){
    const int nb_shunt = static_cast<int>(q_mvar_.size());
    cplx_type tmp;
    int bus_id_me, bus_id_solver;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        // assign diagonal coefficient
        tmp = {p_mw_(shunt_id), my_zero_};
        if(ac) tmp += my_i * q_mvar_(shunt_id);

        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "DataShunt::fillYbus: the shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        if(sn_mva != 1.) tmp /= sn_mva;
        res.push_back(Eigen::Triplet<cplx_type> (bus_id_solver, bus_id_solver, -tmp));
    }
}

void DataShunt::fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver)  // in DC i need that
{
    if(ac) return;  // in AC I do not do that
    // std::cout << " ok i use this function" << std::endl;
    // - bus[:, GS] / baseMVA  # in pandapower
    // yish=gish+jbish -> so g is the MW !
    const int nb_shunt = static_cast<int>(q_mvar_.size());
    int bus_id_me, bus_id_solver;
    for(int shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;
        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            throw std::runtime_error("GridModel::fillSbus: A shunt is connected to a disconnected bus.");
        }
        Sbus.coeffRef(bus_id_solver) -= p_mw_(shunt_id);
    }
}

void DataShunt::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver){
    const int nb_shunt = static_cast<int>(q_mvar_.size());
    //TODO this is no more used!!!! see the other fillYbus
    cplx_type tmp;
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
    //TODO this is no more used!!!! see the other fillYbus
}

void DataShunt::compute_results(const Eigen::Ref<const RealVect> & Va,
                                const Eigen::Ref<const RealVect> & Vm,
                                const Eigen::Ref<const CplxVect> & V,
                                const std::vector<int> & id_grid_to_solver,
                                const RealVect & bus_vn_kv,
                                real_type sn_mva,
                                bool ac)
{
    const int nb_shunt = static_cast<int>(p_mw_.size());
    v_kv_from_vpu(Va, Vm, status_, nb_shunt, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_shunt, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    res_p_ = RealVect::Constant(nb_shunt, my_zero_);
    res_q_ = RealVect::Constant(nb_shunt, my_zero_);
    for(int shunt_id = 0; shunt_id < nb_shunt; ++shunt_id){
        if(!status_[shunt_id]) continue;
        int bus_id_me = bus_id_(shunt_id);
        int bus_solver_id = id_grid_to_solver[bus_id_me];
        if(bus_solver_id == _deactivated_bus_id){
            throw std::runtime_error("DataShunt::compute_results: A shunt is connected to a disconnected bus.");
        }
        cplx_type E = V(bus_solver_id);
        cplx_type y = -my_one_ * (p_mw_(shunt_id) + my_i * q_mvar_(shunt_id)) / sn_mva;
        cplx_type I = y * E;
        I = std::conj(I);
        cplx_type s = E * I;
        res_p_(shunt_id) = std::real(s) * sn_mva;
        if(ac) res_q_(shunt_id) = std::imag(s) * sn_mva;
    }
}

void DataShunt::reset_results(){
    res_p_ = RealVect();  // in MW
    res_q_ = RealVect();  // in MVar
    res_v_ = RealVect();  // in kV
}

void DataShunt::change_p(int shunt_id, real_type new_p, bool & need_reset)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the active value of a disconnected shunt");
    if(p_mw_(shunt_id) != new_p) need_reset = true;
    p_mw_(shunt_id) = new_p;

}

void DataShunt::change_q(int shunt_id, real_type new_q, bool & need_reset)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the reactive value of a disconnected shunt");
    if(q_mvar_(shunt_id) != new_q) need_reset = true;
    q_mvar_(shunt_id) = new_q;
}
