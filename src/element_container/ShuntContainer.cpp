// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "ShuntContainer.h"

#include <iostream>

void ShuntContainer::init(const RealVect & shunt_p_mw,
                          const RealVect & shunt_q_mvar,
                          const Eigen::VectorXi & shunt_bus_id)
{
    int size = static_cast<int>(shunt_p_mw.size());
    GenericContainer::check_size(shunt_p_mw, size, "shunt_p_mw");
    GenericContainer::check_size(shunt_q_mvar, size, "shunt_q_mvar");
    GenericContainer::check_size(shunt_bus_id, size, "shunt_bus_id");

    p_mw_ = shunt_p_mw;
    q_mvar_ = shunt_q_mvar;
    bus_id_ = shunt_bus_id;
    status_ = std::vector<bool>(p_mw_.size(), true); // by default everything is connected
    reset_results();
}

ShuntContainer::StateRes ShuntContainer::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     ShuntContainer::StateRes res(names_, p_mw, q_mvar, bus_id, status);
     return res;
}

void ShuntContainer::set_state(ShuntContainer::StateRes & my_state )
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

void ShuntContainer::fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                              bool ac,
                              const std::vector<int> & id_grid_to_solver,
                              real_type sn_mva) const
{
    if(!ac) return; // no shunt in DC

    const Eigen::Index nb_shunt = static_cast<int>(q_mvar_.size());
    cplx_type tmp;
    int bus_id_me, bus_id_solver;
    for(Eigen::Index shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        // assign diagonal coefficient
        tmp = {p_mw_(shunt_id), -q_mvar_(shunt_id)};

        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "ShuntContainer::fillYbus: the shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        if(sn_mva != 1.) tmp /= sn_mva;
        res.push_back(Eigen::Triplet<cplx_type> (bus_id_solver, bus_id_solver, tmp));
    }
}

void ShuntContainer::fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                                std::vector<Eigen::Triplet<real_type> > & Bpp,
                                const std::vector<int> & id_grid_to_solver,
                                real_type sn_mva,
                                FDPFMethod xb_or_bx) const
{
    const Eigen::Index nb_shunt = static_cast<int>(q_mvar_.size());
    real_type tmp;
    int bus_id_me, bus_id_solver;
    for(Eigen::Index shunt_id=0; shunt_id < nb_shunt; ++shunt_id){
        // i don't do anything if the shunt is disconnected
        if(!status_[shunt_id]) continue;

        bus_id_me = bus_id_(shunt_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            std::ostringstream exc_;
            exc_ << "ShuntContainer::fillBp_Bpp: the shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to a disconnected bus while being connected";
            throw std::runtime_error(exc_.str());
        }
        // assign diagonal coefficient
        tmp = q_mvar_(shunt_id);
        if(sn_mva != 1.) tmp /= sn_mva;
        Bpp.push_back(Eigen::Triplet<real_type> (bus_id_solver, bus_id_solver, tmp));  // -(-tmp) [-tmp for the "correct" value, but then for Bpp i have -(-tmp)]
    }
}

void ShuntContainer::fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const  // in DC i need that
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

void ShuntContainer::fillYbus_spmat(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int> & id_grid_to_solver){
    throw std::runtime_error("ShuntContainer::fillYbus_spmat: should not be used anymore !");
}

void ShuntContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
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
    // res_p_ = RealVect::Constant(nb_shunt, my_zero_);
    // res_q_ = RealVect::Constant(nb_shunt, my_zero_);
    for(int shunt_id = 0; shunt_id < nb_shunt; ++shunt_id){
        if(!status_[shunt_id]) {
            res_p_(shunt_id) = my_zero_;
            res_q_(shunt_id) = my_zero_;
            continue;
        }
        int bus_id_me = bus_id_(shunt_id);
        int bus_solver_id = id_grid_to_solver[bus_id_me];
        if(bus_solver_id == _deactivated_bus_id){
            throw std::runtime_error("ShuntContainer::compute_results: A shunt is connected to a disconnected bus.");
        }
        cplx_type E = V(bus_solver_id);
        cplx_type y = -my_one_ * (p_mw_(shunt_id) + my_i * q_mvar_(shunt_id)) / sn_mva;
        cplx_type I = y * E;
        I = std::conj(I);
        cplx_type s = E * I;
        res_p_(shunt_id) = std::real(s) * sn_mva;
        if(ac) res_q_(shunt_id) = std::imag(s) * sn_mva;
        else res_q_(shunt_id) = my_zero_;
    }
}

void ShuntContainer::reset_results(){
    res_p_ = RealVect(nb());  // in MW
    res_q_ = RealVect(nb());  // in MVar
    res_v_ = RealVect(nb());  // in kV
    res_theta_ = RealVect(nb());  // in deg
}

void ShuntContainer::change_p(int shunt_id, real_type new_p, SolverControl & solver_control)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the active value of a disconnected shunt");
    if(p_mw_(shunt_id) != new_p){
        solver_control.tell_recompute_ybus();
        solver_control.tell_recompute_sbus();  // in dc mode sbus is modified
        p_mw_(shunt_id) = new_p;
    }
}

void ShuntContainer::change_q(int shunt_id, real_type new_q, SolverControl & solver_control)
{
    bool my_status = status_.at(shunt_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the reactive value of a disconnected shunt");
    if(q_mvar_(shunt_id) != new_q){
        solver_control.tell_recompute_ybus();
        q_mvar_(shunt_id) = new_q;
    }
}

void ShuntContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const {
    const int nb_shunt = nb();
    for(int shunt_id = 0; shunt_id < nb_shunt; ++shunt_id)
    {
        if(!status_[shunt_id]) continue;
        const auto my_bus = bus_id_(shunt_id);
        if(my_bus == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "ShuntContainer::reconnect_connected_buses: Shunt with id ";
            exc_ << shunt_id;
            exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_shunt(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[my_bus] = true;  // this bus is connected
    }
}

void ShuntContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
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
