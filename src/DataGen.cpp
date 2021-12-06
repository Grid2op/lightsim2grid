// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataGen.h"
#include <iostream>
#include <sstream>

void DataGen::init(const RealVect & generators_p,
                   const RealVect & generators_v,
                   const RealVect & generators_min_q,
                   const RealVect & generators_max_q,
                   const Eigen::VectorXi & generators_bus_id)
{
    p_mw_ = generators_p;
    vm_pu_ = generators_v;
    bus_id_ = generators_bus_id;
    min_q_ = generators_min_q;
    max_q_ = generators_max_q;
    if(min_q_.size() != max_q_.size())
    {
        std::ostringstream exc_;
        exc_ << "DataGen::init: Impossible to initialize generator with generators_min_q of size ";
        exc_ << min_q_.size();
        exc_ << " and generators_max_q of size ";
        exc_ << max_q_.size();
        exc_ << ". Both should match";
        throw std::runtime_error(exc_.str());
    }
    const int nb_gen = static_cast<int>(min_q_.size());
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if (min_q_(gen_id) > max_q_(gen_id))
        {
            std::ostringstream exc_;
            exc_ << "DataGen::init: Impossible to initialize generator min_q being above max_q for generator ";
            exc_ << gen_id;
            throw std::runtime_error(exc_.str());
        }
    }
    status_ = std::vector<bool>(generators_p.size(), true);
    gen_slackbus_ = std::vector<bool>(generators_p.size(), false);
    gen_slack_weight_ = std::vector<real_type>(generators_p.size(), 0.);
}


DataGen::StateRes DataGen::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> vm_pu(vm_pu_.begin(), vm_pu_.end());
     std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
     std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     std::vector<bool> slack_bus = gen_slackbus_;
     std::vector<real_type> slack_weight = gen_slack_weight_;
     DataGen::StateRes res(p_mw, vm_pu, min_q, max_q, bus_id, status, slack_bus, slack_weight);
     return res;
}
void DataGen::set_state(DataGen::StateRes & my_state )
{
    reset_results();

    std::vector<real_type> & p_mw = std::get<0>(my_state);
    std::vector<real_type> & vm_pu = std::get<1>(my_state);
    std::vector<real_type> & min_q = std::get<2>(my_state);
    std::vector<real_type> & max_q = std::get<3>(my_state);
    std::vector<int> & bus_id = std::get<4>(my_state);
    std::vector<bool> & status = std::get<5>(my_state);
    std::vector<bool> & slack_bus = std::get<6>(my_state);
    std::vector<real_type> & slack_weight = std::get<7>(my_state);
    // TODO check sizes

    // input data
    p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
    vm_pu_ = RealVect::Map(&vm_pu[0], vm_pu.size());
    min_q_ = RealVect::Map(&min_q[0], min_q.size());
    max_q_ = RealVect::Map(&max_q[0], max_q.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
    gen_slackbus_ = slack_bus;
    gen_slack_weight_ = slack_weight;
}

RealVect DataGen::get_slack_weights(Eigen::Index nb_bus_solver, const std::vector<int> & id_grid_to_solver){
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    RealVect res = RealVect::Zero(nb_bus_solver);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;
        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE: only check in debug mode
            std::ostringstream exc_;
            exc_ << "DataGen::get_slack_weights: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        if(gen_slackbus_[gen_id])  res.coeffRef(bus_id_solver) += gen_slack_weight_[gen_id];
    }
    bus_slack_weight_ = res;
    real_type sum_res = res.sum();
    res /= sum_res;
    return res;
}

void DataGen::fillSbus(CplxVect & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    real_type tmp;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "DataGen::fillSbus: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = p_mw_(gen_id);
        Sbus.coeffRef(bus_id_solver) += tmp;
    }
}

void DataGen::fillpv(std::vector<int> & bus_pv,
                     std::vector<bool> & has_bus_been_added,
                     Eigen::VectorXi & slack_bus_id_solver,
                     const std::vector<int> & id_grid_to_solver) const
{
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "DataGen::fillpv: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        // if(bus_id_solver == slack_bus_id_solver) continue;  // slack bus is not PV
        if(is_in_vect(bus_id_solver, slack_bus_id_solver)) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver);
        has_bus_been_added[bus_id_solver] = true;  // don't add it a second time
    }
}

void DataGen::compute_results(const Eigen::Ref<const RealVect> & Va,
                              const Eigen::Ref<const RealVect> & Vm,
                              const Eigen::Ref<const CplxVect> & V,
                              const std::vector<int> & id_grid_to_solver,
                              const RealVect & bus_vn_kv,
                              real_type sn_mva,
                              bool ac)
{
    const int nb_gen = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_gen, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    v_deg_from_va(Va, Vm, status_, nb_gen, bus_id_, id_grid_to_solver, bus_vn_kv, res_theta_);
    res_p_ = p_mw_;
}

void DataGen::reset_results(){
    res_p_ = RealVect();  // in MW
    res_q_ = RealVect();  // in MVar
    res_v_ = RealVect();  // in kV
    res_theta_ = RealVect();  // in deg
    bus_slack_weight_ = RealVect();
}

void DataGen::get_vm_for_dc(RealVect & Vm){
    const int nb_gen = nb();
    int bus_id_me;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;
        bus_id_me = bus_id_(gen_id);
        real_type tmp = vm_pu_(gen_id);
        if(tmp != 0.) Vm(bus_id_me) = tmp;
    }
}

void DataGen::change_p(int gen_id, real_type new_p, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "DataGen::change_p: Impossible to change the active value of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    p_mw_(gen_id) = new_p;
}

void DataGen::change_v(int gen_id, real_type new_v_pu, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "DataGen::change_p: Impossible to change the voltage setpoint of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    vm_pu_(gen_id) = new_v_pu;
}

void DataGen::set_vm(CplxVect & V, const std::vector<int> & id_grid_to_solver)
{
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "DataGen::set_vm: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }

        // scale the input V such that abs(V) = Vm for this generator
        real_type tmp = std::abs(V(bus_id_solver));
        if(tmp == 0.)
        {
            // if it was 0. i force it to 1. (otherwise the rest of the computation would make it O. still)
            V(bus_id_solver) = 1.0;
            tmp = 1.0;
        }
        tmp = 1.0 / tmp;
        tmp *= vm_pu_(gen_id);
        V(bus_id_solver) *= tmp;
    }
}

std::vector<int> DataGen::get_slack_bus_id() const{
    std::vector<int> res;
    const auto nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(gen_slackbus_[gen_id]){
            const auto my_bus = bus_id_(gen_id);
            // do not add twice the same "slack bus"
            if(!is_in_vect(my_bus, res)) res.push_back(my_bus);
        }
    }
    return res;
}

void DataGen::set_p_slack(const RealVect& node_mismatch,
                          const std::vector<int> & id_grid_to_solver)
{
    if(bus_slack_weight_.size() == 0){
        // TODO DEBUG MODE: perform this check only in debug mode
        throw std::runtime_error("Impossible to set the active value of generators for the slack bus");
    }
    const auto nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(!status_[gen_id]) continue;  // nothing to do if gen is disconnected
        if(!gen_slackbus_[gen_id]) continue;  // nothing to do if it's not a slack
        if(gen_slack_weight_[gen_id] == 0.) continue; // nothing to do if no weights are associated to it
        const auto bus_id_me = bus_id_(gen_id);
        const auto bus_id_solver = id_grid_to_solver[bus_id_me];
        // TODO DEBUG MODE: check bus_id_solver >= 0
        // TODO DEBUG MODE: check bus_slack_weight_[bus_id_solver] > 0
        const auto total_contrib_slack = bus_slack_weight_(bus_id_solver);
        const auto my_contrib_slack = gen_slack_weight_[gen_id];
        // now take "my part"
        res_p_(gen_id) += node_mismatch(bus_id_solver) * my_contrib_slack / total_contrib_slack;
    }
}

void DataGen::init_q_vector(int nb_bus_total)  // total number of bus on the grid
{
    const int nb_gen = nb();
    total_q_min_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
    total_q_max_per_bus_ = RealVect::Constant(nb_bus_total, 0.);
    total_gen_per_bus_ = Eigen::VectorXi::Constant(nb_bus_total, 0);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;
        int bus_id = bus_id_(gen_id);
        total_q_min_per_bus_(bus_id) += min_q_(gen_id);
        total_q_max_per_bus_(bus_id) += max_q_(gen_id);
        total_gen_per_bus_(bus_id) += 1;
    }
}

void DataGen::set_q(const RealVect & reactive_mismatch, const std::vector<int> & id_grid_to_solver, bool ac)
{
    const int nb_gen = nb();
    res_q_ = RealVect::Constant(nb_gen, 0.);
    if(!ac) return;  // do not consider Q values in dc mode
    real_type eps_q = 1e-8;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        real_type real_q = 0.;
        if(!status_[gen_id]) continue;  // set at 0 for disconnected generators
        int bus_id = bus_id_(gen_id);
        const auto bus_solver = id_grid_to_solver[bus_id];
        // TODO DEBUG MODE: check that the bus is correct!
        real_type q_to_absorb = reactive_mismatch[bus_solver];
        real_type max_q_me = max_q_(gen_id);
        real_type min_q_me = min_q_(gen_id);
        real_type max_q_bus = total_q_max_per_bus_(bus_id);
        real_type min_q_bus = total_q_min_per_bus_(bus_id);
        int nb_gen_with_me = total_gen_per_bus_(bus_id);
        if(nb_gen_with_me == 1){
            real_q = q_to_absorb;
        }else{
            real_type ratio = (max_q_me - min_q_me + eps_q) / (max_q_bus - min_q_bus + nb_gen_with_me * eps_q) ;
            real_q = q_to_absorb * ratio ;
        }
        res_q_(gen_id) = real_q;
    }
}

