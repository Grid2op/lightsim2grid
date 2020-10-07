// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "DataGen.h"
#include <iostream>

void DataGen::init(const Eigen::VectorXd & generators_p,
                   const Eigen::VectorXd & generators_v,
                   const Eigen::VectorXd & generators_min_q,
                   const Eigen::VectorXd & generators_max_q,
                   const Eigen::VectorXi & generators_bus_id)
{
    p_mw_ = generators_p;
    vm_pu_ = generators_v;
    bus_id_ = generators_bus_id;
    min_q_ = generators_min_q;
    max_q_ = generators_max_q;
    if(min_q_.size() != max_q_.size()) throw std::runtime_error("Impossible to initialize generator with not the same size for min_q and max_q");
    int nb_gen = min_q_.size();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if (min_q_(gen_id) > max_q_(gen_id)) throw std::runtime_error("Impossible to initialize generator min_q being above max_q");
    }
    status_ = std::vector<bool>(generators_p.size(), true);
}


DataGen::StateRes DataGen::get_state() const
{
     std::vector<double> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<double> vm_pu(vm_pu_.begin(), vm_pu_.end());
     std::vector<double> min_q(min_q_.begin(), min_q_.end());
     std::vector<double> max_q(max_q_.begin(), max_q_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     DataGen::StateRes res(p_mw, vm_pu, min_q, max_q, bus_id, status);
     return res;
}
void DataGen::set_state(DataGen::StateRes & my_state )
{
    reset_results();

    std::vector<double> & p_mw = std::get<0>(my_state);
    std::vector<double> & vm_pu = std::get<1>(my_state);
    std::vector<double> & min_q = std::get<2>(my_state);
    std::vector<double> & max_q = std::get<3>(my_state);
    std::vector<int> & bus_id = std::get<4>(my_state);
    std::vector<bool> & status = std::get<5>(my_state);
    // TODO check sizes

    // input data
    p_mw_ = Eigen::VectorXd::Map(&p_mw[0], p_mw.size());
    vm_pu_ = Eigen::VectorXd::Map(&vm_pu[0], vm_pu.size());
    min_q_ = Eigen::VectorXd::Map(&min_q[0], min_q.size());
    max_q_ = Eigen::VectorXd::Map(&max_q[0], max_q.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
}


void DataGen::fillSbus(Eigen::VectorXcd & Sbus, bool ac, const std::vector<int> & id_grid_to_solver){
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    double tmp;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        tmp = p_mw_(gen_id);
        Sbus.coeffRef(bus_id_solver) += tmp;
    }
}

void DataGen::fillpv(std::vector<int> & bus_pv,
                     std::vector<bool> & has_bus_been_added,
                     int slack_bus_id_solver,
                     const std::vector<int> & id_grid_to_solver)
{
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }
        if(bus_id_solver == slack_bus_id_solver) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver);
        has_bus_been_added[bus_id_solver] = true;  // don't add it a second time
    }
}

void DataGen::compute_results(const Eigen::Ref<Eigen::VectorXd> & Va,
                               const Eigen::Ref<Eigen::VectorXd> & Vm,
                               const Eigen::Ref<Eigen::VectorXcd> & V,
                               const std::vector<int> & id_grid_to_solver,
                               const Eigen::VectorXd & bus_vn_kv)
{
    int nb_gen = nb();
    v_kv_from_vpu(Va, Vm, status_, nb_gen, bus_id_, id_grid_to_solver, bus_vn_kv, res_v_);
    res_p_ = p_mw_;
    // res_q_ = q_mvar_;
}

void DataGen::reset_results(){
    res_p_ = Eigen::VectorXd();  // in MW
    res_q_ = Eigen::VectorXd();  // in MVar
    res_v_ = Eigen::VectorXd();  // in kV
}

void DataGen::get_vm_for_dc(Eigen::VectorXd & Vm){
    int nb_gen = nb();
    int bus_id_me;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;
        bus_id_me = bus_id_(gen_id);
        double tmp = vm_pu_(gen_id);
        if(tmp != 0.) Vm(bus_id_me) = tmp;
    }
}

void DataGen::change_p(int gen_id, double new_p, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the active value of a disconnected generator");
    p_mw_(gen_id) = new_p;
}

void DataGen::change_v(int gen_id, double new_v_pu, bool & need_reset)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status) throw std::runtime_error("Impossible to change the voltage setpoint of a disconnected generator");
    vm_pu_(gen_id) = new_v_pu;
}

void DataGen::set_vm(Eigen::VectorXcd & V, const std::vector<int> & id_grid_to_solver)
{
    int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            //TODO improve error message with the gen_id
            throw std::runtime_error("One generator is connected to a disconnected bus.");
        }

        // scale the input V such that abs(V) = Vm for this generator
        double tmp = std::abs(V(bus_id_solver));
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

int DataGen::get_slack_bus_id(int gen_id){
    bool status = status_.at(gen_id);  // also to ensure gen_id is consistent with number of gen
    if(!status) throw std::runtime_error("Generator for slack bus is deactivated");
    int res = bus_id_(gen_id);
    return res;
}

void DataGen::set_p_slack(int slack_bus_id, double p_slack){
    bool status = status_.at(slack_bus_id);  // also to ensure gen_id is consistent with number of gen
    if(!status) throw std::runtime_error("Generator for slack bus is deactivated");
    res_p_(slack_bus_id) = p_slack;
}

void DataGen::init_q_vector(int nb_bus)
{
    int nb_gen = nb();
    total_q_min_per_bus_ = Eigen::VectorXd::Constant(nb_bus, 0.);
    total_q_max_per_bus_ = Eigen::VectorXd::Constant(nb_bus, 0.);
    total_gen_per_bus_ = Eigen::VectorXi::Constant(nb_bus, 0);
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;
        int bus_id = bus_id_(gen_id);
        total_q_min_per_bus_(bus_id) += min_q_(gen_id);
        total_q_max_per_bus_(bus_id) += max_q_(gen_id);
        total_gen_per_bus_(bus_id) += 1;
    }
}

void DataGen::set_q(const std::vector<double> & q_by_bus)
{
    // for(int bus_id = 0; bus_id < q_by_bus.size(); ++bus_id) std::cout << "bus id " << bus_id << " sum q " << q_by_bus[bus_id] << std::endl;
    int nb_gen = nb();
    res_q_ = Eigen::VectorXd::Constant(nb_gen, 0.);
    double eps_q = 1e-8;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        double real_q = 0.;
        if(!status_[gen_id]) continue;  // set at 0 for disconnected generators
        int bus_id = bus_id_(gen_id);
        double q_to_absorb = q_by_bus[bus_id];
        double max_q_me = max_q_(gen_id);
        double min_q_me = min_q_(gen_id);
        double max_q_bus = total_q_max_per_bus_(bus_id);
        double min_q_bus = total_q_min_per_bus_(bus_id);
        int nb_gen_with_me = total_gen_per_bus_(bus_id);
        if(nb_gen_with_me == 1){
            real_q = q_to_absorb;
        }else{
            double ratio = (max_q_me - min_q_me + eps_q) / (max_q_bus - min_q_bus + nb_gen_with_me * eps_q) ;
            real_q = q_to_absorb * ratio ;
        }
        res_q_(gen_id) = real_q;
    }
}

