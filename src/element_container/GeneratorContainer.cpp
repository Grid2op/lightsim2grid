// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GeneratorContainer.h"
#include <iostream>
#include <sstream>

void GeneratorContainer::init(const RealVect & generators_p,
                              const RealVect & generators_v,
                              const RealVect & generators_min_q,
                              const RealVect & generators_max_q,
                              const Eigen::VectorXi & generators_bus_id)
{
    const auto generators_q = RealVect::Zero(generators_p.size());
    const auto voltage_regulator_on = std::vector<bool>(generators_p.size(), true);
    init_full(generators_p, generators_v, generators_q,
              voltage_regulator_on,
              generators_min_q, generators_max_q,
              generators_bus_id);
}

void GeneratorContainer::init_full(const RealVect & generators_p,
                                   const RealVect & generators_v,
                                   const RealVect & generators_q,
                                   const std::vector<bool> & voltage_regulator_on,
                                   const RealVect & generators_min_q,
                                   const RealVect & generators_max_q,
                                   const Eigen::VectorXi & generators_bus_id
                                   )
{
    init_osc(generators_p, generators_q, generators_bus_id, "generators");

    // check the sizes
    int size = static_cast<int>(generators_p.size());
    check_size(generators_v, size, "generators_v");
    check_size(generators_min_q, size, "generators_min_q");
    check_size(generators_max_q, size, "generators_max_q");
    check_size(voltage_regulator_on, size, "voltage_regulator_on");

    // fill the data
    vm_pu_ = generators_v;
    min_q_ = generators_min_q;
    max_q_ = generators_max_q;
    for(int gen_id = 0; gen_id < size; ++gen_id){
        if (min_q_(gen_id) > max_q_(gen_id))
        {
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::init: Impossible to initialize generator min_q being above max_q for generator ";
            exc_ << gen_id;
            throw std::runtime_error(exc_.str());
        }
    }
    gen_slackbus_ = std::vector<bool>(generators_p.size(), false);
    gen_slack_weight_ = std::vector<real_type>(generators_p.size(), 0.);
    turnedoff_gen_pv_ = true;
    voltage_regulator_on_ = voltage_regulator_on;
    reset_results();
}


GeneratorContainer::StateRes GeneratorContainer::get_state() const  // osc : one side container
{
     std::vector<real_type> vm_pu(vm_pu_.begin(), vm_pu_.end());
     std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
     std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
     std::vector<bool> slack_bus = gen_slackbus_;
     std::vector<bool> voltage_regulator_on = voltage_regulator_on_;
     std::vector<real_type> slack_weight = gen_slack_weight_;
     GeneratorContainer::StateRes res(get_osc_state(),  // osc : one side container
                                      turnedoff_gen_pv_,
                                      voltage_regulator_on,
                                      vm_pu,
                                      min_q,
                                      max_q,
                                      slack_bus,
                                      slack_weight);
     return res;
}

void GeneratorContainer::set_state(GeneratorContainer::StateRes & my_state)
{
    OneSideContainer::set_osc_state(std::get<0>(my_state));
    turnedoff_gen_pv_ = std::get<1>(my_state);

    // the generators themelves
    std::vector<bool> & voltage_regulator_on = std::get<2>(my_state);
    std::vector<real_type> & vm_pu = std::get<3>(my_state);
    std::vector<real_type> & min_q = std::get<4>(my_state);
    std::vector<real_type> & max_q = std::get<5>(my_state);
    std::vector<bool> & slack_bus = std::get<6>(my_state);
    std::vector<real_type> & slack_weight = std::get<7>(my_state);

    // check sizes
    const auto size = nb();
    check_size(voltage_regulator_on, size, "voltage_regulator_on");
    check_size(vm_pu, size, "vm_pu");
    check_size(min_q, size, "min_q");
    check_size(max_q, size, "max_q");
    check_size(slack_bus, size, "slack_bus");
    check_size(slack_weight, size, "slack_weight");

    // assign data
    voltage_regulator_on_ = voltage_regulator_on;
    vm_pu_ = RealVect::Map(&vm_pu[0], vm_pu.size());
    min_q_ = RealVect::Map(&min_q[0], min_q.size());
    max_q_ = RealVect::Map(&max_q[0], max_q.size());
    gen_slackbus_ = slack_bus;
    gen_slack_weight_ = slack_weight;
    reset_results();
}

RealVect GeneratorContainer::get_slack_weights_solver(Eigen::Index nb_bus_solver, const std::vector<int> & id_grid_to_solver){
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
            exc_ << "GeneratorContainer::get_slack_weights_solver: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        if(gen_slackbus_[gen_id]) res.coeffRef(bus_id_solver) += gen_slack_weight_[gen_id];
    }
    bus_slack_weight_ = res;
    real_type sum_res = res.sum();
    res /= sum_res;
    return res;
}

void GeneratorContainer::fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const {
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    cplx_type tmp;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the load is disconnected
        if(!status_[gen_id]) continue;

        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::fillSbus: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }
        tmp = {p_mw_(gen_id), 0.};
        if(!voltage_regulator_on_[gen_id]){
            // gen is pq if voltage regulaton is off
            tmp += my_i * q_mvar_(gen_id);
        }
        Sbus.coeffRef(bus_id_solver) += tmp;
    }
}

void GeneratorContainer::fillpv(std::vector<int> & bus_pv,
                                std::vector<bool> & has_bus_been_added,
                                const Eigen::VectorXi & slack_bus_id_solver,
                                const std::vector<int> & id_grid_to_solver) const
{
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        // gen is purposedly not pv
        if (!voltage_regulator_on_[gen_id]) continue;  

        // in this case turned off generators are not pv
        // except the slack that can have a target of 0MW but is still "on"
        // no matter what
        bool gen_pseudo_off = p_mw_(gen_id) == 0.;
        // if (gen_slack_weight_[gen_id] != 0.) gen_pseudo_off = false;  // useless: slack is not PV anyway
        if ((!turnedoff_gen_pv_) && gen_pseudo_off) continue;  
        bus_id_me = bus_id_(gen_id);
        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::fillpv: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to a disconnected bus while being connected to the grid.";
            throw std::runtime_error(exc_.str());
        }

        if(is_in_vect(bus_id_solver, slack_bus_id_solver)) continue;  // slack bus is not PV
        if(has_bus_been_added[bus_id_solver]) continue; // i already added this bus
        bus_pv.push_back(bus_id_solver);
        has_bus_been_added[bus_id_solver] = true;  // don't add it a second time
    }
}

void GeneratorContainer::get_vm_for_dc(RealVect & Vm){
    const int nb_gen = nb();
    int bus_id_me;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;

        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) continue;  // in this case turned off generators are not pv

        bus_id_me = bus_id_(gen_id);
        real_type tmp = vm_pu_(gen_id);
        if(tmp != 0.) Vm(bus_id_me) = tmp;
    }
}

void GeneratorContainer::_change_p(int gen_id, real_type new_p, bool my_status, SolverControl & solver_control)
{
    if(!turnedoff_gen_pv_){
        // if turned off generators (including these with p==0)
        // are not pv, if we change the active generation, it changes
        // the list of pv buses, so I need to refactorize the solver
        // on the other hand, if all generators are pv then I do not need to refactorize in this case
        if((p_mw_(gen_id) == 0. && new_p != 0.) || 
           (p_mw_(gen_id) != 0. && new_p == 0.)){
            solver_control.tell_pv_changed();
           }
    }
}

void GeneratorContainer::change_v(int gen_id, real_type new_v_pu, SolverControl & solver_control)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::change_p: Impossible to change the voltage setpoint of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    change_v_nothrow(gen_id, new_v_pu, solver_control);
}

void GeneratorContainer::change_v_nothrow(int gen_id, real_type new_v_pu, SolverControl & solver_control)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if (vm_pu_(gen_id) != new_v_pu) solver_control.tell_v_changed();
    vm_pu_(gen_id) = new_v_pu;
}

void GeneratorContainer::set_vm(CplxVect & V, const std::vector<int> & id_grid_to_solver) const
{
    const int nb_gen = nb();
    int bus_id_me, bus_id_solver;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        //  i don't do anything if the generator is disconnected
        if(!status_[gen_id]) continue;
        
        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv

        // pseudo off generator (with p == 0)
        bool pseudo_off = p_mw_(gen_id) == 0.;
        // though "pseudo off" a slack is still "PV"
        if(gen_slack_weight_[gen_id] != 0.) pseudo_off = false;  
        if ((!turnedoff_gen_pv_) && pseudo_off) continue;  // in this case turned off generators are not pv

        bus_id_me = bus_id_(gen_id);
        // std::cout << "gen_id " << gen_id << std::endl;
        // std::cout << "id_grid_to_solver.size() " << id_grid_to_solver.size() << std::endl;
        // std::cout << "bus_id_me " << bus_id_me << std::endl;
        // std::cout << "======== " << std::endl;

        bus_id_solver = id_grid_to_solver[bus_id_me];
        if(bus_id_solver == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "GeneratorContainer::set_vm: Generator with id ";
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

Eigen::VectorXi GeneratorContainer::get_slack_bus_id() const{
    std::vector<int> tmp;
    tmp.reserve(gen_slackbus_.size());
    Eigen::VectorXi res;
    const auto nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
        if(gen_slackbus_[gen_id]){
            const auto my_bus = bus_id_(gen_id);
            // do not add twice the same "slack bus"
            if(!is_in_vect(my_bus, tmp)) tmp.push_back(my_bus);
        }
    }
    if(tmp.empty()) throw std::runtime_error("GeneratorContainer::get_slack_bus_id: no generator are tagged slack bus for this grid.");
    res = Eigen::VectorXi::Map(tmp.data(), tmp.size());  // force the copy of the data apparently
    return res;
}

void GeneratorContainer::set_p_slack(const RealVect& node_mismatch,
                                     const std::vector<int> & id_grid_to_solver)
{
    if(bus_slack_weight_.size() == 0){
        // TODO DEBUG MODE: perform this check only in debug mode
        throw std::runtime_error("Generator::set_p_slack: Impossible to set the active value of generators for the slack bus: no known slack (you should haved called Generator::get_slack_weights_solver first)");
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
        res_p_(gen_id) += node_mismatch(bus_id_solver) * my_contrib_slack / total_contrib_slack;
    }
}

void GeneratorContainer::init_q_vector(int nb_bus,
                                       Eigen::VectorXi & total_gen_per_bus,
                                       RealVect & total_q_min_per_bus,
                                       RealVect & total_q_max_per_bus) const
{
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;

        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) continue;  // in this case turned off generators are not pv
        
        int bus_id = bus_id_(gen_id);
        total_q_min_per_bus(bus_id) += min_q_(gen_id);
        total_q_max_per_bus(bus_id) += max_q_(gen_id);
        total_gen_per_bus(bus_id) += 1;
    }
}

void GeneratorContainer::set_q(const RealVect & reactive_mismatch,
                               const std::vector<int> & id_grid_to_solver,
                               bool ac,
                               const Eigen::VectorXi & total_gen_per_bus,
                               const RealVect & total_q_min_per_bus,
                               const RealVect & total_q_max_per_bus)
{
    const int nb_gen = nb();
    if(!ac){
        // do not consider Q values in dc mode
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id) res_q_(gen_id) = 0.;
        return;
    }
    
    real_type eps_q = 1e-8;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]){
            // set at 0 for disconnected generators
            res_q_(gen_id) = 0.;
            continue;  
        }
        real_type real_q = 0.;
        if (!voltage_regulator_on_[gen_id]){
            // gen is purposedly not pv
            res_q_(gen_id) = 0.;
            continue;
        } 
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) {
            // in this case turned off generators are not pv
            res_q_(gen_id) = 0.;
            continue;
        }  

        int bus_id = bus_id_(gen_id);
        const auto bus_solver = id_grid_to_solver[bus_id];
        // TODO DEBUG MODE: check that the bus is correct!
        real_type q_to_absorb = reactive_mismatch[bus_solver];
        real_type max_q_me = max_q_(gen_id);
        real_type min_q_me = min_q_(gen_id);
        real_type max_q_bus = total_q_max_per_bus(bus_id);
        real_type min_q_bus = total_q_min_per_bus(bus_id);
        int nb_gen_with_me = total_gen_per_bus(bus_id);
        if(nb_gen_with_me == 1){
            real_q = q_to_absorb;
        }else{
            real_type ratio = (max_q_me - min_q_me + eps_q) / (max_q_bus - min_q_bus + nb_gen_with_me * eps_q) ;
            real_q = q_to_absorb * ratio ;
        }
        res_q_(gen_id) = real_q;
    }
}

void GeneratorContainer::update_slack_weights(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > could_be_slack,
                                              SolverControl & solver_control)
{
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(could_be_slack(gen_id) && status_[gen_id]){
            // gen is connected and participate to the slack
            if(p_mw_(gen_id) > 0.){
                // gen is properly connected
                if(!gen_slackbus_[gen_id]) solver_control.tell_slack_participate_changed(); // it was not in the slack before, so I need to reset the solver
                add_slackbus(gen_id, p_mw_(gen_id), solver_control);

            }else{
                // gen is now "turned off" (p_mw=0.)
                if(gen_slackbus_[gen_id]) solver_control.tell_slack_participate_changed();  // it was in the slack before, so I need to reset the solver
                remove_slackbus(gen_id, solver_control);
            }
        }else{
            if(gen_slackbus_[gen_id]) solver_control.tell_slack_participate_changed();  // it was in the slack before, I need to reset the solver
            remove_slackbus(gen_id, solver_control);
        }
    }
}
