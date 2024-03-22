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
    int size = static_cast<int>(generators_p.size());
    GenericContainer::check_size(generators_p, size, "generators_p");
    GenericContainer::check_size(generators_v, size, "generators_v");
    GenericContainer::check_size(generators_min_q, size, "generators_min_q");
    GenericContainer::check_size(generators_max_q, size, "generators_max_q");
    GenericContainer::check_size(generators_bus_id, size, "generators_bus_id");

    p_mw_ = generators_p;
    vm_pu_ = generators_v;
    bus_id_ = generators_bus_id;
    min_q_ = generators_min_q;
    max_q_ = generators_max_q;
    if(min_q_.size() != max_q_.size())
    {
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::init: Impossible to initialize generator with generators_min_q of size ";
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
            exc_ << "GeneratorContainer::init: Impossible to initialize generator min_q being above max_q for generator ";
            exc_ << gen_id;
            throw std::runtime_error(exc_.str());
        }
    }
    status_ = std::vector<bool>(generators_p.size(), true);
    gen_slackbus_ = std::vector<bool>(generators_p.size(), false);
    gen_slack_weight_ = std::vector<real_type>(generators_p.size(), 0.);
    turnedoff_gen_pv_ = true;
    voltage_regulator_on_ = std::vector<bool>(generators_p.size(), true);
    q_mvar_ = RealVect::Zero(generators_p.size());
    reset_results();
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
    init(generators_p, generators_v, generators_min_q, generators_max_q, generators_bus_id);
    int size = static_cast<int>(generators_p.size());
    GenericContainer::check_size(generators_q, size, "generators_q");
    GenericContainer::check_size(voltage_regulator_on, size, "voltage_regulator_on");
    voltage_regulator_on_ = voltage_regulator_on;
    q_mvar_ = generators_q;
}


GeneratorContainer::StateRes GeneratorContainer::get_state() const
{
     std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
     std::vector<real_type> vm_pu(vm_pu_.begin(), vm_pu_.end());
     std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
     std::vector<real_type> min_q(min_q_.begin(), min_q_.end());
     std::vector<real_type> max_q(max_q_.begin(), max_q_.end());
     std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
     std::vector<bool> status = status_;
     std::vector<bool> slack_bus = gen_slackbus_;
     std::vector<bool> voltage_regulator_on = voltage_regulator_on_;
     std::vector<real_type> slack_weight = gen_slack_weight_;
     GeneratorContainer::StateRes res(names_, turnedoff_gen_pv_, voltage_regulator_on,
                           p_mw, vm_pu, q_mvar,
                           min_q, max_q, bus_id, status, slack_bus, slack_weight);
     return res;
}

void GeneratorContainer::set_state(GeneratorContainer::StateRes & my_state)
{
    names_ = std::get<0>(my_state);
    turnedoff_gen_pv_ = std::get<1>(my_state);

    // the generators themelves
    std::vector<bool> & voltage_regulator_on = std::get<2>(my_state);
    std::vector<real_type> & p_mw = std::get<3>(my_state);
    std::vector<real_type> & vm_pu = std::get<4>(my_state);
    std::vector<real_type> & q_mvar = std::get<5>(my_state);
    std::vector<real_type> & min_q = std::get<6>(my_state);
    std::vector<real_type> & max_q = std::get<7>(my_state);
    std::vector<int> & bus_id = std::get<8>(my_state);
    std::vector<bool> & status = std::get<9>(my_state);
    std::vector<bool> & slack_bus = std::get<10>(my_state);
    std::vector<real_type> & slack_weight = std::get<11>(my_state);
    // TODO check sizes

    // input data
    voltage_regulator_on_ = voltage_regulator_on;
    p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
    vm_pu_ = RealVect::Map(&vm_pu[0], vm_pu.size());
    q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
    min_q_ = RealVect::Map(&min_q[0], min_q.size());
    max_q_ = RealVect::Map(&max_q[0], max_q.size());
    bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
    status_ = status;
    gen_slackbus_ = slack_bus;
    gen_slack_weight_ = slack_weight;
    reset_results();
}

RealVect GeneratorContainer::get_slack_weights(Eigen::Index nb_bus_solver, const std::vector<int> & id_grid_to_solver){
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
            exc_ << "GeneratorContainer::get_slack_weights: Generator with id ";
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
        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) continue;  // in this case turned off generators are not pv

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

void GeneratorContainer::compute_results(const Eigen::Ref<const RealVect> & Va,
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

void GeneratorContainer::reset_results(){
    res_p_ = RealVect(nb());  // in MW
    res_q_ = RealVect(nb());  // in MVar
    res_v_ = RealVect(nb());  // in kV
    res_theta_ = RealVect(nb());  // in deg
    // bus_slack_weight_ = RealVect();
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

void GeneratorContainer::change_p(int gen_id, real_type new_p, SolverControl & solver_control)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::change_p: Impossible to change the active value of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
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
    if (p_mw_(gen_id) != new_p){
        solver_control.tell_recompute_sbus();
        p_mw_(gen_id) = new_p;
    }
}

void GeneratorContainer::change_q(int gen_id, real_type new_q, SolverControl & solver_control)
{
    bool my_status = status_.at(gen_id); // and this check that load_id is not out of bound
    if(!my_status)
    {
        // TODO DEBUG MODE only this in debug mode
        std::ostringstream exc_;
        exc_ << "GeneratorContainer::change_q: Impossible to change the reactive value of a disconnected generator (check gen. id ";
        exc_ << gen_id;
        exc_ << ")";
        throw std::runtime_error(exc_.str());
    }
    // TODO DEBUG MODE : raise an error if generator is regulating voltage, maybe ? 
    // this would have not effect
    if (q_mvar_(gen_id) != new_q){
        solver_control.tell_recompute_sbus();
        q_mvar_(gen_id) = new_q;
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
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) continue;  // in this case turned off generators are not pv

        bus_id_me = bus_id_(gen_id);
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
        throw std::runtime_error("Generator::set_p_slack: Impossible to set the active value of generators for the slack bus: no known slack (you should haved called Generator::get_slack_weights first)");
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
        // std::cout << "gen_id " << gen_id << " my_contrib_slack " << my_contrib_slack << ", " << total_contrib_slack << " node_mismatch " << node_mismatch(bus_id_solver) << std::endl;
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
    // res_q_ = RealVect::Constant(nb_gen, 0.);
    if(!ac)  
    {
        // do not consider Q values in dc mode
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id) res_q_(gen_id) = 0.;
        return;
    }
    
    real_type eps_q = 1e-8;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        real_type real_q = 0.;
        if(!status_[gen_id]){
            // set at 0 for disconnected generators
            res_q_(gen_id) = 0.;
            continue;  
        }

        if (!voltage_regulator_on_[gen_id]) continue;  // gen is purposedly not pv
        if ((!turnedoff_gen_pv_) && p_mw_(gen_id) == 0.) continue;  // in this case turned off generators are not pv

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

void GeneratorContainer::reconnect_connected_buses(std::vector<bool> & bus_status) const {
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;
        const auto my_bus = bus_id_(gen_id);
        if(my_bus == _deactivated_bus_id){
            // TODO DEBUG MODE only this in debug mode
            std::ostringstream exc_;
            exc_ << "Generator::reconnect_connected_buses: Generator with id ";
            exc_ << gen_id;
            exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_gen(...)` ?.";
            throw std::runtime_error(exc_.str());
        }
        bus_status[my_bus] = true;  // this bus is connected
    }
}

void GeneratorContainer::gen_p_per_bus(std::vector<real_type> & res) const
{
    const int nb_gen = nb();
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;
        const auto my_bus = bus_id_(gen_id);
        if (p_mw_(gen_id) > 0.) res[my_bus] += p_mw_(gen_id);
    }
}

void GeneratorContainer::disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
    const int nb_gen = nb();
    SolverControl unused_solver_control;
    for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
    {
        if(!status_[gen_id]) continue;
        const auto my_bus = bus_id_(gen_id);
        if(!busbar_in_main_component[my_bus]){
            deactivate(gen_id, unused_solver_control);
        }
    }    
}
