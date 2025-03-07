// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GENERATORCONTAINER_H
#define GENERATORCONTAINER_H

#include <iostream>
#include <vector> 

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "OneSideContainer.h"

/**
This class represents the list of all generators.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/gen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/gen.html#electric-model
**/
class GeneratorContainer: public OneSideContainer
{
    public:
    class GenInfo
    {
        public:
            // members
            // TODO add some const here (value should not be changed !) !!!
            int id;  // id of the generator
            std::string name;
            bool connected;
            int bus_id;
            bool is_slack;
            real_type slack_weight;

            bool voltage_regulator_on;
            real_type target_p_mw;
            real_type target_vm_pu;
            real_type target_q_mvar;
            real_type min_q_mvar;
            real_type max_q_mvar;
            bool has_res;
            real_type res_p_mw;
            real_type res_q_mvar;
            real_type res_v_kv;
            real_type res_theta_deg;

            GenInfo(const GeneratorContainer & r_data_gen, int my_id):
            id(-1),
            name(""),
            connected(false),
            bus_id(_deactivated_bus_id),
            is_slack(false),
            slack_weight(-1.0),
            voltage_regulator_on(false),
            target_p_mw(0.),
            target_vm_pu(0.),
            target_q_mvar(0.),
            min_q_mvar(0.),
            max_q_mvar(0.),
            has_res(false),
            res_p_mw(0.),
            res_q_mvar(0.),
            res_v_kv(0.),
            res_theta_deg(0.)
            {
                if((my_id >= 0) & (my_id < r_data_gen.nb()))
                {
                    id = my_id;
                    if(r_data_gen.names_.size()){
                        name = r_data_gen.names_[my_id];
                    }
                    connected = r_data_gen.status_[my_id];
                    if(connected) bus_id = r_data_gen.bus_id_[my_id];
                    is_slack = r_data_gen.gen_slackbus_[my_id];
                    slack_weight = r_data_gen.gen_slack_weight_[my_id];

                    voltage_regulator_on = r_data_gen.voltage_regulator_on_[my_id];
                    target_p_mw = r_data_gen.p_mw_.coeff(my_id);
                    target_vm_pu = r_data_gen.vm_pu_.coeff(my_id);
                    target_q_mvar = r_data_gen.q_mvar_.coeff(my_id);
                    min_q_mvar = r_data_gen.min_q_.coeff(my_id);
                    max_q_mvar = r_data_gen.max_q_.coeff(my_id);

                    has_res = r_data_gen.res_p_.size() > 0;
                    if(has_res)
                    {
                        res_p_mw = r_data_gen.res_p_.coeff(my_id);
                        res_q_mvar = r_data_gen.res_q_.coeff(my_id);
                        res_v_kv = r_data_gen.res_v_.coeff(my_id);
                        res_theta_deg = r_data_gen.res_theta_.coeff(my_id);
                    }
                }
            }
    };
    typedef GenInfo DataInfo;

    private:
    typedef GenericContainerConstIterator<GeneratorContainer> GeneratorConstIterator;

    public:
    typedef std::tuple<
       OneSideContainer::StateRes,
       bool,
       std::vector<bool>, // voltage_regulator_on
       std::vector<real_type>, // vm_pu_
       std::vector<real_type>, // min_q_
       std::vector<real_type>, // max_q_
       std::vector<bool>, // gen_slackbus
       std::vector<real_type> // gen_slack_weight_
       >  StateRes;

    GeneratorContainer():OneSideContainer(), turnedoff_gen_pv_(true){};
    GeneratorContainer(bool turnedoff_gen_pv):OneSideContainer(), turnedoff_gen_pv_(turnedoff_gen_pv) {};

    // TODO add pmin and pmax here !
    void init(const RealVect & generators_p,
              const RealVect & generators_v,
              const RealVect & generators_min_q,
              const RealVect & generators_max_q,
              const Eigen::VectorXi & generators_bus_id
              );

    void init_full(const RealVect & generators_p,
                   const RealVect & generators_v,
                   const RealVect & generators_q,
                   const std::vector<bool> & voltage_regulator_on,
                   const RealVect & generators_min_q,
                   const RealVect & generators_max_q,
                   const Eigen::VectorXi & generators_bus_id
                   );

    // iterator
    typedef GeneratorConstIterator const_iterator_type;
    const_iterator_type begin() const {return GeneratorConstIterator(this, 0); }
    const_iterator_type end() const {return GeneratorConstIterator(this, nb()); }
    GenInfo operator[](int id) const
    {
        if(id < 0)
        {
            throw std::range_error("You cannot ask for a negative generator");
        }
        if(id >= nb())
        {
            throw std::range_error("Generator out of bound. Not enough generator on the grid.");
        }
        return GenInfo(*this, id);
    }

    // pickle
    GeneratorContainer::StateRes get_state() const;
    void set_state(GeneratorContainer::StateRes & my_state );

    // slack handling
    /**
    we suppose that the data are correct (ie gen_id in the proper range, and weight > 0.)
    This is checked in GridModel, and not at this stage
    **/
    void add_slackbus(int gen_id, real_type weight, SolverControl & solver_control){
        // TODO DEBUG MODE
        if(weight <= 0.) throw std::runtime_error("GeneratorContainer::add_slackbus Cannot assign a negative (<=0) weight to the slack bus.");
        if(!gen_slackbus_[gen_id]) solver_control.tell_slack_participate_changed();
        gen_slackbus_[gen_id] = true;
        if(gen_slack_weight_[gen_id] != weight) solver_control.tell_slack_weight_changed();
        gen_slack_weight_[gen_id] = weight;
    }
    void remove_slackbus(int gen_id, SolverControl & solver_control){
        if(gen_slackbus_[gen_id]) solver_control.tell_slack_participate_changed();
        if(gen_slack_weight_[gen_id] != 0.) solver_control.tell_slack_weight_changed();
        gen_slackbus_[gen_id] = false;
        gen_slack_weight_[gen_id] = 0.;
    }
    void remove_all_slackbus(){
        const int nb_gen = nb();
        SolverControl unused_solver_control;
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
        {
            remove_slackbus(gen_id, unused_solver_control);
        }
    }

    // returns only the gen_id with the highest p that is connected to this bus !
    int assign_slack_bus(int slack_bus_id,
                         const std::vector<real_type> & gen_p_per_bus,
                         SolverControl & solver_control){
        const int nb_gen = nb();
        int res_gen_id = -1;
        real_type max_p = -1.;
        for(int gen_id = 0; gen_id < nb_gen; ++gen_id)
        {
            if(!status_[gen_id]) continue;
            if(bus_id_(gen_id) != slack_bus_id) continue;
            const real_type p_mw = p_mw_(gen_id);
            if (p_mw > 0.) add_slackbus(gen_id, p_mw / gen_p_per_bus[slack_bus_id], solver_control);
            if((p_mw > max_p) || (res_gen_id == -1) ){
                res_gen_id = gen_id;
                max_p = p_mw;
            }
        }
        // TODO DEBUG MODE
        if(res_gen_id == -1) throw std::runtime_error("GeneratorContainer::assign_slack_bus No generator connected to the desired buses");
        return res_gen_id;
    }

    /**
    Retrieve the normalized (=sum to 1.000) slack weights for all the buses
    **/
    RealVect get_slack_weights_solver(Eigen::Index nb_bus_solver, const std::vector<int> & id_grid_to_solver);

    Eigen::VectorXi get_slack_bus_id() const;
    void set_p_slack(const RealVect& node_mismatch, const std::vector<int> & id_grid_to_solver);

    // modification
    void turnedoff_no_pv(SolverControl & solver_control){
        solver_control.tell_slack_participate_changed();
        solver_control.tell_slack_weight_changed();
        turnedoff_gen_pv_=false;  // turned off generators are not pv. This is NOT the default.
        }  
    void turnedoff_pv(SolverControl & solver_control){
        solver_control.tell_slack_participate_changed();
        solver_control.tell_slack_weight_changed();
        turnedoff_gen_pv_=true;  // turned off generators are pv. This is the default.
        }  
    bool get_turnedoff_gen_pv() const {return turnedoff_gen_pv_;}
    void update_slack_weights(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > could_be_slack,
                              SolverControl & solver_control);


    real_type get_qmin(int gen_id) {return min_q_.coeff(gen_id);}
    real_type get_qmax(int gen_id) {return max_q_.coeff(gen_id);}

    void change_v(int gen_id, real_type new_v_pu, SolverControl & solver_control);
    void change_v_nothrow(int gen_id, real_type new_v_pu, SolverControl & solver_control);

    virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const;
    virtual void fillpv(std::vector<int>& bus_pv,
                        std::vector<bool> & has_bus_been_added,
                        const Eigen::VectorXi & slack_bus_id_solver,
                        const std::vector<int> & id_grid_to_solver) const;
    void init_q_vector(int nb_bus,
                       Eigen::VectorXi & total_gen_per_bus,
                       RealVect & total_q_min_per_bus,
                       RealVect & total_q_max_per_bus) const; // delta_q_per_gen_

    void set_q(const RealVect & reactive_mismatch,
               const std::vector<int> & id_grid_to_solver,
               bool ac,
               const Eigen::VectorXi & total_gen_per_bus,
               const RealVect & total_q_min_per_bus,
               const RealVect & total_q_max_per_bus);

    void get_vm_for_dc(RealVect & Vm);

    /**
    this functions makes sure that the voltage magnitude of every connected bus is properly used to initialize
    the ac powerflow
    **/
    void set_vm(CplxVect & V, const std::vector<int> & id_grid_to_solver) const;

    void cout_v(){
        for(const auto & el : vm_pu_){
            std::cout << "V " << el << std::endl;
        }
    }


    protected:

    virtual void _deactivate(int gen_id, SolverControl & solver_control) {
        if (status_[gen_id]){
            solver_control.tell_recompute_sbus();
            solver_control.tell_pq_changed();  // bus might now be pq
            if(voltage_regulator_on_[gen_id]) solver_control.tell_v_changed();
            solver_control.tell_pv_changed();
            if(gen_slack_weight_[gen_id] != 0. || gen_slackbus_[gen_id]){
                solver_control.tell_slack_participate_changed();
                solver_control.tell_slack_weight_changed();
            }
        }
        gen_slackbus_[gen_id] = false;
    }

    virtual void _reactivate(int gen_id, SolverControl & solver_control) {
        if(!status_[gen_id]){
            solver_control.tell_recompute_sbus();
            // bus might change between pv / pq depending on the state of the generator
            // TODO speed optim here
            solver_control.tell_pq_changed();
            solver_control.tell_pv_changed(); 

            if(voltage_regulator_on_[gen_id]) solver_control.tell_v_changed();
            if(gen_slack_weight_[gen_id] != 0. || gen_slackbus_[gen_id]){
                solver_control.tell_slack_participate_changed();
                solver_control.tell_slack_weight_changed();
            }
            if(gen_slack_weight_[gen_id] != 0.){
                gen_slackbus_[gen_id] = true;
            }
        }
    }

    virtual void _change_bus(int gen_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        if (new_bus_id != bus_id_[gen_id]){
            if (gen_slack_weight_[gen_id] != 0. || gen_slackbus_[gen_id]) solver_control.has_slack_participate_changed();
            // bus might change between pv / pq depending on the state of the generator
            // TODO speed optim here
            solver_control.tell_pq_changed();
            solver_control.tell_pv_changed(); 
        }
    }

    virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                  const Eigen::Ref<const RealVect> & Vm,
                                  const Eigen::Ref<const CplxVect> & V,
                                  const std::vector<int> & id_grid_to_solver,
                                  const RealVect & bus_vn_kv,
                                  real_type sn_mva,
                                  bool ac){
                                    set_osc_res_p();
                                  }

    virtual void _change_p(int load_id, real_type new_p, bool my_status, SolverControl & solver_control);
    
    protected:
        // physical properties
        RealVect min_q_;
        RealVect max_q_;

        // input data
        std::vector<bool> voltage_regulator_on_;
        RealVect vm_pu_;

        // remember which generators are "slack bus"
        std::vector<bool> gen_slackbus_;  // say for each generator if it's a slack or not
        std::vector<real_type> gen_slack_weight_;

        // intermediate data
        // Eigen::VectorXi total_gen_per_bus_;
        RealVect bus_slack_weight_;  // do not sum to 1., for each node of the grid, say the raw contribution for the generator

        // different parameter of the behaviour of the class
        bool turnedoff_gen_pv_;  // are turned off generators (including one with p=0) pv ?
};

#endif  //GENERATORCONTAINER_H
