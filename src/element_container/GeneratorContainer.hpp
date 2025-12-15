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

#include "Utils.hpp"
#include "OneSideContainer_PQ.hpp"


// class _GeneratorImpl
// {
//     // implementation of the interface
//     protected:
//         void _deactivate(int gen_id, SolverControl & solver_control) {};
//         void _reactivate(int gen_id, SolverControl & solver_control) {};
//         void _change_bus(int gen_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {};
//         void _compute_results(const Eigen::Ref<const RealVect> & Va,
//                                     const Eigen::Ref<const RealVect> & Vm,
//                                     const Eigen::Ref<const CplxVect> & V,
//                                     const std::vector<int> & id_grid_to_solver,
//                                     const RealVect & bus_vn_kv,
//                                     real_type sn_mva,
//                                     bool ac){};
//         void _change_p(int load_id, real_type new_p, bool my_status, SolverControl & solver_control) {};
//         void _reset_results() {};
//         void _change_p(int el_id, real_type new_p, bool my_status, SolverControl & solver_control) {};
//         void _change_q(int el_id, real_type new_p, bool my_status,SolverControl & solver_control) {};

//         // void _deactivate(int gen_id, SolverControl & solver_control) {
//         //     if (status_[gen_id]){
//         //         solver_control.tell_recompute_sbus();
//         //         solver_control.tell_pq_changed();  // bus might now be pq
//         //         if(voltage_regulator_on_[gen_id]) solver_control.tell_v_changed();
//         //         solver_control.tell_pv_changed();
//         //         if(gen_slack_weight_[gen_id] != 0. || gen_slackbus_[gen_id]){
//         //             solver_control.tell_slack_participate_changed();
//         //             solver_control.tell_slack_weight_changed();
//         //         }
//         //     }
//         //     gen_slackbus_[gen_id] = false;
//         // }

//         // void _reactivate(int gen_id, SolverControl & solver_control) {
//         //     if(!status_[gen_id]){
//         //         solver_control.tell_recompute_sbus();
//         //         // bus might change between pv / pq depending on the state of the generator
//         //         // TODO speed optim here
//         //         solver_control.tell_pq_changed();
//         //         solver_control.tell_pv_changed(); 

//         //         if(voltage_regulator_on_[gen_id]) solver_control.tell_v_changed();
//         //         if(gen_slack_weight_[gen_id] != 0. || gen_slackbus_[gen_id]){
//         //             solver_control.tell_slack_participate_changed();
//         //             solver_control.tell_slack_weight_changed();
//         //         }
//         //         if(gen_slack_weight_[gen_id] != 0.){
//         //             gen_slackbus_[gen_id] = true;
//         //         }
//         //     }
//         // }

//         // void _change_bus(int gen_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
//         //     if (new_bus_id != bus_id_[gen_id]){
//         //         if (abs(gen_slack_weight_[gen_id]) > 1e-6 || gen_slackbus_[gen_id]) solver_control.has_slack_participate_changed();
//         //         // bus might change between pv / pq depending on the state of the generator
//         //         // TODO speed optim here
//         //         solver_control.tell_pq_changed();
//         //         solver_control.tell_pv_changed(); 
//         //     }
//         // }

//         // void _compute_results(const Eigen::Ref<const RealVect> & Va,
//         //                             const Eigen::Ref<const RealVect> & Vm,
//         //                             const Eigen::Ref<const CplxVect> & V,
//         //                             const std::vector<int> & id_grid_to_solver,
//         //                             const RealVect & bus_vn_kv,
//         //                             real_type sn_mva,
//         //                             bool ac){
//         //                                 OneSideContainer_PQ<GeneratorContainer>::set_osc_pq_res_p();
//         //                             }

//         // void _change_p(int load_id, real_type new_p, bool my_status, SolverControl & solver_control);
//         // void _reset_results() {};
//         // void _change_p(int el_id, real_type new_p, bool my_status, SolverControl & solver_control) {};
//         // void _change_q(int el_id, real_type new_p, bool my_status,SolverControl & solver_control) {};
// };

// typedef OneSideContainer_PQ<_GeneratorImpl> _BaseGenClass;

/**
This class represents the list of all generators.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/gen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/gen.html#electric-model
**/
class GeneratorContainer: public OneSideContainer_PQ
{

    public:
        class GenInfo : public OneSidePQInfo
        {
            public:
                bool is_slack;
                real_type slack_weight;

                bool voltage_regulator_on;
                real_type target_vm_pu;
                real_type min_q_mvar;
                real_type max_q_mvar;

                GenInfo(const GeneratorContainer & r_data_gen, int my_id):
                OneSidePQInfo(r_data_gen, my_id),
                is_slack(false),
                slack_weight(-1.0),
                voltage_regulator_on(false),
                target_vm_pu(0.),
                min_q_mvar(0.),
                max_q_mvar(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_gen.nb()))
                    {
                        is_slack = r_data_gen.gen_slackbus_[my_id];
                        slack_weight = r_data_gen.gen_slack_weight_[my_id];

                        voltage_regulator_on = r_data_gen.voltage_regulator_on_[my_id];
                        target_vm_pu = r_data_gen.target_vm_pu_.coeff(my_id);
                        min_q_mvar = r_data_gen.min_q_.coeff(my_id);
                        max_q_mvar = r_data_gen.max_q_.coeff(my_id);
                    }
                }
        };
        typedef GenInfo DataInfo;

    private:
    typedef GenericContainerConstIterator<GeneratorContainer> GeneratorConstIterator;


    // iterator
    public:
        GeneratorConstIterator begin() const {return GeneratorConstIterator(this, 0); }
        GeneratorConstIterator end() const {return GeneratorConstIterator(this, nb()); }
        GenInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative load id.");
            }
            if(id >= nb())
            {
                throw std::range_error("Load out of bound. Not enough loads on the grid.");
            }
            return GenInfo(*this, id);
        }

    public:
        typedef std::tuple<
           OneSideContainer_PQ::StateRes,
           bool,                    // turnedoff_gen_pv_
           std::vector<bool>,       // voltage_regulator_on
           std::vector<real_type>,  // target_vm_pu_
           std::vector<real_type>,  // min_q_
           std::vector<real_type>,  // max_q_
           std::vector<bool>,       // gen_slackbus
           std::vector<real_type>   // gen_slack_weight_
        >  StateRes;
        
        GeneratorContainer():OneSideContainer_PQ(), turnedoff_gen_pv_(true){};
        GeneratorContainer(bool turnedoff_gen_pv):OneSideContainer_PQ(), turnedoff_gen_pv_(turnedoff_gen_pv) {};
        
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
                const real_type p_mw = target_p_mw_(gen_id);
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
        void update_slack_weights_by_id(Eigen::Ref<const IntVect> gen_slack_id, SolverControl & solver_control);
        
        
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
            for(const auto & el : target_vm_pu_){
                std::cout << "V " << el << std::endl;
            }
        }

    protected:
        // physical properties
        RealVect min_q_;
        RealVect max_q_;

        // input data
        std::vector<bool> voltage_regulator_on_;
        RealVect target_vm_pu_;

        // remember which generators are "slack bus"
        std::vector<bool> gen_slackbus_;  // say for each generator if it's a slack or not
        std::vector<real_type> gen_slack_weight_;

        // intermediate data
        // Eigen::VectorXi total_gen_per_bus_;
        RealVect bus_slack_weight_;  // do not sum to 1., for each node of the grid, say the raw contribution for the generator

        // different parameter of the behaviour of the class
        bool turnedoff_gen_pv_;  // are turned off generators (including one with p=0) pv ?

    protected:
    virtual void _change_p(int gen_id, real_type new_p, bool my_status, SolverControl & solver_control);

};

#endif  //GENERATORCONTAINER_H
