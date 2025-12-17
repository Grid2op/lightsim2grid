// Copyright (c) 2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DCLINECONTAINER_H
#define DCLINECONTAINER_H

#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "SubstationContainer.hpp"
#include "TwoSidesContainer.hpp"
#include "GeneratorContainer.hpp"


class DCLineContainer;
class DCLineInfo : public TwoSidesContainer<GeneratorContainer>::TwoSidesInfo
{
    public:
        // members
        real_type target_p_1_mw;
        real_type p_2_mw;
        real_type target_vm_1_pu;
        real_type target_vm_2_pu;
        real_type loss_pct;
        real_type loss_mw;
        GenInfo gen_side_1;
        GenInfo gen_side_2;

        DCLineInfo(const DCLineContainer & r_data_dcline, int my_id);
};


class DCLineContainer : public TwoSidesContainer<GeneratorContainer>, public IteratorAdder<DCLineContainer, DCLineInfo>
{
    friend class DCLineInfo;

    public:
        typedef DCLineInfo DataInfo;

    // underlying generators are not pv when powerline is off
    DCLineContainer(){
        SolverControl solver_control_not_used;
        side_1_.turnedoff_no_pv(solver_control_not_used);
        side_2_.turnedoff_no_pv(solver_control_not_used);
    };

        // pickle
        typedef std::tuple<
                TwoSidesContainer<GeneratorContainer>::StateRes,
                std::vector<double>, // loss_percent
                std::vector<double> // loss_mw
                >  StateRes;
        DCLineContainer::StateRes get_state() const;
        void set_state(DCLineContainer::StateRes & my_state);

        // TODO min_p, max_p
        void init(const Eigen::VectorXi & branch_from_id,
                const Eigen::VectorXi & branch_to_id,
                const RealVect & p_mw,
                const RealVect & loss_percent,
                const RealVect & loss_mw,
                const RealVect & vm_or_pu,
                const RealVect & vm_ex_pu,
                const RealVect & min_q_or,
                const RealVect & max_q_or,
                const RealVect & min_q_ex,
                const RealVect & max_q_ex
                );

        // accessor / modifiers
        void deactivate(int dcline_id, SolverControl & solver_control) {  // TODO this in TwoSidesCOntainer !
            _generic_deactivate(dcline_id, status_global_);
            side_1_.deactivate(dcline_id, solver_control);
            side_2_.deactivate(dcline_id, solver_control);
        }
        void reactivate(int dcline_id, SolverControl & solver_control) {  // TODO this in TwoSidesCOntainer !
            _generic_reactivate(dcline_id, status_global_);
            side_1_.reactivate(dcline_id, solver_control);
            side_2_.reactivate(dcline_id, solver_control);
        }

        // for buses only connected through dc line, i don't add them
        // they are not in the same "connected component"
        virtual void reconnect_connected_buses(SubstationContainer & Substation) const {
            // from_gen_.reconnect_connected_buses(bus_status);
            // to_gen_.reconnect_connected_buses(bus_status);
        }

        // for buses only connected through dc line, i don't add them
        // they are not in the same "connected component"
        virtual void get_graph(std::vector<Eigen::Triplet<real_type> > & res) const {};

        // virtual void nb_line_end(std::vector<int> & res) const;
        
        real_type get_qmin_or(int dcline_id) {return side_1_.get_qmin(dcline_id);}
        real_type get_qmax_or(int dcline_id) {return  side_1_.get_qmax(dcline_id);}
        real_type get_qmin_ex(int dcline_id) {return side_2_.get_qmin(dcline_id);}
        real_type get_qmax_ex(int dcline_id) {return  side_2_.get_qmax(dcline_id);}

        real_type get_to_mw(int dcline_id, real_type from_mw){
            // TODO set it to load convention instead of gen convention as in lightsim2grid
            real_type new_p_ext = from_mw >= 0 ? 
                                -(from_mw + loss_mw_(dcline_id)) / (1.0 - 0.01 * loss_percent_(dcline_id)) :
                                -from_mw * (1.0 - 0.01 * loss_percent_(dcline_id)) - loss_mw_(dcline_id)
                                ;
            return new_p_ext;
        }
        void change_p(int dcline_id, real_type new_p, SolverControl & sovler_control){
            side_1_.change_p(dcline_id, -1.0 * new_p, sovler_control);

            side_2_.change_p(dcline_id, -1.0 * get_to_mw(dcline_id, new_p), sovler_control);
        }
        void change_v_side_1(int dcline_id, real_type new_v_pu, SolverControl & sovler_control){
            side_1_.change_v(dcline_id, new_v_pu, sovler_control);
        }
        void change_v_side_2(int dcline_id, real_type new_v_pu, SolverControl & sovler_control){
            side_2_.change_v(dcline_id, new_v_pu, sovler_control);
        }

        // solver stuff
        virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const{
            side_1_.fillSbus(Sbus, id_grid_to_solver, ac);   
            side_2_.fillSbus(Sbus, id_grid_to_solver, ac);   
        } 

        virtual void fillpv(std::vector<int>& bus_pv,
                            std::vector<bool> & has_bus_been_added,
                            const Eigen::VectorXi & slack_bus_id_solver,
                            const std::vector<int> & id_grid_to_solver) const {
            side_1_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_grid_to_solver);   
            side_2_.fillpv(bus_pv, has_bus_been_added, slack_bus_id_solver, id_grid_to_solver);   
        }
        void init_q_vector(int nb_bus,
                        Eigen::VectorXi & total_gen_per_bus,
                        RealVect & total_q_min_per_bus,
                        RealVect & total_q_max_per_bus) const{
            side_1_.init_q_vector(nb_bus, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
            side_2_.init_q_vector(nb_bus, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
        }
        void compute_results(const Eigen::Ref<const RealVect> & Va,
                            const Eigen::Ref<const RealVect> & Vm,
                            const Eigen::Ref<const CplxVect> & V,
                            const std::vector<int> & id_grid_to_solver,
                            const RealVect & bus_vn_kv,
                            real_type sn_mva,
                            bool ac){
            side_1_.compute_results(Va, Vm, V,
                                    id_grid_to_solver,
                                    bus_vn_kv, sn_mva, ac);
            side_2_.compute_results(Va, Vm, V,
                                    id_grid_to_solver,
                                    bus_vn_kv, sn_mva, ac);
        }

        void reset_results(){
            reset_results_tsc();
        }

        void set_q(const RealVect & reactive_mismatch,
                const std::vector<int> & id_grid_to_solver,
                bool ac,
                const Eigen::VectorXi & total_gen_per_bus,
                const RealVect & total_q_min_per_bus,
                const RealVect & total_q_max_per_bus){
            // TODO set it to load convention instead of gen convention as in lightsim2grid
            side_1_.set_q(reactive_mismatch, id_grid_to_solver, ac, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
            side_2_.set_q(reactive_mismatch, id_grid_to_solver, ac, total_gen_per_bus, total_q_min_per_bus, total_q_max_per_bus);
        }
        void get_vm_for_dc(RealVect & Vm){
            side_1_.get_vm_for_dc(Vm);
            side_2_.get_vm_for_dc(Vm);
        }
        void set_vm_or(CplxVect & V, const std::vector<int> & id_grid_to_solver) const{
            side_1_.set_vm(V, id_grid_to_solver);
        }
        void set_vm_ex(CplxVect & V, const std::vector<int> & id_grid_to_solver) const{
            side_2_.set_vm(V, id_grid_to_solver);
        }

        /**
        this functions makes sure that the voltage magnitude of every connected bus is properly used to initialize
        the ac powerflow
        **/
        void set_vm(CplxVect & V, const std::vector<int> & id_grid_to_solver) const{
            side_1_.set_vm(V, id_grid_to_solver);
            side_2_.set_vm(V, id_grid_to_solver);
        }

    protected:
        // it is modeled as 2 generators that are "linked" together
        // see https://pandapower.readthedocs.io/en/v2.0.1/elements/dcline.html#electric-model
        RealVect loss_percent_;
        RealVect loss_mw_;

};

#endif  //DCLINECONTAINER_H