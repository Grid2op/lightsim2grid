// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GRIDMODEL_H
#define GRIDMODEL_H

#include <iostream>
#include <vector>
// #include <set>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI

#include "Utils.h"

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

// import data classes
#include "element_container/GenericContainer.h"
#include "element_container/LineContainer.h"
#include "element_container/ShuntContainer.h"
#include "element_container/TrafoContainer.h"
#include "element_container/LoadContainer.h"
#include "element_container/GeneratorContainer.h"
#include "element_container/SGenContainer.h"
#include "element_container/DCLineContainer.h"

// import newton raphson solvers using different linear algebra solvers
#include "ChooseSolver.h"
// class ChooseSolver;
// enum class SolverType;

//TODO implement a BFS check to make sure the Ymatrix is "connected" [one single component]
class GridModel : public GenericContainer
{
    public:
        typedef Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> IntVectRowMaj;

        typedef std::tuple<
                int, // version major
                int, // version medium
                int, // version minor
                std::vector<int>, // ls_to_orig
                real_type,  // init_vm_pu
                real_type, //sn_mva
                std::vector<real_type>,  // bus_vn_kv
                std::vector<bool>,  // bus_status
                // powerlines
                LineContainer::StateRes ,
                // shunts
                ShuntContainer::StateRes,
                // trafos
                TrafoContainer::StateRes,
                // gens
                GeneratorContainer::StateRes,
                // loads
                LoadContainer::StateRes,
                // static generators
                SGenContainer::StateRes,
                // storage units
                LoadContainer::StateRes,
                //dc lines
                DCLineContainer::StateRes,
                // grid2op specific
                int, // n_sub
                int, // max_nb_bus_per_sub
                std::vector<int>,  // load_pos_topo_vect_
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,  
                std::vector<int>,  // load_to_subid_
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>,
                std::vector<int>
                >  StateRes;

        GridModel():
          timer_last_ac_pf_(0.),
          timer_last_dc_pf_(0.),
          solver_control_(),
          compute_results_(true),
          init_vm_pu_(1.04),
          sn_mva_(1.0),
          max_nb_bus_per_sub_(2){
            _solver.change_solver(SolverType::SparseLU);
            _dc_solver.change_solver(SolverType::DC);
            _solver.set_gridmodel(this);
            _dc_solver.set_gridmodel(this);
            solver_control_.tell_all_changed();
        }
        GridModel(const GridModel & other);
        GridModel copy() const{
            GridModel res(*this);
            return res;
        }

        void set_ls_to_orig(const IntVect & ls_to_orig);  // set both _ls_to_orig and _orig_to_ls
        void set_orig_to_ls(const IntVect & orig_to_ls);  // set both _orig_to_ls and _ls_to_orig
        const IntVect & get_ls_to_orig(void) const {return _ls_to_orig;}
        const IntVect & get_orig_to_ls(void) const {return _orig_to_ls;}
        double timer_last_ac_pf() const {return timer_last_ac_pf_;}
        double timer_last_dc_pf() const {return timer_last_dc_pf_;}

        /**
         * @brief Return the total number of buses (both connected and disconnected)
         * 
         * For the number of connected buses see `nb_bus()`
         * 
         * @return Eigen::Index 
         */
        Eigen::Index total_bus() const {return bus_vn_kv_.size();}

        const std::vector<int> & id_me_to_ac_solver() const {return id_me_to_ac_solver_;}
        const std::vector<int> & id_ac_solver_to_me() const {return id_ac_solver_to_me_;}
        const std::vector<int> & id_me_to_dc_solver() const {return id_me_to_dc_solver_;}
        const std::vector<int> & id_dc_solver_to_me() const {return id_dc_solver_to_me_;}

        // retrieve the underlying data (raw class)
        const GeneratorContainer & get_generators_as_data() const {return generators_;}
        void turnedoff_no_pv(){generators_.turnedoff_no_pv(solver_control_);}  // turned off generators are not pv
        void turnedoff_pv(){generators_.turnedoff_pv(solver_control_);}  // turned off generators are pv
        bool get_turnedoff_gen_pv() {return generators_.get_turnedoff_gen_pv();}
        void update_slack_weights(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > could_be_slack){
            generators_.update_slack_weights(could_be_slack, solver_control_);
        }

        const SGenContainer & get_static_generators_as_data() const {return sgens_;}
        const LoadContainer & get_loads_as_data() const {return loads_;}
        const LineContainer & get_powerlines_as_data() const {return powerlines_;}
        const TrafoContainer & get_trafos_as_data() const {return trafos_;}
        const DCLineContainer & get_dclines_as_data() const {return dc_lines_;}
        Eigen::Ref<const RealVect> get_bus_vn_kv() const {return bus_vn_kv_;}
        std::tuple<int, int> assign_slack_to_most_connected();
        void consider_only_main_component();

        // solver "control"
        void change_solver(const SolverType & type){
            solver_control_.tell_all_changed();
            if(_solver.is_dc(type)) _dc_solver.change_solver(type);
            else _solver.change_solver(type);
        }
        std::vector<SolverType> available_solvers() {return _solver.available_solvers(); }
        SolverType get_solver_type() {return _solver.get_type(); }
        SolverType get_dc_solver_type() {return _dc_solver.get_type(); }
        const ChooseSolver & get_solver() const {return _solver;}
        const ChooseSolver & get_dc_solver() const {return _dc_solver;}

        // do i compute the results (in terms of P,Q,V or loads, generators and flows on lines
        void deactivate_result_computation(){compute_results_=false;}
        void reactivate_result_computation(){compute_results_=true;}

        // All methods to init this data model, all need to be pair unit when applicable
        void init_bus(const RealVect & bus_vn_kv, int nb_line, int nb_trafo);
        void set_init_vm_pu(real_type init_vm_pu) {init_vm_pu_ = init_vm_pu; }
        real_type get_init_vm_pu() {return init_vm_pu_;}
        void set_sn_mva(real_type sn_mva) {sn_mva_ = sn_mva; }
        real_type get_sn_mva() {return sn_mva_;}

        void init_powerlines(const RealVect & branch_r,
                             const RealVect & branch_x,
                             const CplxVect & branch_h,
                             const Eigen::VectorXi & branch_from_id,
                             const Eigen::VectorXi & branch_to_id
                             ){
            powerlines_.init(branch_r, branch_x, branch_h, branch_from_id, branch_to_id);
        }
        void init_powerlines_full(const RealVect & branch_r,
                                  const RealVect & branch_x,
                                  const CplxVect & branch_h_or,
                                  const CplxVect & branch_h_ex,
                                  const Eigen::VectorXi & branch_from_id,
                                  const Eigen::VectorXi & branch_to_id
                             ){
            powerlines_.init(branch_r, branch_x, branch_h_or,
                             branch_h_ex, branch_from_id, 
                             branch_to_id);
        }

        void init_shunt(const RealVect & shunt_p_mw,
                        const RealVect & shunt_q_mvar,
                        const Eigen::VectorXi & shunt_bus_id){
            shunts_.init(shunt_p_mw, shunt_q_mvar, shunt_bus_id);
        }
        void init_trafo(const RealVect & trafo_r,
                        const RealVect & trafo_x,
                        const CplxVect & trafo_b,
                        const RealVect & trafo_tap_step_pct,
                        const RealVect & trafo_tap_pos,
                        const RealVect & trafo_shift_degree,
                        const std::vector<bool> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::VectorXi & trafo_hv_id,
                        const Eigen::VectorXi & trafo_lv_id
                        ){
            trafos_.init(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_shift_degree,
                         trafo_tap_hv, trafo_hv_id, trafo_lv_id);
        }
        void init_generators(const RealVect & generators_p,
                             const RealVect & generators_v,
                             const RealVect & generators_min_q,
                             const RealVect & generators_max_q,
                             const Eigen::VectorXi & generators_bus_id){
            generators_.init(generators_p, generators_v, generators_min_q, generators_max_q, generators_bus_id);
        }
        void init_generators_full(const RealVect & generators_p,
                                  const RealVect & generators_v,
                                  const RealVect & generators_q,
                                  const std::vector<bool> & voltage_regulator_on,
                                  const RealVect & generators_min_q,
                                  const RealVect & generators_max_q,
                                  const Eigen::VectorXi & generators_bus_id){
            generators_.init_full(generators_p, generators_v, generators_q, voltage_regulator_on,
                                  generators_min_q, generators_max_q, generators_bus_id);
        }
        void init_loads(const RealVect & loads_p,
                        const RealVect & loads_q,
                        const Eigen::VectorXi & loads_bus_id){
            loads_.init(loads_p, loads_q, loads_bus_id);
        }
        void init_sgens(const RealVect & sgen_p,
                        const RealVect & sgen_q,
                        const RealVect & sgen_pmin,
                        const RealVect & sgen_pmax,
                        const RealVect & sgen_qmin,
                        const RealVect & sgen_qmax,
                        const Eigen::VectorXi & sgen_bus_id){
            sgens_.init(sgen_p, sgen_q, sgen_pmin, sgen_pmax, sgen_qmin, sgen_qmax, sgen_bus_id);
        }
        void init_storages(const RealVect & storages_p,
                           const RealVect & storages_q,
                           const Eigen::VectorXi & storages_bus_id){
            storages_.init(storages_p, storages_q, storages_bus_id);
        }
        void init_dclines(const Eigen::VectorXi & branch_from_id,
                          const Eigen::VectorXi & branch_to_id,
                          const RealVect & p_mw,
                          const RealVect & loss_percent,
                          const RealVect & loss_mw,
                          const RealVect & vm_or_pu,
                          const RealVect & vm_ex_pu,
                          const RealVect & min_q_or,
                          const RealVect & max_q_or,
                          const RealVect & min_q_ex,
                          const RealVect & max_q_ex){
            dc_lines_.init(branch_from_id, branch_to_id, p_mw,
                           loss_percent, loss_mw, vm_or_pu, vm_ex_pu,
                           min_q_or, max_q_or, min_q_ex, max_q_ex);
        }
        void init_bus_status(){
            const int nb_bus = static_cast<int>(bus_status_.size());
            for(int i = 0; i < nb_bus; ++i) bus_status_[i] = false;
            powerlines_.reconnect_connected_buses(bus_status_);
            shunts_.reconnect_connected_buses(bus_status_);
            trafos_.reconnect_connected_buses(bus_status_);
            generators_.reconnect_connected_buses(bus_status_);
            loads_.reconnect_connected_buses(bus_status_);
            sgens_.reconnect_connected_buses(bus_status_);
            storages_.reconnect_connected_buses(bus_status_);
            dc_lines_.reconnect_connected_buses(bus_status_);
        }

        void add_gen_slackbus(int gen_id, real_type weight);
        void remove_gen_slackbus(int gen_id);

        //pickle
        GridModel::StateRes get_state() const ;
        void set_state(GridModel::StateRes & my_state) ;
        template<class T>
        void check_size(const T& my_state)
        {
            // currently un used
            unsigned int size_th = 6;
            if (my_state.size() != size_th)
            {
                // std::cout << "LightSim::GridModel state size " << my_state.size() << " instead of "<< size_th << std::endl;
                // TODO more explicit error message
                throw std::runtime_error("Invalid state when loading LightSim::GridModel");
            }
        }

        //powerflows
        // control the need to refactorize the topology
        void unset_changes(){
            solver_control_.tell_none_changed();
        }  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_recompute_ybus(){solver_control_.tell_recompute_ybus();}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_recompute_sbus(){solver_control_.tell_recompute_sbus();}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_solver_need_reset(){solver_control_.tell_solver_need_reset();}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_ybus_change_sparsity_pattern(){solver_control_.tell_ybus_change_sparsity_pattern();}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        const SolverControl & get_solver_control() const {return solver_control_;}

        // dc powerflow
        CplxVect dc_pf(const CplxVect & Vinit,
                       int max_iter,  // not used for DC
                       real_type tol  // not used for DC
                       );

        /**
         * @brief Retrieve the PTDF matrice, with bus labeled in the "gridmodel" format.
         * Deactivated buses are represented by a column of 0.
         * 
         * It has the size (nb_line + nb_trafo, total_bus)
         * 
         * @return RealMat 
        */
        RealMat get_ptdf();

        /**
         * @brief Retrieve the LODF matrice.
         * 
         * 
         * It has the size (nb_line + nb_trafo, nb_line + nb_trafo)
         * 
         * @return RealMat 
        */
        RealMat get_lodf();

        /**
         * @brief Retrieve the PTDF matrice, with bus labeled in the "solver" format.
         * 
         * Deactivated buses are represented by a column of 0.
         * 
         * It has the size (nb_line + nb_trafo, nb_bus)
         * NB nb_bus is the number of activated buses !
         * 
         * @return RealMat 
        */
        RealMat get_ptdf_solver();
        
        Eigen::SparseMatrix<real_type> get_Bf_solver();
        Eigen::SparseMatrix<real_type> get_Bf();

        // ac powerflow
        CplxVect ac_pf(const CplxVect & Vinit,
                       int max_iter,
                       real_type tol);

        // check the kirchoff law
        CplxVect check_solution(const CplxVect & V, bool check_q_limits);

        // deactivate a bus. Be careful, if a bus is deactivated, but an element is
        //still connected to it, it will throw an exception
        void deactivate_bus(int bus_id) {
            if(bus_status_[bus_id]){
                // bus was connected, dim of matrix change
                solver_control_.need_reset_solver();
                solver_control_.need_recompute_sbus();
                solver_control_.need_recompute_ybus();
                solver_control_.ybus_change_sparsity_pattern();
                _deactivate(bus_id, bus_status_);
            }
        }
        // if a bus is connected, but isolated, it will make the powerflow diverge
        void reactivate_bus(int bus_id) {
            if(!bus_status_[bus_id]){
                // bus was not connected, dim of matrix change
                solver_control_.need_reset_solver();
                solver_control_.need_recompute_sbus();
                solver_control_.need_recompute_ybus();
                solver_control_.ybus_change_sparsity_pattern();
                _reactivate(bus_id, bus_status_); 
            }
        }
        /**
         * @brief Return the total number of connected buses !
         * 
         * For the total number of buses, see `total_bus()`
         * 
         * @return int 
         */
        int nb_bus() const;
        Eigen::Index nb_powerline() const {return powerlines_.nb();}
        Eigen::Index nb_trafo() const {return trafos_.nb();}

        // read only data accessor
        const LineContainer & get_lines() const {return powerlines_;}
        const DCLineContainer & get_dclines() const {return dc_lines_;}
        const TrafoContainer & get_trafos() const {return trafos_;}
        const GeneratorContainer & get_generators() const {return generators_;}
        const LoadContainer & get_loads() const {return loads_;}
        const LoadContainer & get_storages() const {return storages_;}
        const SGenContainer & get_static_generators() const {return sgens_;}
        const ShuntContainer & get_shunts() const {return shunts_;}
        const std::vector<bool> & get_bus_status() const {return bus_status_;}
        
        void set_line_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, powerlines_.nb(), "set_line_names");
            powerlines_.set_names(names);
        }
        void set_dcline_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, dc_lines_.nb(), "set_dcline_names");
            dc_lines_.set_names(names);
        }
        void set_trafo_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, trafos_.nb(), "set_trafo_names");
            trafos_.set_names(names);
        }
        void set_gen_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, generators_.nb(), "set_gen_names");
            generators_.set_names(names);
        }
        void set_load_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, loads_.nb(), "set_load_names");
            loads_.set_names(names);
        }
        void set_storage_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, storages_.nb(), "set_storage_names");
            storages_.set_names(names);
        }
        void set_sgen_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, sgens_.nb(), "set_sgen_names");
            sgens_.set_names(names);
        }
        void set_shunt_names(const std::vector<std::string> & names){
            GenericContainer::check_size(names, shunts_.nb(), "set_shunt_names");
            shunts_.set_names(names);
        }

        //deactivate a powerline (disconnect it)
        void deactivate_powerline(int powerline_id) {powerlines_.deactivate(powerline_id, solver_control_); }
        void reactivate_powerline(int powerline_id) {powerlines_.reactivate(powerline_id, solver_control_); }
        void change_bus_powerline_or(int powerline_id, int new_bus_id) {powerlines_.change_bus_or(powerline_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_bus_powerline_ex(int powerline_id, int new_bus_id) {powerlines_.change_bus_ex(powerline_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        int get_bus_powerline_or(int powerline_id) {return powerlines_.get_bus_or(powerline_id);}
        int get_bus_powerline_ex(int powerline_id) {return powerlines_.get_bus_ex(powerline_id);}

        //deactivate trafo
        void deactivate_trafo(int trafo_id) {trafos_.deactivate(trafo_id, solver_control_); }
        void reactivate_trafo(int trafo_id) {trafos_.reactivate(trafo_id, solver_control_); }
        void change_bus_trafo_hv(int trafo_id, int new_bus_id) {trafos_.change_bus_hv(trafo_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_bus_trafo_lv(int trafo_id, int new_bus_id) {trafos_.change_bus_lv(trafo_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        int get_bus_trafo_hv(int trafo_id) {return trafos_.get_bus_hv(trafo_id);}
        int get_bus_trafo_lv(int trafo_id) {return trafos_.get_bus_lv(trafo_id);}

        //load
        void deactivate_load(int load_id) {loads_.deactivate(load_id, solver_control_); }
        void reactivate_load(int load_id) {loads_.reactivate(load_id, solver_control_); }
        void change_bus_load(int load_id, int new_bus_id) {loads_.change_bus(load_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_p_load(int load_id, real_type new_p) {loads_.change_p(load_id, new_p, solver_control_); }
        void change_q_load(int load_id, real_type new_q) {loads_.change_q(load_id, new_q, solver_control_); }
        int get_bus_load(int load_id) {return loads_.get_bus(load_id);}

        //generator
        void deactivate_gen(int gen_id) {generators_.deactivate(gen_id, solver_control_); }
        void reactivate_gen(int gen_id) {generators_.reactivate(gen_id, solver_control_); }
        void change_bus_gen(int gen_id, int new_bus_id) {generators_.change_bus(gen_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_p_gen(int gen_id, real_type new_p) {generators_.change_p(gen_id, new_p, solver_control_); }
        void change_v_gen(int gen_id, real_type new_v_pu) {generators_.change_v(gen_id, new_v_pu, solver_control_); }
        int get_bus_gen(int gen_id) {return generators_.get_bus(gen_id);}

        //shunt
        void deactivate_shunt(int shunt_id) {shunts_.deactivate(shunt_id, solver_control_); }
        void reactivate_shunt(int shunt_id) {shunts_.reactivate(shunt_id, solver_control_); }
        void change_bus_shunt(int shunt_id, int new_bus_id) {shunts_.change_bus(shunt_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size()));  }
        void change_p_shunt(int shunt_id, real_type new_p) {shunts_.change_p(shunt_id, new_p, solver_control_); }
        void change_q_shunt(int shunt_id, real_type new_q) {shunts_.change_q(shunt_id, new_q, solver_control_); }
        int get_bus_shunt(int shunt_id) {return shunts_.get_bus(shunt_id);}

        //static gen
        void deactivate_sgen(int sgen_id) {sgens_.deactivate(sgen_id, solver_control_); }
        void reactivate_sgen(int sgen_id) {sgens_.reactivate(sgen_id, solver_control_); }
        void change_bus_sgen(int sgen_id, int new_bus_id) {sgens_.change_bus(sgen_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_p_sgen(int sgen_id, real_type new_p) {sgens_.change_p(sgen_id, new_p, solver_control_); }
        void change_q_sgen(int sgen_id, real_type new_q) {sgens_.change_q(sgen_id, new_q, solver_control_); }
        int get_bus_sgen(int sgen_id) {return sgens_.get_bus(sgen_id);}

        //storage units
        void deactivate_storage(int storage_id) {storages_.deactivate(storage_id, solver_control_); }
        void reactivate_storage(int storage_id) {storages_.reactivate(storage_id, solver_control_); }
        void change_bus_storage(int storage_id, int new_bus_id) {storages_.change_bus(storage_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_p_storage(int storage_id, real_type new_p) {
//            if(new_p == 0.)
//            {
//                std::cout << " i deactivated storage with ID " << storage_id << std::endl;
//                storages_.change_p(storage_id, new_p, need_reset_);
//                deactivate_storage(storage_id);  // requirement from grid2op, might be discussed
//            }else{
//                reactivate_storage(storage_id);  // requirement from grid2op, might be discussed
//                storages_.change_p(storage_id, new_p, need_reset_);
//            }
               storages_.change_p(storage_id, new_p, solver_control_);
            }
        void change_q_storage(int storage_id, real_type new_q) {storages_.change_q(storage_id, new_q, solver_control_); }
        int get_bus_storage(int storage_id) {return storages_.get_bus(storage_id);}

        //deactivate a powerline (disconnect it)
        void deactivate_dcline(int dcline_id) {dc_lines_.deactivate(dcline_id, solver_control_); }
        void reactivate_dcline(int dcline_id) {dc_lines_.reactivate(dcline_id, solver_control_); }
        void change_p_dcline(int dcline_id, real_type new_p) {dc_lines_.change_p(dcline_id, new_p, solver_control_); }
        void change_v_or_dcline(int dcline_id, real_type new_v_pu) {dc_lines_.change_v_or(dcline_id, new_v_pu, solver_control_); }
        void change_v_ex_dcline(int dcline_id, real_type new_v_pu) {dc_lines_.change_v_ex(dcline_id, new_v_pu, solver_control_); }
        void change_bus_dcline_or(int dcline_id, int new_bus_id) {dc_lines_.change_bus_or(dcline_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        void change_bus_dcline_ex(int dcline_id, int new_bus_id) {dc_lines_.change_bus_ex(dcline_id, new_bus_id, solver_control_, static_cast<int>(bus_vn_kv_.size())); }
        int get_bus_dcline_or(int dcline_id) {return dc_lines_.get_bus_or(dcline_id);}
        int get_bus_dcline_ex(int dcline_id) {return dc_lines_.get_bus_ex(dcline_id);}

        // All results access
        tuple3d get_loads_res() const {return loads_.get_res();}
        const std::vector<bool>& get_loads_status() const { return loads_.get_status();}
        tuple3d get_shunts_res() const {return shunts_.get_res();}
        const std::vector<bool>& get_shunts_status() const { return shunts_.get_status();}
        tuple3d get_gen_res() const {return generators_.get_res();}
        const std::vector<bool>& get_gen_status() const { return generators_.get_status();}
        tuple4d get_lineor_res() const {return powerlines_.get_lineor_res();}
        tuple4d get_lineex_res() const {return powerlines_.get_lineex_res();}
        const std::vector<bool>& get_lines_status() const { return powerlines_.get_status();}
        tuple4d get_trafohv_res() const {return trafos_.get_res_hv();}
        tuple4d get_trafolv_res() const {return trafos_.get_res_lv();}
        const std::vector<bool>& get_trafo_status() const { return trafos_.get_status();}
        tuple3d get_storages_res() const {return storages_.get_res();}
        const std::vector<bool>& get_storages_status() const { return storages_.get_status();}
        tuple3d get_sgens_res() const {return sgens_.get_res();}
        const std::vector<bool>& get_sgens_status() const { return sgens_.get_status();}
        tuple3d get_dclineor_res() const {return dc_lines_.get_or_res();}
        tuple3d get_dclineex_res() const {return dc_lines_.get_ex_res();}
        const std::vector<bool>& get_dclines_status() const { return dc_lines_.get_status();}

        Eigen::Ref<const RealVect> get_gen_theta() const  {return generators_.get_theta();}
        Eigen::Ref<const RealVect> get_load_theta() const  {return loads_.get_theta();}
        Eigen::Ref<const RealVect> get_shunt_theta() const  {return shunts_.get_theta();}
        Eigen::Ref<const RealVect> get_storage_theta() const  {return storages_.get_theta();}
        Eigen::Ref<const RealVect> get_lineor_theta() const {return powerlines_.get_theta_or();}
        Eigen::Ref<const RealVect> get_lineex_theta() const {return powerlines_.get_theta_ex();}
        Eigen::Ref<const RealVect> get_trafohv_theta() const {return trafos_.get_theta_hv();}
        Eigen::Ref<const RealVect> get_trafolv_theta() const {return trafos_.get_theta_lv();}
        Eigen::Ref<const RealVect> get_dclineor_theta() const {return dc_lines_.get_theta_or();}
        Eigen::Ref<const RealVect> get_dclineex_theta() const {return dc_lines_.get_theta_ex();}

        Eigen::Ref<const IntVect> get_all_shunt_buses() const {return shunts_.get_buses();}

        // complete results (with theta)
        tuple4d get_loads_res_full() const {return loads_.get_res_full();}
        tuple4d get_shunts_res_full() const {return shunts_.get_res_full();}
        tuple4d get_gen_res_full() const {return generators_.get_res_full();}
        tuple5d get_lineor_res_full() const {return powerlines_.get_res_or_full();}
        tuple5d get_lineex_res_full() const {return powerlines_.get_res_ex_full();}
        tuple5d get_trafohv_res_full() const {return trafos_.get_res_hv_full();}
        tuple5d get_trafolv_res_full() const {return trafos_.get_res_lv_full();}
        tuple4d get_storages_res_full() const {return storages_.get_res_full();}
        tuple4d get_sgens_res_full() const {return sgens_.get_res_full();}
        tuple4d get_dclineor_res_full() const {return dc_lines_.get_res_or_full();}
        tuple4d get_dclineex_res_full() const {return dc_lines_.get_res_ex_full();}

        /**
         * @brief Get the Ybus solver object (AC)
         * 
         * This function allows to retrieve the Ybus passed to the AC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_Ybus` before that)
         * 
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        Eigen::SparseMatrix<cplx_type> get_Ybus_solver(){
            return Ybus_ac_;  // This is copied to python
        }

        /**
         * @brief Get the Ybus solver object (DC)
         * 
         * This function allows to retrieve the Ybus passed to the DC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_dcYbus` before that)
         * 
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        Eigen::SparseMatrix<cplx_type> get_dcYbus_solver(){
            return Ybus_dc_;  // This is copied to python
        }

        /**
         * @brief Get the Sbus solver object (AC)
         * 
         * This function allows to retrieve the Sbus passed to the AC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_Sbus` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        Eigen::Ref<const CplxVect> get_Sbus_solver() const{
            return acSbus_;
        }

        /**
         * @brief Get the Sbus solver object (DC)
         * 
         * This function allows to retrieve the Sbus passed to the DC solver.
         * 
         * It has the "solver bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (nb_bus, nb_bus) (number of active buses)
         * 
         * (new in lightsim2grid 0.9.0, used to be called `get_dcSbus` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        Eigen::Ref<const CplxVect> get_dcSbus_solver() const{
            return dcSbus_;
        }

        /**
         * @brief Get the Ybus gridmodel object (AC powerflow)
         * 
         * This function allows to retrieve the Ybus as seen by the gridmodel.
         * 
         * It may contain empty rows / columns for disconnected buses.
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of active buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Ybus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        const Eigen::SparseMatrix<cplx_type> get_Ybus() const {
            return _relabel_matrix(Ybus_ac_, id_ac_solver_to_me_);
        }

        /**
         * @brief Get the Ybus gridmodel object (DC powerflow)
         * 
         * This function allows to retrieve the Ybus (DC) as seen by the gridmodel.
         * 
         * It may contain empty rows / columns for disconnected buses.
         * 
         * It has the "grimodel bus id" labelling, which is different from the 
         * "gridmodel" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of active buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_dcYbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        const Eigen::SparseMatrix<cplx_type> get_dcYbus() const {
            return _relabel_matrix(Ybus_dc_, id_dc_solver_to_me_);
        }

        /**
         * @brief Get the Sbus solver object (AC)
         * 
         * This function allows to retrieve the Sbus as represented by the "gridmodel"
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "solver" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of total buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Sbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        const CplxVect get_Sbus() const {
            return _relabel_vector(acSbus_, id_ac_solver_to_me_);
        }

        /**
         * @brief Get the Sbus solver object (DC)
         * 
         * This function allows to retrieve the Sbus as represented by the "gridmodel"
         * 
         * It has the "gridmodel bus id" labelling, which is different from the 
         * "solver" bus labelling.
         * 
         * It has the size (total_bus, total_bus) (number of total buses)
         * 
         * (change in lightsim2grid 0.9.0, same as old function is now `get_Sbus_solver` before that)
         * 
         * @return Eigen::Ref<const CplxVect>
         */
        const CplxVect get_dcSbus() const {
            return _relabel_vector(dcSbus_, id_dc_solver_to_me_);
        }

        /**
         * @brief Get the vector (list) of pv buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        Eigen::Ref<const Eigen::VectorXi> get_pv_solver() const{
            return bus_pv_;
        }

        /**
         * @brief Get the vector (list) of pv buses, gridmodel labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return const Eigen::VectorXi
         */
        const Eigen::VectorXi get_pv() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector2(get_pv_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector2(get_pv_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_pv: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the vector (list) of pq buses, solver labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        Eigen::Ref<const Eigen::VectorXi> get_pq_solver() const{
            return bus_pq_;
        }

        /**
         * @brief Get the vector (list) of pq buses, grimodel labelling
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return const Eigen::VectorXi
         */
        const Eigen::VectorXi get_pq() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector2(get_pq_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector2(get_pq_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_pq: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (AC), solver labelling
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        Eigen::Ref<const Eigen::VectorXi> get_slack_ids_solver() const{
            return slack_bus_id_ac_solver_;
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (AC), gridmodel labelling
         * 
         * @return const Eigen::VectorXi 
         */
        const Eigen::VectorXi get_slack_ids() const{
            return _relabel_vector2(slack_bus_id_ac_solver_, id_ac_solver_to_me_);
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (DC), solver labelling
         * 
         * @return Eigen::Ref<const Eigen::VectorXi> 
         */
        Eigen::Ref<const Eigen::VectorXi> get_slack_ids_dc_solver() const{
            return slack_bus_id_dc_solver_;
        }

        /**
         * @brief Get the ids of the buses that participate to the slack (DC), gridmodel labelling
         * 
         * @return const Eigen::VectorXi 
         */
        const Eigen::VectorXi get_slack_ids_dc() const{
            return _relabel_vector2(slack_bus_id_dc_solver_, id_dc_solver_to_me_);
        }

        /**
         * @brief Get the slack weights for each buses (solver labelling)
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        Eigen::Ref<const RealVect> get_slack_weights_solver() const{
            return slack_weights_;
        }

        /**
         * @brief Get the slack weights for each buses (gridmodel labelling)
         * 
         * valid for AC modeling only, TODO in DC
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        const RealVect get_slack_weights() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_slack_weights_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_slack_weights_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_slack_weights: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (complex) voltage angles for each buses (solver labelling)
         * 
         * @return Eigen::Ref<const CplxVect> 
         */
        Eigen::Ref<const CplxVect> get_V_solver() const{
            return _solver.get_V();
        }

        /**
         * @brief Get the (complex) voltage angles for each buses (grimodel labelling)
         * 
         * @return CplxVect
         */
        const CplxVect get_V() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_V_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_V_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_V: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (real) voltage angle for each buses of the grid (solver labelling)
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        Eigen::Ref<const RealVect> get_Va_solver() const{
            return _solver.get_Va();
        }

        /**
         * @brief Get the (real) voltage angle for each buses of the grid (grimodel labelling)
         * 
         * @return const RealVect
         */
        const RealVect get_Va() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_Va_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_Va_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_Va: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        /**
         * @brief Get the (real) voltage magnitude for each buses of the grid (solver labelling)
         * 
         * @return Eigen::Ref<const RealVect> 
         */
        Eigen::Ref<const RealVect> get_Vm_solver() const{
            return _solver.get_Vm();
        }

        /**
         * @brief Get the (real) voltage magnitude for each buses of the grid (grimodel labelling)
         * 
         * @return const RealVect
         */
        const RealVect get_Vm() const{
            if(id_ac_solver_to_me_.size() > 0) return _relabel_vector(get_Vm_solver(), id_ac_solver_to_me_);
            if(id_dc_solver_to_me_.size() > 0) return _relabel_vector(get_Vm_solver(), id_dc_solver_to_me_);
            throw std::runtime_error("GridModel::get_Vm: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
        }

        Eigen::Ref<const Eigen::SparseMatrix<real_type> > get_J_solver() const{
            return _solver.get_J();
        }

        /**
         * @brief Returns the last computed jacobian matrix (solver labelling)
         * 
         * @return Eigen::SparseMatrix<real_type> 
         */
        Eigen::SparseMatrix<real_type> get_J_python_solver() const{
            return _solver.get_J_python();  // This is copied to python
        }
        
        real_type get_computation_time() const{ return _solver.get_computation_time();}
        real_type get_dc_computation_time() const{ return _dc_solver.get_computation_time();}

    // private:
    //     using GenericContainer::update_bus_status;  // to silence clang warnings (overload-virtual)
    // public:
    //     void update_bus_status(int nb_bus_before,
    //                            Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> > active_bus);
        // part dedicated to grid2op backend, optimized for grid2op data representation (for speed)
        // this is not recommended to use it outside of its intended usage within grid2op !
        void update_gens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                           Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_sgens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                           Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_gens_v(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                           Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_loads_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_loads_q(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_topo(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                         Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_storages_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                               Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);

        void set_load_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > load_pos_topo_vect)
        {
            load_pos_topo_vect_.array() = load_pos_topo_vect;
        }
        void set_gen_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > gen_pos_topo_vect)
        {
            gen_pos_topo_vect_.array() = gen_pos_topo_vect;
        }
        void set_line_or_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > line_or_pos_topo_vect)
        {
            line_or_pos_topo_vect_.array() = line_or_pos_topo_vect;
        }
        void set_line_ex_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > line_ex_pos_topo_vect)
        {
            line_ex_pos_topo_vect_.array() = line_ex_pos_topo_vect;
        }
        void set_trafo_hv_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > trafo_hv_pos_topo_vect)
        {
            trafo_hv_pos_topo_vect_.array() = trafo_hv_pos_topo_vect;
        }
        void set_trafo_lv_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > trafo_lv_pos_topo_vect)
        {
            trafo_lv_pos_topo_vect_.array() = trafo_lv_pos_topo_vect;
        }
        void set_storage_pos_topo_vect(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > storage_pos_topo_vect)
        {
            storage_pos_topo_vect_.array() = storage_pos_topo_vect;
        }

        void set_load_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > load_to_subid)
        {
            load_to_subid_.array() = load_to_subid;
        }
        void set_gen_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > gen_to_subid)
        {
            gen_to_subid_.array() = gen_to_subid;
        }
        void set_line_or_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > line_or_to_subid)
        {
            line_or_to_subid_.array() = line_or_to_subid;
        }
        void set_line_ex_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > line_ex_to_subid)
        {
            line_ex_to_subid_.array() = line_ex_to_subid;
        }
        void set_trafo_hv_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > trafo_hv_to_subid)
        {
            trafo_hv_to_subid_.array() = trafo_hv_to_subid;
        }
        void set_trafo_lv_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > trafo_lv_to_subid)
        {
            trafo_lv_to_subid_.array() = trafo_lv_to_subid;
        }
        void set_storage_to_subid(Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > storage_to_subid)
        {
            storage_to_subid_.array() = storage_to_subid;
        }
        void set_n_sub(int n_sub)
        {
            n_sub_ = n_sub;
        }
        void set_max_nb_bus_per_sub(int max_nb_bus_per_sub)
        {
            if(bus_vn_kv_.size() != n_sub_ * max_nb_bus_per_sub){
                std::ostringstream exc_;
                exc_ << "GridModel::set_max_nb_bus_per_sub: ";
                exc_ << "your model counts ";
                exc_ << bus_vn_kv_.size()  << " buses according to `bus_vn_kv_` but ";
                exc_ << n_sub_ * max_nb_bus_per_sub_ << " according to n_sub_ * max_nb_bus_per_sub_.";
                exc_ << "Both should match: either reinit it with another call to `init_bus` or set properly the number of ";
                exc_ << "substations / buses per substations with `set_n_sub` / `set_max_nb_bus_per_sub`";
                throw std::runtime_error(exc_.str());
            }
            max_nb_bus_per_sub_ = max_nb_bus_per_sub;
        }
        int get_max_nb_bus_per_sub() const { return max_nb_bus_per_sub_;}
        
        void fillSbus_other(CplxVect & res, bool ac, const std::vector<int>& id_me_to_solver){
            fillSbus_me(res, ac, id_me_to_solver);
        }

        //for FDPF
    private:
        using GenericContainer::fillBp_Bpp;  // silence clang warning overload-virtual
    public:
        void fillBp_Bpp(Eigen::SparseMatrix<real_type> & Bp, 
                        Eigen::SparseMatrix<real_type> & Bpp, 
                        FDPFMethod xb_or_bx) const;

    private:
        using GenericContainer::fillBf_for_PTDF;  // silence clang warning overload-virtual
    public:
        void fillBf_for_PTDF(Eigen::SparseMatrix<real_type> & Bf, bool transpose=false) const;

        Eigen::SparseMatrix<real_type> debug_get_Bp_python(FDPFMethod xb_or_bx){
            Eigen::SparseMatrix<real_type> Bp;
            Eigen::SparseMatrix<real_type> Bpp;
            fillBp_Bpp(Bp, Bpp, xb_or_bx);
            return Bp;
        }
        Eigen::SparseMatrix<real_type> debug_get_Bpp_python(FDPFMethod xb_or_bx){
            Eigen::SparseMatrix<real_type> Bp;
            Eigen::SparseMatrix<real_type> Bpp;
            fillBp_Bpp(Bp, Bpp, xb_or_bx);
            return Bpp;
        }

    protected:
    // add method to change topology, change ratio of transformers, change

        // compute admittance matrix
        // dc powerflow
        // void init_dcY(Eigen::SparseMatrix<real_type> & dcYbus);

        // ac powerflows
        /**
        computes Ybus_ and Sbus_. It has different flags to have more control on the purpose for this "computation"
        is_ac indicates if you want to perform and AC powerflow or a DC powerflow and reset_solver indicates
        if you will perform a powerflow after it or not. (usually put ``true`` here).
        **/
        CplxVect pre_process_solver(const CplxVect & Vinit,
                                    CplxVect & Sbus,
                                    Eigen::SparseMatrix<cplx_type> & Ybus,
                                    std::vector<int> & id_me_to_solver,
                                    std::vector<int> & id_solver_to_me,
                                    Eigen::VectorXi & slack_bus_id_me,
                                    Eigen::VectorXi & slack_bus_id_solver,
                                    bool is_ac,
                                    const SolverControl & solver_control);

        // init the Ybus matrix (its size, it is filled up elsewhere) and also the 
        // converter from "my bus id" to the "solver bus id" (id_me_to_solver and id_solver_to_me)
        void init_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus,
                       std::vector<int> & id_me_to_solver,
                       std::vector<int>& id_solver_to_me);

        // converts the slack_bus_id from gridmodel ordering into solver ordering
        void init_slack_bus(const CplxVect & Sbus,
                            const std::vector<int> & id_me_to_solver,
                            const std::vector<int>& id_solver_to_me,
                            const Eigen::VectorXi & slack_bus_id_me,
                            Eigen::VectorXi & slack_bus_id_solver
                        );

        /**
         * @brief Build the result matrix (eg Ybus) (labelled using the gridmodel) 
         * from the input same matrix (eg Ybus) but labelled with the solver convention
         * 
         * @param Ybus : solver labelling
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @param relabel_row : whether to relabel also the row id
         * @return Eigen::SparseMatrix<cplx_type> 
         */
        template<typename T>    
        Eigen::SparseMatrix<T> _relabel_matrix(const Eigen::SparseMatrix<T> & Ybus,
                                               const std::vector<int> & id_solver_to_me,
                                               bool relabel_row=true) const {
            // TODO optim : if relabel_row is false, then we can just copy
            // paste the columns easily in the target matrix, which should be
            // way faster than this function.
            typedef typename Eigen::SparseMatrix<T>::StorageIndex index_type;
            if(id_solver_to_me.size() == 0) throw std::runtime_error("GridModel::_relabel_matrix: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Ybus.cols() != nb_bus()) throw std::runtime_error("GridModel::_relabel_matrix: impossible to retrieve the `gridmodel`: the input matrix has not the right number of columns, (.., nb connected bus) expected");
            if(relabel_row & (Ybus.rows() != nb_bus())) throw std::runtime_error("GridModel::_relabel_matrix: impossible to retrieve the `gridmodel`: the input matrix has not the right number of columnd (nb connected bus, ...) expected");
            Eigen::SparseMatrix<T> res(relabel_row ? total_bus() : Ybus.rows(), total_bus());
            res.reserve(Ybus.nonZeros());
            std::vector<Eigen::Triplet<T> > tripletList;
            tripletList.reserve(Ybus.nonZeros());
            const auto n_col = Ybus.cols();
            for (Eigen::Index col_=0; col_ < n_col; ++col_){
                for (typename Eigen::SparseMatrix<T>::InnerIterator it(Ybus, col_); it; ++it)
                {
                    if(relabel_row) tripletList.push_back({static_cast<index_type>(id_solver_to_me[it.row()]),
                                                           static_cast<index_type>(id_solver_to_me[it.col()]),
                                                           it.value()});
                    else tripletList.push_back({static_cast<index_type>(it.row()), 
                                                static_cast<index_type>(id_solver_to_me[it.col()]),
                                                it.value()});
                }
            }
            res.setFromTriplets(tripletList.begin(), tripletList.end());
            res.makeCompressed();
            return res;
        }

        /**
         * @brief Build the Sbus (or any other vector labelled using the gridmodel convention) 
         * from the same vector (input) that uses the solver convention.
         * 
         * TODO copy paste from below, find a better way !
         * 
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect 
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1> > & Sbus,
                                                            const std::vector<int> & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Sbus.size() != nb_bus()) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` input solver has not the right size, expected (nb connected bus, ).");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(total_bus());
            for(auto solver_id = 0; solver_id < Sbus.size(); ++solver_id){
                res[id_solver_to_me[solver_id]] = Sbus[solver_id];
            }
            return res;
        }

        /**
         * @brief Build the Sbus (or any other vector labelled using the gridmodel convention) 
         * from the same vector (input) that uses the solver convention.
         * 
         * TODO copy paste from above, find a better way !
         * 
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect 
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector(const Eigen::Matrix<T, Eigen::Dynamic, 1> & Sbus,
                                                            const std::vector<int> & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            if(Sbus.size() != nb_bus()) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` input solver has not the right size, expected (nb connected bus, ).");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(total_bus());
            for(auto solver_id = 0; solver_id < Sbus.size(); ++solver_id){
                res[id_solver_to_me[solver_id]] = Sbus[solver_id];
            }
            return res;
        }

        /**
         * @brief Build the pv; pq or slack ids (or any other vector labelled using the gridmodel convention) 
         * from the same vector (input) that uses the solver convention.
         * 
         * TODO copy paste from above, find a better way !
         * 
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect 
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector2(const Eigen::Matrix<T, Eigen::Dynamic, 1> & pv_pq_ref_bus,
                                                             const std::vector<int> & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(pv_pq_ref_bus.size());
            Eigen::Index pos_id = 0;
            for(auto el_id : pv_pq_ref_bus){
                res[pos_id] = id_solver_to_me[el_id];
                ++ pos_id;
            }
            return res;
        }

        /**
         * @brief Build the pv; pq or slack ids (or any other vector labelled using the gridmodel convention) 
         * from the same vector (input) that uses the solver convention.
         * 
         * TODO copy paste from above, find a better way !
         * 
         * @param Sbus : Sbus with the solver convention, the one used by the solver
         * @param id_solver_to_me : mapping to convert from the solver id to the gridmodel id
         * @return CplxVect 
         */
        template<class T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> _relabel_vector2(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 1> > & pv_pq_ref_bus,
                                                             const std::vector<int> & id_solver_to_me) const
        {
            if(id_solver_to_me.size() == 0) throw std::runtime_error("GridModel::_relabel_vector: impossible to retrieve the `gridmodel` bus label as it appears no powerflow has run.");
            Eigen::Matrix<T, Eigen::Dynamic, 1> res = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(pv_pq_ref_bus.size());
            Eigen::Index pos_id = 0;
            for(auto el_id : pv_pq_ref_bus){
                res[pos_id] = id_solver_to_me[el_id];
                ++ pos_id;
            }
            return res;
        }

    private:
        using GenericContainer::fillYbus;  // to silence the overload-virtual warning in clang
    protected:
        void fillYbus(Eigen::SparseMatrix<cplx_type> & res, bool ac, const std::vector<int>& id_me_to_solver);
        void fillSbus_me(CplxVect & res, bool ac, const std::vector<int>& id_me_to_solver);
        void fillpv_pq(const std::vector<int>& id_me_to_solver,
                       const std::vector<int>& id_solver_to_me,
                       const Eigen::VectorXi & slack_bus_id_solver,
                       const SolverControl & solver_control);

        // results
        /**process the results from the solver to this instance
        **/
        void process_results(bool conv, CplxVect & res, const CplxVect & Vinit, bool ac,
                             std::vector<int> & id_me_to_solver);

        /**
        Compute the results vector from the Va, Vm post powerflow
        **/
        void compute_results(bool ac);
        /**
        reset the results in case of divergence of the powerflow.
        **/
        void reset_results();

        /**
        reset the solver, and all its results
        **/
        void reset(bool reset_solver, bool reset_ac, bool reset_dc);

        /**
        optimization for grid2op
        **/
        template<class T>
        void update_continuous_values(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                                      Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
                                      T fun)
        {
            for(int el_id = 0; el_id < has_changed.rows(); ++el_id)
            {
                if(has_changed(el_id))
                {
                    (this->*fun)(el_id, static_cast<real_type>(new_values[el_id]));
                }
            }
        }
        template<class CReac, class CChange, class CDeact>
        void update_topo_generic(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
                                 Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
                                 const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> & vect_pos,
                                 const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> & vect_subid,
                                 CReac fun_react,
                                 CChange fun_change,
                                 CDeact fun_deact)
        {
            for(int el_id = 0; el_id < vect_pos.rows(); ++el_id)
            {

                int el_pos = vect_pos(el_id);
                if(! has_changed(el_pos)) continue;
                int new_bus = new_values(el_pos);
                if(new_bus > 0){
                    // new bus is a real bus, so i need to make sure to have it turned on, and then change the bus
                    int sub_id = vect_subid(el_id);
                    int new_bus_backend = sub_id + (new_bus - 1) * n_sub_;
                    bus_status_[new_bus_backend] = true;
                    (this->*fun_react)(el_id); // eg reactivate_load(load_id);
                    (this->*fun_change)(el_id, new_bus_backend); // eg change_bus_load(load_id, new_bus_backend);
                } else{
                    // new bus is negative, we deactivate it
                    (this->*fun_deact)(el_id);// eg deactivate_load(load_id);
                    // bus_status_ is set to "false" in GridModel.update_topo
                    // and a bus is activated if (and only if) one element is connected to it.
                    // I must not set `bus_status_[new_bus_backend] = false;` in this case !
                }
            }
        }

        CplxVect _get_results_back_to_orig_nodes(const CplxVect & res_tmp,
                                                 std::vector<int> & id_me_to_solver,
                                                 int size);

        void check_solution_q_values( CplxVect & res, bool check_q_limits) const;
        void check_solution_q_values_onegen(CplxVect & res, const GeneratorContainer::GenInfo& gen, bool check_q_limits) const;

    protected:
        // memory for the import
        IntVect _ls_to_orig;  // for converter from bus in lightsim2grid index to bus in original file format (*eg* pandapower or pypowsybl)
        IntVect _orig_to_ls;  // for converter from bus in lightsim2grid index to bus in original file format (*eg* pandapower or pypowsybl)

        // member of the grid
        double timer_last_ac_pf_;
        double timer_last_dc_pf_;

        // bool need_reset_solver_;  // some matrices change size, needs to be computed
        // bool need_recompute_sbus_;  // some coeff of sbus changed, need to recompute it
        // bool need_recompute_ybus_;  // some coeff of ybus changed, but not its sparsity pattern
        // bool ybus_change_sparsity_pattern_;  // sparsity pattern of ybus changed (and so are its coeff)
        SolverControl solver_control_;
        bool compute_results_;
        real_type init_vm_pu_;  // default vm initialization, mainly for dc powerflow
        real_type sn_mva_;

        // powersystem representation
        // 1. bus
        RealVect bus_vn_kv_;
        std::vector<bool> bus_status_;  // for each bus, gives its status. true if connected, false otherwise

        // always have the length of the number of buses,
        // id_me_to_model_[id_me] gives -1 if the bus "id_me" is deactivated, or "id_model" if it is activated.
        std::vector<int> id_me_to_ac_solver_;
        // convert the bus id from the model to the bus id of me.
        // it has a variable size, that depends on the number of connected bus. if "id_model" is an id of a bus
        // sent to the solver, then id_model_to_me_[id_model] is the bus id of this model of the grid.
        std::vector<int> id_ac_solver_to_me_;

        std::vector<int> id_me_to_dc_solver_;
        std::vector<int> id_dc_solver_to_me_;

        // 2. powerline
        LineContainer powerlines_;

        // 3. shunt
        ShuntContainer shunts_;

        // 4. transformers
        // have the r, x, h and ratio
        // ratio is computed from the tap, so maybe store tap num and tap_step_pct
        TrafoContainer trafos_;

        // 5. generators
        RealVect total_q_min_per_bus_;
        RealVect total_q_max_per_bus_;
        Eigen::VectorXi total_gen_per_bus_;
        GeneratorContainer generators_;

        // 6. loads
        LoadContainer loads_;

        // 7. static generators (P,Q generators)
        SGenContainer sgens_;

        // 8. storage units
        LoadContainer storages_;

        // 9. hvdc
        DCLineContainer dc_lines_;

        // 10. slack bus
        // std::vector<int> slack_bus_id_;
        Eigen::VectorXi slack_bus_id_ac_me_;  // slack bus id, gridmodel number
        Eigen::VectorXi slack_bus_id_ac_solver_;  // slack bus id, solver number
        Eigen::VectorXi slack_bus_id_dc_me_;
        Eigen::VectorXi slack_bus_id_dc_solver_;
        RealVect slack_weights_;

        // as matrix, for the solver
        Eigen::SparseMatrix<cplx_type> Ybus_ac_;
        Eigen::SparseMatrix<cplx_type> Ybus_dc_;
        CplxVect acSbus_;
        CplxVect dcSbus_;
        Eigen::VectorXi bus_pv_;  // id are the solver internal id and NOT the initial id
        Eigen::VectorXi bus_pq_;  // id are the solver internal id and NOT the initial id

        // TODO have version of the stuff above for the public api, indexed with "me" and not "solver"

        // to solve the newton raphson
        ChooseSolver _solver;
        ChooseSolver _dc_solver;

        // specific grid2op
        int n_sub_;
        int max_nb_bus_per_sub_;
        IntVectRowMaj load_pos_topo_vect_;
        IntVectRowMaj gen_pos_topo_vect_;
        IntVectRowMaj line_or_pos_topo_vect_;
        IntVectRowMaj line_ex_pos_topo_vect_;
        IntVectRowMaj trafo_hv_pos_topo_vect_;
        IntVectRowMaj trafo_lv_pos_topo_vect_;
        IntVectRowMaj storage_pos_topo_vect_;

        IntVectRowMaj load_to_subid_;
        IntVectRowMaj gen_to_subid_;
        IntVectRowMaj line_or_to_subid_;
        IntVectRowMaj line_ex_to_subid_;
        IntVectRowMaj trafo_hv_to_subid_;
        IntVectRowMaj trafo_lv_to_subid_;
        IntVectRowMaj storage_to_subid_;

};

#endif  //GRIDMODEL_H
