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
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <complex>      // std::complex, std::conj
#include <cmath>  // for PI

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

// import data classes
#include "Utils.h"
#include "DataGeneric.h"
#include "DataLine.h"
#include "DataShunt.h"
#include "DataTrafo.h"
#include "DataLoad.h"
#include "DataGen.h"


// import newton raphson solvers using different linear algebra solvers
#include "ChooseSolver.h"

//TODO implement a BFS check to make sure the Ymatrix is "connected" [one single component]
class GridModel : public DataGeneric
{
    public:
        typedef std::tuple<
                std::vector<double>,  // bus_vn_kv
                std::vector<bool>,  // bus_status
                // powerlines
                DataLine::StateRes ,
                // shunts
                DataShunt::StateRes,
                // trafos
                DataTrafo::StateRes,
                // gens
                DataGen::StateRes,
                // loads
                DataLoad::StateRes,
                // slack bus generator id
                int
                >  StateRes;

        GridModel():need_reset_(true){};
        GridModel(const GridModel & other);
        GridModel copy(){
            GridModel res(*this);
            return res;
        }

        // solver "control"
        void change_solver(const SolverType & type){
            need_reset_ = true;
            _solver.change_solver(type);
        }
        std::vector<SolverType> available_solvers() {return _solver.available_solvers(); }

        // All methods to init this data model, all need to be pair unit when applicable
        void init_bus(const Eigen::VectorXd & bus_vn_kv, int nb_line, int nb_trafo);

        void init_powerlines(const Eigen::VectorXd & branch_r,
                             const Eigen::VectorXd & branch_x,
                             const Eigen::VectorXcd & branch_h,
                             const Eigen::VectorXi & branch_from_id,
                             const Eigen::VectorXi & branch_to_id
                             ){
            powerlines_.init(branch_r, branch_x, branch_h, branch_from_id, branch_to_id);
        }
        void init_shunt(const Eigen::VectorXd & shunt_p_mw,
                        const Eigen::VectorXd & shunt_q_mvar,
                        const Eigen::VectorXi & shunt_bus_id){
            shunts_.init(shunt_p_mw, shunt_q_mvar, shunt_bus_id);
        }
        void init_trafo(const Eigen::VectorXd & trafo_r,
                        const Eigen::VectorXd & trafo_x,
                        const Eigen::VectorXcd & trafo_b,
                        const Eigen::VectorXd & trafo_tap_step_pct,
//                        const Eigen::VectorXd & trafo_tap_step_degree,  //TODO handle that too!
                        const Eigen::VectorXd & trafo_tap_pos,
                        const Eigen::Vector<bool, Eigen::Dynamic> & trafo_tap_hv,  // is tap on high voltage (true) or low voltate
                        const Eigen::VectorXi & trafo_hv_id,
                        const Eigen::VectorXi & trafo_lv_id
                        ){
            trafos_.init(trafo_r, trafo_x, trafo_b, trafo_tap_step_pct, trafo_tap_pos, trafo_tap_hv, trafo_hv_id, trafo_lv_id);
        }
        void init_generators(const Eigen::VectorXd & generators_p,
                             const Eigen::VectorXd & generators_v,
                             const Eigen::VectorXd & generators_min_q,
                             const Eigen::VectorXd & generators_max_q,
                             const Eigen::VectorXi & generators_bus_id){
            generators_.init(generators_p, generators_v, generators_min_q, generators_max_q, generators_bus_id);
        }
        void init_loads(const Eigen::VectorXd & loads_p,
                        const Eigen::VectorXd & loads_q,
                        const Eigen::VectorXi & loads_bus_id){
            loads_.init(loads_p, loads_q, loads_bus_id);
        }

        void add_gen_slackbus(int gen_id);

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
                std::cout << "LightSim::GridModel state size " << my_state.size() << " instead of "<< size_th << std::endl;
                // TODO more explicit error message
                throw std::runtime_error("Invalid state when loading LightSim::GridModel");
            }
        }

        //powerflows
        // dc powerflow
        Eigen::VectorXcd dc_pf(const Eigen::VectorXcd & Vinit,
                               int max_iter,  // not used for DC
                               double tol  // not used for DC
                               );

        // ac powerflow
        Eigen::VectorXcd ac_pf(const Eigen::VectorXcd & Vinit,
                               int max_iter,
                               double tol);


        // deactivate a bus. Be careful, if a bus is deactivated, but an element is
        //still connected to it, it will throw an exception
        void deactivate_bus(int bus_id) {_deactivate(bus_id, bus_status_, need_reset_); }
        // if a bus is connected, but isolated, it will make the powerflow diverge
        void reactivate_bus(int bus_id) {_reactivate(bus_id, bus_status_, need_reset_); }
        int nb_bus() const;

        //deactivate a powerline (disconnect it)
        void deactivate_powerline(int powerline_id) {powerlines_.deactivate(powerline_id, need_reset_); }
        void reactivate_powerline(int powerline_id) {powerlines_.reactivate(powerline_id, need_reset_); }
        void change_bus_powerline_or(int powerline_id, int new_bus_id) {powerlines_.change_bus_or(powerline_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        void change_bus_powerline_ex(int powerline_id, int new_bus_id) {powerlines_.change_bus_ex(powerline_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        int get_bus_powerline_or(int powerline_id) {return powerlines_.get_bus_or(powerline_id);}
        int get_bus_powerline_ex(int powerline_id) {return powerlines_.get_bus_ex(powerline_id);}

        //deactivate trafo
        void deactivate_trafo(int trafo_id) {trafos_.deactivate(trafo_id, need_reset_); }
        void reactivate_trafo(int trafo_id) {trafos_.reactivate(trafo_id, need_reset_); }
        void change_bus_trafo_hv(int trafo_id, int new_bus_id) {trafos_.change_bus_hv(trafo_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        void change_bus_trafo_lv(int trafo_id, int new_bus_id) {trafos_.change_bus_lv(trafo_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        int get_bus_trafo_hv(int trafo_id) {return trafos_.get_bus_hv(trafo_id);}
        int get_bus_trafo_lv(int trafo_id) {return trafos_.get_bus_lv(trafo_id);}

        //load
        void deactivate_load(int load_id) {loads_.deactivate(load_id, need_reset_); }
        void reactivate_load(int load_id) {loads_.reactivate(load_id, need_reset_); }
        void change_bus_load(int load_id, int new_bus_id) {loads_.change_bus(load_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        void change_p_load(int load_id, double new_p) {loads_.change_p(load_id, new_p, need_reset_); }
        void change_q_load(int load_id, double new_q) {loads_.change_q(load_id, new_q, need_reset_); }
        int get_bus_load(int load_id) {return loads_.get_bus(load_id);}

        //generator
        void deactivate_gen(int gen_id) {generators_.deactivate(gen_id, need_reset_); }
        void reactivate_gen(int gen_id) {generators_.reactivate(gen_id, need_reset_); }
        void change_bus_gen(int gen_id, int new_bus_id) {generators_.change_bus(gen_id, new_bus_id, need_reset_, bus_vn_kv_.size()); }
        void change_p_gen(int gen_id, double new_p) {generators_.change_p(gen_id, new_p, need_reset_); }
        void change_v_gen(int gen_id, double new_v_pu) {generators_.change_v(gen_id, new_v_pu, need_reset_); }
        int get_bus_gen(int gen_id) {return generators_.get_bus(gen_id);}

        //shunt
        void deactivate_shunt(int shunt_id) {shunts_.deactivate(shunt_id, need_reset_); }
        void reactivate_shunt(int shunt_id) {shunts_.reactivate(shunt_id, need_reset_); }
        void change_bus_shunt(int shunt_id, int new_bus_id) {shunts_.change_bus(shunt_id, new_bus_id, need_reset_, bus_vn_kv_.size());  }
        void change_p_shunt(int shunt_id, double new_p) {shunts_.change_p(shunt_id, new_p, need_reset_); }
        void change_q_shunt(int shunt_id, double new_q) {shunts_.change_q(shunt_id, new_q, need_reset_); }
        int get_bus_shunt(int shunt_id) {return shunts_.get_bus(shunt_id);}

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

        // get some internal information, be cerafull the ID of the buses might not be the same
        // TODO convert it back to this ID, that will make copies, but who really cares ?
        Eigen::SparseMatrix<cdouble> get_Ybus(){
            return Ybus_;
        }
        Eigen::VectorXcd get_Sbus(){
            return Sbus_;
        }
        Eigen::VectorXi get_pv(){
            return bus_pv_;
        }
        Eigen::VectorXi get_pq(){
            return bus_pq_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Va(){
            return _solver.get_Va();
        }
        Eigen::Ref<Eigen::VectorXd> get_Vm(){
            return _solver.get_Vm();
        }
        Eigen::SparseMatrix<double> get_J(){
            return _solver.get_J();
        }

        // part dedicated to grid2op backend, optimized for grid2op data representation (for speed)
        // this is not recommended to use it outside of its intended usage.
        void update_bus_status(int nb_bus_before, Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> > active_bus);
        void update_gens_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                           Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_gens_v(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                           Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_loads_p(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_loads_q(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                            Eigen::Ref<Eigen::Array<float, Eigen::Dynamic, Eigen::RowMajor> > new_values);
        void update_topo(Eigen::Ref<Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > has_changed,
                         Eigen::Ref<Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > new_values);

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
        void set_n_sub(int n_sub)
        {
            n_sub_ = n_sub;
        }

    protected:
    // add method to change topology, change ratio of transformers, change

        // compute admittance matrix
        // dc powerflow
        // void init_dcY(Eigen::SparseMatrix<double> & dcYbus);

        // ac powerflows
        Eigen::VectorXcd pre_process_solver(const Eigen::VectorXcd & Vinit);
        void init_Ybus(Eigen::SparseMatrix<cdouble> & Ybus, Eigen::VectorXcd & Sbus,
                       std::vector<int> & id_me_to_solver, std::vector<int>& id_solver_to_me,
                       int & slack_bus_id_solver);
        void fillYbus(Eigen::SparseMatrix<cdouble> & res, bool ac, const std::vector<int>& id_me_to_solver);
        void fillSbus_me(Eigen::VectorXcd & res, bool ac, const std::vector<int>& id_me_to_solver, int slack_bus_id_solver);
        void fillpv_pq(const std::vector<int>& id_me_to_solver);

        // results
        /**process the results from the solver to this instance
        **/
        void process_results(bool conv, Eigen::VectorXcd & res, const Eigen::VectorXcd & Vinit);

        /**
        Compute the results vector from the Va, Vm post powerflow
        **/
        void compute_results();
        /**
        reset the results in case of divergence of the powerflow.
        **/
        void reset_results();

        /**
        reset the solver, and all its results
        **/
        void reset();

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
                    (this->*fun)(el_id, static_cast<double>(new_values[el_id]));
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
                if(has_changed(el_pos))
                {
                    int new_bus = new_values(el_pos);
                    if(new_bus > 0){
                        // new bus is a real bus, so i need to make sure to have it turned on, and then change the bus
                        int init_bus_me = vect_subid(el_id);
                        int new_bus_backend = new_bus == 1 ? init_bus_me : init_bus_me + n_sub_ ;
                        (this->*fun_react)(el_id); // eg reactivate_load(load_id);
                        (this->*fun_change)(el_id, new_bus_backend); // eg change_bus_load(load_id, new_bus_backend);
                    } else{
                        // new bus is negative, we deactivate it
                        (this->*fun_deact)(el_id);// eg deactivate_load(load_id);
                    }
                }
            }
        }

    protected:
        // member of the grid
        // static const int _deactivated_bus_id;
        bool need_reset_;

        // powersystem representation
        // 1. bus
        Eigen::VectorXd bus_vn_kv_;
        std::vector<bool> bus_status_;

        // always have the length of the number of buses,
        // id_me_to_model_[id_me] gives -1 if the bus "id_me" is deactivated, or "id_model" if it is activated.
        std::vector<int> id_me_to_solver_;
        // convert the bus id from the model to the bus id of me.
        // it has a variable size, that depends on the number of connected bus. if "id_model" is an id of a bus
        // sent to the solver, then id_model_to_me_[id_model] is the bus id of this model of the grid.
        std::vector<int> id_solver_to_me_;

        // 2. powerline
        DataLine powerlines_;

        // 3. shunt
        DataShunt shunts_;

        // 4. transformers
        // have the r, x, h and ratio
        // ratio is computed from the tap, so maybe store tap num and tap_step_pct
        DataTrafo trafos_;

        // 5. generators
        DataGen generators_;

        // 6. loads
        DataLoad loads_;

        // 7. slack bus
        int gen_slackbus_;
        int slack_bus_id_;
        int slack_bus_id_solver_;

        // as matrix, for the solver
        Eigen::SparseMatrix<cdouble> Ybus_;
        Eigen::VectorXcd Sbus_;
        Eigen::VectorXi bus_pv_;  // id are the solver internal id and NOT the initial id
        Eigen::VectorXi bus_pq_;  // id are the solver internal id and NOT the initial id

        // TODO have version of the stuff above for the public api, indexed with "me" and not "solver"

        // to solve the newton raphson
        ChooseSolver _solver;

        // specific grid2op
        int n_sub_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> load_pos_topo_vect_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> gen_pos_topo_vect_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> line_or_pos_topo_vect_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> line_ex_pos_topo_vect_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> trafo_hv_pos_topo_vect_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> trafo_lv_pos_topo_vect_;

        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> load_to_subid_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> gen_to_subid_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> line_or_to_subid_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> line_ex_to_subid_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> trafo_hv_to_subid_;
        Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> trafo_lv_to_subid_;

};

#endif  //GRIDMODEL_H
