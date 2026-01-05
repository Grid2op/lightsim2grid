// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TWO_SIDES_CONTAINER_H
#define TWO_SIDES_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "GenericContainer.hpp"

// TODO other part of the API, like deactivate, reactivate etc.
template<class OneSideType>
class TwoSidesContainer : public GenericContainer
{
    public:
        class TwoSidesInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                std::string name;
                int sub_1_id;
                int sub_2_id;
                int pos_1_topo_vect;
                int pos_2_topo_vect;

                bool connected_global;
                bool connected_1;
                bool connected_2;

                int bus_1_id;
                int bus_2_id;

                bool has_res;
                real_type res_p1_mw;
                real_type res_q1_mvar;
                real_type res_v1_kv;
                real_type res_theta1_deg;
                real_type res_p2_mw;
                real_type res_q2_mvar;
                real_type res_v2_kv;
                real_type res_theta2_deg;

                TwoSidesInfo(const TwoSidesContainer & r_data_two_sides, int my_id) noexcept:
                id(-1),
                name(""),
                sub_1_id(-1),
                sub_2_id(-1),
                pos_1_topo_vect(-1),
                pos_2_topo_vect(-1),
                connected_global(false),
                connected_1(false),
                connected_2(false),
                bus_1_id(_deactivated_bus_id),
                bus_2_id(_deactivated_bus_id),
                has_res(false),
                res_p1_mw(0.),
                res_q1_mvar(0.),
                res_v1_kv(0.),
                res_theta1_deg(0.),
                res_p2_mw(0.),
                res_q2_mvar(0.),
                res_v2_kv(0.),
                res_theta2_deg(0.)
                {
                    if (my_id < 0) return;
                    if (my_id >= r_data_two_sides.nb()) return;
                    id = my_id;

                    if(r_data_two_sides.names_.size()){
                        name = r_data_two_sides.names_[my_id];
                    }

                    connected_global = r_data_two_sides.status_global_[my_id];
                    
                    const auto & side_1_info = r_data_two_sides.side_1_[my_id];
                    const auto & side_2_info = r_data_two_sides.side_2_[my_id];
                    sub_1_id = side_1_info.sub_id;
                    sub_2_id = side_2_info.sub_id;
                    pos_1_topo_vect = side_1_info.pos_topo_vect;
                    pos_2_topo_vect = side_2_info.pos_topo_vect;
                    connected_1 = side_1_info.connected;
                    connected_2 = side_2_info.connected;
                    bus_1_id = side_1_info.bus_id;
                    bus_2_id = side_2_info.bus_id;

                    if(side_1_info.has_res)
                    {
                        res_p1_mw = side_1_info.res_p_mw;
                        res_q1_mvar = side_1_info.res_q_mvar;
                        res_v1_kv = side_1_info.res_v_kv;
                        res_theta1_deg = side_1_info.res_theta_deg;
                    }
                    if(side_2_info.has_res)
                    {
                        res_p2_mw = side_2_info.res_p_mw;
                        res_q2_mvar = side_2_info.res_q_mvar;
                        res_v2_kv = side_2_info.res_v_kv;
                        res_theta2_deg = side_2_info.res_theta_deg;
                    }
                }
        };

    public:
        TwoSidesContainer() noexcept :ignore_status_global_(false), synch_status_both_side_(true){}
        virtual ~TwoSidesContainer() noexcept = default;

        // public generic API
        int nb() const { return side_1_.nb(); }
        GridModelBusId get_bus_side_1(int el_id) const {return side_1_.get_bus(el_id);}
        GridModelBusId get_bus_side_2(int el_id) const {return side_2_.get_bus(el_id);}

        void init_tsc(
            const Eigen::VectorXi & els_bus1_id,
            const Eigen::VectorXi & els_bus2_id,
            const std::string & name_elements
        )  // tsc: two sides container
        {
            auto size = els_bus1_id.size();
            check_size(els_bus2_id, size, name_elements);
            side_1_.init_osc(els_bus1_id);
            side_2_.init_osc(els_bus2_id);
            status_global_ = std::vector<bool>(els_bus1_id.size(), true);
        }

        Eigen::Ref<const GlobalBusIdVect> get_buses_side_1() const {return side_1_.get_buses();}
        Eigen::Ref<const GlobalBusIdVect> get_buses_side_2() const {return side_2_.get_buses();}

        tuple3d get_res_side_1() const {return side_1_.get_res();}
        tuple3d get_res_side_2() const {return side_2_.get_res();}

        tuple4d get_res_full_side_1() const {return side_1_.get_res_full();}
        tuple4d get_res_full_side_2() const {return side_2_.get_res_full();}

        Eigen::Ref<const RealVect> get_theta_side_1() const {return side_1_.get_theta();}
        Eigen::Ref<const RealVect> get_theta_side_2() const {return side_2_.get_theta();}

        const std::vector<bool>& get_status_global() const {return status_global_;}
        const std::vector<bool>& get_status_side_1() const {return side_1_.get_status();}
        const std::vector<bool>& get_status_side_2() const {return side_2_.get_status();}

        Eigen::Ref<const GlobalBusIdVect> get_bus_id_side_1() const {return side_1_.get_bus_id();}
        Eigen::Ref<const GlobalBusIdVect> get_bus_id_side_2() const {return side_2_.get_bus_id();}

        Eigen::Ref<const IntVect> get_bus_id_side_1_numpy() const {return side_1_.get_bus_id_numpy();}
        Eigen::Ref<const IntVect> get_bus_id_side_2_numpy() const {return side_2_.get_bus_id_numpy();}

        void reconnect_connected_buses(SubstationContainer & substation) const{
            side_1_.reconnect_connected_buses(substation);
            side_2_.reconnect_connected_buses(substation);
            // TODO think about status here !
            // Do I do if status_global_, reconnect connected buses of side_1 and side_2
            // (in this case this can do nothing if side_1 or side_2 is not connected)
        }

        virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
            const int nb_el = nb();
            SolverControl unused_solver_control;
            Eigen::Ref<const GlobalBusIdVect> bus_side_1_id_ = get_buses_side_1();
            Eigen::Ref<const GlobalBusIdVect> bus_side_2_id_ = get_buses_side_2();
            for(int i = 0; i < nb_el; ++i){
                if(!status_global_[i]){
                    side_1_.deactivate(i, unused_solver_control);
                    side_2_.deactivate(i, unused_solver_control);
                    continue;
                }
                GlobalBusId bus_side_1 = bus_side_1_id_(i);
                GlobalBusId bus_side_2 = bus_side_2_id_(i);
                if(!busbar_in_main_component[bus_side_1.cast_int()]) side_1_.deactivate(i, unused_solver_control);
                if(!busbar_in_main_component[bus_side_2.cast_int()]) side_2_.deactivate(i, unused_solver_control);
            }
        }
        virtual void nb_line_end(std::vector<int> & res) const{
            const int nb_el = nb();
            for(int el_id = 0; el_id < nb_el; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;

                const GlobalBusId bus_or = get_bus_side_1(el_id);
                if(bus_or.cast_int() != _deactivated_bus_id) res[bus_or.cast_int()] += 1;
                const GlobalBusId bus_ex = get_bus_side_2(el_id);
                if(bus_ex.cast_int() != _deactivated_bus_id) res[bus_ex.cast_int()] += 1;
            }
        }

        void set_pos_topo_vect_side_1(Eigen::Ref<const IntVect> pos_topo_vect)
        {
            side_1_.set_pos_topo_vect(pos_topo_vect);
        }
        void set_pos_topo_vect_side_2(Eigen::Ref<const IntVect> pos_topo_vect)
        {
            side_2_.set_pos_topo_vect(pos_topo_vect);
        }

        void set_subid_side_1(Eigen::Ref<const IntVect> subid)
        {
            side_1_.set_subid(subid);
        }
        void set_subid_side_2(Eigen::Ref<const IntVect> subid)
        {
            side_2_.set_subid(subid);
        }

        void update_topo(
            Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> > & has_changed,
            Eigen::Ref<const Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor> > & new_values,
            SolverControl & solver_control,
            SubstationContainer & substations
        )
        {
            side_1_._check_pos_topo_vect_filled();
            side_2_._check_pos_topo_vect_filled();

            side_1_.update_topo(has_changed, new_values, solver_control, substations);
            side_2_.update_topo(has_changed, new_values, solver_control, substations);

            // set the global status
            int nb_el = nb();
            for(int el_id=0; el_id<nb_el; ++el_id)
            {
                int pos1 = side_1_.pos_topo_vect_(el_id);
                int pos2 = side_2_.pos_topo_vect_(el_id);
                if(has_changed(pos1)){
                    resolve_status(el_id, true, solver_control);
                }
                if(has_changed(pos2)){
                    resolve_status(el_id, false, solver_control);
                }
            }
        }

        // setter (states)
        // methods used within lightsim
        void deactivate(int el_id, SolverControl & solver_control) {
            this->_deactivate(el_id, solver_control);
            _generic_deactivate(el_id, status_global_);
            side_1_.deactivate(el_id, solver_control);
            side_2_.deactivate(el_id, solver_control);
            if(ignore_status_global_) status_global_[el_id] = true;
        }
        void reactivate(int el_id, SolverControl & solver_control) {
            this->_reactivate(el_id, solver_control);
            _generic_reactivate(el_id, status_global_);
            side_1_.reactivate(el_id, solver_control);
            side_2_.reactivate(el_id, solver_control);
        }

        void reset_results_tsc(){
            side_1_.reset_results();
            side_2_.reset_results();
        };

        /**
         * Change the bus on "side 1" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */        
        void change_bus_side_1(int el_id, GridModelBusId new_gridmodel_bus_id, SolverControl & solver_control, const SubstationContainer & substation) {
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 1).");
            side_1_.change_bus(el_id, new_gridmodel_bus_id, solver_control, substation);
            resolve_status(el_id, true, solver_control);
        }
        /**
         * Change the bus on "side 2" of the element el_id.
         * 
         * The bus id is given in the "gridmodel" id, not the "solver id" nor the "local id" **ie** between 0 and `n_busbar_per_sub * n_sub`.
         */  
        void change_bus_side_2(int el_id, GridModelBusId new_gridmodel_bus_id, SolverControl & solver_control, const SubstationContainer & substation) {
            // if(!status_global_[el_id]) throw std::runtime_error("Cannot change the bus of a disconnected element (" + std::to_string(el_id) + ", side 2).");
            side_2_.change_bus(el_id, new_gridmodel_bus_id, solver_control, substation);
            resolve_status(el_id, false, solver_control);
        }

        typedef std::tuple<
            bool,  // ignore_status_global_
            bool,  // synch_status_both_side_
            std::vector<std::string>,
            std::vector<bool>,          // status_global
            typename OneSideType::StateRes, // side_1
            typename OneSideType::StateRes  // side_2
            >  StateRes;

        void resolve_status(int el_id, bool side_1_modif, SolverControl & solver_control){
            OneSideType & side_modified = side_1_modif ? side_1_: side_2_;
            OneSideType & side_to_update = side_1_modif ? side_2_: side_1_;
            if(synch_status_both_side_){
                if(side_modified.get_status(el_id)){
                    // element has been reconnected
                    // I need to reconnect other side
                    side_to_update.reactivate(el_id, solver_control);
                    status_global_[el_id] = true;
                }else{
                    side_to_update.deactivate(el_id, solver_control);
                    status_global_[el_id] = false;
                }
            }
            if(ignore_status_global_) status_global_[el_id] = true;  // always true in this case
            else{
                if(side_modified.get_status(el_id) == side_to_update.get_status(el_id)){
                    status_global_[el_id] = side_modified.get_status(el_id);
                }
            }
        }

        void set_ignore_status_global(bool ignore_status_global){
            ignore_status_global_ = ignore_status_global;
        }
        bool get_ignore_status_global() const{
            return ignore_status_global_;
        }
        void set_synch_status_both_side(bool synch_status_both_side){
            synch_status_both_side_=synch_status_both_side;
        }
        bool get_synch_status_both_side() const{
            return synch_status_both_side_;
        }

    protected:
        bool ignore_status_global_;
        bool synch_status_both_side_;
        
        OneSideType side_1_;
        OneSideType side_2_;

        std::vector<bool> status_global_;

    protected:
        StateRes get_tsc_state() const  // tsc: two sides container
        {
            StateRes res(
                ignore_status_global_,
                synch_status_both_side_,
                names_,
                status_global_,
                side_1_.get_state(),
                side_2_.get_state()
            );
            return res;
        }

        void set_tsc_state(TwoSidesContainer::StateRes & my_state)  // tsc: two sides container
        {
            ignore_status_global_ = std::get<0>(my_state);
            synch_status_both_side_ = std::get<1>(my_state);
            names_ = std::get<2>(my_state);
            status_global_ = std::get<3>(my_state);
            side_1_.set_state(std::get<4>(my_state));
            side_2_.set_state(std::get<5>(my_state));
            auto size = nb();
            check_size(names_, size, "names");
            if(side_1_.nb() != size) throw std::runtime_error("Side_1 do not have the proper size");
            if(side_2_.nb() != size) throw std::runtime_error("Side_2 do not have the proper size");
        }
        virtual void _deactivate(int el_id, SolverControl & solver_control) {}
        virtual void _reactivate(int el_id, SolverControl & solver_control) {}

        // used for example in change_bus_lv(int, int solver_control, int)
        // Eigen::Ref<IntVect> get_buses_not_const_side_1() {return side_1_.get_buses_not_const();}
        // Eigen::Ref<IntVect> get_buses_not_const_side_2() {return side_2_.get_buses_not_const();}


        // DANGER ZONE, for modifiers
        // when it will be fully refactorize, it should disappear
        Eigen::Ref<RealVect> get_res_theta_side_1() {return side_1_.get_res_theta();}
        Eigen::Ref<RealVect> get_res_p_side_1() {return side_1_.get_res_p();}
        Eigen::Ref<RealVect> get_res_q_side_1() {return side_1_.get_res_q();}
        Eigen::Ref<RealVect> get_res_v_side_1() {return side_1_.get_res_v();}
        Eigen::Ref<RealVect> get_res_theta_side_2() {return side_2_.get_res_theta();}
        Eigen::Ref<RealVect> get_res_p_side_2() {return side_2_.get_res_p();}
        Eigen::Ref<RealVect> get_res_q_side_2() {return side_2_.get_res_q();}
        Eigen::Ref<RealVect> get_res_v_side_2() {return side_2_.get_res_v();}

};


#endif  // TWO_SIDES_CONTAINER_H
