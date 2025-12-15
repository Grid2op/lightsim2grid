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
#include "OneSideContainer.hpp"

// TODO other part of the API, like deactivate, reactivate etc.
template<class OneSideType>
class TwoSidesContainer : public GenericContainer
{
    protected:
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

                TwoSidesInfo(const TwoSidesContainer & r_data_two_sides, int my_id):
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
        // public generic API
        int nb() const { return side_1_.nb(); }
        int get_bus_side_1(int el_id) const {return side_1_.get_bus(el_id);}
        int get_bus_side_2(int el_id) const {return side_2_.get_bus(el_id);}

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

        Eigen::Ref<const IntVect> get_buses_side_1() const {return side_1_.get_buses();}
        Eigen::Ref<const IntVect> get_buses_side_2() const {return side_2_.get_buses();}

        tuple3d get_res_side_1() const {return side_1_.get_res();}
        tuple3d get_res_side_2() const {return side_2_.get_res();}

        tuple4d get_res_full_side_1() const {return side_1_.get_res_full();}
        tuple4d get_res_full_side_2() const {return side_2_.get_res_full();}

        Eigen::Ref<const RealVect> get_theta_side_1() const {return side_1_.get_theta();}
        Eigen::Ref<const RealVect> get_theta_side_2() const {return side_2_.get_theta();}

        const std::vector<bool>& get_status_global() const {return status_global_;}
        const std::vector<bool>& get_status_side_1() const {return side_1_.get_status();}
        const std::vector<bool>& get_status_side_2() const {return side_2_.get_status();}

        Eigen::Ref<const Eigen::VectorXi> get_bus_id_side_1() const {return side_1_.get_bus_id();}
        Eigen::Ref<const Eigen::VectorXi> get_bus_id_side_2() const {return side_2_.get_bus_id();}

        void reconnect_connected_buses(Substation & substation) const{
            side_1_.reconnect_connected_buses(substation);
            side_2_.reconnect_connected_buses(substation);
            // TODO think about status here !
            // Do I do if status_global_, reconnect connected buses of side_1 and side_2
            // (in this case this can do nothing if side_1 or side_2 is not connected)
        }

        virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
            const Eigen::Index nb_el = nb();
            SolverControl unused_solver_control;
            Eigen::Ref<const IntVect> bus_side_1_id_ = get_buses_side_1();
            Eigen::Ref<const IntVect> bus_side_2_id_ = get_buses_side_2();
            for(Eigen::Index i = 0; i < nb_el; ++i){
                if(!status_global_[i]) continue;
                auto bus_side_1 = bus_side_1_id_(i);
                auto bus_side_2 = bus_side_2_id_(i);
                if(!busbar_in_main_component[bus_side_1]) side_1_.deactivate(i, unused_solver_control);
                if(!busbar_in_main_component[bus_side_2]) side_2_.deactivate(i, unused_solver_control);
            }
        }
        virtual void nb_line_end(std::vector<int> & res) const{
            const Eigen::Index nb_el = nb();
            for(Eigen::Index el_id = 0; el_id < nb_el; ++el_id){
                // don't do anything if the element is disconnected
                if(!status_global_[el_id]) continue;
                const auto bus_or = get_bus_side_1(el_id);
                const auto bus_ex = get_bus_side_2(el_id);
                res[bus_or] += 1;
                res[bus_ex] += 1;
            }
        }

        void update_bus_status(Substation & substation) const {
            const int nb_ = nb();
            Eigen::Ref<const IntVect> bus_side_1_id_ = get_buses_side_1();
            Eigen::Ref<const IntVect> bus_side_2_id_ = get_buses_side_2();
            const std::vector<bool>& status_side_1_ = get_status_side_1();
            const std::vector<bool>& status_side_2_ = get_status_side_2();
            for(int el_id = 0; el_id < nb_; ++el_id)
            {
                if(!status_global_[el_id]) continue;
                if(status_side_1_[el_id]) substation.reconnect_bus(bus_side_1_id_(el_id));
                if(status_side_2_[el_id]) substation.reconnect_bus(bus_side_2_id_(el_id));
            }
        }   

        void set_pos_topo_vect_side_1(Eigen::Ref<const IntVect> pos_topo_vect)
        {
            side_1_.set_pos_topo_vect(pos_topo_vect);
        }
        void set_ex_pos_topo_vect(Eigen::Ref<const IntVect> pos_topo_vect)
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
            Substation & substations
        )
        {
            side_1_.update_topo(has_changed, new_values, solver_control, substations);
            side_2_.update_topo(has_changed, new_values, solver_control, substations);
        }

        virtual ~TwoSidesContainer() noexcept = default;

        // void compute_results_tsc(const Eigen::Ref<const RealVect> & Va,
        //                          const Eigen::Ref<const RealVect> & Vm,
        //                          const Eigen::Ref<const CplxVect> & V,
        //                          const std::vector<int> & id_grid_to_solver,
        //                          const RealVect & bus_vn_kv,
        //                          real_type sn_mva,
        //                          bool ac)
        // {
        //     side_1_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        //     side_2_.compute_results(Va, Vm, V, id_grid_to_solver, bus_vn_kv, sn_mva, ac);
        // }
        void reset_results_tsc(){
            side_1_.reset_results();
            side_2_.reset_results();
        };

        typedef std::tuple<
            std::vector<bool>,          // status_global
            typename OneSideContainer::StateRes, // side_1
            typename OneSideContainer::StateRes  // side_2
            >  StateRes;

    protected:
        OneSideType side_1_;
        OneSideType side_2_;

        std::vector<bool> status_global_;

    protected:
        StateRes get_tsc_state() const  // tsc: two sides container
        {
            StateRes res(
                status_global_,
                side_1_.get_osc_state(),
                side_2_.get_osc_state()
            );
            return res;
        }

        void set_tsc_state(StateRes & my_state)  // tsc: two sides container
        {
            status_global_ = std::get<0>(my_state);
            side_1_.set_osc_state(std::get<1>(my_state));
            side_2_.set_osc_state(std::get<2>(my_state));
            auto size = nb();
            if(side_1_.nb() != size) throw std::runtime_error("Side_1 do not have the proper size");
            if(side_2_.nb() != size) throw std::runtime_error("Side_2 do not have the proper size");
        }

        // used for example in change_bus_lv(int, int solver_control, int)
        Eigen::Ref<IntVect> get_buses_not_const_side_1() {return side_1_.get_buses_not_const();}
        Eigen::Ref<IntVect> get_buses_not_const_side_2() {return side_2_.get_buses_not_const();}


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
