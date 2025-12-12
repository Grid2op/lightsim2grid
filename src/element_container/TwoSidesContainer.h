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

#include "Utils.h"
#include "OneSideContainer.h"

// TODO other part of the API, like deactivate, reactivate etc.
template<class OneSideType>
class TwoSideContainer : public GenericContainer
{
    public:
        // public generic API
        int nb() const { return side_1_.nb(); }
        int get_bus_side_1(int el_id) {return side_1_.get_bus(el_id);}
        int get_bus_side_2(int el_id) {return side_2_.get_bus(el_id);}

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

        void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
            const Eigen::Index nb_el = nb();
            SolverControl unused_solver_control;
            Eigen::Ref<const IntVect> bus_side_1_id_ = get_buses_side_1();
            Eigen::Ref<const IntVect> bus_side_2_id_ = get_buses_side_2();
            for(Eigen::Index i = 0; i < nb_el; ++i){
                if(!status_global_[i]) continue;
                auto bus_side_1 = bus_side_1_id_(i);
                auto bus_side_2 = bus_side_2_id_(i);
                if(!busbar_in_main_component[bus_side_1] || !busbar_in_main_component[bus_side_2]){
                    // deactivate(i, unused_solver_control);  // TODO (when it compiles)
                }
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

    protected:
        OneSideType side_1_;
        OneSideType side_2_;

        std::vector<bool> status_global_;
};


#endif  // TWO_SIDES_CONTAINER_H
