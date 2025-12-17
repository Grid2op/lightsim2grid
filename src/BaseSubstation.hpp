// Copyright (c) 2025, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASESUBSTATION_H
#define BASESUBSTATION_H

#include <iostream>
#include <vector>
// #include <set>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI

#include "Utils.hpp"
#include "BaseConstants.hpp"
// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"


class Substation : protected BaseConstants
{
    public:

    typedef std::tuple<
       int,  // n_sub_
       int,  // nmax_busbar_per_sub
       std::vector<real_type>, // sub_vn_kv_;
       std::vector<bool>,  // bus_status_;
       std::vector<real_type>,  // bus_vn_kv_;
       std::vector<std::string>  // sub_names_
       > StateRes;
    

    Substation::StateRes get_state() const;
    void set_state(Substation::StateRes & my_state);

    Substation():
        n_sub_(-1), 
        nmax_busbar_per_sub_(-1),
        n_bus_max_(-1){}
        
    Substation(int n_sub, int nmax_busbar_per_sub):
        n_sub_(n_sub),
        nmax_busbar_per_sub_(nmax_busbar_per_sub),
        n_bus_max_(n_sub * nmax_busbar_per_sub),
        sub_vn_kv_(n_sub),
        bus_status_(n_bus_max_, false),
        bus_vn_kv_(n_bus_max_){}

    void reset_bus_status(){
        for(auto i = 0; i < n_bus_max_; ++ i) bus_status_[i] = -1;
    }

    void init_sub(const RealVect & sub_vn_kv){
        for(auto i = 0; i < n_sub_; ++i){
            // store substation vn kv
            sub_vn_kv_[i] = sub_vn_kv[i];
            for(auto j = 0; j < nmax_busbar_per_sub_; ++j){
                // store the buses vn kv, for all buses
                bus_vn_kv_[i + j * n_sub_] = sub_vn_kv[i];
            }
        }

    }

    void init_sub_names(const std::vector<std::string> & sub_names){
        if(sub_names.size() != static_cast<size_t>(n_sub_)){
            throw std::runtime_error("Wrong number of substation when setting their names.");
        }
        sub_names_ = sub_names;
    }
    const std::vector<std::string> & get_sub_names() const {
        return sub_names_;
    }

    // void from_agent_topology(int sub_id, int local_bus_id){
    //     bus_status_[sub_id + local_bus_id * n_sub_] = true;

    // }

    // return the number of substations on the grid
    int nb_sub() const {return n_sub_;}

    // return the maximum number of possible buses on the grid
    unsigned int nb_bus() const {return bus_vn_kv_.size();}
    int nmax_busbar_per_sub() const {return nmax_busbar_per_sub_;}

    Eigen::Ref<const RealVect> get_bus_vn_kv() const {return bus_vn_kv_;}
    bool is_bus_connected(int global_bus_id) const {return bus_status_[global_bus_id];}
    bool is_bus_connected(int sub_id, int local_bus_id) const {return bus_status_[local_to_global(sub_id, local_bus_id)];}

    void init_bus(int n_sub, int nmax_busbar_per_sub, const RealVect & bus_vn_kv)
    {
        if(bus_vn_kv.size() != n_sub * nmax_busbar_per_sub){
            std::ostringstream exc_;
            exc_ << "Substation::init_bus: ";
            exc_ << "your model counts ";
            exc_ << n_sub_  << " substations and ";
            exc_ << nmax_busbar_per_sub_  << " maximum busbars per substation. But you provided a bus_vn_kv with ";
            exc_ << bus_vn_kv.size() << " elements.";
            exc_ << "Both should match.";
            throw std::runtime_error(exc_.str());
        }

        n_sub_ = n_sub;
        nmax_busbar_per_sub_ = nmax_busbar_per_sub;
        bus_vn_kv_ = bus_vn_kv;  // base_kv

        bus_status_ = std::vector<bool>(nb_bus(), true); // by default everything is connected

        // check that a "substation" always has the same vn_kv for all of its buses
        for(int sub_id=0; sub_id < n_sub; sub_id++){
            real_type ref_vn_kv = bus_vn_kv(sub_id);
            for(int bus_id=1; bus_id < nmax_busbar_per_sub; bus_id++)
            {
                if(abs(bus_vn_kv(local_to_global(sub_id, bus_id)) - ref_vn_kv) > _tol_equal_float){
                    const std::string msg = R"mydelimiter(
                    Each bus of each substation must have the same nominal voltage. 
                    Check substation )mydelimiter" + 
                    std::to_string(sub_id) + 
                    " and buses 0 and " + 
                    std::to_string(bus_id)+".";

                    throw std::runtime_error(msg);
                }
            }
        }
    }
    const std::vector<bool> & get_bus_status() const {return bus_status_;}
    /**
    Retrieve the number of connected buses
    **/
    int nb_connected_bus() const
    {
        int res = 0;
        for(const auto & el : bus_status_)
        {
            if(el) ++res;
        }
        return res;
    }
    
    void disconnect_all_buses(){
        for(unsigned int i = 0; i < nb_bus(); ++i) bus_status_[i] = false;
    }
    void reconnect_bus(int global_bus_id){
        bus_status_[global_bus_id] = true;
    }
    void reconnect_bus(int sub_id, int local_bus_id){
        reconnect_bus(local_to_global(sub_id, local_bus_id));
    }
    void disconnect_bus(int global_bus_id){
        bus_status_[global_bus_id] = false;
    }
    void disconnect_bus(int sub_id, int local_bus_id){
        disconnect_bus(local_to_global(sub_id, local_bus_id));
    }

    protected:
        int local_to_global(int sub_id, int local_bus_id) const{
            return sub_id + local_bus_id * n_sub_;
        }

    private:
        int n_sub_;
        int nmax_busbar_per_sub_;
        int n_bus_max_;
        RealVect sub_vn_kv_;
        std::vector<bool> bus_status_;
        RealVect bus_vn_kv_;
        std::vector<std::string> sub_names_;

};

#endif // BASESUBSTATION_H
