// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LOAD_CONTAINER_H
#define LOAD_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "OneSideContainer.h"

/**
This class is a container for all loads on the grid.

The convention used for the generator is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/load.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/load.html#electric-model

NOTE: this class is also used for the storage units! So storage units are modeled as load
which entails that negative storage: the unit is discharging, power is injected in the grid,
positive storage: the unit is charging, power is taken from the grid.
**/
class LoadContainer : public OneSideContainer
{
    // iterators part
    public:
        class LoadInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                std::string name;
                bool connected;
                int bus_id;

                real_type target_p_mw;
                real_type target_q_mvar;
                bool has_res;
                real_type res_p_mw;
                real_type res_q_mvar;
                real_type res_v_kv;
                real_type res_theta_deg;

                LoadInfo(const LoadContainer & r_data_load, int my_id):
                id(-1),
                name(""),
                connected(false),
                bus_id(_deactivated_bus_id),
                target_p_mw(0.),
                target_q_mvar(0.),
                has_res(false),
                res_p_mw(0.),
                res_q_mvar(0.),
                res_v_kv(0.),
                res_theta_deg(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_load.nb()))
                    {
                        id = my_id;
                        if(r_data_load.names_.size()){
                            name = r_data_load.names_[my_id];
                        }
                        connected = r_data_load.status_[my_id];
                        if(connected) bus_id = r_data_load.bus_id_[my_id];

                        target_p_mw = r_data_load.p_mw_.coeff(my_id);
                        target_q_mvar = r_data_load.q_mvar_.coeff(my_id);

                        has_res = r_data_load.res_p_.size() > 0;
                        if(has_res)
                        {
                            res_p_mw = r_data_load.res_p_.coeff(my_id);
                            res_q_mvar = r_data_load.res_q_.coeff(my_id);
                            res_v_kv = r_data_load.res_v_.coeff(my_id);
                            res_theta_deg = r_data_load.res_theta_.coeff(my_id);
                        }
                    }
                }
        };
        typedef LoadInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<LoadContainer> LoadContainerConstIterator;

    public:
        typedef LoadContainerConstIterator const_iterator_type;
        const_iterator_type begin() const {return LoadContainerConstIterator(this, 0); }
        const_iterator_type end() const {return LoadContainerConstIterator(this, nb()); }
        LoadInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative load id.");
            }
            if(id >= nb())
            {
                throw std::range_error("Generator out of bound. Not enough loads on the grid.");
            }
            return LoadInfo(*this, id);
        }

    // regular implementation
    public:
    typedef std::tuple<
       OneSideContainer::StateRes  // state of the base class 
       >  StateRes;

    LoadContainer():OneSideContainer(){};

    // pickle (python)
    LoadContainer::StateRes get_state() const;
    void set_state(LoadContainer::StateRes & my_state);

    void init(const RealVect & load_p_mw,
              const RealVect & load_q_mvar,
              const Eigen::VectorXi & load_bus_id
              )
    {
        init_osc(load_p_mw,
                 load_q_mvar,
                 load_bus_id,
                 "loads");
        reset_results();
    }

    virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const;

    protected:
    virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                  const Eigen::Ref<const RealVect> & Vm,
                                  const Eigen::Ref<const CplxVect> & V,
                                  const std::vector<int> & id_grid_to_solver,
                                  const RealVect & bus_vn_kv,
                                  real_type sn_mva,
                                  bool ac)
                                  {

                                        set_osc_res_p();
                                        set_osc_res_q(ac);
                                  }
};

#endif  //LOAD_CONTAINER_H
