// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SHUNT_CONTAINER_H
#define SHUNT_CONTAINER_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "OneSideContainer.h"

/**
This class is a container for all shunts on the grid.

The convention used for the shunt is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/shunt.html#electric-model
**/
class ShuntContainer : public OneSideContainer
{
    // iterators part
    public:
        class ShuntInfo
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

                ShuntInfo(const ShuntContainer & r_data_shunt, int my_id):
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
                    if((my_id >= 0) & (my_id < r_data_shunt.nb()))
                    {
                        id = my_id;
                        if(r_data_shunt.names_.size()){
                            name = r_data_shunt.names_[my_id];
                        }
                        connected = r_data_shunt.status_[my_id];
                        if(connected)  bus_id = r_data_shunt.bus_id_[my_id];

                        target_p_mw = r_data_shunt.p_mw_.coeff(my_id);
                        target_q_mvar = r_data_shunt.q_mvar_.coeff(my_id);

                        has_res = r_data_shunt.res_p_.size() > 0;
                        if(has_res)
                        {
                            res_p_mw = r_data_shunt.res_p_.coeff(my_id);
                            res_q_mvar = r_data_shunt.res_q_.coeff(my_id);
                            res_v_kv = r_data_shunt.res_v_.coeff(my_id);
                            res_theta_deg = r_data_shunt.res_theta_.coeff(my_id);
                        }
                    }
                }
        };
        typedef ShuntInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<ShuntContainer> ShuntContainerConstIterator;

    public:
        typedef ShuntContainerConstIterator const_iterator_type;
        const_iterator_type begin() const {return ShuntContainerConstIterator(this, 0); }
        const_iterator_type end() const {return ShuntContainerConstIterator(this, nb()); }
        ShuntInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative load id.");
            }
            if(id >= nb())
            {
                throw std::range_error("Generator out of bound. Not enough loads on the grid.");
            }
            return ShuntInfo(*this, id);
        }

    public:
    typedef std::tuple<OneSideContainer::StateRes >  StateRes;

    ShuntContainer():OneSideContainer() {};


    void init(const RealVect & shunt_p_mw,
              const RealVect & shunt_q_mvar,
              const Eigen::VectorXi & shunt_bus_id
              )
    {
        init_osc(shunt_p_mw,
                 shunt_q_mvar,
                 shunt_bus_id,
                 "shunts");
        reset_results();
    }

    // pickle (python)
    ShuntContainer::StateRes get_state() const;
    void set_state(ShuntContainer::StateRes & my_state );
    
    virtual void fillYbus(std::vector<Eigen::Triplet<cplx_type> > & res,
                          bool ac,
                          const std::vector<int> & id_grid_to_solver,
                          real_type sn_mva) const;
    virtual void fillBp_Bpp(std::vector<Eigen::Triplet<real_type> > & Bp,
                            std::vector<Eigen::Triplet<real_type> > & Bpp,
                            const std::vector<int> & id_grid_to_solver,
                            real_type sn_mva,
                            FDPFMethod xb_or_bx) const;
    virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const;  // in DC i need that

    protected:

    virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                  const Eigen::Ref<const RealVect> & Vm,
                                  const Eigen::Ref<const CplxVect> & V,
                                  const std::vector<int> & id_grid_to_solver,
                                  const RealVect & bus_vn_kv,
                                  real_type sn_mva,
                                  bool ac);

    protected:
        // physical properties

        // input data

        //output data

};

#endif  //SHUNT_CONTAINER_H
