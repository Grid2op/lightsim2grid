// Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SGEN_CONTAINER_H
#define SGEN_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "GenericContainer.h"

/**
This class is a container for all static generator (PQ generators) on the grid.
They are given in the generator convention: positive sign for P,Q => the power is produced.

The convention used for the static is the same as in pandapower:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html

and for modeling of the Ybus matrix:
https://pandapower.readthedocs.io/en/latest/elements/sgen.html#electric-model
**/
class SGenContainer: public GenericContainer
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)


    // iterators part
    public:
        class SGenInfo
        {
            public:
                // members
                // TODO add some const here (value should not be changed !) !!!
                int id;  // id of the generator
                std::string name;
                bool connected;
                int bus_id;

                real_type min_q_mvar;
                real_type max_q_mvar;
                real_type min_p_mw;
                real_type max_p_mw;

                real_type target_p_mw;
                real_type target_q_mvar;

                bool has_res;
                real_type res_p_mw;
                real_type res_q_mvar;
                real_type res_v_kv;
                real_type res_theta_deg;

                SGenInfo(const SGenContainer & r_data_sgen, int my_id):
                id(-1),
                name(""),
                connected(false),
                bus_id(-1),
                min_q_mvar(0.),
                max_q_mvar(0.),
                min_p_mw(0.),
                max_p_mw(0.),
                target_p_mw(0.),
                target_q_mvar(0.),
                has_res(false),
                res_p_mw(0.),
                res_q_mvar(0.),
                res_v_kv(0.),
                res_theta_deg(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_sgen.nb()))
                    {
                        id = my_id;
                        if(r_data_sgen.names_.size()){
                            name = r_data_sgen.names_[my_id];
                        }
                        connected = r_data_sgen.status_[my_id];
                        bus_id = r_data_sgen.bus_id_[my_id];

                        min_q_mvar = r_data_sgen.q_min_mvar_(my_id);
                        max_q_mvar = r_data_sgen.q_max_mvar_(my_id);
                        min_p_mw = r_data_sgen.p_min_mw_(my_id);
                        max_p_mw = r_data_sgen.p_max_mw_(my_id);

                        target_p_mw = r_data_sgen.p_mw_.coeff(my_id);
                        target_q_mvar = r_data_sgen.q_mvar_.coeff(my_id);

                        has_res = r_data_sgen.res_p_.size() > 0;
                        if(has_res)
                        {
                            res_p_mw = r_data_sgen.res_p_.coeff(my_id);
                            res_q_mvar = r_data_sgen.res_q_.coeff(my_id);
                            res_v_kv = r_data_sgen.res_v_.coeff(my_id);
                            res_theta_deg = r_data_sgen.res_theta_.coeff(my_id);
                        }
                    }
                }
        };
        typedef SGenInfo DataInfo;

    private:
        typedef GenericContainerConstIterator<SGenContainer> SGenContainerConstIterator;

    public:
        typedef SGenContainerConstIterator const_iterator_type;
        const_iterator_type begin() const {return SGenContainerConstIterator(this, 0); }
        const_iterator_type end() const {return SGenContainerConstIterator(this, nb()); }
        SGenInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative static generator");
            }
            if(id >= nb())
            {
                throw std::range_error("Generator out of bound. Not enough static generators on the grid.");
            }
            return SGenInfo(*this, id);
        }

    public:
    typedef std::tuple<
       std::vector<std::string>,
       std::vector<real_type>, // p_mw
       std::vector<real_type>, // q_mvar
       std::vector<real_type>, // p_min
       std::vector<real_type>, //  p_max
       std::vector<real_type>, //  q_min
       std::vector<real_type>, //  q_max
       std::vector<int>, // bus_id
       std::vector<bool> // status
       >  StateRes;

    SGenContainer() {};

    // pickle (python)
    SGenContainer::StateRes get_state() const;
    void set_state(SGenContainer::StateRes & my_state );


    void init(const RealVect & sgen_p,
              const RealVect & sgen_q,
              const RealVect & sgen_pmin,
              const RealVect & sgen_pmax,
              const RealVect & sgen_qmin,
              const RealVect & sgen_qmax,
              const Eigen::VectorXi & sgen_bus_id
              );

    int nb() const { return static_cast<int>(p_mw_.size()); }

    void deactivate(int sgen_id, SolverControl & solver_control) {
        if(status_[sgen_id]){
            solver_control.tell_recompute_sbus();
        }
        _deactivate(sgen_id, status_);
    }
    void reactivate(int sgen_id, SolverControl & solver_control) {
        if(!status_[sgen_id]){
            solver_control.tell_recompute_sbus();
        }
        _reactivate(sgen_id, status_);
    }
    void change_bus(int sgen_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {_change_bus(sgen_id, new_bus_id, bus_id_, solver_control, nb_bus);}
    int get_bus(int sgen_id) {return _get_bus(sgen_id, status_, bus_id_);}
    void change_p(int sgen_id, real_type new_p, SolverControl & solver_control);
    void change_q(int sgen_id, real_type new_q, SolverControl & solver_control);
    virtual void reconnect_connected_buses(std::vector<bool> & bus_status) const;
    virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component);
    
    virtual void fillSbus(CplxVect & Sbus, const std::vector<int> & id_grid_to_solver, bool ac) const ;
    virtual void gen_p_per_bus(std::vector<real_type> & res) const;
    virtual void update_bus_status(std::vector<bool> & bus_status) const {
        const int nb_ = nb();
        for(int el_id = 0; el_id < nb_; ++el_id)
        {
            if(!status_[el_id]) continue;
            bus_status[bus_id_[el_id]] = true;
        }
    }    

    void compute_results(const Eigen::Ref<const RealVect> & Va,
                         const Eigen::Ref<const RealVect> & Vm,
                         const Eigen::Ref<const CplxVect> & V,
                         const std::vector<int> & id_grid_to_solver,
                         const RealVect & bus_vn_kv,
                         real_type sn_mva,
                         bool ac);
    void reset_results();

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    tuple4d get_res_full() const {return tuple4d(res_p_, res_q_, res_v_, res_theta_);}
    Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
    const std::vector<bool>& get_status() const {return status_;}
    const Eigen::VectorXi & get_bus_id() const {return bus_id_;}

    protected:
        // physical properties
        RealVect p_min_mw_;
        RealVect p_max_mw_;
        RealVect q_min_mvar_;
        RealVect q_max_mvar_;

        // input data
        RealVect p_mw_;
        RealVect q_mvar_;
        Eigen::VectorXi bus_id_;
        std::vector<bool> status_;

        //output data
        RealVect res_p_;  // in MW
        RealVect res_q_;  // in MVar
        RealVect res_v_;  // in kV
        RealVect res_theta_;  // in degree
};

#endif  //SGEN_CONTAINER_H
