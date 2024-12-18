// Copyright (c) 2024, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ONE_SIDE_CONTAINER_H
#define ONE_SIDE_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "GenericContainer.h"


class OneSideContainer : public GenericContainer
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    // regular implementation
    public:
    typedef std::tuple<
       std::vector<std::string>,
       std::vector<real_type>, // p_mw
       std::vector<real_type>, // q_mvar
       std::vector<int>, // bus_id
       std::vector<bool> // status
       >  StateRes;

    OneSideContainer() {};

    // pickle (python)
    OneSideContainer::StateRes get_state() const;
    void set_base_state(OneSideContainer::StateRes & my_state);

    void init_base(const RealVect & els_p,
                   const RealVect & els_q,
                   const Eigen::VectorXi & els_bus_id,
                   const std::string & name_el
                   );

    int nb() const { return static_cast<int>(p_mw_.size()); }
    int get_bus(int load_id) {return _get_bus(load_id, status_, bus_id_);}
    Eigen::Ref<const IntVect> get_buses() const {return bus_id_;}

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    tuple4d get_res_full() const {return tuple4d(res_p_, res_q_, res_v_, res_theta_);}
    
    Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
    const std::vector<bool>& get_status() const {return status_;}
    Eigen::Ref<const Eigen::VectorXi> get_bus_id() const {return bus_id_;}


    virtual void deactivate(int el_id, SolverControl & solver_control) {
        if(status_[el_id]){
            solver_control.tell_recompute_sbus();
        }
        _deactivate(el_id, status_);
    }
    virtual void reactivate(int el_id, SolverControl & solver_control) {
        if(!status_[el_id]){
            solver_control.tell_recompute_sbus();
        }
        _reactivate(el_id, status_);
    }
    virtual void change_bus(int load_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        _change_bus(load_id, new_bus_id, bus_id_, solver_control, nb_bus);
    }
    virtual void change_p(int el_id, real_type new_p, SolverControl & solver_control){
        bool my_status = status_.at(el_id); // and this check that el_id is not out of bound
        if(!my_status)
        {
            std::ostringstream exc_;
            exc_ << "OneSideContainer::change_p: Impossible to change the active value of a disconnected element (check load id ";
            exc_ << el_id;
            exc_ << ")";
            throw std::runtime_error(exc_.str());
        }
        change_p_nothrow(el_id, new_p, solver_control);
    }
    virtual void change_p_nothrow(int load_id, real_type new_p, SolverControl & solver_control)
    {
        if (p_mw_(load_id) != new_p) {
            solver_control.tell_recompute_sbus();
            p_mw_(load_id) = new_p;
        }
    }
    virtual void change_q(int el_id, real_type new_q, SolverControl & solver_control)
    {
        bool my_status = status_.at(el_id); // and this check that el_id is not out of bound
        if(!my_status)
        {
            std::ostringstream exc_;
            exc_ << "OneSideContainer::change_q: Impossible to change the reactive value of a disconnected element (check load id ";
            exc_ << el_id;
            exc_ << ")";
            throw std::runtime_error(exc_.str());
        }
        change_q_nothrow(el_id, new_q, solver_control);
    }
    virtual void change_q_nothrow(int load_id, real_type new_q, SolverControl & solver_control)
    {
        if (q_mvar_(load_id) != new_q) {
            solver_control.tell_recompute_sbus();
            q_mvar_(load_id) = new_q;
        }
    }
    virtual void reconnect_connected_buses(std::vector<bool> & bus_status) const;
    virtual void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component);

    virtual void update_bus_status(std::vector<bool> & bus_status) const {
        const int nb_ = nb();
        for(int el_id = 0; el_id < nb_; ++el_id)
        {
            if(!status_[el_id]) continue;
            bus_status[bus_id_[el_id]] = true;
        }
    }    

    void compute_results_base(const Eigen::Ref<const RealVect> & Va,
                              const Eigen::Ref<const RealVect> & Vm,
                              const Eigen::Ref<const CplxVect> & V,
                              const std::vector<int> & id_grid_to_solver,
                              const RealVect & bus_vn_kv,
                              real_type sn_mva,
                              bool ac);
    virtual void reset_results();

    protected:
        // physical properties

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

#endif  //ONE_SIDE_CONTAINER_H