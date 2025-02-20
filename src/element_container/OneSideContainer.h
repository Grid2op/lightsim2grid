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

// same for all
// - X nb 
// - X get_bus
// - get_buses
// - get_res
// - get_res_full
// - get_theta
// - get_status
// - get_bus_id
// - reconnect_connected_buses
// - update_bus_status
// - gen_p_per_bus

// same public api but need overriden in private api
// - deactivate
// - reactivate
// - change_bus
// - change_p
// - change_q
// - reset_results
// - compute_results

// need to modify in overriden class
// - get_state
// - set_state
// - init


class OneSideContainer : public GenericContainer
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    // regular implementation
    public:
    OneSideContainer() {};

    // public generic API
    int nb() const { return static_cast<int>(p_mw_.size()); }
    int get_bus(int el_id) {return _get_bus(el_id, status_, bus_id_);}
    Eigen::Ref<const IntVect> get_buses() const {return bus_id_;}

    tuple3d get_res() const {return tuple3d(res_p_, res_q_, res_v_);}
    tuple4d get_res_full() const {return tuple4d(res_p_, res_q_, res_v_, res_theta_);}
    
    Eigen::Ref<const RealVect> get_theta() const {return res_theta_;}
    const std::vector<bool>& get_status() const {return status_;}
    Eigen::Ref<const Eigen::VectorXi> get_bus_id() const {return bus_id_;}
    void reconnect_connected_buses(std::vector<bool> & bus_status) const{
        const int nb_els = nb();
        for(int el_id = 0; el_id < nb_els; ++el_id)
        {
            if(!status_[el_id]) continue;
            const auto my_bus = bus_id_(el_id);
            if(my_bus == _deactivated_bus_id){
                // TODO DEBUG MODE only this in debug mode
                std::ostringstream exc_;
                exc_ << "OneSideContainer::reconnect_connected_buses: element with id ";
                exc_ << el_id;
                exc_ << " is connected to bus '-1' (meaning disconnected) while you said it was disconnected. Have you called `gridmodel.deactivate_xxx(...)` ?.";
                throw std::runtime_error(exc_.str());
            }
            bus_status[my_bus] = true;  // this bus is connected
        }
    }
    void disconnect_if_not_in_main_component(std::vector<bool> & busbar_in_main_component){
        const int nb_el = nb();
        SolverControl unused_solver_control;
        for(int el_id = 0; el_id < nb_el; ++el_id)
        {
            if(!status_[el_id]) continue;
            const auto my_bus = bus_id_(el_id);
            if(!busbar_in_main_component[my_bus]){
                deactivate(el_id, unused_solver_control);
            }
        }    
    }
    void update_bus_status(std::vector<bool> & bus_status) const {
        const int nb_ = nb();
        for(int el_id = 0; el_id < nb_; ++el_id)
        {
            if(!status_[el_id]) continue;
            bus_status[bus_id_[el_id]] = true;
        }
    }    

    // base function that can be called
    void gen_p_per_bus(std::vector<real_type> & res) const
    {
        const int nb_gen = nb();
        for(int sgen_id = 0; sgen_id < nb_gen; ++sgen_id)
        {
            if(!status_[sgen_id]) continue;
            const auto my_bus = bus_id_(sgen_id);
            res[my_bus] += p_mw_(sgen_id);
        }
    }

    void deactivate(int el_id, SolverControl & solver_control) {
        if(status_[el_id]){
            solver_control.tell_recompute_sbus();
        }
        this->_deactivate(el_id, solver_control);
        _generic_deactivate(el_id, status_);
    }
    void reactivate(int el_id, SolverControl & solver_control) {
        if(!status_[el_id]){
            solver_control.tell_recompute_sbus();
        }
        this->_reactivate(el_id, solver_control);
        _generic_reactivate(el_id, status_);
    }
    void change_bus(int load_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {
        this->_change_bus(load_id, new_bus_id, solver_control, nb_bus);
        _generic_change_bus(load_id, new_bus_id, bus_id_, solver_control, nb_bus);
    }
    void change_p(int el_id, real_type new_p, SolverControl & solver_control){
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
    void change_p_nothrow(int load_id, real_type new_p, SolverControl & solver_control)
    {
        bool my_status = status_.at(load_id); // and this check that el_id is not out of bound
        this->_change_p(load_id, new_p, my_status, solver_control);
        if (p_mw_(load_id) != new_p) {
            solver_control.tell_recompute_sbus();
            p_mw_(load_id) = new_p;
        }
    }
    void change_q(int el_id, real_type new_q, SolverControl & solver_control)
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
    void change_q_nothrow(int load_id, real_type new_q, SolverControl & solver_control)
    {
        bool my_status = status_.at(load_id); // and this check that el_id is not out of bound
        this->_change_q(load_id, new_q, my_status, solver_control);
        if (q_mvar_(load_id) != new_q) {
            solver_control.tell_recompute_sbus();
            q_mvar_(load_id) = new_q;
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

    protected:

        typedef std::tuple<
        std::vector<std::string>,
        std::vector<real_type>, // p_mw
        std::vector<real_type>, // q_mvar
        std::vector<int>, // bus_id
        std::vector<bool> // status
        >  StateRes;
        OneSideContainer::StateRes get_osc_state() const  // osc: one side element
        {
            std::vector<real_type> p_mw(p_mw_.begin(), p_mw_.end());
            std::vector<real_type> q_mvar(q_mvar_.begin(), q_mvar_.end());
            std::vector<int> bus_id(bus_id_.begin(), bus_id_.end());
            std::vector<bool> status = status_;
            OneSideContainer::StateRes res(names_, p_mw, q_mvar, bus_id, status);
            return res;
        }

        void set_osc_state(OneSideContainer::StateRes & my_state)  // osc: one side element
        {
            // read data
            names_ = std::get<0>(my_state);
            std::vector<real_type> & p_mw = std::get<1>(my_state);
            std::vector<real_type> & q_mvar = std::get<2>(my_state);
            std::vector<int> & bus_id = std::get<3>(my_state);
            std::vector<bool> & status = std::get<4>(my_state);

            // check sizes
            const auto size = p_mw.size();
            if(names_.size() > 0) check_size(names_, size, "names");  // names are optional
            check_size(p_mw, size, "p_mw");
            check_size(q_mvar, size, "q_mvar");
            check_size(bus_id, size, "bus_id");
            check_size(status, size, "status");

            // input data
            p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
            q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
            bus_id_ = Eigen::VectorXi::Map(&bus_id[0], bus_id.size());
            status_ = status;
        }
        
        void init_osc(const RealVect & els_p,
                    const RealVect & els_q,
                    const Eigen::VectorXi & els_bus_id,
                    const std::string & name_el
                    )  // osc: one side element
        {
            int size = static_cast<int>(els_p.size());
            check_size(els_p, size, name_el + "_p");
            check_size(els_q, size, name_el + "_q");
            check_size(els_bus_id, size, name_el + "_bus_id");

            p_mw_ = els_p;
            q_mvar_ = els_q;
            bus_id_ = els_bus_id;
            status_ = std::vector<bool>(els_p.size(), true);
        }

        virtual void _reset_results() {};
        virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<int> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac) {};
        virtual void _deactivate(int el_id, SolverControl & solver_control) {};
        virtual void _reactivate(int el_id, SolverControl & solver_control) {};
        virtual void _change_bus(int load_id, int new_bus_id, SolverControl & solver_control, int nb_bus) {};
        virtual void _change_p(int el_id, real_type new_p, bool my_status, SolverControl & solver_control) {};
        virtual void _change_q(int el_id, real_type new_p, bool my_status,SolverControl & solver_control) {};
        // virtual void _change_v(int el_id, real_type new_p, SolverControl & solver_control) {};

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