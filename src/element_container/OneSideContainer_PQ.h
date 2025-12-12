// Copyright (c) 2024, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ONE_SIDE_CONTAINER_PQ_H
#define ONE_SIDE_CONTAINER_PQ_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.h"
#include "OneSideContainer.h"


/**
 * This class represents a "one side container"
 * with added information about target_p and target_q.
 * 
 * It is used for loads and shunts for example.
 */
class OneSideContainer_PQ : public OneSideContainer
{
    // TODO make a single class for load and shunt and just specialize the part where the
    // TODO powerflow equations are located (when i update the Y matrix)

    protected:
        class OneSidePQInfo: public OneSideContainer::OneSideInfo
        {
            public:
                real_type target_p_mw;
                real_type target_q_mvar;

                OneSidePQInfo(const OneSideContainer_PQ & r_data_pq, int my_id):
                OneSideInfo(r_data_pq, my_id),
                target_p_mw(0.),
                target_q_mvar(0.)
                {
                    if((my_id >= 0) & (my_id < r_data_pq.nb()))
                    {
                        target_p_mw = r_data_pq.target_p_mw_.coeff(my_id);
                        target_q_mvar = r_data_pq.target_q_mvar_.coeff(my_id);
                    }
                }
        };
    
    // regular implementation
    public:
        OneSideContainer_PQ() {};

        // public generic API

        Eigen::Ref<const RealVect> get_target_p() const {return target_p_mw_;}

        // base function that can be called
        void gen_p_per_bus(std::vector<real_type> & res) const
        {
            const int nb_gen = nb();
            for(int sgen_id = 0; sgen_id < nb_gen; ++sgen_id)
            {
                if(!status_[sgen_id]) continue;
                const auto my_bus = bus_id_(sgen_id);
                res[my_bus] += target_p_mw_(sgen_id);
            }
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
            if (target_p_mw_(load_id) != new_p) {
                solver_control.tell_recompute_sbus();
                target_p_mw_(load_id) = new_p;
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
            if (target_q_mvar_(load_id) != new_q) {
                solver_control.tell_recompute_sbus();
                target_q_mvar_(load_id) = new_q;
            }
        }

        typedef std::tuple<
            OneSideContainer::StateRes,
            std::vector<real_type>, // p_mw
            std::vector<real_type> // q_mvar
            >  StateRes;

    protected:
        OneSideContainer_PQ::StateRes get_osc_pq_state() const  // osc: one side element
        {
            std::vector<real_type> target_p_mw(target_p_mw_.begin(), target_p_mw_.end());
            std::vector<real_type> target_q_mvar(target_q_mvar_.begin(), target_q_mvar_.end());
            OneSideContainer_PQ::StateRes res(
                get_osc_state(),
                target_p_mw,
                target_q_mvar);
            return res;
        }

        void set_osc_pq_state(OneSideContainer_PQ::StateRes & my_state)  // osc: one side element
        {
            // read data from my_state
            set_osc_state(std::get<0>(my_state));

            // init target_p and target_q
            std::vector<real_type> & p_mw = std::get<1>(my_state);
            std::vector<real_type> & q_mvar = std::get<2>(my_state);

            // check sizes
            const auto size = nb();
            check_size(p_mw, size, "p_mw");
            check_size(q_mvar, size, "q_mvar");

            // input data
            target_p_mw_ = RealVect::Map(&p_mw[0], p_mw.size());
            target_q_mvar_ = RealVect::Map(&q_mvar[0], q_mvar.size());
        }
        
        void init_osc_pq(const RealVect & els_p,
                         const RealVect & els_q,
                         const Eigen::VectorXi & els_bus_id,
                         const std::string & name_el
                         )  // osc: one side element
        {
            init_osc(els_bus_id);
            int size = nb();
            check_size(els_p, size, name_el + "_p");
            check_size(els_q, size, name_el + "_q");

            target_p_mw_ = els_p;
            target_q_mvar_ = els_q;
        }

        void set_osc_pq_res_p(){
            res_p_ = target_p_mw_;
            set_osc_res_p();
        }

        void set_osc_pq_res_q(bool ac){
            if(ac) res_q_ = target_q_mvar_;
            set_osc_res_q(ac);
        }

    protected:
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

        // data for grid2op compat

        // input data
        RealVect target_p_mw_;
        RealVect target_q_mvar_;

};

#endif  //ONE_SIDE_CONTAINER_PQ_H