// Copyright (c) 2024, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ONE_SIDE_CONTAINER_FORBRANCH_H
#define ONE_SIDE_CONTAINER_FORBRANCH_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer.hpp"


/**
 * This class represents a "one side container". It handles "properly"
 * the "solver_control" when connecting or disconnecting powerlines.
 * 
 * TODO add res_a and other stuff (like ydc_11 and others here instead of in TwoSidesContainer_rxh_A
 */
class OneSideContainer_ForBranch : public OneSideContainer
{
    protected:
        class OneSideForBranchInfo: public OneSideContainer::OneSideInfo
        {
            public:
                OneSideForBranchInfo(const OneSideContainer_ForBranch & r_data_pq, int my_id):
                OneSideInfo(r_data_pq, my_id) {}
        };
    
    /////////////////////////////
    // iterator
    private:
        typedef GenericContainerConstIterator<OneSideContainer_ForBranch> OSCC4BonstIterator;

    public:
        OSCC4BonstIterator begin() const {return OSCC4BonstIterator(this, 0); }
        OSCC4BonstIterator end() const {return OSCC4BonstIterator(this, nb()); }
        OneSideForBranchInfo operator[](int id) const
        {
            if(id < 0)
            {
                throw std::range_error("You cannot ask for a negative load id.");
            }
            if(id >= nb())
            {
                throw std::range_error("Load out of bound. Not enough loads on the grid.");
            }
            return OneSideForBranchInfo(*this, id);
        }
    ////////////////////////////

    // regular implementation
    public:
        OneSideContainer_ForBranch() {};

        // public generic API

        typedef std::tuple<
            OneSideContainer::StateRes
            >  StateRes;

        StateRes get_state() const
        {
            return get_osc_forB_state();
        }

        void set_state(StateRes & state)
        {
            set_osc_forB_state(state);
        }

    protected:
        OneSideContainer_ForBranch::StateRes get_osc_forB_state() const  // osc: one side element
        {
            OneSideContainer_ForBranch::StateRes res(
                get_osc_state());
            return res;
        }

        void set_osc_forB_state(OneSideContainer_ForBranch::StateRes & my_state)  // osc: one side element
        {
            // read data from my_state
            set_osc_state(std::get<0>(my_state));
        }
        
        void init_osc_forB(
            const RealVect & els_p,
            const RealVect & els_q,
            const Eigen::VectorXi & els_bus_id,
            const std::string & name_el
            )  // osc: one side element
        {
            init_osc(els_bus_id);
        }

    protected:
        virtual void _reset_results() {};
        virtual void _compute_results(const Eigen::Ref<const RealVect> & Va,
                                      const Eigen::Ref<const RealVect> & Vm,
                                      const Eigen::Ref<const CplxVect> & V,
                                      const std::vector<SolverBusId> & id_grid_to_solver,
                                      const RealVect & bus_vn_kv,
                                      real_type sn_mva,
                                      bool ac) {};
        virtual void _deactivate(int el_id, SolverControl & solver_control) {
            if(status_[el_id]){
                solver_control.tell_ybus_some_coeffs_zero();
                solver_control.tell_recompute_ybus();
                solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
            }
        };
        virtual void _reactivate(int el_id, SolverControl & solver_control) {
            if(!status_[el_id]){
                solver_control.tell_recompute_ybus();
                solver_control.tell_ybus_change_sparsity_pattern();
                solver_control.tell_dimension_changed();  // if the extremity of the line is alone on a bus, this can happen...
            }
        };
        virtual void _change_bus(int el_id, GridModelBusId new_bus_id, SolverControl & solver_control, int nb_bus) {
            GridModelBusId & bus_me_id = bus_id_(el_id);
            
            if(bus_me_id != new_bus_id) {
                // TODO speed: here the dimension changed only if nothing was connected before
                solver_control.tell_dimension_changed();  // in this case i changed the bus, i need to recompute the jacobian and reset the solver
                
                // TODO speed: sparsity pattern might not change if something is already there  
                solver_control.tell_ybus_change_sparsity_pattern();
                solver_control.tell_recompute_ybus();  // if a bus changed for shunts / line / trafo
            }
        };

    protected:
        // physical properties

        // data for grid2op compat

        // input data

};

#endif  //ONE_SIDE_CONTAINER_FORBRANCH_H