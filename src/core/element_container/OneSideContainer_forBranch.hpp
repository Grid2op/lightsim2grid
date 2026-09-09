// Copyright (c) 2024-2026, RTE (https://www.rte-france.com)
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

namespace ls2g {


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
                OneSideForBranchInfo(const OneSideContainer_ForBranch & r_data_pq, int my_id) noexcept:
                OneSideInfo(r_data_pq, my_id) {}
        };
    
    /////////////////////////////
    // iterator
    private:
        using OSCC4BonstIterator = GenericContainerConstIterator<OneSideContainer_ForBranch>;

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
        OneSideContainer_ForBranch() noexcept = default;
        explicit OneSideContainer_ForBranch(bool /*is_trafo*/) noexcept{};
        ~OneSideContainer_ForBranch() noexcept override = default;

        // public generic API

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            OneSideContainer::StateRes
            > ;

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
            const Eigen::Ref<const RealVect> & /*els_p*/,
            const Eigen::Ref<const RealVect> & /*els_q*/,
            const Eigen::Ref<const Eigen::VectorXi> & els_bus_id,
            const std::string & /*name_el*/
            )  // osc: one side element
        {
            init_osc(els_bus_id);
        }

    protected:
        // A branch end lives in Ybus, not in Sbus: opening / closing / moving it
        // changes the matrix (and, since the end may be alone on its bus, possibly
        // the bus set: `tell_one_el_changed_bus`). The sparsity pattern only ever
        // grows on a reconnection or a move; an opening leaves the coefficients in
        // place (some become 0) -- see AlgoControl for what each flag costs.
        void _on_deactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_ybus_some_coeffs_zero();
            solver_control.tell_recompute_ybus();
            solver_control.tell_one_el_changed_bus();  // if the extremity of the line is alone on a bus, this can happen...
        }
        void _on_reactivate(int /*el_id*/, DualAlgoControl & solver_control) override {
            solver_control.tell_recompute_ybus();
            solver_control.tell_ybus_change_sparsity_pattern();
            solver_control.tell_one_el_changed_bus();  // if the extremity of the line is alone on a bus, this can happen...
        }
        void _on_change_bus(int /*el_id*/, GridModelBusId /*new_bus_id*/, DualAlgoControl & solver_control) override {
            // TODO speed: here the dimension changed only if nothing was connected before
            solver_control.tell_one_el_changed_bus();  // in this case i changed the bus, i need to recompute the jacobian and reset the solver
            // TODO speed: sparsity pattern might not change if something is already there
            solver_control.tell_ybus_change_sparsity_pattern();
            solver_control.tell_recompute_ybus();  // if a bus changed for shunts / line / trafo
        }

    protected:
        // physical properties

        // data for grid2op compat

        // input data

};


} // namespace ls2g

#endif  //ONE_SIDE_CONTAINER_FORBRANCH_H
