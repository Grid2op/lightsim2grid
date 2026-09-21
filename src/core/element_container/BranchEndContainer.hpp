// Copyright (c) 2024-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BRANCH_END_CONTAINER_H
#define BRANCH_END_CONTAINER_H


#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "Utils.hpp"
#include "OneSideContainer.hpp"

namespace ls2g {


/**
 * One END of a line or a transformer: the side type of BranchContainer.
 *
 * A branch end stands in Ybus, not in Sbus, so what it adds to OneSideContainer
 * is the set of flags its own connection changes raise (the `_on_xxx` hooks
 * below). It never counts for the per-bus element counts itself: the branch's
 * global status gates whether an end holds its bus, so TwoSidesContainer does
 * the counting around both ends -- which is why only the `*_no_bus_tracking`
 * mutators are ever called on it.
 *
 * Its StateRes wraps OneSideContainer's in a one-element tuple: that nesting is
 * what every pickle and binary file of a line or transformer carries, so it
 * stays even though the class adds no state of its own.
 */
class BranchEndContainer : public OneSideContainer
{
    public:
        BranchEndContainer() noexcept = default;
        ~BranchEndContainer() noexcept override = default;

        // public generic API

        // /!\ if you change this layout, bump BINARY_FORMAT_VERSION (BinaryArchive.hpp)

        using StateRes = std::tuple<
            OneSideContainer::StateRes
            > ;
        enum StateResIdx {
            OSC_STATE = 0,
            NB_ELEM
        };
        static_assert(std::tuple_size<StateRes>::value == StateResIdx::NB_ELEM,
                      "BranchEndContainer::StateRes and StateResIdx do not match");

        StateRes get_state() const
        {
            return get_osc_forB_state();
        }

        void set_state(StateRes & state)
        {
            set_osc_forB_state(state);
        }

    protected:
        BranchEndContainer::StateRes get_osc_forB_state() const  // osc: one side element
        {
            BranchEndContainer::StateRes res(
                get_osc_state());
            return res;
        }

        void set_osc_forB_state(BranchEndContainer::StateRes & my_state)  // osc: one side element
        {
            // read data from my_state
            set_osc_state(std::get<StateResIdx::OSC_STATE>(my_state));
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

};


} // namespace ls2g

#endif  //BRANCH_END_CONTAINER_H
