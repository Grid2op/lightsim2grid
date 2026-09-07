// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

// out-of-line on purpose: NRSystem.hpp only forward-declares LSGrid (it is
// included BY LSGrid.hpp), the full type is needed to pull the controller data
#include "LSGrid.hpp"

namespace ls2g {

void VoltageControl::update_state(
    const Base                       * /*nr_system_base_ptr*/,
    const LSGrid                     * lsgrid_ptr,
    const EigenRefConstCplxSpMat     & /*Ybus*/,
    const Eigen::Ref<const CplxVect> & /*Sbus*/,
    const Eigen::Ref<const RealVect> & /*slack_weights*/
)
{
    // READ, not re-derived: layer 3 of the plan the grid built into its AC cache
    // during pre_process_solver (see LSGrid::_build_into_cache). Building it here
    // meant walking every generator, SVC and converter station of the grid a second
    // time per solve -- plus the free-Vm slack pass a third -- for an answer that
    // cannot have changed since. It also means the controller list and the pv-pq
    // split it is keyed on are now built from one another rather than from two
    // independent walks of the containers.
    data_.clear();
    if(lsgrid_ptr != nullptr) data_ = lsgrid_ptr->get_ac_voltage_control_plan().controllers();
    my_size_ = data_.n_controllers();
    // per-solve init: the reactive injection state starts at 0 (gen convention)
    q_ = RealVect::Zero(my_size_);
}

} // namespace ls2g
