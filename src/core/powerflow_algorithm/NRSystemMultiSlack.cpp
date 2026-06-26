// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

// out-of-line on purpose: NRSystem.hpp only forward-declares LSGrid (it is
// included BY LSGrid.hpp), the full type is needed to query the slack buses that
// need a free Vm unknown + Q equation.
#include "LSGrid.hpp"

namespace ls2g {

void MultiSlack::update_state(
    const Base *                           /*nr_system_base_ptr*/,
    const LSGrid *                         lsgrid_ptr,
    const Eigen::SparseMatrix<cplx_type>&  /*Ybus*/,
    const CplxVect&                        Sbus,
    const RealVect&                        slack_weights
)
{
    slack_weights_ = slack_weights;
    // initial slack absorbed (see MultiSlackPolicy::initial_slack_absorbed)
    slack_absorbed_ = std::real(Sbus.sum());
    // slack buses not pinned by a local PV generator need a free Vm + Q equation
    // (added in register_in): PQ distributed-slack participants, and remote-voltage
    // / SVC controllers (to which the VoltageControl extension then attaches).
    free_vm_slack_buses_.clear();
    if(lsgrid_ptr != nullptr)
        free_vm_slack_buses_ = lsgrid_ptr->get_free_vm_slack_solver_buses();
}

} // namespace ls2g
