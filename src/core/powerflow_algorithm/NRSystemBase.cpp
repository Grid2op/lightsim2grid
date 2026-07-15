// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

// out-of-line on purpose: NRSystem.hpp only forward-declares LSGrid, the full
// type is needed to query the slack buses that need a free Vm unknown + Q
// equation. Base is the one component present in EVERY NRSystem instantiation
// (single- or multi-slack), so this must live here rather than in MultiSlack:
// a grid with a single slack bus that is not locally voltage-pinned (a PQ
// distributed-slack-style participant, or a remote-voltage / SVC-controlled
// slack) needs the exact same free Vm unknown + Q equation whether solved
// with NR_KLU (MultiSlack present) or NRSing_KLU (no MultiSlack extension).
#include "LSGrid.hpp"

namespace ls2g {

void Base::update_state(
    const LSGrid *                         lsgrid_ptr,
    const Eigen::SparseMatrix<cplx_type>&  /*Ybus*/,
    const Eigen::Ref<const CplxVect> &              /*Sbus*/,
    const Eigen::Ref<const RealVect> &             /*slack_weights*/
)
{
    // Slack buses not pinned by a LOCAL voltage-regulating generator need a
    // free Vm unknown + Q equation (added in register_in), exactly like an
    // ordinary PQ bus. See LSGrid::get_free_vm_slack_solver_buses for the
    // exact criterion (GeneratorContainer::gen_is_local_voltage_controller).
    free_vm_slack_buses_.clear();
    if (lsgrid_ptr != nullptr)
        free_vm_slack_buses_ = lsgrid_ptr->get_free_vm_slack_solver_buses();
}

} // namespace ls2g
