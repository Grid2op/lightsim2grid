// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NRSystem.hpp"

namespace ls2g {

void MultiSlack::update_state(
    const Base *                           /*nr_system_base_ptr*/,
    const LSGrid *                         /*lsgrid_ptr*/,
    const Eigen::SparseMatrix<cplx_type>&  /*Ybus*/,
    const CplxVect&                        Sbus,
    Eigen::Ref<const RealVect>             slack_weights
)
{
    slack_weights_ = slack_weights;
    // initial slack absorbed (see MultiSlackPolicy::initial_slack_absorbed)
    slack_absorbed_ = std::real(Sbus.sum());
}

} // namespace ls2g
