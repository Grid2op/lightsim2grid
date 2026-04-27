// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SCALING_POLICIES_H
#define SCALING_POLICIES_H

namespace ls2g {

/**
 * Runtime step-scaling strategies for the Newton-Raphson loop.
 *
 *   NoScaling       — full Newton step (alpha = 1), zero overhead.
 *   MaxVoltageChange — clamp step so max|dVa| <= max_dVa and max|dVm| <= max_dVm.
 *   LineSearch       — Armijo backtracking line search.
 *   Iwamoto          — Iwamoto optimal multiplier step.
 */
enum class ScalingPolicyType : int {
    NoScaling        = 0,
    MaxVoltageChange = 1,
    LineSearch       = 2,
    Iwamoto          = 3
};

} // namespace ls2g

#endif // SCALING_POLICIES_H
