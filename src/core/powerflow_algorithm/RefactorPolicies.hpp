// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef REFACTOR_POLICIES_H
#define REFACTOR_POLICIES_H

namespace ls2g {

/**
 * Runtime refactorization strategies for the Newton-Raphson loop.
 *
 *   AlwaysRefactor — rebuild and refactorize J every iteration (default).
 *   EveryN         — refactorize every N iterations; update values only in between
 *                    (chord method variant).
 *   Chord          — build J once on the first iteration, reuse factorization for
 *                    all subsequent iterations (pure chord / Shamanskii with N=inf).
 */
enum class RefactorPolicyType : int {
    AlwaysRefactor = 0,
    EveryN         = 1,
    Chord          = 2
};

} // namespace ls2g

#endif // REFACTOR_POLICIES_H
