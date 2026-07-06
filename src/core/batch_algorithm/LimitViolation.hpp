// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LIMITVIOLATION_H
#define LIMITVIOLATION_H

#include "Utils.hpp"
#include <string>

namespace ls2g {

// the kind of element on which a limit was violated
enum class LS2G_API ViolationElementType : int {
    BUS = 0,
    LINE = 1,
    TRAFO = 2,
    GRID = 3  // the whole grid / contingency, not a specific element (see LimitViolationType::DIVERGENCE)
};

// the kind of limit that was violated
enum class LS2G_API LimitViolationType : int {
    LOW_VOLTAGE = 0,
    HIGH_VOLTAGE = 1,
    CURRENT = 2,
    NOT_SIMULATED = 3,  // a pre-check (graph connectivity) skipped the contingency: the solver was never invoked
    DIVERGENCE = 4  // the solver was invoked (for the contingency, or the pre-contingency "n" case) but did not converge
};

// a single limit violation, as detected by ContingencyAnalysis (see compute_limit_violations)
struct LS2G_API LimitViolation {
    ViolationElementType element_type;
    int element_id;  // grid-model bus id for BUS ; local (0-based, own type) line / trafo id otherwise ; unused (-1) for GRID
    int side;  // 1 or 2 for LINE / TRAFO ; unused (0) for BUS / GRID
    LimitViolationType violation_type;
    real_type value;  // value reached ; unused (NaN) for NOT_SIMULATED / DIVERGENCE
    real_type limit;  // limit that was violated ; unused (NaN) for NOT_SIMULATED / DIVERGENCE
    // element name (LINE / TRAFO only, from LSGrid::set_line_names / set_trafo_names) ; empty
    // string if the grid never had names set for this element type, or for BUS / GRID (no
    // per-bus name exists in LSGrid, only per-substation ones)
    std::string name{};
};

} // namespace ls2g

#endif  // LIMITVIOLATION_H
