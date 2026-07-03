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

namespace ls2g {

// the kind of element on which a limit was violated
enum class ViolationElementType : int {
    BUS = 0,
    LINE = 1,
    TRAFO = 2
};

// the kind of limit that was violated
enum class LimitViolationType : int {
    LOW_VOLTAGE = 0,
    HIGH_VOLTAGE = 1,
    CURRENT = 2
};

// a single limit violation, as detected by ContingencyAnalysis (see compute_limit_violations)
struct LimitViolation {
    ViolationElementType element_type;
    int element_id;  // grid-model bus id for BUS ; local (0-based, own type) line / trafo id otherwise
    int side;  // 1 or 2 for LINE / TRAFO ; unused (0) for BUS
    LimitViolationType violation_type;
    real_type value;  // value reached
    real_type limit;  // limit that was violated
};

} // namespace ls2g

#endif  // LIMITVIOLATION_H
