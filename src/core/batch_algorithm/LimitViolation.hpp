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
    GRID = 3,  // the whole grid / contingency, not a specific element (see LimitViolationType::DIVERGENCE)
    GENERATOR = 4  // see LimitViolationType::LOW_Q / HIGH_Q and compute_gen_q_violations
};

// the kind of limit that was violated
enum class LS2G_API LimitViolationType : int {
    LOW_VOLTAGE = 0,
    HIGH_VOLTAGE = 1,
    CURRENT = 2,
    NOT_SIMULATED = 3,  // a pre-check (graph connectivity) skipped the contingency: the solver was never invoked
    DIVERGENCE = 4,  // the solver was invoked (for the contingency, or the pre-contingency "n" case) but did not converge
    // a voltage-regulating generator's reactive output left [min_q_mvar, max_q_mvar]: what
    // OpenLoadFlow's `ReactiveLimits` outer loop would act on. Reported (see
    // compute_gen_q_violations), never enforced -- lightsim2grid does not switch the bus
    // PV -> PQ and re-solve.
    LOW_Q = 5,
    HIGH_Q = 6
};

// a single limit violation, as detected by ContingencyAnalysis (see compute_limit_violations)
struct LS2G_API LimitViolation {
    ViolationElementType element_type;
    // grid-model bus id for BUS ; local (0-based, own type) line / trafo / generator id
    // otherwise ; unused (-1) for GRID
    int element_id;
    int side;  // 1 or 2 for LINE / TRAFO ; unused (0) for BUS / GENERATOR / GRID
    LimitViolationType violation_type;
    // value reached (MVAr for LOW_Q / HIGH_Q) ; unused (NaN) for NOT_SIMULATED / DIVERGENCE
    real_type value;
    // limit that was violated (min_q_mvar for LOW_Q, max_q_mvar for HIGH_Q) ; unused (NaN)
    // for NOT_SIMULATED / DIVERGENCE
    real_type limit;
    // element name: LINE / TRAFO / GENERATOR (from LSGrid::set_line_names /
    // set_trafo_names / set_gen_names) or, for BUS, the name of the *substation* the
    // violating bus belongs to (from LSGrid::set_substation_names) -- there is no per-bus
    // name in LSGrid, only per-substation ones. Empty string if the grid never had names
    // set for the relevant kind, or for GRID.
    std::string name{};
};

} // namespace ls2g

#endif  // LIMITVIOLATION_H
