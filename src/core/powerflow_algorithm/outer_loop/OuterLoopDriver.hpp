// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_LOOP_DRIVER_H
#define OUTER_LOOP_DRIVER_H

#include <vector>

#include "Utils.hpp"

namespace ls2g {

class LSGrid;

/**
 * Parameters of the outer-loop driver itself (not of any loop), with OpenLoadFlow's names.
 * The defaults are those of the pypowsybl build the comparisons are made with, see
 * docs/dev_notes/outer_loops_fixed_sparsity.md.
 */
struct OuterLoopDriverParams
{
    int max_outer_iterations = 30;                         ///< maxOuterLoopIterations
    bool voltage_remote_control_robust_mode = true;        ///< voltageRemoteControlRobustMode
    real_type min_realistic_voltage = 0.8;                 ///< minRealisticVoltage, pu
    real_type max_realistic_voltage = 1.2;                 ///< maxRealisticVoltage, pu
    real_type min_nominal_voltage_realistic_check = 180.;  ///< minNominalVoltageRealisticVoltageCheck, kV
};

/**
 * OpenLoadFlow's isStateUnrealistic: whether a bus whose voltage magnitude is an unknown of
 * the Newton (`vm_unknown[bus]`, solver numbering) and whose nominal voltage is at least
 * `min_nominal_voltage_realistic_check` lies outside the realistic band. Buses held at a
 * set-point are not looked at, as in OpenLoadFlow (only its BUS_V variables are).
 */
LS2G_API bool is_state_unrealistic(const LSGrid & grid,
                                   const CplxVect & V,
                                   const std::vector<bool> & vm_unknown,
                                   const OuterLoopDriverParams & params);

}  // namespace ls2g

#endif  // OUTER_LOOP_DRIVER_H
