// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef OUTER_CONTROL_LAYOUT_H
#define OUTER_CONTROL_LAYOUT_H

#include <vector>

#include "Utils.hpp"

namespace ls2g {

/**
 * Where the continuous controls of the outer loops sit in the Newton's system, read-only: the
 * transformers of BranchControl (PhaseControl, TransformerVoltageControl) and the shunt
 * controllers of ShuntControl (ShuntVoltageControl), with the rows and columns they own in J.
 * For whoever reproduces that system elsewhere (another solver of the same equations): the
 * structure is the one of the last rebuild, so it holds for every solve until the next one.
 *
 * Solver bus ids, J row / column indices (-1 where not owned). The values (shift, ratio,
 * susceptance, which form a row holds) are per solve and not here; the taps' blocks follow
 * from the grid (TrafoContainer::pi_block_at, shift_sign, ratio_exponent).
 */
struct OuterControlLayout {
    /// a transformer BranchControl handles
    struct Branch {
        int trafo = -1;              ///< grid id
        bool phase = false;          ///< a phase shifter (its tap may move)
        bool shift_column = false;   ///< ... solving for its shift: a column and a row
        bool ratio = false;          ///< in a ratio group (its ratio tap may move)
        bool ratio_column = false;   ///< ... of a solved one: a ratio column
        int bus1_solver = -1;
        int bus2_solver = -1;
        int shift_col = -1;          ///< the shift unknown
        int shift_row = -1;          ///< its row (alpha = target, or P = target)
        int ratio_col = -1;          ///< the ratio unknown
    };
    /// the transformers of a solved ratio group (regulating one bus)
    struct RatioGroup {
        int bus_solver = -1;         ///< the regulated bus
        real_type target_vm = 1.;    ///< pu
        std::vector<int> members;    ///< indices into `branches`, in the group's order
        int vm_col = -1;             ///< the regulated bus' Vm column
        std::vector<int> rows;       ///< one per member, same order
    };
    /// a shunt controller of ShuntControl: the regulating shunts of one bus
    struct Shunt {
        int bus_solver = -1;
        std::vector<int> shunts;     ///< grid ids
        bool column = false;         ///< of a solved group: a susceptance column
        int col = -1;                ///< the susceptance unknown
        int q_row = -1;              ///< its bus' Q row (dQ/dB)
    };
    /// the shunt controllers of a solved group (regulating one bus)
    struct ShuntGroup {
        int bus_solver = -1;
        real_type target_vm = 1.;
        std::vector<int> members;    ///< indices into `shunts`
        int vm_col = -1;
        std::vector<int> rows;
    };

    std::vector<Branch> branches;
    std::vector<RatioGroup> ratio_groups;
    std::vector<Shunt> shunts;
    std::vector<ShuntGroup> shunt_groups;
};

}  // namespace ls2g

#endif  // OUTER_CONTROL_LAYOUT_H
