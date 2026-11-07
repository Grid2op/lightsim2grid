// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ADJOINT_SOLVE_INFO_H
#define ADJOINT_SOLVE_INFO_H

#include <vector>

#include "Utils.hpp"

namespace ls2g {

/**
What the adjoint of one solve needs to know about it beyond its Jacobian (see
BaseAlgo::get_adjoint_solve_info).

The gradient of a loss through a converged powerflow is read off the system of its
last Newton solve: J^T lambda = xbar, then one product with dF/du per input. The
Jacobian alone does not say which of its rows are not equations of the physics any
more, nor how an input that is not a right-hand side enters a row. That is this:

- `pinned_buses`: the solver buses whose Q row the solve set to the identity (a bus
  still PV that was reserved a magnitude unknown, see NRSystem::set_pv_pinned_buses).
  Their multiplier belongs to no reactive balance and is dropped.
- `vset_lost_groups`: the voltage-control groups whose voltage row no longer holds
  their set-point (every controller held at an output, or all of them stranded): a
  generator of such a group has no gen_v gradient there.
- `vm_col_*`: for each bus asked for, dF/d|V_b| by Jacobian row -- the column J holds
  for a bus whose magnitude is an unknown, and the one it would hold for a bus whose
  magnitude is fixed (a PV bus, the slack). It is what a set-point fixing |V_b| moves,
  evaluated on the solve's own admittance matrix (taps and sections the loops settled
  included). Compressed: bus `vm_col_bus[k]` owns the entries
  `[vm_col_start[k], vm_col_start[k + 1])` of `vm_col_row` / `vm_col_val`.
- `gen_p_terms`: what the gen_p gradient owes beyond lambda at the unit's P row (pu
  per MW): a unit frozen at a limit of its capability curve, whose reactive output
  follows its active set-point (`row` the row that output enters, `coeff`
  d(its specified injection)/d(target P)); a unit whose set-point reaches nothing (the
  distributed slack holds it at a bound), which takes its P row back off (`coeff`
  -1 / sn_mva).
**/
struct AdjointSolveInfo
{
    struct GenPTerm {
        int gen_id;
        int row;           // the Jacobian row the term reads lambda at
        real_type coeff;   // pu per MW
    };

    std::vector<int> pinned_buses;
    std::vector<int> vset_lost_groups;

    std::vector<int> vm_col_bus;
    std::vector<int> vm_col_start;     // one more entry than vm_col_bus
    std::vector<int> vm_col_row;
    std::vector<real_type> vm_col_val;

    std::vector<GenPTerm> gen_p_terms;

    void clear() {
        pinned_buses.clear();
        vset_lost_groups.clear();
        vm_col_bus.clear();
        vm_col_start.clear();
        vm_col_row.clear();
        vm_col_val.clear();
        gen_p_terms.clear();
    }
};

}  // namespace ls2g

#endif  // ADJOINT_SOLVE_INFO_H
