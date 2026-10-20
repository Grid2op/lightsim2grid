// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef POWERFLOW_ADJOINT_H
#define POWERFLOW_ADJOINT_H

#include <vector>

#include "Utils.hpp"
#include "powerflow_algorithm/AdjointSolveInfo.hpp"

namespace ls2g {

class LSGrid;
struct SolverBusLayout;
class VoltageControlPlan;

/**
The part of reverse-mode differentiation that is about ONE converged solve, shared by a
batch (BaseBatchSweep, once per row) and a single grid (LSGrid::solve_JT).

Once lambda solves J^T lambda = xbar, an injection's gradient is lambda at its row: it is
a right-hand side. What is left is the inputs that are not, and the rows that are not
equations of the physics any more -- and that is a property of the solve, described by
its AdjointSolveInfo (BaseAlgo::get_adjoint_solve_info), not of who ran it:

- the identity rows, whose multiplier is dropped;
- gen_v, a set-point fixing a bus magnitude (directly, or through the identity row of a
  bus pinned PV), or a voltage-control group's v_set;
- the gen_p of a unit frozen at a limit of its capability curve, whose reactive output
  then follows its active set-point, and of one the distributed slack holds at a bound,
  which its set-point no longer moves.
**/
namespace power_flow_adjoint {

/**
Which part of the gen_v gradient each generator has, the same for every solve of a grid
labelling.

- `target_bus`: the solver bus whose magnitude its set-point fixes, -1 where it fixes
  none -- a generator disconnected, not regulating, regulating a bus a voltage-control
  group holds (see `vc_group`), or one whose bus an SVC or an hvdc converter station
  holds (their set-point no input moves, so the generator's is not free either). A bus
  with a magnitude unknown is a target too: whether its set-point reaches the solve is
  then per solve (pinned PV, or switched PQ), see share_in_solve.
- `vc_group`: the voltage-control group whose v_set its set-point is, -1 otherwise.
- `share`: how much of its bus' (or group's) derivative it carries -- 1/n where n tie
  -- see BaseBatchSweep::get_gen_v_share for why a share and not a derivative each.
- `vm_buses`: the distinct target buses, what BaseAlgo::get_adjoint_solve_info is asked
  the dF/d|V| column of.
**/
struct GenVLayout
{
    IntVect target_bus;
    IntVect vc_group;
    RealVect share;
    std::vector<int> vm_buses;
};

GenVLayout gen_v_layout(const LSGrid & grid, const SolverBusLayout & layout, const VoltageControlPlan & plan);

/** The Jacobian rows of one solve that are identities: both rows of each bus of `masked`
(null: none), and the Q row of each bus the solve pinned PV. Appended to `out`. */
void identity_rows(const AdjointSolveInfo & info, const std::vector<int> * masked,
                   const IntVect & p_row, const IntVect & q_row, std::vector<int> & out);

/** Whether generator `g`'s set-point reaches this solve, as a share of the derivative:
`layout.share[g]` if so, 0 where it reaches nothing (a target bus switched PQ, a group
whose voltage row lost its set-point). `vm_col` is the solve's bus -> magnitude column. */
real_type share_in_solve(const GenVLayout & layout, const AdjointSolveInfo & info,
                         const IntVect & vm_col, Eigen::Index g);

/** The indirect half of one solve's gen_v gradient, -lambda^T dF/dv per generator, already
weighted by share_in_solve (the caller weights the direct half the same way). `group_v_row`
is the J row of each group's voltage constraint. `out` has one entry per generator. */
void gen_v_indirect(const GenVLayout & layout, const AdjointSolveInfo & info,
                    const IntVect & vm_col, const IntVect & group_v_row,
                    const real_type * lambda, real_type * out);

/** What a solve adds to the gen_p gradient beyond lambda at the P row (AdjointSolveInfo::
gen_p_terms): the units whose reactive output follows their target P, the ones whose target
reaches nothing. Added to `out`, one entry per generator. */
void gen_p_indirect(const AdjointSolveInfo & info, const real_type * lambda, real_type * out, Eigen::Index nb_gen);

/** A solver-bus-keyed map re-keyed on grid buses, -1 where the bus is not solved. */
IntVect bus_map_to_grid(const SolverBusLayout & layout, const IntVect & solver_keyed);

}  // namespace power_flow_adjoint

}  // namespace ls2g

#endif  // POWERFLOW_ADJOINT_H
