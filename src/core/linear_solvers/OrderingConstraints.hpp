// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef ORDERING_CONSTRAINTS_H
#define ORDERING_CONSTRAINTS_H

#include <vector>

namespace ls2g {

/**
What a Newton system knows of its Jacobian that a fill-reducing ordering does not: some
indices should be eliminated after some others.

The case it exists for: a row the outer loops switch by value between two forms, each with
its single entry on a different column -- a voltage row, "Vm(reg) = target" or "Q_c = q"; a
ratio-group row, "Vm(bus) = target" or "ratio = held". A pivot order chosen in one form puts
the pivot on an entry the other form zeroes, and every switch then costs a refactorization
failure and a fresh factorization. When the Vm columns of the buses involved are eliminated
first, their pivots on the buses' Q rows (never switched) fill the row's other entry in, and
one pivot order serves both forms (see docs/dev_notes/speed_optimisation.md).

An index is a row and a column at once: the ordering is symmetric, and these rows sit on the
diagonal of the column they pivot on. A lead index is never a follow one (the leads are Vm
columns, the follows controller unknowns and custom rows), so that every constraint can be
honoured off one AMD order (see constrained_amd_order).
**/
struct OrderingConstraint
{
    std::vector<int> follow;  // eliminated after...
    std::vector<int> lead;    // ...every one of these
};

using OrderingConstraints = std::vector<OrderingConstraint>;

namespace detail {
    // a linear solver that takes ordering constraints gets the Newton system's, which are only
    // built then; any other one is left as it is (C++14 detection, as for the rgrowth threshold)
    template<class Solver, class System>
    auto set_ordering_constraints_impl(Solver & solver, const System & system, int)
        -> decltype(solver.set_ordering_constraints(system.ordering_constraints()), void()) {
        solver.set_ordering_constraints(system.ordering_constraints());
    }
    template<class Solver, class System>
    void set_ordering_constraints_impl(Solver &, const System &, long) {}
}  // namespace detail

} // namespace ls2g

#endif // ORDERING_CONSTRAINTS_H
