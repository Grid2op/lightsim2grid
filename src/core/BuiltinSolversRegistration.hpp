// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BUILTIN_SOLVERS_REGISTRATION_H
#define BUILTIN_SOLVERS_REGISTRATION_H

namespace ls2g {

class SolverRegistry;

// Register every in-tree solver into `reg`.
// Call this once from the pybind11 module init before any GridModel is created.
// To add a new built-in solver: add one register_solver(...) call here, guarded
// by the same #ifdef as the solver class itself.
void register_builtin_solvers(SolverRegistry& reg);


} // namespace ls2g

#endif // BUILTIN_SOLVERS_REGISTRATION_H