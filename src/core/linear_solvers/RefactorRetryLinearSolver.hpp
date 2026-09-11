// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef REFACTOR_RETRY_LINEAR_SOLVER_H
#define REFACTOR_RETRY_LINEAR_SOLVER_H

#include "Utils.hpp"
#include "CustTimer.hpp"
#include "LinearSolverStats.hpp"
#include "LinearSolverPolicy.hpp"

namespace ls2g {

/**
LinearSolverPolicy<LinearSolver> with its refactorize-failure fallback switched on from
construction (see LinearSolverPolicy::set_refactor_fallback): if refactorize() fails, a
full numeric factorize() (reusing the symbolic analysis the underlying solver already
holds) is tried before giving up. A defensive measure recommended by SuiteSparse's own
docs for KLU's klu_refactor/klu_factor pair, generalized to any LinearSolver exposing a
real factorize/refactorize distinction (KLU, CKTSO -- for SparseLU/NICSLU, factorize() and
refactorize() are the same call, so the fallback is a harmless no-op there).

The fallback itself lives in the policy, behind a switch, so that an algorithm which is
NOT one of the NRRefactorRetry_* family can still turn it on when it knows its own edits
may move a pivot -- the batch algorithms do, when they mask buses (BaseAlgo::
set_refactor_fallback). This class only sets the default: the NRRefactorRetry_* names
are "the fallback, always on".

`final`: matches the project's existing convention (e.g. KLULinearSolver, NRAlgo) --
this is meant to be used only as a concrete, leaf LinearSolver type, never derived from
further.
**/
template<class LinearSolver>
class RefactorRetryLinearSolver final : public LinearSolverPolicy<LinearSolver>
{
    public:
        RefactorRetryLinearSolver() noexcept { this->set_refactor_fallback(true); }
};

} // namespace ls2g

#endif // REFACTOR_RETRY_LINEAR_SOLVER_H
