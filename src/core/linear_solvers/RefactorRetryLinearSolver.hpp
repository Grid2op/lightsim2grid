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
Same counting/timing behavior as LinearSolverPolicy<LinearSolver>, plus: if refactorize()
fails, fall back to a full factorize() (reusing whatever symbolic state the underlying
solver already holds) before giving up. This is a defensive measure recommended by
SuiteSparse's own docs for KLU's klu_refactor/klu_factor pair, generalized here to any
LinearSolver exposing a real factorize/refactorize distinction (KLU, CKTSO -- for
SparseLU/NICSLU, factorize() and refactorize() are the same call, so the fallback is a
harmless no-op there).

`final`: matches the project's existing convention (e.g. KLULinearSolver, NRAlgo) --
this is meant to be used only as a concrete, leaf LinearSolver type, never derived from
further. It does NOT enable any virtual-dispatch optimization (there is nothing virtual
here to devirtualize): refactorize() below hides (does not override) the non-virtual
base method of the same name, resolved entirely at compile time since this class is only
ever used as a template parameter, never accessed through a LinearSolverPolicy<>&/* base
reference.
**/
template<class LinearSolver>
class RefactorRetryLinearSolver final : public LinearSolverPolicy<LinearSolver>
{
    public:
        ErrorType refactorize(const EigenRefConstRealSpMat & J) {
            ++this->stats_.nb_refactorize;
            auto timer = CustTimer();
            ErrorType res = this->inner_.refactorize(J);
            this->stats_.timer_refactor_ += timer.duration();
            if (res != ErrorType::NoError) {
                ++this->stats_.nb_refactorize_failed;
                ++this->stats_.nb_fallback_factorize;
                auto timer_f = CustTimer();
                res = this->inner_.factorize(J);
                this->stats_.timer_factor_ += timer_f.duration();
                ++this->stats_.nb_factorize;
                if (res != ErrorType::NoError) ++this->stats_.nb_fallback_factorize_failed;
            }
            return res;
        }
};

} // namespace ls2g

#endif // REFACTOR_RETRY_LINEAR_SOLVER_H
