// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LINEAR_SOLVER_POLICY_H
#define LINEAR_SOLVER_POLICY_H

#include "Utils.hpp"
#include "CustTimer.hpp"
#include "LinearSolverStats.hpp"

namespace ls2g {

/**
Default linear-solver policy: a transparent pass-through wrapper around `LinearSolver`
that counts and times every call (reset/analyze/factorize/refactorize/solve) into a
LinearSolverStats. This is the wrapper every built-in NR/DC/FDPF solver uses by
default (see Solvers.hpp): NRAlgo/BaseDCAlgo/BaseFDPFAlgo no longer keep their own
timer_factor_/timer_refactor_/timer_initialize_ bookkeeping, they read it from here
instead (single source of truth).

Deliberately NOT virtual: `_linear_solver`'s static type is always known at compile time
by the algorithm that holds it (it's a template parameter, never accessed through a base
pointer/reference), so plain method-hiding is correct here and avoids an unnecessary
vtable / indirect-call overhead that a virtual interface would add to every single
solver call for no actual polymorphism ever exercised.
**/
template<class LinearSolver>
class LinearSolverPolicy
{
    public:
        LinearSolverPolicy() noexcept = default;
        ~LinearSolverPolicy() noexcept = default;

        // can this linear solver solve problem where RHS is a matrix
        static const bool CAN_SOLVE_MAT;

        ErrorType reset() {
            ++stats_.nb_reset;
            return inner_.reset();
        }

        ErrorType analyze(const EigenRefConstRealSpMat & J) {
            ++stats_.nb_analyze;
            auto timer = CustTimer();
            ErrorType res = inner_.analyze(J);
            stats_.timer_initialize_ += timer.duration();
            return res;
        }

        ErrorType factorize(const EigenRefConstRealSpMat & J) {
            ++stats_.nb_factorize;
            auto timer = CustTimer();
            ErrorType res = inner_.factorize(J);
            stats_.timer_factor_ += timer.duration();
            return res;
        }

        ErrorType refactorize(const EigenRefConstRealSpMat & J) {
            ++stats_.nb_refactorize;
            auto timer = CustTimer();
            ErrorType res = inner_.refactorize(J);
            stats_.timer_refactor_ += timer.duration();
            return res;
        }

        ErrorType solve(Eigen::Ref<RealVect> b) {
            ++stats_.nb_solve;
            auto timer = CustTimer();
            ErrorType res = inner_.solve(b);
            stats_.timer_solve_ += timer.duration();
            return res;
        }

        const LinearSolverStats & get_linear_solver_stats() const noexcept { return stats_; }

        // Called from the owning algorithm's reset_timer() (itself invoked at the start
        // of every compute_pf/compute_pf_dc): zeroes only the timer_* fields, so
        // get_timers_jacobian() keeps reporting "last call only" like it always has.
        // Counters (nb_*) are untouched -- see detail::reset_stats_timers_impl.
        void reset_stats_timers() noexcept {
            stats_.timer_initialize_ = 0.;
            stats_.timer_factor_     = 0.;
            stats_.timer_refactor_   = 0.;
            stats_.timer_solve_      = 0.;
        }

    protected:
        LinearSolver inner_;
        LinearSolverStats stats_;

    private:
        // no copy allowed (matches the concrete solver classes' own convention)
        LinearSolverPolicy(const LinearSolverPolicy&) = delete;
        LinearSolverPolicy(LinearSolverPolicy&&) = delete;
        LinearSolverPolicy & operator=(LinearSolverPolicy&&) = delete;
        LinearSolverPolicy & operator=(const LinearSolverPolicy&) = delete;
};

template<class LinearSolver>
const bool LinearSolverPolicy<LinearSolver>::CAN_SOLVE_MAT = LinearSolver::CAN_SOLVE_MAT;

} // namespace ls2g

#endif // LINEAR_SOLVER_POLICY_H
