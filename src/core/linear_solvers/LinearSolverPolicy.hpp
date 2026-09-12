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
        static constexpr bool CAN_SOLVE_MAT = LinearSolver::CAN_SOLVE_MAT;

        // can this linear solver solve J^T x = b using the factorization of J itself,
        // without ever forming J^T (see solve_transpose)
        static constexpr bool CAN_SOLVE_TRANSPOSE = LinearSolver::CAN_SOLVE_TRANSPOSE;

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
            if (res != ErrorType::NoError) {
                ++stats_.nb_refactorize_failed;
                if (refactor_fallback_) {
                    // A refactorize reuses the pivot sequence of the last factorize.
                    // KLU halts on a pivot that is now exactly zero -- which is what a
                    // value-level edit that changes a bus's role does (a masked bus,
                    // a PV bus released to PQ, a stranded controller's row repurposed
                    // into "Q_c = 0"): the entry the base matrix pivoted at is a zero
                    // in the edited one. The symbolic analysis (ordering, pattern) is
                    // still right; only the pivots must be chosen again, and that is
                    // exactly a numeric factorize. SparseLU never gets here: its
                    // refactorize IS a factorize.
                    ++stats_.nb_fallback_factorize;
                    auto timer_f = CustTimer();
                    res = inner_.factorize(J);
                    stats_.timer_factor_ += timer_f.duration();
                    ++stats_.nb_factorize;
                    if (res != ErrorType::NoError) ++stats_.nb_fallback_factorize_failed;
                }
            }
            return res;
        }

        // Whether a failed refactorize() falls back to a numeric factorize() before
        // reporting the failure (see there). Off by default: a plain algorithm
        // reports the failure, so a systematic one stays visible in the stats.
        // RefactorRetryLinearSolver is this policy with the switch on from
        // construction; the batch algorithms turn it on for whatever algorithm they
        // run when they mask buses or switch PV / PQ (BaseAlgo::set_refactor_fallback).
        void set_refactor_fallback(bool val) noexcept { refactor_fallback_ = val; }
        bool refactor_fallback() const noexcept { return refactor_fallback_; }

        ErrorType solve(Eigen::Ref<RealVect> b) {
            ++stats_.nb_solve;
            auto timer = CustTimer();
            ErrorType res = inner_.solve(b);
            stats_.timer_solve_ += timer.duration();
            return res;
        }

        // Solves J^T x = b reusing the factorization of J that factorize() /
        // refactorize() produced -- no transposed copy of J, no second analyze, no
        // second numeric factorization. Only meaningful where CAN_SOLVE_TRANSPOSE is
        // true; the solvers where it is false return ErrorType::NotImplemented rather
        // than silently solving the wrong system, so a caller that cannot check the
        // flag at compile time can check the return value instead.
        ErrorType solve_transpose(Eigen::Ref<RealVect> b) {
            ++stats_.nb_solve_transpose;
            auto timer = CustTimer();
            ErrorType res = inner_.solve_transpose(b);
            stats_.timer_solve_transpose_ += timer.duration();
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
            stats_.timer_solve_transpose_ = 0.;
        }

    protected:
        LinearSolver inner_;
        LinearSolverStats stats_;
        bool refactor_fallback_ = false;

    private:
        // no copy allowed (matches the concrete solver classes' own convention)
        LinearSolverPolicy(const LinearSolverPolicy&) = delete;
        LinearSolverPolicy(LinearSolverPolicy&&) = delete;
        LinearSolverPolicy & operator=(LinearSolverPolicy&&) = delete;
        LinearSolverPolicy & operator=(const LinearSolverPolicy&) = delete;
};

} // namespace ls2g

#endif // LINEAR_SOLVER_POLICY_H
