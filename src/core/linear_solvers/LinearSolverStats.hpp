// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LINEAR_SOLVER_STATS_H
#define LINEAR_SOLVER_STATS_H

#include <cmath>
#include <cstddef>
#include <limits>

namespace ls2g {

/**
Per-call counters and timings for a linear solver, as tracked by LinearSolverPolicy<>
(and its RefactorRetryLinearSolver<> derivative). All algorithms wrapping a LinearSolver
(NRAlgo, BaseDCAlgo, BaseFDPFAlgo) source their own timer_factor_/timer_refactor_/
timer_initialize_/timer_solve_ reporting from an instance of this struct instead of
keeping their own separate CustTimer-based bookkeeping.

Field names for the timers deliberately match TimerJac's (see BaseAlgo.hpp) so callers
can copy them across directly.
**/
struct LinearSolverStats {
    std::size_t nb_reset = 0;
    std::size_t nb_analyze = 0;
    std::size_t nb_factorize = 0;                  // includes fallback-triggered factorize() calls
    std::size_t nb_refactorize = 0;                // refactorize() attempts
    std::size_t nb_refactorize_failed = 0;         // of those, how many failed
    std::size_t nb_fallback_factorize = 0;         // subset of nb_factorize triggered by a failed refactorize
    std::size_t nb_fallback_factorize_failed = 0;  // of those, how many also failed

    std::size_t nb_solve = 0;
    std::size_t nb_solve_transpose = 0;            // solve_transpose() calls (the adjoint / VJP path)

    double timer_initialize_ = 0.;  // time spent in analyze()
    double timer_factor_     = 0.;  // time spent in factorize() (incl. fallback factors)
    double timer_refactor_   = 0.;  // time spent in refactorize() itself (not the fallback factor)
    double timer_solve_      = 0.;  // time spent in solve()
    double timer_solve_transpose_ = 0.;  // time spent in solve_transpose()

    // the smallest reciprocal pivot growth of the refactorizations accepted since the timers
    // were last reset (so, like them, the last solve's); NaN when none was measured. Only
    // measured with the refactor fallback on, by a solver that can tell (KLU, KLU static pivot)
    double min_rgrowth_ = std::numeric_limits<double>::quiet_NaN();
};

namespace detail {
    // min_rgrowth_ of two stats taken together: NaN (nothing measured) only if both are
    inline void merge_min_rgrowth(double & into, double from) noexcept {
        if (!std::isnan(from) && (std::isnan(into) || from < into)) into = from;
    }

    // The field lists of LinearSolverStats, kept next to it: a field added to the struct is
    // added to the one of these it belongs to.

    // what the factorizations cost, summed: analyses, factorizations, refactorizations, their
    // failures and fallbacks, their timers, and the smallest pivot growth of the two
    inline void add_factorization_stats(LinearSolverStats & into, const LinearSolverStats & from) noexcept {
        into.nb_analyze += from.nb_analyze;
        into.nb_factorize += from.nb_factorize;
        into.nb_refactorize += from.nb_refactorize;
        into.nb_refactorize_failed += from.nb_refactorize_failed;
        into.nb_fallback_factorize += from.nb_fallback_factorize;
        into.nb_fallback_factorize_failed += from.nb_fallback_factorize_failed;
        into.timer_initialize_ += from.timer_initialize_;
        into.timer_factor_ += from.timer_factor_;
        into.timer_refactor_ += from.timer_refactor_;
        merge_min_rgrowth(into.min_rgrowth_, from.min_rgrowth_);
    }
    // a forward solver's own work, summed: its resets and its solve() calls
    inline void add_solve_stats(LinearSolverStats & into, const LinearSolverStats & from) noexcept {
        into.nb_reset += from.nb_reset;
        into.nb_solve += from.nb_solve;
        into.timer_solve_ += from.timer_solve_;
    }
    // an adjoint solver's own work, summed: its solve_transpose() calls
    inline void add_solve_transpose_stats(LinearSolverStats & into, const LinearSolverStats & from) noexcept {
        into.nb_solve_transpose += from.nb_solve_transpose;
        into.timer_solve_transpose_ += from.timer_solve_transpose_;
    }
    // `now` with every counter replaced by how much it grew since `before` (the counters cover
    // a solver's lifetime); the timers and min_rgrowth_ are `now`'s, which cover its last solve
    inline LinearSolverStats counters_since(const LinearSolverStats & now, const LinearSolverStats & before) noexcept {
        LinearSolverStats res = now;
        res.nb_reset -= before.nb_reset;
        res.nb_analyze -= before.nb_analyze;
        res.nb_factorize -= before.nb_factorize;
        res.nb_refactorize -= before.nb_refactorize;
        res.nb_refactorize_failed -= before.nb_refactorize_failed;
        res.nb_fallback_factorize -= before.nb_fallback_factorize;
        res.nb_fallback_factorize_failed -= before.nb_fallback_factorize_failed;
        res.nb_solve -= before.nb_solve;
        res.nb_solve_transpose -= before.nb_solve_transpose;
        return res;
    }

    // C++14-safe detection idiom (deliberately not `if constexpr`, which is C++17: this
    // project falls back to cxx_std_14 on older compilers, see CMakeLists.txt). Prefers
    // the real accessor when the wrapped LinearSolver provides one (LinearSolverPolicy<>
    // and derivatives), else returns an all-zero struct for plain, unwrapped solvers.
    template<class T>
    auto get_stats_impl(const T& solver, int) -> decltype(solver.get_linear_solver_stats()) {
        return solver.get_linear_solver_stats();
    }
    template<class T>
    LinearSolverStats get_stats_impl(const T&, long) {
        return LinearSolverStats{};
    }

    // Same idiom for resetting only the *timer* fields (see LinearSolverPolicy::
    // reset_stats_timers): a no-op for plain, unwrapped solvers. Counters (nb_*) are
    // deliberately left alone -- they accumulate over the algorithm's whole lifetime,
    // not per powerflow call, so a systematic (vs. occasional) fallback shows up even
    // if nobody inspects stats after every single solve.
    template<class T>
    auto reset_stats_timers_impl(T& solver, int) -> decltype(solver.reset_stats_timers(), void()) {
        solver.reset_stats_timers();
    }
    template<class T>
    void reset_stats_timers_impl(T&, long) {}
} // namespace detail

} // namespace ls2g

#endif // LINEAR_SOLVER_STATS_H
