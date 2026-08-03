// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef LINEAR_SOLVER_STATS_H
#define LINEAR_SOLVER_STATS_H

#include <cstddef>

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

    double timer_initialize_ = 0.;  // time spent in analyze()
    double timer_factor_     = 0.;  // time spent in factorize() (incl. fallback factors)
    double timer_refactor_   = 0.;  // time spent in refactorize() itself (not the fallback factor)
    double timer_solve_      = 0.;  // time spent in solve()
};

namespace detail {
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
