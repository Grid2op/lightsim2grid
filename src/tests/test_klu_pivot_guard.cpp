// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// A KLU refactorization keeps the pivot order of the last factorization and only reports
// an exactly-zero pivot. A pivot that became tiny but not zero gives a garbage step without
// any error. The pivot-growth guard (KLULinearSolver::set_pivot_growth_guard, switched on by
// the algorithm for a Jacobian with voltage-control rows, effective with the refactor
// fallback) turns that case into a failed refactorization, which the fallback answers with a
// numeric factorize: the pivots are chosen again, the symbolic analysis is kept.

#ifdef KLU_SOLVER_AVAILABLE

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Eigen/Dense"

#include "Utils.hpp"
#include "linear_solvers/LinearSolverPolicy.hpp"
#include "linear_solvers/KLUSolver.hpp"

using Catch::Approx;
using ls2g::ErrorType;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

using SpMat = Eigen::SparseMatrix<real_type>;
using Solver = ls2g::LinearSolverPolicy<ls2g::KLULinearSolver>;

// full 2x2 pattern; KLU keeps the diagonal pivots of `healthy` on a refactorization
SpMat make_matrix(real_type a00)
{
    SpMat res(2, 2);
    res.insert(0, 0) = a00;
    res.insert(0, 1) = 1.;
    res.insert(1, 0) = 1.;
    res.insert(1, 1) = 3.;
    res.makeCompressed();
    return res;
}

RealVect solve_with(Solver & solver, const SpMat & J, const RealVect & b)
{
    RealVect x = b;
    REQUIRE(solver.solve(x) == ErrorType::NoError);
    // the answer is judged against J itself, not against another solver
    return J * x - b;
}

}  // namespace

TEST_CASE("a refactorization with a collapsed pivot is re-pivoted, not trusted", "[klu][refactorize]")
{
    const SpMat healthy = make_matrix(4.);
    const SpMat collapsed = make_matrix(1e-14);   // same pattern, the (0, 0) pivot is tiny but not zero
    const RealVect b = (RealVect(2) << 1., 2.).finished();

    SECTION("without the guard KLU accepts it, as before")
    {
        Solver solver;
        REQUIRE(solver.analyze(healthy) == ErrorType::NoError);
        REQUIRE(solver.factorize(healthy) == ErrorType::NoError);
        REQUIRE(solver.refactorize(collapsed) == ErrorType::NoError);
        REQUIRE(solver.get_linear_solver_stats().nb_refactorize_failed == 0);
    }

    SECTION("with the guard and the fallback it is re-pivoted: one analyze, an exact answer")
    {
        Solver solver;
        solver.set_refactor_fallback(true);
        solver.set_pivot_growth_guard(true);
        REQUIRE(solver.analyze(healthy) == ErrorType::NoError);
        REQUIRE(solver.factorize(healthy) == ErrorType::NoError);
        REQUIRE(solver.refactorize(collapsed) == ErrorType::NoError);
        const auto stats = solver.get_linear_solver_stats();
        CHECK(stats.nb_analyze == 1);
        CHECK(stats.nb_refactorize_failed == 1);
        CHECK(stats.nb_fallback_factorize == 1);
        CHECK(solve_with(solver, collapsed, b).norm() == Approx(0.).margin(1e-9));
    }

    SECTION("a healthy change of values is refactorized, nothing falls back")
    {
        Solver solver;
        solver.set_refactor_fallback(true);
        solver.set_pivot_growth_guard(true);
        REQUIRE(solver.analyze(healthy) == ErrorType::NoError);
        REQUIRE(solver.factorize(healthy) == ErrorType::NoError);
        const SpMat moved = make_matrix(5.);
        REQUIRE(solver.refactorize(moved) == ErrorType::NoError);
        CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 0);
        CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 0);
        CHECK(solve_with(solver, moved, b).norm() == Approx(0.).margin(1e-12));
    }

    SECTION("the guard asked for without the fallback does nothing: nothing would catch it")
    {
        Solver solver;
        solver.set_pivot_growth_guard(true);
        REQUIRE(solver.analyze(healthy) == ErrorType::NoError);
        REQUIRE(solver.factorize(healthy) == ErrorType::NoError);
        REQUIRE(solver.refactorize(collapsed) == ErrorType::NoError);
        CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 0);
    }
}

#endif  // KLU_SOLVER_AVAILABLE
