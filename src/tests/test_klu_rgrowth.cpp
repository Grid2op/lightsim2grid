// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// A KLU refactorization reuses the pivots of the last factorization: values that make one
// of them tiny (but not zero) go through klu_refactor, with factors too inaccurate to
// trust. With the refactor fallback on, the reciprocal pivot growth is checked and such a
// refactorization is redone as a factorization (new pivots), as powsybl-math-native does;
// with it off, nothing changes. C++14 only.

#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Utils.hpp"
#include "linear_solvers/LinearSolverPolicy.hpp"
#ifdef KLU_SOLVER_AVAILABLE
#include "linear_solvers/KLUSolver.hpp"
#endif

#ifdef KLU_SOLVER_AVAILABLE

using Catch::Approx;
using ls2g::ErrorType;
using ls2g::KLULinearSolver;
using ls2g::LinearSolverPolicy;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

using SpMat = Eigen::SparseMatrix<real_type>;

SpMat make(real_type a00, real_type a01, real_type a10, real_type a11)
{
    SpMat res(2, 2);
    res.insert(0, 0) = a00;
    res.insert(0, 1) = a01;
    res.insert(1, 0) = a10;
    res.insert(1, 1) = a11;
    res.makeCompressed();
    return res;
}

// a diagonally dominant first matrix: the diagonal is pivoted on
const SpMat FIRST = make(2., 1., 1., 2.);
// the same pattern, the (0, 0) pivot now tiny
const SpMat TINY_PIVOT = make(1e-14, 1., 1., 1.);

}  // namespace

TEST_CASE("fallback off: a tiny reused pivot is not checked", "[linear_solver][klu]")
{
    LinearSolverPolicy<KLULinearSolver> solver;
    REQUIRE(solver.analyze(FIRST) == ErrorType::NoError);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    CHECK(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 0);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 0);
}

TEST_CASE("fallback on: a tiny reused pivot is factorized again", "[linear_solver][klu]")
{
    LinearSolverPolicy<KLULinearSolver> solver;
    solver.set_refactor_fallback(true);
    REQUIRE(solver.analyze(FIRST) == ErrorType::NoError);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);

    // well-conditioned new values: a plain refactorization
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 0);

    REQUIRE(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 1);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize_failed == 0);

    // and the factors are the good ones: x = (1, 2) solves TINY_PIVOT x = b
    RealVect b(2);
    b << 1e-14 + 2., 3.;
    REQUIRE(solver.solve(b) == ErrorType::NoError);
    CHECK(b(0) == Approx(1.).margin(1e-12));
    CHECK(b(1) == Approx(2.).margin(1e-12));

    // switched off again, no check
    solver.set_refactor_fallback(false);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    CHECK(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 1);
}

TEST_CASE("min_rgrowth: the smallest growth of the accepted refactorizations", "[linear_solver][klu]")
{
    LinearSolverPolicy<KLULinearSolver> solver;
    REQUIRE(solver.analyze(FIRST) == ErrorType::NoError);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    // the fallback off: nothing measured
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(std::isnan(solver.get_linear_solver_stats().min_rgrowth_));

    solver.set_refactor_fallback(true);
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    const double first = solver.get_linear_solver_stats().min_rgrowth_;
    CHECK(first > 0.);

    // a refactorization rejected on its growth is not one whose pivots were used
    REQUIRE(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    REQUIRE(solver.get_linear_solver_stats().nb_fallback_factorize == 1);
    CHECK(solver.get_linear_solver_stats().min_rgrowth_ == first);

    // and it covers what the timers cover
    solver.reset_stats_timers();
    CHECK(std::isnan(solver.get_linear_solver_stats().min_rgrowth_));
}

TEST_CASE("the fallback's growth threshold can be set", "[linear_solver][klu]")
{
    LinearSolverPolicy<KLULinearSolver> solver;
    CHECK(solver.fallback_rgrowth_threshold() == LinearSolverPolicy<KLULinearSolver>::RGROWTH_THRESHOLD);
    solver.set_refactor_fallback(true);
    REQUIRE(solver.analyze(FIRST) == ErrorType::NoError);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);

    // above any growth: every refactorization falls back
    solver.set_fallback_rgrowth_threshold(1e300);
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);

    // 0: only what klu_refactor reports, so a tiny pivot goes through
    solver.set_fallback_rgrowth_threshold(0.);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    REQUIRE(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);

    // set while the fallback is off, applied when it is turned on
    solver.set_refactor_fallback(false);
    solver.set_fallback_rgrowth_threshold(1e300);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);
    solver.set_refactor_fallback(true);
    REQUIRE(solver.refactorize(make(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 2);
}

TEST_CASE("the KLU solver alone checks nothing by default", "[linear_solver][klu]")
{
    KLULinearSolver solver;
    CHECK(solver.rgrowth_threshold() == 0.);
    REQUIRE(solver.analyze(FIRST) == ErrorType::NoError);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    CHECK(solver.refactorize(TINY_PIVOT) == ErrorType::NoError);
    solver.set_rgrowth_threshold(1e-10);
    REQUIRE(solver.factorize(FIRST) == ErrorType::NoError);
    CHECK(solver.refactorize(TINY_PIVOT) == ErrorType::SolverReFactor);
}

#endif  // KLU_SOLVER_AVAILABLE
