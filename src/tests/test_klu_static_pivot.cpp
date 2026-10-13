// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// KLUStaticPivotLinearSolver: KLU's factorization, then refactorizations on its frozen pivot
// order and nonzero pattern. They must give the solution klu_refactor gives (to round-off),
// solve the transposed system out of the same factors, and refuse a frozen pivot that turned
// zero (always) or tiny (above an rgrowth threshold), which LinearSolverPolicy's fallback
// then factorizes again. C++14 only.

#include <cmath>
#include <limits>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "Utils.hpp"
#include "linear_solvers/LinearSolverPolicy.hpp"
#ifdef KLU_SOLVER_AVAILABLE
#include "linear_solvers/KLUSolver.hpp"
#include "linear_solvers/KLUStaticPivotSolver.hpp"
#endif

#ifdef KLU_SOLVER_AVAILABLE

using ls2g::ErrorType;
using ls2g::KLULinearSolver;
using ls2g::KLUStaticPivotLinearSolver;
using ls2g::LinearSolverPolicy;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

using SpMat = Eigen::SparseMatrix<real_type>;

SpMat make2(real_type a00, real_type a01, real_type a10, real_type a11)
{
    SpMat res(2, 2);
    res.insert(0, 0) = a00;
    res.insert(0, 1) = a01;
    res.insert(1, 0) = a10;
    res.insert(1, 1) = a11;
    res.makeCompressed();
    return res;
}

// The pattern of a power-flow Jacobian: structurally symmetric, a diagonal, a ring of
// neighbours and a few long-range links. Unsymmetric values, not diagonally dominant, so
// that KLU's partial pivoting has choices to make.
std::vector<std::pair<int, int> > jacobian_like_pattern(int n)
{
    std::set<std::pair<int, int> > entries;
    for(int i = 0; i < n; ++i){
        entries.insert({i, i});
        const int j = (i + 1) % n;
        entries.insert({i, j});
        entries.insert({j, i});
    }
    for(int i = 0; i + 7 < n; i += 5){
        entries.insert({i, i + 7});
        entries.insert({i + 7, i});
    }
    return std::vector<std::pair<int, int> >(entries.begin(), entries.end());
}

SpMat with_values(int n, const std::vector<std::pair<int, int> > & pattern, unsigned seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<real_type> val(-1., 1.);
    std::vector<Eigen::Triplet<real_type> > triplets;
    for(const auto & e : pattern){
        real_type v = val(gen);
        if(e.first == e.second) v += 1.5;
        triplets.emplace_back(e.first, e.second, v);
    }
    SpMat res(n, n);
    res.setFromTriplets(triplets.begin(), triplets.end());
    res.makeCompressed();
    return res;
}

RealVect rhs(int n)
{
    RealVect b(n);
    for(int i = 0; i < n; ++i) b(i) = 1. + 0.1 * i;
    return b;
}

real_type residual(const SpMat & A, const RealVect & x, const RealVect & b)
{
    return (A * x - b).cwiseAbs().maxCoeff();
}

const int N = 40;

}  // namespace

TEST_CASE("a static-pivot refactorization solves what klu_refactor solves", "[linear_solver][klu][static_pivot]")
{
    const auto pattern = jacobian_like_pattern(N);
    const SpMat first = with_values(N, pattern, 1);

    KLUStaticPivotLinearSolver solver;
    KLULinearSolver klu;
    REQUIRE(solver.analyze(first) == ErrorType::NoError);
    REQUIRE(solver.factorize(first) == ErrorType::NoError);
    REQUIRE(klu.analyze(first) == ErrorType::NoError);
    REQUIRE(klu.factorize(first) == ErrorType::NoError);

    // right after the factorization: KLU's own factors
    RealVect x = rhs(N);
    REQUIRE(solver.solve(x) == ErrorType::NoError);
    CHECK(residual(first, x, rhs(N)) < 1e-12);

    // several refactorizations in a row, on the same pattern and close values (a Newton)
    for(unsigned k = 0; k < 5; ++k){
        SpMat next = first;
        const SpMat delta = with_values(N, pattern, 100 + k);
        for(int p = 0; p < next.nonZeros(); ++p) next.valuePtr()[p] += 0.05 * delta.valuePtr()[p];

        REQUIRE(solver.refactorize(next) == ErrorType::NoError);
        REQUIRE(klu.refactorize(next) == ErrorType::NoError);
        RealVect x_sp = rhs(N);
        RealVect x_klu = rhs(N);
        REQUIRE(solver.solve(x_sp) == ErrorType::NoError);
        REQUIRE(klu.solve(x_klu) == ErrorType::NoError);
        CHECK(residual(next, x_sp, rhs(N)) < 1e-12);
        CHECK((x_sp - x_klu).cwiseAbs().maxCoeff() < 1e-12);

        // and against a dense LU with partial pivoting
        const Eigen::MatrixXd dense(next);
        const RealVect x_dense = dense.partialPivLu().solve(rhs(N));
        CHECK((x_sp - x_dense).cwiseAbs().maxCoeff() < 1e-12);
    }
}

TEST_CASE("the transposed system is solved out of the same factors", "[linear_solver][klu][static_pivot]")
{
    const auto pattern = jacobian_like_pattern(N);
    const SpMat first = with_values(N, pattern, 2);
    const SpMat next = with_values(N, pattern, 3);
    const SpMat next_t = SpMat(next.transpose());

    KLUStaticPivotLinearSolver solver;
    REQUIRE(solver.analyze(first) == ErrorType::NoError);
    REQUIRE(solver.factorize(first) == ErrorType::NoError);

    // on KLU's factors
    RealVect x = rhs(N);
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    CHECK(residual(SpMat(first.transpose()), x, rhs(N)) < 1e-12);

    // on the solver's own
    REQUIRE(solver.refactorize(next) == ErrorType::NoError);
    x = rhs(N);
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    CHECK(residual(next_t, x, rhs(N)) < 1e-11);

    KLULinearSolver klu;
    REQUIRE(klu.analyze(first) == ErrorType::NoError);
    REQUIRE(klu.factorize(first) == ErrorType::NoError);
    REQUIRE(klu.refactorize(next) == ErrorType::NoError);
    RealVect x_klu = rhs(N);
    REQUIRE(klu.solve_transpose(x_klu) == ErrorType::NoError);
    CHECK((x - x_klu).cwiseAbs().maxCoeff() < 1e-11);
}

TEST_CASE("a zero frozen pivot is always rejected", "[linear_solver][klu][static_pivot]")
{
    // whatever the symmetric ordering, the first pivot is a diagonal entry
    KLUStaticPivotLinearSolver solver;
    REQUIRE(solver.rgrowth_threshold() == 0.);
    REQUIRE(solver.analyze(make2(2., 1., 1., 2.)) == ErrorType::NoError);
    REQUIRE(solver.factorize(make2(2., 1., 1., 2.)) == ErrorType::NoError);
    CHECK(solver.refactorize(make2(0., 1., 1., 0.)) == ErrorType::SolverReFactor);
    // nothing to solve with after a rejection
    RealVect b(2);
    b << 1., 2.;
    CHECK(solver.solve(b) == ErrorType::SolverSolve);

    // the work vector was left clean: a good refactorization still goes through
    REQUIRE(solver.refactorize(make2(3., 1., 1., 3.)) == ErrorType::NoError);
    b << 4., 4.;
    REQUIRE(solver.solve(b) == ErrorType::NoError);
    CHECK(std::abs(b(0) - 1.) < 1e-14);
    CHECK(std::abs(b(1) - 1.) < 1e-14);

    // a non-finite pivot too
    const real_type nan = std::numeric_limits<real_type>::quiet_NaN();
    CHECK(solver.refactorize(make2(nan, 1., 1., nan)) == ErrorType::SolverReFactor);
}

TEST_CASE("a tiny frozen pivot is rejected above an rgrowth threshold", "[linear_solver][klu][static_pivot]")
{
    const SpMat first = make2(2., 1., 1., 2.);
    const SpMat tiny = make2(1e-14, 1., 1., 1e-14);

    KLUStaticPivotLinearSolver solver;
    REQUIRE(solver.analyze(first) == ErrorType::NoError);
    REQUIRE(solver.factorize(first) == ErrorType::NoError);
    // no threshold: only a zero pivot is refused, as klu_refactor
    CHECK(solver.refactorize(tiny) == ErrorType::NoError);

    solver.set_rgrowth_threshold(1e-10);
    REQUIRE(solver.factorize(first) == ErrorType::NoError);
    REQUIRE(solver.refactorize(make2(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.last_rgrowth() > 0.1);
    CHECK(solver.refactorize(tiny) == ErrorType::SolverReFactor);
    CHECK(solver.last_rgrowth() < 1e-10);

    // KLU's own klu_rgrowth refuses it as well
    KLULinearSolver klu;
    klu.set_rgrowth_threshold(1e-10);
    REQUIRE(klu.analyze(first) == ErrorType::NoError);
    REQUIRE(klu.factorize(first) == ErrorType::NoError);
    CHECK(klu.refactorize(tiny) == ErrorType::SolverReFactor);
}

TEST_CASE("the fallback factorizes a rejected refactorization again", "[linear_solver][klu][static_pivot]")
{
    LinearSolverPolicy<KLUStaticPivotLinearSolver> solver;
    solver.set_refactor_fallback(true);
    const SpMat first = make2(2., 1., 1., 2.);
    REQUIRE(solver.analyze(first) == ErrorType::NoError);
    REQUIRE(solver.factorize(first) == ErrorType::NoError);
    REQUIRE(solver.refactorize(make2(3., 1., 1., 3.)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 0);

    // a zero pivot: new pivots (off the diagonal), and the solution is right
    const SpMat swapped = make2(0., 1., 1., 0.);
    REQUIRE(solver.refactorize(swapped) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 1);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);
    RealVect b(2);
    b << 2., 1.;
    REQUIRE(solver.solve(b) == ErrorType::NoError);
    CHECK(std::abs(b(0) - 1.) < 1e-14);
    CHECK(std::abs(b(1) - 2.) < 1e-14);

    // and the plan follows the new pivots: a refactorization on them goes through
    const SpMat near_swapped = make2(1e-3, 2., 3., 1e-3);
    REQUIRE(solver.refactorize(near_swapped) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_fallback_factorize == 1);
    b << 1e-3 * 1. + 2. * 2., 3. * 1. + 1e-3 * 2.;
    REQUIRE(solver.solve(b) == ErrorType::NoError);
    CHECK(std::abs(b(0) - 1.) < 1e-12);
    CHECK(std::abs(b(1) - 2.) < 1e-12);
}

TEST_CASE("a factorization that keeps the pivots keeps the plan", "[linear_solver][klu][static_pivot]")
{
    const auto pattern = jacobian_like_pattern(N);
    const SpMat first = with_values(N, pattern, 4);
    const SpMat second = with_values(N, pattern, 5);

    KLUStaticPivotLinearSolver solver;
    REQUIRE(solver.analyze(first) == ErrorType::NoError);
    // factorize / refactorize / factorize (same values, so same pivots) / refactorize ...
    for(int round = 0; round < 3; ++round){
        REQUIRE(solver.factorize(first) == ErrorType::NoError);
        REQUIRE(solver.refactorize(second) == ErrorType::NoError);
        RealVect x = rhs(N);
        REQUIRE(solver.solve(x) == ErrorType::NoError);
        CHECK(residual(second, x, rhs(N)) < 1e-11);
    }
    // ... and a factorization of other values, which may pivot elsewhere
    REQUIRE(solver.factorize(second) == ErrorType::NoError);
    REQUIRE(solver.refactorize(first) == ErrorType::NoError);
    RealVect x = rhs(N);
    REQUIRE(solver.solve(x) == ErrorType::NoError);
    CHECK(residual(first, x, rhs(N)) < 1e-11);

    // reset, then a fresh analysis
    REQUIRE(solver.reset() == ErrorType::NoError);
    RealVect y = rhs(N);
    CHECK(solver.solve(y) == ErrorType::SolverSolve);
    REQUIRE(solver.analyze(second) == ErrorType::NoError);
    REQUIRE(solver.factorize(second) == ErrorType::NoError);
    REQUIRE(solver.refactorize(first) == ErrorType::NoError);
    y = rhs(N);
    REQUIRE(solver.solve(y) == ErrorType::NoError);
    CHECK(residual(first, y, rhs(N)) < 1e-11);
}

#endif  // KLU_SOLVER_AVAILABLE
