// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// solve_transpose() answers J^T x = b out of the factorization of J itself -- no
// transposed copy, no second analyze, no second numeric factorization. The reference
// here is a dense LU of the explicitly transposed matrix, and every matrix below is
// deliberately asymmetric in BOTH its pattern and its values: a solve_transpose that
// quietly solved J x = b instead would return a different vector and fail.

#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Eigen/Dense"

#include "Utils.hpp"
#include "linear_solvers/LinearSolverPolicy.hpp"
#include "linear_solvers/SparseLUSolver.hpp"
#ifdef KLU_SOLVER_AVAILABLE
#include "linear_solvers/KLUSolver.hpp"
#endif

using Catch::Approx;
using ls2g::ErrorType;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

using SpMat = Eigen::SparseMatrix<real_type>;

// A deterministic, strongly asymmetric, diagonally dominant matrix: an entry (i, j)
// exists without (j, i) existing, so J^T does not even share J's pattern. `shift`
// perturbs the values without touching the pattern -- that is the batch-adjoint case
// (one analyze, then refactorize per row).
SpMat make_asym_matrix(int n, real_type shift = 0.)
{
    std::vector<Eigen::Triplet<real_type> > coeffs;
    for(int i = 0; i < n; ++i){
        coeffs.push_back(Eigen::Triplet<real_type>(i, i, 10. + i + shift));
        // one sub-diagonal, two super-diagonals: no (j, i) mirror for most of them
        if(i + 1 < n) coeffs.push_back(Eigen::Triplet<real_type>(i, i + 1, -1. - 0.25 * i + shift));
        if(i + 3 < n) coeffs.push_back(Eigen::Triplet<real_type>(i, i + 3, 0.5 + 0.125 * i));
        if(i - 2 >= 0) coeffs.push_back(Eigen::Triplet<real_type>(i, i - 2, 2. + 0.5 * i - shift));
    }
    // a couple of far-off entries, on one side only, to break any residual symmetry
    coeffs.push_back(Eigen::Triplet<real_type>(0, n - 1, 3.5));
    coeffs.push_back(Eigen::Triplet<real_type>(n - 2, 1, -2.25));

    SpMat res(n, n);
    res.setFromTriplets(coeffs.begin(), coeffs.end());
    res.makeCompressed();
    return res;
}

RealVect make_rhs(int n)
{
    RealVect b(n);
    for(int i = 0; i < n; ++i) b(i) = 1. + 0.75 * i - (i % 3);
    return b;
}

// x such that J^T x = b, by a dense LU of the explicitly transposed matrix
RealVect dense_transposed_solve(const SpMat & J, const RealVect & b)
{
    const Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> dense = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>(J).transpose();
    return dense.fullPivLu().solve(b);
}

void check_close(const RealVect & got, const RealVect & expected, real_type tol = 1e-10)
{
    REQUIRE(got.size() == expected.size());
    for(Eigen::Index i = 0; i < got.size(); ++i){
        REQUIRE(got(i) == Approx(expected(i)).margin(tol));
    }
}

}  // namespace


TEST_CASE("SparseLU solves the transposed system out of the factorization of J")
{
    const int n = 12;
    const SpMat J = make_asym_matrix(n);
    const RealVect b = make_rhs(n);

    ls2g::SparseLULinearSolver solver;
    REQUIRE(solver.analyze(J) == ErrorType::NoError);
    REQUIRE(solver.factorize(J) == ErrorType::NoError);

    RealVect x = b;
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    check_close(x, dense_transposed_solve(J, b));

    // ... and the residual of the transposed system is what it claims to be
    const RealVect residual = SpMat(J.transpose()) * x - b;
    REQUIRE(residual.lpNorm<Eigen::Infinity>() == Approx(0.).margin(1e-10));

    // the untransposed solve still solves J x = b, and gives a DIFFERENT vector --
    // without this the test above would pass on a solve_transpose that ignores the
    // transposition entirely
    RealVect y = b;
    REQUIRE(solver.solve(y) == ErrorType::NoError);
    REQUIRE((y - x).lpNorm<Eigen::Infinity>() > 1e-3);
    REQUIRE((J * y - b).lpNorm<Eigen::Infinity>() == Approx(0.).margin(1e-10));
}


TEST_CASE("a transposed solve after a refactorize uses the new values")
{
    const int n = 12;
    const SpMat J = make_asym_matrix(n);
    const SpMat J2 = make_asym_matrix(n, 0.75);   // same pattern, different values
    REQUIRE(J.nonZeros() == J2.nonZeros());
    const RealVect b = make_rhs(n);

    ls2g::SparseLULinearSolver solver;
    REQUIRE(solver.analyze(J) == ErrorType::NoError);
    REQUIRE(solver.factorize(J) == ErrorType::NoError);
    REQUIRE(solver.refactorize(J2) == ErrorType::NoError);

    RealVect x = b;
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    check_close(x, dense_transposed_solve(J2, b));
}


TEST_CASE("LinearSolverPolicy times and counts the transposed solves separately")
{
    const int n = 12;
    const SpMat J = make_asym_matrix(n);
    const RealVect b = make_rhs(n);

    ls2g::LinearSolverPolicy<ls2g::SparseLULinearSolver> solver;
    // the flags are constexpr: asserted at compile time, and deliberately NOT through
    // REQUIRE, whose expression decomposition binds a reference to its operands -- that
    // would odr-use a constexpr static member the headers no longer define out of line
    // (needed before C++17, and this project still builds at C++14).
    static_assert(ls2g::LinearSolverPolicy<ls2g::SparseLULinearSolver>::CAN_SOLVE_TRANSPOSE ==
                  ls2g::SparseLULinearSolver::CAN_SOLVE_TRANSPOSE,
                  "the policy must mirror the capability of the solver it wraps");

    REQUIRE(solver.analyze(J) == ErrorType::NoError);
    REQUIRE(solver.factorize(J) == ErrorType::NoError);

    RealVect x = b;
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    check_close(x, dense_transposed_solve(J, b));

    RealVect y = b;
    REQUIRE(solver.solve(y) == ErrorType::NoError);

    const ls2g::LinearSolverStats & stats = solver.get_linear_solver_stats();
    REQUIRE(stats.nb_solve_transpose == 1);
    REQUIRE(stats.nb_solve == 1);   // a transposed solve is not counted as a plain one

    solver.reset_stats_timers();
    REQUIRE(stats.timer_solve_transpose_ == 0.);
}


#ifdef KLU_SOLVER_AVAILABLE
TEST_CASE("KLU solves the transposed system out of the factorization of J")
{
    const int n = 12;
    const SpMat J = make_asym_matrix(n);
    const RealVect b = make_rhs(n);

    static_assert(ls2g::KLULinearSolver::CAN_SOLVE_TRANSPOSE, "KLU has klu_tsolve");

    ls2g::KLULinearSolver solver;
    REQUIRE(solver.analyze(J) == ErrorType::NoError);
    REQUIRE(solver.factorize(J) == ErrorType::NoError);

    RealVect x = b;
    REQUIRE(solver.solve_transpose(x) == ErrorType::NoError);
    check_close(x, dense_transposed_solve(J, b));

    RealVect y = b;
    REQUIRE(solver.solve(y) == ErrorType::NoError);
    REQUIRE((y - x).lpNorm<Eigen::Infinity>() > 1e-3);
}


TEST_CASE("KLU and SparseLU agree on the transposed solve, before and after a refactorize")
{
    const int n = 20;
    const RealVect b = make_rhs(n);

    ls2g::KLULinearSolver klu;
    ls2g::SparseLULinearSolver splu;

    const SpMat J = make_asym_matrix(n);
    REQUIRE(klu.analyze(J) == ErrorType::NoError);
    REQUIRE(klu.factorize(J) == ErrorType::NoError);
    REQUIRE(splu.analyze(J) == ErrorType::NoError);
    REQUIRE(splu.factorize(J) == ErrorType::NoError);

    RealVect x_klu = b;
    RealVect x_splu = b;
    REQUIRE(klu.solve_transpose(x_klu) == ErrorType::NoError);
    REQUIRE(splu.solve_transpose(x_splu) == ErrorType::NoError);
    check_close(x_klu, x_splu);

    // the adjoint of a batch does exactly this: one analyze, then per row new values
    // into the same pattern, refactorize, transposed solve
    const SpMat J2 = make_asym_matrix(n, -1.5);
    REQUIRE(klu.refactorize(J2) == ErrorType::NoError);
    REQUIRE(splu.refactorize(J2) == ErrorType::NoError);

    x_klu = b;
    x_splu = b;
    REQUIRE(klu.solve_transpose(x_klu) == ErrorType::NoError);
    REQUIRE(splu.solve_transpose(x_splu) == ErrorType::NoError);
    check_close(x_klu, x_splu);
    check_close(x_klu, dense_transposed_solve(J2, b));
}
#endif  // KLU_SOLVER_AVAILABLE
