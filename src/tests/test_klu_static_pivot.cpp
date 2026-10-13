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
// then factorizes again. Given ordering constraints, its analysis orders with AMD then moves
// each constrained index after the ones it must follow, which is what keeps a pivot for a row
// switched by value between two forms. C++14 only.

#include <algorithm>
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
#include "linear_solvers/OrderingConstraints.hpp"
#endif

#ifdef KLU_SOLVER_AVAILABLE

using ls2g::ErrorType;
using ls2g::KLULinearSolver;
using ls2g::KLUStaticPivotLinearSolver;
using ls2g::LinearSolverPolicy;
using ls2g::OrderingConstraint;
using ls2g::OrderingConstraints;
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

namespace {

std::vector<int> positions(const std::vector<int> & perm)
{
    std::vector<int> pos(perm.size(), -1);
    for(std::size_t k = 0; k < perm.size(); ++k) pos[perm[k]] = static_cast<int>(k);
    return pos;
}

// A ring of n "bus" indices, plus index n: the Q unknown of a controller at bus `ctrl` (an
// entry in that bus' row) and the voltage row switched by value between its two forms,
// "Vm(ctrl) = target" (an entry on column `ctrl`, the one on its own column a stored zero)
// and "Q = q" (the other way round). Both forms share one pattern.
SpMat with_switched_row(int n, int ctrl, bool voltage_form)
{
    const SpMat ring = with_values(n, jacobian_like_pattern(n), 7);
    std::vector<Eigen::Triplet<real_type> > triplets;
    for(int col = 0; col < ring.outerSize(); ++col){
        for(SpMat::InnerIterator it(ring, col); it; ++it) triplets.emplace_back(it.row(), it.col(), it.value());
    }
    triplets.emplace_back(ctrl, n, -1.);
    triplets.emplace_back(n, ctrl, voltage_form ? 1. : 0.);
    triplets.emplace_back(n, n, voltage_form ? 0. : 1.);
    SpMat res(n + 1, n + 1);
    res.setFromTriplets(triplets.begin(), triplets.end());
    res.makeCompressed();
    return res;
}

}  // namespace

TEST_CASE("a constrained ordering is AMD's with the follow indices moved after their leads", "[linear_solver][klu][static_pivot][ordering]")
{
    const SpMat A = with_values(N, jacobian_like_pattern(N), 6);
    std::vector<int> amd_perm(N);
    double control[AMD_CONTROL];
    double info[AMD_INFO];
    ls2g::amd_defaults(control);
    REQUIRE(ls2g::amd_order(N, A.outerIndexPtr(), A.innerIndexPtr(), amd_perm.data(), control, info) == AMD_OK);

    // no constraint: AMD's order
    std::vector<int> perm;
    REQUIRE(ls2g::constrained_amd_order(N, A.outerIndexPtr(), A.innerIndexPtr(), OrderingConstraints(), perm));
    CHECK(perm == amd_perm);

    // the indices AMD places first and last, so that some constraints have to move things
    OrderingConstraints constraints;
    constraints.push_back(OrderingConstraint{{amd_perm[0], amd_perm[1]}, {amd_perm[N - 1], amd_perm[N / 2]}});
    constraints.push_back(OrderingConstraint{{amd_perm[2]}, {amd_perm[3]}});
    constraints.push_back(OrderingConstraint{{amd_perm[N - 2]}, {amd_perm[4]}});  // already after
    REQUIRE(ls2g::constrained_amd_order(N, A.outerIndexPtr(), A.innerIndexPtr(), constraints, perm));

    // a permutation, each follow index after every one of its leads
    std::vector<int> sorted = perm;
    std::sort(sorted.begin(), sorted.end());
    for(int i = 0; i < N; ++i) REQUIRE(sorted[i] == i);
    const std::vector<int> pos = positions(perm);
    for(const OrderingConstraint & c : constraints){
        for(int f : c.follow){
            for(int l : c.lead) CHECK(pos[f] > pos[l]);
        }
    }
    // a moved index lands right after its last lead
    CHECK(pos[amd_perm[0]] == pos[amd_perm[N - 1]] + 1);
    CHECK(pos[amd_perm[1]] == pos[amd_perm[N - 1]] + 2);
    CHECK(pos[amd_perm[2]] == pos[amd_perm[3]] + 1);
    // and everything else keeps AMD's relative order (the one already after its lead included)
    std::vector<int> rest_amd, rest;
    const std::set<int> moved = {amd_perm[0], amd_perm[1], amd_perm[2]};
    for(int i : amd_perm) if(!moved.count(i)) rest_amd.push_back(i);
    for(int i : perm) if(!moved.count(i)) rest.push_back(i);
    CHECK(rest == rest_amd);
}

TEST_CASE("ordered after its leads, a row switched by value keeps a pivot in both its forms", "[linear_solver][klu][static_pivot][ordering]")
{
    const int ctrl = 5;
    const SpMat q_form = with_switched_row(N, ctrl, false);
    const SpMat v_form = with_switched_row(N, ctrl, true);
    REQUIRE(q_form.nonZeros() == v_form.nonZeros());
    const RealVect b = rhs(N + 1);

    // AMD alone eliminates the switched index first (it has one neighbour): factorized in
    // its Q form, it pivots on its own diagonal, which the voltage form zeroes
    {
        KLUStaticPivotLinearSolver solver;
        REQUIRE(solver.analyze(q_form) == ErrorType::NoError);
        REQUIRE(solver.factorize(q_form) == ErrorType::NoError);
        CHECK(solver.refactorize(v_form) == ErrorType::SolverReFactor);
    }

    // after the Vm column of the bus, whose pivot fills its diagonal in: both forms go through,
    // whichever the factorization saw
    OrderingConstraints constraints;
    constraints.push_back(OrderingConstraint{{N}, {ctrl}});
    for(int first_q = 0; first_q < 2; ++first_q){
        const SpMat & first = first_q ? q_form : v_form;
        const SpMat & other = first_q ? v_form : q_form;
        KLUStaticPivotLinearSolver solver;
        solver.set_ordering_constraints(constraints);
        REQUIRE(solver.analyze(first) == ErrorType::NoError);
        REQUIRE(solver.factorize(first) == ErrorType::NoError);
        for(int round = 0; round < 2; ++round){
            for(const SpMat * A : {&other, &first}){
                REQUIRE(solver.refactorize(*A) == ErrorType::NoError);
                RealVect x = b;
                REQUIRE(solver.solve(x) == ErrorType::NoError);
                CHECK(residual(*A, x, b) < 1e-11);
            }
        }
    }

    // reset forgets them
    KLUStaticPivotLinearSolver solver;
    solver.set_ordering_constraints(constraints);
    REQUIRE(solver.reset() == ErrorType::NoError);
    CHECK(solver.ordering_constraints().empty());
}

TEST_CASE("the policy forwards ordering constraints to a solver that takes them", "[linear_solver][klu][static_pivot][ordering]")
{
    const int ctrl = 5;
    LinearSolverPolicy<KLUStaticPivotLinearSolver> solver;
    OrderingConstraints constraints;
    constraints.push_back(OrderingConstraint{{N}, {ctrl}});
    solver.set_ordering_constraints(constraints);
    const SpMat q_form = with_switched_row(N, ctrl, false);
    REQUIRE(solver.analyze(q_form) == ErrorType::NoError);
    REQUIRE(solver.factorize(q_form) == ErrorType::NoError);
    REQUIRE(solver.refactorize(with_switched_row(N, ctrl, true)) == ErrorType::NoError);
    CHECK(solver.get_linear_solver_stats().nb_refactorize_failed == 0);
}

#endif  // KLU_SOLVER_AVAILABLE
