// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Two levels, because the adjoint can be wrong in two independent ways.
//
// BatchAdjoint on its own is a batch of "same pattern, different values" transposed
// solves: the reference there is a dense LU of each row's explicitly transposed
// matrix, which owes nothing to powerflows.
//
// Wired into a sweep, what matters is a different claim: that lambda really is the
// gradient. The reference there is the sweep itself -- a central finite difference of
// a scalar loss with respect to an injection, which knows nothing about Jacobians. If
// the implicit-function sign were flipped, the ledger read backwards, or the Jacobian
// captured one Newton iterate before the solution, this is what would say so.

#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "Eigen/Dense"

#include "LSGrid.hpp"
#include "Solvers.hpp"
#include "batch_algorithm/BatchAdjoint.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::BatchAdjoint;
using ls2g::CplxVect;
using ls2g::InjectionSweep;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMatRM = BatchAdjoint::RealMatRM;
using SpMat = Eigen::SparseMatrix<real_type>;
using DenseMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic>;

// ============================================================ part A: the solves

// asymmetric in pattern AND values (so J^T is a genuinely different matrix), with a
// `shift` that moves the values without touching the pattern -- one batch row each
SpMat make_asym_matrix(int n, real_type shift)
{
    std::vector<Eigen::Triplet<real_type> > coeffs;
    for(int i = 0; i < n; ++i){
        coeffs.push_back(Eigen::Triplet<real_type>(i, i, 10. + i + shift));
        if(i + 1 < n) coeffs.push_back(Eigen::Triplet<real_type>(i, i + 1, -1. - 0.25 * i + shift));
        if(i + 3 < n) coeffs.push_back(Eigen::Triplet<real_type>(i, i + 3, 0.5 + 0.125 * i));
        if(i - 2 >= 0) coeffs.push_back(Eigen::Triplet<real_type>(i, i - 2, 2. + 0.5 * i - shift));
    }
    coeffs.push_back(Eigen::Triplet<real_type>(0, n - 1, 3.5));
    coeffs.push_back(Eigen::Triplet<real_type>(n - 2, 1, -2.25));
    SpMat res(n, n);
    res.setFromTriplets(coeffs.begin(), coeffs.end());
    res.makeCompressed();
    return res;
}

RealVect dense_transposed_solve(const SpMat & J, const RealVect & b)
{
    const DenseMat dense = DenseMat(J).transpose();
    return dense.fullPivLu().solve(b);
}

real_type row_shift(int row) { return 0.1 * static_cast<real_type>((row * 7) % 13); }

}  // namespace


TEST_CASE("the batch adjoint solves every row against its own Jacobian")
{
    const int n = 14;
    const int nb_rows = 9;

    BatchAdjoint adj;
    adj.allocate(nb_rows, make_asym_matrix(n, 0.));
    REQUIRE(adj.dim_J() == n);
    REQUIRE(adj.nb_rows() == nb_rows);
    REQUIRE(adj.memory_bytes() == static_cast<std::size_t>(nb_rows) * static_cast<std::size_t>(adj.nnz_J()) * sizeof(real_type));

    RealMatRM xbar(nb_rows, n);
    for(int i = 0; i < nb_rows; ++i){
        adj.store_row(i, make_asym_matrix(n, row_shift(i)));
        for(int j = 0; j < n; ++j) xbar(i, j) = 1. + 0.5 * j - (i % 3);
    }

    std::vector<char> row_ok;
    const RealMatRM lambda = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 1, row_ok);

    for(int i = 0; i < nb_rows; ++i){
        REQUIRE(row_ok[static_cast<size_t>(i)] == 1);
        const RealVect expected = dense_transposed_solve(make_asym_matrix(n, row_shift(i)), xbar.row(i).transpose());
        for(int j = 0; j < n; ++j) REQUIRE(lambda(i, j) == Approx(expected(j)).margin(1e-9));
    }

    // ... and a row really is solved against ITS matrix, not against some other row's
    const RealVect wrong = dense_transposed_solve(make_asym_matrix(n, row_shift(0)), xbar.row(3).transpose());
    REQUIRE((lambda.row(3).transpose() - wrong).lpNorm<Eigen::Infinity>() > 1e-6);
}


TEST_CASE("the batch adjoint gives the same answer however many threads it uses")
{
    const int n = 14;
    const int nb_rows = 23;   // deliberately not a multiple of the thread counts below

    BatchAdjoint adj;
    adj.allocate(nb_rows, make_asym_matrix(n, 0.));
    RealMatRM xbar(nb_rows, n);
    for(int i = 0; i < nb_rows; ++i){
        adj.store_row(i, make_asym_matrix(n, row_shift(i)));
        for(int j = 0; j < n; ++j) xbar(i, j) = 0.3 * j - 0.7 * i;
    }

    std::vector<char> ok1, ok4;
    const RealMatRM l1 = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 1, ok1);
    const RealMatRM l4 = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 4, ok4);
    REQUIRE(ok1 == ok4);
    for(int i = 0; i < nb_rows; ++i){
        for(int j = 0; j < n; ++j) REQUIRE(l4(i, j) == Approx(l1(i, j)).margin(1e-10));
    }
}


TEST_CASE("a row the batch never solved has no adjoint, and says so")
{
    const int n = 10;
    const int nb_rows = 5;
    BatchAdjoint adj;
    adj.allocate(nb_rows, make_asym_matrix(n, 0.));
    for(int i = 0; i < nb_rows; ++i){
        if(i == 2) continue;   // as a diverging / skipped row leaves it
        adj.store_row(i, make_asym_matrix(n, row_shift(i)));
    }
    REQUIRE_FALSE(adj.has_row(2));

    const RealMatRM xbar = RealMatRM::Constant(nb_rows, n, 1.);
    std::vector<char> row_ok;
    const RealMatRM lambda = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 1, row_ok);

    REQUIRE(row_ok[2] == 0);
    REQUIRE(lambda.row(2).lpNorm<Eigen::Infinity>() == Approx(0.).margin(0.));
    for(int i = 0; i < nb_rows; ++i){
        if(i == 2) continue;
        REQUIRE(row_ok[static_cast<size_t>(i)] == 1);
        REQUIRE(lambda.row(i).lpNorm<Eigen::Infinity>() > 0.);
    }
}


TEST_CASE("the multiplier of an identity-pinned equation is dropped")
{
    const int n = 12;
    const int nb_rows = 4;
    BatchAdjoint adj;
    adj.allocate(nb_rows, make_asym_matrix(n, 0.));
    for(int i = 0; i < nb_rows; ++i) adj.store_row(i, make_asym_matrix(n, row_shift(i)));

    const RealMatRM xbar = RealMatRM::Constant(nb_rows, n, 1.5);
    std::vector<std::vector<int> > identity_rows(nb_rows);
    identity_rows[1].push_back(4);
    identity_rows[1].push_back(7);
    identity_rows[3].push_back(0);

    std::vector<char> ok_free, ok_pinned;
    const RealMatRM free_ = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 1, ok_free);
    const RealMatRM pinned = adj.solve_JT(xbar, identity_rows, 1, ok_pinned);

    for(int i = 0; i < nb_rows; ++i){
        for(int j = 0; j < n; ++j){
            const bool dropped = (i == 1 && (j == 4 || j == 7)) || (i == 3 && j == 0);
            if(dropped) REQUIRE(pinned(i, j) == Approx(0.).margin(0.));
            else REQUIRE(pinned(i, j) == Approx(free_(i, j)).margin(1e-12));
        }
    }
}


TEST_CASE("the batch adjoint refuses a Jacobian whose pattern is not the batch's")
{
    const int n = 10;
    BatchAdjoint adj;
    adj.allocate(3, make_asym_matrix(n, 0.));

    SpMat other(n, n);                       // a different pattern entirely
    std::vector<Eigen::Triplet<real_type> > coeffs;
    for(int i = 0; i < n; ++i) coeffs.push_back(Eigen::Triplet<real_type>(i, i, 1.));
    other.setFromTriplets(coeffs.begin(), coeffs.end());
    other.makeCompressed();

    REQUIRE_THROWS(adj.store_row(0, other));
    REQUIRE_THROWS(adj.store_row(7, make_asym_matrix(n, 0.)));   // out of range
}


TEST_CASE("several cotangents per row cost solves, not factorizations")
{
    const int n = 11;
    const int nb_rows = 6;
    const int k = 3;
    BatchAdjoint adj;
    adj.allocate(nb_rows, make_asym_matrix(n, 0.));
    for(int i = 0; i < nb_rows; ++i) adj.store_row(i, make_asym_matrix(n, row_shift(i)));

    RealMatRM xbar(nb_rows, k * n);
    for(int i = 0; i < nb_rows; ++i){
        for(int d = 0; d < k; ++d){
            for(int j = 0; j < n; ++j) xbar(i, d * n + j) = 1. + d - 0.3 * j + 0.1 * i;
        }
    }

    std::vector<char> row_ok;
    const RealMatRM lambda = adj.solve_JT(xbar, std::vector<std::vector<int> >(), 1, row_ok);
    REQUIRE(lambda.cols() == k * n);

    for(int i = 0; i < nb_rows; ++i){
        const SpMat J = make_asym_matrix(n, row_shift(i));
        for(int d = 0; d < k; ++d){
            const RealVect expected = dense_transposed_solve(J, xbar.row(i).segment(d * n, n).transpose());
            for(int j = 0; j < n; ++j) REQUIRE(lambda(i, d * n + j) == Approx(expected(j)).margin(1e-9));
        }
    }

    // one analysis and one numeric factorization per row, k triangular solves per row
    const ls2g::LinearSolverStats & stats = adj.get_linear_solver_stats();
    REQUIRE(stats.nb_analyze == 1);
    REQUIRE(stats.nb_solve_transpose == static_cast<std::size_t>(nb_rows * k));
    REQUIRE(stats.nb_factorize + stats.nb_refactorize == static_cast<std::size_t>(nb_rows));
}


// ======================================================= part B: it is the gradient

namespace {

const int NB_BUS = 4;
const int NB_GEN = 2;
const int LOAD_BUS = NB_BUS - 1;
const int GEN_BUS = 1;
const real_type SN_MVA = 100.;

// the 4-bus radial feeder of test_injection_sweep.cpp: slack generator at bus 0, PV
// generator at bus 1, one load at bus 3
LSGrid make_grid()
{
    LSGrid g;
    g.set_sn_mva(SN_MVA);
    g.set_init_vm_pu(1.0);
    g.init_bus(static_cast<unsigned int>(NB_BUS), 1, RealVect::Constant(NB_BUS, 138.), 0, 0);

    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for(int i = 0; i < NB_BUS - 1; ++i){ f(i) = i; t(i) = i + 1; }
    g.init_powerlines(RealVect::Constant(NB_BUS - 1, 0.01), RealVect::Constant(NB_BUS - 1, 0.1),
                      CplxVect::Zero(NB_BUS - 1), f, t);

    RealVect lp(1), lq(1); lp << 50.; lq << 10.;
    Eigen::VectorXi lb(1); lb << LOAD_BUS;
    g.init_loads(lp, lq, lb);

    RealVect p(NB_GEN), v(NB_GEN), q(NB_GEN), qmin(NB_GEN), qmax(NB_GEN);
    Eigen::VectorXi b(NB_GEN);
    p << 0., 20.;
    v << 1.02, 1.01;
    q << 0., 0.;
    qmin << -1000., -1000.;
    qmax << 1000., 1000.;
    b << 0, GEN_BUS;
    g.init_generators_full(p, v, q, std::vector<bool>{true, true}, qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
    return g;
}

const int NB_STEPS = 5;

struct Inputs
{
    RealMatRM gen_p{NB_STEPS, NB_GEN};
    RealMatRM sgen_p{NB_STEPS, 0};
    RealMatRM load_p{NB_STEPS, 1};
    RealMatRM load_q{NB_STEPS, 1};

    Inputs()
    {
        for(int i = 0; i < NB_STEPS; ++i){
            const real_type jitter = static_cast<real_type>((i * 37) % 11) / 10.;
            gen_p(i, 0) = 0.;
            gen_p(i, 1) = 10. + 15. * jitter;
            load_p(i, 0) = 30. + 40. * jitter;
            load_q(i, 0) = 5. + 10. * jitter;
        }
    }
};

// A loss with a genuine dependence on BOTH the angle and the magnitude of every bus
// voltage: L = sum over rows and buses of (wr * Re V + wi * Im V). Deliberately not
// |V|^2, whose angle cotangent is identically zero and would leave half of the ledger
// untested. Its cotangent (torch's convention for a real loss of a complex tensor,
// dL/dRe V + i dL/dIm V) is the constant wr + i wi.
real_type loss_weight_re(int bus) { return 1. + 0.5 * bus; }
real_type loss_weight_im(int bus) { return -0.75 + 0.25 * bus; }

real_type loss_of(const Eigen::Ref<const Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > & V)
{
    real_type res = 0.;
    for(Eigen::Index i = 0; i < V.rows(); ++i){
        for(Eigen::Index b = 0; b < V.cols(); ++b){
            res += loss_weight_re(static_cast<int>(b)) * V(i, b).real() +
                   loss_weight_im(static_cast<int>(b)) * V(i, b).imag();
        }
    }
    return res;
}

real_type run_loss(const Inputs & in)
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    const int status = sweep.compute_Vs(in.gen_p, in.sgen_p, in.load_p, in.load_q,
                                        CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-11);
    REQUIRE(status == 1);
    return loss_of(sweep.get_voltages());
}

}  // namespace


TEST_CASE("the adjoint of a sweep is the gradient a finite difference measures")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.set_keep_jacobian(true);
    REQUIRE(sweep.get_keep_jacobian());

    const Inputs in;
    REQUIRE(sweep.compute_Vs(in.gen_p, in.sgen_p, in.load_p, in.load_q,
                             CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-11) == 1);
    REQUIRE(sweep.adjoint_memory_bytes() > 0);

    const int dim_J = sweep.dim_J();
    REQUIRE(dim_J > 0);
    const IntVect theta_col = sweep.get_theta_col_of_bus();
    const IntVect vm_col = sweep.get_vm_col_of_bus();
    const IntVect p_row = sweep.get_p_row_of_bus();
    const IntVect q_row = sweep.get_q_row_of_bus();
    REQUIRE(theta_col.size() == NB_BUS);

    // project the loss's cotangent onto the Newton-Raphson unknowns, bus by bus
    const auto V = sweep.get_voltages();
    RealMatRM xbar = RealMatRM::Zero(NB_STEPS, dim_J);
    for(int i = 0; i < NB_STEPS; ++i){
        for(int b = 0; b < NB_BUS; ++b){
            const cplx_type gV(loss_weight_re(b), loss_weight_im(b));
            const cplx_type v = V(i, b);
            if(theta_col[b] >= 0) xbar(i, theta_col[b]) = (std::conj(v) * gV).imag();
            if(vm_col[b] >= 0) xbar(i, vm_col[b]) = (std::conj(v / std::abs(v)) * gV).real();
        }
    }

    const RealMatRM lambda = sweep.solve_JT(xbar);
    REQUIRE(lambda.rows() == NB_STEPS);
    REQUIRE(lambda.cols() == dim_J);
    for(int i = 0; i < NB_STEPS; ++i) REQUIRE(sweep.adjoint_row_ok()[static_cast<size_t>(i)] == 1);

    // lambda is the gradient with respect to the per-unit bus injection; an element's
    // own gradient is that, times how it enters its bus (a load subtracts, a generator
    // adds, both scaled by sn_mva)
    REQUIRE(p_row[LOAD_BUS] >= 0);
    REQUIRE(q_row[LOAD_BUS] >= 0);
    REQUIRE(p_row[GEN_BUS] >= 0);

    const real_type delta = 1e-4;   // MW / MVAr
    for(int i = 0; i < NB_STEPS; ++i){
        // ---- load P
        {
            Inputs plus = in, minus = in;
            plus.load_p(i, 0) += delta;
            minus.load_p(i, 0) -= delta;
            const real_type fd = (run_loss(plus) - run_loss(minus)) / (2. * delta);
            const real_type adjoint = -lambda(i, p_row[LOAD_BUS]) / SN_MVA;
            REQUIRE(adjoint == Approx(fd).margin(1e-7));
        }
        // ---- load Q
        {
            Inputs plus = in, minus = in;
            plus.load_q(i, 0) += delta;
            minus.load_q(i, 0) -= delta;
            const real_type fd = (run_loss(plus) - run_loss(minus)) / (2. * delta);
            const real_type adjoint = -lambda(i, q_row[LOAD_BUS]) / SN_MVA;
            REQUIRE(adjoint == Approx(fd).margin(1e-7));
        }
        // ---- generator P, at a PV bus (the other sign, and another corner of the ledger)
        {
            Inputs plus = in, minus = in;
            plus.gen_p(i, 1) += delta;
            minus.gen_p(i, 1) -= delta;
            const real_type fd = (run_loss(plus) - run_loss(minus)) / (2. * delta);
            const real_type adjoint = lambda(i, p_row[GEN_BUS]) / SN_MVA;
            REQUIRE(adjoint == Approx(fd).margin(1e-7));
        }
    }
}


TEST_CASE("refreshing the Jacobian at the solution moves it off the last iterate")
{
    // Newton-Raphson tests its residual BEFORE deciding it needs a new Jacobian, so
    // the J it stops on was built one iterate earlier. refresh_J_at_solution() is what
    // re-evaluates it at the voltage actually returned -- the point the adjoint's
    // gradient is exact at. The effect is of the order of the last Newton step, so a
    // loose tolerance is what makes it visible: at the 1e-11 the gradient test above
    // uses, the two Jacobians agree to well past its margin, which is precisely why
    // that test cannot stand in for this one.
    LSGrid grid = make_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    const CplxVect v_conv = grid.ac_pf(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-8);
    REQUIRE(v_conv.size() > 0);

    ls2g::NR_SparseLU algo;
    algo.set_lsgrid(&grid);
    const Eigen::SparseMatrix<cplx_type> Ybus = grid.get_Ybus_solver();
    const CplxVect Sbus = grid.get_Sbus_solver();
    const IntVect pv = grid.get_pv_solver_numpy();
    const IntVect pq = grid.get_pq_solver_numpy();
    const IntVect slack_ids = grid.get_slack_ids_solver_numpy();
    const RealVect slack_weights = grid.get_slack_weights_solver();

    // a tolerance loose enough that the last step taken is not negligible
    CplxVect V = CplxVect::Constant(Ybus.rows(), cplx_type(1., 0.));
    REQUIRE(algo.compute_pf(Ybus, V, Sbus, slack_ids, slack_weights, pv, pq, 30, 1e-4));

    const Eigen::SparseMatrix<real_type> J_last_iterate = algo.get_J();
    algo.refresh_J_at_solution();
    const Eigen::SparseMatrix<real_type> J_at_solution = algo.get_J();

    // same pattern -- refreshing values is all it does, which is what lets the batch
    // keep one symbolic factorization
    REQUIRE(J_at_solution.rows() == J_last_iterate.rows());
    REQUIRE(J_at_solution.nonZeros() == J_last_iterate.nonZeros());

    real_type biggest_change = 0.;
    for(Eigen::Index k = 0; k < J_last_iterate.nonZeros(); ++k){
        biggest_change = std::max(biggest_change,
                                  std::abs(J_at_solution.valuePtr()[k] - J_last_iterate.valuePtr()[k]));
    }
    REQUIRE(biggest_change > 1e-9);   // it really did re-evaluate

    // ... and it is a function of the converged state alone: refreshing again changes
    // nothing, so nothing here depends on how many times it is called
    algo.refresh_J_at_solution();
    const Eigen::SparseMatrix<real_type> J_again = algo.get_J();
    for(Eigen::Index k = 0; k < J_again.nonZeros(); ++k){
        REQUIRE(J_again.valuePtr()[k] == Approx(J_at_solution.valuePtr()[k]).margin(1e-14));
    }
}


TEST_CASE("a sweep keeps no Jacobian unless it was asked to")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    const Inputs in;
    REQUIRE(sweep.compute_Vs(in.gen_p, in.sgen_p, in.load_p, in.load_q,
                             CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-11) == 1);

    REQUIRE_FALSE(sweep.get_keep_jacobian());
    REQUIRE(sweep.adjoint_memory_bytes() == 0);
    REQUIRE_THROWS(sweep.solve_JT(RealMatRM::Zero(NB_STEPS, 1)));
}


TEST_CASE("a sweep refuses to keep a Jacobian no algorithm builds")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::DC_SparseLU);
    sweep.set_keep_jacobian(true);

    // through compute(), not the compute_Vs wrapper used everywhere else here: that
    // one turns EVERY std::exception into its historical -1 return (see its body), so
    // it would hide this error rather than report it
    const Inputs in;
    sweep.modify_gen_p(in.gen_p);
    sweep.modify_sgen_p(in.sgen_p);
    sweep.modify_load_p(in.load_p);
    sweep.modify_load_q(in.load_q);
    REQUIRE_THROWS(sweep.compute(CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-11));
}
