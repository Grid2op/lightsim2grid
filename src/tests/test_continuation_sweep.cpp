// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ContinuationSweep is a batch algorithm (a sibling of the four BaseBatchSweep
// instantiations, see batch_algorithm/ContinuationSweep.hpp): it traces the solution
// curve from the grid's own injections to a target state and stops at the nose.
//
// The oracle here is deliberately not another continuation: every traced point must be
// a solution of the ORDINARY powerflow at the injections that point claims, which is
// checked by re-solving with LSGrid::ac_pf. The rest pins what the class exists for --
// one symbolic factorization for the whole curve -- and the two failure modes that
// would otherwise produce a plausible-looking curve instead of an error: a direction
// whose sign is inverted (the voltages would RISE under load) and a zero direction.

#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/ContinuationSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContinuationSweep;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

const int NB_BUS = 4;
const int NB_GEN = 2;

// 4-bus radial feeder 0-1-2-3, one load at bus 3, a slack generator at bus 0 and a PV
// generator at bus 1. sn_mva deliberately != 1, so a direction built in MW that forgot
// to divide by it would be off by a factor of 100.
LSGrid make_grid()
{
    LSGrid g;
    g.set_sn_mva(100.);
    g.set_init_vm_pu(1.0);
    const RealVect vn = RealVect::Constant(NB_BUS, 138.);
    g.init_bus(static_cast<unsigned int>(NB_BUS), 1, vn, 0, 0);

    const RealVect r = RealVect::Constant(NB_BUS - 1, 0.01);
    const RealVect x = RealVect::Constant(NB_BUS - 1, 0.1);
    const CplxVect h = CplxVect::Zero(NB_BUS - 1);
    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for (int i = 0; i < NB_BUS - 1; ++i) { f(i) = i; t(i) = i + 1; }
    g.init_powerlines(r, x, h, f, t);

    RealVect lp(1), lq(1); lp << 50.; lq << 10.;
    Eigen::VectorXi lb(1); lb << NB_BUS - 1;
    g.init_loads(lp, lq, lb);

    RealVect p(NB_GEN), v(NB_GEN), q(NB_GEN), qmin(NB_GEN), qmax(NB_GEN);
    Eigen::VectorXi b(NB_GEN);
    p << 0., 20.;
    v << 1.02, 1.01;
    q << 0., 0.;
    qmin << -1000., -1000.;
    qmax << 1000., 1000.;
    b << 0, 1;
    g.init_generators_full(p, v, q, std::vector<bool>{true, true}, qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
    return g;
}

CplxVect flat(const LSGrid & g)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(g.total_bus()), cplx_type(1., 0.));
}

// the load's target active power at a given lambda, for a run that doubles it
real_type load_p_at(real_type lam) { return 50. * (1. + lam); }
real_type load_q_at(real_type lam) { return 10. * (1. + lam); }

// Doubles the single load's P and Q, generation held fixed. A sweep cannot be returned
// by value (BaseBatchSolverSynch is deliberately neither copyable nor movable -- the
// solvers hold a pointer to its own _grid_model member), so this configures one in place.
void setup_sweep(ContinuationSweep & sweep)
{
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    RealVect lp(1), lq(1);
    lp << load_p_at(1.);
    lq << load_q_at(1.);
    sweep.set_target_load_p(lp);
    sweep.set_target_load_q(lq);
}

}  // namespace

TEST_CASE("a continuation traces points that are genuine powerflow solutions", "[cpf]")
{
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    setup_sweep(sweep);
    sweep.set_stop_at_lam(1.0);
    sweep.compute(flat(grid), 20, 1e-10);

    REQUIRE(sweep.get_status() == 1);
    REQUIRE(sweep.nb_points() > 5);
    REQUIRE(sweep.get_lam()(0) == Approx(0.).margin(1e-14));
    REQUIRE(sweep.get_lam()(sweep.nb_points() - 1) == Approx(1.).margin(1e-12));

    // every traced point, re-solved with an ordinary powerflow at its own injections
    for (Eigen::Index i = 0; i < sweep.nb_points(); ++i) {
        const real_type lam = sweep.get_lam()(i);
        LSGrid ref = make_grid();
        ref.change_p_load(0, load_p_at(lam));
        ref.change_q_load(0, load_q_at(lam));
        const CplxVect V = ref.ac_pf(flat(ref), 30, 1e-11);
        REQUIRE(V.size() > 0);
        for (Eigen::Index b = 0; b < V.size(); ++b) {
            REQUIRE(std::abs(V(b) - sweep.get_voltages()(i, b)) < 1e-8);
        }
    }
}

TEST_CASE("the whole curve costs one symbolic factorization", "[cpf]")
{
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    setup_sweep(sweep);
    sweep.compute(flat(grid), 20, 1e-10);

    REQUIRE(sweep.nb_points() > 10);
    // ONE analyze however many points were traced -- the reason a continuation belongs
    // in the batch layer rather than in a loop around LSGrid::ac_pf.
    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
    REQUIRE(sweep.get_linear_solver_stats().nb_refactorize > sweep.nb_points());
}

TEST_CASE("loading the grid lowers its voltages", "[cpf]")
{
    // The sign of the load direction. A load is a NEGATIVE injection (LoadContainer::
    // fillSbus subtracts it), so a direction that increases the load must decrease the
    // voltages; inverting it traces the grid UNLOADING, converging happily all the way.
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    setup_sweep(sweep);
    sweep.compute(flat(grid), 20, 1e-10);
    REQUIRE(sweep.nb_points() > 2);

    // the direction itself: one load of +50 MW at the end bus, ie -0.5 pu of injection
    const CplxVect & dir = sweep.get_direction_solver();
    REQUIRE(std::real(dir(NB_BUS - 1)) == Approx(-0.5));
    REQUIRE(std::imag(dir(NB_BUS - 1)) == Approx(-0.1));

    const Eigen::Index last = sweep.nb_points() - 1;
    const real_type vm_first = std::abs(sweep.get_voltages()(0, NB_BUS - 1));
    const real_type vm_last = std::abs(sweep.get_voltages()(last, NB_BUS - 1));
    REQUIRE(vm_last < vm_first);
    for (Eigen::Index i = 1; i < sweep.nb_points(); ++i) {
        REQUIRE(std::abs(sweep.get_voltages()(i, NB_BUS - 1)) <=
                std::abs(sweep.get_voltages()(i - 1, NB_BUS - 1)) + 1e-9);
    }
}

TEST_CASE("the tangent's lambda component collapses at the nose", "[cpf]")
{
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    setup_sweep(sweep);
    sweep.compute(flat(grid), 20, 1e-10);

    REQUIRE(sweep.get_status() == 1);  // stopped at the nose, not on the step cap
    // strictly positive throughout -- this parameterisation cannot make it change sign,
    // it only tends to zero (see ContinuationSweep.hpp)
    for (Eigen::Index i = 0; i + 1 < sweep.nb_points(); ++i) {
        REQUIRE(sweep.get_tangent_lam()(i) > 0.);
    }
    const Eigen::Index before_last = sweep.nb_points() - 2;
    REQUIRE(sweep.get_tangent_lam()(before_last) < sweep.get_tangent_lam()(0));
}

TEST_CASE("a zero direction is refused rather than traced", "[cpf]")
{
    // Not a degenerate curve but a meaningless run: J . z = 0 gives z = 0, the tangent
    // normalises to (0, ..., 1), and lambda would march to its target reporting success
    // for a continuation that continued nothing.
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    RealVect same(1);
    same << 50.;  // exactly the grid's own value
    sweep.set_target_load_p(same);
    REQUIRE_THROWS_AS(sweep.compute(flat(grid), 20, 1e-10), std::runtime_error);

    // and so is setting no target at all
    ContinuationSweep untargeted(grid);
    untargeted.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE_THROWS_AS(untargeted.compute(flat(grid), 20, 1e-10), std::runtime_error);
}

TEST_CASE("a direction the slack absorbs entirely is refused", "[cpf]")
{
    // The companion of the zero-direction guard, and the reason slack machines are NOT
    // excluded from gen_steering the way MATPOWER excludes them from a target case: at a
    // SINGLE slack bus the machine's target_p is inert, so a direction that only scales
    // it is non-zero, produces a non-zero tangent (that bus does have a P equation, paired
    // with the slack-absorbed unknown), and yet moves no voltage whatsoever -- the slack
    // absorption cancels it one for one. Left alone, lambda marches to its target over a
    // curve on which nothing happens, and the run reports success.
    LSGrid grid = make_grid();
    grid.change_p_gen(0, 30.);  // the slack machine, given a non-zero setpoint to scale
    ContinuationSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    RealVect gp(NB_GEN);
    gp << 60., 20.;  // ONLY the slack machine moves
    sweep.set_target_gen_p(gp);
    REQUIRE_THROWS_AS(sweep.compute(flat(grid), 20, 1e-10), std::runtime_error);
}

TEST_CASE("a single slack machine's setpoint does not change the solution", "[cpf]")
{
    // What makes the test above true, stated on its own: with one slack bus the machine's
    // target_p is an input the powerflow ignores -- its output is whatever balances the
    // grid. (With a DISTRIBUTED slack this stops holding, which is exactly why excluding
    // slack machines from the steering would be wrong there.)
    LSGrid a = make_grid();
    a.change_p_gen(0, 30.);
    const CplxVect Va = a.ac_pf(flat(a), 30, 1e-11);

    LSGrid b = make_grid();
    b.change_p_gen(0, 60.);
    const CplxVect Vb = b.ac_pf(flat(b), 30, 1e-11);

    REQUIRE(Va.size() == Vb.size());
    for (Eigen::Index i = 0; i < Va.size(); ++i) REQUIRE(Va(i) == Vb(i));
}

TEST_CASE("a continuation refuses an algorithm with no Jacobian", "[cpf]")
{
    LSGrid grid = make_grid();

    ContinuationSweep gauss(grid);
    setup_sweep(gauss);
    gauss.change_algorithm(AlgorithmType::GaussSeidel);
    // change_algorithm() clears the object, so the target has to be set again
    RealVect lp(1); lp << load_p_at(1.);
    gauss.set_target_load_p(lp);
    REQUIRE_THROWS_AS(gauss.compute(flat(grid), 20, 1e-10), std::runtime_error);

    ContinuationSweep dc(grid);
    setup_sweep(dc);
    dc.change_algorithm(AlgorithmType::DC_SparseLU);
    dc.set_target_load_p(lp);
    REQUIRE_THROWS_AS(dc.compute(flat(grid), 20, 1e-10), std::runtime_error);
}

TEST_CASE("a continuation cannot be split over threads", "[cpf]")
{
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    REQUIRE_FALSE(sweep.supports_multithread());
    REQUIRE_THROWS_AS(sweep.set_nb_thread(4), std::runtime_error);
    sweep.set_nb_thread(1);  // still allowed
    REQUIRE(sweep.get_nb_thread() == 1);
}

TEST_CASE("stop_at_lam lands exactly on the requested lambda", "[cpf]")
{
    LSGrid grid = make_grid();
    ContinuationSweep sweep(grid);
    setup_sweep(sweep);
    sweep.set_stop_at_lam(0.4);
    sweep.compute(flat(grid), 20, 1e-10);

    REQUIRE(sweep.get_status() == 1);
    REQUIRE(sweep.get_lam_max() == Approx(0.4).margin(1e-12));

    LSGrid ref = make_grid();
    ref.change_p_load(0, load_p_at(0.4));
    ref.change_q_load(0, load_q_at(0.4));
    const CplxVect V = ref.ac_pf(flat(ref), 30, 1e-11);
    REQUIRE(V.size() > 0);
    const Eigen::Index last = sweep.nb_points() - 1;
    for (Eigen::Index b = 0; b < V.size(); ++b) {
        REQUIRE(std::abs(V(b) - sweep.get_voltages()(last, b)) < 1e-8);
    }
}

TEST_CASE("an adaptive step reaches the same nose", "[cpf]")
{
    LSGrid grid = make_grid();

    ContinuationSweep fixed(grid);
    setup_sweep(fixed);
    fixed.compute(flat(grid), 20, 1e-10);

    ContinuationSweep adapt(grid);
    setup_sweep(adapt);
    adapt.set_adapt_step(true);
    adapt.compute(flat(grid), 20, 1e-10);

    REQUIRE(adapt.get_status() == 1);
    REQUIRE(adapt.get_lam_max() == Approx(fixed.get_lam_max()).epsilon(1e-3));
    REQUIRE(adapt.nb_points() < fixed.nb_points());
}

TEST_CASE("an exact tangent reaches the same nose", "[cpf]")
{
    LSGrid grid = make_grid();

    ContinuationSweep loose(grid);
    setup_sweep(loose);
    loose.compute(flat(grid), 20, 1e-10);

    ContinuationSweep exact(grid);
    setup_sweep(exact);
    exact.set_exact_tangent(true);
    exact.compute(flat(grid), 20, 1e-10);

    REQUIRE(exact.get_status() == 1);
    REQUIRE(exact.get_lam_max() == Approx(loose.get_lam_max()).epsilon(1e-3));
    REQUIRE(exact.get_linear_solver_stats().nb_analyze == 1);
}
