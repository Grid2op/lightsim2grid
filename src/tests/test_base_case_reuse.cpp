// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Keeping a batch's base case between two compute() calls is only worth anything if
// the second call gives exactly what it would have given without the cache. Every
// test here is that comparison: the same batch run twice on one object, against the
// same batch run on a fresh one. What differs is what happened in between -- a
// contingency registered, the thread count changed, the grid solved by somebody else
// -- and each of those either invalidates the cache or must not affect the answer.
//
// The cache is three nested levels (see BaseBatchSolverSynch's block comment):
// clear_grid_results() (L1) -> clear_batch_inputs() (L2) -> clear_batch_outputs()
// (L3). "The base case was kept" is exactly "neither L1 nor L2 was dropped since the
// last compute()", which is what base_case_was_reused() reports.

#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::InjectionSweep;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using CplxMat = Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int NB_BUS = 6;
const int NB_GEN = 2;
const int NB_STEPS = 7;

// a 6-bus radial feeder with two loads, a slack generator and a PV generator
LSGrid make_grid()
{
    LSGrid g;
    g.set_sn_mva(100.);
    g.set_init_vm_pu(1.0);
    g.init_bus(static_cast<unsigned int>(NB_BUS), 1, RealVect::Constant(NB_BUS, 138.), 0, 0);

    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for(int i = 0; i < NB_BUS - 1; ++i){ f(i) = i; t(i) = i + 1; }
    g.init_powerlines(RealVect::Constant(NB_BUS - 1, 0.01), RealVect::Constant(NB_BUS - 1, 0.1),
                      CplxVect::Zero(NB_BUS - 1), f, t);

    RealVect lp(2), lq(2); lp << 40., 30.; lq << 8., 6.;
    Eigen::VectorXi lb(2); lb << 3, NB_BUS - 1;
    g.init_loads(lp, lq, lb);

    RealVect p(NB_GEN), v(NB_GEN), q(NB_GEN), qmin(NB_GEN), qmax(NB_GEN);
    Eigen::VectorXi b(NB_GEN);
    p << 0., 25.;
    v << 1.02, 1.01;
    q << 0., 0.;
    qmin << -1000., -1000.;
    qmax << 1000., 1000.;
    b << 0, 1;
    g.init_generators_full(p, v, q, std::vector<bool>{true, true}, qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
    return g;
}

CplxVect flat() { return CplxVect::Constant(NB_BUS, cplx_type(1., 0.)); }

struct Inputs
{
    RealMat gen_p{NB_STEPS, NB_GEN};
    RealMat sgen_p{NB_STEPS, 0};
    RealMat load_p{NB_STEPS, 2};
    RealMat load_q{NB_STEPS, 2};

    explicit Inputs(real_type offset = 0.)
    {
        for(int i = 0; i < NB_STEPS; ++i){
            const real_type jitter = static_cast<real_type>((i * 37) % 11) / 10. + offset;
            gen_p(i, 0) = 0.;
            gen_p(i, 1) = 10. + 15. * jitter;
            load_p(i, 0) = 30. + 20. * jitter;
            load_p(i, 1) = 20. + 15. * jitter;
            load_q(i, 0) = 5. + 5. * jitter;
            load_q(i, 1) = 4. + 4. * jitter;
        }
    }
};

void run(InjectionSweep & sweep, const Inputs & in)
{
    sweep.modify_gen_p(in.gen_p);
    sweep.modify_sgen_p(in.sgen_p);
    sweep.modify_load_p(in.load_p);
    sweep.modify_load_q(in.load_q);
    sweep.compute(flat(), 30, 1e-11);
    REQUIRE(sweep.get_status() == 1);
}

void check_same(const CplxMat & got, const CplxMat & expected)
{
    REQUIRE(got.rows() == expected.rows());
    REQUIRE(got.cols() == expected.cols());
    for(Eigen::Index i = 0; i < got.rows(); ++i){
        for(Eigen::Index b = 0; b < got.cols(); ++b){
            REQUIRE(got(i, b).real() == Approx(expected(i, b).real()).margin(1e-11));
            REQUIRE(got(i, b).imag() == Approx(expected(i, b).imag()).margin(1e-11));
        }
    }
}

// the same batch, on an object that has never computed anything
CplxMat fresh_result(const Inputs & in, const CplxVect & v_start = flat())
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    sweep.modify_gen_p(in.gen_p);
    sweep.modify_sgen_p(in.sgen_p);
    sweep.modify_load_p(in.load_p);
    sweep.modify_load_q(in.load_q);
    sweep.compute(v_start, 30, 1e-11);
    REQUIRE(sweep.get_status() == 1);
    return sweep.get_voltages();
}

}  // namespace


TEST_CASE("a kept base case gives what a fresh one would have")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(sweep.get_reuse_base_case());          // on by default

    const Inputs first;
    run(sweep, first);
    REQUIRE_FALSE(sweep.base_case_was_reused());   // nothing to reuse yet
    check_same(sweep.get_voltages(), fresh_result(first));

    // ... and again, with different injections: this one keeps the base case
    const Inputs second(0.7);
    run(sweep, second);
    REQUIRE(sweep.base_case_was_reused());
    check_same(sweep.get_voltages(), fresh_result(second));
}


TEST_CASE("keeping the base case is what stops the second batch analyzing again")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    const Inputs first;
    const Inputs second(0.7);

    run(sweep, first);
    const std::size_t analyze_after_first = sweep.get_linear_solver_stats().nb_analyze;
    REQUIRE(analyze_after_first >= 1);

    run(sweep, second);
    // the whole point: the symbolic factorization of the first call is still the one
    // being used, so the second call never asked for another
    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == analyze_after_first);
    REQUIRE(sweep.get_linear_solver_stats().nb_refactorize > 0);

    // turned off, every call pays for its own again
    LSGrid grid2 = make_grid();
    InjectionSweep no_reuse(grid2);
    no_reuse.change_algorithm(AlgorithmType::NR_SparseLU);
    no_reuse.set_reuse_base_case(false);
    run(no_reuse, first);
    const std::size_t analyze_off = no_reuse.get_linear_solver_stats().nb_analyze;
    run(no_reuse, second);
    REQUIRE_FALSE(no_reuse.base_case_was_reused());
    REQUIRE(no_reuse.get_linear_solver_stats().nb_analyze > analyze_off);
    check_same(no_reuse.get_voltages(), fresh_result(second));
}


TEST_CASE("a starting voltage is never the kept part of a base case")
{
    // The base case is kept; the voltage a call starts from is not -- it is mapped onto
    // the kept labelling every time. That distinction is observable: the reference slack
    // holds the angle it was started at, so a batch started somewhere else converges to
    // a rotated solution, and a cached start would silently give the FIRST call's.
    const CplxVect elsewhere = CplxVect::Constant(NB_BUS, cplx_type(0.95, 0.05));
    const Inputs in;

    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);
    run(sweep, in);                                   // builds the base case, flat start
    const CplxMat from_flat = sweep.get_voltages();

    sweep.modify_gen_p(in.gen_p);
    sweep.modify_sgen_p(in.sgen_p);
    sweep.modify_load_p(in.load_p);
    sweep.modify_load_q(in.load_q);
    sweep.compute(elsewhere, 30, 1e-11);              // ... kept, but started elsewhere
    REQUIRE(sweep.get_status() == 1);
    REQUIRE(sweep.base_case_was_reused());

    check_same(sweep.get_voltages(), fresh_result(in, elsewhere));

    // and the two starts really do give different answers, so the check above is not
    // passing for want of anything to notice
    REQUIRE((sweep.get_voltages() - from_flat).cwiseAbs().maxCoeff() > 1e-6);
}


TEST_CASE("a different number of simulations rebuilds the base case")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    const Inputs in;
    run(sweep, in);

    // a smaller batch: the row count is part of what a base case was prepared for
    sweep.clear();
    Inputs small;
    RealMat gen_p = small.gen_p.topRows(3);
    RealMat sgen_p = small.sgen_p.topRows(3);
    RealMat load_p = small.load_p.topRows(3);
    RealMat load_q = small.load_q.topRows(3);
    sweep.modify_gen_p(gen_p);
    sweep.modify_sgen_p(sgen_p);
    sweep.modify_load_p(load_p);
    sweep.modify_load_q(load_q);
    sweep.compute(flat(), 30, 1e-11);
    REQUIRE(sweep.get_status() == 1);
    REQUIRE_FALSE(sweep.base_case_was_reused());
    REQUIRE(sweep.get_voltages().rows() == 3);
}


TEST_CASE("changing the thread count rebuilds the base case")
{
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    const Inputs in;
    run(sweep, in);
    const CplxMat one_thread = sweep.get_voltages();

    sweep.set_nb_thread(3);
    run(sweep, in);
    REQUIRE_FALSE(sweep.base_case_was_reused());   // the workers are not the ones it was built with
    check_same(sweep.get_voltages(), one_thread);

    run(sweep, in);
    REQUIRE(sweep.base_case_was_reused());         // ... and now they are
    check_same(sweep.get_voltages(), one_thread);
}


TEST_CASE("registering a contingency rebuilds the base case")
{
    LSGrid grid = make_grid();
    ContingencyAnalysis analysis(grid);
    analysis.change_algorithm(AlgorithmType::NR_SparseLU);

    analysis.add_n1(1);
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE_FALSE(analysis.base_case_was_reused());
    const CplxMat with_one = analysis.get_voltages();

    // a second contingency: what the graph walk settled no longer describes this batch
    analysis.add_n1(2);
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE_FALSE(analysis.base_case_was_reused());

    LSGrid grid2 = make_grid();
    ContingencyAnalysis fresh(grid2);
    fresh.change_algorithm(AlgorithmType::NR_SparseLU);
    fresh.add_n1(1);
    fresh.add_n1(2);
    fresh.compute(flat(), 30, 1e-11);
    check_same(analysis.get_voltages(), fresh.get_voltages());

    // re-running the same set keeps it
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE(analysis.base_case_was_reused());
    check_same(analysis.get_voltages(), fresh.get_voltages());
}


TEST_CASE("dropping computed results keeps the base case")
{
    // Tightening a violation threshold changes what is CHECKED of a row, and nothing
    // the base case is made of: it is an L3 modifier, and L3 leaves the algorithm alone.
    //
    // It did not always: clear_results_only() used to reset the algorithm, taking the
    // ledger, the sparsity and the factorization with it, and a base case kept across
    // that is not a stale answer but a solve against a default-constructed system --
    // which is exactly how the limit-violation suite segfaulted. The levels make that
    // combination unrepresentable: the algorithm belongs to L2, so nothing can drop it
    // without also dropping the base case that would have used it.
    LSGrid grid = make_grid();
    ContingencyAnalysis analysis(grid);
    analysis.change_algorithm(AlgorithmType::NR_SparseLU);
    analysis.add_n1(1);
    analysis.compute(flat(), 30, 1e-11);
    const CplxMat before = analysis.get_voltages();

    analysis.set_violation_threshold(0.5);        // tightened -> clear_results_only()
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE(analysis.base_case_was_reused());
    check_same(analysis.get_voltages(), before);

    // the same through the public entry point a caller can reach directly
    analysis.clear_results_only();
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE(analysis.base_case_was_reused());
    check_same(analysis.get_voltages(), before);
}


TEST_CASE("the three cache levels nest")
{
    // L3 drops the results; L2 also drops the batch inputs and the algorithm they
    // configured; L1 also drops what was read off the grid. Each is observable: the
    // results by their row count, L2 by base_case_was_reused() and by a fresh symbolic
    // analysis, and every one of them by the answer still being right afterwards.
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    const Inputs in;
    run(sweep, in);
    const CplxMat expected = sweep.get_voltages();
    check_same(expected, fresh_result(in));

    // ---- L3: the results, and only the results
    sweep.clear_batch_outputs();
    REQUIRE(sweep.get_voltages().rows() == 0);
    const std::size_t analyze_before_l3 = sweep.get_linear_solver_stats().nb_analyze;
    run(sweep, in);
    REQUIRE(sweep.base_case_was_reused());
    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == analyze_before_l3);
    check_same(sweep.get_voltages(), expected);

    // ---- L2: + the batch inputs and the algorithm
    sweep.clear_batch_inputs();
    REQUIRE(sweep.get_voltages().rows() == 0);   // L3 went with it
    const std::size_t analyze_before_l2 = sweep.get_linear_solver_stats().nb_analyze;
    run(sweep, in);
    REQUIRE_FALSE(sweep.base_case_was_reused());
    REQUIRE(sweep.get_linear_solver_stats().nb_analyze > analyze_before_l2);
    check_same(sweep.get_voltages(), expected);

    // ---- L1: + what was read off the grid
    sweep.clear_grid_results();
    REQUIRE(sweep.get_voltages().rows() == 0);   // L2 and L3 went with it
    run(sweep, in);
    REQUIRE_FALSE(sweep.base_case_was_reused());
    check_same(sweep.get_voltages(), expected);
}


TEST_CASE("a cache level is not a registration")
{
    // Dropping a level says "what I built from this is stale", never "forget what you
    // were asked to compute". Only clear() does the latter -- a different axis, which
    // is why it is not a fourth level.
    LSGrid grid = make_grid();
    ContingencyAnalysis analysis(grid);
    analysis.change_algorithm(AlgorithmType::NR_SparseLU);
    analysis.add_n1(1);
    analysis.add_n1(2);
    REQUIRE(analysis.my_defaults().size() == 2);

    analysis.compute(flat(), 30, 1e-11);
    const CplxMat expected = analysis.get_voltages();

    analysis.clear_grid_results();                 // the top of the hierarchy
    REQUIRE(analysis.my_defaults().size() == 2);   // ... and the set is still registered
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE_FALSE(analysis.base_case_was_reused());
    check_same(analysis.get_voltages(), expected);

    analysis.clear();
    REQUIRE(analysis.my_defaults().empty());
}


TEST_CASE("switching algorithm rebuilds everything but keeps the contingencies")
{
    // change_algorithm() is L1 on purpose: a different algorithm may be a different
    // family (AC <-> DC), and then nothing read off the grid survives. It is not a
    // clear(): the contingencies were registered by the caller, not derived from the
    // algorithm.
    LSGrid grid = make_grid();
    ContingencyAnalysis analysis(grid);
    analysis.change_algorithm(AlgorithmType::NR_SparseLU);
    analysis.add_n1(1);
    analysis.compute(flat(), 30, 1e-11);

    analysis.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(analysis.my_defaults().size() == 1);
    analysis.compute(flat(), 30, 1e-11);
    REQUIRE_FALSE(analysis.base_case_was_reused());

    LSGrid grid2 = make_grid();
    ContingencyAnalysis fresh(grid2);
    fresh.change_algorithm(AlgorithmType::NR_SparseLU);
    fresh.add_n1(1);
    fresh.compute(flat(), 30, 1e-11);
    check_same(analysis.get_voltages(), fresh.get_voltages());
}


TEST_CASE("the seed a row starts from is part of the base case")
{
    // set_init_from_n_powerflow decides whether every row starts from the "n" solve's
    // answer or from the caller's own Vinit -- L2, and it has to be, because the kept
    // base case is where that answer is kept.
    LSGrid grid = make_grid();
    InjectionSweep sweep(grid);
    sweep.change_algorithm(AlgorithmType::NR_SparseLU);

    const Inputs in;
    run(sweep, in);
    const CplxMat from_vinit = sweep.get_voltages();

    sweep.set_init_from_n_powerflow(true);
    run(sweep, in);
    REQUIRE_FALSE(sweep.base_case_was_reused());
    // same root either way (a converged powerflow does not depend on its seed here),
    // which is the point: what must not happen is reading a seed that was never stored
    check_same(sweep.get_voltages(), from_vinit);

    run(sweep, in);
    REQUIRE(sweep.base_case_was_reused());
    check_same(sweep.get_voltages(), from_vinit);
}
