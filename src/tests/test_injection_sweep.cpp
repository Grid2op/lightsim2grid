// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// TimeSeries and InjectionSweep are two instantiations of the same template
// (batch_algorithm/BaseInjectionSweep.hpp): same inputs, same outputs, same interface. They differ
// only in how each step is initialized -- TimeSeries warm starts step i with the solution
// of step i-1, InjectionSweep restarts every step from the same seed.
//
// The tests below pin what follows from that difference, and nothing else:
//   * both must converge to the SAME solution (the seed moves the starting point, not the
//     root of a well conditioned system);
//   * InjectionSweep's results must not depend on the order the steps are given in, nor on
//     how many threads they were split over -- exactly the two properties a chained
//     TimeSeries cannot offer, which is why it rejects nb_thread > 1.

#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseInjectionSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::InjectionSweep;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::TimeSeries;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using CplxMat = Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int NB_BUS = 4;
const int NB_GEN = 2;
const int NB_STEPS = 12;

// 4-bus radial feeder 0-1-2-3 (r = 0.01, x = 0.1 pu), one load at bus 3, a slack generator
// at bus 0 and a PV generator at bus 1. sn_mva deliberately != 1.
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

// One injection "scenario" per row. Deliberately NOT a smooth ramp: consecutive rows jump
// around, so warm starting from the previous one is neither obviously good nor obviously
// harmless -- which is the situation InjectionSweep exists for.
struct Scenarios
{
    RealMat gen_p{NB_STEPS, NB_GEN};
    RealMat sgen_p{NB_STEPS, 0};
    RealMat load_p{NB_STEPS, 1};
    RealMat load_q{NB_STEPS, 1};

    Scenarios()
    {
        for (int i = 0; i < NB_STEPS; ++i) {
            // pseudo random but fully deterministic
            const real_type jitter = static_cast<real_type>((i * 37) % 11) / 10.;
            gen_p(i, 0) = 0.;
            gen_p(i, 1) = 10. + 15. * jitter;
            load_p(i, 0) = 30. + 40. * jitter;
            load_q(i, 0) = 5. + 10. * jitter;
        }
    }

    // the same scenarios in a different order
    Scenarios permuted(const std::vector<int> & perm) const
    {
        Scenarios res;
        for (int i = 0; i < NB_STEPS; ++i) {
            res.gen_p.row(i) = gen_p.row(perm[i]);
            res.load_p.row(i) = load_p.row(perm[i]);
            res.load_q.row(i) = load_q.row(perm[i]);
        }
        return res;
    }
};

template<class BatchAlgo>
CplxMat run(BatchAlgo & algo, const Scenarios & sc)
{
    const int status = algo.compute_Vs(sc.gen_p, sc.sgen_p, sc.load_p, sc.load_q,
                                       CplxVect::Constant(NB_BUS, cplx_type(1., 0.)),
                                       30, 1e-11);
    REQUIRE(status == 1);
    return algo.get_voltages();
}

// a fixed, non trivial permutation of [0, NB_STEPS)
std::vector<int> shuffled()
{
    std::vector<int> perm(NB_STEPS);
    for (int i = 0; i < NB_STEPS; ++i) perm[i] = (i * 5 + 3) % NB_STEPS;  // 5 is coprime with 12
    return perm;
}

}  // namespace

TEST_CASE("InjectionSweep converges to the same solution as TimeSeries", "[injection_sweep]")
{
    const Scenarios sc;

    SECTION("ac")
    {
        LSGrid g_ts = make_grid();
        LSGrid g_is = make_grid();
        g_ts.change_algorithm(AlgorithmType::NR_SparseLU);
        g_is.change_algorithm(AlgorithmType::NR_SparseLU);

        TimeSeries ts(g_ts);
        InjectionSweep is(g_is);
        const CplxMat V_ts = run(ts, sc);
        const CplxMat V_is = run(is, sc);

        REQUIRE(V_ts.rows() == NB_STEPS);
        REQUIRE((V_is - V_ts).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-9));
        // and they really are different problems from step to step
        REQUIRE(V_is.col(NB_BUS - 1).cwiseAbs().maxCoeff() -
                V_is.col(NB_BUS - 1).cwiseAbs().minCoeff() > 1e-3);
    }

    SECTION("dc")
    {
        LSGrid g_ts = make_grid();
        LSGrid g_is = make_grid();
        g_ts.change_algorithm(AlgorithmType::DC_SparseLU);
        g_is.change_algorithm(AlgorithmType::DC_SparseLU);

        TimeSeries ts(g_ts);
        InjectionSweep is(g_is);
        ts.change_algorithm(AlgorithmType::DC_SparseLU);
        is.change_algorithm(AlgorithmType::DC_SparseLU);
        const CplxMat V_ts = run(ts, sc);
        const CplxMat V_is = run(is, sc);
        REQUIRE((V_is - V_ts).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-9));
    }
}

TEST_CASE("InjectionSweep results do not depend on the order of the steps", "[injection_sweep]")
{
    // The defining invariant: because every step restarts from the same seed, permuting
    // the input rows permutes the output rows and changes nothing else -- bit for bit.
    const Scenarios sc;
    const std::vector<int> perm = shuffled();

    LSGrid g1 = make_grid();
    LSGrid g2 = make_grid();
    g1.change_algorithm(AlgorithmType::NR_SparseLU);
    g2.change_algorithm(AlgorithmType::NR_SparseLU);

    InjectionSweep is1(g1);
    InjectionSweep is2(g2);
    const CplxMat V = run(is1, sc);
    const CplxMat V_perm = run(is2, sc.permuted(perm));

    for (int i = 0; i < NB_STEPS; ++i) {
        REQUIRE(V_perm.row(i) == V.row(perm[i]));
    }
}

TEST_CASE("InjectionSweep results do not depend on the number of threads", "[injection_sweep]")
{
    const Scenarios sc;
    LSGrid g_ref = make_grid();
    g_ref.change_algorithm(AlgorithmType::NR_SparseLU);
    InjectionSweep ref(g_ref);
    const CplxMat V_ref = run(ref, sc);
    REQUIRE(ref.thread_init_time() == 0.);  // nothing to build when nb_thread == 1

    for (int nb_thread : {2, 3, 4, 8}) {
        LSGrid g = make_grid();
        g.change_algorithm(AlgorithmType::NR_SparseLU);
        InjectionSweep is(g);
        is.set_nb_thread(nb_thread);
        REQUIRE(is.get_nb_thread() == nb_thread);
        const CplxMat V = run(is, sc);
        REQUIRE(V == V_ref);  // bit for bit: the threads solve independent problems
        REQUIRE(is.nb_solved() == NB_STEPS);
    }
}

TEST_CASE("only a batch whose steps are independent can be multi-threaded", "[injection_sweep]")
{
    LSGrid g_ts = make_grid();
    LSGrid g_is = make_grid();
    TimeSeries ts(g_ts);
    InjectionSweep is(g_is);

    REQUIRE_FALSE(ts.supports_multithread());
    REQUIRE(is.supports_multithread());
    REQUIRE(ts.get_nb_thread() == 1);
    REQUIRE(is.get_nb_thread() == 1);

    SECTION("TimeSeries rejects more than one thread")
    {
        REQUIRE_THROWS_AS(ts.set_nb_thread(2), std::runtime_error);
        REQUIRE(ts.get_nb_thread() == 1);  // a rejected assignment changes nothing
        // 1 and below stay legal (they mean "sequential", which is what it does)
        REQUIRE_NOTHROW(ts.set_nb_thread(1));
        REQUIRE_NOTHROW(ts.set_nb_thread(0));
        REQUIRE(ts.get_nb_thread() == 1);
    }

    SECTION("values below 1 are clamped rather than rejected")
    {
        is.set_nb_thread(0);
        REQUIRE(is.get_nb_thread() == 1);
        is.set_nb_thread(-5);
        REQUIRE(is.get_nb_thread() == 1);
    }
}

TEST_CASE("InjectionSweep can seed every step with the n powerflow", "[injection_sweep]")
{
    // init_from_n_powerflow picks WHICH voltage is used as the seed; for InjectionSweep it
    // applies to every step (like ContingencyAnalysis), not only the first one. Either way
    // a well conditioned system converges to the same root.
    const Scenarios sc;
    LSGrid g_v = make_grid();
    LSGrid g_n = make_grid();
    g_v.change_algorithm(AlgorithmType::NR_SparseLU);
    g_n.change_algorithm(AlgorithmType::NR_SparseLU);

    InjectionSweep from_v(g_v);
    InjectionSweep from_n(g_n);
    REQUIRE_FALSE(from_v.get_init_from_n_powerflow());  // the default
    from_n.set_init_from_n_powerflow(true);
    REQUIRE(from_n.get_init_from_n_powerflow());

    const CplxMat V_v = run(from_v, sc);
    const CplxMat V_n = run(from_n, sc);
    REQUIRE((V_n - V_v).cwiseAbs().maxCoeff() == Approx(0.).margin(1e-9));

    // ... and the seed is computed once, before the threads are spawned, so threading it
    // changes nothing either
    LSGrid g_th = make_grid();
    g_th.change_algorithm(AlgorithmType::NR_SparseLU);
    InjectionSweep threaded(g_th);
    threaded.set_init_from_n_powerflow(true);
    threaded.set_nb_thread(4);
    REQUIRE(run(threaded, sc) == V_n);
}

TEST_CASE("a multi-threaded InjectionSweep reports errors like the sequential one", "[injection_sweep]")
{
    // a malformed input must surface as an exception on the caller's thread, not crash a
    // worker: run_step_range captures it and compute_Vs rethrows the first one.
    Scenarios sc;
    LSGrid g = make_grid();
    g.change_algorithm(AlgorithmType::NR_SparseLU);
    InjectionSweep is(g);
    is.set_nb_thread(4);

    const RealMat wrong_shape = RealMat::Zero(NB_STEPS, NB_GEN + 3);
    REQUIRE_THROWS_AS(is.compute_Vs(wrong_shape, sc.sgen_p, sc.load_p, sc.load_q,
                                    CplxVect::Constant(NB_BUS, cplx_type(1., 0.)), 30, 1e-11),
                      std::runtime_error);
}
