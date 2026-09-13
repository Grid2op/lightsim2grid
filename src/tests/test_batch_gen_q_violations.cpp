// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// `compute_gen_q_violations`: the batch-side reactive-limit check of the voltage
// regulating generators (see batch_algorithm/GenQCheck.hpp).
//
// The reactive output it reports is RE-DERIVED from the algorithm's per-bus mismatch and
// its controller list, because a batch publishes voltages and nothing else. Every test
// here therefore pins that re-derivation against the one thing that cannot be wrong: what
// a single-shot LSGrid::ac_pf publishes on the very same grid (GeneratorContainer's
// res_q_, through LSGrid::get_gen_res). If the two disagree it is the re-derivation that
// is wrong -- the sharing rule (LSGrid::_split_q_residual_per_bus), the controller
// write-back (LSGrid::_write_back_controller_q), or the choice between them.
//
// A generator's reactive LIMITS are also its sharing key (span = qmax - qmin) and can
// only be set at init_generators time, so a test that needs a violation builds the grid
// with the limits it wants from the start and compares against an ac_pf of that same
// grid -- never against one with different ranges.

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::TimeSeries;
using ls2g::ViolationElementType;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type V_SET = 1.05;
const real_type LOAD_P = 80.;
const real_type LOAD_Q = 60.;
const int NB_BUS = 4;
const real_type WIDE_Q = 1000.;

struct GenSpec
{
    int bus;
    real_type vset;
    real_type p;
    real_type min_q;
    real_type max_q;
    int regulated_bus = -1;  // -1: regulates its own bus (the classical PV path)
};

// The 4-bus radial feeder 0-1-2-3 (r = 0.01, x = 0.1 pu, sn_mva = 100) shared by
// test_batch_voltage_control.cpp and test_scenario_sweep_violations.cpp, with the load on
// bus 3. `meshed` adds a fourth line 0--3, so that a single line outage leaves the load
// connected and the contingency rows actually converge.
// Generator 0 is the slack.
LSGrid make_grid(const std::vector<GenSpec> & gens, bool meshed = false)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const RealVect bus_vn_kv = RealVect::Constant(NB_BUS, 138.);
    grid.init_bus(static_cast<unsigned int>(NB_BUS), 1, bus_vn_kv, 0, 0);

    const int n_line = meshed ? NB_BUS : NB_BUS - 1;
    const RealVect branch_r = RealVect::Constant(n_line, 0.01);
    const RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < NB_BUS - 1; ++i) { from_id(i) = i; to_id(i) = i + 1; }
    if (meshed) { from_id(NB_BUS - 1) = 0; to_id(NB_BUS - 1) = NB_BUS - 1; }
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << LOAD_P;
    load_q << LOAD_Q;
    Eigen::VectorXi load_bus(1);
    load_bus << NB_BUS - 1;
    grid.init_loads(load_p, load_q, load_bus);

    const int nb_gen = static_cast<int>(gens.size());
    RealVect gen_p(nb_gen), gen_v(nb_gen), gen_min_q(nb_gen), gen_max_q(nb_gen);
    Eigen::VectorXi gen_bus(nb_gen);
    for (int k = 0; k < nb_gen; ++k) {
        gen_p(k) = gens[k].p;
        gen_v(k) = gens[k].vset;
        gen_min_q(k) = gens[k].min_q;
        gen_max_q(k) = gens[k].max_q;
        gen_bus(k) = gens[k].bus;
    }
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    for (int k = 0; k < nb_gen; ++k) {
        if (gens[k].regulated_bus >= 0) grid.set_gen_regulated_bus(k, gens[k].regulated_bus);
    }
    grid.tell_solver_need_reset();
    return grid;
}

// the slack generator, with limits wide enough never to be reported
GenSpec slack_gen() { return GenSpec{0, 1.02, 0., -WIDE_Q, WIDE_Q, -1}; }

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

// every generator's reactive output as a single-shot ac_pf publishes it: the reference
// this whole file is written against
RealVect reference_gen_q(LSGrid & grid)
{
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    grid.ac_pf(flat_start(grid), 30, 1e-11);
    RealVect res(std::get<1>(grid.get_gen_res()));
    return res;
}

// the violation reported on `gen_id`, or nullptr
const LimitViolation * find_gen(const std::vector<LimitViolation> & viols, int gen_id)
{
    for (std::size_t k = 0; k < viols.size(); ++k) {
        if (viols[k].element_type == ViolationElementType::GENERATOR && viols[k].element_id == gen_id) {
            return &viols[k];
        }
    }
    return nullptr;
}

// one TimeSeries row whose injection is the grid's own state, with the check on. A batch
// is neither copyable nor movable (and it works on a private copy of the grid taken at
// construction), so it is built by the caller and configured here.
void setup_one_row(TimeSeries & ts)
{
    ts.set_compute_gen_q_violations(true);
    ts.set_gen_q_violation_tol_mvar(0.);
    RealMat load_p(1, 1);
    load_p << LOAD_P;
    ts.modify_load_p(load_p);
    RealMat load_q(1, 1);
    load_q << LOAD_Q;
    ts.modify_load_q(load_q);
}

}  // namespace

TEST_CASE("a batch row reports the reactive output a single ac_pf publishes", "[batch][gen_q]")
{
    SECTION("wide limits: nothing reported")
    {
        std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -WIDE_Q, WIDE_Q, -1}};
        LSGrid grid = make_grid(gens);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);
        REQUIRE(ts.converged_mask()[0] == 1);
        REQUIRE(ts.get_gen_q_violations().size() == 1);
        CHECK(ts.get_gen_q_violations()[0].empty());
        CHECK(ts.get_gen_q_violations_n().empty());
    }

    SECTION("a limit the solved value leaves: reported, with the value ac_pf publishes")
    {
        // +/- 10 MVAr on the machine at bus 1, which holds bus 1 at 1.05 pu against a
        // 80 MW / 60 MVAr load two lines away: it takes far more than 10 MVAr to do it.
        std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
        LSGrid ref_grid = make_grid(gens);
        const RealVect q_ref = reference_gen_q(ref_grid);
        REQUIRE(q_ref.size() == 2);
        REQUIRE(q_ref(1) > 10.);  // the case is only meaningful if the limit IS left

        LSGrid grid = make_grid(gens);
        std::vector<std::string> names{"slack_gen", "gen_bus1"};
        grid.set_gen_names(names);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        ts.compute(flat_start(grid), 30, 1e-11);

        REQUIRE(ts.converged_mask()[0] == 1);
        const std::vector<LimitViolation> & viols = ts.get_gen_q_violations()[0];
        REQUIRE(viols.size() == 1);
        CHECK(viols[0].element_type == ViolationElementType::GENERATOR);
        CHECK(viols[0].element_id == 1);
        CHECK(viols[0].side == 0);
        CHECK(viols[0].violation_type == LimitViolationType::HIGH_Q);
        CHECK(viols[0].value == Approx(q_ref(1)).margin(1e-6));
        CHECK(viols[0].limit == Approx(10.));
        CHECK(viols[0].name == "gen_bus1");
        // the base ("n") case solves that same grid, so it reports the same thing
        REQUIRE(ts.get_gen_q_violations_n().size() == 1);
        CHECK(ts.get_gen_q_violations_n()[0].value == Approx(q_ref(1)).margin(1e-6));
    }

    SECTION("the tolerance is what keeps a machine resting on its limit quiet")
    {
        std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
        LSGrid ref_grid = make_grid(gens);
        const real_type q1 = reference_gen_q(ref_grid)(1);

        LSGrid grid = make_grid(gens);
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        TimeSeries ts(grid);
        setup_one_row(ts);
        // a tolerance wider than the overshoot hides it; anything smaller does not
        ts.set_gen_q_violation_tol_mvar(q1 - 10. + 1.);
        ts.compute(flat_start(grid), 30, 1e-11);
        CHECK(ts.get_gen_q_violations()[0].empty());
    }
}

TEST_CASE("two generators on one bus split the residual the way LSGrid does", "[batch][gen_q]")
{
    // gen 1 and gen 2 both regulate bus 1 with reactive ranges in a 1:2 ratio, so each
    // takes a DIFFERENT share of that bus' reactive residual (LSGrid's share is
    // proportional to qmax - qmin). A batch that split the residual equally, or gave one
    // machine all of it, would fail here and pass every single-generator test.
    std::vector<GenSpec> gens{slack_gen(),
                              GenSpec{1, V_SET, 10., -20., 20., -1},
                              GenSpec{1, V_SET, 10., -40., 40., -1}};
    LSGrid ref_grid = make_grid(gens);
    const RealVect q_ref = reference_gen_q(ref_grid);
    REQUIRE(q_ref.size() == 3);
    REQUIRE(q_ref(1) > 20.);   // both machines leave their range
    REQUIRE(q_ref(2) > 40.);
    // ... and unequally, which is what makes the test discriminating
    REQUIRE(q_ref(2) == Approx(2. * q_ref(1)).epsilon(1e-6));

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);

    const std::vector<LimitViolation> & viols = ts.get_gen_q_violations()[0];
    REQUIRE(viols.size() == 2);
    for (std::size_t k = 0; k < viols.size(); ++k) {
        const int gen_id = viols[k].element_id;
        REQUIRE(gen_id >= 1);
        CHECK(viols[k].violation_type == LimitViolationType::HIGH_Q);
        CHECK(viols[k].value == Approx(q_ref(gen_id)).margin(1e-6));
    }
}

TEST_CASE("a remotely regulating generator's reactive output comes from the controller list",
          "[batch][gen_q][vctrl]")
{
    // gen 1 stands on bus 1 and regulates bus 3: its reactive output is a Jacobian
    // unknown of the VoltageControl extension, and its own bus carries no reactive
    // residual of its own. Reading a residual there would report ~0 instead of the real
    // output, so this fails outright if the controller path is missing.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., NB_BUS - 1}};
    LSGrid ref_grid = make_grid(gens);
    const RealVect q_ref = reference_gen_q(ref_grid);
    REQUIRE(std::abs(q_ref(1)) > 10.);

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.converged_mask()[0] == 1);

    const LimitViolation * viol = find_gen(ts.get_gen_q_violations()[0], 1);
    REQUIRE(viol != nullptr);
    CHECK(viol->value == Approx(q_ref(1)).margin(1e-6));
    CHECK(viol->violation_type == (q_ref(1) > 0. ? LimitViolationType::HIGH_Q
                                                 : LimitViolationType::LOW_Q));
}

TEST_CASE("the generator reactive-limit check is opt in, and says so when it is off",
          "[batch][gen_q]")
{
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    CHECK_FALSE(ts.get_compute_gen_q_violations());
    CHECK_THROWS_AS(ts.get_gen_q_violations(), std::runtime_error);
    CHECK_THROWS_AS(ts.get_gen_q_violations_n(), std::runtime_error);
    CHECK(ts.get_gen_q_violation_tol_mvar() == Approx(1e-4));
    CHECK_THROWS_AS(ts.set_gen_q_violation_tol_mvar(-1.), std::runtime_error);
}

TEST_CASE("the generator reactive-limit check refuses an algorithm that cannot feed it",
          "[batch][gen_q]")
{
    // Every built-in AC family publishes its per-bus mismatch (NR, fast-decoupled AND
    // Gauss-Seidel -- see BaseAlgo::FILLS_BUS_MISMATCH and its overrides), so the only
    // rejections reachable from here are DC and a plugin that does not opt in. Gauss
    // Seidel is checked below to AGREE with the NR reference instead.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};

    // the algorithm a batch runs is its OWN: inherited from the grid at construction,
    // and changed afterwards through the BATCH, never through the grid (see
    // BaseBatchSolverSynch's constructor)
    SECTION("DC: no reactive power at all")
    {
        LSGrid grid = make_grid(gens);
        TimeSeries ts(grid);
        ts.change_algorithm(AlgorithmType::DC_SparseLU);
        ts.set_compute_gen_q_violations(true);
        RealMat load_p(1, 1);
        load_p << LOAD_P;
        ts.modify_load_p(load_p);
        CHECK_THROWS_AS(ts.compute(flat_start(grid), 30, 1e-11), std::runtime_error);
    }
}

TEST_CASE("a second compute() that reuses the base case reports the same thing",
          "[batch][gen_q]")
{
    // `reuse_base_case` (on by default) skips the "n" solve of the second call, which
    // leaves the member algorithm holding the mismatch of the LAST ROW of the first call.
    // Deriving the base case's report from that would report a row as the base case; both
    // calls must report the very same thing.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    ts.set_compute_gen_q_violations(true);
    ts.set_gen_q_violation_tol_mvar(0.);
    // two rows with DIFFERENT loads, so the last row's reactive output is not the base
    // case's and a stale read is visible
    RealMat load_p(2, 1);
    load_p << LOAD_P, 0.5 * LOAD_P;
    ts.modify_load_p(load_p);
    RealMat load_q(2, 1);
    load_q << LOAD_Q, 0.5 * LOAD_Q;
    ts.modify_load_q(load_q);

    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.get_compute_gen_q_violations());
    REQUIRE(ts.get_gen_q_violations_n().size() == 1);
    const real_type q_n_first = ts.get_gen_q_violations_n()[0].value;
    REQUIRE(ts.converged_mask()[1] == 1);
    REQUIRE(ts.get_gen_q_violations().size() == 2);
    REQUIRE(ts.get_gen_q_violations()[1].size() == 1);
    const real_type q_row1_first = ts.get_gen_q_violations()[1][0].value;
    REQUIRE(std::abs(q_row1_first - q_n_first) > 1.);  // the rows differ from the base case

    ts.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(ts.base_case_was_reused());
    REQUIRE(ts.get_gen_q_violations_n().size() == 1);
    CHECK(ts.get_gen_q_violations_n()[0].value == Approx(q_n_first));
    REQUIRE(ts.get_gen_q_violations()[1].size() == 1);
    CHECK(ts.get_gen_q_violations()[1][0].value == Approx(q_row1_first));
}

TEST_CASE("Gauss-Seidel reports the same reactive output as Newton-Raphson", "[batch][gen_q]")
{
    // the reactive output is read off the ALGORITHM's mismatch, so it has to be right for
    // every family that publishes one -- not just for the Newton-Raphson the rest of this
    // file runs on.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
    LSGrid ref_grid = make_grid(gens);
    const real_type q1 = reference_gen_q(ref_grid)(1);

    LSGrid grid = make_grid(gens);
    TimeSeries ts(grid);
    ts.change_algorithm(AlgorithmType::GaussSeidel);
    setup_one_row(ts);
    ts.compute(flat_start(grid), 10000, 1e-9);
    REQUIRE(ts.converged_mask()[0] == 1);
    const LimitViolation * viol = find_gen(ts.get_gen_q_violations()[0], 1);
    REQUIRE(viol != nullptr);
    CHECK(viol->value == Approx(q1).margin(1e-4));
}

TEST_CASE("a contingency row reports its own reactive output, not the base case's",
          "[batch][gen_q][contingency]")
{
    // meshed: line 2 (bus2--bus3) can go without islanding the load, which now reaches
    // bus 3 through line 3 (bus0--bus3). The outage moves the reactive flows, so the base
    // case and the contingency row must report DIFFERENT values -- each matching the
    // ac_pf of the corresponding grid.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
    LSGrid ref_n = make_grid(gens, /*meshed=*/true);
    const RealVect q_n = reference_gen_q(ref_n);
    LSGrid ref_c = make_grid(gens, /*meshed=*/true);
    ref_c.deactivate_powerline(2);
    const RealVect q_c = reference_gen_q(ref_c);
    REQUIRE(q_n(1) > 10.);
    REQUIRE(q_c(1) > 10.);
    REQUIRE(std::abs(q_c(1) - q_n(1)) > 1e-3);  // the outage really does change it

    LSGrid grid = make_grid(gens, /*meshed=*/true);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ContingencyAnalysis ca(grid);
    ca.set_compute_gen_q_violations(true);
    ca.set_gen_q_violation_tol_mvar(0.);
    ca.add_n1(2);
    ca.compute(flat_start(grid), 30, 1e-11);

    REQUIRE(ca.converged_mask()[0] == 1);
    const LimitViolation * viol_n = find_gen(ca.get_gen_q_violations_n(), 1);
    REQUIRE(viol_n != nullptr);
    CHECK(viol_n->value == Approx(q_n(1)).margin(1e-6));
    CHECK(viol_n->limit == Approx(10.));

    REQUIRE(ca.get_gen_q_violations().size() == 1);
    const LimitViolation * viol_c = find_gen(ca.get_gen_q_violations()[0], 1);
    REQUIRE(viol_c != nullptr);
    CHECK(viol_c->value == Approx(q_c(1)).margin(1e-6));
}

TEST_CASE("a row that was never simulated reports nothing at all", "[batch][gen_q][contingency]")
{
    // line 0 (bus0--bus1) carries the whole radial feeder: taking it out islands
    // everything past bus 0, so the row is skipped before the solver. A skipped row must
    // come back empty rather than with a stale or zero-voltage answer.
    std::vector<GenSpec> gens{slack_gen(), GenSpec{1, V_SET, 10., -10., 10., -1}};
    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ContingencyAnalysis ca(grid);
    ca.set_compute_gen_q_violations(true);
    ca.add_n1(0);
    ca.compute(flat_start(grid), 30, 1e-11);

    CHECK(ca.converged_mask()[0] == 0);
    REQUIRE(ca.get_gen_q_violations().size() == 1);
    CHECK(ca.get_gen_q_violations()[0].empty());
    // ... while the base case, which did converge, still reports on its own
    CHECK(find_gen(ca.get_gen_q_violations_n(), 1) != nullptr);
}

TEST_CASE("a generator a row disconnects is never reported, and the other takes its share",
          "[batch][gen_q][scenario_sweep]")
{
    // two machines regulating bus 1; row 0 keeps both, row 1 disconnects gen 2. The one
    // that stays must then take the WHOLE reactive residual of the bus -- which is what a
    // grid with gen 2 deactivated publishes.
    std::vector<GenSpec> gens{slack_gen(),
                              GenSpec{1, V_SET, 10., -20., 20., -1},
                              GenSpec{1, V_SET, 10., -40., 40., -1}};
    LSGrid ref_both = make_grid(gens);
    const RealVect q_both = reference_gen_q(ref_both);
    LSGrid ref_alone = make_grid(gens);
    ref_alone.deactivate_gen(2);
    const RealVect q_alone = reference_gen_q(ref_alone);
    REQUIRE(q_both(1) > 20.);
    REQUIRE(q_alone(1) > 20.);
    REQUIRE(std::abs(q_alone(1) - q_both(1)) > 1e-3);  // gen 1 alone carries more

    LSGrid grid = make_grid(gens);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ScenarioSweep sweep(grid);
    sweep.set_compute_gen_q_violations(true);
    sweep.set_gen_q_violation_tol_mvar(0.);
    RealMat load_p(2, 1);
    load_p << LOAD_P, LOAD_P;
    sweep.modify_load_p(load_p);
    RealMat load_q(2, 1);
    load_q << LOAD_Q, LOAD_Q;
    sweep.modify_load_q(load_q);
    BoolMat gen_off(2, 3);
    gen_off << false, false, false,
               false, false, true;  // row 1 disconnects gen 2
    sweep.set_contingency_gens(gen_off);

    sweep.compute(flat_start(grid), 30, 1e-11);
    REQUIRE(sweep.converged_mask()[0] == 1);
    REQUIRE(sweep.converged_mask()[1] == 1);

    const LimitViolation * row0_g1 = find_gen(sweep.get_gen_q_violations()[0], 1);
    REQUIRE(row0_g1 != nullptr);
    CHECK(row0_g1->value == Approx(q_both(1)).margin(1e-6));

    const LimitViolation * row1_g1 = find_gen(sweep.get_gen_q_violations()[1], 1);
    REQUIRE(row1_g1 != nullptr);
    CHECK(row1_g1->value == Approx(q_alone(1)).margin(1e-6));
    // the disconnected machine produces nothing and is never reported
    CHECK(find_gen(sweep.get_gen_q_violations()[1], 2) == nullptr);
}
