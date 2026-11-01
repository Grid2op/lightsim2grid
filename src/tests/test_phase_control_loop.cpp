// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's PhaseControl outer loop (PhaseControlLoop + the BranchControl extension): an
// active power controller solved for its shift then rounded, a current limiter moved tap by
// tap, on one symbolic analysis, the inputs untouched. C++14 only.

#include <cmath>
#include <memory>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "powerflow_algorithm/outer_loop/PhaseControlLoop.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::OuterLoopStatus;
using ls2g::PhaseControlLoop;
using ls2g::RealVect;
using ls2g::RegulationMode;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// a triangle 0-1-2: lines 0-1 and 1-2, a transformer 0-2 (the phase shifter), loads at 1 and
// 2, the slack generator at 0. The phase tap changer runs from -10 to 10 degrees in `n_steps`
// positions, its reactance growing by `x_pct_per_deg` % per degree, at 0 degree.
LSGrid make_grid(int n_steps, real_type x_pct_per_deg = 0.5)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi line_from(2), line_to(2);
    line_from << 0, 1;
    line_to << 1, 2;
    grid.init_powerlines(RealVect::Constant(2, 0.01), RealVect::Constant(2, 0.08), CplxVect::Zero(2), line_from, line_to);
    RealVect trafo_r(1), trafo_x(1), trafo_ratio(1), trafo_shift(1);
    trafo_r << 0.005;
    trafo_x << 0.06;
    trafo_ratio << 1.;
    trafo_shift << 0.;
    Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
    trafo_bus1 << 0;
    trafo_bus2 << 2;
    grid.init_trafo(trafo_r, trafo_x, CplxVect::Zero(1), trafo_ratio, trafo_shift, {false}, trafo_bus1, trafo_bus2, false);
    RealVect load_p(2), load_q(2);
    load_p << 60., 90.;
    load_q << 20., 30.;
    Eigen::VectorXi load_bus(2);
    load_bus << 1, 2;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);

    std::vector<real_type> rho, alpha, r, x, g, b;
    for (int k = 0; k < n_steps; ++k) {
        const real_type a = -10. + 20. * k / (n_steps - 1);
        rho.push_back(1.);
        alpha.push_back(a);
        r.push_back(0.);
        x.push_back(x_pct_per_deg * a);
        g.push_back(0.);
        b.push_back(0.);
    }
    grid.set_trafo_phase_tap_changer(0, 0, (n_steps - 1) / 2, rho, alpha, r, x, g, b);
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

void use_phase_control(LSGrid & grid)
{
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<PhaseControlLoop>());
}

}  // namespace

TEST_CASE("an active power phase shifter is solved for its shift, then rounded", "[outer][phase_control]")
{
    // its flow without control. The reactance does not move with the tap here: OpenLoadFlow
    // solves the shift with the starting tap's, so a rounded tap of another reactance lands
    // off the target (test_phase_control.py compares that case with OpenLoadFlow)
    LSGrid ref = make_grid(2001, 0.);
    ref.change_algorithm("NRSing_SparseLU");
    REQUIRE(ref.ac_pf(flat(ref), 30, 1e-12).size() == 3);
    const real_type p0 = std::get<0>(ref.get_trafo_res1())(0);

    // a target 5 MW above, on a fine table: rounded to within a tap's worth of it
    LSGrid grid = make_grid(2001, 0.);
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
    use_phase_control(grid);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    REQUIRE(stats.loop_iterations.size() == 1);
    CHECK(stats.loop_iterations[0].second == 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    const ls2g::TrafoInfo pst = grid.get_trafos()[0];
    CHECK(pst.phase_tap_position == 1000);  // the input stays
    CHECK(pst.res_phase_tap_position != 1000);
    CHECK(pst.res_p1_mw == Approx(p0 + 5.).margin(0.05));

    // the same as the transformer moved to that tap
    LSGrid moved = make_grid(2001, 0.);
    moved.change_trafo_phase_tap(0, pst.res_phase_tap_position);
    moved.change_algorithm("NRSing_SparseLU");
    REQUIRE(moved.ac_pf(flat(moved), 30, 1e-12).size() == 3);
    CHECK(moved.get_trafos()[0].res_p1_mw == Approx(pst.res_p1_mw).margin(1e-9));
}

TEST_CASE("a current limiter moves tap by tap until below its limit", "[outer][phase_control]")
{
    LSGrid ref = make_grid(21);
    ref.change_algorithm("NRSing_SparseLU");
    REQUIRE(ref.ac_pf(flat(ref), 30, 1e-12).size() == 3);
    const real_type i0 = std::get<3>(ref.get_trafo_res1())(0) * 1000.;  // A

    LSGrid grid = make_grid(21);
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::CURRENT_LIMITER, true, 0.85 * i0, 0., 1);
    use_phase_control(grid);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    const ls2g::TrafoInfo pst = grid.get_trafos()[0];
    CHECK(pst.res_a1_ka * 1000. < 0.85 * i0);
    CHECK(pst.res_phase_tap_position != 10);

    // one tap back, it was above
    LSGrid back = make_grid(21);
    const int step = pst.res_phase_tap_position < 10 ? 1 : -1;
    back.change_trafo_phase_tap(0, pst.res_phase_tap_position + step);
    back.change_algorithm("NRSing_SparseLU");
    REQUIRE(back.ac_pf(flat(back), 30, 1e-12).size() == 3);
    CHECK(back.get_trafos()[0].res_a1_ka * 1000. > 0.85 * i0);
}

TEST_CASE("a phase shifter is reported where the loop would act", "[outer][phase_control]")
{
    LSGrid grid = make_grid(21);
    grid.change_algorithm("NRSing_SparseLU");
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    const real_type p0 = std::get<0>(grid.get_trafo_res1())(0);
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 1., 1);
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<PhaseControlLoop>());
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    std::vector<ls2g::LimitViolation> found;
    for (const auto & v : grid.get_physical_violations(true, 1e-3, 0.)) {
        if (v.violation_type == LimitViolationType::PHASE_CONTROL_P) found.push_back(v);
    }
    REQUIRE(found.size() == 1);
    CHECK(found[0].value == Approx(p0).margin(1e-9));
    CHECK(found[0].limit == p0 + 5.);
    CHECK(found[0].category() == ls2g::ViolationCategory::CONTROL);
    // within its deadband: nothing
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 6., 1);
    found.clear();
    for (const auto & v : grid.get_physical_violations(true, 1e-3, 0.)) {
        if (v.violation_type == LimitViolationType::PHASE_CONTROL_P) found.push_back(v);
    }
    CHECK(found.empty());
}

TEST_CASE("a phase regulation set between two solves is the one the next solve holds", "[outer][phase_control]")
{
    // what the outer loops reserve in the Newton is built from the regulation data: a phase
    // shifter that starts regulating, or a new target, must reach the next solve
    LSGrid ref = make_grid(2001, 0.);
    ref.change_algorithm("NRSing_SparseLU");
    REQUIRE(ref.ac_pf(flat(ref), 30, 1e-12).size() == 3);
    const real_type p0 = std::get<0>(ref.get_trafo_res1())(0);

    LSGrid expected = make_grid(2001, 0.);
    expected.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
    use_phase_control(expected);
    const CplxVect V_expected = expected.ac_pf(flat(expected), 30, 1e-12);
    REQUIRE(V_expected.size() == 3);
    const int pos_expected = expected.get_trafos()[0].res_phase_tap_position;
    REQUIRE(pos_expected != 1000);

    SECTION("a phase shifter that starts regulating") {
        LSGrid grid = make_grid(2001, 0.);
        use_phase_control(grid);
        REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
        const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
        REQUIRE(V.size() == 3);
        CHECK(grid.get_trafos()[0].res_phase_tap_position == pos_expected);
        CHECK((V - V_expected).cwiseAbs().maxCoeff() < 1e-9);
    }
    SECTION("a new target") {
        LSGrid grid = make_grid(2001, 0.);
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 - 5., 0., 1);
        use_phase_control(grid);
        REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
        const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
        REQUIRE(V.size() == 3);
        CHECK(grid.get_trafos()[0].res_phase_tap_position == pos_expected);
        CHECK((V - V_expected).cwiseAbs().maxCoeff() < 1e-9);
    }
}

namespace {

int nb_phase_reports(const std::vector<LimitViolation> & out)
{
    int nb = 0;
    for (const LimitViolation & v : out) {
        if (v.violation_type != LimitViolationType::PHASE_CONTROL_P &&
            v.violation_type != LimitViolationType::PHASE_LIMITER_CURRENT) continue;
        INFO("reported: " << v.value << " against " << v.limit);
        ++nb;
    }
    return nb;
}

}  // namespace

TEST_CASE("the report after a solve reads the flow at the tap the loop moved to", "[outer][phase_control]")
{
    LSGrid ref = make_grid(2001, 0.);
    ref.change_algorithm("NRSing_SparseLU");
    REQUIRE(ref.ac_pf(flat(ref), 30, 1e-12).size() == 3);
    const real_type p0 = std::get<0>(ref.get_trafo_res1())(0);

    // the loop meets its target to within a tap's worth: a solve that did that has nothing
    // to report, at the tap it ends on (the grid's own tap is 5 MW away)
    LSGrid grid = make_grid(2001, 0.);
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
    use_phase_control(grid);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    REQUIRE(grid.get_trafos()[0].res_phase_tap_position != 1000);
    REQUIRE(grid.get_trafos()[0].res_p1_mw == Approx(p0 + 5.).margin(0.05));

    SECTION("the grid's report") {
        CHECK(nb_phase_reports(grid.get_physical_violations(true, 0.1, 1e-4)) == 0);
    }
    SECTION("a batch row's, and its \"n\" report") {
        ls2g::ContingencyAnalysis ca(grid);
        ca.set_compute_physical_violations(true);
        ca.set_physical_violation_tol_mva(0.1);
        ca.add_nk(std::vector<int>());
        ca.compute(flat(grid), 30, 1e-12);
        REQUIRE(ca.converged_mask()[0]);
        CHECK(nb_phase_reports(ca.get_physical_violations()[0]) == 0);
        CHECK(nb_phase_reports(ca.get_physical_violations_n()) == 0);
    }
}

TEST_CASE("the taps a solve left on transformers replaced since are refused, not read", "[outer][phase_control]")
{
    // the checks after a solve read the taps its loops left, one per transformer of that
    // solve; init_trafo then replaces the transformers (none here). They used to read those
    // taps in the tap tables init_trafo had freed (a heap-use-after-free, then an id check
    // raising after the read): the message shows the readback is now refused before any read
    LSGrid ref = make_grid(21);
    ref.change_algorithm("NRSing_SparseLU");
    REQUIRE(ref.ac_pf(flat(ref), 30, 1e-12).size() == 3);
    const real_type p0 = std::get<0>(ref.get_trafo_res1())(0);

    LSGrid grid = make_grid(21);
    grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p0 + 5., 0., 1);
    use_phase_control(grid);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    std::vector<int> taps;
    grid.get_algo().get_outer_phase_tap(taps);
    REQUIRE(taps.size() == 1);
    CHECK_NOTHROW(grid.get_violations(1., true, 0.));
    CHECK_NOTHROW(grid.get_physical_violations(true, 1e-4, 1e-4));

    grid.init_trafo(RealVect(0), RealVect(0), CplxVect(0), RealVect(0), RealVect(0), std::vector<bool>(),
                    Eigen::VectorXi(0), Eigen::VectorXi(0), false);
    REQUIRE(grid.get_trafos().nb() == 0);
    CHECK_THROWS_MATCHES(grid.get_violations(1., true, 0.), std::out_of_range,
                         MessageMatches(ContainsSubstring("were the transformers replaced")));
    CHECK_THROWS_MATCHES(grid.get_physical_violations(true, 1e-4, 1e-4), std::out_of_range,
                         MessageMatches(ContainsSubstring("were the transformers replaced")));
}

TEST_CASE("after a DC solve, a phase shifter's flow is the DC one", "[outer][phase_control]")
{
    // the DC flow through the phase shifter at its tap
    LSGrid ref = make_grid(21);
    REQUIRE(ref.dc_pf(flat(ref), 1, 1e-12).size() == 3);
    const real_type p_dc = std::get<0>(ref.get_trafo_res1())(0);

    LSGrid grid = make_grid(21);
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<PhaseControlLoop>());
    SECTION("at its target: nothing") {
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p_dc, 0., 1);
        REQUIRE(grid.dc_pf(flat(grid), 1, 1e-12).size() == 3);
        CHECK(nb_phase_reports(grid.get_physical_violations(false, 1e-3, 1e-4)) == 0);
    }
    SECTION("away from it: the DC flow") {
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::ACTIVE_POWER, true, p_dc + 5., 0., 1);
        REQUIRE(grid.dc_pf(flat(grid), 1, 1e-12).size() == 3);
        const std::vector<LimitViolation> out = grid.get_physical_violations(false, 1e-3, 1e-4);
        REQUIRE(nb_phase_reports(out) == 1);
        for (const LimitViolation & v : out) {
            if (v.violation_type != LimitViolationType::PHASE_CONTROL_P) continue;
            CHECK(v.value == Approx(p_dc).margin(1e-9));
        }
    }
    SECTION("a current limiter: no current in DC") {
        grid.set_trafo_phase_tap_regulation(0, RegulationMode::CURRENT_LIMITER, true, 1., 0., 1);
        REQUIRE(grid.dc_pf(flat(grid), 1, 1e-12).size() == 3);
        CHECK(nb_phase_reports(grid.get_physical_violations(false, 1e-3, 1e-4)) == 0);
    }
}
