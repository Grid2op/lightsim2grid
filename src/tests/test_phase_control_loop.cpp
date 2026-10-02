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

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/PhaseControlLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
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
