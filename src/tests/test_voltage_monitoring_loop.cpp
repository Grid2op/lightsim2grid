// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's VoltageMonitoring outer loop (VoltageMonitoringLoop): an idle standby SVC
// (a voltage monitor) held in its voltage-control group at Q = 0, switched on at the
// automaton's low / high set-point when the voltage it monitors leaves the thresholds; and
// the automaton's b0, a shunt carried by the SVC. C++14 only.

#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/ReactiveLimitsLoop.hpp"
#include "powerflow_algorithm/outer_loop/VoltageMonitoringLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::real_type;
using ls2g::VoltageMonitoringLoop;

namespace {

const real_type NAN_ = std::numeric_limits<real_type>::quiet_NaN();

// buses 0-1-2-3 in a line, a heavy load at bus 3 pulling the voltages down, the slack
// generator at bus 0 (1.0 pu), an SVC at bus 2 regulating `reg_bus`, of mode `mode`
// (0 off, 1 voltage) at `target_vm` pu
LSGrid make_grid(int mode = 0, real_type target_vm = 1., int reg_bus = 2, real_type b_range = 5.)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
    Eigen::VectorXi from_id(3), to_id(3);
    from_id << 0, 1, 2;
    to_id << 1, 2, 3;
    grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.08),
                         CplxVect::Zero(3), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 80.;
    load_q << 40.;
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);

    RealVect target(1), q(1), slope(1), b_min(1), b_max(1);
    target << target_vm;
    q << 0.;
    slope << 0.;
    b_min << -b_range;
    b_max << b_range;
    Eigen::VectorXi reg(1), bus(1);
    reg << reg_bus;
    bus << 2;
    grid.init_svcs({mode}, target, q, slope, b_min, b_max, reg, bus);
    return grid;
}

// flag the SVC standby, thresholds [low, high], set-points 0.97 / 1.03 pu
void set_monitor(LSGrid & grid, real_type low, real_type high)
{
    RealVect lo(1), hi(1), lo_t(1), hi_t(1);
    lo << low;
    hi << high;
    lo_t << 0.97;
    hi_t << 1.03;
    grid.set_svc_standby({true}, lo, hi, lo_t, hi_t);
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

CplxVect solve_outer(LSGrid & grid)
{
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<VoltageMonitoringLoop>());
    return grid.ac_pf(flat_start(grid), 30, 1e-10);
}

real_type svc_q(const LSGrid & grid)
{
    return std::get<1>(grid.get_svcs().get_res())(0);
}

}  // namespace

TEST_CASE("a monitor below its low threshold is switched on at the low set-point", "[outer][svc]")
{
    LSGrid grid = make_grid();
    set_monitor(grid, 0.99, 1.10);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    REQUIRE(stats.loop_iterations.size() == 1);
    CHECK(stats.loop_iterations[0].second == 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    CHECK(std::abs(V(2)) == Approx(0.97).margin(1e-9));
    CHECK(svc_q(grid) > 1.);  // it produces to hold the bus up
    // the SVC itself is untouched: still off, still standby
    CHECK(grid.get_svcs().get_regulation_mode(0) == 0);

    // the same as the SVC regulating at that set-point from the start
    LSGrid ref = make_grid(/*mode=*/1, /*target_vm=*/0.97);
    ref.change_algorithm("NRSing_SparseLU");
    const CplxVect V_ref = ref.ac_pf(flat_start(ref), 30, 1e-10);
    REQUIRE(V_ref.size() == 4);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);
    CHECK(svc_q(grid) == Approx(svc_q(ref)).margin(1e-6));

    // a next solve starts again from the idle SVC
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().loop_iterations[0].second == 1);
}

TEST_CASE("a monitor above its high threshold is switched on at the high set-point", "[outer][svc]")
{
    LSGrid grid = make_grid();
    set_monitor(grid, 0.50, 0.60);  // the idle voltage is well above both
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(std::abs(V(2)) == Approx(1.03).margin(1e-9));
}

TEST_CASE("a monitor inside its thresholds stays idle", "[outer][svc]")
{
    LSGrid grid = make_grid();
    set_monitor(grid, 0.80, 1.10);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().nb_outer_iterations == 0);
    CHECK(svc_q(grid) == Approx(0.).margin(1e-9));
    // held at 0 is the grid with the SVC off
    LSGrid ref = make_grid();
    ref.change_algorithm("NRSing_SparseLU");
    const CplxVect V_ref = ref.ac_pf(flat_start(ref), 30, 1e-10);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);
}

TEST_CASE("only a monitor of its own bus makes the loop run", "[outer][svc]")
{
    LSGrid grid = make_grid(/*mode=*/0, /*target_vm=*/1., /*reg_bus=*/3);
    set_monitor(grid, 0.99, 1.10);
    REQUIRE(solve_outer(grid).size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().loop_iterations.empty());
}

TEST_CASE("detection: a plain solve reports the switch", "[outer][svc]")
{
    LSGrid grid = make_grid();
    set_monitor(grid, 0.99, 1.10);
    grid.change_algorithm("NRSing_SparseLU");
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 4);
    int found = 0;
    for (const LimitViolation & v : grid.get_physical_violations()) {
        if (v.violation_type == LimitViolationType::LOW_VOLTAGE_SVC_STANDBY) ++found;
    }
    CHECK(found == 1);
}

TEST_CASE("b0 is a shunt carried by the SVC", "[svc]")
{
    // the SVC off with b0, against the same grid with that susceptance as a shunt
    LSGrid grid = make_grid();
    RealVect b0(1);
    b0 << 0.2;  // pu
    grid.set_svc_b0(b0);
    grid.change_algorithm("NRSing_SparseLU");
    const CplxVect V = grid.ac_pf(flat_start(grid), 30, 1e-10);
    REQUIRE(V.size() == 4);

    LSGrid ref = make_grid();
    RealVect sh_p(1), sh_q(1);
    sh_p << 0.;
    sh_q << -20.;  // MVar consumed at 1 pu: -b0 * sn_mva
    Eigen::VectorXi sh_bus(1);
    sh_bus << 2;
    ref.init_shunt(sh_p, sh_q, sh_bus);
    ref.change_algorithm("NRSing_SparseLU");
    const CplxVect V_ref = ref.ac_pf(flat_start(ref), 30, 1e-10);
    REQUIRE(V_ref.size() == 4);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-10);
    // its output is the SVC's, b0.V^2 produced
    CHECK(svc_q(grid) == Approx(20. * std::norm(V(2))).margin(1e-8));
}

TEST_CASE("a monitor switched on is held at its limit by ReactiveLimits", "[outer][svc][reactive_limits]")
{
    // switched on at 0.97 pu, but a 0.05 pu susceptance cannot hold it there
    LSGrid grid = make_grid(0, 1., 2, /*b_range=*/0.05);
    set_monitor(grid, 0.99, 1.10);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<VoltageMonitoringLoop>());
    grid.add_outer_loop(std::make_shared<ls2g::ReactiveLimitsLoop>());
    const CplxVect V = grid.ac_pf(flat_start(grid), 30, 1e-10);
    REQUIRE(V.size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    CHECK(std::abs(V(2)) < 0.97);
    // held at b_max V², V that of the solve which switched it (the value is frozen)
    CHECK(svc_q(grid) == Approx(0.05 * std::norm(V(2)) * 100.).margin(1e-3));
}
