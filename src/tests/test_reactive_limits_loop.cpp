// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's ReactiveLimits outer loop (ReactiveLimitsLoop) on buses holding their own
// voltage: frozen PQ at a reactive limit, the slack bus included, the strongest PV bus kept,
// on one symbolic analysis. C++14 only.

#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/ReactiveLimitsLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::ReactiveLimitsLoop;
using ls2g::real_type;

namespace {

// buses 0-1-2-3 in a line, a load at bus 3; the slack generator at bus 0 (1.02 pu) and a
// generator at bus 2 (`v2` pu), with the given reactive ranges
LSGrid make_grid(real_type v2, real_type min_q0, real_type max_q0, real_type min_q2, real_type max_q2,
                 bool regulating2 = true, real_type q2 = 0.)
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
    load_p << 60.;
    load_q << 40.;
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(2), gen_v(2), q(2), gen_min_q(2), gen_max_q(2);
    gen_p << 0., 20.;
    gen_v << 1.02, v2;
    q << 0., q2;
    gen_min_q << min_q0, min_q2;
    gen_max_q << max_q0, max_q2;
    Eigen::VectorXi gen_bus(2);
    gen_bus << 0, 2;
    grid.init_generators_full(gen_p, gen_v, q, {true, regulating2}, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

CplxVect solve_outer(LSGrid & grid)
{
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ReactiveLimitsLoop>());
    return grid.ac_pf(flat_start(grid), 30, 1e-10);
}

CplxVect solve_plain(LSGrid & grid)
{
    grid.change_algorithm("NRSing_SparseLU");
    return grid.ac_pf(flat_start(grid), 30, 1e-10);
}

real_type gen_q(const LSGrid & grid, int gen_id)
{
    return std::get<1>(grid.get_gen_res())(gen_id);
}

}  // namespace

TEST_CASE("a unit beyond its max is frozen there", "[outer][reactive_limits]")
{
    // gen 2 would have to produce more than 10 MVar to hold 1.05 pu
    LSGrid grid = make_grid(1.05, -500., 500., -10., 10.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    REQUIRE(stats.loop_iterations.size() == 1);
    CHECK(stats.loop_iterations[0].second == 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    CHECK(gen_q(grid, 1) == Approx(10.).margin(1e-6));
    CHECK(std::abs(V(2)) < 1.05);

    // the same as the unit injecting 10 MVar
    LSGrid ref = make_grid(1.05, -500., 500., -10., 10., /*regulating2=*/false, /*q2=*/10.);
    const CplxVect V_ref = solve_plain(ref);
    REQUIRE(V_ref.size() == 4);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);

    // a next solve starts again from the PV bus
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().loop_iterations[0].second == 1);
}

TEST_CASE("a unit below its min is frozen there", "[outer][reactive_limits]")
{
    // holding 0.85 pu, below what the bus would sit at, gen 2 would have to absorb more
    // than 5 MVar
    LSGrid grid = make_grid(0.85, -500., 500., -5., 100.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(gen_q(grid, 1) == Approx(-5.).margin(1e-6));
    CHECK(std::abs(V(2)) > 0.85);
}

TEST_CASE("inside its limits nothing moves", "[outer][reactive_limits]")
{
    LSGrid grid = make_grid(1.02, -500., 500., -500., 500.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().nb_outer_iterations == 0);
    LSGrid ref = make_grid(1.02, -500., 500., -500., 500.);
    CHECK((V - solve_plain(ref)).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("the slack bus is frozen like any other", "[outer][reactive_limits]")
{
    // the slack unit cannot produce what holding 1.02 pu takes; gen 2 can
    LSGrid grid = make_grid(1.0, -5., 5., -500., 500.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(gen_q(grid, 0) == Approx(5.).margin(1e-6));
    CHECK(std::abs(V(0)) != Approx(1.02).margin(1e-6));
    CHECK(std::arg(V(0)) == Approx(0.).margin(1e-12));  // still the angle reference
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
}

TEST_CASE("the strongest bus stays PV when every one would switch", "[outer][reactive_limits]")
{
    // neither unit can hold its bus: the one with the highest target P (gen 2) stays PV
    LSGrid grid = make_grid(1.05, -1., 1., -1., 1.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 4);
    CHECK(gen_q(grid, 0) == Approx(-1.).margin(1e-6));  // the slack unit has to absorb
    CHECK(std::abs(V(2)) == Approx(1.05).margin(1e-9));
}

namespace {

// buses 0-1-2-3-4 in a line, a load at bus 4; the slack generator at bus 0 (1.02 pu), and
// generators at buses 1 and 3 regulating bus 2 remotely at 1.04 pu (when `regulating`),
// with the given reactive ranges, or injecting `q` when not
LSGrid make_group_grid(real_type max_q1, real_type max_q3, bool regulating1 = true, real_type q1 = 0.,
                       bool with_gen3 = true, real_type min_q1 = std::numeric_limits<real_type>::quiet_NaN())
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(5, 1, RealVect::Constant(5, 138.), 0, 0);
    Eigen::VectorXi from_id(4), to_id(4);
    from_id << 0, 1, 2, 3;
    to_id << 1, 2, 3, 4;
    grid.init_powerlines(RealVect::Constant(4, 0.01), RealVect::Constant(4, 0.08),
                         CplxVect::Zero(4), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 60.;
    load_q << 50.;
    Eigen::VectorXi load_bus(1);
    load_bus << 4;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(3), gen_v(3), q(3), gen_min_q(3), gen_max_q(3);
    gen_p << 0., 10., 10.;
    gen_v << 1.02, 1.04, 1.04;
    q << 0., q1, 0.;
    gen_min_q << -500., std::isnan(min_q1) ? -max_q1 : min_q1, -max_q3;
    gen_max_q << 500., max_q1, max_q3;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 3;
    grid.init_generators_full(gen_p, gen_v, q, {true, regulating1, with_gen3}, gen_min_q, gen_max_q, gen_bus);
    grid.set_gen_regulated_bus(1, 2);
    grid.set_gen_regulated_bus(2, 2);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

}  // namespace

TEST_CASE("a remote controller beyond its max is held there", "[outer][reactive_limits]")
{
    LSGrid grid = make_group_grid(/*max_q1=*/5., /*max_q3=*/500., true, 0., /*with_gen3=*/false);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 5);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    CHECK(gen_q(grid, 1) == Approx(5.).margin(1e-6));
    CHECK(std::abs(V(2)) < 1.04);

    LSGrid ref = make_group_grid(5., 500., /*regulating1=*/false, /*q1=*/5., false);
    const CplxVect V_ref = solve_plain(ref);
    REQUIRE(V_ref.size() == 5);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);
}

TEST_CASE("the first controller of a group is held, the other one regulates on", "[outer][reactive_limits]")
{
    // the group shares by reactive range: gen 1 (its first controller) has a wide one, so it
    // is asked a large share, but it can produce 5 MVar only
    LSGrid grid = make_group_grid(/*max_q1=*/5., /*max_q3=*/500., true, 0., true, /*min_q1=*/-500.);
    const CplxVect V = solve_outer(grid);
    REQUIRE(V.size() == 5);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    CHECK(gen_q(grid, 1) == Approx(5.).margin(1e-6));
    CHECK(std::abs(V(2)) == Approx(1.04).margin(1e-9));  // gen 3 still holds it

    // the same as gen 1 injecting 5 MVar, gen 3 regulating alone
    LSGrid ref = make_group_grid(5., 500., /*regulating1=*/false, /*q1=*/5., true, -500.);
    const CplxVect V_ref = solve_plain(ref);
    REQUIRE(V_ref.size() == 5);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);
    CHECK(gen_q(grid, 2) == Approx(gen_q(ref, 2)).margin(1e-6));
}
