// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OuterRowView / OuterGridView: a batch row as the outer loops see it -- the grid's data with
// the row's overrides on top -- and a loop reading a row through it exactly as it reads a grid
// where the row's change is real. C++14 only.

#include <cmath>
#include <limits>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/BaseOuterLoop.hpp"
#include "powerflow_algorithm/outer_loop/DistributedSlackLoop.hpp"
#include "powerflow_algorithm/outer_loop/OuterRowView.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::DistributedSlackLoop;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::OuterContext;
using ls2g::OuterGridView;
using ls2g::OuterRowView;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

const real_type NaN = std::numeric_limits<real_type>::quiet_NaN();

// buses 0-1-2-3 in a line, a 200 MW load at bus 3; generators at buses 0 (the slack, 50 MW),
// 1 (40 MW) and 2 (60 MW), sharing the slack with the given weights
LSGrid make_grid(const std::vector<real_type> & weights)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
    Eigen::VectorXi from_id(3), to_id(3);
    from_id << 0, 1, 2;
    to_id << 1, 2, 3;
    grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.05),
                         CplxVect::Zero(3), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 200.;
    load_q << 20.;
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    gen_p << 50., 40., 60.;
    gen_v << 1.02, 1.01, 1.01;
    gen_min_q << -500., -500., -500.;
    gen_max_q << 500., 500., 500.;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 2;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    RealVect w(3);
    for (int i = 0; i < 3; ++i) w(i) = weights[i];
    grid.set_gen_can_participate_slack({weights[0] > 0., weights[1] > 0., weights[2] > 0.}, w);
    grid.set_gen_p_limits(RealVect::Zero(3), RealVect::Constant(3, 1000.));
    grid.change_algorithm("NRSing_SparseLU");
    return grid;
}

CplxVect solve(LSGrid & grid)
{
    return grid.ac_pf(CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.)), 30, 1e-10);
}

int solver_bus(const LSGrid & grid, int bus_me)
{
    return grid.id_me_to_ac_solver()[bus_me].cast_int();
}

// what a loop detects after the grid's last AC solve, the row `row` shows (null: the grid)
std::vector<LimitViolation> detect(const LSGrid & grid, const ls2g::BaseOuterLoop & loop, const OuterRowView * row)
{
    const CplxVect V = grid.get_algo().get_V();
    const RealVect Va = grid.get_algo().get_Va();
    const CplxVect mismatch = grid.get_algo().get_bus_mismatch();
    const RealVect ctrl_q = grid.get_algo().get_controller_q();
    OuterRowView as_solved;
    if (row == nullptr) row = &as_solved;
    OuterContext ctx;
    ctx.grid = &grid;
    ctx.V = &V;
    ctx.Va = &Va;
    ctx.bus_mismatch = &mismatch;
    ctx.controller_q = &ctrl_q;
    ctx.slack_bus = grid.get_slack_ids_solver()[0].cast_int();
    ctx.row = row;
    std::vector<LimitViolation> out;
    loop.detect(ctx, out);
    return out;
}

int count(const std::vector<LimitViolation> & out, LimitViolationType type)
{
    int res = 0;
    for (const LimitViolation & v : out) res += v.violation_type == type ? 1 : 0;
    return res;
}

}  // namespace

TEST_CASE("without a row, the view is the grid", "[outer_loop][row_view]")
{
    LSGrid grid = make_grid({1., 1., 1.});
    REQUIRE(solve(grid).size() == 4);
    const OuterGridView view(grid, nullptr);
    CHECK(&view.layout() == &grid.get_ac_layout());
    CHECK(view.masked() == nullptr);
    for (int g = 0; g < 3; ++g) {
        CHECK(view.gen_on(g));
        CHECK(view.gen_bus_solver(g) == solver_bus(grid, g));
        CHECK(view.gen_bus_grid(g) == g);
        CHECK(view.gen_is_local_voltage_controller(g));
    }
    CHECK(view.gen_target_p(1) == Approx(40.));
    CHECK(view.gen_target_vm(2) == Approx(1.01));
    CHECK(view.load_target_q(0) == Approx(20.));
    CHECK(view.load_bus_solver(0) == solver_bus(grid, 3));
    int b1, b2;
    REQUIRE(view.line_closed(1, b1, b2));
    CHECK(b1 == solver_bus(grid, 1));
    CHECK(b2 == solver_bus(grid, 2));
    const std::vector<real_type> loops = {NaN, 45., NaN};
    std::vector<real_type> scratch;
    CHECK(&view.gen_target_p_overrides(loops, scratch) == &loops);  // nothing to merge
}

TEST_CASE("a row's overrides win over the grid, element by element", "[outer_loop][row_view]")
{
    LSGrid grid = make_grid({1., 1., 1.});
    REQUIRE(solve(grid).size() == 4);
    OuterRowView row;
    std::vector<int> masked = {solver_bus(grid, 3)};
    row.masked = &masked;
    row.gens.off = {0, 1, 0};                           // generator 1 disconnected
    row.gens.bus_solver = {OuterRowView::Units::GRID_BUS, OuterRowView::Units::GRID_BUS,
                           solver_bus(grid, 3)};        // generator 2 moved to bus 3
    row.gens.target_p = {NaN, NaN, 70.};
    row.loads.target_q = {35.};
    row.lines.off = {0, 1, 0};                          // line 1 open
    row.lines.bus1_solver = {OuterRowView::Units::GRID_BUS, OuterRowView::Units::GRID_BUS, -1};  // line 2 open on side 1
    const OuterGridView view(grid, &row);

    CHECK(&view.layout() == &grid.get_ac_layout());     // no layout of its own: the grid's
    CHECK(view.is_masked(solver_bus(grid, 3)));
    CHECK_FALSE(view.is_masked(solver_bus(grid, 0)));

    CHECK(view.gen_on(0));
    CHECK_FALSE(view.gen_on(1));
    CHECK(view.gen_bus_solver(1) == -1);
    CHECK(view.gen_bus_grid(1) == -1);
    CHECK_FALSE(view.gen_is_local_voltage_controller(1));
    CHECK(view.gen_bus_solver(2) == solver_bus(grid, 3));
    CHECK(view.gen_bus_grid(2) == 3);
    CHECK(view.gen_target_p(2) == Approx(70.));
    CHECK(view.gen_target_p(0) == Approx(50.));        // NaN: the grid's
    CHECK(view.load_target_q(0) == Approx(35.));

    int b1, b2;
    CHECK(view.line_closed(0, b1, b2));
    CHECK_FALSE(view.line_closed(1, b1, b2));
    CHECK_FALSE(view.line_closed(2, b1, b2));

    // the loops' moves win over the row's
    const std::vector<real_type> loops = {55., NaN, NaN};
    std::vector<real_type> scratch;
    const std::vector<real_type> & merged = view.gen_target_p_overrides(loops, scratch);
    REQUIRE(merged.size() == 3);
    CHECK(merged[0] == Approx(55.));
    CHECK(std::isnan(merged[1]));
    CHECK(merged[2] == Approx(70.));
}

TEST_CASE("a placed generator is on in its row, even where the grid has it off", "[outer_loop][row_view]")
{
    LSGrid grid = make_grid({1., 1., 1.});
    REQUIRE(solve(grid).size() == 4);
    grid.deactivate_gen(1);
    OuterRowView row;
    row.gens.bus_solver = {OuterRowView::Units::GRID_BUS, solver_bus(grid, 2), OuterRowView::Units::GRID_BUS};
    const OuterGridView view(grid, &row);
    CHECK(view.gen_on(1));
    CHECK(view.gen_bus_solver(1) == solver_bus(grid, 2));
    CHECK(view.gen_is_local_voltage_controller(1));
}

TEST_CASE("a loop detects on a row as on the grid where the row's change is real", "[outer_loop][row_view]")
{
    // generator 1 is the only unit off the slack bus sharing the slack: without it, the slack
    // bus' own unit would take what it already absorbs, and detection says nothing
    DistributedSlackLoop loop;
    {
        LSGrid grid = make_grid({1., 1., 0.});
        REQUIRE(solve(grid).size() == 4);
        CHECK(count(detect(grid, loop, nullptr), LimitViolationType::SLACK_MISMATCH) == 1);

        OuterRowView row;
        row.gens.off = {0, 1, 0};
        CHECK(count(detect(grid, loop, &row), LimitViolationType::SLACK_MISMATCH) == 0);
    }
    {
        LSGrid grid = make_grid({1., 1., 0.});
        grid.deactivate_gen(1);
        REQUIRE(solve(grid).size() == 4);
        CHECK(count(detect(grid, loop, nullptr), LimitViolationType::SLACK_MISMATCH) == 0);
    }
}
