// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's DistributedSlack outer loop (DistributedSlackLoop) in the NROuter_*
// algorithms: the slack bus' mismatch shared on the units flagged "can participate in the
// slack", bounded by their limits, and the same trigger in detection mode. C++14 only.

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "powerflow_algorithm/outer_loop/DistributedSlackLoop.hpp"
#include "powerflow_algorithm/outer_loop/ReactiveLimitsLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::DistributedSlackLoop;
using ls2g::ErrorType;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::ReactiveLimitsLoop;
using ls2g::real_type;

namespace {

// buses 0-1-2-3 in a line, a 200 MW load at bus 3; generators at buses 0 (the slack, 50 MW),
// 1 (40 MW) and 2 (60 MW), every one flagged "can participate in the slack" with the given
// weights, between 0 and its `max_p` MW
LSGrid make_grid(const std::vector<real_type> & weights,
                 const std::vector<real_type> & max_p = {1000., 1000., 1000.})
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
    grid.set_gen_can_participate_slack({true, true, true}, w);
    RealVect p_max(3);
    for (int i = 0; i < 3; ++i) p_max(i) = max_p[i];
    grid.set_gen_p_limits(RealVect::Zero(3), p_max);
    return grid;
}

CplxVect solve(LSGrid & grid)
{
    return grid.ac_pf(CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.)), 30, 1e-10);
}

void use_loop(LSGrid & grid, real_type threshold_mw = 1e-6, bool fail_on_residue = true)
{
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    DistributedSlackLoop::Params params;
    params.slack_bus_p_max_mismatch_mw = threshold_mw;
    params.fail_on_residue = fail_on_residue;
    grid.add_outer_loop(std::make_shared<DistributedSlackLoop>(params));
}

RealVect gen_p(const LSGrid & grid)
{
    const auto res_p = std::get<0>(grid.get_generators().get_res());
    RealVect res(res_p.size());
    for (int i = 0; i < res.size(); ++i) res(i) = static_cast<real_type>(res_p(i));
    return res;
}

}  // namespace

TEST_CASE("the loop shares the slack as the in-Newton distributed slack does", "[outer_loop][distributed_slack]")
{
    const std::vector<real_type> weights{1., 2., 3.};
    // the reference: lightsim2grid's in-Newton distributed slack, with the same weights
    LSGrid ref = make_grid(weights);
    ref.add_gen_slackbus(1, 2.);
    ref.add_gen_slackbus(2, 3.);
    ref.change_algorithm("NR_SparseLU");
    const CplxVect V_ref = solve(ref);
    REQUIRE(V_ref.size() == 4);

    LSGrid grid = make_grid(weights);
    use_loop(grid);
    const CplxVect V = solve(grid);
    REQUIRE(V.size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    CHECK(stats.nb_outer_iterations >= 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);

    // the loop stops once less than OpenLoadFlow's residue (1e-3 MW) is left
    for (int i = 0; i < 4; ++i) CHECK(std::abs(V(i) - V_ref(i)) < 1e-5);
    const RealVect p = gen_p(grid), p_ref = gen_p(ref);
    for (int i = 0; i < 3; ++i) CHECK(p(i) == Approx(p_ref(i)).margin(2e-3));
    // the shares follow the weights (compared between the two units that are not the slack
    // generator, which also publishes the leftover below the residue)
    CHECK((p(2) - 60.) == Approx(1.5 * (p(1) - 40.)).margin(1e-6));
}

TEST_CASE("a unit at its limit leaves the distribution", "[outer_loop][distributed_slack]")
{
    // generation 150 MW for a 200 MW load: more than 50 MW to share once the losses are in;
    // capped at 70 MW, the third unit (key 3 of 6) stops after taking 10 MW
    LSGrid grid = make_grid({1., 2., 3.}, {1000., 1000., 70.});
    use_loop(grid);
    REQUIRE(solve(grid).size() == 4);
    const RealVect p = gen_p(grid);
    CHECK(p(2) == Approx(70.));
    // the other two shared the rest, in the ratio of their keys (the slack generator also
    // publishes the leftover below OpenLoadFlow's residue, hence the margin)
    CHECK((p(1) - 40.) == Approx(2. * (p(0) - 50.)).margin(3e-3));
    CHECK((p(0) - 50.) + (p(1) - 40.) > 40.);
}

TEST_CASE("every unit at its limit: the solve fails, or the slack bus keeps the rest", "[outer_loop][distributed_slack]")
{
    SECTION("fail_on_residue (OpenLoadFlow's FAIL)") {
        LSGrid grid = make_grid({1., 2., 3.}, {65., 65., 65.});  // 195 MW at most for a 200 MW load
        use_loop(grid, 1e-6, true);
        CHECK(solve(grid).size() == 0);
        CHECK(grid.get_algo().get_error() == ErrorType::OuterLoopFailed);
        CHECK(grid.get_algo().get_outer_loop_stats().failed_loop == "DistributedSlack");
    }
    SECTION("leave it on the slack bus") {
        LSGrid grid = make_grid({1., 2., 3.}, {65., 65., 65.});
        use_loop(grid, 1e-6, false);
        REQUIRE(solve(grid).size() == 4);
        const RealVect p = gen_p(grid);
        CHECK(p(1) == Approx(65.));
        CHECK(p(2) == Approx(65.));
        CHECK(p(0) > 65.);  // the slack generator: its limit plus what is left
    }
}

TEST_CASE("below the threshold nothing moves", "[outer_loop][distributed_slack]")
{
    LSGrid grid = make_grid({1., 2., 3.});
    use_loop(grid, 1e6);  // a mismatch can never be that large here
    REQUIRE(solve(grid).size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().nb_outer_iterations == 0);
    const RealVect p = gen_p(grid);
    CHECK(p(1) == Approx(40.));
    CHECK(p(2) == Approx(60.));
}

TEST_CASE("detection mode reports the slack mismatch the loop would share", "[outer_loop][distributed_slack]")
{
    for (const char * algo : {"NRSing_SparseLU", "NR_SparseLU"}) {
        LSGrid grid = make_grid({1., 2., 3.});
        grid.change_algorithm(algo);
        // the default loop list holds a DistributedSlack loop (threshold 1 MW)
        REQUIRE(solve(grid).size() == 4);
        const std::vector<LimitViolation> out = grid.get_physical_violations(true, 1e-4, 1e-4);
        int nb = 0;
        for (const LimitViolation & v : out) {
            if (v.violation_type != LimitViolationType::SLACK_MISMATCH) continue;
            ++nb;
            // what the slack generator injected beyond its target
            CHECK(v.value == Approx(gen_p(grid)(0) - 50.).margin(1e-6));
            CHECK(v.limit == Approx(1.));
        }
        CHECK(nb == 1);
    }
}

TEST_CASE("without any unit flagged, the loop has nothing to do", "[outer_loop][distributed_slack]")
{
    // no distributed slack set up: OpenLoadFlow with distributedSlack off, the single slack kept
    LSGrid grid = make_grid({1., 2., 3.});
    grid.set_gen_can_participate_slack({false, false, false}, RealVect::Zero(3));
    grid.change_algorithm("NROuter_SparseLU");
    grid.reset_outer_loops();  // the default list, DistributedSlack in it
    REQUIRE(solve(grid).size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    // left out by is_needed (the other loops of the default list may run)
    for (const auto & it : stats.loop_iterations) CHECK(it.first != "DistributedSlack");
    CHECK(gen_p(grid)(1) == Approx(40.));
}

TEST_CASE("only the slack unit flagged: the loop leaves it past its limit, no failure", "[outer_loop][distributed_slack]")
{
    // the slack generator keeps its flag (it is in the slack) but has nothing to share with:
    // OpenLoadFlow with distributedSlack off, which converges with it beyond its max_p
    LSGrid grid = make_grid({1., 2., 3.}, {80., 1000., 1000.});
    grid.set_gen_can_participate_slack({false, false, false}, RealVect::Zero(3));
    REQUIRE(grid.get_generators().get_can_participate_slack(0));
    use_loop(grid, 1e-6, true);
    REQUIRE(solve(grid).size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    for (const auto & it : stats.loop_iterations) CHECK(it.first != "DistributedSlack");
    const RealVect p = gen_p(grid);
    CHECK(p(0) > 80.);
    CHECK(p(1) == Approx(40.));
    CHECK(p(2) == Approx(60.));
}

namespace {

// generator 1 has a capability curve, its fixed limits are the curve's at its target P
// (40 MW: +-30 MVar), and a voltage it cannot hold: the loops move it to a higher P, where
// its range is wider, and freeze it at that limit
LSGrid make_curve_grid()
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
    RealVect gen_p0(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    gen_p0 << 50., 40., 60.;
    gen_v << 1.02, 1.06, 1.01;
    gen_min_q << -500., -30., -500.;
    gen_max_q << 500., 30., 500.;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 2;
    grid.init_generators(gen_p0, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_can_participate_slack({true, true, true}, RealVect::Ones(3));
    grid.set_gen_p_limits(RealVect::Zero(3), RealVect::Constant(3, 1000.));
    Eigen::VectorXi curve_gen(2);
    curve_gen << 1, 1;
    RealVect curve_p(2), curve_min_q(2), curve_max_q(2);
    curve_p << 0., 100.;
    curve_min_q << -10., -60.;
    curve_max_q << 10., 60.;
    grid.set_gen_capability_curves(curve_gen, curve_p, curve_min_q, curve_max_q);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<DistributedSlackLoop>());
    grid.add_outer_loop(std::make_shared<ReactiveLimitsLoop>());
    return grid;
}

void check_no_reactive_violation(const std::vector<LimitViolation> & out)
{
    for (const LimitViolation & v : out) {
        INFO("violation " << static_cast<int>(v.violation_type) << " on " << v.element_id
             << ": " << v.value << " against " << v.limit);
        CHECK_FALSE((v.violation_type == LimitViolationType::HIGH_Q || v.violation_type == LimitViolationType::LOW_Q));
    }
}

}  // namespace

TEST_CASE("the reactive limits a solve's report reads are at the target P the loops gave", "[outer_loop][distributed_slack][reactive_limits]")
{
    LSGrid grid = make_curve_grid();
    REQUIRE(solve(grid).size() == 4);
    REQUIRE(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    const real_type p1 = gen_p(grid)(1);
    REQUIRE(p1 > 45.);
    const real_type q_max = 10. + 50. * p1 / 100.;
    const real_type q1 = std::get<1>(grid.get_generators().get_res())(1);
    REQUIRE(q1 == Approx(q_max).margin(1e-4));
    REQUIRE(q1 > 30. + 1.);

    SECTION("the grid's report") {
        // what the solve did is within the limits it read: nothing for that bus
        check_no_reactive_violation(grid.get_physical_violations(true, 1e-4, 1e-4));
    }
    SECTION("a batch row's") {
        ls2g::ContingencyAnalysis ca(grid);
        ca.set_compute_physical_violations(true);
        ca.add_nk(std::vector<int>());
        ca.compute(CplxVect::Constant(grid.total_bus(), ls2g::cplx_type(1., 0.)), 30, 1e-10);
        REQUIRE(ca.converged_mask()[0]);
        REQUIRE(ca.get_physical_violations().size() == 1);
        check_no_reactive_violation(ca.get_physical_violations()[0]);
        check_no_reactive_violation(ca.get_physical_violations_n());
    }
}
