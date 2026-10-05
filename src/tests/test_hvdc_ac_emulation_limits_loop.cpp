// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's AcHvdcAcEmulationLimits outer loop (HvdcAcEmulationLimitsLoop): an
// angle-droop hvdc line saturated at the limit of the direction its droop flow exceeds,
// released when that flow comes back inside, saturated on the other side when it reversed.
// C++14 only.

#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/HvdcAcEmulationLimitsLoop.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::HvdcAcEmulationLimitsLoop;
using ls2g::LSGrid;
using ls2g::OuterContext;
using ls2g::OuterLoopStatus;
using ls2g::OuterState;
using ls2g::RealVect;
using ls2g::real_type;

namespace {

const real_type P0_MW = 30.;
const real_type K_MW_PER_DEG = 400.;

// buses 0-1-2-3 in a line, an 80 MW load at bus 3, the slack generator at bus 0, and a
// lossless angle-droop hvdc line from bus 1 to bus 3 (p = P0 + K * (theta1 - theta2))
LSGrid make_grid(real_type pmax12, real_type pmax21, bool droop = true)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(4, 1, RealVect::Constant(4, 138.), 0, 0);
    Eigen::VectorXi from_id(3), to_id(3);
    from_id << 0, 1, 2;
    to_id << 1, 2, 3;
    grid.init_powerlines(RealVect::Constant(3, 0.01), RealVect::Constant(3, 0.1),
                         CplxVect::Zero(3), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 80.;
    load_q << 60.;
    Eigen::VectorXi load_bus(1);
    load_bus << 3;
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

    Eigen::VectorXi bus1(1), bus2(1);
    bus1 << 1;
    bus2 << 3;
    const std::vector<int> type{0}, mode{0};
    const std::vector<bool> vreg{false}, droop_on{droop};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm(1), q_lim(1), pf(1), p0(1), k(1), p12(1), p21(1);
    vm << 1.;
    q_lim << 1000.;
    pf << 1.;
    p0 << P0_MW;
    k << K_MW_PER_DEG;
    p12 << pmax12;
    p21 << pmax21;
    grid.init_hvdc_lines(bus1, bus2, type, type, zero, zero, vreg, vreg, vm, vm, zero, zero,
                         -q_lim, q_lim, -q_lim, q_lim, pf, pf, mode, zero, zero, zero, droop_on,
                         p0, k, p12, p21);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

// drives the loop by hand: the angles of the hvdc ends give the droop flow it sees
struct Driver
{
    LSGrid grid;
    HvdcAcEmulationLimitsLoop loop;
    CplxVect Sbus;
    OuterState state;
    ls2g::OuterControls controls;
    RealVect Va;
    int b1, b2;

    explicit Driver(real_type pmax12, real_type pmax21) : grid(make_grid(pmax12, pmax21))
    {
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 4);  // builds the bus maps
        b1 = grid.id_me_to_ac_solver()[1].cast_int();
        b2 = grid.id_me_to_ac_solver()[3].cast_int();
        Sbus = CplxVect::Zero(4);
        state.Sbus = &Sbus;
        state.Sbus_init = &Sbus;
        Va = RealVect::Zero(4);
        OuterContext ctx = context();
        ls2g::OuterDeclaration decl;
        loop.declare(ctx, decl);  // reserves the line's regime
        REQUIRE(loop.is_needed(ctx));
        loop.initialize(ctx);
    }

    OuterContext context()
    {
        OuterContext ctx;
        ctx.grid = &grid;
        ctx.Va = &Va;
        ctx.state = &state;
        ctx.controls = &controls;
        return ctx;
    }

    // the regime the loop set for the line
    int regime() const { return controls.hvdc_regime(0)->regime(); }

    // the loop's check with a droop flow of `p_mw` (side 1 -> side 2)
    OuterLoopStatus check(real_type p_mw)
    {
        Va.setZero();
        Va(b1) = (p_mw - P0_MW) / K_MW_PER_DEG * M_PI / 180.;
        OuterContext ctx = context();
        return loop.check(ctx);
    }
};

}  // namespace

TEST_CASE("the regime follows OpenLoadFlow's rules", "[outer][hvdc]")
{
    Driver d(/*pmax12=*/50., /*pmax21=*/40.);
    CHECK(d.check(45.) == OuterLoopStatus::STABLE);  // inside: nothing to do
    CHECK(d.regime() == ls2g::HvdcRegimeControl::KEEP);

    CHECK(d.check(70.) == OuterLoopStatus::UNSTABLE);  // above pmax12: saturated 1 -> 2
    CHECK(d.regime() == 1);
    CHECK(d.check(70.) == OuterLoopStatus::STABLE);    // still beyond it: kept

    CHECK(d.check(45.) == OuterLoopStatus::UNSTABLE);  // back strictly inside: released
    CHECK(d.regime() == 0);

    CHECK(d.check(70.) == OuterLoopStatus::UNSTABLE);
    CHECK(d.check(-60.) == OuterLoopStatus::UNSTABLE);  // reversed beyond pmax21
    CHECK(d.regime() == -1);
    CHECK(d.check(-60.) == OuterLoopStatus::STABLE);
    CHECK(d.check(-30.) == OuterLoopStatus::UNSTABLE);  // back inside the reverse limit
    CHECK(d.regime() == 0);

    // the grid's own regime is never written
    CHECK(d.grid.get_status_droop_hvdc(0) == 0);
}

TEST_CASE("NROuter saturates the line, with one analysis", "[outer][hvdc]")
{
    LSGrid grid = make_grid(/*pmax12=*/5., /*pmax21=*/500.);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<HvdcAcEmulationLimitsLoop>());
    const CplxVect V = grid.ac_pf(flat_start(grid), 30, 1e-10);
    REQUIRE(V.size() == 4);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    REQUIRE(stats.loop_iterations.size() == 1);
    CHECK(stats.loop_iterations[0].second == 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);

    // published as saturated (lossless line), the line's own regime untouched
    CHECK(std::get<0>(grid.get_dclines().get_res_side_1())(0) == Approx(-5.));
    CHECK(std::get<0>(grid.get_dclines().get_res_side_2())(0) == Approx(5.));
    CHECK(grid.get_status_droop_hvdc(0) == 0);

    // the same solve as a single-slack Newton with the line saturated by hand
    LSGrid ref = make_grid(5., 500.);
    ref.change_algorithm(AlgorithmType::NRSing_SparseLU);
    ref.set_status_droop_hvdc(0, 1);
    const CplxVect V_ref = ref.ac_pf(flat_start(ref), 30, 1e-10);
    REQUIRE(V_ref.size() == 4);
    CHECK((V - V_ref).cwiseAbs().maxCoeff() < 1e-9);

    // a next solve starts again from the linear regime
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 4);
    CHECK(grid.get_algo().get_outer_loop_stats().loop_iterations[0].second == 1);
}

TEST_CASE("not needed without a line in AC emulation", "[outer][hvdc]")
{
    LSGrid off = make_grid(5., 500., /*droop=*/false);
    off.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(off.ac_pf(flat_start(off), 30, 1e-10).size() == 4);
    OuterContext ctx;
    ctx.grid = &off;
    HvdcAcEmulationLimitsLoop loop;
    CHECK_FALSE(loop.is_needed(ctx));

    // a line the caller froze at a limit is not in AC emulation either
    LSGrid frozen = make_grid(5., 500.);
    frozen.set_status_droop_hvdc(0, 1);
    frozen.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(frozen.ac_pf(flat_start(frozen), 30, 1e-10).size() == 4);
    ctx.grid = &frozen;
    CHECK_FALSE(loop.is_needed(ctx));
}
