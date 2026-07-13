// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Basic tests of the LSGrid main API, driven entirely from C++ (no python,
// no grid2op, no external grid file): a 3-bus radial grid is built through
// the init_* methods and solved with the always-available Eigen SparseLU
// algorithms (the LSGrid constructor defaults, no KLU needed). Covers the
// AC / DC powerflow contract (converged => per-bus V, diverged => empty
// vector), physically-checkable results (power balance, DC angles), copy /
// get_state / save_binary round trips, setpoint changes, element
// deactivation, and the documented error paths.

#include <cmath>
#include <complex>
#include <stdexcept>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "test_helpers.hpp"

using Catch::Approx;
using ls2g::LSGrid;
using ls2g::CplxVect;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// Slack generator (1.02 pu) at bus 0, two identical lines 0-1 and 1-2
// (x = 0.1 pu each), one 50 MW / 10 MVar load at bus 2. sn_mva = 100, so the
// load is 0.5 pu and (with r = 0) the exact DC angles are 0 / -0.05 / -0.10.
LSGrid make_three_bus_grid(real_type r_pu = 0.01)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    RealVect bus_vn_kv(3);
    bus_vn_kv << 138., 138., 138.;
    grid.init_bus(3, 1, bus_vn_kv, 0, 0);  // nb_line / nb_trafo args are unused

    RealVect branch_r(2), branch_x(2);
    branch_r << r_pu, r_pu;
    branch_x << 0.1, 0.1;
    const CplxVect branch_h = CplxVect::Zero(2);
    Eigen::VectorXi from_id(2), to_id(2);
    from_id << 0, 1;
    to_id << 1, 2;
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 10.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
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
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

CplxVect solve_ac(LSGrid & grid)
{
    return grid.ac_pf(flat_start(grid), 20, 1e-8);
}

}  // anonymous namespace

TEST_CASE("AC powerflow on a hand-built 3-bus grid converges and balances", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    const CplxVect V = solve_ac(grid);

    // converged <=> non-empty per-bus voltage vector
    REQUIRE(V.size() == 3);
    CHECK(grid.get_algo().converged());
    CHECK(grid.get_algo().get_error() == ls2g::ErrorType::NoError);
    CHECK(grid.get_algo().get_nb_iter() >= 1);

    // slack holds its voltage setpoint; magnitude drops along the feeder
    CHECK(std::abs(V(0)) == Approx(1.02).epsilon(1e-8));
    CHECK(std::abs(V(1)) < std::abs(V(0)));
    CHECK(std::abs(V(2)) < std::abs(V(1)));

    // the load draws exactly its setpoint
    const auto load_res = grid.get_loads_res();  // (p_mw, q_mvar, v_kv)
    CHECK(std::get<0>(load_res)(0) == Approx(50.).epsilon(1e-8));
    CHECK(std::get<1>(load_res)(0) == Approx(10.).epsilon(1e-8));

    // power balance: the slack covers the load plus the (positive) line
    // losses, and its injection is what enters line 0 at bus 0
    const auto line1 = grid.get_line_res1();  // "from" side (p_mw, q_mvar, v_kv, a)
    const auto line2 = grid.get_line_res2();  // "to" side
    const real_type loss0 = std::get<0>(line1)(0) + std::get<0>(line2)(0);
    const real_type loss1 = std::get<0>(line1)(1) + std::get<0>(line2)(1);
    CHECK(loss0 > 0.);
    CHECK(loss1 > 0.);
    const auto gen_res = grid.get_gen_res();
    const real_type slack_p = std::get<0>(gen_res)(0);
    CHECK(slack_p == Approx(50. + loss0 + loss1).epsilon(1e-6));
    CHECK(slack_p == Approx(std::get<0>(line1)(0)).epsilon(1e-6));
}

TEST_CASE("DC powerflow reproduces the analytic angles of a lossless feeder", "[LSGrid]")
{
    // r = 0 so theta_k = -P * x * k holds exactly: 0.5 pu through x = 0.1 pu
    // per line gives 0 / -0.05 / -0.10 rad
    LSGrid grid = make_three_bus_grid(0.);
    const CplxVect V = grid.dc_pf(flat_start(grid), 1, 1e-8);

    REQUIRE(V.size() == 3);
    CHECK(std::arg(V(0)) == Approx(0.).margin(1e-10));
    CHECK(std::arg(V(1)) == Approx(-0.05).epsilon(1e-8));
    CHECK(std::arg(V(2)) == Approx(-0.10).epsilon(1e-8));

    // both lines carry the full 50 MW, losslessly
    const auto line1 = grid.get_line_res1();
    CHECK(std::get<0>(line1)(0) == Approx(50.).epsilon(1e-8));
    CHECK(std::get<0>(line1)(1) == Approx(50.).epsilon(1e-8));
}

TEST_CASE("voltage getters follow the powerflow", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();

    SECTION("before any powerflow they throw") {
        CHECK_THROWS_AS(grid.get_V(), std::runtime_error);
        CHECK_THROWS_AS(grid.get_Va(), std::runtime_error);
        CHECK_THROWS_AS(grid.get_Vm(), std::runtime_error);
    }

    SECTION("after a powerflow they match the returned vector") {
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 3);
        const CplxVect V_get = grid.get_V();
        const RealVect Vm = grid.get_Vm();
        const RealVect Va = grid.get_Va();
        REQUIRE(V_get.size() == 3);
        for (Eigen::Index i = 0; i < 3; ++i) {
            CHECK(std::abs(V_get(i) - V(i)) < 1e-10);
            CHECK(Vm(i) == Approx(std::abs(V(i))).epsilon(1e-10));
            CHECK(Va(i) == Approx(std::arg(V(i))).margin(1e-10));
        }
    }
}

TEST_CASE("a diverging AC powerflow returns an empty vector, not garbage", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    // 100 GW through a 0.1 pu line cannot converge
    grid.change_p_load(0, 1e5);
    const CplxVect V = solve_ac(grid);
    CHECK(V.size() == 0);
    CHECK_FALSE(grid.get_algo().converged());
    CHECK(grid.get_algo().get_error() != ls2g::ErrorType::NoError);
}

TEST_CASE("changing setpoints takes effect on the next solve", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    REQUIRE(solve_ac(grid).size() == 3);
    const real_type slack_p_50 = std::get<0>(grid.get_gen_res())(0);

    grid.change_p_load(0, 80.);
    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::get<0>(grid.get_loads_res())(0) == Approx(80.).epsilon(1e-8));
    CHECK(std::get<0>(grid.get_gen_res())(0) > slack_p_50);

    grid.change_v_gen(0, 1.05);
    const CplxVect V = solve_ac(grid);
    REQUIRE(V.size() == 3);
    CHECK(std::abs(V(0)) == Approx(1.05).epsilon(1e-8));
}

TEST_CASE("deactivating and reactivating a load", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();

    grid.deactivate_load(0);
    REQUIRE(solve_ac(grid).size() == 3);
    // nothing consumed: the slack only covers (near-zero) losses
    CHECK(std::get<0>(grid.get_loads_res())(0) == Approx(0.).margin(1e-8));
    CHECK(std::get<0>(grid.get_gen_res())(0) == Approx(0.).margin(1e-6));

    grid.reactivate_load(0);
    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::get<0>(grid.get_loads_res())(0) == Approx(50.).epsilon(1e-8));
}

TEST_CASE("copies are independent and solve to the same state", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    LSGrid other = grid.copy();

    const CplxVect V = solve_ac(grid);
    const CplxVect V_other = solve_ac(other);
    REQUIRE(V.size() == 3);
    REQUIRE(V_other.size() == 3);
    CHECK((V - V_other).norm() < 1e-10);

    // mutating the copy must not touch the original
    other.change_p_load(0, 80.);
    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::get<0>(grid.get_loads_res())(0) == Approx(50.).epsilon(1e-8));
}

TEST_CASE("get_state / set_state round-trips the whole grid", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    LSGrid::StateRes state = grid.get_state();

    LSGrid restored;
    restored.set_state(state);

    const CplxVect V = solve_ac(grid);
    const CplxVect V_restored = solve_ac(restored);
    REQUIRE(V.size() == 3);
    REQUIRE(V_restored.size() == 3);
    CHECK((V - V_restored).norm() < 1e-10);
    CHECK(restored.get_sn_mva() == Approx(100.));
}

TEST_CASE("save_binary / load_binary round-trips the whole grid", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();
    ls2g_test::TempFile file;
    grid.save_binary(file.str());

    LSGrid loaded = LSGrid::load_binary(file.str());
    const CplxVect V = solve_ac(grid);
    const CplxVect V_loaded = solve_ac(loaded);
    REQUIRE(V.size() == 3);
    REQUIRE(V_loaded.size() == 3);
    CHECK((V - V_loaded).norm() < 1e-10);
}

TEST_CASE("documented error paths", "[LSGrid]")
{
    LSGrid grid = make_three_bus_grid();

    SECTION("Vinit of the wrong size is rejected") {
        CHECK_THROWS_AS(grid.ac_pf(CplxVect::Constant(2, cplx_type(1., 0.)), 20, 1e-8),
                        std::runtime_error);
        CHECK_THROWS_AS(grid.dc_pf(CplxVect::Constant(7, cplx_type(1., 0.)), 1, 1e-8),
                        std::runtime_error);
    }
    SECTION("slack registration validates the generator id and weight") {
        CHECK_THROWS_AS(grid.add_gen_slackbus(12, 1.), std::runtime_error);
        CHECK_THROWS_AS(grid.add_gen_slackbus(-1, 1.), std::runtime_error);
        CHECK_THROWS_AS(grid.add_gen_slackbus(0, 0.), std::runtime_error);
    }
    SECTION("element setters validate the element id") {
        CHECK_THROWS_AS(grid.change_p_load(3, 10.), std::out_of_range);
        CHECK_THROWS_AS(grid.change_v_gen(3, 1.), std::out_of_range);
        CHECK_THROWS_AS(grid.deactivate_load(3), std::out_of_range);
    }
}
