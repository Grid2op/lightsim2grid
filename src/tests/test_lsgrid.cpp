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
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "case_exotic_elements.hpp"
#include "test_helpers.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::LSGrid;
using ls2g::CplxVect;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// Three 138 kV buses, two identical lines 0-1 and 1-2 (x = 0.1 pu each), one
// 50 MW / 10 MVar load at bus 2. sn_mva = 100, so the load is 0.5 pu and
// (with r = 0) the exact DC angles are 0 / -0.05 / -0.10. No generator yet.
LSGrid make_three_bus_skeleton(real_type r_pu)
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
    return grid;
}

// the skeleton plus a single slack generator (1.02 pu) at bus 0
LSGrid make_three_bus_grid(real_type r_pu = 0.01)
{
    LSGrid grid = make_three_bus_skeleton(r_pu);
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

// the skeleton plus generators at bus 0 (1.02 pu) and bus 1 (gen1_v_pu).
// NO slack is registered: each test declares the slack setup it needs.
LSGrid make_three_bus_grid_two_gens(real_type gen1_v_pu = 1.02)
{
    LSGrid grid = make_three_bus_skeleton(0.01);
    RealVect gen_p(2), gen_v(2), gen_min_q(2), gen_max_q(2);
    gen_p << 0., 0.;
    gen_v << 1.02, gen1_v_pu;
    gen_min_q << -1000., -1000.;
    gen_max_q << 1000., 1000.;
    Eigen::VectorXi gen_bus(2);
    gen_bus << 0, 1;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
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
    // the exotic-elements case (IEEE14 + SVC + storage + 3 HVDC lines + a
    // phase-shifting transformer, see src/tests/case_exotic_elements.hpp) rather
    // than the 3-bus skeleton, so this round trip actually exercises every
    // element type LSGrid knows how to serialize, not just lines/loads/a slack gen.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    ls2g_test::TempFile file;
    grid.save_binary(file.str());

    LSGrid loaded = LSGrid::load_binary(file.str());
    const CplxVect V = solve_ac(grid);
    const CplxVect V_loaded = solve_ac(loaded);
    REQUIRE(V.size() == 14);
    REQUIRE(V_loaded.size() == 14);
    CHECK((V - V_loaded).norm() < 1e-10);
}

TEST_CASE("a capacitive shunt raises the local voltage", "[LSGrid][shunt]")
{
    LSGrid base = make_three_bus_grid();
    REQUIRE(solve_ac(base).size() == 3);
    const real_type vm2_without = std::abs(base.get_V()(2));

    LSGrid grid = make_three_bus_grid();
    // pandapower convention: q_mvar is what the shunt consumes at v = 1 pu,
    // so -25 MVar is a 25 MVar capacitor bank
    RealVect shunt_p(1), shunt_q(1);
    shunt_p << 0.;
    shunt_q << -25.;
    Eigen::VectorXi shunt_bus(1);
    shunt_bus << 2;
    grid.init_shunt(shunt_p, shunt_q, shunt_bus);

    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::abs(grid.get_V()(2)) > vm2_without);
    // the shunt injects reactive power (generator sign convention on results)
    CHECK(std::get<1>(grid.get_shunts_res())(0) < 0.);
}

TEST_CASE("a charging storage unit behaves as an extra load", "[LSGrid][storage]")
{
    LSGrid grid = make_three_bus_grid();
    // load convention: positive target_p = charging (consumes from the grid)
    RealVect storage_p(1), storage_q(1);
    storage_p << 10.;
    storage_q << 0.;
    Eigen::VectorXi storage_bus(1);
    storage_bus << 1;
    grid.init_storages(storage_p, storage_q, storage_bus);

    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::get<0>(grid.get_storages_res())(0) == Approx(10.).epsilon(1e-8));
    // the slack now covers the load, the storage and the (positive) losses
    CHECK(std::get<0>(grid.get_gen_res())(0) > 60.);
}

TEST_CASE("static var compensators", "[LSGrid][svc]")
{
    LSGrid base = make_three_bus_grid();
    REQUIRE(solve_ac(base).size() == 3);
    const CplxVect V_without = base.get_V();

    // SvcContainer::RegulationMode: OFF = 0, VOLTAGE = 1, REACTIVE_POWER = 2
    RealVect slope(1), b_min(1), b_max(1);
    slope << 0.;
    b_min << -100.;
    b_max << 100.;
    Eigen::VectorXi reg_bus(1), svc_bus(1);
    reg_bus << 2;
    svc_bus << 2;

    SECTION("VOLTAGE mode holds the regulated bus at its setpoint") {
        LSGrid grid = make_three_bus_grid();
        RealVect target_vm(1), q_set(1);
        target_vm << 1.03;
        q_set << 0.;
        grid.init_svcs({1}, target_vm, q_set, slope, b_min, b_max, reg_bus, svc_bus);
        grid.tell_solver_need_reset();
        REQUIRE(solve_ac(grid).size() == 3);
        CHECK(std::abs(grid.get_V()(2)) == Approx(1.03).epsilon(1e-8));
        // raising the bus above its natural voltage takes reactive injection
        CHECK(std::get<1>(grid.get_svcs().get_res())(0) > 0.);
    }
    SECTION("REACTIVE_POWER mode injects exactly its setpoint") {
        LSGrid grid = make_three_bus_grid();
        RealVect target_vm(1), q_set(1);
        target_vm << 1.0;
        q_set << 25.;
        grid.init_svcs({2}, target_vm, q_set, slope, b_min, b_max, reg_bus, svc_bus);
        grid.tell_solver_need_reset();
        REQUIRE(solve_ac(grid).size() == 3);
        CHECK(std::get<1>(grid.get_svcs().get_res())(0) == Approx(25.).epsilon(1e-8));
        CHECK(std::abs(grid.get_V()(2)) > std::abs(V_without(2)));
    }
    SECTION("OFF mode changes nothing") {
        LSGrid grid = make_three_bus_grid();
        RealVect target_vm(1), q_set(1);
        target_vm << 1.03;
        q_set << 25.;
        grid.init_svcs({0}, target_vm, q_set, slope, b_min, b_max, reg_bus, svc_bus);
        grid.tell_solver_need_reset();
        REQUIRE(solve_ac(grid).size() == 3);
        CHECK((grid.get_V() - V_without).norm() < 1e-10);
    }
}

TEST_CASE("HVDC links between two converter stations", "[LSGrid][hvdc]")
{
    // one hvdc line between bus 1 (side 1, rectifier) and bus 2 (side 2);
    // helper filling the many init_hvdc_lines arguments
    const auto add_hvdc = [](LSGrid & grid, int type2, bool vreg2_on,
                             real_type vm2, real_type power_factor2,
                             real_type p_setpoint = 20.,
                             bool droop = false, real_type droop_p0 = 0.,
                             real_type droop_k_per_deg = 0.) {
        Eigen::VectorXi bus1(1), bus2(1);
        bus1 << 1;
        bus2 << 2;
        // ConverterStationContainer::ConverterType: VSC = 0, LCC = 1
        const std::vector<int> type1{0};
        // HvdcLineContainer::ConvertersMode: SIDE_1_RECTIFIER = 0
        const std::vector<int> mode{0};
        const std::vector<bool> vreg1{false}, vreg2{vreg2_on}, droop_on{droop};
        const RealVect zero = RealVect::Zero(1);
        RealVect vm1_pu(1), vm2_pu(1), q_min(1), q_max(1);
        RealVect pf1(1), pf2(1), p_set(1), p_max(1), p0(1), k_deg(1);
        vm1_pu << 1.0;
        vm2_pu << vm2;
        q_min << -1e6;
        q_max << 1e6;
        pf1 << 1.;
        pf2 << power_factor2;
        p_set << p_setpoint;
        p_max << 300.;
        p0 << droop_p0;
        k_deg << droop_k_per_deg;
        grid.init_hvdc_lines(bus1, bus2, type1, {type2},
                             zero, zero,           // no station losses
                             vreg1, vreg2, vm1_pu, vm2_pu,
                             zero, zero,           // q setpoints
                             q_min, q_max, q_min, q_max,
                             pf1, pf2, mode, p_set,
                             zero, zero,           // r_ohm, nominal_v (no line loss)
                             droop_on, p0, k_deg, p_max, p_max);
        grid.tell_solver_need_reset();
    };

    SECTION("a lossless VSC-VSC link moves its power setpoint") {
        LSGrid grid = make_three_bus_grid();
        add_hvdc(grid, /*type2=VSC*/ 0, false, 1.0, 1.);
        REQUIRE(solve_ac(grid).size() == 3);
        // generator sign convention: the rectifier consumes at bus 1, the
        // inverter injects the (lossless) 20 MW at bus 2
        CHECK(std::get<0>(grid.get_dcline_res1())(0) == Approx(-20.).epsilon(1e-8));
        CHECK(std::get<0>(grid.get_dcline_res2())(0) == Approx(20.).epsilon(1e-8));
    }
    SECTION("a voltage-regulating VSC station holds its bus") {
        LSGrid grid = make_three_bus_grid();
        add_hvdc(grid, /*type2=VSC*/ 0, true, 1.03, 1.);
        REQUIRE(solve_ac(grid).size() == 3);
        CHECK(std::abs(grid.get_V()(2)) == Approx(1.03).epsilon(1e-8));
    }
    SECTION("an LCC station consumes reactive per its power factor") {
        LSGrid grid = make_three_bus_grid();
        add_hvdc(grid, /*type2=LCC*/ 1, false, 1.0, 0.9);
        REQUIRE(solve_ac(grid).size() == 3);
        const real_type p2 = std::get<0>(grid.get_dcline_res2())(0);
        const real_type q2 = std::get<1>(grid.get_dcline_res2())(0);
        CHECK(p2 == Approx(20.).epsilon(1e-8));
        CHECK(q2 == Approx(-std::abs(p2) * std::tan(std::acos(0.9))).epsilon(1e-6));
    }
    SECTION("a VSC-VSC link with angle droop follows the droop law") {
        LSGrid grid = make_three_bus_grid();
        const real_type p0 = 10., k_per_deg = 2.;
        add_hvdc(grid, /*type2=VSC*/ 0, false, 1.0, 1., /*p_setpoint=*/0.,
                 /*droop=*/true, p0, k_per_deg);
        const CplxVect V = solve_ac(grid);
        REQUIRE(V.size() == 3);
        // linear regime: transferred power (1 -> 2) is p0 + k * (theta1 - theta2),
        // with k converted from MW/deg to MW/rad
        const real_type k_per_rad = k_per_deg * 180. / std::acos(real_type(-1.));
        const real_type expected = p0 + k_per_rad * (std::arg(V(1)) - std::arg(V(2)));
        CHECK(std::get<0>(grid.get_dcline_res2())(0) == Approx(expected).epsilon(1e-6));
        CHECK(std::get<0>(grid.get_dcline_res1())(0) == Approx(-expected).epsilon(1e-6));
    }
    SECTION("angle droop cannot be combined with an LCC station") {
        LSGrid grid = make_three_bus_grid();
        REQUIRE_THROWS_MATCHES(
            add_hvdc(grid, /*type2=LCC*/ 1, false, 1.0, 0.9, /*p_setpoint=*/0.,
                     /*droop=*/true, 10., 2.),
            std::runtime_error,
            MessageMatches(ContainsSubstring("only available for VSC - VSC")));
    }
}

TEST_CASE("distributed slack splits the mismatch by weight", "[LSGrid][slack]")
{
    LSGrid grid = make_three_bus_grid_two_gens();
    grid.add_gen_slackbus(0, 0.7);
    grid.add_gen_slackbus(1, 0.3);

    REQUIRE(solve_ac(grid).size() == 3);
    // both setpoints are 0 MW, so each gen's whole output is its slack share,
    // and shares are proportional to the weights
    const auto gen_res = grid.get_gen_res();
    const real_type p0 = std::get<0>(gen_res)(0);
    const real_type p1 = std::get<0>(gen_res)(1);
    CHECK(p0 > 0.);
    CHECK(p1 > 0.);
    CHECK(p0 / 0.7 == Approx(p1 / 0.3).epsilon(1e-6));
}

TEST_CASE("remote voltage control pins the regulated bus", "[LSGrid][vctrl]")
{
    // the PV generator at bus 1 (target 1.04 pu) regulates bus 2 instead of
    // its own bus (bordered VoltageControl extension, AC NR only)
    LSGrid grid = make_three_bus_grid_two_gens(1.04);
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 2);
    grid.tell_solver_need_reset();

    REQUIRE(solve_ac(grid).size() == 3);
    CHECK(std::abs(grid.get_V()(2)) == Approx(1.04).epsilon(1e-8));

    // out-of-range regulated bus is rejected
    CHECK_THROWS_AS(grid.set_gen_regulated_bus(1, 12), std::out_of_range);
}

TEST_CASE("a local and a remote voltage controller on one bus is unfeasible", "[LSGrid][vctrl]")
{
    // two generators on bus 1: gen 1 regulates its own bus (ordinary PV
    // pinning) while gen 2 regulates bus 2 remotely. Gen 2's own bus then
    // has no reactive equation left to trade against the remote target, so
    // the solve must reject the state instead of producing garbage.
    LSGrid grid = make_three_bus_skeleton(0.01);
    RealVect gen_p(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    gen_p << 0., 0., 0.;
    gen_v << 1.02, 1.02, 1.04;
    gen_min_q << -1000., -1000., -1000.;
    gen_max_q << 1000., 1000., 1000.;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 1;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(2, 2);
    grid.tell_solver_need_reset();

    REQUIRE_THROWS_MATCHES(
        solve_ac(grid), std::runtime_error,
        MessageMatches(ContainsSubstring("no reactive (Q) equation")));
}

TEST_CASE("transformers: tap ratio and phase shifter", "[LSGrid][trafo]")
{
    // buses at 138 kV, slack gen (1.02 pu) at bus 0, 50 MW / 10 MVar load at
    // bus 2. `loop` false: line 0-1 plus a trafo 1-2 (radial, for the ratio
    // effect); true: lines 0-1 and 1-2 plus a trafo 0-2 closing a loop (a
    // phase shifter only redirects flow inside a loop).
    const auto make_trafo_grid = [](bool loop, real_type ratio, real_type shift_deg) {
        LSGrid grid;
        grid.set_sn_mva(100.);
        grid.set_init_vm_pu(1.0);
        RealVect bus_vn_kv(3);
        bus_vn_kv << 138., 138., 138.;
        grid.init_bus(3, 1, bus_vn_kv, 0, 0);

        const int nb_line = loop ? 2 : 1;
        RealVect line_r(nb_line), line_x(nb_line);
        Eigen::VectorXi line_from(nb_line), line_to(nb_line);
        if (loop) {
            line_r << 0.01, 0.01;
            line_x << 0.1, 0.1;
            line_from << 0, 1;
            line_to << 1, 2;
        } else {
            line_r << 0.01;
            line_x << 0.1;
            line_from << 0;
            line_to << 1;
        }
        grid.init_powerlines(line_r, line_x, CplxVect::Zero(nb_line), line_from, line_to);

        RealVect trafo_r(1), trafo_x(1), trafo_ratio(1), trafo_shift(1);
        trafo_r << 0.01;
        trafo_x << 0.2;
        trafo_ratio << ratio;
        trafo_shift << shift_deg;  // degrees at init time
        Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
        trafo_bus1 << (loop ? 0 : 1);
        trafo_bus2 << 2;
        grid.init_trafo(trafo_r, trafo_x, CplxVect::Zero(1), trafo_ratio, trafo_shift,
                        {true}, trafo_bus1, trafo_bus2, true);

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
    };

    SECTION("the tap ratio scales the downstream voltage") {
        LSGrid grid = make_trafo_grid(false, 1., 0.);
        REQUIRE(solve_ac(grid).size() == 3);
        const real_type vm2_ratio1 = std::abs(grid.get_V()(2));

        grid.change_ratio_trafo(0, 1.05);
        const CplxVect V_changed = solve_ac(grid);
        REQUIRE(V_changed.size() == 3);
        // raising the tap on the hv side lowers the lv-side voltage, by
        // roughly the ratio (not exactly: the load's voltage dependence)
        CHECK(std::abs(V_changed(2)) < vm2_ratio1);
        CHECK(std::abs(V_changed(2)) == Approx(vm2_ratio1 / 1.05).epsilon(0.02));

        // changing the ratio on a solved grid == initializing with it
        LSGrid fresh = make_trafo_grid(false, 1.05, 0.);
        const CplxVect V_fresh = solve_ac(fresh);
        REQUIRE(V_fresh.size() == 3);
        CHECK((V_changed - V_fresh).norm() < 1e-8);
    }

    SECTION("the phase shifter redirects the loop flow") {
        LSGrid grid = make_trafo_grid(true, 1., 0.);
        REQUIRE(solve_ac(grid).size() == 3);
        const real_type p_trafo_0 = std::get<0>(grid.get_trafo_res1())(0);

        grid.change_shift_trafo_deg(0, 5.);
        REQUIRE(solve_ac(grid).size() == 3);
        const real_type p_trafo_plus = std::get<0>(grid.get_trafo_res1())(0);
        const CplxVect V_plus = grid.get_V();

        grid.change_shift_trafo_deg(0, -5.);
        REQUIRE(solve_ac(grid).size() == 3);
        const real_type p_trafo_minus = std::get<0>(grid.get_trafo_res1())(0);

        // opposite shifts move the trafo flow in opposite directions, and
        // the load is served throughout
        CHECK(p_trafo_plus != Approx(p_trafo_0).margin(1.));
        CHECK((p_trafo_plus - p_trafo_0) * (p_trafo_minus - p_trafo_0) < 0.);
        CHECK(std::get<0>(grid.get_loads_res())(0) == Approx(50.).epsilon(1e-8));

        // changing the shift on a solved grid == initializing with it
        LSGrid fresh = make_trafo_grid(true, 1., 5.);
        const CplxVect V_fresh = solve_ac(fresh);
        REQUIRE(V_fresh.size() == 3);
        CHECK((V_plus - V_fresh).norm() < 1e-8);
    }
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
