// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Exercises `case_exotic_elements.hpp` (an AUTO-GENERATED, dependency-free C++
// reproduction of lightsim2grid.tests._exotic_elements_fixture.build_exotic_elements_grid():
// IEEE14 + SVC + storage + 3 HVDC lines + a phase-shifting transformer) purely from C++,
// no pypowsybl / Python involved. The expected bus voltages below were captured once from
// the Python side (lightsim2grid.tests._exotic_elements_case_ls.build_exotic_elements_case_grid,
// which is bit-identical to the pypowsybl-built reference -- see
// lightsim2grid/tests/test_exotic_elements_case_ls.py) with a flat start and the same
// (20, 1e-7) ac_pf tolerance used here.

#include <cmath>
#include <complex>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "case_exotic_elements.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

TEST_CASE("exotic elements case grid: AC powerflow converges", "[LSGrid][exotic_elements]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    const CplxVect v0 = CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
    const CplxVect V = grid.ac_pf(v0, 20, 1e-7);

    REQUIRE(V.size() == 14);
    CHECK(grid.get_algo().converged());
    CHECK(grid.get_algo().get_error() == ls2g::ErrorType::NoError);

    // power balance sanity check: every load is served (nothing lost to a botched
    // conversion), and the slack generator's output is strictly positive
    const real_type total_load_p = std::get<0>(grid.get_loads_res()).sum();
    CHECK(total_load_p > 0.);
    const real_type slack_p = std::get<0>(grid.get_gen_res())(0);
    CHECK(slack_p > 0.);
}

TEST_CASE("exotic elements case grid: bus voltages match the Python (pypowsybl) reference",
          "[LSGrid][exotic_elements]")
{
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    const CplxVect v0 = CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
    const CplxVect V = grid.ac_pf(v0, 20, 1e-7);
    REQUIRE(V.size() == 14);

    // captured from lightsim2grid.tests._exotic_elements_case_ls.build_exotic_elements_case_grid()
    // (bit-identical to the pypowsybl-built reference, see test_exotic_elements_case_ls.py)
    CplxVect expected(14);
    expected << std::complex<double>(1.0600000000000001, 0),
                std::complex<double>(0.96799351594580019, -0.25097520412759489),
                std::complex<double>(0.97102700718869095, -0.23896977070368952),
                std::complex<double>(0.97196619148095098, -0.23512065544740082),
                std::complex<double>(0.99485520142802564, -0.26337707742414879),
                std::complex<double>(0.9594907122432127, -0.27213173889703429),
                std::complex<double>(1.041138716753307, -0.089750612684686787),
                std::complex<double>(0.98545996908740219, -0.22128861092757629),
                std::complex<double>(0.99403667576525778, -0.17572183506772299),
                std::complex<double>(1.0027959169054035, -0.15194114497901523),
                std::complex<double>(1.0363593510923492, -0.26619409347963563),
                std::complex<double>(1.0077294791673204, -0.22951146148151697),
                std::complex<double>(1.0627848391631813, -0.24205037832007373),
                std::complex<double>(0.96878618486527368, -0.24789781768742494);

    // tolerance matches ac_pf's own convergence tolerance (1e-7)
    for (Eigen::Index i = 0; i < expected.size(); ++i) {
        CHECK(V(i).real() == Approx(expected(i).real()).margin(1e-7));
        CHECK(V(i).imag() == Approx(expected(i).imag()).margin(1e-7));
    }
}

TEST_CASE("exotic elements case grid: per-element SVC / HVDC / storage results match the Python reference",
          "[LSGrid][exotic_elements]")
{
    // captured from lightsim2grid.tests._exotic_elements_case_ls.build_exotic_elements_case_grid()
    // (grid.get_svcs() / grid.get_dclines() / grid.get_storages(), same solve as the TEST_CASE above)
    // -- this exercises the 3 HVDC lines (LCC, VSC+droop, VSC no-droop), the SVC, and the storage
    // unit individually, rather than only their combined effect on bus voltages.
    LSGrid grid = ls2g_test::make_exotic_elements_grid();
    const CplxVect v0 = CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
    const CplxVect V = grid.ac_pf(v0, 20, 1e-7);
    REQUIRE(V.size() == 14);

    // SVC1: VOLTAGE mode, target 1.0 pu -- res_p_mw is always 0 (SVCs are reactive-only),
    // res_q_mvar is the reactive injection needed to hold the target
    const auto svc_res = grid.get_svcs().get_res();  // (p_mw, q_mvar, v_kv)
    REQUIRE(std::get<0>(svc_res).size() == 1);
    CHECK(std::get<0>(svc_res)(0) == Approx(0.).margin(1e-7));
    CHECK(std::get<1>(svc_res)(0) == Approx(3.1973936014341695).margin(1e-6));

    // 3 HVDC lines, in the order they were captured: [0] HVDC_LCC, [1] HVDC_VSC_DROOP,
    // [2] HVDC_VSC_NODROOP (see init_hvdc_lines' type1/type2/droop_enabled arrays)
    const auto dc_res1 = grid.get_dcline_res1();  // side 1 (p_mw, q_mvar, v_kv)
    const auto dc_res2 = grid.get_dcline_res2();  // side 2
    REQUIRE(std::get<0>(dc_res1).size() == 3);

    RealVect expected_p1(3), expected_q1(3), expected_p2(3), expected_q2(3);
    expected_p1 << -1, 0.27049315600354495, -2;
    expected_q1 << -0.48432217236497815, 9.4080026411983049, -34.899112612469736;
    expected_p2 << 0.97140321016521081, -0.27707090360682907, 1.9293708416040272;
    expected_q2 << -0.47047211298952835, -31.247308364599498, -43.520804752267765;

    for (Eigen::Index i = 0; i < 3; ++i) {
        CHECK(std::get<0>(dc_res1)(i) == Approx(expected_p1(i)).margin(1e-6));
        CHECK(std::get<1>(dc_res1)(i) == Approx(expected_q1(i)).margin(1e-6));
        CHECK(std::get<0>(dc_res2)(i) == Approx(expected_p2(i)).margin(1e-6));
        CHECK(std::get<1>(dc_res2)(i) == Approx(expected_q2(i)).margin(1e-6));
    }

    // BATT1: fixed setpoint (not solved for), load-sign convention -- negative means
    // it is producing (IIDM battery target_p=5.0/target_q=1.0, generator convention,
    // negated on conversion, see _aux_add_storage.py)
    const auto storage_res = grid.get_storages_res();  // (p_mw, q_mvar, v_kv)
    REQUIRE(std::get<0>(storage_res).size() == 1);
    CHECK(std::get<0>(storage_res)(0) == Approx(-5.).margin(1e-7));
    CHECK(std::get<1>(storage_res)(0) == Approx(-1.).margin(1e-7));
}
