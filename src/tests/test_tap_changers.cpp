// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The tap changers of the transformers and the sections of the shunts: the pi model taken at
// the taps (as OpenLoadFlow's Transformers.getTapCharacteristics), moving them, rounding to
// the closest tap, and their state. C++14 only.

#include <cmath>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::RegulationMode;
using ls2g::TrafoContainer;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// buses 0-1-2: a line 0-1, a transformer 1-2 (r, x, b total, ratio, shift in degree), a load
// and a shunt (q_shunt MVar at 1 pu) at bus 2, the slack generator at bus 0
LSGrid make_grid(real_type r, real_type x, real_type b, real_type ratio, real_type shift_deg, real_type q_shunt = 0.)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi line_from(1), line_to(1);
    line_from << 0;
    line_to << 1;
    grid.init_powerlines(RealVect::Constant(1, 0.01), RealVect::Constant(1, 0.1), CplxVect::Zero(1), line_from, line_to);
    RealVect trafo_r(1), trafo_x(1), trafo_ratio(1), trafo_shift(1);
    CplxVect trafo_b(1);
    trafo_r << r;
    trafo_x << x;
    trafo_b << cplx_type(0., b);
    trafo_ratio << ratio;
    trafo_shift << shift_deg;
    Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
    trafo_bus1 << 1;
    trafo_bus2 << 2;
    grid.init_trafo(trafo_r, trafo_x, trafo_b, trafo_ratio, trafo_shift, {true}, trafo_bus1, trafo_bus2, false);
    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 20.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect shunt_p(1), shunt_q(1);
    shunt_p << 0.;
    shunt_q << q_shunt;
    Eigen::VectorXi shunt_bus(1);
    shunt_bus << 2;
    grid.init_shunt(shunt_p, shunt_q, shunt_bus);
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

CplxVect solve(LSGrid & grid)
{
    return grid.ac_pf(CplxVect::Constant(3, {1.0, 0.}), 30, 1e-12);
}

const real_type R0 = 0.01, X0 = 0.2, B0 = 0.04;

// a ratio tap changer: positions 0 .. 2, the transformer at position 1
void add_ratio_taps(LSGrid & grid)
{
    grid.set_trafo_ratio_tap_changer(0, 0, 1, {0.95, 1.05, 1.1}, {0., 10., 20.}, {0., 5., 10.},
                                     {0., 0., 0.}, {0., 50., 100.});
}

}  // namespace

TEST_CASE("a ratio tap changer's step gives the pi model", "[tap_changer]")
{
    // given at its current tap: ratio 1.05, the neutral impedance
    LSGrid grid = make_grid(R0, X0, B0, 1.05, 0.);
    add_ratio_taps(grid);
    const ls2g::TrafoInfo info = grid.get_trafos()[0];
    CHECK(info.has_ratio_tap_changer);
    CHECK(info.ratio_tap_position == 1);
    CHECK(info.ratio_low_tap == 0);
    CHECK(info.ratio_high_tap == 2);
    CHECK(!info.has_phase_tap_changer);
    CHECK(info.ratio == Approx(1.05).epsilon(1e-15));
    CHECK(info.r_pu == Approx(R0 * 1.1).epsilon(1e-15));
    CHECK(info.x_pu == Approx(X0 * 1.05).epsilon(1e-15));

    // the same solve as the transformer given with the corrected values
    LSGrid ref = make_grid(R0 * 1.1, X0 * 1.05, B0 * 1.5, 1.05, 0.);
    const CplxVect V = solve(grid);
    REQUIRE(V.size() == 3);
    CHECK((V - solve(ref)).cwiseAbs().maxCoeff() < 1e-12);

    // moved: the ratio and the impedance follow, the Jacobian keeps its pattern
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    REQUIRE(solve(grid).size() == 3);
    grid.change_trafo_ratio_tap(0, 2);
    CHECK(grid.get_trafos()[0].ratio == Approx(1.1 / 1.05 * 1.05).epsilon(1e-15));
    const CplxVect V2 = solve(grid);
    REQUIRE(V2.size() == 3);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    LSGrid ref2 = make_grid(R0 * 1.2, X0 * 1.1, B0 * 2., 1.1, 0.);
    CHECK((V2 - solve(ref2)).cwiseAbs().maxCoeff() < 1e-12);

    // a position outside the table
    CHECK_THROWS(grid.change_trafo_ratio_tap(0, 3));
    CHECK_THROWS(grid.change_trafo_phase_tap(0, 0));  // no phase changer
}

TEST_CASE("a phase tap changer's step gives the shift", "[tap_changer]")
{
    LSGrid grid = make_grid(R0, X0, 0., 1., 10.);
    grid.set_trafo_phase_tap_changer(0, -1, 1, {1., 1., 1.}, {-10., 0., 10.}, {0., 0., 0.}, {-20., 0., 20.},
                                     {0., 0., 0.}, {0., 0., 0.});
    CHECK(grid.get_trafos()[0].shift_rad == Approx(10. / 180. * M_PI).epsilon(1e-15));
    CHECK(grid.get_trafos()[0].x_pu == Approx(X0 * 1.2).epsilon(1e-15));
    grid.change_trafo_phase_tap(0, -1);
    CHECK(grid.get_trafos()[0].shift_rad == Approx(-10. / 180. * M_PI).epsilon(1e-15));
    LSGrid ref = make_grid(R0, X0 * 0.8, 0., 1., -10.);
    const CplxVect V = solve(grid);
    REQUIRE(V.size() == 3);
    CHECK((V - solve(ref)).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("the closest tap keeps the current one unless another is strictly closer", "[tap_changer]")
{
    LSGrid grid = make_grid(R0, X0, 0., 1., 0.);
    grid.set_trafo_phase_tap_changer(0, 0, 1, {1., 1., 1.}, {-10., 0., 10.}, {0., 0., 0.}, {0., 0., 0.},
                                     {0., 0., 0.}, {0., 0., 0.});
    const real_type deg = M_PI / 180.;
    CHECK(grid.closest_trafo_phase_tap(0, 5. * deg) == 1);   // a tie with position 2
    CHECK(grid.closest_trafo_phase_tap(0, 6. * deg) == 2);
    CHECK(grid.closest_trafo_phase_tap(0, -7. * deg) == 0);
    CHECK_THROWS(grid.closest_trafo_ratio_tap(0, 1.));  // no ratio changer

    // a ratio, the phase changer's rho included
    LSGrid grid2 = make_grid(R0, X0, 0., 1.05, 0.);
    add_ratio_taps(grid2);
    CHECK(grid2.closest_trafo_ratio_tap(0, 0.96) == 0);
    CHECK(grid2.closest_trafo_ratio_tap(0, 1.09) == 2);
}

TEST_CASE("a shift between two taps keeps the shift-dependent impedance", "[tap_changer]")
{
    // set_trafo_shift_dependent_rx: at a tap the step decides, between taps the correction
    // follows the shift, relative to the tap's
    LSGrid grid = make_grid(R0, X0, 0., 1., 0.);
    grid.set_trafo_shift_dependent_rx(true, {{-10. / 180. * M_PI, 0., 10. / 180. * M_PI}}, {{-20., 0., 20.}});
    grid.set_trafo_phase_tap_changer(0, 0, 1, {1., 1., 1.}, {-10., 0., 10.}, {-20., 0., 20.}, {-20., 0., 20.},
                                     {0., 0., 0.}, {0., 0., 0.});
    CHECK(grid.get_trafos()[0].x_pu == Approx(X0).epsilon(1e-15));
    grid.change_shift_trafo(0, 5. / 180. * M_PI);
    CHECK(grid.get_trafos()[0].x_pu == Approx(X0 * 1.1).epsilon(1e-12));
    grid.change_trafo_phase_tap(0, 2);
    CHECK(grid.get_trafos()[0].x_pu == Approx(X0 * 1.2).epsilon(1e-15));
}

TEST_CASE("tap changers and sections survive a copy and a state", "[tap_changer]")
{
    LSGrid grid = make_grid(R0, X0, B0, 1.05, 0., 10.);
    add_ratio_taps(grid);
    grid.set_trafo_ratio_tap_regulation(0, true, 1.01, 0.005, 2);
    grid.set_shunt_sections(0, 1, {0., 0., 0.}, {10., 25., 45.});
    grid.set_shunt_section_regulation(0, true, 1.0, 0.01, 2);
    grid.change_trafo_ratio_tap(0, 0);

    LSGrid copy(grid);
    const ls2g::TrafoInfo info = copy.get_trafos()[0];
    CHECK(info.ratio_tap_position == 0);
    CHECK(info.ratio == Approx(0.95).epsilon(1e-15));
    CHECK(info.ratio_regulation_mode == RegulationMode::VOLTAGE);
    CHECK(info.ratio_regulating);
    CHECK(info.ratio_target == 1.01);
    CHECK(info.ratio_deadband == 0.005);
    CHECK(info.ratio_regulated == 2);
    const ls2g::ShuntInfo shunt = copy.get_shunts()[0];
    CHECK(shunt.has_sections);
    CHECK(shunt.section_count == 1);
    CHECK(shunt.max_section_count == 3);
    CHECK(shunt.regulating);
    CHECK(shunt.regulated_bus == 2);

    // the state of the container, and a move back from it
    TrafoContainer trafos;
    TrafoContainer::StateRes state = grid.get_trafos().get_state();
    trafos.set_state(state);
    CHECK(trafos[0].r_pu == Approx(R0).epsilon(1e-15));
    CHECK(trafos.get_tap_changers(false).position(0) == 0);
    // a position outside the table is refused
    std::get<0>(std::get<TrafoContainer::RATIO_TAPS>(state))[0] = 5;
    CHECK_THROWS(trafos.set_state(state));
    CHECK_THROWS(grid.set_trafo_ratio_tap_regulation(0, true, 1., 0., 3));  // no such bus
}

TEST_CASE("a shunt's sections give its reactive power", "[tap_changer][shunt]")
{
    LSGrid grid = make_grid(R0, X0, 0., 1., 0.);
    grid.set_shunt_sections(0, 2, {0., 0., 0.}, {10., 25., 45.});
    CHECK(grid.get_shunts()[0].target_q_mvar == 25.);
    LSGrid ref = make_grid(R0, X0, 0., 1., 0., 25.);
    const CplxVect V = solve(grid);
    REQUIRE(V.size() == 3);
    CHECK((V - solve(ref)).cwiseAbs().maxCoeff() < 1e-12);

    grid.change_shunt_section_count(0, 0);
    CHECK(grid.get_shunts()[0].target_q_mvar == 0.);
    grid.change_shunt_section_count(0, 3);
    CHECK(grid.get_shunts()[0].target_q_mvar == 45.);
    LSGrid ref3 = make_grid(R0, X0, 0., 1., 0., 45.);
    const CplxVect V3 = solve(grid);
    REQUIRE(V3.size() == 3);
    CHECK((V3 - solve(ref3)).cwiseAbs().maxCoeff() < 1e-12);
    CHECK_THROWS(grid.change_shunt_section_count(0, 4));
    CHECK_THROWS(grid.change_shunt_section_count(0, -1));
}
