// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's TransformerVoltageControl outer loop (TransformerVoltageControlLoop + the
// BranchControl extension's ratios): a transformer solved for its ratio, then rounded to a tap,
// on one symbolic analysis, the inputs untouched. C++14 only.

#include <cmath>
#include <memory>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/TransformerVoltageControlLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolationType;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::TransformerVoltageControlLoop;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// bus 0 (138 kV, the slack generator), a transformer 0 -> 1 (20 kV) with a ratio tap changer
// (rho 0.9 .. 1.1 in `n_steps` positions, at 1.0) regulating bus 1, a line 1 -> 2 and a load at 2
LSGrid make_grid(int n_steps, real_type target_vm, real_type deadband_pu, bool regulating = true)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    RealVect vn(3);
    vn << 138., 20., 20.;
    grid.init_bus(3, 1, vn, 0, 0);
    Eigen::VectorXi line_from(1), line_to(1);
    line_from << 1;
    line_to << 2;
    grid.init_powerlines(RealVect::Constant(1, 0.02), RealVect::Constant(1, 0.1), CplxVect::Zero(1), line_from, line_to);
    RealVect trafo_r(1), trafo_x(1), trafo_ratio(1), trafo_shift(1);
    trafo_r << 0.005;
    trafo_x << 0.08;
    trafo_ratio << 1.;
    trafo_shift << 0.;
    Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
    trafo_bus1 << 0;
    trafo_bus2 << 1;
    grid.init_trafo(trafo_r, trafo_x, CplxVect::Zero(1), trafo_ratio, trafo_shift, {false}, trafo_bus1, trafo_bus2, false);
    RealVect load_p(1), load_q(1);
    load_p << 40.;
    load_q << 20.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.0;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);

    std::vector<real_type> rho, zero(static_cast<std::size_t>(n_steps), 0.);
    for (int k = 0; k < n_steps; ++k) rho.push_back(0.9 + 0.2 * k / (n_steps - 1));
    grid.set_trafo_ratio_tap_changer(0, 0, (n_steps - 1) / 2, rho, zero, zero, zero, zero);
    grid.set_trafo_ratio_tap_regulation(0, regulating, target_vm, deadband_pu, 1);
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

real_type vm_at_tap(int n_steps, int position)
{
    LSGrid grid = make_grid(n_steps, 1., 0., false);
    grid.change_trafo_ratio_tap(0, position);
    grid.change_algorithm("NRSing_SparseLU");
    const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
    REQUIRE(V.size() == 3);
    return std::abs(V(1));
}

}  // namespace

TEST_CASE("a transformer is solved for its ratio, then rounded to the closest tap", "[outer][transformer_voltage_control]")
{
    const int n = 2001;
    LSGrid grid = make_grid(n, 1.0, 0.002);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
    const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
    REQUIRE(V.size() == 3);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    const ls2g::TrafoInfo trafo = grid.get_trafos()[0];
    CHECK(trafo.ratio_tap_position == (n - 1) / 2);  // the input stays
    CHECK(trafo.res_ratio_tap_position != (n - 1) / 2);
    // on a fine table, the voltage lands within a tap's worth of its target
    CHECK(std::abs(V(1)) == Approx(1.0).margin(2e-4));
    // the same as the transformer moved to that tap
    CHECK(std::abs(V(1)) == Approx(vm_at_tap(n, trafo.res_ratio_tap_position)).margin(1e-10));
}

TEST_CASE("a target out of reach ends at the extreme tap", "[outer][transformer_voltage_control]")
{
    LSGrid grid = make_grid(21, 1.3, 0.002);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    CHECK(grid.get_algo().get_outer_loop_stats().status == OuterLoopStatus::STABLE);
    // the tap is on side 2: a higher voltage there takes the lowest ratio
    const int pos = grid.get_trafos()[0].res_ratio_tap_position;
    CHECK((pos == 0 || pos == 20));
    CHECK(vm_at_tap(21, pos) >= vm_at_tap(21, pos == 0 ? 20 : 0));
}

TEST_CASE("within its deadband nothing moves, and nothing is reported", "[outer][transformer_voltage_control]")
{
    LSGrid ref = make_grid(21, 1.0, 0., false);
    ref.change_algorithm("NRSing_SparseLU");
    const CplxVect V0 = ref.ac_pf(flat(ref), 30, 1e-12);
    REQUIRE(V0.size() == 3);
    const real_type v1 = std::abs(V0(1));

    // a deadband around the voltage it has
    LSGrid grid = make_grid(21, v1 + 0.004, 0.01);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    CHECK(grid.get_algo().get_outer_loop_stats().nb_outer_iterations == 0);
    CHECK(grid.get_trafos()[0].res_ratio_tap_position == 10);
    std::size_t found = 0;
    for (const auto & v : grid.get_physical_violations(true, 1e-3, 1e-6)) {
        if (v.violation_type == LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND) ++found;
    }
    CHECK(found == 0);

    // outside it: reported on the regulated bus, a CONTROL
    LSGrid off = make_grid(21, v1 + 0.02, 0.01);
    off.change_algorithm("NRSing_SparseLU");
    off.clear_outer_loops();
    off.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
    REQUIRE(off.ac_pf(flat(off), 30, 1e-12).size() == 3);
    std::vector<ls2g::LimitViolation> viols;
    for (const auto & v : off.get_physical_violations(true, 1e-3, 1e-6)) {
        if (v.violation_type == LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND) viols.push_back(v);
    }
    REQUIRE(viols.size() == 1);
    CHECK(viols[0].element_id == 1);
    CHECK(viols[0].value == Approx(v1 * 20.).epsilon(1e-9));
    CHECK(viols[0].limit == Approx((v1 + 0.02) * 20.).epsilon(1e-12));
    CHECK(viols[0].category() == ls2g::ViolationCategory::CONTROL);
}
