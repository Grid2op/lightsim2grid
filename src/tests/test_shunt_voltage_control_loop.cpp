// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// OpenLoadFlow's ShuntVoltageControl outer loop (ShuntVoltageControlLoop + the ShuntControl
// extension): a bus' shunts solved for their susceptance, then rounded to sections, on one
// symbolic analysis, the inputs untouched. C++14 only.

#include <cmath>
#include <memory>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "powerflow_algorithm/outer_loop/ShuntVoltageControlLoop.hpp"
#include "powerflow_algorithm/outer_loop/TransformerVoltageControlLoop.hpp"

using Catch::Approx;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::OuterLoopStatus;
using ls2g::RealVect;
using ls2g::ShuntVoltageControlLoop;
using ls2g::TransformerVoltageControlLoop;
using ls2g::real_type;

namespace {

// buses 0-1-2 in a line, a load at 2, the slack generator at 0; a shunt at bus 2 with
// `n` sections of 5 MVar (producing), `on` of them on, regulating bus 2 at `target` pu
LSGrid make_grid(int n, int on, real_type target, bool regulating = true)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi from_id(2), to_id(2);
    from_id << 0, 1;
    to_id << 1, 2;
    grid.init_powerlines(RealVect::Constant(2, 0.01), RealVect::Constant(2, 0.1), CplxVect::Zero(2), from_id, to_id);
    RealVect load_p(1), load_q(1);
    load_p << 60.;
    load_q << 40.;
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(load_p, load_q, load_bus);
    RealVect shunt_p(1), shunt_q(1);
    shunt_p << 0.;
    shunt_q << 0.;
    Eigen::VectorXi shunt_bus(1);
    shunt_bus << 2;
    grid.init_shunt(shunt_p, shunt_q, shunt_bus);
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.0;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    std::vector<real_type> p(static_cast<std::size_t>(n), 0.), q;
    for (int k = 1; k <= n; ++k) q.push_back(-5. * k);  // MVar at 1 pu, absorbing positive
    grid.set_shunt_sections(0, on, p, q);
    grid.set_shunt_section_regulation(0, regulating, target, 0., 2);
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

real_type vm_with(int n, int on)
{
    LSGrid grid = make_grid(n, on, 1., false);
    grid.change_algorithm("NRSing_SparseLU");
    const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
    REQUIRE(V.size() == 3);
    return std::abs(V(2));
}

}  // namespace

TEST_CASE("a shunt is solved for its susceptance, then rounded to the closest section", "[outer][shunt_voltage_control]")
{
    const int n = 40;
    LSGrid grid = make_grid(n, 0, 0.99);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ShuntVoltageControlLoop>());
    const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
    REQUIRE(V.size() == 3);
    const ls2g::OuterLoopStats stats = grid.get_algo().get_outer_loop_stats();
    CHECK(stats.status == OuterLoopStatus::STABLE);
    REQUIRE(stats.loop_iterations.size() == 1);
    CHECK(stats.loop_iterations[0].second == 1);
    CHECK(grid.get_algo().get_linear_solver_stats().nb_analyze == 1);
    const ls2g::ShuntInfo shunt = grid.get_shunts()[0];
    CHECK(shunt.section_count == 0);  // the input stays
    REQUIRE(shunt.res_section_count > 0);
    CHECK(shunt.res_q_mvar < 0.);  // producing
    // the same as the shunt switched to that count, and the closest one to the target
    CHECK(std::abs(V(2)) == Approx(vm_with(n, shunt.res_section_count)).margin(1e-10));
    CHECK(std::abs(std::abs(V(2)) - 0.99) <= std::abs(vm_with(n, shunt.res_section_count - 1) - 0.99) + 1e-12);
    CHECK(std::abs(std::abs(V(2)) - 0.99) <= std::abs(vm_with(n, shunt.res_section_count + 1) - 0.99) + 1e-12);
}

TEST_CASE("a non-regulating shunt does not run the loop", "[outer][shunt_voltage_control]")
{
    LSGrid grid = make_grid(10, 2, 0.99, false);
    grid.change_algorithm("NROuter_SparseLU");
    grid.clear_outer_loops();
    grid.add_outer_loop(std::make_shared<ShuntVoltageControlLoop>());
    const CplxVect V = grid.ac_pf(flat(grid), 30, 1e-12);
    REQUIRE(V.size() == 3);
    CHECK(grid.get_algo().get_outer_loop_stats().loop_iterations.empty());
    CHECK(grid.get_shunts()[0].res_section_count == 2);
    CHECK(std::abs(V(2)) == Approx(vm_with(10, 2)).margin(1e-12));
}

namespace {

// bus 0 (138 kV, the slack generator), a transformer 0 -> 1 (20 kV) whose ratio tap changer
// regulates bus 1 at 1.02 pu, a line 1 -> 2 and a load at 2; a shunt at bus 1 regulating the
// same bus at 1.05 pu. Neither target is met at the grid's tap and sections.
LSGrid make_shared_bus_grid()
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
    RealVect shunt_p(1), shunt_q(1);
    shunt_p << 0.;
    shunt_q << 0.;
    Eigen::VectorXi shunt_bus(1);
    shunt_bus << 1;
    grid.init_shunt(shunt_p, shunt_q, shunt_bus);
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.0;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);

    const int n_steps = 21;
    std::vector<real_type> rho, zero(static_cast<std::size_t>(n_steps), 0.);
    for (int k = 0; k < n_steps; ++k) rho.push_back(0.9 + 0.2 * k / (n_steps - 1));
    grid.set_trafo_ratio_tap_changer(0, 0, (n_steps - 1) / 2, rho, zero, zero, zero, zero);
    grid.set_trafo_ratio_tap_regulation(0, true, 1.02, 0., 1);
    std::vector<real_type> p(10, 0.), q;
    for (int k = 1; k <= 10; ++k) q.push_back(-5. * k);
    grid.set_shunt_sections(0, 0, p, q);
    grid.set_shunt_section_regulation(0, true, 1.05, 0., 1);
    return grid;
}

int count_of(const std::vector<ls2g::LimitViolation> & viols, ls2g::LimitViolationType type)
{
    int res = 0;
    for (const auto & v : viols) res += (v.violation_type == type) ? 1 : 0;
    return res;
}

}  // namespace

TEST_CASE("detection: a transformer voltage control listed first hides the shunt one of its bus",
          "[outer][shunt_voltage_control][detection]")
{
    // the solve leaves such a bus to the transformers (ShuntVoltageControlLoop::_declare): the
    // detection reports what the loops would act on, so the same
    LSGrid grid = make_shared_bus_grid();
    grid.change_algorithm("NRSing_SparseLU");
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-12).size() == 3);
    SECTION("the transformer's loop first: only it is reported") {
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
        grid.add_outer_loop(std::make_shared<ShuntVoltageControlLoop>());
        const std::vector<ls2g::LimitViolation> viols = grid.get_physical_violations();
        CHECK(count_of(viols, ls2g::LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND) == 1);
        CHECK(count_of(viols, ls2g::LimitViolationType::SHUNT_VOLTAGE_CONTROL) == 0);
    }
    SECTION("the shunt's loop first: both are") {
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<ShuntVoltageControlLoop>());
        grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
        const std::vector<ls2g::LimitViolation> viols = grid.get_physical_violations();
        CHECK(count_of(viols, ls2g::LimitViolationType::TRANSFORMER_VOLTAGE_DEADBAND) == 1);
        CHECK(count_of(viols, ls2g::LimitViolationType::SHUNT_VOLTAGE_CONTROL) == 1);
    }
    SECTION("without the transformer's loop, the shunt is reported") {
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<ShuntVoltageControlLoop>());
        CHECK(count_of(grid.get_physical_violations(), ls2g::LimitViolationType::SHUNT_VOLTAGE_CONTROL) == 1);
    }
}
