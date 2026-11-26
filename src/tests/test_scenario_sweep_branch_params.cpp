// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ScenarioSweep::set_line_parameters / set_trafo_parameters: a row solves as the grid would
// with its parameters (a transformer's being the neutral ones, its tap changer's corrections
// applied on top), in one symbolic analysis; a row changing a transformer whose taps an outer
// loop regulates is left out of the batch. C++14 only.

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "powerflow_algorithm/outer_loop/TransformerVoltageControlLoop.hpp"

using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::TransformerVoltageControlLoop;
using ls2g::cplx_type;
using ls2g::real_type;
using RealMat = ls2g::YbusPolicy::ParamRealMat;
using CplxMat = ls2g::YbusPolicy::ParamCplxMat;

namespace {

const real_type NaN = std::numeric_limits<real_type>::quiet_NaN();
const real_type R_TR = 0.005, X_TR = 0.08;

// bus 0 (138 kV, the slack generator), a transformer 0 -> 1 (20 kV) with a ratio tap changer
// (5 positions, at the middle one, whose step corrects r and x) regulating bus 1 when
// `regulating`, a line 1 -> 2 and a load at 2
LSGrid make_grid(bool regulating)
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
    trafo_r << R_TR;
    trafo_x << X_TR;
    trafo_ratio << 1.;
    trafo_shift << 0.;
    Eigen::VectorXi trafo_bus1(1), trafo_bus2(1);
    trafo_bus1 << 0;
    trafo_bus2 << 1;
    grid.init_trafo(trafo_r, trafo_x, CplxVect::Zero(1), trafo_ratio, trafo_shift, {false}, trafo_bus1, trafo_bus2, false);
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(RealVect::Constant(1, 40.), RealVect::Constant(1, 20.), load_bus);
    Eigen::VectorXi gen_bus(1);
    gen_bus << 0;
    grid.init_generators(RealVect::Zero(1), RealVect::Constant(1, 1.0), RealVect::Constant(1, -1000.),
                         RealVect::Constant(1, 1000.), gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.set_trafo_ratio_tap_changer(0, 0, 2, {0.9, 0.95, 1., 1.05, 1.1}, {-20., -10., 10., 20., 30.},
                                     {-10., -5., 10., 15., 20.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.});
    grid.set_trafo_ratio_tap_regulation(0, regulating, 1.0, 0.002, 1);
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

}  // namespace

TEST_CASE("a row's branch parameters solve as the grid given them, in one analysis", "[batch][scenario_sweep][branch_params]")
{
    LSGrid grid = make_grid(false);
    grid.change_algorithm("NRSing_SparseLU");
    const int nb_rows = 4;
    ScenarioSweep sweep(grid);
    IntVect line_ids(1), trafo_ids(1);
    line_ids << 0;
    trafo_ids << 0;
    // row 0: nothing; row 1: the line's x; row 2: the transformer's neutral x and its ratio;
    // row 3: both, and the line's shunt admittance
    RealMat line_r = RealMat::Constant(nb_rows, 1, NaN), line_x = RealMat::Constant(nb_rows, 1, NaN);
    CplxMat line_h1 = CplxMat::Constant(nb_rows, 1, {NaN, NaN}), line_h2 = line_h1;
    line_x(1, 0) = 0.2;
    line_x(3, 0) = 0.05;
    line_h1(3, 0) = cplx_type(0., 0.01);
    RealMat tr_r = RealMat::Constant(nb_rows, 1, NaN), tr_x = tr_r, tr_ratio = tr_r, tr_shift = tr_r;
    CplxMat tr_h1 = CplxMat::Constant(nb_rows, 1, {NaN, NaN}), tr_h2 = tr_h1;
    tr_x(2, 0) = 1.5 * X_TR;
    tr_ratio(2, 0) = 1.02;
    tr_x(3, 0) = 0.5 * X_TR;
    tr_shift(3, 0) = 0.05;
    sweep.set_line_parameters(line_ids, line_r, line_x, line_h1, line_h2);
    sweep.set_trafo_parameters(trafo_ids, tr_r, tr_x, tr_h1, tr_h2, tr_ratio, tr_shift);
    sweep.compute(flat(grid), 30, 1e-12);
    REQUIRE(sweep.get_voltages().rows() == nb_rows);
    for (const auto & stats : sweep.get_linear_solver_stats_per_algo()) CHECK(stats.nb_analyze <= 1);

    for (int row = 0; row < nb_rows; ++row) {
        INFO("row " << row);
        REQUIRE(sweep.converged_mask()[static_cast<std::size_t>(row)]);
        // the oracle: the grid given the row's values the way LSGrid takes them (the
        // transformer's x a NEUTRAL value: its tap's correction applies on top)
        LSGrid ref = grid;
        RealVect r(1), x(1);
        CplxVect h1(1), h2(1);
        r << 0.02;
        x << (std::isnan(line_x(row, 0)) ? 0.1 : line_x(row, 0));
        h1 << (std::isnan(std::real(line_h1(row, 0))) ? cplx_type(0., 0.) : line_h1(row, 0));
        h2 << cplx_type(0., 0.);
        ref.update_powerlines_parameters(r, x, h1, h2);
        RealVect tr(1), tx(1);
        CplxVect th(1);
        tr << R_TR;
        tx << (std::isnan(tr_x(row, 0)) ? X_TR : tr_x(row, 0));
        th << cplx_type(0., 0.);
        ref.update_trafos_parameters_full(tr, tx, th, th);
        if (!std::isnan(tr_ratio(row, 0))) ref.change_ratio_trafo(0, tr_ratio(row, 0));
        if (!std::isnan(tr_shift(row, 0))) ref.change_shift_trafo(0, tr_shift(row, 0));
        const CplxVect V = ref.ac_pf(flat(ref), 30, 1e-12);
        REQUIRE(V.size() == 3);
        CHECK((sweep.get_voltages().row(row).transpose() - V).cwiseAbs().maxCoeff() < 1e-10);
    }
    // the transformer's own x at its tap is the corrected one: the row's neutral value was
    // corrected too (a row giving the CORRECTED value would differ)
    CHECK(grid.get_trafos()[0].x_pu != X_TR);
}

TEST_CASE("a row's branch parameters are checked when given", "[batch][scenario_sweep][branch_params]")
{
    LSGrid grid = make_grid(false);
    ScenarioSweep sweep(grid);
    IntVect ids(1);
    ids << 0;
    const RealMat nan2 = RealMat::Constant(2, 1, NaN);
    const CplxMat cnan2 = CplxMat::Constant(2, 1, {NaN, NaN});
    IntVect bad_id(1);
    bad_id << 1;
    CHECK_THROWS_AS(sweep.set_line_parameters(bad_id, nan2, nan2, cnan2, cnan2), std::invalid_argument);
    IntVect twice(2);
    twice << 0, 0;
    CHECK_THROWS_AS(sweep.set_line_parameters(twice, RealMat::Constant(2, 2, NaN), RealMat::Constant(2, 2, NaN),
                                              CplxMat::Constant(2, 2, {NaN, NaN}), CplxMat::Constant(2, 2, {NaN, NaN})),
                    std::invalid_argument);
    CHECK_THROWS_AS(sweep.set_line_parameters(ids, nan2, RealMat::Constant(3, 1, NaN), cnan2, cnan2), std::invalid_argument);
    RealMat zero = RealMat::Zero(2, 1);
    CHECK_THROWS_AS(sweep.set_line_parameters(ids, zero, zero, cnan2, cnan2), std::invalid_argument);
    RealMat inf = RealMat::Constant(2, 1, std::numeric_limits<real_type>::infinity());
    CHECK_THROWS_AS(sweep.set_line_parameters(ids, nan2, inf, cnan2, cnan2), std::invalid_argument);
    CHECK_THROWS_AS(sweep.set_trafo_parameters(ids, nan2, nan2, cnan2, cnan2, RealMat::Constant(2, 1, -1.), nan2),
                    std::invalid_argument);
    // nothing was registered: another count of rows is still free
    const RealMat nan3 = RealMat::Constant(3, 1, NaN);
    const CplxMat cnan3 = CplxMat::Constant(3, 1, {NaN, NaN});
    CHECK_NOTHROW(sweep.set_line_parameters(ids, nan3, nan3, cnan3, cnan3));
    CHECK_THROWS(sweep.set_trafo_parameters(ids, nan2, nan2, cnan2, cnan2, nan2, nan2));
}

TEST_CASE("a row changing a transformer whose taps a loop regulates is left out of the batch",
          "[batch][scenario_sweep][branch_params][outer]")
{
    for (bool regulating : {true, false}) {
        INFO("regulating = " << regulating);
        LSGrid grid = make_grid(regulating);
        grid.change_algorithm("NROuter_SparseLU");
        grid.clear_outer_loops();
        grid.add_outer_loop(std::make_shared<TransformerVoltageControlLoop>());
        ScenarioSweep sweep(grid);
        IntVect line_ids(1), trafo_ids(1);
        line_ids << 0;
        trafo_ids << 0;
        // row 0 changes the line only, row 1 the transformer
        RealMat line_x(2, 1);
        line_x << 0.2, NaN;
        const RealMat nan2 = RealMat::Constant(2, 1, NaN);
        const CplxMat cnan2 = CplxMat::Constant(2, 1, {NaN, NaN});
        RealMat tr_x(2, 1);
        tr_x << NaN, 1.2 * X_TR;
        sweep.set_line_parameters(line_ids, nan2, line_x, cnan2, cnan2);
        sweep.set_trafo_parameters(trafo_ids, nan2, tr_x, cnan2, cnan2, nan2, nan2);
        sweep.compute(flat(grid), 30, 1e-10);
        const std::vector<std::string> & reasons = sweep.get_batch_unsupported_reasons();
        REQUIRE(reasons.size() == 2);
        CHECK(reasons[0].empty());
        CHECK(sweep.converged_mask()[0]);
        if (regulating) {
            CHECK(reasons[1].find("TransformerVoltageControl") != std::string::npos);
            CHECK_FALSE(sweep.converged_mask()[1]);
        } else {
            CHECK(reasons[1].empty());
            CHECK(sweep.converged_mask()[1]);
        }
    }
}
