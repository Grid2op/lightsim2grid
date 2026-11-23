// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ScenarioSweep::set_gen_slack_weights / set_storage_slack_weights: a row solves, and splits
// the slack onto its generators, as the grid would with its weights (add_gen_slackbus /
// remove_gen_slackbus / add_storage_slackbus), AC and DC, in one symbolic analysis; what the
// batch cannot do with them is refused (the bounded redistribution with them is checked
// against the one-off path in test_ScenarioSweep_slack_weights.py). C++14 only.

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::cplx_type;
using ls2g::real_type;
using WeightMat = ls2g::YbusPolicy::ParamRealMat;

namespace {

const real_type NaN = std::numeric_limits<real_type>::quiet_NaN();

// a meshed 3-bus grid: generator 0 on bus 0 and generator 1 on bus 1 share the slack (weights
// 1 and 2) with a storage unit on bus 1 (weight 1); generator 2 on bus 2 takes no share; a
// load on bus 2
LSGrid make_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(3, 1, RealVect::Constant(3, 138.), 0, 0);
    Eigen::VectorXi line_from(3), line_to(3);
    line_from << 0, 1, 0;
    line_to << 1, 2, 2;
    RealVect line_r(3), line_x(3);
    line_r << 0.01, 0.02, 0.015;
    line_x << 0.1, 0.12, 0.08;
    grid.init_powerlines(line_r, line_x, CplxVect::Zero(3), line_from, line_to);
    Eigen::VectorXi load_bus(1);
    load_bus << 2;
    grid.init_loads(RealVect::Constant(1, 120.), RealVect::Constant(1, 30.), load_bus);
    RealVect gen_p(3), gen_v(3);
    gen_p << 30., 40., 20.;
    gen_v << 1.02, 1.01, 1.0;
    Eigen::VectorXi gen_bus(3);
    gen_bus << 0, 1, 2;
    grid.init_generators(gen_p, gen_v, RealVect::Constant(3, -1000.), RealVect::Constant(3, 1000.), gen_bus);
    Eigen::VectorXi sto_bus(1);
    sto_bus << 1;
    grid.init_storages(RealVect::Constant(1, -5.), RealVect::Zero(1), sto_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.add_gen_slackbus(1, 2.);
    grid.add_storage_slackbus(0, 1.);
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

// the largest gap between two voltage vectors once a common rotation is taken out: the batch
// keeps the angle reference it picked for the whole sweep, the grid takes its first slack bus
real_type gap_up_to_rotation(const CplxVect & a, const CplxVect & b)
{
    const cplx_type rot = b.dot(a);  // conj(b) . a
    const cplx_type unit = std::abs(rot) > 0. ? rot / std::abs(rot) : cplx_type(1., 0.);
    return (a - b * unit).cwiseAbs().maxCoeff();
}

// row `row` of the weights given to the units `ids` (NaN: the grid's), applied to `grid` the
// way LSGrid takes them
void apply_row(LSGrid & grid, bool storage, const IntVect & ids, const WeightMat & w, int row)
{
    for (Eigen::Index k = 0; k < ids.size(); ++k) {
        const real_type val = w(row, k);
        if (std::isnan(val)) continue;
        if (storage) {
            if (val > 0.) grid.add_storage_slackbus(ids(k), val);
            else grid.remove_storage_slackbus(ids(k));
        } else {
            if (val > 0.) grid.add_gen_slackbus(ids(k), val);
            else grid.remove_gen_slackbus(ids(k));
        }
    }
}

}  // namespace

TEST_CASE("a row's slack weights solve as the grid given them, in one analysis", "[batch][scenario_sweep][slack_weights]")
{
    for (bool ac : {true, false}) {
        INFO("ac = " << ac);
        LSGrid grid = make_grid();
        const AlgorithmType algo = ac ? AlgorithmType::NR_SparseLU : AlgorithmType::DC_SparseLU;
        grid.change_algorithm(algo);
        const int nb_rows = 5;
        ScenarioSweep sweep(grid);
        sweep.change_algorithm(algo);
        sweep.set_compute_gen_results(true);
        IntVect gen_ids(2), sto_ids(1);
        gen_ids << 0, 1;
        sto_ids << 0;
        // row 0: nothing; row 1: generator 1's weight; row 2: generator 0 out of the slack
        // (the reference bus keeps the angle, with no share); row 3: the storage unit out,
        // generator 1 lighter; row 4: the storage unit heavier
        WeightMat gen_w = WeightMat::Constant(nb_rows, 2, NaN), sto_w = WeightMat::Constant(nb_rows, 1, NaN);
        gen_w(1, 1) = 10.;
        gen_w(2, 0) = 0.;
        gen_w(3, 1) = 0.5;
        sto_w(3, 0) = 0.;
        sto_w(4, 0) = 6.;
        sweep.set_gen_slack_weights(gen_ids, gen_w);
        sweep.set_storage_slack_weights(sto_ids, sto_w);
        sweep.compute(flat(grid), 30, 1e-11);
        REQUIRE(sweep.get_voltages().rows() == nb_rows);
        for (const auto & stats : sweep.get_linear_solver_stats_per_algo()) CHECK(stats.nb_analyze <= 1);

        const auto & gen_res = sweep.get_gen_results();
        RealVect first_row_p;
        for (int row = 0; row < nb_rows; ++row) {
            INFO("row " << row);
            REQUIRE(sweep.converged_mask()[static_cast<std::size_t>(row)]);
            LSGrid ref = grid;
            apply_row(ref, false, gen_ids, gen_w, row);
            apply_row(ref, true, sto_ids, sto_w, row);
            const CplxVect V = ac ? ref.ac_pf(flat(ref), 30, 1e-11) : ref.dc_pf(flat(ref), 30, 1e-11);
            REQUIRE(V.size() == 3);
            CHECK(gap_up_to_rotation(sweep.get_voltages().row(row).transpose(), V) < 1e-9);
            const RealVect ref_p(std::get<0>(ref.get_gen_res()));
            RealVect row_p(3);
            for (int g = 0; g < 3; ++g) row_p(g) = gen_res(row, 2 * g);
            CHECK((row_p - ref_p).cwiseAbs().maxCoeff() < 1e-7);
            if (row == 0) first_row_p = row_p;
            else CHECK((row_p - first_row_p).cwiseAbs().maxCoeff() > 1e-3);  // the weights did act
        }
    }
}

TEST_CASE("a row's slack weights are checked when given and at compute", "[batch][scenario_sweep][slack_weights]")
{
    LSGrid grid = make_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    IntVect ids(1);
    ids << 0;
    {
        ScenarioSweep sweep(grid);
        IntVect bad_id(1);
        bad_id << 3;
        CHECK_THROWS_AS(sweep.set_gen_slack_weights(bad_id, WeightMat::Constant(2, 1, 1.)), std::invalid_argument);
        IntVect twice(2);
        twice << 0, 0;
        CHECK_THROWS_AS(sweep.set_gen_slack_weights(twice, WeightMat::Constant(2, 2, 1.)), std::invalid_argument);
        CHECK_THROWS_AS(sweep.set_gen_slack_weights(ids, WeightMat::Constant(2, 2, 1.)), std::invalid_argument);
        CHECK_THROWS_AS(sweep.set_gen_slack_weights(ids, WeightMat::Constant(2, 1, -1.)), std::invalid_argument);
        CHECK_THROWS_AS(sweep.set_gen_slack_weights(ids, WeightMat::Constant(2, 1, std::numeric_limits<real_type>::infinity())),
                        std::invalid_argument);
        // nothing was registered: another count of rows is still free
        CHECK_NOTHROW(sweep.set_gen_slack_weights(ids, WeightMat::Constant(3, 1, NaN)));
        CHECK_THROWS(sweep.set_storage_slack_weights(ids, WeightMat::Constant(2, 1, NaN)));
    }
    SECTION("a weight on a bus outside the grid's slack") {
        // generator 2 stands on bus 2, where no unit shares the grid's slack
        ScenarioSweep sweep(grid);
        IntVect gen2(1);
        gen2 << 2;
        WeightMat w(2, 1);
        w << NaN, 1.;
        sweep.set_gen_slack_weights(gen2, w);
        CHECK_THROWS_AS(sweep.compute(flat(grid), 30, 1e-10), std::invalid_argument);
        // 0 changes nothing to the slack set: accepted
        w << NaN, 0.;
        sweep.set_gen_slack_weights(gen2, w);
        CHECK_NOTHROW(sweep.compute(flat(grid), 30, 1e-10));
    }
    SECTION("with redistribute_slack: accepted, its pre-pass shares with the row's weights") {
        ScenarioSweep sweep(grid);
        sweep.set_redistribute_slack(true);
        sweep.set_gen_slack_weights(ids, WeightMat::Constant(2, 1, 3.));
        CHECK_NOTHROW(sweep.compute(flat(grid), 30, 1e-10));
        CHECK(sweep.converged_mask()[0]);
        CHECK(sweep.converged_mask()[1]);
    }
}
