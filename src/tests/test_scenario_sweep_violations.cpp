// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Covers the two features added to ScenarioSweep on top of what test_batch_voltage_control.cpp
// and test_injection_sweep.cpp already exercise: `handle_disconnected_grid` (main-component-only
// masking, shared with ContingencyAnalysis) and inline limit-violation checking
// (`compute_limit_violations` / `get_violations` / `get_violations_n`). See
// batch_algorithm/BaseBatchSweep.hpp for the SFINAE gating this relies on, and
// batch_algorithm/YbusPolicy.hpp for `branch_ids_for_row`, tested directly below.

#include <complex>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "batch_algorithm/YbusPolicy.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::Coeff;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::LimitViolationType;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::ViolationElementType;
using ls2g::YbusPolicy;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type V_SET = 1.04;
const real_type LOAD_P = 50.;
const real_type LOAD_Q = 10.;
const int NB_BUS = 4;

// same 4-bus radial feeder (0-1-2-3, r=0.01, x=0.1 pu, sn_mva=100) as
// test_batch_voltage_control.cpp's make_skeleton()/make_remote_gen_grid() -- not
// shared across files in this codebase, so duplicated here rather than reaching into
// another test file's anonymous namespace.
LSGrid make_skeleton()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const RealVect bus_vn_kv = RealVect::Constant(NB_BUS, 138.);
    grid.init_bus(static_cast<unsigned int>(NB_BUS), 1, bus_vn_kv, 0, 0);

    const int n_line = NB_BUS - 1;
    const RealVect branch_r = RealVect::Constant(n_line, 0.01);
    const RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < n_line; ++i) { from_id(i) = i; to_id(i) = i + 1; }
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << LOAD_P;
    load_q << LOAD_Q;
    Eigen::VectorXi load_bus(1);
    load_bus << NB_BUS - 1;
    grid.init_loads(load_p, load_q, load_bus);
    return grid;
}

// slack generator at bus 0, plus a generator at bus 1 REMOTELY regulating bus 3 (the
// load bus) at V_SET -- disconnecting line 2 (bus2--bus3) strands the regulated bus.
LSGrid make_remote_gen_grid()
{
    LSGrid grid = make_skeleton();
    RealVect gen_p(2), gen_v(2), gen_min_q(2), gen_max_q(2);
    Eigen::VectorXi gen_bus(2);
    gen_p << 0., 0.;
    gen_v << 1.02, V_SET;
    gen_min_q << -1000., -1000.;
    gen_max_q << 1000., 1000.;
    gen_bus << 0, 1;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);

    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

}  // namespace

TEST_CASE("a contingency stranding a regulated bus is reported through ScenarioSweep, not silently wrong",
          "[batch][scenario_sweep][vctrl]")
{
    // mirrors test_batch_voltage_control.cpp's identically-named ContingencyAnalysis
    // test: line 2 is bus2--bus3, and bus 3 is the REGULATED bus. In the default mode
    // the row is skipped outright (zero voltages, one NOT_SIMULATED sentinel
    // violation). In "handle disconnected grid" mode the row comes back as a
    // DIVERGENCE (masking pins Vm(3), which conflicts with VoltageControl's own
    // bordered row against that same bus) -- not a good answer, but an honest one:
    // what must never happen is a converged row with the regulation quietly dropped.
    for (bool handle_disconnected : {false, true}) {
        INFO("handle_disconnected_grid = " << handle_disconnected);
        LSGrid grid = make_remote_gen_grid();
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        ScenarioSweep sweep(grid);
        // compute_limit_violations must be set BEFORE handle_disconnected_grid (and
        // before any injection/contingency is registered): toggling it clear()s the
        // whole object, by design -- see set_compute_limit_violations's own comment.
        sweep.set_compute_limit_violations(true);
        sweep.set_handle_disconnected_grid(handle_disconnected);

        RealMat load_p(1, 1);
        load_p << LOAD_P;
        sweep.modify_load_p(load_p);
        BoolMat line_mask(1, NB_BUS - 1);
        line_mask << false, false, true;  // disconnect line 2 (bus2--bus3)
        sweep.set_contingency_lines(line_mask);

        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_voltages().rows() == 1);
        REQUIRE(sweep.get_violations().size() == 1);
        REQUIRE(sweep.get_violations()[0].size() == 1);
        CHECK(sweep.get_violations()[0][0].element_type == ViolationElementType::GRID);
        const auto expected_reason = handle_disconnected ? LimitViolationType::DIVERGENCE
                                                          : LimitViolationType::NOT_SIMULATED;
        CHECK(sweep.get_violations()[0][0].violation_type == expected_reason);
        for (int b = 0; b < NB_BUS; ++b) CHECK(std::abs(sweep.get_voltages()(0, b)) == 0.);
    }
}

namespace {
// slack generator at bus 0 only -- NO remote voltage control, so a masked solve of
// the main component has nothing extra to conflict with the masking. On the 0-1-2-3
// radial feeder, disconnecting line 2 (bus2--bus3, the line right before the load)
// strands only a load-only leaf: the slack (bus 0) is never stranded by this cut, so
// the main component always converges normally under masking.
LSGrid make_plain_slack_grid()
{
    LSGrid grid = make_skeleton();
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    Eigen::VectorXi gen_bus(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.tell_solver_need_reset();
    return grid;
}
}  // namespace

TEST_CASE("ScenarioSweep + handle_disconnected_grid matches ContingencyAnalysis on the same contingency",
          "[batch][scenario_sweep][vctrl]")
{
    // a contingency that strands only a load-only leaf (never the slack) -- both
    // should converge, and agree with each other bit for bit (both funnel through
    // the same _select_ref_slack_and_masks(), representation-agnostic over
    // li_defaults vs. line_mask/trafo_mask).
    LSGrid grid = make_plain_slack_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);

    ContingencyAnalysis ca(grid, /*compute_limit_violations=*/true);
    ca.set_handle_disconnected_grid(true);
    ca.add_n1(2);  // disconnect line 2 (bus2--bus3)
    ca.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(ca.converged()[0] == 1);

    ScenarioSweep sweep(grid);
    sweep.set_compute_limit_violations(true);
    sweep.set_handle_disconnected_grid(true);
    RealMat load_p(1, 1);
    load_p << LOAD_P;
    sweep.modify_load_p(load_p);
    BoolMat line_mask(1, NB_BUS - 1);
    line_mask << false, false, true;  // disconnect line 2
    sweep.set_contingency_lines(line_mask);
    sweep.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);

    for (int b = 0; b < NB_BUS; ++b) {
        CHECK(sweep.get_voltages()(0, b).real() == Approx(ca.get_voltages()(0, b).real()));
        CHECK(sweep.get_voltages()(0, b).imag() == Approx(ca.get_voltages()(0, b).imag()));
    }
}

TEST_CASE("YbusPolicy::Contingency::branch_ids_for_row matches init_li_coeffs_from_masks's own scan",
          "[batch][scenario_sweep][ybus_policy]")
{
    YbusPolicy::Contingency policy;
    const size_t n_line = 3;
    // 2 rows: row 0 disconnects line 1 only; row 1 disconnects line 0 and trafo 1.
    BoolMat line_mask(2, 3);
    line_mask << false, true, false,
                 true, false, false;
    BoolMat trafo_mask(2, 2);
    trafo_mask << false, false,
                  false, true;
    policy.line_mask = line_mask;
    policy.trafo_mask = trafo_mask;

    const auto row0 = policy.branch_ids_for_row(0, n_line);
    REQUIRE(row0.size() == 1);
    CHECK(row0[0] == 1);

    const auto row1 = policy.branch_ids_for_row(1, n_line);
    REQUIRE(row1.size() == 2);
    CHECK(row1[0] == 0);
    CHECK(row1[1] == static_cast<int>(n_line) + 1);  // trafo 1 -> branch id n_line+1
}

namespace {
// 4-bus RING (0-1, 1-2, 2-3, 3-0), slack at bus 0, single load at bus 3 -- unlike
// the radial feeder used by the other tests in this file, every line here is
// redundant (2-edge-connected): disconnecting any ONE line still leaves every bus
// reachable from the slack, so a plain (non-masked) single-line contingency always
// converges normally. This is deliberately a different grid from make_skeleton()'s
// chain, where every line is a bridge and this test's premise (no
// handle_disconnected_grid, expect every row to converge) would not hold.
LSGrid make_ring_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    grid.init_bus(static_cast<unsigned int>(NB_BUS), 1, RealVect::Constant(NB_BUS, 138.), 0, 0);

    const int n_line = NB_BUS;
    const RealVect branch_r = RealVect::Constant(n_line, 0.01);
    const RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    for (int i = 0; i < n_line; ++i) { from_id(i) = i; to_id(i) = (i + 1) % NB_BUS; }
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(1), load_q(1);
    load_p << LOAD_P;
    load_q << LOAD_Q;
    Eigen::VectorXi load_bus(1);
    load_bus << NB_BUS - 1;
    grid.init_loads(load_p, load_q, load_bus);

    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    Eigen::VectorXi gen_bus(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    gen_bus << 0;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    return grid;
}
}  // namespace

TEST_CASE("ScenarioSweep multi-row masks: different rows disconnect different branches",
          "[batch][scenario_sweep]")
{
    // 3 rows, each disconnecting a different single line of the 4-line ring --
    // every row must independently match the equivalent single-contingency
    // ContingencyAnalysis result (exercises _row_skip_branch_ids's per-row
    // correctness, not just a single-row case).
    LSGrid grid = make_ring_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);

    ScenarioSweep sweep(grid);
    RealMat load_p(3, 1);
    load_p << LOAD_P, LOAD_P, LOAD_P;
    sweep.modify_load_p(load_p);
    BoolMat line_mask = BoolMat::Constant(3, NB_BUS, false);
    line_mask(0, 0) = true;  // row 0: disconnect line 0
    line_mask(1, 1) = true;  // row 1: disconnect line 1
    // row 2: nothing disconnected
    sweep.set_contingency_lines(line_mask);
    sweep.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);

    for (int row = 0; row < 2; ++row) {
        ContingencyAnalysis ca(grid);
        ca.add_n1(row);
        ca.compute(flat_start(grid), 30, 1e-10);
        for (int b = 0; b < NB_BUS; ++b) {
            CHECK(sweep.get_voltages()(row, b).real() == Approx(ca.get_voltages()(0, b).real()));
            CHECK(sweep.get_voltages()(row, b).imag() == Approx(ca.get_voltages()(0, b).imag()));
        }
    }
    // row 2 (nothing disconnected) must match a plain single-shot ac_pf
    const CplxVect V_ref = grid.ac_pf(flat_start(grid), 30, 1e-10);
    for (int b = 0; b < NB_BUS; ++b) {
        CHECK(sweep.get_voltages()(2, b).real() == Approx(V_ref(b).real()));
        CHECK(sweep.get_voltages()(2, b).imag() == Approx(V_ref(b).imag()));
    }
}

TEST_CASE("a mid-batch islanding contingency does not abandon the rows after it (no handle_disconnected_grid)",
          "[batch][scenario_sweep]")
{
    // Regression test for _run_range's early-return: without handle_disconnected_grid,
    // row 0 disconnects line 2 (bus2--bus3), which strands the load-only leaf bus3 --
    // a genuinely disconnected grid, so row 0 must NOT converge. Rows 1 and 2 register
    // no contingency at all and must still be solved normally: _run_range must
    // `continue` past row 0, not `return` and silently leave rows 1-2 as zero voltages
    // with an empty (mis-readable-as-"converged, nothing found") violation list.
    LSGrid grid = make_plain_slack_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);

    ScenarioSweep sweep(grid);
    sweep.set_compute_limit_violations(true);  // set first: see set_compute_limit_violations's note
    RealMat load_p(3, 1);
    load_p << LOAD_P, LOAD_P, LOAD_P;
    sweep.modify_load_p(load_p);
    BoolMat line_mask = BoolMat::Constant(3, NB_BUS - 1, false);
    line_mask(0, 2) = true;  // row 0: disconnect line 2 (bus2--bus3), strands bus3
    // rows 1, 2: nothing disconnected
    sweep.set_contingency_lines(line_mask);
    sweep.compute(flat_start(grid), 30, 1e-10);

    REQUIRE(sweep.get_status() == 0);  // row 0 failed, so the batch as a whole is not clean

    // row 0: genuinely not simulated / stranded -- zero voltages, one GRID sentinel
    REQUIRE(sweep.get_violations()[0].size() == 1);
    CHECK(sweep.get_violations()[0][0].element_type == ViolationElementType::GRID);
    for (int b = 0; b < NB_BUS; ++b) CHECK(std::abs(sweep.get_voltages()(0, b)) == 0.);

    // rows 1 and 2 must still have been solved -- not left over from row 0's abort --
    // and must match a plain single-shot ac_pf on the unmodified grid exactly.
    const CplxVect V_ref = grid.ac_pf(flat_start(grid), 30, 1e-10);
    for (int row = 1; row < 3; ++row) {
        CHECK(sweep.get_violations()[row].empty());
        bool any_nonzero = false;
        for (int b = 0; b < NB_BUS; ++b) {
            CHECK(sweep.get_voltages()(row, b).real() == Approx(V_ref(b).real()));
            CHECK(sweep.get_voltages()(row, b).imag() == Approx(V_ref(b).imag()));
            if (std::abs(sweep.get_voltages()(row, b)) > 0.) any_nonzero = true;
        }
        CHECK(any_nonzero);
    }
}

TEST_CASE("compute_flows/compute_power_flows reject a mask re-set after compute() without a fresh compute()",
          "[batch][scenario_sweep]")
{
    // Regression test: unlike ContingencyAnalysis (whose own
    // _maybe_check_results_match_defaults compares the registered contingency COUNT),
    // ScenarioSweep's row-count lock means set_contingency_lines can only ever replace
    // a row's CONTENT, never its count -- so re-setting the same-shaped mask after
    // compute() must still be caught, via _results_stale_ instead of a row-count
    // comparison.
    LSGrid grid = make_ring_grid();
    grid.change_algorithm(AlgorithmType::NR_SparseLU);

    ScenarioSweep sweep(grid);
    RealMat load_p(1, 1);
    load_p << LOAD_P;
    sweep.modify_load_p(load_p);
    BoolMat line_mask = BoolMat::Constant(1, NB_BUS, false);
    line_mask(0, 0) = true;  // disconnect line 0
    sweep.set_contingency_lines(line_mask);
    sweep.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);
    CHECK_NOTHROW(sweep.compute_flows());
    CHECK_NOTHROW(sweep.compute_power_flows());

    // same shape, different content -- the row-count lock alone would never notice
    BoolMat line_mask2 = BoolMat::Constant(1, NB_BUS, false);
    line_mask2(0, 1) = true;  // disconnect line 1 instead
    sweep.set_contingency_lines(line_mask2);
    CHECK_THROWS_AS(sweep.compute_flows(), std::runtime_error);
    CHECK_THROWS_AS(sweep.compute_power_flows(), std::runtime_error);

    // a fresh compute() clears the staleness again
    sweep.compute(flat_start(grid), 30, 1e-10);
    CHECK_NOTHROW(sweep.compute_flows());
    CHECK_NOTHROW(sweep.compute_power_flows());
}
