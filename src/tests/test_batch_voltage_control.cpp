// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The NR extensions (Base's free-Vm slack unknowns, Hvdc's angle droop,
// VoltageControl's remote regulation) are not built from the arguments handed to
// the solver: they are pulled out of the LSGrid itself, through the `lsgrid_ptr`
// the algorithm holds, by LSGrid::get_free_vm_slack_solver_buses /
// fill_hvdc_droop_solver_data / fill_voltage_control_solver_data. Those three read
// the grid's *member* solver-side labelling (id_me_to_ac_solver_,
// id_ac_solver_to_me_, slack_bus_id_ac_solver_).
//
// LSGrid::ac_pf passes those very members to pre_process_solver, so they are always
// current there. The batch algorithms (TimeSeries / ContingencyAnalysis /
// SecurityAnalysis, through BaseBatchSolverSynch) do not: they own local mapping
// vectors and hand those over instead, on a private COPY of the grid whose copy
// constructor reset() the members to empty. A member left empty is indistinguishable
// from "this grid has no such controller", so the extension quietly builds nothing
// and the batch returns a converged, plausible, WRONG answer -- remote voltage
// control simply not applied.
//
// Every test here therefore compares a batch solve against the single-shot ac_pf of
// the same grid: they must agree to solver tolerance.

#include <cmath>
#include <complex>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/ContingencyAnalysis.hpp"
#include "batch_algorithm/BaseInjectionSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::GlobalBusId;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::TimeSeries;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type V_SET = 1.04;
const real_type LOAD_P = 50.;
const real_type LOAD_Q = 10.;
const int NB_BUS = 4;

// 4-bus radial feeder 0-1-2-3 (r = 0.01, x = 0.1 pu per line), one
// 50 MW / 10 MVAr load at bus 3, sn_mva = 100. Generation is added per grid.
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

// slack generator at bus 0, plus a generator at bus 1 REMOTELY regulating bus 3
// (the load bus) at V_SET -- the VoltageControl extension.
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

// slack generator at bus 0, plus TWO generators (buses 1 and 2, deliberately
// different reactive ranges) both remotely regulating bus 3 at the same setpoint.
// One control group of two controllers: two Q columns against one voltage row plus
// one sharing row, splitting Q by the (qmax - qmin) key. Sharing a regulated bus is
// supported -- what is not is regulating a bus whose magnitude is already
// determined (the slack, or a PV bus pinned by a LOCAL regulator).
LSGrid make_shared_reg_bus_grid()
{
    LSGrid grid = make_skeleton();
    RealVect gen_p(3), gen_v(3), gen_min_q(3), gen_max_q(3);
    Eigen::VectorXi gen_bus(3);
    gen_p << 0., 0., 0.;
    gen_v << 1.02, V_SET, V_SET;
    gen_min_q << -1000., -1000., -250.;
    gen_max_q << 1000., 1000., 250.;
    gen_bus << 0, 1, 2;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);

    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    grid.set_gen_regulated_bus(2, 3);
    grid.tell_solver_need_reset();
    return grid;
}

// slack generator at bus 0, plus a voltage-mode SVC at bus 2 regulating bus 3
LSGrid make_remote_svc_grid(real_type slope_pu)
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

    RealVect target_vm(1), q_set(1), slope(1), b_min(1), b_max(1);
    target_vm << V_SET;
    q_set << 0.;
    slope << slope_pu;
    b_min << -100.;
    b_max << 100.;
    Eigen::VectorXi reg(1), bus(1);
    reg << 3;
    bus << 2;
    // SvcContainer::RegulationMode: VOLTAGE = 1
    grid.init_svcs({1}, target_vm, q_set, slope, b_min, b_max, reg, bus);
    grid.tell_solver_need_reset();
    return grid;
}

// NO voltage control at all: two DISTRIBUTED-SLACK generators, gen 0 (bus 0) a
// local voltage regulator, gen 1 (bus 1) with voltage_regulator_on == false. Bus 1
// is a slack participant that no local generator voltage-pins, so
// LSGrid::get_free_vm_slack_solver_buses must give it a free Vm unknown + Q
// equation -- the `Base` extension, which reads slack_bus_id_ac_solver_. Without
// it bus 1 stays frozen at init_vm_pu.
LSGrid make_free_vm_slack_grid()
{
    LSGrid grid = make_skeleton();
    RealVect gen_p(2), gen_v(2), gen_q(2), gen_min_q(2), gen_max_q(2);
    Eigen::VectorXi gen_bus(2);
    gen_p << 0., 0.;
    gen_v << 1.02, 1.02;
    gen_q << 0., 0.;
    gen_min_q << -1000., -1000.;
    gen_max_q << 1000., 1000.;
    gen_bus << 0, 1;
    const std::vector<bool> vreg_on{true, false};
    grid.init_generators_full(gen_p, gen_v, gen_q, vreg_on, gen_min_q, gen_max_q, gen_bus);

    grid.add_gen_slackbus(0, 0.5);
    grid.add_gen_slackbus(1, 0.5);
    grid.tell_solver_need_reset();
    return grid;
}

// slack generator at bus 0 and an angle-droop hvdc line bus 1 -- bus 3. The Hvdc
// extension reads the FORWARD map only, the one that was already kept in sync, so
// this is the control that must have been right all along.
LSGrid make_droop_hvdc_grid()
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

    Eigen::VectorXi bus1(1), bus2(1);
    bus1 << 1;
    bus2 << 3;
    // ConverterStationContainer::ConverterType: VSC = 0
    // HvdcLineContainer::ConvertersMode: SIDE_1_RECTIFIER = 0
    const std::vector<int> type1{0}, type2{0}, mode{0};
    const std::vector<bool> vreg1{false}, vreg2{false}, droop_on{true};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm1_pu(1), vm2_pu(1), q_min(1), q_max(1);
    RealVect pf1(1), pf2(1), p_set(1), p_max(1), p0(1), k_deg(1);
    vm1_pu << 1.0;
    vm2_pu << 1.0;
    q_min << -1e6;
    q_max << 1e6;
    pf1 << 1.;
    pf2 << 1.;
    p_set << 0.;
    p_max << 300.;
    p0 << 10.;
    k_deg << 5.;
    grid.init_hvdc_lines(bus1, bus2, type1, type2, zero, zero,
                         vreg1, vreg2, vm1_pu, vm2_pu, zero, zero,
                         q_min, q_max, q_min, q_max, pf1, pf2, mode, p_set,
                         zero, zero, droop_on, p0, k_deg, p_max, p_max);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

CplxVect single_shot(LSGrid & grid)
{
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    const CplxVect V = grid.ac_pf(flat_start(grid), 30, 1e-10);
    REQUIRE(V.size() == NB_BUS);
    return V;
}

// one TimeSeries step carrying the grid's own (unchanged) injections
CplxVect one_step_timeseries(LSGrid & grid, int nb_gen)
{
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(grid);
    RealMat gen_p(1, nb_gen), sgen_p(1, 0), load_p(1, 1), load_q(1, 1);
    gen_p.setZero();
    load_p << LOAD_P;
    load_q << LOAD_Q;
    REQUIRE(ts.compute_Vs(gen_p, sgen_p, load_p, load_q, flat_start(grid), 30, 1e-10) == 1);
    return ts.get_voltages().row(0);
}

// one ContingencyAnalysis "contingency" that disconnects nothing: the base case,
// solved through the whole ContingencyAnalysis machinery
CplxVect empty_contingency(LSGrid & grid)
{
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ContingencyAnalysis ca(grid);
    ca.add_nk(std::vector<int>());
    ca.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(ca.get_voltages().rows() == 1);
    return ca.get_voltages().row(0);
}

void check_same_voltages(const CplxVect & batch, const CplxVect & ref)
{
    REQUIRE(batch.size() == ref.size());
    for (int b = 0; b < ref.size(); ++b) {
        INFO("bus " << b);
        CHECK(std::abs(batch(b)) == Approx(std::abs(ref(b))).epsilon(1e-8));
        CHECK(std::arg(batch(b)) == Approx(std::arg(ref(b))).margin(1e-9));
    }
}

}  // namespace


TEST_CASE("a remote-regulating generator holds its bus through TimeSeries", "[batch][vctrl]")
{
    LSGrid ref_grid = make_remote_gen_grid();
    const CplxVect ref = single_shot(ref_grid);
    // the point of the setup: the single shot really does regulate bus 3
    REQUIRE(std::abs(ref(3)) == Approx(V_SET).epsilon(1e-8));

    LSGrid grid = make_remote_gen_grid();
    check_same_voltages(one_step_timeseries(grid, 2), ref);
}

TEST_CASE("a remote-regulating generator holds its bus through ContingencyAnalysis", "[batch][vctrl]")
{
    LSGrid ref_grid = make_remote_gen_grid();
    const CplxVect ref = single_shot(ref_grid);
    REQUIRE(std::abs(ref(3)) == Approx(V_SET).epsilon(1e-8));

    LSGrid grid = make_remote_gen_grid();
    check_same_voltages(empty_contingency(grid), ref);
}

TEST_CASE("a batch on an already-solved grid regulates too", "[batch][vctrl]")
{
    // BaseBatchSolverSynch copies the grid, and LSGrid's copy constructor
    // reset()s the solver-side labelling -- so pre-solving the source grid does
    // NOT hand the batch a usable mapping. It has to be rebuilt by
    // pre_process_solver and written back, on every path.
    LSGrid ref_grid = make_remote_gen_grid();
    const CplxVect ref = single_shot(ref_grid);

    LSGrid grid = make_remote_gen_grid();
    single_shot(grid);   // solved BEFORE the batch copies it
    check_same_voltages(one_step_timeseries(grid, 2), ref);
}

TEST_CASE("two generators sharing one regulated bus share it in a batch too", "[batch][vctrl]")
{
    // Sharing is the supported multi-controller case (test_powerflow_algorithm.cpp
    // pins the Q split itself, single-shot). Here: the whole group has to survive the
    // batch path, which before the mapping sync got no controllers at all.
    LSGrid ref_grid = make_shared_reg_bus_grid();
    const CplxVect ref = single_shot(ref_grid);
    REQUIRE(std::abs(ref(3)) == Approx(V_SET).epsilon(1e-8));
    // and the sharing really is active: both controllers inject, in proportion to
    // their (qmax - qmin) keys, 2000 MVAr against 500 MVAr
    const RealVect q = ref_grid.get_algo().get_controller_q();
    const IntVect elem = ref_grid.get_algo().get_controller_elem_id();
    REQUIRE(q.size() == 2);
    real_type q1 = 0., q2 = 0.;
    for (int c = 0; c < q.size(); ++c) {
        if (elem(c) == 1) q1 = q(c);
        if (elem(c) == 2) q2 = q(c);
    }
    REQUIRE(std::abs(q1) > 1e-6);
    CHECK(q1 / 2000. == Approx(q2 / 500.).epsilon(1e-8));

    LSGrid ts_grid = make_shared_reg_bus_grid();
    check_same_voltages(one_step_timeseries(ts_grid, 3), ref);

    LSGrid ca_grid = make_shared_reg_bus_grid();
    check_same_voltages(empty_contingency(ca_grid), ref);
}

TEST_CASE("a voltage-mode SVC regulates through the batch algorithms", "[batch][vctrl][svc]")
{
    const real_type slope_pu = GENERATE(static_cast<real_type>(0.), static_cast<real_type>(0.02));
    LSGrid ref_grid = make_remote_svc_grid(slope_pu);
    const CplxVect ref = single_shot(ref_grid);
    // a flat SVC lands exactly on the setpoint, a sloped one strictly below it
    // (V = v_set - s.Q with Q > 0) -- either way it must be raised well above the
    // ~0.96 pu the same grid settles at with the SVC ignored
    REQUIRE(std::abs(ref(3)) > 1.);

    LSGrid ts_grid = make_remote_svc_grid(slope_pu);
    check_same_voltages(one_step_timeseries(ts_grid, 1), ref);

    LSGrid ca_grid = make_remote_svc_grid(slope_pu);
    check_same_voltages(empty_contingency(ca_grid), ref);
}

TEST_CASE("a free-Vm distributed-slack participant keeps its Vm unknown in a batch", "[batch][slack]")
{
    // Same root cause as the voltage-control tests, different extension: `Base`
    // reads slack_bus_id_ac_solver_ through get_free_vm_slack_solver_buses().
    LSGrid ref_grid = make_free_vm_slack_grid();
    const CplxVect ref = single_shot(ref_grid);
    // bus 1 is NOT pinned: it must have moved off the 1.0 pu flat start
    REQUIRE(std::abs(std::abs(ref(1)) - 1.0) > 1e-4);

    LSGrid ts_grid = make_free_vm_slack_grid();
    check_same_voltages(one_step_timeseries(ts_grid, 2), ref);

    LSGrid ca_grid = make_free_vm_slack_grid();
    check_same_voltages(empty_contingency(ca_grid), ref);
}

TEST_CASE("an hvdc angle droop survives the batch algorithms", "[batch][hvdc]")
{
    // The control: the Hvdc extension only ever needed the FORWARD map, which was
    // already synced, so this passed before the reverse map / slack sync too.
    LSGrid ref_grid = make_droop_hvdc_grid();
    const CplxVect ref = single_shot(ref_grid);

    LSGrid ts_grid = make_droop_hvdc_grid();
    check_same_voltages(one_step_timeseries(ts_grid, 1), ref);

    LSGrid ca_grid = make_droop_hvdc_grid();
    check_same_voltages(empty_contingency(ca_grid), ref);
}

TEST_CASE("a contingency stranding a regulated bus is reported, not silently wrong", "[batch][vctrl]")
{
    // Line 2 is bus2--bus3, and bus 3 is the REGULATED bus: removing it strands
    // the very bus the generator controls.
    //
    // In the default mode such a contingency is skipped outright (zero voltages).
    // In "handle disconnected grid" mode the stranded bus is masked -- its P/Q
    // rows are replaced by identity rows -- but masking is deliberately a
    // value-only change that must not touch the J pattern, so the VoltageControl
    // extension keeps its bordered row `Vm(3) - v_set = 0` against a bus whose Vm
    // the mask pins. That system has no solution and the contingency comes back as
    // a DIVERGENCE. Not a good answer, but an honest one: what must never happen
    // is a converged row computed with the regulation quietly dropped.
    for (bool handle_disconnected : {false, true}) {
        INFO("handle_disconnected_grid = " << handle_disconnected);
        LSGrid grid = make_remote_gen_grid();
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        ContingencyAnalysis ca(grid, /*compute_limit_violations=*/true);
        ca.set_handle_disconnected_grid(handle_disconnected);
        ca.add_n1(2);
        ca.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(ca.get_voltages().rows() == 1);
        REQUIRE(ca.converged().size() == 1);
        CHECK(ca.converged()[0] == 0);
        for (int b = 0; b < NB_BUS; ++b) CHECK(std::abs(ca.get_voltages()(0, b)) == 0.);
    }
}

TEST_CASE("stale solver state on the source grid cannot leak into a batch", "[batch][stale]")
{
    // The nastiest shape of the mapping bug this file is about is not an EMPTY
    // member but a stale NON-empty one, sized for a grid that no longer exists:
    //
    //   1. ac_pf runs                -> the AC maps are built for N buses
    //   2. the grid changes so that a bus leaves the solver
    //                                -> the AC maps are now wrong, and still N long
    //   3. dc_pf runs                -> the DC maps are rebuilt at N-1; AC stays stale
    //   4. a batch is built and runs  -> it inherits the grid's AC algorithm
    //                                   (LSGrid::change_algorithm routes a DC type to
    //                                   the separate _dc_algo slot, so asking the GRID
    //                                   for DC never switches a batch over)
    //
    // A stale map is worse than an empty one: an empty one reads as "nothing to do",
    // while a too-long one indexes past the end of the current solver vectors. What
    // makes this safe is that a batch never reads the source grid's solver state --
    // BaseBatchSolverSynch copies the grid (and LSGrid's copy constructor reset()s all
    // of it), then prepare_solver_input_base calls tell_all_changed() before
    // pre_process_solver, rebuilding every mapping from scratch on that copy.
    //
    // This case is deliberately cheap and allocation-light so it also runs under the
    // valgrind pass of the C++ suite and the ASan/UBSan + Eigen-assertion builds,
    // where an out-of-bounds read that stays inside the heap would still be caught.
    const int nb_bus = 5;
    auto build = [](LSGrid & grid) {
        grid.set_sn_mva(100.);
        grid.set_init_vm_pu(1.0);
        grid.init_bus(static_cast<unsigned int>(nb_bus), 1, RealVect::Constant(nb_bus, 138.), 0, 0);
        // lines 0:0-1, 1:1-2, 2:2-3, 3:1-4 -- bus 4 hangs off bus 1 by line 3 alone
        Eigen::VectorXi f(4), t(4);
        f << 0, 1, 2, 1;
        t << 1, 2, 3, 4;
        grid.init_powerlines(RealVect::Constant(4, 0.01), RealVect::Constant(4, 0.1),
                             CplxVect::Zero(4), f, t);
        RealVect lp(2), lq(2);
        lp << LOAD_P, 5.;
        lq << LOAD_Q, 1.;
        Eigen::VectorXi lb(2);
        lb << 3, 4;                      // load 1 is the only element on bus 4
        grid.init_loads(lp, lq, lb);

        RealVect gp(2), gv(2), gqn(2), gqx(2);
        Eigen::VectorXi gb(2);
        gp << 0., 0.;
        gv << 1.02, V_SET;
        gqn << -1000., -1000.;
        gqx << 1000., 1000.;
        gb << 0, 1;
        grid.init_generators(gp, gv, gqn, gqx, gb);
        grid.add_gen_slackbus(0, 1.);
        grid.set_gen_regulated_bus(1, 3);   // remote control: needs id_ac_solver_to_me_
        grid.tell_solver_need_reset();
    };
    // takes bus 4 out of the solver: open its only line and its only load
    auto strand_bus_4 = [](LSGrid & grid) {
        grid.deactivate_powerline(3);
        grid.deactivate_load(1);
        grid.deactivate_bus(GlobalBusId(4));
    };
    const CplxVect V0 = CplxVect::Constant(nb_bus, cplx_type(1., 0.));

    LSGrid grid;
    build(grid);
    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(grid.ac_pf(V0, 30, 1e-11).size() == nb_bus);
    // the AC maps now describe the FULL grid
    REQUIRE(static_cast<int>(grid.id_ac_solver_to_me().size()) == nb_bus);

    strand_bus_4(grid);
    grid.change_algorithm(AlgorithmType::DC_SparseLU);
    REQUIRE(grid.dc_pf(V0, 30, 1e-11).size() == nb_bus);
    // the DC maps followed the change; the AC ones are still the old, too-long ones.
    // If they were not, this test would be vacuous.
    REQUIRE(static_cast<int>(grid.id_dc_solver_to_me().size()) == nb_bus - 1);
    REQUIRE(static_cast<int>(grid.id_ac_solver_to_me().size()) == nb_bus);

    // the batch inherits the AC algorithm and must nonetheless reproduce a clean
    // single-shot ac_pf of the stranded grid
    TimeSeries ts(grid);
    RealMat gen_p(1, 2), sgen_p(1, 0), load_p(1, 2), load_q(1, 2);
    gen_p << 0., 0.;
    load_p << LOAD_P, 5.;
    load_q << LOAD_Q, 1.;
    REQUIRE(ts.compute_Vs(gen_p, sgen_p, load_p, load_q, V0, 30, 1e-11) == 1);
    const CplxVect Vts = ts.get_voltages().row(0);

    LSGrid ref;
    build(ref);
    strand_bus_4(ref);
    ref.change_algorithm(AlgorithmType::NR_SparseLU);
    const CplxVect Vref = ref.ac_pf(V0, 30, 1e-11);
    REQUIRE(Vref.size() == nb_bus);
    CHECK(std::abs(Vref(3)) == Approx(V_SET).epsilon(1e-8));  // control still applied

    for (int b = 0; b < nb_bus; ++b) {
        if (b == 4) continue;  // stranded: the batch reports 0, ac_pf reports init_vm_pu
        INFO("bus " << b);
        CHECK(std::abs(Vts(b)) == Approx(std::abs(Vref(b))).epsilon(1e-9));
        CHECK(std::arg(Vts(b)) == Approx(std::arg(Vref(b))).margin(1e-10));
    }
}
