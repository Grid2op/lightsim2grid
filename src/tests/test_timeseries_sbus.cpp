// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// TimeSeries::compute_Vs does NOT use the injection vector the gridmodel builds. It
// rebuilds its own, per step, out of the four matrices the caller supplies:
// generator and static-generator ACTIVE power, and load active + reactive power.
//
// LSGrid::fillSbus_me stamps rather more than that -- storage units,
// REACTIVE_POWER-mode SVCs, the reactive setpoint of a generator whose voltage
// regulation is off, a static generator's reactive setpoint, the hvdc injections and
// (dc only) the phase-shifter term. None of those is a per-step input, so all of them
// are constant across the steps and have to be taken from the grid. Only the hvdc
// term ever was, so everything else used to be silently dropped: the powerflow
// converged on an injection the caller never asked for, off by as much as 0.08 pu of
// voltage on the four-bus feeder below.
//
// Every test here solves the same grid twice -- once single-shot, once as a one-step
// TimeSeries carrying the grid's own injections -- and requires the voltages to
// agree. That comparison is the real invariant: with the per-step matrices set to the
// grid's own target values, a TimeSeries step must reproduce ac_pf / dc_pf exactly.

#include <cmath>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::TimeSeries;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const real_type LOAD_P = 50.;
const real_type LOAD_Q = 10.;
const int NB_BUS = 4;

// 4-bus radial feeder 0-1-2-3 (r = 0.01, x = 0.1 pu), 50 MW / 10 MVAr at bus 3,
// sn_mva = 100 -- deliberately != 1 so the per-unit conversion is exercised too.
LSGrid make_skeleton()
{
    LSGrid g;
    g.set_sn_mva(100.);
    g.set_init_vm_pu(1.0);
    const RealVect vn = RealVect::Constant(NB_BUS, 138.);
    g.init_bus(static_cast<unsigned int>(NB_BUS), 1, vn, 0, 0);
    const RealVect r = RealVect::Constant(NB_BUS - 1, 0.01);
    const RealVect x = RealVect::Constant(NB_BUS - 1, 0.1);
    const CplxVect h = CplxVect::Zero(NB_BUS - 1);
    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for (int i = 0; i < NB_BUS - 1; ++i) { f(i) = i; t(i) = i + 1; }
    g.init_powerlines(r, x, h, f, t);
    RealVect lp(1), lq(1); lp << LOAD_P; lq << LOAD_Q;
    Eigen::VectorXi lb(1); lb << NB_BUS - 1;
    g.init_loads(lp, lq, lb);
    return g;
}

// a slack generator at bus 0 and, optionally, a second generator at bus 1 with
// voltage_regulator_on == false -- which is what makes its target_q part of Sbus
void add_gens(LSGrid & g, bool with_q_mode_gen)
{
    const int nb = with_q_mode_gen ? 2 : 1;
    RealVect p(nb), v(nb), q(nb), qmin(nb), qmax(nb);
    Eigen::VectorXi b(nb);
    p(0) = 0.; v(0) = 1.02; q(0) = 0.; qmin(0) = -1000.; qmax(0) = 1000.; b(0) = 0;
    if (with_q_mode_gen) {
        p(1) = 20.; v(1) = 1.0; q(1) = 30.;
        qmin(1) = -1000.; qmax(1) = 1000.; b(1) = 1;
    }
    const std::vector<bool> vreg = with_q_mode_gen ? std::vector<bool>{true, false}
                                                   : std::vector<bool>{true};
    g.init_generators_full(p, v, q, vreg, qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
}

// a slack generator at bus 0 plus a genuine voltage-regulating (PV) generator at bus 2,
// whose target vm_pu is `pv_target_vm` -- used by the modify_gen_v tests below, which
// need a real PV bus (unlike add_gens' optional second generator, which is Q-mode /
// voltage_regulator_on == false, ie NOT a PV bus).
void add_slack_and_pv_gens(LSGrid & g, real_type pv_target_vm)
{
    RealVect p(2), v(2), q(2), qmin(2), qmax(2);
    Eigen::VectorXi b(2);
    p(0) = 0.; v(0) = 1.02; q(0) = 0.; qmin(0) = -1000.; qmax(0) = 1000.; b(0) = 0;
    p(1) = 20.; v(1) = pv_target_vm; q(1) = 0.; qmin(1) = -1000.; qmax(1) = 1000.; b(1) = 2;
    const std::vector<bool> vreg{true, true};
    g.init_generators_full(p, v, q, vreg, qmin, qmax, b);
    g.add_gen_slackbus(0, 1.);
}

// a triangular (2-connected) 3-bus mesh -- bus 0 (slack), bus 1 (load), bus 2 (PV
// generator) -- so that disconnecting any single line still leaves every bus
// connected. Used by the ScenarioSweep modify_gen_v test below, which needs
// handle_disconnected_grid's masked path (_run_range_masked) to actually run without
// islanding the PV bus away from the slack.
LSGrid make_mesh_skeleton()
{
    LSGrid g;
    g.set_sn_mva(100.);
    g.set_init_vm_pu(1.0);
    const RealVect vn = RealVect::Constant(3, 138.);
    g.init_bus(3u, 1u, vn, 0, 0);
    const RealVect r = RealVect::Constant(3, 0.01);
    const RealVect x = RealVect::Constant(3, 0.1);
    const CplxVect h = CplxVect::Zero(3);
    Eigen::VectorXi f(3), t(3);
    f(0) = 0; t(0) = 1;  // bus0 - bus1
    f(1) = 1; t(1) = 2;  // bus1 - bus2
    f(2) = 0; t(2) = 2;  // bus0 - bus2 (the redundant edge)
    g.init_powerlines(r, x, h, f, t);
    RealVect lp(1), lq(1); lp << LOAD_P; lq << LOAD_Q;
    Eigen::VectorXi lb(1); lb << 1;
    g.init_loads(lp, lq, lb);
    add_slack_and_pv_gens(g, 1.02);
    return g;
}

void add_sgen(LSGrid & g, int bus, real_type p_mw, real_type q_mvar)
{
    RealVect sp(1), sq(1), pmin(1), pmax(1), qmin(1), qmax(1);
    sp << p_mw; sq << q_mvar;
    pmin << -1000.; pmax << 1000.; qmin << -1000.; qmax << 1000.;
    Eigen::VectorXi sb(1); sb << bus;
    g.init_sgens(sp, sq, pmin, pmax, qmin, qmax, sb);
}

void add_storage(LSGrid & g, int bus, real_type p_mw, real_type q_mvar)
{
    RealVect sp(1), sq(1);
    sp << p_mw; sq << q_mvar;
    Eigen::VectorXi sb(1); sb << bus;
    g.init_storages(sp, sq, sb);
}

void add_q_mode_svc(LSGrid & g, int bus, real_type q_set_mvar)
{
    RealVect tv(1), qs(1), sl(1), bmin(1), bmax(1);
    tv << 1.0; qs << q_set_mvar; sl << 0.; bmin << -1000.; bmax << 1000.;
    Eigen::VectorXi reg(1), b(1); reg << bus; b << bus;
    // SvcContainer::RegulationMode: REACTIVE_POWER = 2 (the mode that stamps Sbus;
    // VOLTAGE mode is solved by the VoltageControl extension instead)
    g.init_svcs({2}, tv, qs, sl, bmin, bmax, reg, b);
}

// a phase-shifting transformer in parallel with line 1-2. In DC its shift enters the
// injection through TrafoContainer::hack_Sbus_for_dc_phase_shifter, applied by
// fillSbus_me AFTER the per-unit division.
void add_phase_shifter(LSGrid & g, int bus1, int bus2, real_type shift_deg)
{
    RealVect tr(1), tx(1), ratio(1), shift(1);
    CplxVect tb(1);
    tr << 0.01; tx << 0.1; tb << cplx_type(0., 0.); ratio << 1.; shift << shift_deg;
    const std::vector<bool> tap_hv{true};
    Eigen::VectorXi b1(1), b2(1); b1 << bus1; b2 << bus2;
    g.init_trafo(tr, tx, tb, ratio, shift, tap_hv, b1, b2, false);
}

CplxVect flat(const LSGrid & g)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(g.total_bus()), cplx_type(1., 0.));
}

// Solve `ref` single-shot and `ts_grid` (an identical grid) as a one-step TimeSeries
// fed the grid's OWN injections, then require the two to agree. `ac` picks ac_pf/dc_pf
// and, through the grid's algorithm, the batch path used.
void check_timeseries_matches_single_shot(LSGrid & ref, LSGrid & ts_grid, bool ac = true)
{
    const AlgorithmType algo = ac ? AlgorithmType::NR_SparseLU : AlgorithmType::DC_SparseLU;
    ref.change_algorithm(algo);
    const CplxVect V = ac ? ref.ac_pf(flat(ref), 30, 1e-11) : ref.dc_pf(flat(ref), 30, 1e-11);
    REQUIRE(V.size() == NB_BUS);

    ts_grid.change_algorithm(algo);
    TimeSeries ts(ts_grid);
    // The batch classes inherit the grid's *AC* algorithm (BaseBatchSolverSynch's
    // constructor reads get_algo(), and LSGrid::change_algorithm routes a DC type to
    // the separate _dc_algo slot), so asking the GRID for DC is not enough -- the
    // batch has to be switched over explicitly, or it silently runs AC.
    if (!ac) ts.change_algorithm(algo);
    const int nb_gen = static_cast<int>(ref.get_generators_as_data().nb());
    const int nb_sgen = static_cast<int>(ref.get_static_generators_as_data().nb());
    const int nb_load = static_cast<int>(ref.get_loads_as_data().nb());
    RealMat gen_p(1, nb_gen), sgen_p(1, nb_sgen), load_p(1, nb_load), load_q(1, nb_load);
    gen_p.row(0) = ref.get_gen_target_p().transpose();
    if (nb_sgen > 0) sgen_p.row(0) = ref.get_sgen_target_p().transpose();
    load_p.row(0) = ref.get_load_target_p().transpose();
    load_q.row(0) = ref.get_loads_as_data().get_target_q().transpose();

    REQUIRE(ts.compute_Vs(gen_p, sgen_p, load_p, load_q, flat(ts_grid), 30, 1e-11) == 1);
    const CplxVect Vts = ts.get_voltages().row(0);
    REQUIRE(Vts.size() == NB_BUS);
    for (int b = 0; b < NB_BUS; ++b) {
        INFO("bus " << b << (ac ? " (ac)" : " (dc)"));
        CHECK(std::abs(Vts(b)) == Approx(std::abs(V(b))).epsilon(1e-9));
        CHECK(std::arg(Vts(b)) == Approx(std::arg(V(b))).margin(1e-10));
    }
}

}  // namespace


TEST_CASE("a TimeSeries step reproduces the single shot with only loads and generators",
          "[tssbus]")
{
    // the control: everything in this grid IS one of the four per-step matrices, so
    // it passed even before the constant remainder was taken from the grid
    LSGrid a = make_skeleton(); add_gens(a, false); a.tell_solver_need_reset();
    LSGrid b = make_skeleton(); add_gens(b, false); b.tell_solver_need_reset();
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a TimeSeries step keeps the reactive setpoint of a Q-mode generator", "[tssbus]")
{
    // voltage_regulator_on == false, so GeneratorContainer::fillSbus stamps
    // `target_q` -- which no per-step matrix carries
    LSGrid a = make_skeleton(); add_gens(a, true); a.tell_solver_need_reset();
    LSGrid b = make_skeleton(); add_gens(b, true); b.tell_solver_need_reset();
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a TimeSeries step keeps a static generator's reactive setpoint", "[tssbus]")
{
    // sgen_p carries the active power; the reactive setpoint has no matrix at all
    auto build = [](LSGrid & g) {
        add_gens(g, false);
        add_sgen(g, 2, 10., 25.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a TimeSeries step keeps storage units", "[tssbus]")
{
    // storages have no per-step matrix whatsoever: active AND reactive were dropped
    auto build = [](LSGrid & g) {
        add_gens(g, false);
        add_storage(g, 2, 15., 20.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a TimeSeries step keeps a REACTIVE_POWER-mode SVC", "[tssbus]")
{
    // the SVC mode that DOES stamp Sbus (VOLTAGE mode is a solver unknown instead,
    // and is legitimately absent from it)
    auto build = [](LSGrid & g) {
        add_gens(g, false);
        add_q_mode_svc(g, 2, 40.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a TimeSeries step keeps all of them at once", "[tssbus]")
{
    // all four together, so a fix that happened to handle them only one at a time
    // (or double-counted one) would show up here
    auto build = [](LSGrid & g) {
        add_gens(g, true);
        add_sgen(g, 2, 10., 25.);
        add_storage(g, 2, 15., 20.);
        add_q_mode_svc(g, 1, 40.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b);
}

TEST_CASE("a DC TimeSeries step reproduces dc_pf", "[tssbus][dc]")
{
    // DC consumes the real part of the very same per-step matrix
    // (BaseBatchSolverSynch::compute_one_powerflow falls back to Pbus_ only when no
    // Sbus is supplied, which is not the TimeSeries case), so a storage unit's active
    // power was dropped in DC too.
    auto build = [](LSGrid & g) {
        add_gens(g, false);
        add_storage(g, 2, 15., 20.);
        add_sgen(g, 1, 10., 25.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b, /*ac=*/false);
}

TEST_CASE("a DC TimeSeries step keeps the phase-shifter injection", "[tssbus][dc]")
{
    // TrafoContainer::hack_Sbus_for_dc_phase_shifter is applied by fillSbus_me AFTER
    // the per-unit division, i.e. deliberately in physical units -- so it also
    // exercises that the derived remainder does not re-scale what it must not.
    auto build = [](LSGrid & g) {
        add_gens(g, false);
        add_phase_shifter(g, 1, 2, 5.);
        g.tell_solver_need_reset();
    };
    LSGrid a = make_skeleton(); build(a);
    LSGrid b = make_skeleton(); build(b);
    check_timeseries_matches_single_shot(a, b, /*ac=*/false);
}

TEST_CASE("several TimeSeries steps still track the single shot", "[tssbus]")
{
    // the remainder is constant across steps, so it must be added to EVERY row, not
    // just the first -- and the varying matrices must still vary
    auto build = [](LSGrid & g) {
        add_gens(g, true);
        add_storage(g, 2, 15., 20.);
        add_q_mode_svc(g, 1, 40.);
        g.tell_solver_need_reset();
    };
    LSGrid ts_grid = make_skeleton(); build(ts_grid);
    ts_grid.change_algorithm(AlgorithmType::NR_SparseLU);

    const std::vector<real_type> scales{1.0, 1.15};
    TimeSeries ts(ts_grid);
    const int nb_steps = static_cast<int>(scales.size());
    RealMat gen_p(nb_steps, 2), sgen_p(nb_steps, 0), load_p(nb_steps, 1), load_q(nb_steps, 1);
    for (int s = 0; s < nb_steps; ++s) {
        gen_p(s, 0) = 0.;
        gen_p(s, 1) = 20. * scales[s];
        load_p(s, 0) = LOAD_P * scales[s];
        load_q(s, 0) = LOAD_Q * scales[s];
    }
    REQUIRE(ts.compute_Vs(gen_p, sgen_p, load_p, load_q, flat(ts_grid), 30, 1e-11) == 1);

    for (int s = 0; s < nb_steps; ++s) {
        INFO("step " << s);
        LSGrid ref = make_skeleton(); build(ref);
        ref.change_algorithm(AlgorithmType::NR_SparseLU);
        ref.change_p_gen(1, gen_p(s, 1));
        ref.change_p_load(0, load_p(s, 0));
        ref.change_q_load(0, load_q(s, 0));
        const CplxVect V = ref.ac_pf(flat(ref), 30, 1e-11);
        REQUIRE(V.size() == NB_BUS);
        const CplxVect Vts = ts.get_voltages().row(s);
        for (int b = 0; b < NB_BUS; ++b) {
            INFO("bus " << b);
            CHECK(std::abs(Vts(b)) == Approx(std::abs(V(b))).epsilon(1e-9));
        }
    }
}

// --- modify_gen_v: per-step generator voltage-magnitude target ------------------

TEST_CASE("a TimeSeries step with an unset modify_gen_v keeps the grid's own PV target",
          "[tssbus][genv]")
{
    // non-regression: modify_gen_v is deliberately never called here, so the PV
    // bus's converged |V| must still match the grid's own target_vm_pu, exactly as
    // before this setter existed (BaseBatchSolverSynch::_finish_preprocessing's
    // once-only set_vm() call is untouched).
    LSGrid ref = make_skeleton(); add_slack_and_pv_gens(ref, 1.05);
    ref.change_algorithm(AlgorithmType::NR_SparseLU);
    const CplxVect V = ref.ac_pf(flat(ref), 30, 1e-11);
    REQUIRE(V.size() == NB_BUS);

    LSGrid ts_grid = make_skeleton(); add_slack_and_pv_gens(ts_grid, 1.05);
    ts_grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(ts_grid);
    RealMat gen_p(1, 2), load_p(1, 1), load_q(1, 1);
    gen_p(0, 0) = 0.; gen_p(0, 1) = 20.;
    load_p(0, 0) = LOAD_P; load_q(0, 0) = LOAD_Q;
    ts.modify_gen_p(gen_p);
    ts.modify_load_p(load_p);
    ts.modify_load_q(load_q);
    ts.compute(flat(ts_grid), 30, 1e-11);
    REQUIRE(ts.get_status() == 1);
    const CplxVect Vts = ts.get_voltages().row(0);
    for (int b = 0; b < NB_BUS; ++b) {
        INFO("bus " << b);
        CHECK(std::abs(Vts(b)) == Approx(std::abs(V(b))).epsilon(1e-9));
    }
}

TEST_CASE("a TimeSeries step honors a per-row modify_gen_v target", "[tssbus][genv]")
{
    LSGrid ts_grid = make_skeleton(); add_slack_and_pv_gens(ts_grid, 1.02);
    ts_grid.change_algorithm(AlgorithmType::NR_SparseLU);
    TimeSeries ts(ts_grid);

    const std::vector<real_type> targets{1.00, 1.03, 0.97};
    const int nb_steps = static_cast<int>(targets.size());
    RealMat gen_p(nb_steps, 2), load_p(nb_steps, 1), load_q(nb_steps, 1), gen_v(nb_steps, 2);
    for (int s = 0; s < nb_steps; ++s) {
        gen_p(s, 0) = 0.; gen_p(s, 1) = 20.;
        load_p(s, 0) = LOAD_P; load_q(s, 0) = LOAD_Q;
        gen_v(s, 0) = 1.02;         // slack: kept constant across rows
        gen_v(s, 1) = targets[s];  // the PV generator at bus 2, varied per row
    }
    ts.modify_gen_p(gen_p);
    ts.modify_load_p(load_p);
    ts.modify_load_q(load_q);
    ts.modify_gen_v(gen_v);
    ts.compute(flat(ts_grid), 30, 1e-11);
    REQUIRE(ts.get_status() == 1);
    for (int s = 0; s < nb_steps; ++s) {
        INFO("step " << s);
        CHECK(std::abs(ts.get_voltages().row(s)(2)) == Approx(targets[s]).epsilon(1e-8));
    }
}

TEST_CASE("a ScenarioSweep row honors a per-row modify_gen_v target through the masked path",
          "[tssbus][genv]")
{
    // exercises _run_range_masked (BaseBatchSweep.hpp), the ONLY per-step solve path
    // that does not go through _run_one_step -- handle_disconnected_grid=true routes
    // every row through it regardless of whether that row's contingency actually
    // strands anything. make_mesh_skeleton is 2-connected so no row ever islands the
    // PV bus away from the slack; this isolates the modify_gen_v call site itself.
    LSGrid ts_grid = make_mesh_skeleton();
    ts_grid.change_algorithm(AlgorithmType::NR_SparseLU);
    ScenarioSweep sw(ts_grid);
    sw.set_handle_disconnected_grid(true);

    const std::vector<real_type> targets{1.00, 1.03, 0.97};
    const int nb_steps = static_cast<int>(targets.size());
    RealMat gen_p(nb_steps, 2), load_p(nb_steps, 1), load_q(nb_steps, 1), gen_v(nb_steps, 2);
    ScenarioSweep::BoolMat line_mask = ScenarioSweep::BoolMat::Zero(nb_steps, 3);
    for (int s = 0; s < nb_steps; ++s) {
        gen_p(s, 0) = 0.; gen_p(s, 1) = 20.;
        load_p(s, 0) = LOAD_P; load_q(s, 0) = LOAD_Q;
        gen_v(s, 0) = 1.02;
        gen_v(s, 1) = targets[s];
        line_mask(s, s % 3) = true;  // disconnect a different (redundant) line each row
    }
    sw.modify_gen_p(gen_p);
    sw.modify_load_p(load_p);
    sw.modify_load_q(load_q);
    sw.modify_gen_v(gen_v);
    sw.set_contingency_lines(line_mask);
    sw.compute(flat(ts_grid), 30, 1e-11);
    REQUIRE(sw.get_status() == 1);
    for (int s = 0; s < nb_steps; ++s) {
        INFO("step " << s);
        CHECK(std::abs(sw.get_voltages().row(s)(2)) == Approx(targets[s]).epsilon(1e-8));
    }
}
