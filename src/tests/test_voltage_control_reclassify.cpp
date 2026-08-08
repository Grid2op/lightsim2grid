// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Which mechanism sets a bus' voltage magnitude: the classical PV treatment, or a
// bordered VoltageControl group?
//
// The rule (LSGrid::get_group_controlled_buses, applied in LSGrid::fillpv_pq): a bus
// keeps the PV treatment as long as the only things regulating it stand ON it --
// generators, or voltage-regulating hvdc converter stations. As soon as an ACTIVE
// REMOTE regulator -- or any voltage-mode SVC, which is always a group controller --
// aims at the bus, the bus keeps its own Vm unknown and every regulator of it, the
// local ones included, becomes a member of the group. Generators and converter
// stations share on the same key (a reactive range in MVAr); an SVC's key is a
// susceptance range, which is why an SVC must still be alone in its group.
//
// The bug this pins: fillpv only knew how to keep a controller's OWN bus out of PV.
// A local regulator standing on somebody else's regulated bus still claimed it as
// PV, which left the group's voltage row with no Vm column to write into -- a
// structurally empty row, hence a singular Jacobian -- so the configuration had to
// be refused outright even though it is well posed (the local machine belongs in the
// group, and the sharing row then fixes the reactive split).

#include <cmath>
#include <set>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

const real_type V_SET = 1.03;
const real_type LOAD_P = 50.;
const real_type LOAD_Q = 10.;
const int NB_BUS = 4;

// 4-bus radial feeder 0-1-2-3, 50 MW / 10 MVAr at bus 3, sn_mva = 100
LSGrid make_skeleton()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    const RealVect vn = RealVect::Constant(NB_BUS, 138.);
    grid.init_bus(static_cast<unsigned int>(NB_BUS), 1, vn, 0, 0);
    const RealVect r = RealVect::Constant(NB_BUS - 1, 0.01);
    const RealVect x = RealVect::Constant(NB_BUS - 1, 0.1);
    const CplxVect h = CplxVect::Zero(NB_BUS - 1);
    Eigen::VectorXi f(NB_BUS - 1), t(NB_BUS - 1);
    for (int i = 0; i < NB_BUS - 1; ++i) { f(i) = i; t(i) = i + 1; }
    grid.init_powerlines(r, x, h, f, t);
    RealVect lp(1), lq(1); lp << LOAD_P; lq << LOAD_Q;
    Eigen::VectorXi lb(1); lb << NB_BUS - 1;
    grid.init_loads(lp, lq, lb);
    return grid;
}

void add_gens(LSGrid & grid,
              const std::vector<int> & buses,
              const std::vector<real_type> & v_pu,
              const std::vector<real_type> & q_range_mvar,
              const std::vector<bool> & vreg_on = std::vector<bool>())
{
    const int nb = static_cast<int>(buses.size());
    RealVect p(nb), v(nb), q(nb), qmin(nb), qmax(nb);
    Eigen::VectorXi b(nb);
    for (int i = 0; i < nb; ++i) {
        p(i) = 0.; v(i) = v_pu[i]; q(i) = 0.;
        qmin(i) = -0.5 * q_range_mvar[i];
        qmax(i) =  0.5 * q_range_mvar[i];
        b(i) = buses[i];
    }
    if (vreg_on.empty()) grid.init_generators(p, v, qmin, qmax, b);
    else                 grid.init_generators_full(p, v, q, vreg_on, qmin, qmax, b);
}

void add_voltage_svc(LSGrid & grid, int svc_bus, int reg_bus,
                     real_type target_vm_pu, real_type slope_pu)
{
    RealVect tv(1), qs(1), sl(1), bmin(1), bmax(1);
    tv << target_vm_pu; qs << 0.; sl << slope_pu; bmin << -100.; bmax << 100.;
    Eigen::VectorXi reg(1), bus(1);
    reg << reg_bus; bus << svc_bus;
    grid.init_svcs({1}, tv, qs, sl, bmin, bmax, reg, bus);   // RegulationMode::VOLTAGE
}

// One VSC-VSC hvdc line between `bus_1` and `bus_2`, fixed active setpoint (no
// droop). `vreg_1` / `vreg_2` say whether each station regulates the magnitude of
// its own bus, at `target_vm_pu`, with +/- q_range/2 MVAr of reactive capability.
void add_vreg_hvdc(LSGrid & grid, int bus_1, int bus_2,
                   bool vreg_1, bool vreg_2,
                   real_type target_vm_pu, real_type q_range_mvar,
                   real_type p_setpoint_mw = 0.)
{
    Eigen::VectorXi b1(1), b2(1);
    b1 << bus_1;
    b2 << bus_2;
    const std::vector<int> type1{0}, type2{0}, mode{0};   // VSC / VSC, side 1 rectifier
    const std::vector<bool> vr1{vreg_1}, vr2{vreg_2}, droop_off{false};
    const RealVect zero = RealVect::Zero(1);
    RealVect vm1(1), vm2(1), qmin(1), qmax(1), pf1(1), pf2(1), pset(1), pmax(1);
    vm1 << target_vm_pu;
    vm2 << target_vm_pu;
    qmin << -0.5 * q_range_mvar;
    qmax << 0.5 * q_range_mvar;
    pf1 << 1.;
    pf2 << 1.;
    pset << p_setpoint_mw;
    pmax << 300.;
    grid.init_hvdc_lines(b1, b2, type1, type2, zero, zero,
                         vr1, vr2, vm1, vm2, zero, zero,
                         qmin, qmax, qmin, qmax, pf1, pf2, mode, pset,
                         zero, zero, droop_off, zero, zero, pmax, pmax);
    grid.tell_solver_need_reset();
}

CplxVect flat(const LSGrid & g)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(g.total_bus()), cplx_type(1., 0.));
}

CplxVect solve(LSGrid & g)
{
    g.change_algorithm(AlgorithmType::NR_SparseLU);
    return g.ac_pf(flat(g), 30, 1e-10);
}

// per-controller reactive output (pu), keyed by generator id
real_type gen_ctrl_q(const LSGrid & g, int gen_id)
{
    const RealVect q = g.get_algo().get_controller_q();
    const IntVect kind = g.get_algo().get_controller_kind();
    const IntVect elem = g.get_algo().get_controller_elem_id();
    for (int c = 0; c < q.size(); ++c)
        if (kind(c) == ls2g::VoltageControlSolverData::GEN && elem(c) == gen_id) return q(c);
    return std::nan("");
}

// reactive output (pu) of the controller of the given kind and element id
real_type ctrl_q_of(const LSGrid & g, int want_kind, int want_elem)
{
    const RealVect q = g.get_algo().get_controller_q();
    const IntVect kind = g.get_algo().get_controller_kind();
    const IntVect elem = g.get_algo().get_controller_elem_id();
    for (int c = 0; c < q.size(); ++c)
        if (kind(c) == want_kind && elem(c) == want_elem) return q(c);
    return std::nan("");
}

}  // namespace


TEST_CASE("a locally regulated bus stays PV unless a group claims it", "[reclassify]")
{
    // the classification rule itself, on grid ids, before any solve
    SECTION("one local generator: PV, no group") {
        LSGrid g = make_skeleton();
        add_gens(g, {0, 1}, {1.02, V_SET}, {2000., 2000.});
        g.add_gen_slackbus(0, 1.);
        CHECK(g.get_group_controlled_buses().empty());
    }
    SECTION("SEVERAL local generators on one bus: still PV, still no group") {
        // the common case that must not change: two machines on bus 1, both
        // regulating it locally, share Q through the per-bus redistribution
        LSGrid g = make_skeleton();
        add_gens(g, {0, 1, 1}, {1.02, V_SET, V_SET}, {2000., 2000., 500.});
        g.add_gen_slackbus(0, 1.);
        CHECK(g.get_group_controlled_buses().empty());
    }
    SECTION("a remote regulator claims its target") {
        LSGrid g = make_skeleton();
        add_gens(g, {0, 1}, {1.02, V_SET}, {2000., 2000.});
        g.add_gen_slackbus(0, 1.);
        g.set_gen_regulated_bus(1, 3);
        const std::set<int> grp = g.get_group_controlled_buses();
        REQUIRE(grp.size() == 1);
        CHECK(*grp.begin() == 3);
    }
    SECTION("a voltage-mode SVC always claims its target") {
        LSGrid g = make_skeleton();
        add_gens(g, {0}, {1.02}, {2000.});
        g.add_gen_slackbus(0, 1.);
        add_voltage_svc(g, 2, 3, V_SET, 0.);
        const std::set<int> grp = g.get_group_controlled_buses();
        REQUIRE(grp.size() == 1);
        CHECK(*grp.begin() == 3);
    }
}

TEST_CASE("several local generators on one bus keep the classical PV treatment", "[reclassify]")
{
    // regression guard for the rule above: this is the ordinary configuration, it
    // must solve exactly as it did before any reclassification existed
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1, 1}, {1.02, V_SET, V_SET}, {2000., 2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.tell_solver_need_reset();
    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(1)) == Approx(V_SET).epsilon(1e-8));
    // no VoltageControl group was built: the PV path handled it
    CHECK(g.get_algo().get_controller_q().size() == 0);
}

TEST_CASE("a remote regulator may target a locally regulated bus", "[reclassify]")
{
    // THE BUG. gen 2 sits on bus 2 and regulates it locally; gen 1 (bus 1) aims at
    // bus 2 remotely with the same setpoint. Before the reclassification bus 2 was
    // taken as PV and the whole grid was refused.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1, 2}, {1.02, V_SET, V_SET}, {2000., 2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 2);
    g.tell_solver_need_reset();

    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(2)) == Approx(V_SET).epsilon(1e-8));

    // both machines are group members, and the sharing row split their reactive
    // output proportionally to (qmax - qmin): 2000 MVAr against 500 MVAr
    const real_type q_remote = gen_ctrl_q(g, 1);
    const real_type q_local  = gen_ctrl_q(g, 2);
    REQUIRE(std::isfinite(q_remote));
    REQUIRE(std::isfinite(q_local));
    REQUIRE(std::abs(q_remote) > 1e-6);
    CHECK(q_remote / 2000. == Approx(q_local / 500.).epsilon(1e-8));
}

TEST_CASE("a remote regulator disagreeing with the local one is rejected on setpoints", "[reclassify]")
{
    // same shape, conflicting setpoints: now that both are in one group this is
    // caught by the group's own setpoint check, which names the real problem
    // instead of complaining about the bus classification
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1, 2}, {1.02, V_SET + 0.02, V_SET}, {2000., 2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 2);
    g.tell_solver_need_reset();

    REQUIRE_THROWS_WITH(solve(g), Catch::Matchers::ContainsSubstring("conflicting voltage setpoints"));
}

TEST_CASE("a remote regulator may target a locally pinned slack bus", "[reclassify]")
{
    // the slack's own generator pins it, so before the reclassification the slack
    // had no Vm unknown and a remote controller aiming at it was refused. Now the
    // slack generator joins the group and MultiSlack grants the free Vm.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1}, {V_SET, V_SET}, {2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 0);   // -> the slack bus
    g.tell_solver_need_reset();

    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(0)) == Approx(V_SET).epsilon(1e-8));
    const real_type q_remote = gen_ctrl_q(g, 1);
    const real_type q_slack  = gen_ctrl_q(g, 0);
    REQUIRE(std::isfinite(q_remote));
    REQUIRE(std::isfinite(q_slack));
    CHECK(q_remote / 500. == Approx(q_slack / 2000.).epsilon(1e-8));
}

TEST_CASE("a remote regulator may target an unpinned slack bus", "[reclassify]")
{
    // the slack generator has voltage_regulator_on == false, so nothing pins the
    // slack: MultiSlack gives it a free Vm + Q equation and the remote controller
    // can hold it. Previously refused because a slack bus is never in `bus_pq_`,
    // the regulated-bus check having no `has_free_q` escape hatch.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1}, {1.02, V_SET}, {2000., 2000.}, {false, true});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 0);
    g.tell_solver_need_reset();

    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(0)) == Approx(V_SET).epsilon(1e-8));
}

TEST_CASE("a remote regulator on a bus pinned by another machine is still rejected", "[reclassify]")
{
    // The configuration that genuinely has NO solution, and must stay refused: gen 2
    // sits on bus 1 and pins it locally; gen 1, ALSO on bus 1, aims remotely at bus
    // 3. Whatever reactive gen 1 injects at bus 1 is absorbed one-for-one by gen 2
    // holding Vm(bus 1) with free Q, so gen 1 has no influence whatsoever on bus 3.
    // Nothing about bus classification can fix that.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1, 1}, {1.02, V_SET, V_SET}, {2000., 2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 3);
    g.tell_solver_need_reset();

    REQUIRE_THROWS_WITH(solve(g),
                        Catch::Matchers::ContainsSubstring("OWN bus has no reactive (Q) equation"));
}

TEST_CASE("an SVC-regulated bus is taken off the PV path too", "[reclassify][svc]")
{
    // An SVC is always a group controller, so its regulated bus must leave the PV
    // path even when a local generator sits on it. Sharing an SVC with another
    // controller is a separate v1 restriction (the sharing key of an SVC is a
    // susceptance range, not a reactive range), so the outcome is still a refusal --
    // but now it is the SVC-sharing restriction that says so, not a bus-class error.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 2}, {1.02, V_SET}, {2000., 500.});
    g.add_gen_slackbus(0, 1.);
    add_voltage_svc(g, 1, 2, V_SET, 0.);   // SVC at bus 1 regulating bus 2
    g.tell_solver_need_reset();

    CHECK(g.get_group_controlled_buses().count(2) == 1);
    REQUIRE_THROWS_WITH(solve(g),
                        Catch::Matchers::ContainsSubstring("shares a regulated bus with other controllers"));
}

TEST_CASE("an SVC alone on a locally-free bus still works", "[reclassify][svc]")
{
    // the control: no generator on the regulated bus, so the SVC is alone in its
    // group and the reclassification changes nothing
    LSGrid g = make_skeleton();
    add_gens(g, {0}, {1.02}, {2000.});
    g.add_gen_slackbus(0, 1.);
    add_voltage_svc(g, 2, 3, V_SET, 0.);
    g.tell_solver_need_reset();

    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(3)) == Approx(V_SET).epsilon(1e-8));
}

TEST_CASE("a voltage-regulating hvdc station keeps the PV path when nothing else claims its bus",
          "[reclassify][hvdc]")
{
    // the control: an hvdc station regulating its own bus, with no group anywhere.
    // It pins the bus exactly as before -- no controller is built for it.
    LSGrid g = make_skeleton();
    add_gens(g, {0}, {1.02}, {2000.});
    g.add_gen_slackbus(0, 1.);
    add_vreg_hvdc(g, /*bus_1=*/1, /*bus_2=*/3, /*vreg_1=*/true, /*vreg_2=*/false,
                  V_SET, /*q_range_mvar=*/600.);

    CHECK(g.get_group_controlled_buses().empty());
    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(1)) == Approx(V_SET).epsilon(1e-8));
    CHECK(g.get_algo().get_controller_q().size() == 0);  // PV path, no group
}

TEST_CASE("a voltage-regulating hvdc station joins the group that claims its bus",
          "[reclassify][hvdc]")
{
    // hvdc side 1 sits on bus 1 and regulates it; gen 1 (bus 2) aims at bus 1
    // remotely with the same setpoint. Bus 1 leaves the PV path, and the station is
    // enrolled alongside the generator instead of having its regulation dropped.
    LSGrid g = make_skeleton();
    add_gens(g, {0, 2}, {1.02, V_SET}, {2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 1);          // gen 1 -> bus 1, where the station is
    add_vreg_hvdc(g, 1, 3, true, false, V_SET, /*q_range_mvar=*/2000.);

    CHECK(g.get_group_controlled_buses().count(1) == 1);
    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(1)) == Approx(V_SET).epsilon(1e-8));

    // two controllers in the group: the generator and the side-1 station, sharing
    // reactive power in proportion to their MVAr ranges (500 against 2000)
    REQUIRE(g.get_algo().get_controller_q().size() == 2);
    const real_type q_gen = ctrl_q_of(g, ls2g::VoltageControlSolverData::GEN, 1);
    const real_type q_sta = ctrl_q_of(g, ls2g::VoltageControlSolverData::HVDC_SIDE_1, 0);
    REQUIRE(std::isfinite(q_gen));
    REQUIRE(std::isfinite(q_sta));
    REQUIRE(std::abs(q_gen) > 1e-6);
    CHECK(q_gen / 500. == Approx(q_sta / 2000.).epsilon(1e-8));

    // and the write-back reached the station's own result, in MVAr
    const real_type res_q_mvar = g.get_dclines()[0].station_side_1.res_q_mvar;
    CHECK(res_q_mvar == Approx(q_sta * 100.).epsilon(1e-8));   // sn_mva == 100
}

TEST_CASE("the far hvdc station may be the one in the group", "[reclassify][hvdc]")
{
    // same, on side 2, so the HVDC_SIDE_2 kind and its write-back are exercised too
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1}, {1.02, V_SET}, {2000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(1, 2);          // gen 1 -> bus 2, the side-2 station bus
    add_vreg_hvdc(g, 1, 2, false, true, V_SET, /*q_range_mvar=*/2000.);

    CHECK(g.get_group_controlled_buses().count(2) == 1);
    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(2)) == Approx(V_SET).epsilon(1e-8));

    REQUIRE(g.get_algo().get_controller_q().size() == 2);
    const real_type q_sta = ctrl_q_of(g, ls2g::VoltageControlSolverData::HVDC_SIDE_2, 0);
    REQUIRE(std::isfinite(q_sta));
    CHECK(g.get_dclines()[0].station_side_2.res_q_mvar == Approx(q_sta * 100.).epsilon(1e-8));
}

TEST_CASE("a generator, an hvdc station and one regulated bus share it three ways",
          "[reclassify][hvdc]")
{
    // a local generator AND an hvdc station on bus 1, plus a remote generator aiming
    // at bus 1: one group of three, all keyed on MVAr ranges
    LSGrid g = make_skeleton();
    add_gens(g, {0, 1, 2}, {1.02, V_SET, V_SET}, {2000., 1000., 500.});
    g.add_gen_slackbus(0, 1.);
    g.set_gen_regulated_bus(2, 1);          // gen 2 (bus 2) -> bus 1
    add_vreg_hvdc(g, 1, 3, true, false, V_SET, /*q_range_mvar=*/2000.);

    const CplxVect V = solve(g);
    REQUIRE(V.size() == NB_BUS);
    CHECK(std::abs(V(1)) == Approx(V_SET).epsilon(1e-8));

    REQUIRE(g.get_algo().get_controller_q().size() == 3);
    const real_type q_local  = ctrl_q_of(g, ls2g::VoltageControlSolverData::GEN, 1);
    const real_type q_remote = ctrl_q_of(g, ls2g::VoltageControlSolverData::GEN, 2);
    const real_type q_sta    = ctrl_q_of(g, ls2g::VoltageControlSolverData::HVDC_SIDE_1, 0);
    REQUIRE(std::isfinite(q_local));
    REQUIRE(std::isfinite(q_remote));
    REQUIRE(std::isfinite(q_sta));
    REQUIRE(std::abs(q_remote) > 1e-6);
    CHECK(q_local / 1000. == Approx(q_remote / 500.).epsilon(1e-8));
    CHECK(q_sta / 2000.   == Approx(q_remote / 500.).epsilon(1e-8));
}
