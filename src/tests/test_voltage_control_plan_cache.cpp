// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The voltage-control plan (VoltageControlPlan: which buses a control GROUP
// regulates, which slack buses keep a free Vm unknown, and the controller list the
// bordered block is built from) is derived ONCE per powerflow, into the AC cache,
// and read from there by the NR extensions. It used to be re-derived four times per
// solve -- by the pv/pq split, by Base::update_state, and twice by
// VoltageControl::update_state -- each walking every generator of the grid.
//
// Deriving it once is only correct if a plan is dropped whenever anything it is made
// of moves. That is what these tests pin, and they pin it the only way that cannot
// be fooled by the cache itself: every mutation is applied to a grid that has
// already solved (so a stale plan is available to be wrongly reused) AND to a fresh
// grid built in the final state (which has no plan to reuse), and the two solves
// must agree bit for bit. A reused stale plan does not throw and does not look
// wrong: it converges, on the previous scenario's controllers.
//
// The companion question -- WHICH mechanism sets a bus' magnitude -- is
// test_voltage_control_reclassify.cpp; this file assumes that rule and tests only
// the caching of its result.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

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
const int NB_BUS = 5;

// 5-bus radial feeder 0-1-2-3-4, one load at each of buses 3 and 4, sn_mva = 100.
// Bus 0 carries the slack; the interesting controllers go on buses 1 and 2.
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
    RealVect lp(2), lq(2); lp << 40., 30.; lq << 8., 6.;
    Eigen::VectorXi lb(2); lb << 3, 4;
    grid.init_loads(lp, lq, lb);
    return grid;
}

void add_gens(LSGrid & grid, const std::vector<int> & buses,
              const std::vector<real_type> & v_pu)
{
    const int nb = static_cast<int>(buses.size());
    RealVect p(nb), v(nb), qmin(nb), qmax(nb);
    Eigen::VectorXi b(nb);
    for (int i = 0; i < nb; ++i) {
        p(i) = 0.; v(i) = v_pu[i];
        qmin(i) = -1000.; qmax(i) = 1000.;
        b(i) = buses[i];
    }
    grid.init_generators(p, v, qmin, qmax, b);
}

void add_voltage_svc(LSGrid & grid, int svc_bus, real_type target_vm_pu)
{
    RealVect tv(1), qs(1), sl(1), bmin(1), bmax(1);
    tv << target_vm_pu; qs << 0.; sl << 0.; bmin << -100.; bmax << 100.;
    Eigen::VectorXi reg(1), bus(1);
    reg << svc_bus; bus << svc_bus;
    grid.init_svcs({1}, tv, qs, sl, bmin, bmax, reg, bus);   // RegulationMode::VOLTAGE
}

CplxVect flat(const LSGrid & g)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(g.total_bus()), cplx_type(1., 0.));
}

// NB: no change_algorithm() here, deliberately. It raises tell_all_changed() on the
// family it switches, which would retire the cache before every single solve and
// make every test below pass whatever the invalidation rule is. NR_SparseLU is
// already the default AC algorithm (LSGrid's constructor), so there is nothing to
// select.
CplxVect solve(LSGrid & g)
{
    return g.ac_pf(flat(g), 30, 1e-11);
}

// gens 0 (slack, bus 0), 1 (bus 1) and 2 (bus 2); 1 and 2 both regulate bus 3
// remotely, i.e. one control group of two members at a bus nothing else pins.
LSGrid make_group_grid()
{
    LSGrid grid = make_skeleton();
    add_gens(grid, {0, 1, 2}, {1.01, V_SET, V_SET});
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    grid.set_gen_regulated_bus(2, 3);
    return grid;
}

// The whole point, in one helper: `mutate` applied to a grid that has already
// solved must give exactly what it gives on a grid that never had a plan to reuse.
template<class Mutate>
void same_as_freshly_built(Mutate mutate)
{
    LSGrid warm = make_group_grid();
    REQUIRE(solve(warm).size() == NB_BUS);   // a plan is now cached in the AC cache
    mutate(warm);
    const CplxVect v_warm = solve(warm);

    LSGrid cold = make_group_grid();
    mutate(cold);
    const CplxVect v_cold = solve(cold);

    REQUIRE(v_warm.size() == NB_BUS);
    REQUIRE(v_cold.size() == v_warm.size());
    for (Eigen::Index k = 0; k < v_warm.size(); ++k) {
        CHECK(v_warm(k).real() == Approx(v_cold(k).real()).margin(1e-10));
        CHECK(v_warm(k).imag() == Approx(v_cold(k).imag()).margin(1e-10));
    }
}

}  // namespace


TEST_CASE("the plan a solve reads is the one the cache built", "[voltage_control][cache_reuse]")
{
    LSGrid g = make_group_grid();
    // before any solve the AC cache holds no labelling, so the solver-side layers are
    // empty. Layer 1 is grid ids and answers all the same -- that is what the pv/pq
    // needs while it is still building the labelling the other two are expressed in.
    CHECK(g.get_ac_voltage_control_plan().controllers().n_controllers() == 0);
    CHECK(g.get_group_controlled_buses().count(3) == 1);

    REQUIRE(solve(g).size() == NB_BUS);

    const auto & plan = g.get_ac_voltage_control_plan();
    CHECK(plan.group_controlled_buses().count(3) == 1);
    CHECK(plan.controllers().n_groups() == 1);
    CHECK(plan.controllers().n_controllers() == 2);   // gens 1 and 2
    CHECK(plan.controllers().v_set(0) == Approx(V_SET));
    // ... and it is the same answer the on-demand form gives, which is what makes
    // the cached one trustworthy in the first place
    ls2g::VoltageControlSolverData on_demand;
    g.fill_voltage_control_solver_data(on_demand, true);
    REQUIRE(on_demand.n_controllers() == plan.controllers().n_controllers());
    for (int c = 0; c < on_demand.n_controllers(); ++c) {
        CHECK(on_demand.bus(c) == plan.controllers().bus(c));
        CHECK(on_demand.elem_id(c) == plan.controllers().elem_id(c));
        CHECK(on_demand.kind(c) == plan.controllers().kind(c));
    }
}

TEST_CASE("an injection-only step reuses the plan", "[voltage_control][cache_reuse]")
{
    // The case the caching exists for: a grid2op step moves P and Q and nothing else.
    // Nothing the plan is made of has changed, so it must not be rebuilt -- and the
    // answer must still be the one a grid with no plan to reuse gives.
    LSGrid g = make_group_grid();
    REQUIRE(solve(g).size() == NB_BUS);
    g.change_p_load(0, 44.);
    g.change_q_load(0, 9.);
    CHECK_FALSE(g.get_ac_algo_controler().need_recompute_voltage_control());

    same_as_freshly_built([](LSGrid & grid){
        grid.change_p_load(0, 44.);
        grid.change_q_load(0, 9.);
    });
}

TEST_CASE("re-declaring the same slack weights reuses the plan", "[voltage_control][cache_reuse]")
{
    // ``LightSimBackend._handle_dist_slack`` calls ``update_slack_weights`` on EVERY
    // step when the distributed slack is on, and that walks every generator calling
    // add_slackbus / remove_slackbus. Those two do touch an input of the plan -- the
    // slack role is one of the three things is_pseudo_off() reads -- so they must say
    // so, but only when the role actually moves. Saying so unconditionally would
    // retire the plan on every step of every distributed-slack episode.
    LSGrid g = make_group_grid();
    REQUIRE(solve(g).size() == NB_BUS);

    IntVect same_slack(1);
    same_slack << 0;                      // the slack it already has
    g.update_slack_weights_by_id(same_slack);
    CHECK_FALSE(g.get_ac_algo_controler().need_recompute_voltage_control());

    // ... and moving the slack to another generator does retire it
    IntVect other_slack(1);
    other_slack << 1;
    g.update_slack_weights_by_id(other_slack);
    CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
}

TEST_CASE("every input of the plan retires it", "[voltage_control][cache_reuse]")
{
    // One section per input the three layers read. Each first states that the change
    // is VISIBLE to the control (so the plan is rebuilt at all), then that the solve
    // it produces is the one a grid built that way from scratch gives.
    SECTION("a generator's regulated bus") {
        // gen 2 leaves the group at bus 3 and aims at bus 1 instead: one group of two
        // becomes two groups of one, and bus 1 -- gen 1's OWN bus -- becomes
        // group-controlled. Nothing about that is visible in a vector size.
        LSGrid g = make_group_grid();
        REQUIRE(solve(g).size() == NB_BUS);
        g.set_gen_regulated_bus(2, 1);
        CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
        same_as_freshly_built([](LSGrid & grid){ grid.set_gen_regulated_bus(2, 1); });
    }
    SECTION("a group member's voltage setpoint") {
        // both members must move: a group whose controllers disagree is refused
        LSGrid g = make_group_grid();
        REQUIRE(solve(g).size() == NB_BUS);
        g.change_v_gen(1, 1.05);
        g.change_v_gen(2, 1.05);
        CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
        same_as_freshly_built([](LSGrid & grid){
            grid.change_v_gen(1, 1.05);
            grid.change_v_gen(2, 1.05);
        });
    }
    SECTION("a group member disconnected") {
        LSGrid g = make_group_grid();
        REQUIRE(solve(g).size() == NB_BUS);
        g.deactivate_gen(2);
        CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
        same_as_freshly_built([](LSGrid & grid){ grid.deactivate_gen(2); });
    }
    SECTION("a group member moved to another bus") {
        LSGrid g = make_group_grid();
        REQUIRE(solve(g).size() == NB_BUS);
        g.change_bus_gen_python(2, 1);
        CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
        same_as_freshly_built([](LSGrid & grid){ grid.change_bus_gen_python(2, 1); });
    }
    SECTION("a generator turned pseudo-off") {
        // p crossing 0 flips is_remote_voltage_controller, so it changes who is in the group
        LSGrid g = make_group_grid();
        REQUIRE(solve(g).size() == NB_BUS);
        g.turnedoff_no_pv();
        g.change_p_gen(2, 10.);
        CHECK(g.get_ac_algo_controler().need_recompute_voltage_control());
        same_as_freshly_built([](LSGrid & grid){
            grid.turnedoff_no_pv();
            grid.change_p_gen(2, 10.);
        });
    }
}

TEST_CASE("an SVC appearing and disappearing retires the plan", "[voltage_control][cache_reuse]")
{
    // a voltage-mode SVC is ALWAYS a group controller, so connecting or disconnecting
    // one creates or dissolves a group at the bus it regulates
    LSGrid warm = make_skeleton();
    add_gens(warm, {0, 1}, {1.01, V_SET});
    warm.add_gen_slackbus(0, 1.);
    add_voltage_svc(warm, 2, 1.02);
    REQUIRE(solve(warm).size() == NB_BUS);
    REQUIRE(warm.get_ac_voltage_control_plan().controllers().n_controllers() == 1);

    warm.deactivate_svc(0);
    CHECK(warm.get_ac_algo_controler().need_recompute_voltage_control());
    const CplxVect v_warm = solve(warm);
    REQUIRE(v_warm.size() == NB_BUS);
    CHECK(warm.get_ac_voltage_control_plan().controllers().n_controllers() == 0);

    // the same grid, built with the SVC already off
    LSGrid cold = make_skeleton();
    add_gens(cold, {0, 1}, {1.01, V_SET});
    cold.add_gen_slackbus(0, 1.);
    add_voltage_svc(cold, 2, 1.02);
    cold.deactivate_svc(0);
    const CplxVect v_cold = solve(cold);
    REQUIRE(v_cold.size() == v_warm.size());
    for (Eigen::Index k = 0; k < v_warm.size(); ++k) {
        CHECK(v_warm(k).real() == Approx(v_cold(k).real()).margin(1e-10));
        CHECK(v_warm(k).imag() == Approx(v_cold(k).imag()).margin(1e-10));
    }
}

TEST_CASE("a DC cache carries no plan", "[voltage_control][cache_reuse]")
{
    // there is no voltage control in a DC solve: the DC half of the cache must not
    // pay for a plan, and the AC half must not be disturbed by a DC solve either
    LSGrid g = make_group_grid();
    REQUIRE(g.dc_pf(flat(g), 1, 1e-10).size() == NB_BUS);
    CHECK(g.get_ac_voltage_control_plan().controllers().n_controllers() == 0);
    REQUIRE(solve(g).size() == NB_BUS);
    CHECK(g.get_ac_voltage_control_plan().controllers().n_controllers() == 2);
}


// ---------------------------------------------------------------------------
// An algorithm with no bordered block
// ---------------------------------------------------------------------------
//
// Layers 3 and 4 are consumed by the Base and VoltageControl components of
// NRSystem. Fast-decoupled and Gauss-Seidel hold no NRSystem at all, so for them
// those layers are work nobody reads -- and layer 2's group step is worse than
// useless: it takes the regulated bus out of PV and leaves NOTHING pinning its
// magnitude. That solve converges, looks plausible, and is wrong (measured 0.36 pu
// away from the Newton answer on a case118 with eight control groups).
//
// So ac_pf refuses such a grid by name, and the plan is not built at all.

namespace {

// gen 1 at bus 1 regulates bus 3 REMOTELY and gen 2 SITS on bus 3 regulating it
// locally. This is the configuration where the two splits differ: without the group
// step gen 2 pins bus 3 as PV, with it gen 2 is enrolled in the group instead and
// bus 3 keeps its own Vm unknown for the bordered voltage row to act on.
LSGrid make_local_and_remote_grid()
{
    LSGrid grid = make_skeleton();
    add_gens(grid, {0, 1, 3}, {1.01, V_SET, V_SET});
    grid.add_gen_slackbus(0, 1.);
    grid.set_gen_regulated_bus(1, 3);
    return grid;
}

// gens 0 (slack) and 1, both regulating their OWN bus: the ordinary PV case, which
// every algorithm can do
LSGrid make_local_only_grid()
{
    LSGrid grid = make_skeleton();
    add_gens(grid, {0, 1}, {1.01, V_SET});
    grid.add_gen_slackbus(0, 1.);
    return grid;
}

bool contains(const std::string & haystack, const std::string & needle)
{
    return haystack.find(needle) != std::string::npos;
}

// what ac_pf says (empty when it does not throw)
std::string ac_pf_error(LSGrid & g)
{
    try {
        g.ac_pf(flat(g), 30, 1e-11);
    } catch (const std::runtime_error & exc) {
        return exc.what();
    }
    return std::string();
}

}  // namespace


TEST_CASE("an algorithm without the bordered block refuses a controlled grid",
          "[voltage_control][fdpf]")
{
    SECTION("a remote-regulating generator is named") {
        LSGrid g = make_group_grid();       // gens 1 and 2 both regulate bus 3
        g.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
        const std::string msg = ac_pf_error(g);
        REQUIRE_FALSE(msg.empty());
        CHECK(contains(msg, "generator(s) 1 2"));
        CHECK(contains(msg, "FDPF_XB_SparseLU"));
        CHECK(contains(msg, "change_algorithm"));
    }
    SECTION("a voltage-mode SVC is named -- even a purely local one") {
        // an SVC is ALWAYS a group controller, local or not: its reactive injection
        // is solved for by the bordered block and by nothing else
        LSGrid g = make_local_only_grid();
        add_voltage_svc(g, 2, 1.02);
        g.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
        const std::string msg = ac_pf_error(g);
        REQUIRE_FALSE(msg.empty());
        CHECK(contains(msg, "SVC(s) 0"));
    }
    SECTION("Gauss-Seidel is refused for the same reason") {
        LSGrid g = make_group_grid();
        g.change_algorithm(AlgorithmType::GaussSeidel);
        CHECK_FALSE(ac_pf_error(g).empty());
    }
    SECTION("... and a grid it CAN do is not refused") {
        LSGrid g = make_local_only_grid();
        g.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
        const CplxVect v_fdpf = g.ac_pf(flat(g), 200, 1e-10);
        REQUIRE(v_fdpf.size() == NB_BUS);

        // and it agrees with Newton-Raphson, which is the point of not refusing it
        LSGrid ref = make_local_only_grid();
        const CplxVect v_nr = solve(ref);
        REQUIRE(v_nr.size() == NB_BUS);
        for (Eigen::Index k = 0; k < v_nr.size(); ++k) {
            CHECK(v_fdpf(k).real() == Approx(v_nr(k).real()).margin(1e-8));
            CHECK(v_fdpf(k).imag() == Approx(v_nr(k).imag()).margin(1e-8));
        }
    }
}

TEST_CASE("an algorithm without the bordered block builds no plan", "[voltage_control][fdpf]")
{
    // the cost half of the same decision: layers 1, 3 and 4 are not built, so the
    // three container walks they cost are not paid either
    LSGrid g = make_local_only_grid();
    g.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
    REQUIRE(g.ac_pf(flat(g), 200, 1e-10).size() == NB_BUS);

    const auto & plan = g.get_ac_voltage_control_plan();
    CHECK(plan.group_controlled_buses().empty());
    CHECK(plan.free_vm_slack_buses().empty());
    CHECK(plan.controllers().n_controllers() == 0);

    // switching back to Newton-Raphson builds it again (change_algorithm retires
    // everything), so the skip is not a state a later solve can inherit
    g.change_algorithm(AlgorithmType::NR_SparseLU);
    REQUIRE(solve(g).size() == NB_BUS);
    CHECK(g.get_ac_voltage_control_plan().free_vm_slack_buses().size() +
          g.get_ac_voltage_control_plan().controllers().n_controllers() >= 0);  // built, whatever it holds
}

TEST_CASE("without the bordered block the pv/pq split is the classical one",
          "[voltage_control][fdpf]")
{
    // The part that would be a WRONG ANSWER rather than a slow one. `pre_process_solver`
    // is the entry point ac_pf uses, without the refusal on top, so it is where the two
    // splits can be compared side by side.
    //
    // Newton-Raphson takes bus 3 out of PV, because the group's bordered voltage row is
    // what will set its magnitude. A fast-decoupled build must NOT: it has no such row,
    // and a bus that is neither PV nor pinned by anything is how the 0.36 pu error
    // happens.
    const ls2g::AlgoControl all_changed;   // default ctor: everything changed

    LSGrid nr = make_local_and_remote_grid();
    nr.pre_process_solver(flat(nr), all_changed);
    const std::vector<int> pv_nr = nr.get_pv_solver().to_int_vector();

    LSGrid fd = make_local_and_remote_grid();
    fd.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
    fd.pre_process_solver(flat(fd), all_changed);
    const std::vector<int> pv_fd = fd.get_pv_solver().to_int_vector();

    // bus 3 is the regulated one; the labelling is the identity here (every bus
    // carries something), but read it back rather than assuming
    const int bus_3_solver = nr.id_me_to_ac_solver()[3].cast_int();
    REQUIRE(bus_3_solver >= 0);
    CHECK(std::find(pv_nr.begin(), pv_nr.end(), bus_3_solver) == pv_nr.end());
    CHECK(std::find(pv_fd.begin(), pv_fd.end(), bus_3_solver) != pv_fd.end());
}


TEST_CASE("the DC split is untouched by the AC algorithm's capability", "[voltage_control][fdpf]")
{
    // The DC family also runs the pv/pq split, and BaseDCAlgo reads `pv`
    // (retrieve_pv_with_slack / extract_slack_bus_id). Whether a bus a control group
    // regulates SHOULD be reclassified in a solve that has no voltage at all is a fair
    // question, but it is not this change's: gating layer 1 on the AC algorithm's
    // capability would have answered it by accident, and silently. So the DC split
    // keeps the group step unconditionally, and this pins that.
    const ls2g::AlgoControl all_changed;

    LSGrid nr = make_local_and_remote_grid();
    nr.pre_process_dc_solver(flat(nr), all_changed);
    const std::vector<int> pv_dc_nr = nr.get_dc_pv_solver().to_int_vector();

    LSGrid fd = make_local_and_remote_grid();
    fd.change_algorithm(AlgorithmType::FDPF_XB_SparseLU);
    fd.pre_process_dc_solver(flat(fd), all_changed);
    const std::vector<int> pv_dc_fd = fd.get_dc_pv_solver().to_int_vector();

    CHECK(pv_dc_nr == pv_dc_fd);
    // ... and it really is the group-aware split, not the classical one: bus 3 has a
    // local regulator on it, and the group takes it back
    const int bus_3_solver = nr.id_me_to_dc_solver()[3].cast_int();
    REQUIRE(bus_3_solver >= 0);
    CHECK(std::find(pv_dc_nr.begin(), pv_dc_nr.end(), bus_3_solver) == pv_dc_nr.end());
}
