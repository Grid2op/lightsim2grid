// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The projection of the switches onto the elements (LSGrid::project_switches,
// set_switch_open, update_switches): a switch moves, the labels of its
// substation are recomputed, and every terminal there goes where its component
// says -- through the same mutators a topology-vector action uses, so the bus
// counts, the change flags and the powerflow are the ones that action would
// have produced.

#include <cmath>
#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::AlgoControl;
using ls2g::CplxVect;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::SwitchKind;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using BoolArr = Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor>;
using IntArr = Eigen::Array<int, Eigen::Dynamic, Eigen::RowMajor>;

IntVect ivec(std::initializer_list<int> values)
{
    IntVect res(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for(const int v : values) res(i++) = v;
    return res;
}

// Three substations of two busbars (buses s and s+3), a slack generator in
// substation 0, lines 0-1, 1-2 and a second 0-1, a load in substation 2; wired
// for update_topo (pos_topo_vect: load 0, gen 1, line ends 1 at 2-4, ends 2 at
// 5-7). The detailed topology, switches numbered grid-wide in this order:
//   sub 0 (6 nodes): sections A (node 0) and B (node 1)
//     0: coupler A-B, open        1: gen (node 2) -> A, closed
//     2: gen -> B, open           3: line 0 bay (node 3) -> A, closed
//     4: line 0 bay -> B, open    5: internal connection bay <-> terminal (node 4)
//     6: line 2 end 1 (node 5) -> A, closed (its only way in: A stays live)
//   sub 1 (4 nodes): one section (node 0)
//     7: line 0 end 2 (node 1) breaker, closed   8: line 1 end 1 (node 2) breaker, closed
//     9: line 2 end 2 (node 3) breaker, closed
//   sub 2 (3 nodes): one section (node 0)
//     10: load (node 1) breaker, closed          11: line 1 end 2 (node 2) breaker, closed
// In this state every element is exactly where the plain grid put it: gen,
// line 0 and line 2 on bus 0, line 1 on 1-2, load on bus 2.
LSGrid make_switch_grid()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    RealVect bus_vn_kv = RealVect::Constant(6, 138.);
    grid.init_bus(3, 2, bus_vn_kv, 0, 0);
    RealVect branch_r(3), branch_x(3);
    branch_r << 0.01, 0.01, 0.02;
    branch_x << 0.1, 0.1, 0.2;
    grid.init_powerlines(branch_r, branch_x, CplxVect::Zero(3), ivec({0, 1, 0}), ivec({1, 2, 1}));
    RealVect load_p(1), load_q(1);
    load_p << 50.;
    load_q << 10.;
    grid.init_loads(load_p, load_q, ivec({2}));
    RealVect gen_p(1), gen_v(1), gen_min_q(1), gen_max_q(1);
    gen_p << 0.;
    gen_v << 1.02;
    gen_min_q << -1000.;
    gen_max_q << 1000.;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, ivec({0}));
    grid.add_gen_slackbus(0, 1.);

    grid.set_load_to_subid(ivec({2}));
    grid.set_gen_to_subid(ivec({0}));
    grid.set_line_to_sub1_id(ivec({0, 1, 0}));
    grid.set_line_to_sub2_id(ivec({1, 2, 1}));
    grid.set_load_pos_topo_vect(ivec({0}));
    grid.set_gen_pos_topo_vect(ivec({1}));
    grid.set_line_pos1_topo_vect(ivec({2, 3, 4}));
    grid.set_line_pos2_topo_vect(ivec({5, 6, 7}));

    const int B = static_cast<int>(SwitchKind::BREAKER);
    const int D = static_cast<int>(SwitchKind::DISCONNECTOR);
    const int I = static_cast<int>(SwitchKind::INTERNAL_CONNECTION);
    grid.init_detailed_topology(
        ivec({6, 4, 3}),
        ivec({0, 0, 1, 2}), ivec({0, 1, 0, 0}),
        ivec({0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2}),
        ivec({0, 2, 2, 3, 3, 4, 5, 1, 2, 3, 1, 2}),
        ivec({1, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0}),
        {B, D, D, D, D, I, D, B, B, B, B, B},
        {true, false, true, false, true, false, false, false, false, false, false, false},
        {true, false, false, false, false, false, false, true, true, true, true, true});
    grid.set_load_to_node_id(ivec({1}));
    grid.set_gen_to_node_id(ivec({2}));
    grid.set_line_to_node1_id(ivec({4, 2, 5}));
    grid.set_line_to_node2_id(ivec({1, 2, 3}));
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(grid.get_init_vm_pu(), 0.));
}

void solve_and_settle(LSGrid & grid)
{
    REQUIRE(grid.dc_pf(flat_start(grid), 30, 1e-10).size() > 0);
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() > 0);
    grid.unset_changes();
}

std::vector<bool> flags_of(const AlgoControl & c)
{
    return {c.has_dimension_changed(), c.has_pv_changed(), c.has_pq_changed(),
            c.has_slack_participate_changed(), c.need_reset_solver(), c.need_recompute_sbus(),
            c.need_recompute_ybus(), c.has_slack_weight_changed(), c.has_v_changed(),
            c.has_ybus_some_coeffs_zero(), c.has_one_el_changed_bus(), c.has_voltage_control_changed()};
}

void require_same_flags(const LSGrid & a, const LSGrid & b)
{
    CHECK(flags_of(a.get_ac_algo_controler()) == flags_of(b.get_ac_algo_controler()));
    CHECK(flags_of(a.get_dc_algo_controler()) == flags_of(b.get_dc_algo_controler()));
}

void require_counts_exact(LSGrid & grid, const char * what)
{
    INFO("after: " << what);
    const std::vector<std::size_t> incremental = grid.get_substations().get_nb_elements_per_bus();
    CHECK(grid.get_substations().connected_bus_count_is_exact());
    grid.recompute_bus_element_counts();
    const std::vector<std::size_t> & from_scratch = grid.get_substations().get_nb_elements_per_bus();
    REQUIRE(incremental.size() == from_scratch.size());
    for(std::size_t b = 0; b < from_scratch.size(); ++b){
        INFO("bus " << b);
        CHECK(incremental[b] == from_scratch[b]);
    }
}

// some entries of the grid2op topology vector, the rest untouched
void set_topo(LSGrid & grid, std::initializer_list<int> positions, std::initializer_list<int> local_buses)
{
    BoolArr changed = BoolArr::Constant(8, false);
    IntArr values = IntArr::Zero(8);
    auto it_bus = local_buses.begin();
    for(const int pos : positions){ changed(pos) = true; values(pos) = *it_bus++; }
    grid.update_topo(changed, values);
}

void set_switches(LSGrid & grid, std::initializer_list<int> ids, std::initializer_list<bool> open)
{
    BoolArr changed = BoolArr::Constant(12, false);
    BoolArr values = BoolArr::Constant(12, false);
    auto it_open = open.begin();
    for(const int id : ids){ changed(id) = true; values(id) = *it_open++; }
    grid.update_switches(changed, values);
}

real_type max_abs_diff(const CplxVect & a, const CplxVect & b)
{
    REQUIRE(a.size() == b.size());
    return (a - b).cwiseAbs().maxCoeff();
}

}  // anonymous namespace

TEST_CASE("projecting a grid already where its switches say changes nothing", "[LSGrid][switches]")
{
    LSGrid grid = make_switch_grid();
    solve_and_settle(grid);
    grid.project_switches();
    CHECK(grid.get_ac_algo_controler().nothing_changed());
    CHECK(grid.get_dc_algo_controler().nothing_changed());
    CHECK(grid.get_generators()[0].bus_id == 0);
    CHECK(grid.get_lines()[0].bus_1_id == 0);
    CHECK(grid.get_lines()[0].bus_2_id == 1);
    CHECK(grid.get_lines()[1].bus_2_id == 2);
    CHECK(grid.get_lines()[2].bus_1_id == 0);
    CHECK(grid.get_lines()[2].bus_2_id == 1);
    CHECK(grid.get_loads()[0].bus_id == 2);
    require_counts_exact(grid, "a no-op projection");
}

TEST_CASE("moving a bay to the other section is the same as the topology-vector action", "[LSGrid][switches]")
{
    // the generator's and line 0's disconnectors to A open, to B close: both now
    // stand on section B, which becomes bus 3 (substation 0, local bus 2); section
    // A keeps line 2's end and stays bus 0
    LSGrid by_switch = make_switch_grid();
    LSGrid by_topo = make_switch_grid();
    solve_and_settle(by_switch);
    solve_and_settle(by_topo);

    set_switches(by_switch, {1, 2, 3, 4}, {true, false, true, false});
    set_topo(by_topo, {1, 2}, {2, 2});

    CHECK(by_switch.get_generators()[0].bus_id == 3);
    CHECK(by_switch.get_generators()[0].connected);
    CHECK(by_switch.get_lines()[0].bus_1_id == 3);
    CHECK(by_switch.get_lines()[2].bus_1_id == 0);
    CHECK(by_switch.get_substation_topology(0).nb_buses() == 2);
    CHECK(by_switch.get_node_bus()(0) == 0);
    CHECK(by_switch.get_node_bus()(1) == 3);
    CHECK(by_switch.get_busbar_sections()[0].bus_id == 0);
    CHECK(by_switch.get_busbar_sections()[1].bus_id == 3);
    CHECK(by_topo.get_generators()[0].bus_id == 3);
    CHECK(by_topo.get_lines()[0].bus_1_id == 3);
    require_same_flags(by_switch, by_topo);
    require_counts_exact(by_switch, "bay moved by switches");

    // and the two grids solve to the same voltages
    const CplxVect v_switch = by_switch.ac_pf(flat_start(by_switch), 30, 1e-10);
    const CplxVect v_topo = by_topo.ac_pf(flat_start(by_topo), 30, 1e-10);
    REQUIRE(v_switch.size() > 0);
    REQUIRE(v_topo.size() > 0);
    CHECK(max_abs_diff(v_switch, v_topo) < 1e-12);
    // bus 3 is live, and bus 0 still is (line 2's end)
    CHECK(by_switch.get_substations().is_bus_connected(ls2g::GridModelBusId(3)));
    CHECK(by_switch.get_substations().is_bus_connected(ls2g::GridModelBusId(0)));

    SECTION("closing the coupler merges the sections back into one bus")
    {
        by_switch.unset_changes();
        CHECK(by_switch.set_switch_open(0, false));
        CHECK(by_switch.get_substation_topology(0).nb_buses() == 1);
        // section A comes first in the numbering: the merged bus is bus 0
        CHECK(by_switch.get_generators()[0].bus_id == 0);
        CHECK(by_switch.get_lines()[0].bus_1_id == 0);
        CHECK(by_switch.get_busbar_sections()[1].bus_id == 0);
        CHECK(by_switch.get_ac_algo_controler().has_dimension_changed());  // bus 3 left
        require_counts_exact(by_switch, "coupler closed");
        CHECK(by_switch.ac_pf(flat_start(by_switch), 30, 1e-10).size() > 0);
    }
}

TEST_CASE("opening a feeder's breaker disconnects it, like the topology-vector action", "[LSGrid][switches]")
{
    LSGrid by_switch = make_switch_grid();
    LSGrid by_topo = make_switch_grid();
    solve_and_settle(by_switch);
    solve_and_settle(by_topo);

    CHECK(by_switch.set_switch_open(10, true));
    set_topo(by_topo, {0}, {-1});

    CHECK_FALSE(by_switch.get_loads()[0].connected);
    CHECK(by_switch.get_loads()[0].bus_id == -1);
    // its section keeps line 1's end: still a bus
    CHECK(by_switch.get_busbar_sections()[3].connected);
    CHECK(by_switch.get_lines()[1].bus_2_id == 2);
    require_same_flags(by_switch, by_topo);
    require_counts_exact(by_switch, "load breaker opened");
    CHECK(max_abs_diff(by_switch.ac_pf(flat_start(by_switch), 30, 1e-10),
                       by_topo.ac_pf(flat_start(by_topo), 30, 1e-10)) < 1e-12);

    SECTION("closing it again brings the load back onto its bus")
    {
        by_switch.unset_changes();
        CHECK(by_switch.set_switch_open(10, false));
        CHECK(by_switch.get_loads()[0].connected);
        CHECK(by_switch.get_loads()[0].bus_id == 2);
        require_counts_exact(by_switch, "load breaker closed");
    }
    SECTION("a switch already in that position does not move")
    {
        by_switch.unset_changes();
        CHECK_FALSE(by_switch.set_switch_open(10, true));
        CHECK(by_switch.get_ac_algo_controler().nothing_changed());
    }
}

TEST_CASE("opening a line end's breaker: synched ends drop the line, unsynched leave it half-open", "[LSGrid][switches]")
{
    SECTION("synched (the default): both ends off, like the topology vector")
    {
        LSGrid by_switch = make_switch_grid();
        LSGrid by_topo = make_switch_grid();
        solve_and_settle(by_switch);
        solve_and_settle(by_topo);
        CHECK(by_switch.set_switch_open(7, true));
        set_topo(by_topo, {5}, {-1});
        CHECK_FALSE(by_switch.get_lines()[0].connected_global);
        CHECK_FALSE(by_switch.get_lines()[0].connected_1);
        CHECK_FALSE(by_switch.get_lines()[0].connected_2);
        require_same_flags(by_switch, by_topo);
        require_counts_exact(by_switch, "synched line end opened");
        CHECK(max_abs_diff(by_switch.ac_pf(flat_start(by_switch), 30, 1e-10),
                           by_topo.ac_pf(flat_start(by_topo), 30, 1e-10)) < 1e-12);
        by_switch.unset_changes();
        CHECK(by_switch.set_switch_open(7, false));
        CHECK(by_switch.get_lines()[0].connected_global);
        CHECK(by_switch.get_lines()[0].bus_1_id == 0);
        CHECK(by_switch.get_lines()[0].bus_2_id == 1);
        require_counts_exact(by_switch, "synched line end closed again");
        CHECK(by_switch.ac_pf(flat_start(by_switch), 30, 1e-10).size() > 0);
    }
    SECTION("unsynched: the line stays on with its far end open")
    {
        LSGrid grid = make_switch_grid();
        grid.set_synch_status_both_side(false);
        solve_and_settle(grid);
        CHECK(grid.set_switch_open(7, true));
        CHECK(grid.get_lines()[0].connected_global);
        CHECK(grid.get_lines()[0].connected_1);
        CHECK_FALSE(grid.get_lines()[0].connected_2);
        CHECK(grid.get_lines()[0].bus_1_id == 0);
        require_counts_exact(grid, "unsynched line end opened");
    }
}

TEST_CASE("only the substations whose switches moved are projected", "[LSGrid][switches]")
{
    LSGrid grid = make_switch_grid();
    solve_and_settle(grid);
    // an escape hatch in substation 0: the generator taken out by hand
    grid.deactivate_gen(0);
    REQUIRE_FALSE(grid.get_generators()[0].connected);
    // a switch move in substation 2 leaves it alone
    CHECK(grid.set_switch_open(10, true));
    CHECK_FALSE(grid.get_generators()[0].connected);
    CHECK_FALSE(grid.get_loads()[0].connected);
    require_counts_exact(grid, "switch elsewhere");
    // a projection of the whole grid puts the generator back where its switches say
    grid.project_switches();
    CHECK(grid.get_generators()[0].connected);
    CHECK(grid.get_generators()[0].bus_id == 0);
    CHECK_FALSE(grid.get_loads()[0].connected);  // its breaker is still open
    require_counts_exact(grid, "whole projection");
}

TEST_CASE("the switch positions and the projected elements survive copy and state", "[LSGrid][switches][serialization]")
{
    LSGrid grid = make_switch_grid();
    set_switches(grid, {1, 2, 3, 4}, {true, false, true, false});
    REQUIRE(grid.get_generators()[0].bus_id == 3);
    const IntVect node_bus = grid.get_node_bus();

    SECTION("copy")
    {
        LSGrid copied = grid.copy();
        CHECK(copied.get_node_bus() == node_bus);
        CHECK(copied.get_generators()[0].bus_id == 3);
        CHECK(copied.get_switches()[1].open);
        CHECK_FALSE(copied.get_switches()[2].open);
        // and they keep behaving the same afterwards
        CHECK(copied.set_switch_open(0, false));
        CHECK(grid.set_switch_open(0, false));
        CHECK(copied.get_node_bus() == grid.get_node_bus());
        CHECK(copied.get_generators()[0].bus_id == grid.get_generators()[0].bus_id);
    }
    SECTION("get_state / set_state")
    {
        LSGrid::StateRes state = grid.get_state();
        LSGrid restored;
        restored.set_state(state);
        // a restored grid starts with its bus counts disarmed (recounted before any
        // use); arm them here so the exactness check below has something to compare
        restored.recompute_bus_element_counts();
        CHECK(restored.get_node_bus() == node_bus);
        CHECK(restored.get_generators()[0].bus_id == 3);
        CHECK(restored.get_switches()[1].open);
        CHECK(restored.set_switch_open(0, false));
        CHECK(restored.get_generators()[0].bus_id == 0);
        require_counts_exact(restored, "restored then switched");
        CHECK(restored.ac_pf(flat_start(restored), 30, 1e-10).size() > 0);
    }
}

TEST_CASE("the busbar sections read their bus' voltage after a powerflow", "[LSGrid][switches]")
{
    LSGrid grid = make_switch_grid();
    CHECK_FALSE(grid.get_busbar_sections()[0].has_res);
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() > 0);
    const RealVect vm = grid.get_Vm();
    const RealVect va = grid.get_Va();
    const auto sections = grid.get_busbar_sections();
    REQUIRE(sections[0].has_res);
    CHECK(sections[0].bus_id == 0);
    CHECK(std::abs(sections[0].res_v_kv - vm(0) * 138.) < 1e-9);
    CHECK(std::abs(sections[0].res_theta_deg - va(0) * ls2g::BaseConstants::my_180_pi_) < 1e-9);
    CHECK(sections[2].bus_id == 1);
    CHECK(std::abs(sections[2].res_v_kv - vm(1) * 138.) < 1e-9);
    // section B is in no valid component: no voltage to report
    CHECK_FALSE(sections[1].connected);
    CHECK(std::isnan(sections[1].res_v_kv));
}

TEST_CASE("update_switches refuses what it cannot apply", "[LSGrid][switches]")
{
    SECTION("no detailed topology")
    {
        LSGrid grid;
        grid.set_sn_mva(100.);
        grid.set_init_vm_pu(1.0);
        RealVect vn = RealVect::Constant(2, 138.);
        grid.init_bus(2, 1, vn, 0, 0);
        BoolArr none(0);
        CHECK_THROWS_MATCHES(grid.update_switches(none, none), std::runtime_error,
                             MessageMatches(ContainsSubstring("no detailed topology")));
        CHECK_THROWS_AS(grid.set_switch_open(0, true), std::runtime_error);
        CHECK_THROWS_AS(grid.project_switches(), std::runtime_error);
    }
    SECTION("arrays of the wrong length")
    {
        LSGrid grid = make_switch_grid();
        BoolArr short_arr = BoolArr::Constant(3, false);
        CHECK_THROWS_MATCHES(grid.update_switches(short_arr, short_arr), std::runtime_error,
                             MessageMatches(ContainsSubstring("one entry per switch")));
    }
    SECTION("a bad switch id, and an internal connection")
    {
        LSGrid grid = make_switch_grid();
        solve_and_settle(grid);
        CHECK_THROWS_AS(grid.set_switch_open(12, true), std::out_of_range);
        CHECK_THROWS_MATCHES(grid.set_switch_open(5, true), std::runtime_error,
                             MessageMatches(ContainsSubstring("internal connection")));
        CHECK(grid.get_ac_algo_controler().nothing_changed());
        CHECK(grid.get_lines()[0].bus_1_id == 0);
    }
}
