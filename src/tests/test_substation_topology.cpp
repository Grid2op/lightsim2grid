// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The detailed topology (switches inside each substation) as DATA: what a
// SubstationTopology declares and refuses, how SubstationContainer numbers the
// substations' local ids grid-wide, how the node ids ride on the element
// containers, and that all of it survives a state / binary round trip -- while a
// grid that has no detailed topology is exactly what it was.
//
// Nothing here labels a component or moves an element: that is the projection,
// tested separately.

#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "SubstationContainer.hpp"
#include "SubstationTopology.hpp"
#include "case_exotic_elements.hpp"
#include "test_helpers.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::BusbarSectionInfo;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::SubstationContainer;
using ls2g::SubstationInfo;
using ls2g::SubstationTopology;
using ls2g::SwitchInfo;
using ls2g::SwitchKind;
using ls2g::Terminal;
using ls2g::TerminalKind;
using ls2g_test::TempFile;

namespace {

IntVect ivec(std::initializer_list<int> values)
{
    IntVect res(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for(const int v : values) res(i++) = v;
    return res;
}

// Substation 0 of the fixture below, on its own: two busbar sections (nodes 0 and
// 1) with an open coupler, a generator bay (node 2, disconnector to section A), a
// line bay whose bay node (3) can reach either section and whose terminal (node
// 4) sits behind an internal connection.
SubstationTopology make_sub0_topology()
{
    SubstationTopology topo;
    topo.init(5,
              ivec({0, 1}),
              ivec({0, 2, 3, 3, 4}),
              ivec({1, 0, 0, 1, 3}),
              {SwitchKind::BREAKER, SwitchKind::DISCONNECTOR, SwitchKind::DISCONNECTOR,
               SwitchKind::DISCONNECTOR, SwitchKind::INTERNAL_CONNECTION},
              {true, false, false, true, false},
              {true, false, false, false, false});
    return topo;
}

// Three substations of two busbars each (6 buses), two lines 0-1 and 1-2, a load
// at substation 2 and a slack generator at substation 0 -- the layout of
// test_bus_element_count's two-busbar grid, with substation ids set on every
// terminal so a node id has something to be local to.
LSGrid make_grid_without_topology()
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    RealVect bus_vn_kv(6);
    bus_vn_kv << 138., 138., 138., 138., 138., 138.;
    grid.init_bus(3, 2, bus_vn_kv, 0, 0);

    RealVect branch_r(2), branch_x(2);
    branch_r << 0.01, 0.01;
    branch_x << 0.1, 0.1;
    const ls2g::CplxVect branch_h = ls2g::CplxVect::Zero(2);
    grid.init_powerlines(branch_r, branch_x, branch_h, ivec({0, 1}), ivec({1, 2}));

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
    grid.set_line_to_sub1_id(ivec({0, 1}));
    grid.set_line_to_sub2_id(ivec({1, 2}));
    return grid;
}

// The same grid with its detailed topology declared:
//   sub 0: make_sub0_topology (gen on node 2, line 0 side 1 on node 4)
//   sub 1: one section (node 0), line 0 side 2 on node 1 and line 1 side 1 on
//          node 2, each behind a closed breaker
//   sub 2: one section (node 0), the load on node 1 with NO switch at all, line 1
//          side 2 on node 2 behind a closed breaker
// Busbar sections and switches are given sorted by substation, so their position
// is their grid-wide id: sections 0-1 in sub 0, 2 in sub 1, 3 in sub 2; switches
// 0-4 in sub 0, 5-6 in sub 1, 7 in sub 2.
LSGrid make_grid_with_topology()
{
    LSGrid grid = make_grid_without_topology();
    grid.init_detailed_topology(
        ivec({5, 3, 3}),
        ivec({0, 0, 1, 2}), ivec({0, 1, 0, 0}),
        ivec({0, 0, 0, 0, 0, 1, 1, 2}),
        ivec({0, 2, 3, 3, 4, 1, 2, 2}),
        ivec({1, 0, 0, 1, 3, 0, 0, 0}),
        {static_cast<int>(SwitchKind::BREAKER), static_cast<int>(SwitchKind::DISCONNECTOR),
         static_cast<int>(SwitchKind::DISCONNECTOR), static_cast<int>(SwitchKind::DISCONNECTOR),
         static_cast<int>(SwitchKind::INTERNAL_CONNECTION), static_cast<int>(SwitchKind::BREAKER),
         static_cast<int>(SwitchKind::BREAKER), static_cast<int>(SwitchKind::BREAKER)},
        {true, false, false, true, false, false, false, false},
        {true, false, false, false, false, true, true, true});
    grid.set_load_to_node_id(ivec({1}));
    grid.set_gen_to_node_id(ivec({2}));
    grid.set_line_to_node1_id(ivec({4, 2}));
    grid.set_line_to_node2_id(ivec({1, 2}));
    return grid;
}

bool same_terminal(const Terminal & term, TerminalKind kind, int el_id, int node)
{
    return (term.kind == kind) && (term.el_id == el_id) && (term.node == node);
}

// The topology half of two grids' states: the substation container (which owns
// the switches) and the node id of every terminal. NB a whole LSGrid::StateRes
// never compares equal to its own round trip: it carries NaN sentinels (voltage
// limits, active power limits), and NaN != NaN.
bool same_topology_state(const LSGrid & a, const LSGrid & b)
{
    if(!(a.get_substations().get_state() == b.get_substations().get_state())) return false;
    if(a.get_loads().nb() != b.get_loads().nb()) return false;
    for(int i = 0; i < a.get_loads().nb(); ++i) if(a.get_loads()[i].node_id != b.get_loads()[i].node_id) return false;
    if(a.get_generators().nb() != b.get_generators().nb()) return false;
    for(int i = 0; i < a.get_generators().nb(); ++i) if(a.get_generators()[i].node_id != b.get_generators()[i].node_id) return false;
    if(a.get_lines().nb() != b.get_lines().nb()) return false;
    for(int i = 0; i < a.get_lines().nb(); ++i){
        if(a.get_lines()[i].node_1_id != b.get_lines()[i].node_1_id) return false;
        if(a.get_lines()[i].node_2_id != b.get_lines()[i].node_2_id) return false;
    }
    return true;
}

}  // anonymous namespace

// ---- SubstationTopology on its own ------------------------------------------

TEST_CASE("a SubstationTopology declares nodes, busbar sections and switches", "[SubstationTopology]")
{
    const SubstationTopology topo = make_sub0_topology();
    CHECK(topo.nb_nodes() == 5);
    CHECK(topo.nb_busbar_sections() == 2);
    CHECK(topo.nb_switches() == 5);
    CHECK(topo.bbs_node(1) == 1);
    CHECK(topo.sw_node1(4) == 4);
    CHECK(topo.sw_node2(4) == 3);
    CHECK(topo.sw_kind(4) == SwitchKind::INTERNAL_CONNECTION);
    CHECK(topo.is_open(0));
    CHECK_FALSE(topo.is_open(1));
    CHECK(topo.is_retained(0));
    CHECK_FALSE(topo.is_retained(1));
    // no names set: empty, not an error
    CHECK(topo.sw_name(0) == "");
    CHECK(topo.bbs_name(0) == "");
    CHECK(topo.terminals().empty());
    // out-of-range ids are refused, never read
    CHECK_THROWS_AS(topo.bbs_node(2), std::out_of_range);
    CHECK_THROWS_AS(topo.sw_node1(5), std::out_of_range);
    CHECK_THROWS_AS(topo.sw_kind(-1), std::out_of_range);
}

TEST_CASE("a SubstationTopology round-trips through its state", "[SubstationTopology][serialization]")
{
    SubstationTopology topo = make_sub0_topology();
    topo.set_sw_names({"coupler", "gen_disc", "line_disc_A", "line_disc_B", "line_ic"});
    topo.set_bbs_names({"A", "B"});
    SubstationTopology::StateRes state = topo.get_state();
    SubstationTopology restored;
    restored.set_state(state);
    CHECK(restored.get_state() == topo.get_state());
    CHECK(restored.sw_name(2) == "line_disc_A");
    CHECK(restored.bbs_name(1) == "B");
    CHECK(restored.sw_kind(4) == SwitchKind::INTERNAL_CONNECTION);

    SECTION("an empty topology (a substation with nothing to describe) round-trips too")
    {
        const SubstationTopology empty;
        SubstationTopology::StateRes empty_state = empty.get_state();
        SubstationTopology restored_empty;
        restored_empty.set_state(empty_state);
        CHECK(restored_empty.nb_nodes() == 0);
        CHECK(restored_empty.get_state() == empty.get_state());
    }
}

TEST_CASE("a SubstationTopology refuses an inconsistent declaration", "[SubstationTopology]")
{
    SubstationTopology topo;
    const std::vector<SwitchKind> one_breaker{SwitchKind::BREAKER};

    SECTION("a busbar section on a node that does not exist")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({2}), IntVect(), IntVect(), {}, {}, {}),
                             std::out_of_range, MessageMatches(ContainsSubstring("busbar section 0")));
    }
    SECTION("two busbar sections on one node")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({1, 1}), IntVect(), IntVect(), {}, {}, {}),
                             std::runtime_error, MessageMatches(ContainsSubstring("two busbar sections")));
    }
    SECTION("a switch to a node that does not exist")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({0}), ivec({0}), ivec({5}), one_breaker, {false}, {false}),
                             std::out_of_range, MessageMatches(ContainsSubstring("switch (side 2) 0")));
    }
    SECTION("a switch joining a node to itself")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({0}), ivec({1}), ivec({1}), one_breaker, {false}, {false}),
                             std::runtime_error, MessageMatches(ContainsSubstring("to itself")));
    }
    SECTION("an open internal connection")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({0}), ivec({0}), ivec({1}),
                                       {SwitchKind::INTERNAL_CONNECTION}, {true}, {false}),
                             std::runtime_error, MessageMatches(ContainsSubstring("always closed")));
    }
    SECTION("switch vectors of different lengths")
    {
        CHECK_THROWS_MATCHES(topo.init(2, ivec({0}), ivec({0}), ivec({1}), one_breaker, {false, true}, {false}),
                             std::runtime_error, MessageMatches(ContainsSubstring("'sw_open'")));
    }
    SECTION("a negative node count")
    {
        CHECK_THROWS_AS(topo.init(-1, IntVect(), IntVect(), IntVect(), {}, {}, {}), std::runtime_error);
    }
    SECTION("the wrong number of names")
    {
        topo = make_sub0_topology();
        CHECK_THROWS_AS(topo.set_sw_names({"only one"}), std::runtime_error);
        CHECK_THROWS_AS(topo.set_bbs_names({"a", "b", "c"}), std::runtime_error);
    }
    SECTION("a poisoned state: unknown switch kind")
    {
        SubstationTopology::StateRes state = make_sub0_topology().get_state();
        std::get<SubstationTopology::SW_KIND>(state)[1] = 99;
        CHECK_THROWS_MATCHES(topo.set_state(state), std::runtime_error,
                             MessageMatches(ContainsSubstring("unknown kind")));
        // and a rejected state leaves the object untouched
        CHECK(topo.nb_nodes() == 0);
    }
    SECTION("a poisoned state: fewer nodes than the switches use")
    {
        SubstationTopology::StateRes state = make_sub0_topology().get_state();
        std::get<SubstationTopology::NB_NODES>(state) = 3;
        CHECK_THROWS_AS(topo.set_state(state), std::out_of_range);
    }
    SECTION("a terminal on a node that does not exist")
    {
        topo = make_sub0_topology();
        CHECK_THROWS_AS(topo.add_terminal(TerminalKind::LOAD, 0, 5), std::out_of_range);
        CHECK(topo.terminals().empty());
    }
}

// ---- SubstationContainer: one topology per substation, numbered grid-wide ----

TEST_CASE("SubstationContainer numbers the substations' local ids grid-wide", "[SubstationContainer][detailed_topology]")
{
    const LSGrid grid = make_grid_with_topology();
    const SubstationContainer & subs = grid.get_substations();
    REQUIRE(subs.has_detailed_topology());
    CHECK(grid.has_detailed_topology());

    CHECK(subs.nb_nodes() == 5 + 3 + 3);
    CHECK(subs.nb_switches() == 8);
    CHECK(subs.nb_busbar_sections() == 4);
    CHECK(subs.first_node(0) == 0);
    CHECK(subs.first_node(1) == 5);
    CHECK(subs.first_node(2) == 8);
    CHECK(subs.first_switch(1) == 5);
    CHECK(subs.first_switch(2) == 7);
    CHECK(subs.first_busbar_section(1) == 2);
    CHECK(subs.first_busbar_section(2) == 3);
    // a grid-wide switch id resolves to (substation, local id)
    CHECK(subs.switch_sub(4) == 0);
    CHECK(subs.switch_local(4) == 4);
    CHECK(subs.switch_sub(5) == 1);
    CHECK(subs.switch_local(5) == 0);
    CHECK(subs.switch_sub(7) == 2);
    CHECK(subs.switch_local(7) == 0);
    CHECK(subs.bbs_sub(3) == 2);
    CHECK(subs.bbs_local(3) == 0);
    CHECK(subs.bbs_sub(1) == 0);
    CHECK(subs.bbs_local(1) == 1);
    CHECK_THROWS_AS(subs.switch_sub(8), std::out_of_range);
    CHECK_THROWS_AS(subs.bbs_sub(-1), std::out_of_range);
    CHECK_THROWS_AS(subs.topology(3), std::out_of_range);

    SECTION("the substation, switch and busbar-section infos carry it")
    {
        const SubstationInfo sub1 = subs[1];
        CHECK(sub1.nb_nodes == 3);
        CHECK(sub1.nb_switches == 2);
        CHECK(sub1.nb_busbar_sections == 1);
        CHECK(sub1.first_node == 5);
        CHECK(sub1.first_switch == 5);
        CHECK(sub1.first_busbar_section == 2);

        const auto switches = grid.get_switches();
        REQUIRE(switches.nb() == 8);
        const SwitchInfo sw5 = switches[5];
        CHECK(sw5.id == 5);
        CHECK(sw5.sub_id == 1);
        CHECK(sw5.local_id == 0);
        CHECK(sw5.node1 == 1);
        CHECK(sw5.node2 == 0);
        CHECK(sw5.kind == SwitchKind::BREAKER);
        CHECK_FALSE(sw5.open);
        CHECK(sw5.retained);
        const SwitchInfo sw4 = switches[4];
        CHECK(sw4.sub_id == 0);
        CHECK(sw4.local_id == 4);
        CHECK(sw4.kind == SwitchKind::INTERNAL_CONNECTION);
        CHECK(sw4.name == "");
        // iterable, one snapshot per switch
        int count = 0;
        for(const SwitchInfo & sw : switches){ CHECK(sw.id == count); ++count; }
        CHECK(count == 8);

        const auto sections = grid.get_busbar_sections();
        REQUIRE(sections.nb() == 4);
        const BusbarSectionInfo bbs1 = sections[1];
        CHECK(bbs1.id == 1);
        CHECK(bbs1.sub_id == 0);
        CHECK(bbs1.local_id == 1);
        CHECK(bbs1.node == 1);
        const BusbarSectionInfo bbs3 = sections[3];
        CHECK(bbs3.sub_id == 2);
        CHECK(bbs3.local_id == 0);
        CHECK(bbs3.node == 0);
    }

    SECTION("names are given in the grid-wide order and land in the right substation")
    {
        LSGrid named = make_grid_with_topology();
        named.set_switch_names({"s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7"});
        named.set_busbar_section_names({"A", "B", "C", "D"});
        CHECK(named.get_switches()[6].name == "s6");
        CHECK(named.get_substation_topology(1).sw_name(1) == "s6");
        CHECK(named.get_busbar_sections()[3].name == "D");
        CHECK(named.get_substation_topology(2).bbs_name(0) == "D");
        CHECK_THROWS_AS(named.set_switch_names({"too", "few"}), std::runtime_error);
        CHECK_THROWS_AS(named.set_busbar_section_names({"too", "few"}), std::runtime_error);
    }
}

TEST_CASE("LSGrid::init_detailed_topology refuses what it cannot number", "[LSGrid][detailed_topology]")
{
    LSGrid grid = make_grid_without_topology();
    const std::vector<int> two_breakers{static_cast<int>(SwitchKind::BREAKER), static_cast<int>(SwitchKind::BREAKER)};

    SECTION("switches not sorted by substation (their position is their grid-wide id)")
    {
        CHECK_THROWS_MATCHES(grid.init_detailed_topology(ivec({2, 2, 2}), ivec({0, 1, 2}), ivec({0, 0, 0}),
                                                         ivec({1, 0}), ivec({0, 0}), ivec({1, 1}),
                                                         two_breakers, {false, false}, {true, true}),
                             std::runtime_error, MessageMatches(ContainsSubstring("sorted by substation")));
        CHECK_FALSE(grid.has_detailed_topology());
    }
    SECTION("a switch in a substation that does not exist")
    {
        CHECK_THROWS_AS(grid.init_detailed_topology(ivec({2, 2, 2}), ivec({0, 1, 2}), ivec({0, 0, 0}),
                                                    ivec({0, 3}), ivec({0, 0}), ivec({1, 1}),
                                                    two_breakers, {false, false}, {true, true}),
                        std::out_of_range);
    }
    SECTION("one node count per substation")
    {
        CHECK_THROWS_AS(grid.init_detailed_topology(ivec({2, 2}), ivec({0, 1, 2}), ivec({0, 0, 0}),
                                                    IntVect(), IntVect(), IntVect(), {}, {}, {}),
                        std::runtime_error);
    }
    SECTION("an inconsistent substation is named")
    {
        // sub 1 has 2 nodes but its section stands on node 2
        CHECK_THROWS_MATCHES(grid.init_detailed_topology(ivec({2, 2, 2}), ivec({0, 1, 2}), ivec({0, 2, 0}),
                                                         IntVect(), IntVect(), IntVect(), {}, {}, {}),
                             std::runtime_error, MessageMatches(ContainsSubstring("substation 1")));
        CHECK_FALSE(grid.has_detailed_topology());
    }
    SECTION("an unknown switch kind")
    {
        CHECK_THROWS_MATCHES(grid.init_detailed_topology(ivec({2, 2, 2}), ivec({0, 1, 2}), ivec({0, 0, 0}),
                                                         ivec({0}), ivec({0}), ivec({1}),
                                                         {42}, {false}, {true}),
                             std::runtime_error, MessageMatches(ContainsSubstring("unknown kind")));
    }
    SECTION("before the substations exist")
    {
        LSGrid empty;
        CHECK_THROWS_AS(empty.init_detailed_topology(IntVect(), IntVect(), IntVect(), IntVect(), IntVect(),
                                                     IntVect(), {}, {}, {}),
                        std::runtime_error);
    }
}

// ---- the node ids ride on the element containers -----------------------------

TEST_CASE("the terminals of each substation are indexed from the elements' node ids", "[LSGrid][detailed_topology]")
{
    const LSGrid grid = make_grid_with_topology();

    // per substation, in container order: loads, gens, ..., line side 1, line side 2
    const std::vector<Terminal> & t0 = grid.get_substation_topology(0).terminals();
    REQUIRE(t0.size() == 2);
    CHECK(same_terminal(t0[0], TerminalKind::GEN, 0, 2));
    CHECK(same_terminal(t0[1], TerminalKind::LINE_1, 0, 4));

    const std::vector<Terminal> & t1 = grid.get_substation_topology(1).terminals();
    REQUIRE(t1.size() == 2);
    CHECK(same_terminal(t1[0], TerminalKind::LINE_1, 1, 2));
    CHECK(same_terminal(t1[1], TerminalKind::LINE_2, 0, 1));

    const std::vector<Terminal> & t2 = grid.get_substation_topology(2).terminals();
    REQUIRE(t2.size() == 2);
    CHECK(same_terminal(t2[0], TerminalKind::LOAD, 0, 1));
    CHECK(same_terminal(t2[1], TerminalKind::LINE_2, 1, 2));

    // and the infos expose them
    CHECK(grid.get_loads()[0].node_id == 1);
    CHECK(grid.get_generators()[0].node_id == 2);
    CHECK(grid.get_lines()[0].node_1_id == 4);
    CHECK(grid.get_lines()[0].node_2_id == 1);
    CHECK(grid.get_lines()[1].node_1_id == 2);
    CHECK(grid.get_lines()[1].node_2_id == 2);

    // consistent: check_grid accepts it
    CHECK_NOTHROW(grid.check_grid());

    SECTION("a terminal the topology does not describe (-1) is simply absent")
    {
        LSGrid partial = make_grid_with_topology();
        partial.set_load_to_node_id(ivec({-1}));
        const std::vector<Terminal> & t2b = partial.get_substation_topology(2).terminals();
        REQUIRE(t2b.size() == 1);
        CHECK(same_terminal(t2b[0], TerminalKind::LINE_2, 1, 2));
        CHECK(partial.get_loads()[0].node_id == -1);
        CHECK_NOTHROW(partial.check_grid());
    }
    SECTION("the substation ids can be set after the node ids: the index follows")
    {
        LSGrid moved = make_grid_with_topology();
        // the load now belongs to substation 1 (node 1 exists there too)
        moved.set_load_to_subid(ivec({1}));
        CHECK(moved.get_substation_topology(2).terminals().size() == 1);
        REQUIRE(moved.get_substation_topology(1).terminals().size() == 3);
        CHECK(same_terminal(moved.get_substation_topology(1).terminals()[0], TerminalKind::LOAD, 0, 1));
    }
    SECTION("a node the substation does not have is refused at once")
    {
        LSGrid bad = make_grid_with_topology();
        CHECK_THROWS_MATCHES(bad.set_load_to_node_id(ivec({3})), std::out_of_range,
                             MessageMatches(ContainsSubstring("node 3")));
        // below -1 is never a node id
        CHECK_THROWS_AS(bad.set_gen_to_node_id(ivec({-2})), std::out_of_range);
    }
    SECTION("a node id without a substation id is refused")
    {
        LSGrid no_sub = make_grid_without_topology();
        no_sub.init_detailed_topology(ivec({2, 2, 2}), ivec({0, 1, 2}), ivec({0, 0, 0}),
                                      IntVect(), IntVect(), IntVect(), {}, {}, {});
        // the static generators never got a substation id in this fixture (none exist),
        // so use a fresh grid whose loads have no subid
        LSGrid raw;
        raw.set_sn_mva(100.);
        raw.set_init_vm_pu(1.0);
        RealVect bus_vn_kv(2);
        bus_vn_kv << 138., 138.;
        raw.init_bus(2, 1, bus_vn_kv, 0, 0);
        RealVect load_p(1), load_q(1);
        load_p << 1.;
        load_q << 0.;
        raw.init_loads(load_p, load_q, ivec({1}));
        raw.init_detailed_topology(ivec({2, 2}), ivec({0, 1}), ivec({0, 0}),
                                   IntVect(), IntVect(), IntVect(), {}, {}, {});
        CHECK_THROWS_MATCHES(raw.set_load_to_node_id(ivec({1})), std::runtime_error,
                             MessageMatches(ContainsSubstring("substation ids first")));
    }
}

TEST_CASE("check_grid catches a node id a poisoned state made inconsistent", "[LSGrid][check_grid][detailed_topology]")
{
    // sub 2 keeps its busbar section on node 0 and has no switch, so shrinking it
    // to one node leaves the SubstationTopology consistent on its own -- but the
    // load stands on node 1, which no longer exists. Only the grid can see that
    // (here as soon as it re-indexes the terminals, before check_grid even runs).
    LSGrid grid = make_grid_without_topology();
    grid.init_detailed_topology(ivec({5, 3, 2}), ivec({0, 0, 1, 2}), ivec({0, 1, 0, 0}),
                                IntVect(), IntVect(), IntVect(), {}, {}, {});
    grid.set_load_to_node_id(ivec({1}));
    REQUIRE_NOTHROW(grid.check_grid());
    LSGrid::StateRes state = grid.get_state();
    std::vector<SubstationTopology::StateRes> & topos = std::get<7>(std::get<LSGrid::SUBSTATION_ID>(state));
    std::get<SubstationTopology::NB_NODES>(topos[2]) = 1;
    LSGrid restored;
    CHECK_THROWS_MATCHES(restored.set_state(state), std::out_of_range,
                         MessageMatches(ContainsSubstring("node 1, out of range [0, 1)")));
}

// ---- state, binary, copy ------------------------------------------------------

TEST_CASE("the detailed topology survives the state, binary and copy round trips", "[LSGrid][detailed_topology][serialization]")
{
    LSGrid grid = make_grid_with_topology();
    grid.set_switch_names({"s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7"});
    grid.set_busbar_section_names({"A", "B", "C", "D"});

    SECTION("get_state / set_state")
    {
        LSGrid::StateRes state = grid.get_state();
        LSGrid restored;
        restored.set_state(state);
        CHECK(same_topology_state(restored, grid));
        CHECK(restored.has_detailed_topology());
        CHECK(restored.get_switches()[6].name == "s6");
        CHECK(restored.get_switches()[0].open);
        CHECK(restored.get_substation_topology(0).terminals().size() == 2);
        CHECK(restored.get_lines()[0].node_1_id == 4);
    }
    SECTION("save_binary / load_binary")
    {
        TempFile file;
        grid.save_binary(file.str());
        const LSGrid loaded = LSGrid::load_binary(file.str());
        CHECK(same_topology_state(loaded, grid));
        CHECK(loaded.get_substations().nb_switches() == 8);
        CHECK(loaded.get_substation_topology(2).terminals().size() == 2);
    }
    SECTION("the substation container alone, through its own binary file")
    {
        TempFile file;
        grid.get_substations().save_binary(file.str());
        const SubstationContainer loaded = SubstationContainer::load_binary(file.str());
        CHECK(loaded.get_state() == grid.get_substations().get_state());
        CHECK(loaded.first_switch(2) == 7);
    }
    SECTION("copy")
    {
        const LSGrid copied = grid.copy();
        CHECK(same_topology_state(copied, grid));
        CHECK(copied.get_busbar_sections()[2].name == "C");
        CHECK(copied.get_substation_topology(1).terminals().size() == 2);
    }
}

// ---- and a grid without one is exactly what it was ----------------------------

TEST_CASE("a grid without detailed topology reports none and round-trips unchanged", "[LSGrid][detailed_topology]")
{
    SECTION("the exotic-elements fixture")
    {
        const LSGrid grid = ls2g_test::make_exotic_elements_grid();
        CHECK_FALSE(grid.has_detailed_topology());
        CHECK_FALSE(grid.get_substations().has_detailed_topology());
        CHECK(grid.get_substations().nb_nodes() == 0);
        CHECK(grid.get_switches().nb() == 0);
        CHECK(grid.get_busbar_sections().nb() == 0);
        const SubstationInfo sub0 = grid.get_substations()[0];
        CHECK(sub0.nb_nodes == 0);
        CHECK(sub0.nb_switches == 0);
        CHECK(sub0.first_node == -1);
        CHECK(grid.get_loads()[0].node_id == -1);
        CHECK(grid.get_lines()[0].node_1_id == -1);
        CHECK_THROWS_AS(grid.get_substation_topology(0), std::runtime_error);

        LSGrid::StateRes state = grid.get_state();
        CHECK(std::get<7>(std::get<LSGrid::SUBSTATION_ID>(state)).empty());
        LSGrid restored;
        restored.set_state(state);
        CHECK(same_topology_state(restored, grid));
        CHECK_FALSE(restored.has_detailed_topology());

        TempFile file;
        grid.save_binary(file.str());
        CHECK(same_topology_state(LSGrid::load_binary(file.str()), grid));
    }
    SECTION("the two-busbar fixture, with substation ids but no node ids")
    {
        const LSGrid grid = make_grid_without_topology();
        CHECK_FALSE(grid.has_detailed_topology());
        CHECK_NOTHROW(grid.check_grid());
        CHECK(grid.get_generators()[0].node_id == -1);
    }
}
