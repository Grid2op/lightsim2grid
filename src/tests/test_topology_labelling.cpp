// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Labelling the components of a substation's closed-switch graph as electrical
// buses (SubstationTopology::label): the validity rule -- pypowsybl's, checked
// row by row against what `Network.get_buses()` does on pypowsybl 1.15 -- the
// numbering, the capacity check, and the grid-wide view LSGrid gives of it.
//
// Every case here is a pure function of the declared switches and terminals:
// nothing moves an element. That push is the projection, tested separately.

#include <stdexcept>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "LSGrid.hpp"
#include "SubstationTopology.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using ls2g::IntVect;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::SubstationTopology;
using ls2g::SwitchKind;
using ls2g::TerminalKind;

namespace {

IntVect ivec(std::initializer_list<int> values)
{
    IntVect res(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for(const int v : values) res(i++) = v;
    return res;
}

constexpr int NO_LIMIT = 1000;  // a capacity no test here reaches

// One busbar section on node 0 and `nb_extra` bare nodes, all joined to it by
// closed breakers (node k <-> node 0 is switch k-1). Terminals are added by the
// test, on whichever node it wants.
SubstationTopology make_star(int nb_extra)
{
    IntVect node1(nb_extra), node2(nb_extra);
    std::vector<SwitchKind> kinds(static_cast<std::size_t>(nb_extra), SwitchKind::BREAKER);
    std::vector<bool> open(static_cast<std::size_t>(nb_extra), false);
    std::vector<bool> retained(static_cast<std::size_t>(nb_extra), true);
    for(int k = 0; k < nb_extra; ++k){
        node1(k) = k + 1;
        node2(k) = 0;
    }
    SubstationTopology topo;
    topo.init(nb_extra + 1, ivec({0}), node1, node2, kinds, open, retained);
    return topo;
}

// Two nodes joined by one closed breaker and NO busbar section: the component
// is whatever the two terminals put on it make of it.
SubstationTopology make_pair_no_bbs()
{
    SubstationTopology topo;
    topo.init(2, IntVect(), ivec({0}), ivec({1}), {SwitchKind::BREAKER}, {false}, {true});
    return topo;
}

}  // anonymous namespace

// ---- the validity rule, row by row -------------------------------------------

TEST_CASE("an isolated busbar section is not a bus", "[SubstationTopology][label]")
{
    SubstationTopology topo = make_star(0);
    topo.label(NO_LIMIT);
    CHECK(topo.labels_ready());
    CHECK(topo.nb_buses() == 0);
    CHECK(topo.node_bus(0) == -1);
    CHECK(topo.bbs_bus(0) == -1);
}

TEST_CASE("a busbar section with one feeder is a bus", "[SubstationTopology][label]")
{
    SECTION("a load")
    {
        SubstationTopology topo = make_star(1);
        topo.add_terminal(TerminalKind::LOAD, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
        CHECK(topo.node_bus(0) == 1);
        CHECK(topo.node_bus(1) == 1);
    }
    SECTION("a line end")
    {
        SubstationTopology topo = make_star(1);
        topo.add_terminal(TerminalKind::LINE_2, 3, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
        CHECK(topo.bbs_bus(0) == 1);
    }
    SECTION("but not when the feeder's breaker is open")
    {
        SubstationTopology topo = make_star(1);
        topo.add_terminal(TerminalKind::LOAD, 0, 1);
        REQUIRE(topo.set_open(0, true));
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
        CHECK(topo.node_bus(0) == -1);
        CHECK(topo.node_bus(1) == -1);
    }
}

TEST_CASE("without a busbar section, a branch end with another feeder is a bus", "[SubstationTopology][label]")
{
    SECTION("line end + line end")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::LINE_1, 0, 0);
        topo.add_terminal(TerminalKind::LINE_2, 1, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
        CHECK(topo.node_bus(0) == 1);
        CHECK(topo.node_bus(1) == 1);
    }
    SECTION("line end + load")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::LINE_1, 0, 0);
        topo.add_terminal(TerminalKind::LOAD, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
    }
    SECTION("transformer end + shunt")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::TRAFO_2, 0, 0);
        topo.add_terminal(TerminalKind::SHUNT, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
    }
    SECTION("an HVDC converter station counts as a branch: station + load")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::HVDC_1, 0, 0);
        topo.add_terminal(TerminalKind::LOAD, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 1);
    }
}

TEST_CASE("without a busbar section, a lone branch end or injections only are not a bus", "[SubstationTopology][label]")
{
    SECTION("a single line end (its breaker open): pypowsybl disconnects it")
    {
        SubstationTopology topo = make_star(1);
        topo.add_terminal(TerminalKind::LINE_1, 0, 1);
        topo.set_open(0, true);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
        CHECK(topo.node_bus(1) == -1);
    }
    SECTION("load + generator")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::LOAD, 0, 0);
        topo.add_terminal(TerminalKind::GEN, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
    }
    SECTION("load + load")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::LOAD, 0, 0);
        topo.add_terminal(TerminalKind::LOAD, 1, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
    }
    SECTION("SVC + load")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::SVC, 0, 0);
        topo.add_terminal(TerminalKind::LOAD, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
    }
    SECTION("storage + static generator")
    {
        SubstationTopology topo = make_pair_no_bbs();
        topo.add_terminal(TerminalKind::STORAGE, 0, 0);
        topo.add_terminal(TerminalKind::SGEN, 0, 1);
        topo.label(NO_LIMIT);
        CHECK(topo.nb_buses() == 0);
    }
}

// ---- the graph: what joins, what separates ------------------------------------

TEST_CASE("an internal connection joins like a closed switch and cannot be opened", "[SubstationTopology][label]")
{
    // bbs on node 0; node 1 reachable through an internal connection; node 2
    // through a closed disconnector from node 1
    SubstationTopology topo;
    topo.init(3, ivec({0}), ivec({1, 2}), ivec({0, 1}),
              {SwitchKind::INTERNAL_CONNECTION, SwitchKind::DISCONNECTOR}, {false, false}, {false, false});
    topo.add_terminal(TerminalKind::LOAD, 0, 2);
    topo.label(NO_LIMIT);
    CHECK(topo.nb_buses() == 1);
    CHECK(topo.node_bus(2) == 1);
    CHECK_THROWS_MATCHES(topo.set_open(0, true), std::runtime_error,
                         MessageMatches(ContainsSubstring("internal connection")));
    // and refusing changed nothing
    CHECK_FALSE(topo.is_open(0));
    CHECK(topo.labels_ready());
}

TEST_CASE("an open coupler splits two sections into two buses, a closed one merges them", "[SubstationTopology][label]")
{
    // sections A (node 0) and B (node 1), coupler 0-1 (switch 0), a load on node 2
    // joined to A (switch 1), a line end on node 3 joined to B (switch 2)
    SubstationTopology topo;
    topo.init(4, ivec({0, 1}), ivec({0, 2, 3}), ivec({1, 0, 1}),
              {SwitchKind::BREAKER, SwitchKind::BREAKER, SwitchKind::BREAKER},
              {true, false, false}, {true, true, true});
    topo.add_terminal(TerminalKind::LOAD, 0, 2);
    topo.add_terminal(TerminalKind::LINE_1, 0, 3);

    topo.label(NO_LIMIT);
    CHECK(topo.nb_buses() == 2);
    CHECK(topo.bbs_bus(0) == 1);
    CHECK(topo.bbs_bus(1) == 2);
    CHECK(topo.node_bus(2) == 1);
    CHECK(topo.node_bus(3) == 2);

    REQUIRE(topo.set_open(0, false));
    CHECK_FALSE(topo.labels_ready());  // stale until relabelled
    CHECK_THROWS_AS(topo.node_bus(0), std::runtime_error);
    topo.label(NO_LIMIT);
    CHECK(topo.nb_buses() == 1);
    for(int node = 0; node < 4; ++node) CHECK(topo.node_bus(node) == 1);

    // setting a switch to the position it already has changes nothing
    CHECK_FALSE(topo.set_open(0, false));
    CHECK(topo.labels_ready());
}

// ---- the numbering ------------------------------------------------------------

TEST_CASE("buses are numbered busbar sections first, in section order, then by lowest terminal node", "[SubstationTopology][label]")
{
    // sections declared in the order node 5 (section 0), node 0 (section 1); a
    // section-less line-line component on nodes 2-3; the section on node 5 holds a
    // load on node 6, the one on node 0 a load on node 1. Node 4 is bare.
    SubstationTopology topo;
    topo.init(7, ivec({5, 0}),
              ivec({6, 1, 2}), ivec({5, 0, 3}),
              {SwitchKind::BREAKER, SwitchKind::BREAKER, SwitchKind::BREAKER},
              {false, false, false}, {true, true, true});
    topo.add_terminal(TerminalKind::LOAD, 0, 6);
    topo.add_terminal(TerminalKind::LOAD, 1, 1);
    topo.add_terminal(TerminalKind::LINE_1, 0, 2);
    topo.add_terminal(TerminalKind::LINE_2, 1, 3);
    topo.label(NO_LIMIT);
    CHECK(topo.nb_buses() == 3);
    // section 0 (node 5) is bus 1 although its node id is larger
    CHECK(topo.bbs_bus(0) == 1);
    CHECK(topo.node_bus(6) == 1);
    CHECK(topo.bbs_bus(1) == 2);
    CHECK(topo.node_bus(1) == 2);
    // the section-less component comes last
    CHECK(topo.node_bus(2) == 3);
    CHECK(topo.node_bus(3) == 3);
    CHECK(topo.node_bus(4) == -1);

    SECTION("the labels are a function of the switch positions only")
    {
        // reach the same positions along another history: open, relabel, close
        SubstationTopology other = topo;
        other.set_open(1, true);
        other.label(NO_LIMIT);
        CHECK(other.nb_buses() == 2);
        other.set_open(1, false);
        other.label(NO_LIMIT);
        CHECK(other.node_bus() == topo.node_bus());
        CHECK(other.nb_buses() == topo.nb_buses());
    }
}

TEST_CASE("more buses than the layout holds is refused before the labels move", "[SubstationTopology][label]")
{
    // two sections, each with a load, coupler open: two buses
    SubstationTopology topo;
    topo.init(4, ivec({0, 1}), ivec({0, 2, 3}), ivec({1, 0, 1}),
              {SwitchKind::BREAKER, SwitchKind::BREAKER, SwitchKind::BREAKER},
              {false, false, false}, {true, true, true});
    topo.add_terminal(TerminalKind::LOAD, 0, 2);
    topo.add_terminal(TerminalKind::LOAD, 1, 3);
    topo.label(2);
    REQUIRE(topo.nb_buses() == 1);  // coupler closed: one bus

    topo.set_open(0, true);
    CHECK_THROWS_MATCHES(topo.label(1, 7), std::runtime_error,
                         MessageMatches(ContainsSubstring("substation 7") && ContainsSubstring("2 electrical buses")));
    // the labels are the old ones, and say so
    CHECK_FALSE(topo.labels_ready());
    CHECK(topo.nb_buses() == 1);
    // with the room, it works
    topo.label(2);
    CHECK(topo.nb_buses() == 2);
}

TEST_CASE("an empty substation labels to nothing", "[SubstationTopology][label]")
{
    SubstationTopology topo;
    topo.label(1);
    CHECK(topo.labels_ready());
    CHECK(topo.nb_buses() == 0);
    CHECK(topo.node_bus().size() == 0);
}

// ---- the grid-wide view -------------------------------------------------------

namespace {

// The fixture of test_substation_topology.cpp: three substations of two busbars
// each, with the detailed topology declared and every terminal placed.
//   sub 0: sections A (node 0) and B (node 1), coupler open; the generator (node 2)
//          and the line-0 bay (node 3, terminal node 4 behind an internal
//          connection) both reach A -> {0,2,3,4} is bus 1, B alone is nothing
//   sub 1: one section, both line ends behind closed breakers -> one bus
//   sub 2: one section with line 1's end -> bus 1; the load on node 1 with no
//          switch at all -> nothing
LSGrid make_labelled_grid()
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
    grid.init_powerlines(branch_r, branch_x, ls2g::CplxVect::Zero(2), ivec({0, 1}), ivec({1, 2}));
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

}  // anonymous namespace

TEST_CASE("LSGrid labels every substation as soon as the terminals are known", "[LSGrid][label]")
{
    const LSGrid grid = make_labelled_grid();
    for(int sub_id = 0; sub_id < 3; ++sub_id) CHECK(grid.get_substation_topology(sub_id).labels_ready());
    CHECK(grid.get_substation_topology(0).nb_buses() == 1);
    CHECK(grid.get_substation_topology(1).nb_buses() == 1);
    CHECK(grid.get_substation_topology(2).nb_buses() == 1);

    // by grid-wide node id, in the grid's bus numbering (bus = sub + (local-1) * n_sub)
    const IntVect node_bus = grid.get_node_bus();
    REQUIRE(node_bus.size() == 11);
    CHECK(node_bus == ivec({0, -1, 0, 0, 0, 1, 1, 1, 2, -1, 2}));

    // the busbar sections know their bus
    const auto sections = grid.get_busbar_sections();
    CHECK(sections[0].connected);
    CHECK(sections[0].bus_id == 0);
    CHECK_FALSE(sections[1].connected);
    CHECK(sections[1].bus_id == -1);
    CHECK(sections[2].bus_id == 1);
    CHECK(sections[3].bus_id == 2);

    SECTION("moving a terminal relabels: the load reaches its section, a second bus appears elsewhere")
    {
        LSGrid moved = make_labelled_grid();
        // the load now stands on section B of substation 0 (node 1): B becomes bus 2 there
        moved.set_load_to_subid(ivec({0}));
        moved.set_load_to_node_id(ivec({1}));
        CHECK(moved.get_substation_topology(0).nb_buses() == 2);
        CHECK(moved.get_node_bus()(1) == 3);  // sub 0, local bus 2 -> 0 + 1 * 3
        // and substation 2 lost its load: its section keeps line 1's end, still a bus
        CHECK(moved.get_substation_topology(2).nb_buses() == 1);
        CHECK(moved.get_node_bus()(9) == -1);
    }
    SECTION("the labels come back with the state")
    {
        LSGrid::StateRes state = grid.get_state();
        LSGrid restored;
        restored.set_state(state);
        CHECK(restored.get_node_bus() == node_bus);
        CHECK(restored.get_busbar_sections()[3].bus_id == 2);
        const LSGrid copied = grid.copy();
        CHECK(copied.get_node_bus() == node_bus);
    }
    SECTION("a substation the layout cannot hold is refused when its terminals are indexed")
    {
        // one busbar per substation: substation 0 with the load on section B needs two
        LSGrid tight;
        tight.set_sn_mva(100.);
        tight.set_init_vm_pu(1.0);
        RealVect vn(2);
        vn << 138., 138.;
        tight.init_bus(2, 1, vn, 0, 0);
        RealVect p(2), q(2);
        p << 1., 1.;
        q << 0., 0.;
        tight.init_loads(p, q, ivec({0, 1}));
        tight.set_load_to_subid(ivec({0, 0}));
        tight.init_detailed_topology(ivec({4, 0}), ivec({0, 0}), ivec({0, 1}),
                                     ivec({0, 0}), ivec({2, 3}), ivec({0, 1}),
                                     {static_cast<int>(SwitchKind::BREAKER), static_cast<int>(SwitchKind::BREAKER)},
                                     {false, false}, {true, true});
        CHECK_THROWS_MATCHES(tight.set_load_to_node_id(ivec({2, 3})), std::runtime_error,
                             MessageMatches(ContainsSubstring("2 electrical buses") && ContainsSubstring("at most 1")));
    }
}

TEST_CASE("get_node_bus is empty without detailed topology", "[LSGrid][label]")
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);
    RealVect vn(2);
    vn << 138., 138.;
    grid.init_bus(2, 1, vn, 0, 0);
    CHECK(grid.get_node_bus().size() == 0);
}
