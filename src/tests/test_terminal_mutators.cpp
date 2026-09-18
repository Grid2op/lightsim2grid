// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// The two terminal mutators the switch projection drives the elements through:
// OneSideContainer::set_terminal ("this terminal goes on that bus, or off") and
// TwoSidesContainer::set_terminals (both ends of a branch at once, the global
// status following). On standalone containers, with a SubstationContainer
// counting the buses and a DualAlgoControl collecting the flags -- no LSGrid,
// no powerflow.
//
// What is pinned: the per-bus element counts stay exact through every
// transition, `tell_dimension_changed` is raised exactly when a bus crosses
// 0 <-> 1, a no-op call raises nothing, and the global status of a branch
// follows its ends by the same rule the topology vector uses.

#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "SubstationContainer.hpp"
#include "Utils.hpp"
#include "element_container/LineContainer.hpp"
#include "element_container/LoadContainer.hpp"

using ls2g::CplxVect;
using ls2g::DualAlgoControl;
using ls2g::GridModelBusId;
using ls2g::LineContainer;
using ls2g::LoadContainer;
using ls2g::RealVect;
using ls2g::SubstationContainer;

namespace {

const GridModelBusId OFF(ls2g::GenericContainer::_deactivated_bus_id);

Eigen::VectorXi ivec(std::initializer_list<int> values)
{
    Eigen::VectorXi res(static_cast<Eigen::Index>(values.size()));
    Eigen::Index i = 0;
    for(const int v : values) res(i++) = v;
    return res;
}

// Three substations of two busbars (6 buses: substation s owns buses s and
// s+3), two loads on buses 0 and 2, two lines 0-1 and 1-2. Counts established
// from the elements, flags cleared: from here every flag is raised by the test.
struct Fixture
{
    SubstationContainer subs;
    LoadContainer loads;
    LineContainer lines;
    DualAlgoControl ctrl;

    Fixture()
    {
        RealVect vn = RealVect::Constant(6, 138.);
        subs.init_bus(3, 2, vn);
        RealVect p(2), q(2);
        p << 50., 20.;
        q << 10., 5.;
        loads.init(p, q, ivec({0, 2}));
        RealVect r(2), x(2);
        r << 0.01, 0.01;
        x << 0.1, 0.1;
        lines.init(r, x, CplxVect::Zero(2), ivec({0, 1}), ivec({1, 2}));
        recount();
        clear_flags();
    }

    void recount()
    {
        // armed first: a contribution to disarmed counts is ignored (see
        // SubstationContainer::bus_gained_element), the order
        // LSGrid::recompute_bus_element_counts uses
        subs.reset_bus_element_counts();
        subs.mark_bus_counts_ready();
        bool unused = false;
        for(int i = 0; i < loads.nb(); ++i) loads.contribute_to_buses(i, subs, +1, unused);
        for(int i = 0; i < lines.nb(); ++i) lines.contribute_to_buses(i, subs, +1, unused);
        subs.recount_connected_buses();
    }
    void clear_flags()
    {
        ctrl.ac_algo_controler().tell_none_changed();
        ctrl.dc_algo_controler().tell_none_changed();
    }
    bool dimension_changed() const
    {
        return ctrl.ac_algo_controler().has_dimension_changed() || ctrl.dc_algo_controler().has_dimension_changed();
    }
    bool nothing_changed() const
    {
        return ctrl.ac_algo_controler().nothing_changed() && ctrl.dc_algo_controler().nothing_changed();
    }
    // the incremental counts against a fresh count from the elements
    void require_counts_exact(const char * what)
    {
        INFO("after: " << what);
        const std::vector<std::size_t> incremental = subs.get_nb_elements_per_bus();
        CHECK(subs.connected_bus_count_is_exact());
        recount();
        const std::vector<std::size_t> & from_scratch = subs.get_nb_elements_per_bus();
        REQUIRE(incremental.size() == from_scratch.size());
        for(std::size_t b = 0; b < from_scratch.size(); ++b){
            INFO("bus " << b);
            CHECK(incremental[b] == from_scratch[b]);
        }
    }
    std::size_t count(int bus) const { return subs.get_nb_elements_per_bus()[static_cast<std::size_t>(bus)]; }
};

}  // anonymous namespace

// ---- set_terminal -------------------------------------------------------------

TEST_CASE("set_terminal moves, drops and restores a terminal with exact counts", "[set_terminal][bus_count]")
{
    Fixture f;
    REQUIRE(f.count(0) == 2);  // load 0 + line 0's first end
    REQUIRE(f.count(3) == 0);

    SECTION("a move to an empty bus of the same substation makes that bus exist")
    {
        CHECK(f.loads.set_terminal(0, GridModelBusId(3), f.ctrl, f.subs));
        CHECK(f.loads.get_bus(0) == GridModelBusId(3));
        CHECK(f.count(0) == 1);
        CHECK(f.count(3) == 1);
        CHECK(f.dimension_changed());  // bus 3: 0 -> 1
        f.require_counts_exact("move to an empty bus");
    }
    SECTION("a move to a bus that already exists changes no dimension")
    {
        CHECK(f.loads.set_terminal(0, GridModelBusId(1), f.ctrl, f.subs));
        CHECK(f.count(0) == 1);
        CHECK(f.count(1) == 3);
        CHECK_FALSE(f.dimension_changed());
        CHECK_FALSE(f.nothing_changed());  // but the Sbus moved
        f.require_counts_exact("move to a live bus");
    }
    SECTION("off, then back on the same bus")
    {
        CHECK(f.loads.set_terminal(1, OFF, f.ctrl, f.subs));
        CHECK_FALSE(f.loads.get_status(1));
        CHECK(f.count(2) == 1);  // line 1's second end stays
        CHECK_FALSE(f.dimension_changed());
        f.require_counts_exact("off");
        f.clear_flags();
        CHECK(f.loads.set_terminal(1, GridModelBusId(2), f.ctrl, f.subs));
        CHECK(f.loads.get_status(1));
        CHECK(f.loads.get_bus(1) == GridModelBusId(2));
        CHECK(f.count(2) == 2);
        f.require_counts_exact("back on");
    }
    SECTION("off from a bus it held alone, then back on another: two crossings")
    {
        REQUIRE(f.loads.set_terminal(0, GridModelBusId(3), f.ctrl, f.subs));
        f.clear_flags();
        CHECK(f.loads.set_terminal(0, OFF, f.ctrl, f.subs));
        CHECK(f.count(3) == 0);
        CHECK(f.dimension_changed());  // bus 3: 1 -> 0
        f.require_counts_exact("off from a lone bus");
        f.clear_flags();
        // reactivate AND move, in one bracket: the last bus (3) is never re-held
        CHECK(f.loads.set_terminal(0, GridModelBusId(4), f.ctrl, f.subs));
        CHECK(f.count(3) == 0);
        CHECK(f.count(4) == 1);
        CHECK(f.dimension_changed());  // bus 4: 0 -> 1
        f.require_counts_exact("back on elsewhere");
    }
}

TEST_CASE("set_terminal on a terminal already there is a no-op and raises nothing", "[set_terminal][bus_count]")
{
    Fixture f;
    CHECK_FALSE(f.loads.set_terminal(0, GridModelBusId(0), f.ctrl, f.subs));
    CHECK(f.nothing_changed());
    REQUIRE(f.loads.set_terminal(1, OFF, f.ctrl, f.subs));
    f.clear_flags();
    CHECK_FALSE(f.loads.set_terminal(1, OFF, f.ctrl, f.subs));
    CHECK(f.nothing_changed());
    f.require_counts_exact("no-ops");
}

TEST_CASE("set_terminal refuses a bad id before touching the counts", "[set_terminal][bus_count]")
{
    Fixture f;
    CHECK_THROWS_AS(f.loads.set_terminal(2, GridModelBusId(0), f.ctrl, f.subs), std::out_of_range);
    CHECK_THROWS_AS(f.loads.set_terminal(-1, GridModelBusId(0), f.ctrl, f.subs), std::out_of_range);
    CHECK_THROWS_AS(f.loads.set_terminal(0, GridModelBusId(6), f.ctrl, f.subs), std::out_of_range);
    CHECK_THROWS_AS(f.loads.set_terminal(0, GridModelBusId(-2), f.ctrl, f.subs), std::out_of_range);
    CHECK(f.nothing_changed());
    CHECK(f.count(0) == 2);
    CHECK(f.loads.get_bus(0) == GridModelBusId(0));
    f.require_counts_exact("refused calls");
}

// ---- set_terminals ------------------------------------------------------------

TEST_CASE("set_terminals moves both ends of a branch in one bracket", "[set_terminals][bus_count]")
{
    Fixture f;
    REQUIRE(f.lines.get_synch_status_both_side());
    // line 0 goes from 0-1 to 3-4: two empty buses appear, two live ones thin out
    CHECK(f.lines.set_terminals(0, GridModelBusId(3), GridModelBusId(4), f.ctrl, f.subs));
    CHECK(f.lines.get_bus_side_1(0) == GridModelBusId(3));
    CHECK(f.lines.get_bus_side_2(0) == GridModelBusId(4));
    CHECK(f.lines.get_status_global()[0]);
    CHECK(f.count(0) == 1);
    CHECK(f.count(1) == 1);
    CHECK(f.count(3) == 1);
    CHECK(f.count(4) == 1);
    CHECK(f.dimension_changed());
    f.require_counts_exact("both ends moved");
    f.clear_flags();
    // and back: the same call again is a no-op
    CHECK(f.lines.set_terminals(0, GridModelBusId(0), GridModelBusId(1), f.ctrl, f.subs));
    f.clear_flags();
    CHECK_FALSE(f.lines.set_terminals(0, GridModelBusId(0), GridModelBusId(1), f.ctrl, f.subs));
    CHECK(f.nothing_changed());
    f.require_counts_exact("back and no-op");
}

TEST_CASE("with synched ends, one end off takes the whole branch down", "[set_terminals][bus_count]")
{
    Fixture f;
    CHECK(f.lines.set_terminals(1, GridModelBusId(1), OFF, f.ctrl, f.subs));
    CHECK_FALSE(f.lines.get_status_global()[1]);
    CHECK_FALSE(f.lines.get_connected_side_1(1));
    CHECK_FALSE(f.lines.get_connected_side_2(1));
    CHECK(f.count(1) == 1);  // line 0's end 2 only: line 1's end 1 is gone
    CHECK(f.count(2) == 1);  // the load alone
    CHECK_FALSE(f.dimension_changed());
    f.require_counts_exact("synched branch down");
    f.clear_flags();
    // both ends on again: global back on, both ends where they are told
    CHECK(f.lines.set_terminals(1, GridModelBusId(4), GridModelBusId(2), f.ctrl, f.subs));
    CHECK(f.lines.get_status_global()[1]);
    CHECK(f.lines.get_bus_side_1(1) == GridModelBusId(4));
    CHECK(f.lines.get_bus_side_2(1) == GridModelBusId(2));
    CHECK(f.count(4) == 1);
    CHECK(f.dimension_changed());  // bus 4: 0 -> 1
    f.require_counts_exact("synched branch back");
    f.clear_flags();
    // both off: the same as one off, and repeating it is a no-op
    CHECK(f.lines.set_terminals(1, OFF, OFF, f.ctrl, f.subs));
    CHECK(f.count(4) == 0);
    CHECK(f.dimension_changed());
    f.clear_flags();
    CHECK_FALSE(f.lines.set_terminals(1, OFF, OFF, f.ctrl, f.subs));
    CHECK_FALSE(f.lines.set_terminals(1, GridModelBusId(4), OFF, f.ctrl, f.subs));  // also "both off"
    CHECK(f.nothing_changed());
    f.require_counts_exact("both off");
}

TEST_CASE("without synched ends, a branch can be half-open and is on while one end is", "[set_terminals][bus_count]")
{
    Fixture f;
    f.lines.set_synch_status_both_side(false);
    CHECK(f.lines.set_terminals(1, GridModelBusId(1), OFF, f.ctrl, f.subs));
    CHECK(f.lines.get_status_global()[1]);
    CHECK(f.lines.get_connected_side_1(1));
    CHECK_FALSE(f.lines.get_connected_side_2(1));
    CHECK(f.count(1) == 2);  // line 0's end 2 and line 1's end 1
    CHECK(f.count(2) == 1);  // the load alone
    f.require_counts_exact("half-open");
    f.clear_flags();
    // the other end alone
    CHECK(f.lines.set_terminals(1, OFF, GridModelBusId(5), f.ctrl, f.subs));
    CHECK(f.lines.get_status_global()[1]);
    CHECK_FALSE(f.lines.get_connected_side_1(1));
    CHECK(f.lines.get_connected_side_2(1));
    CHECK(f.count(1) == 1);
    CHECK(f.count(5) == 1);
    CHECK(f.dimension_changed());
    f.require_counts_exact("other end alone");
    f.clear_flags();
    // both off: now the branch is off
    CHECK(f.lines.set_terminals(1, OFF, OFF, f.ctrl, f.subs));
    CHECK_FALSE(f.lines.get_status_global()[1]);
    CHECK(f.count(5) == 0);
    f.require_counts_exact("both off, unsynched");
    f.clear_flags();
    // and back on with both ends
    CHECK(f.lines.set_terminals(1, GridModelBusId(1), GridModelBusId(2), f.ctrl, f.subs));
    CHECK(f.lines.get_status_global()[1]);
    CHECK(f.count(1) == 2);
    CHECK(f.count(2) == 2);
    f.require_counts_exact("back on, unsynched");
}

TEST_CASE("with the global gate ignored, the branch stays on whatever its ends do", "[set_terminals][bus_count]")
{
    Fixture f;
    f.lines.set_synch_status_both_side(false);
    f.lines.set_ignore_status_global(true);
    CHECK(f.lines.set_terminals(1, OFF, OFF, f.ctrl, f.subs));
    CHECK(f.lines.get_status_global()[1]);
    CHECK_FALSE(f.lines.get_connected_side_1(1));
    CHECK_FALSE(f.lines.get_connected_side_2(1));
    f.require_counts_exact("ends off, gate ignored");
    f.clear_flags();
    CHECK_FALSE(f.lines.set_terminals(1, OFF, OFF, f.ctrl, f.subs));
    CHECK(f.nothing_changed());
}

TEST_CASE("set_terminals refuses a bad id before touching the counts", "[set_terminals][bus_count]")
{
    Fixture f;
    CHECK_THROWS_AS(f.lines.set_terminals(2, GridModelBusId(0), GridModelBusId(1), f.ctrl, f.subs), std::out_of_range);
    CHECK_THROWS_AS(f.lines.set_terminals(0, GridModelBusId(0), GridModelBusId(6), f.ctrl, f.subs), std::out_of_range);
    CHECK_THROWS_AS(f.lines.set_terminals(0, GridModelBusId(-3), GridModelBusId(1), f.ctrl, f.subs), std::out_of_range);
    CHECK(f.nothing_changed());
    CHECK(f.lines.get_bus_side_2(0) == GridModelBusId(1));
    f.require_counts_exact("refused calls");
}
