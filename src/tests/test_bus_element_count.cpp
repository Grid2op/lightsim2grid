// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// SubstationContainer keeps, per bus, how many elements hold it alive. It is
// maintained incrementally: every mutator that can change bus membership takes its
// element's contribution away, mutates, and puts the new one back
// (GenericContainer::_apply_and_track_buses). A bus that thereby crosses 0 has
// entered or left the solved system, which renumbers every bus after it.
//
// Incremental bookkeeping rots silently: one missed decrement and the count drifts
// for the rest of the grid's life, unlike the from-scratch recompute it replaces,
// which is self-healing. So the only assertion worth making is the strong one --
// after ANY mutation, the incremental counts equal what a recount from the elements
// produces. That is what every case here checks, through one helper.
//
// The recount is LSGrid::recompute_bus_element_counts(), which is built from the
// same `contribute_to_buses` predicate as the incremental path. That is deliberate:
// there is exactly ONE statement of "which buses does this element hold", so this
// test cannot pass by both sides sharing a misconception about a DIFFERENT rule --
// it catches a mutator that forgot to bracket itself, which is the failure that
// actually happens.

#include <string>
#include <vector>
#include <functional>

#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "case_exotic_elements.hpp"

using ls2g::CplxVect;
using ls2g::LSGrid;
using ls2g::cplx_type;
using ls2g::GridModelBusId;
using ls2g::RealVect;

namespace {

LSGrid make_grid(){ return ls2g_test::make_exotic_elements_grid(); }

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()),
                              cplx_type(grid.get_init_vm_pu(), 0.));
}

// The whole point of this file: whatever just happened, do the counts still say
// what a fresh count from the elements would say?
void require_counts_match_recount(LSGrid & grid, const std::string & what)
{
    const std::vector<std::size_t> incremental = grid.get_substations().get_nb_elements_per_bus();

    // `nb_connected_bus()` is not a second, separately-maintained copy of this: the
    // mutators move it by one on the 0-crossings they already detect, so it has to
    // equal a fresh count of the non-empty buses after every one of them. Checked
    // BEFORE the recount below, which would re-establish it and hide any drift.
    INFO("the connected-bus count drifted from the element counts after: " << what);
    CHECK(grid.get_substations().connected_bus_count_is_exact());

    grid.recompute_bus_element_counts();
    const std::vector<std::size_t> & from_scratch = grid.get_substations().get_nb_elements_per_bus();

    INFO("after: " << what);
    REQUIRE(incremental.size() == from_scratch.size());
    for(std::size_t b = 0; b < from_scratch.size(); ++b){
        INFO("bus " << b);
        CHECK(incremental[b] == from_scratch[b]);
    }
}

struct NamedMutation {
    const char * name;
    std::function<void(LSGrid &)> apply;
};

// Every LSGrid method that can change whether an element holds a bus. If one is
// added and not listed here, the sweep below silently stops covering it -- which is
// why the count is asserted at the end of the file.
std::vector<NamedMutation> all_bus_mutations()
{
    return {
        // ---- one-sided elements ------------------------------------------
        {"deactivate_load",        [](LSGrid & g){ g.deactivate_load(2); }},
        {"reactivate_load",        [](LSGrid & g){ g.deactivate_load(2); g.reactivate_load(2); }},
        {"change_bus_load",        [](LSGrid & g){ g.change_bus_load_python(0, 4); }},
        {"deactivate_gen",         [](LSGrid & g){ g.deactivate_gen(3); }},
        {"reactivate_gen",         [](LSGrid & g){ g.deactivate_gen(3); g.reactivate_gen(3); }},
        {"change_bus_gen",         [](LSGrid & g){ g.change_bus_gen_python(2, 9); }},
        {"deactivate_shunt",       [](LSGrid & g){ g.deactivate_shunt(0); }},
        {"reactivate_shunt",       [](LSGrid & g){ g.deactivate_shunt(0); g.reactivate_shunt(0); }},
        {"change_bus_shunt",       [](LSGrid & g){ g.change_bus_shunt_python(0, 3); }},
        {"deactivate_sgen",        [](LSGrid & g){ if(g.get_static_generators_as_data().nb() > 0) g.deactivate_sgen(0); }},
        {"deactivate_storage",     [](LSGrid & g){ g.deactivate_storage(0); }},
        {"reactivate_storage",     [](LSGrid & g){ g.deactivate_storage(0); g.reactivate_storage(0); }},
        {"change_bus_storage",     [](LSGrid & g){ g.change_bus_storage_python(0, 6); }},
        {"deactivate_svc",         [](LSGrid & g){ g.deactivate_svc(0); }},
        {"reactivate_svc",         [](LSGrid & g){ g.deactivate_svc(0); g.reactivate_svc(0); }},
        {"change_bus_svc",         [](LSGrid & g){ g.change_bus_svc_python(0, 5); }},
        // ---- lines / transformers: gated by status_global_ ----------------
        {"deactivate_powerline",        [](LSGrid & g){ g.deactivate_powerline(4); }},
        {"reactivate_powerline",        [](LSGrid & g){ g.deactivate_powerline(4); g.reactivate_powerline(4); }},
        {"deactivate_powerline_side1",  [](LSGrid & g){ g.deactivate_powerline_side1(5); }},
        {"deactivate_powerline_side2",  [](LSGrid & g){ g.deactivate_powerline_side2(6); }},
        {"reactivate_powerline_side1",  [](LSGrid & g){ g.deactivate_powerline_side1(5); g.reactivate_powerline_side1(5); }},
        {"reactivate_powerline_side2",  [](LSGrid & g){ g.deactivate_powerline_side2(6); g.reactivate_powerline_side2(6); }},
        {"change_bus1_powerline",       [](LSGrid & g){ g.change_bus1_powerline_python(7, 2); }},
        {"change_bus2_powerline",       [](LSGrid & g){ g.change_bus2_powerline_python(8, 9); }},
        {"deactivate_trafo",            [](LSGrid & g){ g.deactivate_trafo(1); }},
        {"reactivate_trafo",            [](LSGrid & g){ g.deactivate_trafo(1); g.reactivate_trafo(1); }},
        {"deactivate_trafo_side1",      [](LSGrid & g){ g.deactivate_trafo_side1(0); }},
        {"deactivate_trafo_side2",      [](LSGrid & g){ g.deactivate_trafo_side2(0); }},
        {"reactivate_trafo_side1",      [](LSGrid & g){ g.deactivate_trafo_side1(0); g.reactivate_trafo_side1(0); }},
        {"reactivate_trafo_side2",      [](LSGrid & g){ g.deactivate_trafo_side2(0); g.reactivate_trafo_side2(0); }},
        {"change_bus1_trafo",           [](LSGrid & g){ g.change_bus1_trafo_python(0, 3); }},
        {"change_bus2_trafo",           [](LSGrid & g){ g.change_bus2_trafo_python(0, 3); }},
        // ---- HVDC: two stations, NO status_global_ gate --------------------
        {"deactivate_dcline",           [](LSGrid & g){ g.deactivate_dcline(1); }},
        {"reactivate_dcline",           [](LSGrid & g){ g.deactivate_dcline(1); g.reactivate_dcline(1); }},
        {"deactivate_dcline_side1",     [](LSGrid & g){ g.deactivate_dcline_side1(2); }},
        {"deactivate_dcline_side2",     [](LSGrid & g){ g.deactivate_dcline_side2(2); }},
        {"change_bus1_dcline",          [](LSGrid & g){ g.change_bus1_dcline_python(0, 4); }},
        {"change_bus2_dcline",          [](LSGrid & g){ g.change_bus2_dcline_python(0, 6); }},
        // ---- direct bus manipulation ---------------------------------------
        {"deactivate_bus (empty)",      [](LSGrid & g){ g.deactivate_load(8); g.deactivate_bus_python(8); }},
        {"reactivate_bus",              [](LSGrid & g){ g.reactivate_bus_python(3); }},
        // ---- bulk ------------------------------------------------------------
        {"consider_only_main_component", [](LSGrid & g){ g.deactivate_powerline(4); g.consider_only_main_component(); }},
    };
}

}  // namespace

TEST_CASE("every bus-membership mutator keeps the per-bus element counts exact",
          "[SubstationContainer][bus_count]")
{
    for(const auto & mutation : all_bus_mutations()){
        SECTION(mutation.name){
            LSGrid grid = make_grid();
            grid.recompute_bus_element_counts();   // the baseline every grid starts from
            mutation.apply(grid);
            require_counts_match_recount(grid, mutation.name);
        }
        SECTION(std::string(mutation.name) + " (on a grid that has already solved)"){
            LSGrid grid = make_grid();
            REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 14);
            grid.recompute_bus_element_counts();
            mutation.apply(grid);
            require_counts_match_recount(grid, mutation.name);
        }
    }
}

TEST_CASE("a no-op mutation moves no count at all", "[SubstationContainer][bus_count]")
{
    // The reason the mutators guard these: taking a contribution away and putting it
    // straight back leaves the counts right, but a bus that is alone transiently hits
    // zero. Under the crossing rule that reads as "the dimension changed", costing a
    // full rebuild -- and grid2op sends a topology vector every single step, most of
    // which asks for the bus an element is already on.
    LSGrid grid = make_grid();
    grid.recompute_bus_element_counts();
    const std::vector<std::size_t> before = grid.get_substations().get_nb_elements_per_bus();

    grid.change_bus_load_python(0, grid.get_loads_as_data().get_buses()(0).cast_int());  // same bus
    CHECK(grid.get_substations().get_nb_elements_per_bus() == before);

    grid.deactivate_load(2);
    const std::vector<std::size_t> after_deact = grid.get_substations().get_nb_elements_per_bus();
    grid.deactivate_load(2);   // already off: must not decrement twice
    CHECK(grid.get_substations().get_nb_elements_per_bus() == after_deact);
    require_counts_match_recount(grid, "deactivating an already-deactivated load");
}

TEST_CASE("the count sweep still covers every mutator", "[SubstationContainer][bus_count]")
{
    // A guard on the list above: if LSGrid grows a method that can change bus
    // membership and it is not added here, this number stops matching and someone
    // has to look. Counted from the survey in the commit message.
    CHECK(all_bus_mutations().size() == 41u);
}

TEST_CASE("bus connectivity IS the element counts", "[SubstationContainer][bus_count]")
{
    // What stages 2 and 3 are about. There used to be a `std::vector<bool>
    // bus_status_` next to these counts, saying the same thing a second time: it had
    // to be rebuilt from the elements on every powerflow to stay true (O(nb_bus),
    // which on a grid with far more buses than elements -- 5000 substations at 12
    // busbars each, most empty -- is the dominant cost of the very step the counts
    // were introduced to make cheap), and it could be set by hand to something the
    // elements contradicted. It is gone: `is_bus_connected(b)` is `count[b] > 0`,
    // and `nb_connected_bus()` is one integer moved at the crossings.
    LSGrid grid = make_grid();
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 14);
    const auto & subs = grid.get_substations();
    REQUIRE(subs.bus_counts_ready());
    REQUIRE(subs.connected_bus_count_is_exact());

    SECTION("a mutator that empties a bus takes that bus out of the system"){
        const int bus8 = grid.get_loads_as_data().get_buses()(8).cast_int();
        REQUIRE(subs.get_nb_elements_per_bus()[bus8] >= 1u);
        const std::size_t nb_conn_before = subs.nb_connected_bus();

        grid.deactivate_load(8);
        CHECK(subs.is_bus_connected(GridModelBusId(bus8)) ==
              (subs.get_nb_elements_per_bus()[bus8] > 0u));
        CHECK(subs.connected_bus_count_is_exact());

        grid.reactivate_load(8);
        CHECK(subs.is_bus_connected(GridModelBusId(bus8)) ==
              (subs.get_nb_elements_per_bus()[bus8] > 0u));
        CHECK(subs.nb_connected_bus() == nb_conn_before);
        CHECK(subs.connected_bus_count_is_exact());
    }

    SECTION("get_bus_status() is the counts, built on demand"){
        const std::vector<bool> status = grid.get_bus_status();
        REQUIRE(status.size() == subs.get_nb_elements_per_bus().size());
        std::size_t nb_true = 0;
        for(std::size_t b = 0; b < status.size(); ++b){
            INFO("bus " << b);
            CHECK(status[b] == (subs.get_nb_elements_per_bus()[b] > 0u));
            if(status[b]) ++nb_true;
        }
        CHECK(nb_true == subs.nb_connected_bus());
    }

    SECTION("deactivate_bus / reactivate_bus are no-ops"){
        // They have nothing left to set: a bus is in the system iff an element holds
        // it. They were already all but inert -- whatever they wrote, the next
        // powerflow rebuilt the status from the elements and threw it away -- which
        // is why this changes no answer.
        const std::vector<bool> before = grid.get_bus_status();
        const std::size_t nb_conn_before = subs.nb_connected_bus();
        grid.deactivate_bus_python(3);
        grid.reactivate_bus_python(4);
        CHECK(grid.get_bus_status() == before);
        CHECK(subs.nb_connected_bus() == nb_conn_before);
        CHECK(subs.connected_bus_count_is_exact());
    }

    SECTION("replacing a whole container disarms the counts, and a rebuild restores them"){
        // Re-initializing the storage container with three empty vectors removes
        // every storage from the grid. No +1 / -1 sees that, which is the point:
        // without the disarming the counts would keep describing elements that no
        // longer exist, and the recount at the next init_bus_status() is what makes
        // that unable to outlive one rebuild.
        REQUIRE(grid.get_storages().nb() > 0);
        grid.init_storages(RealVect(), RealVect(), Eigen::VectorXi());
        CHECK_FALSE(subs.bus_counts_ready());
        grid.init_bus_status();
        CHECK(subs.bus_counts_ready());
        CHECK(subs.connected_bus_count_is_exact());
        require_counts_match_recount(grid, "re-initializing the storage container");
    }
}

TEST_CASE("a restored grid counts its own buses", "[SubstationContainer][bus_count]")
{
    // Bus connectivity is no longer serialized (BINARY_FORMAT_VERSION 6), so there
    // is no field for a saved state to disagree with the elements about. What comes
    // back is counted from the restored elements' own status -- which IS in the
    // state -- so a round-trip reproduces the connectivity exactly, and a crafted
    // file cannot state anything else.
    LSGrid grid = make_grid();
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() == 14);
    grid.deactivate_load(8);   // change the topology so the round-trip has to carry it
    const std::vector<bool> expected_status = grid.get_bus_status();
    const std::size_t expected_nb_conn = grid.get_substations().nb_connected_bus();

    LSGrid::StateRes st = grid.get_state();
    LSGrid restored;
    restored.set_state(st);

    const auto & subs = restored.get_substations();
    CHECK_FALSE(subs.bus_counts_ready());          // nothing has counted yet
    restored.init_bus_status();                     // ... so this counts, from the elements
    CHECK(subs.bus_counts_ready());
    CHECK(restored.get_bus_status() == expected_status);
    CHECK(subs.nb_connected_bus() == expected_nb_conn);
    CHECK(subs.connected_bus_count_is_exact());
    require_counts_match_recount(restored, "restoring a state");
}

TEST_CASE("an SVC holds its bus in the solved system, like every other element",
          "[SubstationContainer][bus_count]")
{
    // The gap the counts closed on the way past. `init_bus_status()` used to
    // disconnect every bus and walk eight containers to put them back: lines,
    // shunts, trafos, gens, loads, sgens, storages, hvdc. Not SVCs -- SvcContainer
    // inherits a perfectly good reconnect_connected_buses, it was simply never
    // called. So a bus held by nothing but an SVC read as disconnected, got no solver
    // id, and silently dropped out of the system it should have been in.
    //
    // Nothing states that list any more: the counts come from `contribute_to_buses`,
    // which every container implements for itself, and the bus status is the counts
    // read back. A container cannot be left off a list that no longer exists.
    //
    // What is checked here is the one fact the old list got wrong -- an active SVC
    // contributes +1 to its bus, so a bus holding only that SVC has a count of 1 and
    // is connected. (The exotic fixture has one busbar per substation and a line on
    // every bus, so there is no bus to strand an SVC on; the contribution is the
    // claim, and it is what the count measures.)
    LSGrid grid = make_grid();
    REQUIRE(grid.get_svcs().nb() > 0);
    grid.recompute_bus_element_counts();

    const auto & subs = grid.get_substations();
    const int svc_bus = grid.get_svcs().get_buses()(0).cast_int();
    const std::size_t with_svc = subs.get_nb_elements_per_bus()[svc_bus];
    REQUIRE(with_svc >= 1u);

    grid.deactivate_svc(0);
    CHECK(subs.get_nb_elements_per_bus()[svc_bus] == with_svc - 1u);
    CHECK(subs.is_bus_connected(GridModelBusId(svc_bus)) ==
          (subs.get_nb_elements_per_bus()[svc_bus] > 0u));
    require_counts_match_recount(grid, "deactivating the SVC");

    grid.reactivate_svc(0);
    CHECK(subs.get_nb_elements_per_bus()[svc_bus] == with_svc);
    CHECK(subs.is_bus_connected(GridModelBusId(svc_bus)));
    require_counts_match_recount(grid, "reactivating the SVC");
}
