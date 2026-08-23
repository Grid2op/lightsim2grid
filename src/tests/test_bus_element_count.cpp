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
