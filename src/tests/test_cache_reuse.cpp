// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Cache reuse: a powerflow reuses what the previous one of the SAME family built
// (bus labelling, Ybus / Sbus, the pv-pq split, the slack weights) and only
// re-stamps what the grid reports as modified since. On by default since 1.0.0,
// per family, and no longer requiring anything from the caller.
//
// Three things must hold, and this file tests each of them:
//
//  1. IT IS INVISIBLE. Reusing must never change an answer. The bulk of this file
//     is a completeness sweep: every mutating method of LSGrid is applied to a
//     grid that has already solved (and therefore has a live cache), and the
//     result compared to the same mutation on a grid that rebuilds everything. A
//     setter that forgets to invalidate what it touched shows up here as a
//     numerical difference -- which is the failure mode automatic caching makes
//     dangerous, and the reason this sweep exists.
//
//  2. IT IS SAFE. A "nothing changed" claim may cost time, never memory safety.
//     `unset_changes()` used to be an unchecked assertion, and three sequences of
//     documented API calls turned it into an out-of-bounds read. Release wheels
//     are -O3 -DNDEBUG, so those were segfaults, not errors. The sections below
//     reproduce them; they also run under the valgrind / ASan+UBSan passes of the
//     cpp_unit_tests workflow, which is what proves the reuse decision is made
//     before anything is indexed.
//
//  3. IT IS CONTROLLABLE. allow_*_cache_reuse / prevent_*_cache_reuse do what
//     they say, per family, and the two families never write over each other.

#include <complex>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "case_exotic_elements.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::GridModelBusId;
using ls2g::LSGrid;
using ls2g::RealVect;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

// IEEE14 + SVC + storage + 3 HVDC lines + a phase-shifting transformer: every
// container the injection / admittance assembly can touch is non-empty.
LSGrid make_grid(){ return ls2g_test::make_exotic_elements_grid(); }

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()),
                              cplx_type(grid.get_init_vm_pu(), 0.));
}

CplxVect solve_ac(LSGrid & grid){ return grid.ac_pf(flat_start(grid), 30, 1e-10); }
CplxVect solve_dc(LSGrid & grid){ return grid.dc_pf(flat_start(grid), 1, 1e-10); }

// every mutating entry point of LSGrid that can move a number in the powerflow.
// `apply` is run once, on a grid built by make_grid().
struct NamedMutation {
    const char * name;
    std::function<void(LSGrid &)> apply;
};

std::vector<NamedMutation> all_mutations()
{
    return {
        // --- injections -----------------------------------------------------
        {"change_p_load",      [](LSGrid & g){ g.change_p_load(3, 42.); }},
        {"change_q_load",      [](LSGrid & g){ g.change_q_load(3, 17.); }},
        {"change_p_gen",       [](LSGrid & g){ g.change_p_gen(1, 33.); }},
        {"change_q_gen",       [](LSGrid & g){ g.change_q_gen(1, 5.); }},
        {"change_v_gen",       [](LSGrid & g){ g.change_v_gen(1, 1.03); }},
        {"change_p_storage",   [](LSGrid & g){ g.change_p_storage(0, 7.); }},
        {"change_q_storage",   [](LSGrid & g){ g.change_q_storage(0, 3.); }},
        {"change_p_shunt",     [](LSGrid & g){ g.change_p_shunt(0, 2.); }},
        {"change_q_shunt",     [](LSGrid & g){ g.change_q_shunt(0, -12.); }},
        {"change_p_dcline",    [](LSGrid & g){ g.change_p_dcline(0, 11.); }},
        {"change_v1_dcline",   [](LSGrid & g){ g.change_v1_dcline(0, 1.04); }},
        {"change_v2_dcline",   [](LSGrid & g){ g.change_v2_dcline(0, 1.01); }},
        // --- connectivity ---------------------------------------------------
        {"deactivate_load",    [](LSGrid & g){ g.deactivate_load(2); }},
        {"deactivate_gen",     [](LSGrid & g){ g.deactivate_gen(3); }},
        {"deactivate_shunt",   [](LSGrid & g){ g.deactivate_shunt(0); }},
        {"deactivate_storage", [](LSGrid & g){ g.deactivate_storage(0); }},
        {"deactivate_svc",     [](LSGrid & g){ g.deactivate_svc(0); }},
        {"deactivate_powerline",       [](LSGrid & g){ g.deactivate_powerline(4); }},
        {"deactivate_powerline_side1", [](LSGrid & g){ g.deactivate_powerline_side1(5); }},
        {"deactivate_powerline_side2", [](LSGrid & g){ g.deactivate_powerline_side2(6); }},
        {"deactivate_trafo",           [](LSGrid & g){ g.deactivate_trafo(1); }},
        {"deactivate_trafo_side1",     [](LSGrid & g){ g.deactivate_trafo_side1(0); }},
        {"deactivate_dcline",          [](LSGrid & g){ g.deactivate_dcline(1); }},
        {"deactivate_dcline_side1",    [](LSGrid & g){ g.deactivate_dcline_side1(2); }},
        {"deactivate then reactivate a line",
            [](LSGrid & g){ g.deactivate_powerline(4); g.reactivate_powerline(4); }},
        // --- topology (element moved to another bus) ------------------------
        {"change_bus_load",        [](LSGrid & g){ g.change_bus_load_python(0, 4); }},
        {"change_bus_gen",         [](LSGrid & g){ g.change_bus_gen_python(2, 9); }},
        {"change_bus_shunt",       [](LSGrid & g){ g.change_bus_shunt_python(0, 3); }},
        {"change_bus_storage",     [](LSGrid & g){ g.change_bus_storage_python(0, 6); }},
        {"change_bus1_powerline",  [](LSGrid & g){ g.change_bus1_powerline_python(7, 2); }},
        {"change_bus2_powerline",  [](LSGrid & g){ g.change_bus2_powerline_python(8, 9); }},
        {"change_bus1_trafo",      [](LSGrid & g){ g.change_bus1_trafo_python(0, 3); }},
        // --- branch characteristics -----------------------------------------
        {"change_ratio_trafo",     [](LSGrid & g){ g.change_ratio_trafo(0, 1.05); }},
        {"change_shift_trafo_deg", [](LSGrid & g){ g.change_shift_trafo_deg(2, 4.); }},
        // --- slack ----------------------------------------------------------
        {"add_gen_slackbus",       [](LSGrid & g){ g.add_gen_slackbus(1, 1.); }},
        {"add then remove a slack gen",
            [](LSGrid & g){ g.add_gen_slackbus(1, 1.); g.remove_gen_slackbus(1); }},
        {"update_slack_weights_by_id", [](LSGrid & g){
            ls2g::IntVect ids(2); ids << 0, 1; g.update_slack_weights_by_id(ids); }},
        // --- global settings -------------------------------------------------
        {"turnedoff_no_pv",        [](LSGrid & g){ g.turnedoff_no_pv(); }},
        {"set_reference_slack_bus",    [](LSGrid & g){ g.add_gen_slackbus(1, 1.); g.set_reference_slack_bus(6); }},
        {"consider_only_main_component",   [](LSGrid & g){ g.deactivate_powerline(4); g.consider_only_main_component(); }},
        {"set_gen_regulated_bus",      [](LSGrid & g){ g.set_gen_regulated_bus(3, 9); }},
        {"deactivate_bus",             [](LSGrid & g){ g.deactivate_load(8); g.deactivate_bus_python(8); }},
        {"change_algorithm (ac)",  [](LSGrid & g){ g.change_algorithm(AlgorithmType::NRSing_SparseLU); }},
        {"set_sn_mva",             [](LSGrid & g){ g.set_sn_mva(50.); }},
        {"set_init_vm_pu",         [](LSGrid & g){ g.set_init_vm_pu(1.01); }},
    };
}

// mutations that only make sense (or only move a number) in AC
bool is_ac_only(const std::string & name)
{
    return name == "change_q_load" || name == "change_q_gen" || name == "change_v_gen" ||
           name == "change_q_storage" || name == "change_p_shunt" || name == "change_q_shunt" ||
           name == "change_bus_shunt" || name == "deactivate_shunt" || name == "deactivate_svc" ||
           name == "change_v1_dcline" || name == "change_v2_dcline" ||
           name == "turnedoff_no_pv" || name == "change_algorithm (ac)";
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// 1. reuse is invisible: every setter invalidates what it touched
// ---------------------------------------------------------------------------

TEST_CASE("every grid modification survives cache reuse", "[LSGrid][cache_reuse]")
{
    // Reference for each mutation: the same mutation applied to a grid that has
    // never solved and is not allowed to reuse anything, so it assembles the
    // whole problem from the containers. If a setter forgets to invalidate the
    // part of the cache it just invalidated in fact, the "cached" grid solves a
    // slightly different problem and the two voltages differ.
    for(const auto & mutation : all_mutations()){
        SECTION(std::string("ac: ") + mutation.name){
            LSGrid cached = make_grid();
            REQUIRE(solve_ac(cached).size() == 14);  // builds and marks the AC cache
            mutation.apply(cached);
            const CplxVect v_cached = solve_ac(cached);

            LSGrid rebuilt = make_grid();
            rebuilt.allow_cache_reuse(false);
            mutation.apply(rebuilt);
            const CplxVect v_rebuilt = solve_ac(rebuilt);

            REQUIRE(v_cached.size() == v_rebuilt.size());
            CHECK((v_cached - v_rebuilt).norm() < 1e-9);
        }
        if(is_ac_only(mutation.name)) continue;
        SECTION(std::string("dc: ") + mutation.name){
            LSGrid cached = make_grid();
            REQUIRE(solve_dc(cached).size() == 14);
            mutation.apply(cached);
            const CplxVect v_cached = solve_dc(cached);

            LSGrid rebuilt = make_grid();
            rebuilt.allow_cache_reuse(false);
            mutation.apply(rebuilt);
            const CplxVect v_rebuilt = solve_dc(rebuilt);

            REQUIRE(v_cached.size() == v_rebuilt.size());
            CHECK((v_cached - v_rebuilt).norm() < 1e-9);
        }
    }
}

TEST_CASE("a long run of modifications stays identical with and without reuse", "[LSGrid][cache_reuse]")
{
    // the sweep above applies one mutation to a freshly-solved grid; here the
    // cache is carried across many solves, so an invalidation that is merely
    // *late* (rather than missing) also shows up.
    LSGrid cached = make_grid();
    LSGrid rebuilt = make_grid();
    rebuilt.allow_cache_reuse(false);

    for(int step = 0; step < 12; ++step){
        const auto mutate = [step](LSGrid & g){
            g.change_p_load(step % 11, 10. + step);
            g.change_q_load((step + 3) % 11, 2. + 0.5 * step);
            if(step == 3) g.deactivate_powerline(4);
            if(step == 5) g.change_bus_gen_python(2, 9);
            if(step == 7) g.reactivate_powerline(4);
            if(step == 9) g.change_v_gen(1, 1.035);
        };
        mutate(cached);
        mutate(rebuilt);

        // alternate the two families: each must stay on its own cache
        const CplxVect v_cached = (step % 2 == 0) ? solve_ac(cached) : solve_dc(cached);
        const CplxVect v_rebuilt = (step % 2 == 0) ? solve_ac(rebuilt) : solve_dc(rebuilt);
        REQUIRE(v_cached.size() == v_rebuilt.size());
        CHECK((v_cached - v_rebuilt).norm() < 1e-9);
    }
}

// ---------------------------------------------------------------------------
// 2. reuse is safe: a wrong "nothing changed" costs time, never memory safety
// ---------------------------------------------------------------------------

TEST_CASE("a 'nothing changed' claim never makes a powerflow read data that was never built",
          "[LSGrid][cache_reuse][unset_changes]")
{
    // `unset_changes()` asserted "the cached solver-side data matches the grid",
    // for BOTH families at once, and nothing checked it. Every sequence here
    // segfaulted before the cache-consistency guard in _pre_process_solver_impl.
    // It is still reachable today: `unset_changes()` is kept for backward
    // compatibility, and does still mark when automatic reuse is off.
    SECTION("on a grid that was never solved at all"){
        LSGrid grid = make_grid();
        grid.allow_cache_reuse(false);  // so unset_changes() is not a no-op
        grid.unset_changes();           // nothing has been built yet: the claim is false
        REQUIRE(solve_ac(grid).size() == 14);
    }
    SECTION("AC after a DC solve: the AC family has nothing cached"){
        LSGrid grid = make_grid();
        grid.allow_cache_reuse(false);
        REQUIRE(solve_dc(grid).size() == 14);
        grid.unset_changes();
        const CplxVect v = solve_ac(grid);
        REQUIRE(v.size() == 14);

        LSGrid ref = make_grid();
        CHECK((v - solve_ac(ref)).norm() < 1e-9);
    }
    SECTION("DC after an AC solve: the DC family has nothing cached"){
        LSGrid grid = make_grid();
        grid.allow_cache_reuse(false);
        REQUIRE(solve_ac(grid).size() == 14);
        grid.unset_changes();
        const CplxVect v = solve_dc(grid);
        REQUIRE(v.size() == 14);

        LSGrid ref = make_grid();
        CHECK((v - solve_dc(ref)).norm() < 1e-9);
    }
    SECTION("alternating families, marking both after each solve"){
        LSGrid grid = make_grid();
        grid.allow_cache_reuse(false);
        LSGrid ref_ac = make_grid();
        LSGrid ref_dc = make_grid();
        const CplxVect v_ac_ref = solve_ac(ref_ac);
        const CplxVect v_dc_ref = solve_dc(ref_dc);
        for(int round = 0; round < 3; ++round){
            CHECK((solve_ac(grid) - v_ac_ref).norm() < 1e-9);
            grid.unset_changes();
            CHECK((solve_dc(grid) - v_dc_ref).norm() < 1e-9);
            grid.unset_changes();
        }
    }
    SECTION("set_state cannot leave a stale cache behind"){
        // a same-sized grid restored into an LSGrid that already solved: the
        // structural checks alone would not notice, so set_state invalidates.
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);  // grid now carries a live AC cache
        REQUIRE_FALSE(grid.get_ac_algo_controler().need_reset_solver());

        LSGrid other = make_grid();
        other.deactivate_powerline(4);
        other.change_p_load(3, 77.);
        LSGrid::StateRes state = other.get_state();
        grid.set_state(state);

        // cold, whatever this instance had cached before
        CHECK(grid.get_ac_algo_controler().need_reset_solver());
        CHECK(grid.get_dc_algo_controler().need_reset_solver());

        const CplxVect v = solve_ac(grid);
        REQUIRE(v.size() == 14);
        LSGrid ref = make_grid();
        ref.allow_cache_reuse(false);
        ref.deactivate_powerline(4);
        ref.change_p_load(3, 77.);
        CHECK((v - solve_ac(ref)).norm() < 1e-9);
    }
}

TEST_CASE("unset_changes never records a claim the cache cannot back",
          "[LSGrid][cache_reuse][unset_changes]")
{
    // `unset_changes()` marks BOTH families -- that is its historical contract and
    // code written before 1.0.0 relies on it. It is also what used to make it
    // dangerous: a family that had never solved got marked "in sync" next to one
    // that had, and its next powerflow took the "nothing to rebuild" path with an
    // empty bus labelling, indexing it with bus ids in the hundreds. Under
    // -O3 -DNDEBUG (what the wheels ship) that is a segfault, not an error.
    //
    // It still marks both families. What changed is that each one is verified first,
    // so the claim is only ever RECORDED where it is true; a family that cannot back
    // it is retired instead and rebuilds on its next powerflow. The powerflow path
    // no longer re-checks, which is the whole point -- so these sections are what
    // stands between `unset_changes()` and an out-of-bounds read.
    LSGrid ref_ac = make_grid();
    LSGrid ref_dc = make_grid();
    const CplxVect v_ac_ref = solve_ac(ref_ac);
    const CplxVect v_dc_ref = solve_dc(ref_dc);

    SECTION("a family that never solved is retired, not marked"){
        LSGrid grid = make_grid();
        grid.unset_changes();   // neither family has built anything

        CHECK(grid.get_ac_algo_controler().need_reset_solver());
        CHECK(grid.get_dc_algo_controler().need_reset_solver());
        CHECK((solve_ac(grid) - v_ac_ref).norm() < 1e-9);
        CHECK((solve_dc(grid) - v_dc_ref).norm() < 1e-9);
    }
    SECTION("the family that did solve is marked, the other is retired"){
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);   // AC has a live cache, DC has nothing
        grid.unset_changes();

        CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());  // true claim, recorded
        CHECK(grid.get_dc_algo_controler().need_reset_solver());        // false claim, refused
        CHECK((solve_dc(grid) - v_dc_ref).norm() < 1e-9);
        CHECK((solve_ac(grid) - v_ac_ref).norm() < 1e-9);
    }
    SECTION("and the other way round"){
        LSGrid grid = make_grid();
        REQUIRE(solve_dc(grid).size() == 14);
        grid.unset_changes();

        CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());
        CHECK(grid.get_ac_algo_controler().need_reset_solver());
        CHECK((solve_ac(grid) - v_ac_ref).norm() < 1e-9);
    }
    SECTION("a pending topology change is not papered over"){
        // The dangerous direction, and the reason unset_changes() refreshes the bus
        // status before comparing: called with a real change outstanding, "forget
        // past changes" would drop it and solve the OLD system with the new grid.
        // The connectivity comparison is what refuses.
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);
        grid.deactivate_powerline(4);
        grid.unset_changes();
        CHECK(grid.get_ac_algo_controler().need_reset_solver());

        LSGrid ref = make_grid();
        ref.deactivate_powerline(4);
        CHECK((solve_ac(grid) - solve_ac(ref)).norm() < 1e-9);
    }
    SECTION("a genuinely up-to-date cache is left alone, however often it is claimed"){
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);
        for(int i = 0; i < 3; ++i){
            grid.unset_changes();
            CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());
            CHECK((solve_ac(grid) - v_ac_ref).norm() < 1e-9);
        }
    }
    SECTION("turning one family's reuse off does not make the other unsafe"){
        // the sequence that used to segfault: both families marked, only one built
        LSGrid grid = make_grid();
        grid.allow_ac_cache_reuse(false);
        grid.unset_changes();
        CHECK((solve_dc(grid) - v_dc_ref).norm() < 1e-9);

        LSGrid grid2 = make_grid();
        grid2.allow_dc_cache_reuse(false);
        grid2.unset_changes();
        CHECK((solve_ac(grid2) - v_ac_ref).norm() < 1e-9);

        LSGrid grid3 = make_grid();
        REQUIRE(solve_ac(grid3).size() == 14);
        grid3.allow_ac_cache_reuse(false);
        grid3.unset_changes();
        CHECK((solve_dc(grid3) - v_dc_ref).norm() < 1e-9);
    }
}

// ---------------------------------------------------------------------------
// 3b. a cache is one object, and it belongs to one owner
// ---------------------------------------------------------------------------

TEST_CASE("a build into someone else's cache never leaves this grid solving a mixture",
          "[LSGrid][cache_reuse][batch]")
{
    // `pre_process_solver` serves two kinds of caller: this grid (ac_pf / dc_pf /
    // check_solution pass ac_cache_ / dc_cache_) and a foreign builder -- the batch
    // algorithms, which pass a cache they own.
    //
    // For the foreign one, the grid's own cache is still an output: the NR
    // extensions are not built from what the solver was handed, they are pulled out
    // of the LSGrid through `lsgrid_ptr`, and they read the labelling AND the pv-pq
    // split (see fill_voltage_control_solver_data's use of bus_pq). So the build has
    // to publish there. What it must never do is leave that publication looking like
    // a cache: the labelling and the split would be the caller's while the matrix
    // and the injections are still this grid's, which is the right size, passes
    // every structural check, converges, and is wrong.
    //
    // Retiring the snapshot is what prevents it, and that is what these sections
    // pin. Note what is NOT tested here: that the nine parts cannot be handed over
    // separately in the first place. That is not a runtime property any more, it is
    // the type -- SolverSideCache is one object, so there is nothing to mismatch.
    SECTION("AC"){
        LSGrid grid = make_grid();
        const CplxVect v_before = solve_ac(grid);
        REQUIRE(v_before.size() == 14);
        REQUIRE_FALSE(grid.get_ac_algo_controler().need_reset_solver());  // cache is live

        // exactly what a batch does: its own cache, its own control
        ls2g::AcSolverCache foreign;
        const ls2g::AlgoControl foreign_control;  // default ctor: everything changed
        grid.build_solver_input(flat_start(grid), foreign, foreign_control);

        // it built the WHOLE cache into the caller's object ...
        CHECK(foreign.mat.rows() == 14);
        CHECK(foreign.mat.cols() == 14);
        CHECK(foreign.inj.size() == 14);
        CHECK(foreign.slack_weights.size() == 14);
        CHECK(foreign.bus_pq.size() > 0);
        CHECK(foreign.id_solver_to_me.size() == 14);
        // ... published the labelling and the split the NR extensions read back ...
        CHECK(grid.get_ac_pq_solver().size() == foreign.bus_pq.size());
        CHECK(grid.get_ac_pv_solver().size() == foreign.bus_pv.size());
        // ... and retired this grid's own cache rather than leaving it a mixture
        CHECK(grid.get_ac_algo_controler().need_reset_solver());

        // so the grid's own next powerflow is unaffected by the detour
        CHECK((solve_ac(grid) - v_before).norm() < 1e-9);
    }
    SECTION("DC"){
        LSGrid grid = make_grid();
        const CplxVect v_before = solve_dc(grid);
        REQUIRE(v_before.size() == 14);
        REQUIRE_FALSE(grid.get_dc_algo_controler().need_reset_solver());

        ls2g::DcSolverCache foreign;
        const ls2g::AlgoControl foreign_control;
        grid.build_dc_solver_input(flat_start(grid), foreign, foreign_control);

        CHECK(foreign.mat.rows() == 14);
        CHECK(foreign.inj.size() == 14);
        CHECK(grid.get_dc_pq_solver().size() == foreign.bus_pq.size());
        CHECK(grid.get_dc_algo_controler().need_reset_solver());

        CHECK((solve_dc(grid) - v_before).norm() < 1e-9);
    }
    SECTION("a foreign build ignores a 'nothing changed' control"){
        // The own / foreign distinction used to be a runtime flag inside ONE
        // function (`own_cache = _is_own_cache(cache)`); it is now which of two
        // functions you call, and `build_solver_input` hardcodes a full rebuild
        // instead of consulting the reuse policy.
        //
        // This is the section that would catch that being undone. `solver_control`,
        // the connectivity snapshot and the change flags all describe THIS grid and
        // say nothing about the caller's cache -- so honouring a "nothing changed"
        // claim about a foreign cache means re-stamping only part of it and leaving
        // the rest from whatever was there before. Here "before" is nothing at all,
        // so a foreign build that believed the claim would hand back an empty
        // labelling and an empty matrix; with a non-empty one it would be worse than
        // empty, because every size would still look right.
        LSGrid grid = make_grid();
        const CplxVect v_before = solve_ac(grid);
        REQUIRE(v_before.size() == 14);

        ls2g::AcSolverCache foreign;   // nothing has ever been built into it
        ls2g::AlgoControl lying_control;
        lying_control.tell_none_changed();   // ... and a control that says so is fine
        grid.build_solver_input(flat_start(grid), foreign, lying_control);

        // built in full regardless of what the control claimed
        CHECK(foreign.mat.rows() == 14);
        CHECK(foreign.mat.cols() == 14);
        CHECK(foreign.inj.size() == 14);
        CHECK(foreign.slack_weights.size() == 14);
        CHECK(foreign.bus_pq.size() > 0);
        CHECK(foreign.id_solver_to_me.size() == 14);
        CHECK(foreign.id_me_to_solver.size() == grid.total_bus());

        // and the grid's own next powerflow is still unaffected by the detour
        CHECK((solve_ac(grid) - v_before).norm() < 1e-9);
    }
    SECTION("one family's foreign build leaves the other family's cache alone"){
        LSGrid grid = make_grid();
        const CplxVect v_ac_before = solve_ac(grid);
        const CplxVect v_dc_before = solve_dc(grid);

        ls2g::DcSolverCache foreign;
        const ls2g::AlgoControl foreign_control;
        grid.build_dc_solver_input(flat_start(grid), foreign, foreign_control);

        CHECK(grid.get_dc_algo_controler().need_reset_solver());
        CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());  // AC untouched
        CHECK((solve_ac(grid) - v_ac_before).norm() < 1e-9);
        CHECK((solve_dc(grid) - v_dc_before).norm() < 1e-9);
    }
}

// ---------------------------------------------------------------------------
// 4. nothing cached ever crosses a serialization boundary
// ---------------------------------------------------------------------------

TEST_CASE("a deserialized grid always starts with a cold cache", "[LSGrid][cache_reuse][serialization]")
{
    // The solver caches are derived state: a second copy of something the elements
    // already determine. Restoring one from a file would mean trusting a copy that
    // cannot be checked against the elements it claims to describe -- check_grid()
    // can validate an index, it cannot validate "this really is the admittance
    // matrix of this grid" -- so none of it is serialized, and the bus-connectivity
    // snapshot that IS in the layout is not believed on the way back in either.
    const auto check_cold = [](const LSGrid & g){
        // both families ask for a full rebuild ...
        CHECK(g.get_ac_algo_controler().need_reset_solver());
        CHECK(g.get_dc_algo_controler().need_reset_solver());
        // ... and nothing solver-side came back from the file
        CHECK(g.get_ac_pv_solver().size() == 0);
        CHECK(g.get_dc_pv_solver().size() == 0);
        CHECK(g.get_ac_pq_solver().size() == 0);
        CHECK(g.get_dc_pq_solver().size() == 0);
        CHECK(g.get_ac_slack_weights_solver().size() == 0);
        CHECK(g.get_dc_slack_weights_solver().size() == 0);
        CHECK(g.id_me_to_ac_solver().size() == 0);
        CHECK(g.id_me_to_dc_solver().size() == 0);
    };

    SECTION("through get_state / set_state (the pickle path)"){
        LSGrid source = make_grid();
        REQUIRE(solve_ac(source).size() == 14);   // source carries a live AC cache
        REQUIRE(solve_dc(source).size() == 14);   // ... and a live DC one
        LSGrid::StateRes state = source.get_state();

        // restored INTO a grid that has already solved: a fresh LSGrid is cold to
        // begin with, so it would not tell us whether set_state invalidates
        LSGrid restored = make_grid();
        REQUIRE(solve_ac(restored).size() == 14);
        REQUIRE(solve_dc(restored).size() == 14);
        REQUIRE_FALSE(restored.get_ac_algo_controler().need_reset_solver());
        restored.set_state(state);
        check_cold(restored);

        // and it still solves to the same answer as a grid that never serialized
        LSGrid ref = make_grid();
        ref.allow_cache_reuse(false);
        CHECK((solve_ac(restored) - solve_ac(ref)).norm() < 1e-9);
        CHECK((solve_dc(restored) - solve_dc(ref)).norm() < 1e-9);
    }

    SECTION("through save_binary / load_binary"){
        LSGrid source = make_grid();
        REQUIRE(solve_ac(source).size() == 14);
        REQUIRE(solve_dc(source).size() == 14);

        const std::string path = "test_cache_reuse_cold_load.lsb";
        source.save_binary(path);
        LSGrid restored = LSGrid::load_binary(path);
        std::remove(path.c_str());
        check_cold(restored);

        LSGrid ref = make_grid();
        ref.allow_cache_reuse(false);
        CHECK((solve_ac(restored) - solve_ac(ref)).norm() < 1e-9);
        CHECK((solve_dc(restored) - solve_dc(ref)).norm() < 1e-9);
    }

    SECTION("there is no cache metadata left in StateRes to hand-edit"){
        // This section used to poison StateRes' bus-connectivity photograph -- the
        // one piece of cache metadata a serialized grid carried -- and check that a
        // crafted file still could not authorize any reuse. That field is gone: a
        // bus entering or leaving the solved system is reported by
        // SubstationContainer's element counts as it happens, so nothing about the
        // cache needs to survive a save any more. The attack surface was removed
        // rather than defended, which is why there is nothing to poison here.
        //
        // What remains testable is the property it was protecting: a restored grid
        // is cold and answers exactly like one that never cached. Element status
        // (which the counts are derived from) IS still serialized, so a faithful
        // round trip must not smuggle warmth in with it.
        LSGrid source = make_grid();
        REQUIRE(solve_ac(source).size() == 14);
        source.deactivate_powerline(4);          // a real connectivity change
        REQUIRE(solve_ac(source).size() == 14);
        LSGrid::StateRes state = source.get_state();

        LSGrid restored = make_grid();
        REQUIRE(solve_ac(restored).size() == 14);   // warm, and for a DIFFERENT topology
        restored.set_state(state);
        check_cold(restored);

        LSGrid ref = make_grid();
        ref.allow_cache_reuse(false);
        ref.deactivate_powerline(4);
        CHECK((solve_ac(restored) - solve_ac(ref)).norm() < 1e-9);
    }
}

// ---------------------------------------------------------------------------
// 3. reuse is controllable
// ---------------------------------------------------------------------------

TEST_CASE("cache reuse is on by default and marks only the family that ran", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    CHECK(grid.get_allow_cache_reuse());
    CHECK(grid.get_allow_ac_cache_reuse());
    CHECK(grid.get_allow_dc_cache_reuse());

    // before any powerflow both families ask for a full rebuild
    CHECK(grid.get_ac_algo_controler().need_reset_solver());
    CHECK(grid.get_dc_algo_controler().need_reset_solver());

    REQUIRE(solve_ac(grid).size() == 14);
    // ... and after an AC solve, ONLY the AC family is marked in sync
    CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());
    CHECK_FALSE(grid.get_ac_algo_controler().has_dimension_changed());
    CHECK(grid.get_dc_algo_controler().need_reset_solver());

    REQUIRE(solve_dc(grid).size() == 14);
    CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());
    CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());

    // a modification invalidates both (it is a modification of the grid, which
    // both families model)
    grid.change_p_load(0, 25.);
    CHECK(grid.get_ac_algo_controler().need_recompute_sbus());
    CHECK(grid.get_dc_algo_controler().need_recompute_sbus());
}

TEST_CASE("the two families keep separate solver-side data", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    // nothing built yet
    CHECK(grid.get_ac_pv_solver().size() == 0);
    CHECK(grid.get_dc_pv_solver().size() == 0);

    REQUIRE(solve_ac(grid).size() == 14);
    CHECK(grid.get_ac_pv_solver().size() > 0);
    CHECK(grid.get_ac_slack_weights_solver().size() > 0);
    CHECK(grid.get_dc_pv_solver().size() == 0);          // untouched by the AC solve
    CHECK(grid.get_dc_slack_weights_solver().size() == 0);

    const auto ac_pv = grid.get_ac_pv_solver();
    const auto ac_pq = grid.get_ac_pq_solver();
    const RealVect ac_sw = grid.get_ac_slack_weights_solver();

    REQUIRE(solve_dc(grid).size() == 14);
    CHECK(grid.get_dc_pv_solver().size() > 0);
    CHECK(grid.get_dc_slack_weights_solver().size() > 0);

    // the DC solve must not have written over the AC family's copies
    REQUIRE(grid.get_ac_pv_solver().size() == ac_pv.size());
    REQUIRE(grid.get_ac_pq_solver().size() == ac_pq.size());
    REQUIRE(grid.get_ac_slack_weights_solver().size() == ac_sw.size());
    for(std::size_t k = 0; k < ac_pv.size(); ++k) CHECK(grid.get_ac_pv_solver()(k) == ac_pv(k));
    for(std::size_t k = 0; k < ac_pq.size(); ++k) CHECK(grid.get_ac_pq_solver()(k) == ac_pq(k));
    CHECK((RealVect(grid.get_ac_slack_weights_solver()) - ac_sw).norm() == Approx(0.).margin(1e-14));
}

TEST_CASE("prevent_*_cache_reuse forces the next powerflow to rebuild", "[LSGrid][cache_reuse]")
{
    SECTION("per family"){
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);
        REQUIRE(solve_dc(grid).size() == 14);
        CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());
        CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());

        grid.prevent_ac_cache_reuse();
        CHECK(grid.get_ac_algo_controler().need_reset_solver());
        CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());  // DC untouched

        grid.prevent_dc_cache_reuse();
        CHECK(grid.get_dc_algo_controler().need_reset_solver());

        // and the rebuilt answers are the same ones
        LSGrid ref = make_grid();
        ref.allow_cache_reuse(false);
        CHECK((solve_ac(grid) - solve_ac(ref)).norm() < 1e-9);
        CHECK((solve_dc(grid) - solve_dc(ref)).norm() < 1e-9);
    }
    SECTION("both, and under its historical name"){
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);
        REQUIRE(solve_dc(grid).size() == 14);
        grid.tell_solver_need_reset();  // == prevent_cache_reuse()
        CHECK(grid.get_ac_algo_controler().need_reset_solver());
        CHECK(grid.get_dc_algo_controler().need_reset_solver());
        CHECK(solve_ac(grid).size() == 14);
    }
    SECTION("it is a one-shot invalidation, not a mode"){
        LSGrid grid = make_grid();
        REQUIRE(solve_ac(grid).size() == 14);
        grid.prevent_cache_reuse();
        REQUIRE(solve_ac(grid).size() == 14);
        CHECK(grid.get_allow_cache_reuse());                             // still enabled
        CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());   // and marked again
    }
}

TEST_CASE("allow_cache_reuse(false) keeps a family rebuilding for good", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    grid.allow_ac_cache_reuse(false);
    CHECK_FALSE(grid.get_allow_ac_cache_reuse());
    CHECK(grid.get_allow_dc_cache_reuse());
    CHECK_FALSE(grid.get_allow_cache_reuse());  // "both" is false as soon as one is

    for(int i = 0; i < 3; ++i){
        REQUIRE(solve_ac(grid).size() == 14);
        // never marked in sync, however many times it solves
        CHECK(grid.get_ac_algo_controler().need_reset_solver());
    }
    // the DC family is unaffected by the AC switch
    REQUIRE(solve_dc(grid).size() == 14);
    CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());

    // ... and unset_changes() cannot re-enable what the flag forbids
    grid.unset_changes();
    LSGrid ref = make_grid();
    ref.allow_cache_reuse(false);
    CHECK((solve_ac(grid) - solve_ac(ref)).norm() < 1e-9);
}

TEST_CASE("unset_changes is a no-op when both families cache automatically", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    REQUIRE(solve_ac(grid).size() == 14);
    grid.prevent_cache_reuse();
    REQUIRE(grid.get_ac_algo_controler().need_reset_solver());
    grid.unset_changes();  // both families are automatic: nothing must happen
    CHECK(grid.get_ac_algo_controler().need_reset_solver());
    CHECK(grid.get_dc_algo_controler().need_reset_solver());
}

TEST_CASE("changing the algorithm invalidates that family's cache", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    REQUIRE(solve_ac(grid).size() == 14);
    REQUIRE(solve_dc(grid).size() == 14);
    CHECK_FALSE(grid.get_ac_algo_controler().need_reset_solver());
    CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());

    grid.change_algorithm(AlgorithmType::NRSing_SparseLU);  // an AC solver
    CHECK(grid.get_ac_algo_controler().need_reset_solver());
    CHECK_FALSE(grid.get_dc_algo_controler().need_reset_solver());  // DC keeps its cache

    const CplxVect v = solve_ac(grid);
    REQUIRE(v.size() == 14);
    LSGrid ref = make_grid();
    ref.allow_cache_reuse(false);
    ref.change_algorithm(AlgorithmType::NRSing_SparseLU);
    CHECK((v - solve_ac(ref)).norm() < 1e-9);
}

TEST_CASE("a divergence keeps the cached data but resets the algorithm", "[LSGrid][cache_reuse]")
{
    // Divergence is a numerical failure, not a data one: Ybus / Sbus and the
    // labelling still describe the grid, and are kept. What is thrown away is the
    // algorithm's own state (half-converged iterate, factorization of a system it
    // gave up on), and the next solve rebuilds those from the cached matrices.
    LSGrid grid = make_grid();
    REQUIRE(solve_ac(grid).size() == 14);

    // make it diverge: one iteration is not enough from a flat start
    const CplxVect diverged = grid.ac_pf(flat_start(grid), 1, 1e-12);
    CHECK(diverged.size() == 0);
    CHECK_FALSE(grid.get_algo().converged());

    // the same grid, solved properly again, must reach the reference solution
    const CplxVect v = solve_ac(grid);
    REQUIRE(v.size() == 14);
    LSGrid ref = make_grid();
    ref.allow_cache_reuse(false);
    CHECK((v - solve_ac(ref)).norm() < 1e-9);

    // and so must a grid modified while it was in the diverged state
    grid.change_p_load(2, 55.);
    const CplxVect v2 = solve_ac(grid);
    REQUIRE(v2.size() == 14);
    LSGrid ref2 = make_grid();
    ref2.allow_cache_reuse(false);
    ref2.change_p_load(2, 55.);
    CHECK((v2 - solve_ac(ref2)).norm() < 1e-9);
}

TEST_CASE("a divergence recovers under every built-in AC algorithm", "[LSGrid][cache_reuse]")
{
    // Each algorithm keeps derived state of its own -- the Jacobian sparsity for
    // Newton-Raphson, Bp / Bpp for the fast-decoupled ones -- which reset()
    // discards. After a divergence the next powerflow must therefore rebuild those
    // from the (kept) cached matrices rather than reuse them.
    // check_solution() is thrown into the sequence on purpose: it runs the same
    // pre-processing without ever calling the algorithm, so it must not consume
    // the "rebuild your internals" the diverged solve asked for.
    // The fast-decoupled solvers do not support the hvdc angle droop that line 1
    // of this grid has enabled: disconnect it, on the test grid and on its
    // reference alike.
    const auto disable_droop = [](LSGrid & g){ g.deactivate_dcline(1); };
    for(const auto algo : {AlgorithmType::NR_SparseLU,
                           AlgorithmType::NRSing_SparseLU,
                           AlgorithmType::FDPF_XB_SparseLU,
                           AlgorithmType::FDPF_BX_SparseLU}){
        LSGrid grid = make_grid();
        disable_droop(grid);
        grid.change_algorithm(algo);
        REQUIRE(grid.ac_pf(flat_start(grid), 100, 1e-8).size() == 14);

        CHECK(grid.ac_pf(flat_start(grid), 1, 1e-12).size() == 0);  // diverges

        const CplxVect probe = grid.check_solution(flat_start(grid), false);
        CHECK(probe.size() == 14);

        const CplxVect v = grid.ac_pf(flat_start(grid), 100, 1e-8);
        REQUIRE(v.size() == 14);
        LSGrid ref = make_grid();
        disable_droop(ref);
        ref.allow_cache_reuse(false);
        ref.change_algorithm(algo);
        CHECK((v - ref.ac_pf(flat_start(ref), 100, 1e-8)).norm() < 1e-7);
    }
}

TEST_CASE("a copy of a grid does not inherit its cache", "[LSGrid][cache_reuse]")
{
    LSGrid grid = make_grid();
    REQUIRE(solve_ac(grid).size() == 14);
    grid.change_p_load(1, 61.);

    LSGrid copy(grid);   // the copy constructor reset()s the solver state
    CHECK(copy.get_ac_algo_controler().need_reset_solver());
    CHECK(copy.get_allow_cache_reuse());  // ... but the settings travel with it

    const CplxVect v = solve_ac(copy);
    REQUIRE(v.size() == 14);
    LSGrid ref = make_grid();
    ref.allow_cache_reuse(false);
    ref.change_p_load(1, 61.);
    CHECK((v - solve_ac(ref)).norm() < 1e-9);
}
