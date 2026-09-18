// Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ScenarioSweep's topological axis (set_topo_actions): one TopoAction per row, played
// as value edits on the sweep's fixed solver layout. The oracle is a one-off ac_pf on a
// copy of the grid with the action really applied. This first version plays
// disconnections only; what it refuses is pinned here too.

#include <complex>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "light_env/topo_action.hpp"
#include "batch_algorithm/LimitViolation.hpp"

using Catch::Approx;
using ls2g::AlgorithmType;
using ls2g::CplxVect;
using ls2g::ElementType;
using ls2g::LSGrid;
using ls2g::LimitViolation;
using ls2g::LimitViolationType;
using ls2g::ViolationElementType;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::TopoAction;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int N_SUB = 4;
const int N_BUSBAR = 2;

// a 4-substation feeder 0-1-2-3 with a second line between 1 and 2, a slack generator
// at 0, a PV generator and a load at 2, a load at 3; two busbars per substation (the
// elements all on busbar 1) so that a "move" can be expressed
LSGrid make_grid(real_type q_lim = 1000.)
{
    LSGrid grid;
    grid.set_sn_mva(100.);
    grid.set_init_vm_pu(1.0);

    const RealVect bus_vn_kv = RealVect::Constant(N_SUB * N_BUSBAR, 138.);
    grid.init_bus(static_cast<unsigned int>(N_SUB), static_cast<unsigned int>(N_BUSBAR), bus_vn_kv, 0, 0);

    const int n_line = 4;
    RealVect branch_r = RealVect::Constant(n_line, 0.01);
    RealVect branch_x = RealVect::Constant(n_line, 0.1);
    const CplxVect branch_h = CplxVect::Zero(n_line);
    Eigen::VectorXi from_id(n_line), to_id(n_line);
    from_id << 0, 1, 2, 1;
    to_id << 1, 2, 3, 2;
    branch_x(3) = 0.15;
    grid.init_powerlines(branch_r, branch_x, branch_h, from_id, to_id);

    RealVect load_p(2), load_q(2);
    load_p << 10., 50.;
    load_q << 2., 10.;
    Eigen::VectorXi load_bus(2);
    load_bus << 2, 3;
    grid.init_loads(load_p, load_q, load_bus);

    // generator 2, at 3, does not regulate voltage (a PQ injection) and is OFF in the
    // base grid: what a row may put back
    RealVect gen_p(3), gen_v(3), gen_q(3), gen_min_q(3), gen_max_q(3);
    Eigen::VectorXi gen_bus(3);
    gen_p << 0., 20., 15.;
    gen_v << 1.04, 1.02, 1.0;
    gen_q << 0., 0., 5.;
    gen_min_q << -q_lim, -q_lim, -q_lim;
    gen_max_q << q_lim, q_lim, q_lim;
    gen_bus << 0, 2, 3;
    const std::vector<bool> regulates = {true, true, false};
    grid.init_generators_full(gen_p, gen_v, gen_q, regulates, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    grid.deactivate_gen(2);

    // what a TopoAction needs: the substation of every element and the busbar count
    Eigen::VectorXi load_sub(2), gen_sub(3);
    load_sub << 2, 3;
    gen_sub << 0, 2, 3;
    grid.set_load_to_subid(load_sub);
    grid.set_gen_to_subid(gen_sub);
    grid.set_line_to_sub1_id(from_id);
    grid.set_line_to_sub2_id(to_id);
    grid.set_max_nb_bus_per_sub(N_BUSBAR);

    grid.change_algorithm(AlgorithmType::NR_SparseLU);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat_start(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), {1.0, 0.});
}

CplxVect reference(const LSGrid & base, const TopoAction & action)
{
    LSGrid grid(base);
    TopoAction act = action;
    act.check_validity(grid);
    act.apply_to_gridmodel(grid);
    return grid.ac_pf(flat_start(grid), 30, 1e-10);
}

// the reference, and the buses it solved (a busbar the row uses is one of them)
CplxVect reference(const LSGrid & base, const TopoAction & action, std::vector<int> & solved)
{
    LSGrid grid(base);
    TopoAction act = action;
    act.check_validity(grid);
    act.apply_to_gridmodel(grid);
    const CplxVect res = grid.ac_pf(flat_start(grid), 30, 1e-10);
    solved.clear();
    for(const auto & b : grid.id_ac_solver_to_me().to_int_vector()) solved.push_back(b);
    return res;
}

void require_row_matches(const ScenarioSweep & sweep, int row, const LSGrid & base, const TopoAction & action)
{
    std::vector<int> solved;
    const CplxVect ref = reference(base, action, solved);
    REQUIRE(ref.size() > 0);
    REQUIRE(!solved.empty());
    REQUIRE(sweep.converged_mask()[static_cast<size_t>(row)]);
    const auto & Vs = sweep.get_voltages();
    for(int b : solved){
        INFO("row " << row << ", bus " << b);
        REQUIRE(std::abs(Vs(row, b) - ref(b)) == Approx(0.).margin(1e-8));
    }
}

}  // namespace

TEST_CASE("a row's disconnections through a TopoAction match a one-off powerflow", "[batch][scenario_sweep][topo]")
{
    LSGrid grid = make_grid();
    REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() > 0);

    std::vector<TopoAction> actions(6);
    actions[1].set_line_status(3, -1);                        // the second 1-2 line, by status
    actions[2].add_element(ElementType::line_ex, 3, -1);      // the same, by one end
    actions[3].add_element(ElementType::load, 0, -1);         // the load at 2
    actions[4].add_element(ElementType::gen, 1, -1);          // the PV generator at 2: bus 2 turns PQ
    actions[5].set_line_status(3, -1);                        // all three at once
    actions[5].add_element(ElementType::load, 0, -1);
    actions[5].add_element(ElementType::gen, 1, -1);

    ScenarioSweep sweep(grid);
    sweep.set_topo_actions(actions);
    sweep.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);
    for(int row = 0; row < 6; ++row) require_row_matches(sweep, row, grid, actions[static_cast<size_t>(row)]);

    // one symbolic analysis for the whole sweep
    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);

    // row 4 really released bus 2: its magnitude left the setpoint, row 0 held it
    REQUIRE(std::abs(sweep.get_voltages()(0, 2)) == Approx(1.02).margin(1e-8));
    REQUIRE(std::abs(sweep.get_voltages()(4, 2)) != Approx(1.02).margin(1e-4));

    // the disconnected branches, masks and actions together
    REQUIRE(sweep.get_row_disconnected_branches(0).empty());
    REQUIRE(sweep.get_row_disconnected_branches(1) == std::vector<int>{3});
    REQUIRE(sweep.get_row_disconnected_branches(2) == std::vector<int>{3});
    REQUIRE(sweep.get_row_disconnected_branches(3).empty());

    // ... and the flow of a disconnected branch is 0
    const auto & amps = sweep.compute_flows();
    REQUIRE(amps(1, 3) == 0.);
    REQUIRE(amps(0, 3) != 0.);
}

TEST_CASE("an action and a mask on one row: both apply, the same element twice is refused",
          "[batch][scenario_sweep][topo]")
{
    LSGrid grid = make_grid();
    std::vector<TopoAction> actions(2);
    actions[0].add_element(ElementType::load, 0, -1);
    actions[1].set_line_status(3, -1);
    BoolMat line_mask = BoolMat::Constant(2, 4, false);
    line_mask(0, 3) = true;

    ScenarioSweep sweep(grid);
    sweep.set_contingency_lines(line_mask);
    sweep.set_topo_actions(actions);
    sweep.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);
    TopoAction both;
    both.add_element(ElementType::load, 0, -1);
    both.set_line_status(3, -1);
    require_row_matches(sweep, 0, grid, both);
    require_row_matches(sweep, 1, grid, actions[1]);
    REQUIRE(sweep.get_row_disconnected_branches(0) == std::vector<int>{3});

    // line 3 both masked and in the action of row 1
    line_mask(1, 3) = true;
    sweep.set_contingency_lines(line_mask);
    REQUIRE_THROWS_AS(sweep.compute(flat_start(grid), 30, 1e-10), std::runtime_error);
}

TEST_CASE("what this version refuses: moves, reconnections, DC, invalid actions",
          "[batch][scenario_sweep][topo]")
{
    LSGrid grid = make_grid();

    SECTION("moving the slack generator is refused at compute") {
        std::vector<TopoAction> actions(2);
        actions[1].add_element(ElementType::gen, 0, 2);
        ScenarioSweep sweep(grid);
        sweep.set_topo_actions(actions);
        REQUIRE_THROWS_AS(sweep.compute(flat_start(grid), 30, 1e-10), std::runtime_error);
    }
    SECTION("a no-op reconnection of a connected line is a plain row") {
        std::vector<TopoAction> actions(1);
        actions[0].set_line_status(3, 1);
        ScenarioSweep plain(grid);
        plain.set_topo_actions(actions);   // line 3 is connected: nothing to do
        plain.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(plain.get_status() == 1);
        require_row_matches(plain, 0, grid, TopoAction());
    }
    SECTION("the DC algorithm refuses a row with an action, plays a do-nothing one") {
        std::vector<TopoAction> actions(2);
        actions[1].set_line_status(3, -1);
        ScenarioSweep sweep(grid);
        sweep.change_algorithm(AlgorithmType::DC_SparseLU);
        sweep.set_topo_actions(actions);
        REQUIRE_THROWS_AS(sweep.compute(flat_start(grid), 30, 1e-10), std::runtime_error);
        sweep.set_topo_actions(std::vector<TopoAction>(2));
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
    }
    SECTION("an invalid action is refused when registered, and nothing is registered") {
        std::vector<TopoAction> actions(2);
        actions[1].add_element(ElementType::load, 18, -1);
        ScenarioSweep sweep(grid);
        REQUIRE_THROWS_AS(sweep.set_topo_actions(actions), std::invalid_argument);
        actions[1] = TopoAction();
        actions[1].add_element(ElementType::load, 0, -2);
        REQUIRE_THROWS_AS(sweep.set_topo_actions(actions), std::invalid_argument);
        // the row count is still free
        REQUIRE_NOTHROW(sweep.set_topo_actions(std::vector<TopoAction>(3)));
    }
}

TEST_CASE("a base-off generator reactivated on its bus: PQ -> PV at constant sparsity", "[batch][scenario_sweep][topo]")
{
    SECTION("a regulating generator: its bus turns PV, held at the set-point") {
        LSGrid grid = make_grid();
        grid.deactivate_gen(1);   // bus 2 is PQ in this base grid
        REQUIRE(grid.ac_pf(flat_start(grid), 30, 1e-10).size() > 0);

        std::vector<TopoAction> actions(4);
        actions[0].add_element(ElementType::gen, 1, 1);        // back on busbar 1 of substation 2
        actions[2].add_element(ElementType::gen, 1, 1);
        actions[2].set_line_status(3, -1);                     // ... with a line out
        actions[3].add_element(ElementType::gen, 1, 1);
        actions[3].add_element(ElementType::load, 0, -1);      // ... and the load of its bus out

        ScenarioSweep sweep(grid);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 4; ++row) require_row_matches(sweep, row, grid, actions[static_cast<size_t>(row)]);
        REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
        REQUIRE(std::abs(sweep.get_voltages()(0, 2)) == Approx(1.02).margin(1e-8));
        REQUIRE(std::abs(sweep.get_voltages()(1, 2)) != Approx(1.02).margin(1e-4));

        // a per-row set-point (modify_gen_v) is what the bus is held at
        RealMat gen_v(2, 3);
        gen_v << 1.04, 1.03, 1.0,
                 1.04, 1.01, 1.0;
        ScenarioSweep with_v(grid);
        with_v.modify_gen_v(gen_v);
        with_v.set_topo_actions(std::vector<TopoAction>(2, actions[0]));
        with_v.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(with_v.get_status() == 1);
        for(int row = 0; row < 2; ++row){
            LSGrid ref_grid(grid);
            ref_grid.change_v_gen(1, gen_v(row, 1));
            require_row_matches(with_v, row, ref_grid, actions[0]);
            REQUIRE(std::abs(with_v.get_voltages()(row, 2)) == Approx(gen_v(row, 1)).margin(1e-8));
        }
    }
    SECTION("a non-regulating generator: an injection, its bus stays PQ") {
        LSGrid grid = make_grid();   // generator 2 (PQ, at bus 3) is off in it
        std::vector<TopoAction> actions(3);
        actions[0].add_element(ElementType::gen, 2, 1);
        actions[2].add_element(ElementType::gen, 2, 1);
        actions[2].add_element(ElementType::gen, 1, -1);       // ... while bus 2 loses its PV generator

        ScenarioSweep sweep(grid);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 3; ++row) require_row_matches(sweep, row, grid, actions[static_cast<size_t>(row)]);
        REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
        // the injection really reached the row: bus 3 is not where it is without it
        REQUIRE(std::abs(sweep.get_voltages()(0, 3) - sweep.get_voltages()(1, 3)) > 1e-4);
    }
    SECTION("refused: a slack participant") {
        LSGrid grid = make_grid();
        grid.add_gen_slackbus(1, 0.5);
        grid.deactivate_gen(1);
        std::vector<TopoAction> actions(1);
        actions[0].add_element(ElementType::gen, 1, 1);
        ScenarioSweep sweep(grid);
        sweep.set_topo_actions(actions);
        REQUIRE_THROWS_AS(sweep.compute(flat_start(grid), 30, 1e-10), std::runtime_error);
    }
}

TEST_CASE("elements moved between busbars: a bus created or merged, at constant sparsity", "[batch][scenario_sweep][topo]")
{
    // substation 2 holds load 0, generator 1, the extremity of lines 1 and 3 (both
    // from 1) and the origin of line 2 (to 3); busbar 2 of substation 2 is bus 6
    LSGrid grid = make_grid();
    const int bus_2b = 2 + N_SUB;

    SECTION("splits, merges, a dangling bus, a generator moved") {
        std::vector<TopoAction> actions(6);
        actions[0].add_element(ElementType::load, 0, 2);          // load 0 and line 3 on busbar 2
        actions[0].add_element(ElementType::line_ex, 3, 2);
        // actions[1]: a plain row (merged back)
        actions[2].add_element(ElementType::line_ex, 3, 2);       // only a line end: a bus with one branch
        actions[3].add_element(ElementType::gen, 1, 2);           // the generator and a line: busbar 2 turns PV, 1 PQ
        actions[3].add_element(ElementType::line_ex, 3, 2);
        actions[4].add_element(ElementType::load, 0, 2);          // a split and a disconnection together
        actions[4].add_element(ElementType::line_ex, 3, 2);       // (line 1 off would island 2-3: the
        actions[4].add_element(ElementType::gen, 1, -1);          // generator goes instead)
        actions[5].add_element(ElementType::load, 0, 2);          // a load alone: an island of one, masked
        ScenarioSweep sweep(grid);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 5; ++row) require_row_matches(sweep, row, grid, actions[static_cast<size_t>(row)]);
        REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
        const auto & Vs = sweep.get_voltages();
        REQUIRE(std::abs(Vs(0, bus_2b)) > 0.);
        REQUIRE(std::abs(Vs(1, bus_2b)) == 0.);
        REQUIRE(std::abs(Vs(3, bus_2b)) == Approx(1.02).margin(1e-8));
        REQUIRE(std::abs(Vs(3, 2)) != Approx(1.02).margin(1e-4));
        // the island of one: masked, the row solved without the load
        REQUIRE(sweep.converged_mask()[5]);
        REQUIRE(std::abs(Vs(5, bus_2b)) == 0.);
        TopoAction load_off;
        load_off.add_element(ElementType::load, 0, -1);
        require_row_matches(sweep, 5, grid, load_off);
        // the flows of a moved line are read with the row's buses: the reference agrees
        {
            LSGrid ref(grid);
            TopoAction act = actions[0];
            act.check_validity(ref);
            act.apply_to_gridmodel(ref);
            REQUIRE(ref.ac_pf(flat_start(ref), 30, 1e-10).size() > 0);
            const auto & amps = sweep.compute_flows();
            const auto res = ref.get_line_res1();
            for(int l = 0; l < 4; ++l){
                INFO("line " << l);
                REQUIRE(amps(0, l) == Approx(std::get<3>(res)(l)).margin(1e-8));
            }
        }
    }
    SECTION("a branch disconnected in the base grid, reconnected: where it was, or on a new busbar") {
        LSGrid off = make_grid();
        off.deactivate_powerline(3);
        std::vector<TopoAction> actions(4);
        actions[0].set_line_status(3, 1);
        actions[1].add_element(ElementType::line_or, 3, 1);
        actions[1].add_element(ElementType::line_ex, 3, 1);
        // actions[2]: plain
        actions[3].add_element(ElementType::line_ex, 3, 2);       // on busbar 2 of substation 2, with the load
        actions[3].add_element(ElementType::load, 0, 2);
        ScenarioSweep sweep(off);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(off), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 4; ++row) require_row_matches(sweep, row, off, actions[static_cast<size_t>(row)]);
        REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
        const auto & amps = sweep.compute_flows();
        REQUIRE(amps(0, 3) != 0.);
        REQUIRE(amps(2, 3) == 0.);
    }
    SECTION("four threads agree with one") {
        std::vector<TopoAction> actions(4);
        actions[0].add_element(ElementType::load, 0, 2);
        actions[0].add_element(ElementType::line_ex, 3, 2);
        actions[2].add_element(ElementType::gen, 1, 2);
        actions[2].add_element(ElementType::line_ex, 3, 2);
        actions[3].add_element(ElementType::gen, 1, -1);
        ScenarioSweep one(grid);
        one.set_topo_actions(actions);
        one.compute(flat_start(grid), 30, 1e-10);
        ScenarioSweep four(grid);
        four.set_nb_thread(4);
        four.set_topo_actions(actions);
        four.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(one.get_status() == 1);
        REQUIRE(four.get_status() == 1);
        for(int row = 0; row < 4; ++row){
            for(int b = 0; b < N_SUB * N_BUSBAR; ++b){
                REQUIRE(std::abs(one.get_voltages()(row, b) - four.get_voltages()(row, b)) == Approx(0.).margin(1e-10));
            }
        }
    }
}

namespace {
// the physical violations a plain sweep (no action) reports on the grid rewired for real
std::vector<LimitViolation> physical_reference(const LSGrid & base, const TopoAction & action)
{
    LSGrid grid(base);
    TopoAction act = action;
    act.check_validity(grid);
    act.apply_to_gridmodel(grid);
    ScenarioSweep plain(grid);
    plain.set_compute_physical_violations(true);
    plain.set_physical_violation_tol_mva(0.);
    RealMat load_p(1, 2);
    load_p << 10., 50.;
    plain.modify_load_p(load_p);
    plain.compute(flat_start(grid), 30, 1e-10);
    REQUIRE(plain.get_status() == 1);
    return plain.get_physical_violations()[0];
}
void require_same_violations(const std::vector<LimitViolation> & got, const std::vector<LimitViolation> & ref)
{
    REQUIRE(got.size() == ref.size());
    for(const LimitViolation & r : ref){
        bool found = false;
        for(const LimitViolation & g : got){
            if(g.element_type != r.element_type || g.element_id != r.element_id || g.violation_type != r.violation_type) continue;
            INFO("bus " << r.element_id);
            REQUIRE(g.value == Approx(r.value).margin(1e-6));
            REQUIRE(g.limit == Approx(r.limit).margin(1e-9));
            found = true;
        }
        REQUIRE(found);
    }
}
}  // namespace

TEST_CASE("the reactive-capability check follows a generator the row moves or reactivates", "[batch][scenario_sweep][topo][physical]")
{
    // +/- 0.5 MVAr per machine: every regulated bus violates, so what is reported is
    // exactly which buses hold a generator in the row, and how much they ask
    SECTION("a generator moved to a new busbar with a line") {
        LSGrid grid = make_grid(0.5);
        std::vector<TopoAction> actions(3);
        actions[0].add_element(ElementType::gen, 1, 2);
        actions[0].add_element(ElementType::line_ex, 3, 2);
        actions[2].add_element(ElementType::gen, 1, -1);
        ScenarioSweep sweep(grid);
        sweep.set_compute_physical_violations(true);
        sweep.set_physical_violation_tol_mva(0.);
        RealMat load_p(3, 2);
        load_p << 10., 50., 10., 50., 10., 50.;
        sweep.modify_load_p(load_p);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 3; ++row){
            INFO("row " << row);
            require_same_violations(sweep.get_physical_violations()[static_cast<size_t>(row)],
                                    physical_reference(grid, actions[static_cast<size_t>(row)]));
        }
        // the moved generator is reported on its new bus (6), no longer on bus 2
        bool at_new = false, at_old = false;
        for(const LimitViolation & v : sweep.get_physical_violations()[0]){
            if(v.element_type != ViolationElementType::BUS) continue;
            if(v.element_id == 2 + N_SUB) at_new = true;
            if(v.element_id == 2) at_old = true;
        }
        REQUIRE(at_new);
        REQUIRE(!at_old);
    }
    SECTION("a generator reactivated on its bus") {
        LSGrid grid = make_grid(0.5);
        grid.deactivate_gen(1);
        std::vector<TopoAction> actions(2);
        actions[0].add_element(ElementType::gen, 1, 1);
        ScenarioSweep sweep(grid);
        sweep.set_compute_physical_violations(true);
        sweep.set_physical_violation_tol_mva(0.);
        RealMat load_p(2, 2);
        load_p << 10., 50., 10., 50.;
        sweep.modify_load_p(load_p);
        sweep.set_topo_actions(actions);
        sweep.compute(flat_start(grid), 30, 1e-10);
        REQUIRE(sweep.get_status() == 1);
        for(int row = 0; row < 2; ++row){
            INFO("row " << row);
            require_same_violations(sweep.get_physical_violations()[static_cast<size_t>(row)],
                                    physical_reference(grid, actions[static_cast<size_t>(row)]));
        }
        bool at_bus_2 = false;
        for(const LimitViolation & v : sweep.get_physical_violations()[0]) if(v.element_id == 2) at_bus_2 = true;
        REQUIRE(at_bus_2);
    }
}
