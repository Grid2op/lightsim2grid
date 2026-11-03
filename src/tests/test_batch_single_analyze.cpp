// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ============================================================================================
// THE POINT OF THE BATCH CLASSES -- AND WHY THEY ARE FAST
// ============================================================================================
//
// TimeSeries, InjectionSweep, ContingencyAnalysis and ScenarioSweep exist for one reason:
// the sparsity pattern of the Jacobian is the SAME for every row of the batch. So the
// expensive symbolic work -- the fill-reducing ordering and the symbolic factorization,
// `analyze` -- is paid ONCE for the whole batch, and every row after it only pays a numeric
// refactorization (`refactorize`), which is a fraction of the price.
//
// That is what makes a batch of thousands of contingencies or time steps cheap, and it is
// the contract of these classes: `nb_analyze == 1` per algorithm, whatever the number of
// rows, whatever happens inside a row.
//
// What a row IS allowed to do: a numeric `factorize` (choosing its pivots again, symbolic
// analysis kept) -- that is what the refactorize fallback does when a row's matrix no longer
// suits the pivots of the last factorization. That is cheap enough and sometimes necessary.
//
// What a row must NEVER do: trigger another `analyze`. A fix that "solves" a failing row by
// re-analyzing it is not a fix, it throws the speed of the batch classes away (a sweep that
// re-analyzes every row is as slow as solving every row from scratch). If one of the tests
// below fails, do not relax it: find what changed the sparsity pattern, or what asked for a
// new analysis, and remove that.
//
// With nb_thread > 1 each worker owns its algorithm and analyzes once, and so does the
// member algorithm that solves the base ("n") case: the figure that must be 1 is the one of
// EACH algorithm, which is why the threaded cases look at get_linear_solver_stats_per_algo()
// and not at the sum.
//
// Every class is run on two grids: a plain one, and one where a generator regulates a REMOTE
// bus. The latter adds bordered voltage-control rows to the Jacobian, which is where rows
// need a numeric re-factorization (see test_klu_pivot_guard.cpp): the case where the
// temptation to re-analyze is.
// ============================================================================================

#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"
#include "light_env/topo_action.hpp"

using ls2g::AlgorithmType;
using ls2g::ContingencyAnalysis;
using ls2g::CplxVect;
using ls2g::InjectionSweep;
using ls2g::LSGrid;
using ls2g::LinearSolverStats;
using ls2g::RealVect;
using ls2g::ScenarioSweep;
using ls2g::TimeSeries;
using ls2g::TopoAction;
using ls2g::cplx_type;
using ls2g::real_type;

namespace {

using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int N_SUB = 4;
const int N_BUSBAR = 2;
// enough rows that "one analyze for the batch" and "one analyze per row" cannot be confused
const int NB_ROWS = 12;

// a 4-substation feeder 0-1-2-3 with a second line between 1 and 2, a slack generator at 0,
// a PV generator and a load at 2, a load at 3. With `remote`, the generator at 2 regulates
// the voltage of bus 3 instead of its own.
LSGrid make_grid(bool remote)
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

    RealVect gen_p(2), gen_v(2), gen_min_q(2), gen_max_q(2);
    Eigen::VectorXi gen_bus(2);
    gen_p << 0., 20.;
    gen_v << 1.04, 1.02;
    gen_min_q << -1000., -1000.;
    gen_max_q << 1000., 1000.;
    gen_bus << 0, 2;
    grid.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus);
    grid.add_gen_slackbus(0, 1.);
    if(remote) grid.set_gen_regulated_bus(1, 3);

    // what a TopoAction needs: the substation of every element and the busbar count
    Eigen::VectorXi load_sub(2), gen_sub(2);
    load_sub << 2, 3;
    gen_sub << 0, 2;
    grid.set_load_to_subid(load_sub);
    grid.set_gen_to_subid(gen_sub);
    grid.set_line_to_sub1_id(from_id);
    grid.set_line_to_sub2_id(to_id);
    grid.set_max_nb_bus_per_sub(N_BUSBAR);

    grid.change_algorithm(AlgorithmType::NR_KLU);
    grid.tell_solver_need_reset();
    return grid;
}

CplxVect flat(const LSGrid & grid)
{
    return CplxVect::Constant(static_cast<Eigen::Index>(grid.total_bus()), cplx_type(1., 0.));
}

// the injections of NB_ROWS rows, moving ~ a few % from one row to the next
struct Rows
{
    RealMat gen_p{NB_ROWS, 2};
    RealMat sgen_p{NB_ROWS, 0};
    RealMat load_p{NB_ROWS, 2};
    RealMat load_q{NB_ROWS, 2};
    Rows()
    {
        for(int i = 0; i < NB_ROWS; ++i){
            const real_type j = static_cast<real_type>((i * 37) % 11) / 10.;
            gen_p(i, 0) = 0.;
            gen_p(i, 1) = 20. + 5. * j;
            load_p(i, 0) = 10. + 4. * j;
            load_p(i, 1) = 50. + 8. * j;
            load_q(i, 0) = 2. + j;
            load_q(i, 1) = 10. + 2. * j;
        }
    }
};

// the contract. `n_algo` algorithms were used (1 single-threaded; the member one plus one
// per worker otherwise): each of them analyzed at most once, and the batch analyzed at all.
void require_one_analyze_per_algorithm(const std::vector<LinearSolverStats> & per_algo)
{
    REQUIRE_FALSE(per_algo.empty());
    std::size_t total = 0;
    for(std::size_t i = 0; i < per_algo.size(); ++i){
        INFO("algorithm #" << i << " (0 is the member one)");
        // NEVER more than one: a row must not ask for a new analysis
        REQUIRE(per_algo[i].nb_analyze <= 1);
        total += per_algo[i].nb_analyze;
    }
    REQUIRE(total >= 1);
}

}  // namespace


TEST_CASE("TimeSeries analyzes once for the whole batch", "[batch][single_analyze]")
{
    const bool remote = GENERATE(false, true);
    INFO("remote regulator: " << remote);
    LSGrid grid = make_grid(remote);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-10).size() > 0);

    TimeSeries ts(grid);
    const Rows rows;
    REQUIRE(ts.compute_Vs(rows.gen_p, rows.sgen_p, rows.load_p, rows.load_q, flat(grid), 30, 1e-10) == 1);

    REQUIRE(ts.get_linear_solver_stats().nb_analyze == 1);
    REQUIRE(ts.get_linear_solver_stats().nb_refactorize > 0);   // the rows did refactorize
    require_one_analyze_per_algorithm(ts.get_linear_solver_stats_per_algo());
}


TEST_CASE("InjectionSweep analyzes once for the whole batch", "[batch][single_analyze]")
{
    const bool remote = GENERATE(false, true);
    INFO("remote regulator: " << remote);
    LSGrid grid = make_grid(remote);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-10).size() > 0);

    InjectionSweep sweep(grid);
    const Rows rows;
    sweep.modify_gen_p(rows.gen_p);
    sweep.modify_sgen_p(rows.sgen_p);
    sweep.modify_load_p(rows.load_p);
    sweep.modify_load_q(rows.load_q);
    sweep.compute(flat(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);

    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
    REQUIRE(sweep.get_linear_solver_stats().nb_refactorize > 0);
    require_one_analyze_per_algorithm(sweep.get_linear_solver_stats_per_algo());
}


TEST_CASE("ContingencyAnalysis analyzes once for the whole batch, whatever the thread count",
          "[batch][single_analyze]")
{
    const bool remote = GENERATE(false, true);
    const int nb_thread = GENERATE(1, 2);
    INFO("remote regulator: " << remote << ", nb_thread: " << nb_thread);
    LSGrid grid = make_grid(remote);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-10).size() > 0);

    ContingencyAnalysis ca(grid);
    ca.set_nb_thread(nb_thread);
    ca.add_all_n1();
    ca.compute(flat(grid), 30, 1e-10);

    if(nb_thread == 1) REQUIRE(ca.get_linear_solver_stats().nb_analyze == 1);
    require_one_analyze_per_algorithm(ca.get_linear_solver_stats_per_algo());
    // the contingencies were solved by refactorizing (an islanding one is skipped, not solved)
    REQUIRE(ca.get_linear_solver_stats().nb_refactorize > 0);
}


TEST_CASE("ScenarioSweep analyzes once for the whole batch", "[batch][single_analyze]")
{
    const bool remote = GENERATE(false, true);
    INFO("remote regulator: " << remote);
    LSGrid grid = make_grid(remote);
    REQUIRE(grid.ac_pf(flat(grid), 30, 1e-10).size() > 0);

    // row 0 is the base case, the others alternate between the two parallel 1-2 lines (1 and
    // 3): the only outages that island nothing
    std::vector<TopoAction> actions(NB_ROWS);
    for(int row = 1; row < NB_ROWS; ++row){
        actions[static_cast<std::size_t>(row)].set_line_status(row % 2 == 1 ? 1 : 3, -1);
    }

    ScenarioSweep sweep(grid);
    sweep.set_topo_actions(actions);
    sweep.compute(flat(grid), 30, 1e-10);
    REQUIRE(sweep.get_status() == 1);

    REQUIRE(sweep.get_linear_solver_stats().nb_analyze == 1);
    REQUIRE(sweep.get_linear_solver_stats().nb_refactorize > 0);
    require_one_analyze_per_algorithm(sweep.get_linear_solver_stats_per_algo());
}
