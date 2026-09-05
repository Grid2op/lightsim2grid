// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

/**
 * Instruction-count audit of the CACHED powerflow path -- a powerflow that
 * follows an already successful one on the same grid, which is what every
 * grid2op step after the first actually runs.
 *
 * Run under callgrind. The point of the driver is that callgrind only ever
 * COLLECTS the powerflow calls themselves: loading the grid, the first (cold)
 * solve that fills the cache, and the mutation applied between two steps are
 * all executed with collection off, so the numbers that come out are the cost
 * of `ac_pf` (and `dc_pf`) alone, per solve.
 *
 *   valgrind --tool=callgrind --instr-atstart=no --collect-atstart=no \
 *            ./profile_cached_pf <grid.lsb> <phase> [nb_solves] [algo]
 *
 * Phases (each is one callgrind run, so the numbers are never mixed):
 *   cold       one solve on a grid that has never solved: the full build.
 *   idem       repeated solves with NOTHING changed in between -- the floor:
 *              what the library spends when the answer is already known.
 *   inj        repeated solves with every load's P and Q moved ~2% in between:
 *              the ordinary grid2op step (injections change, topology does not).
 *   inj_nores  `inj` with result computation switched off: the difference
 *              against `inj` is exactly what compute_results() costs.
 *   dcac       `inj`, but running dc_pf then ac_pf as LightSimBackend does by
 *              default (initdc=True).
 *   topo       repeated solves with one line opened / closed in between: the
 *              cache is invalidated (Ybus sparsity changes), i.e. the cost the
 *              cached path is supposed to be avoiding.
 *   nocache    `inj` with allow_ac_cache_reuse(false): the same solves with the
 *              cache switched off, which is the "what does the cache buy today"
 *              reference.
 */

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "AlgoConfig.hpp"
#include "LSGrid.hpp"

#ifdef LS2G_HAS_CALLGRIND
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION do {} while (0)
#define CALLGRIND_STOP_INSTRUMENTATION  do {} while (0)
#define CALLGRIND_TOGGLE_COLLECT        do {} while (0)
#define CALLGRIND_ZERO_STATS            do {} while (0)
#endif

using ls2g::LSGrid;
using ls2g::CplxVect;
using ls2g::RealVect;
using ls2g::real_type;
using ls2g::AlgorithmType;

namespace {

const int MAX_ITER = 10;
const real_type TOL = 1e-8;

// Move every load ~2% around its base value, in a way that changes on every
// step and is bounded (so the grid stays solvable over a long run). This is
// what a grid2op step does to the grid between two powerflows: the injections
// move, the topology does not.
void perturb_injections(LSGrid & grid,
                        const RealVect & p0,
                        const RealVect & q0,
                        int step)
{
    const int nb_load = static_cast<int>(p0.size());
    for(int load_id = 0; load_id < nb_load; ++load_id){
        // deterministic, cheap, and different for each (load, step) pair
        const real_type phase = static_cast<real_type>((load_id * 7919 + step * 104729) % 1000) / 1000.;
        const real_type factor = 1. + 0.02 * (phase - 0.5) * 2.;
        grid.change_p_load(load_id, p0(load_id) * factor);
        grid.change_q_load(load_id, q0(load_id) * factor);
    }
}

// The Newton loop refactorizes J on every iteration by default
// (RefactorPolicyType::AlwaysRefactor). The `EveryN` and `Chord` policies reuse
// a factorization across iterations; measuring the same solves under each is
// how the audit puts a number on what the refactorization actually costs.
// Set through AlgoConfig rather than a typed setter because AlgorithmSelector
// only exposes the algorithm through the BaseAlgo interface.
void set_refactor_policy(LSGrid & grid, const std::string & policy)
{
    if(policy == "always") return;
    ls2g::AlgoConfig cfg = grid.get_ac_algo_config();
    if(cfg.int_params.size() < 4){
        throw std::runtime_error("the selected AC algorithm has no refactorization policy");
    }
    if(policy == "chord"){
        cfg.int_params[1] = 2;   // RefactorPolicyType::Chord
    } else if(policy.rfind("every", 0) == 0){
        cfg.int_params[1] = 1;   // RefactorPolicyType::EveryN
        cfg.int_params[3] = std::atoi(policy.c_str() + 5);
    } else {
        throw std::runtime_error("unknown refactor policy '" + policy +
                                 "' (expected always, chord or everyN)");
    }
    grid.set_ac_algo_config(cfg);
}

// "KLU" / "SparseLU" are shorthands that pick the AC and the DC member of the
// same linear-solver family; anything else is passed to the registry by name
// (e.g. "NRSing_KLU", "NRRefactorRetry_KLU"), with the DC side left as it is.
void select_algo(LSGrid & grid, const std::string & name)
{
    if(name == "KLU"){
        grid.change_algorithm(AlgorithmType::NR_KLU);
        grid.change_algorithm(AlgorithmType::DC_KLU);
    } else if(name == "SparseLU"){
        grid.change_algorithm(AlgorithmType::NR_SparseLU);
        grid.change_algorithm(AlgorithmType::DC_SparseLU);
    } else {
        grid.change_algorithm(name);
    }
}

}  // namespace


int main(int argc, char ** argv)
{
    if(argc < 3){
        std::cerr << "usage: " << argv[0]
                  << " <grid.lsb> <cold|idem|inj|inj_nores|dcac|topo|nocache>"
                     " [nb_solves] [KLU|SparseLU|<registry name>] [always|chord|everyN]\n";
        return 2;
    }
    const std::string path = argv[1];
    const std::string phase = argv[2];
    const int nb_solves = (argc > 3) ? std::atoi(argv[3]) : 10;
    const std::string algo_name = (argc > 4) ? argv[4] : "KLU";
    const std::string refactor = (argc > 5) ? argv[5] : "always";

    LSGrid grid = LSGrid::load_binary(path);
    try {
        select_algo(grid, algo_name);
    } catch (const std::exception & exc){
        std::cerr << "cannot select " << algo_name << ": " << exc.what() << "\n";
        return 3;
    }

    try {
        set_refactor_policy(grid, refactor);
    } catch (const std::exception & exc){
        std::cerr << "cannot set the refactorization policy: " << exc.what() << "\n";
        return 3;
    }

    const int nb_bus = static_cast<int>(grid.total_bus());
    const RealVect p0 = grid.get_loads().get_target_p();
    const RealVect q0 = grid.get_loads().get_target_q();

    CplxVect V0 = CplxVect::Constant(nb_bus, grid.get_init_vm_pu());

    if(phase == "cold"){
        // Nothing to warm up: the ONE solve we measure is the first one.
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
        const CplxVect V = grid.ac_pf(V0, MAX_ITER, TOL);
        CALLGRIND_TOGGLE_COLLECT;
        CALLGRIND_STOP_INSTRUMENTATION;
        if(V.size() == 0){ std::cerr << "cold solve diverged\n"; return 4; }
        std::cout << "cold: 1 solve, " << grid.get_algo().get_nb_iter() << " iterations\n";
        return 0;
    }

    if(phase == "nocache") grid.allow_ac_cache_reuse(false);
    if(phase == "inj_nores") grid.deactivate_result_computation();

    // ---- warm up: the solve that fills the cache, never collected -----------
    CplxVect V = grid.ac_pf(V0, MAX_ITER, TOL);
    if(V.size() == 0){ std::cerr << "warm-up solve diverged\n"; return 4; }

    // A second uncollected solve, so that what we measure is really the steady
    // state (the first cached solve still refactorizes from a cold KLU object).
    V = grid.ac_pf(V, MAX_ITER, TOL);
    if(V.size() == 0){ std::cerr << "second warm-up solve diverged\n"; return 4; }

    // `topo` needs a line whose opening leaves the grid solvable -- on the
    // pegase cases most lines are radial feeders and opening them islands a
    // bus. Found here, with collection still off, so the scan costs nothing.
    int topo_line_id = -1;
    if(phase == "topo"){
        const int nb_line = static_cast<int>(grid.nb_powerline());
        const int scan_max = (nb_line < 200) ? nb_line : 200;
        for(int line_id = 0; line_id < scan_max; ++line_id){
            if(!grid.get_lines_status()[line_id]) continue;
            grid.deactivate_powerline(line_id);
            const CplxVect Vtry = grid.ac_pf(V, MAX_ITER, TOL);
            grid.reactivate_powerline(line_id);
            const CplxVect Vback = grid.ac_pf(V, MAX_ITER, TOL);
            if(Vtry.size() != 0 && Vback.size() != 0){ topo_line_id = line_id; break; }
        }
        if(topo_line_id < 0){
            std::cerr << "topo: no line found whose opening keeps the grid solvable\n";
            return 5;
        }
        std::cout << "topo: toggling line " << topo_line_id << "\n";
        V = grid.ac_pf(V, MAX_ITER, TOL);
    }

    long total_iter = 0;
    int nb_ok = 0;
    const auto t_start = std::chrono::steady_clock::now();

    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_ZERO_STATS;
    for(int step = 0; step < nb_solves; ++step){
        // ---- what happens BETWEEN two powerflows: never collected ----------
        if(phase == "inj" || phase == "inj_nores" || phase == "dcac" ||
           phase == "nocache"){
            perturb_injections(grid, p0, q0, step);
        } else if (phase == "topo"){
            // open a line on even steps, close it back on odd ones: the sparsity
            // pattern of Ybus changes, which is what retires the cache
            if(step % 2 == 0) grid.deactivate_powerline(topo_line_id);
            else grid.reactivate_powerline(topo_line_id);
        } else if (phase != "idem"){
            std::cerr << "unknown phase '" << phase << "'\n";
            return 2;
        }

        // ---- the powerflow itself: this is what is collected ---------------
        if(phase == "dcac"){
            CplxVect Vdc_init = CplxVect::Constant(nb_bus, 1.);
            CALLGRIND_TOGGLE_COLLECT;
            grid.deactivate_result_computation();
            const CplxVect Vdc = grid.dc_pf(Vdc_init, MAX_ITER, TOL);
            grid.reactivate_result_computation();
            const CplxVect Vac = (Vdc.size() != 0) ? grid.ac_pf(Vdc, MAX_ITER, TOL)
                                                   : CplxVect();
            CALLGRIND_TOGGLE_COLLECT;
            if(Vac.size() == 0){ std::cerr << "dcac diverged at step " << step << "\n"; return 4; }
            V = Vac;
        } else {
            CALLGRIND_TOGGLE_COLLECT;
            const CplxVect Vnew = grid.ac_pf(V, MAX_ITER, TOL);
            CALLGRIND_TOGGLE_COLLECT;
            if(Vnew.size() == 0){ std::cerr << "diverged at step " << step << "\n"; return 4; }
            V = Vnew;
        }
        total_iter += grid.get_algo().get_nb_iter();
        ++nb_ok;
    }
    CALLGRIND_STOP_INSTRUMENTATION;

    const auto t_end = std::chrono::steady_clock::now();
    const double secs = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << phase << ": " << nb_ok << " solves, "
              << (static_cast<double>(total_iter) / nb_ok) << " NR iterations/solve, "
              << (secs / nb_ok * 1e3) << " ms/solve (wall, meaningless under valgrind)\n";
    return 0;
}
