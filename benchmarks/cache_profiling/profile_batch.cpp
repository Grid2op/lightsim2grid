// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

/**
 * Instruction-count audit of the BATCH algorithms -- TimeSeries and
 * ContingencyAnalysis -- the sibling of profile_cached_pf.cpp for the code that
 * runs one solve per row on a fixed Jacobian sparsity.
 *
 * Same discipline: callgrind only ever COLLECTS the call under audit. Loading the
 * grid, the plain powerflow that seeds the batch, registering the rows and reading
 * the results back all run with collection off.
 *
 *   valgrind --tool=callgrind --instr-atstart=no --collect-atstart=no \
 *            ./profile_batch <grid.lsb> <phase> [nb_rows] [algo] [always] [trace_file]
 *
 * The positional layout is that of profile_cached_pf so the same scripts drive
 * both; the 5th argument (a refactorization policy there) is accepted and ignored.
 *
 * Phases:
 *   ts_ac / ts_dc      TimeSeries::compute() over nb_rows rows whose loads move
 *                      ~2% from row to row (the shape of a chronics file): the
 *                      per-row cost of a solve plus everything the batch does
 *                      around it -- Sbus assembly, result storage.
 *   ts_flows           after such a compute(), compute_flows() + compute_power_flows()
 *                      alone: reading nb_rows x nb_bus voltages back into
 *                      nb_rows x nb_branch flows.
 *   ca_ac / ca_dc      ContingencyAnalysis::compute() over the first nb_rows N-1
 *                      contingencies: the Ybus edit, the connectivity check, the
 *                      solve, the result storage.
 *   ca_ac_mask / ca_dc_mask
 *                      idem with handle_disconnected_grid: the islanding
 *                      contingencies are solved on the surviving component instead
 *                      of being skipped, which costs a component labelling per
 *                      contingency up front.
 *   ca_flows           the flows of a ca_ac run, as ts_flows.
 *   ca_construct       what a grid2op loop pays per step when it builds a fresh
 *                      ContingencyAnalysis from an already-solved grid: the copy of
 *                      the grid, 8 contingencies, compute(). One construction per
 *                      "row"; nothing is shared between two of them.
 *
 * The "per solve" figure the scripts print is per ROW here (per construction for
 * ca_construct). Running a compute phase with nb_rows = 1 and again with N gives,
 * by difference, the fixed cost of a compute() (the rebuild of the solver input
 * and the "n" solve) apart from the marginal cost of a row -- see
 * run_profile_batch.sh, which does exactly that.
 *
 * With a trace_file the batch's answer is written with 17 significant digits, one
 * "step" per row: the row's complex voltages (grid numbering) for the compute
 * phases, its (amps, MW) pairs for the flow phases, so compare_traces.py reads
 * both. A batch does not report per-row iteration counts, so the "iter" field of
 * a step carries the row's converged flag instead.
 */

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <sys/resource.h>

#include "LSGrid.hpp"
#include "batch_algorithm/BaseBatchSweep.hpp"

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
using ls2g::TimeSeries;
using ls2g::ContingencyAnalysis;

namespace {

const real_type TOL = 1e-8;
const int MAX_ITER = 10;
const int NB_CONT_CONSTRUCT = 8;

using RealMat = ls2g::BaseBatchSolverSynch::RealMat;
using CplxMat = ls2g::BaseBatchSolverSynch::CplxMat;

// same as profile_cached_pf: "KLU" / "SparseLU" pick the AC and DC member of a
// linear-solver family, anything else goes to the registry by name.
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

// A batch inherits the AC algorithm of the grid it was built from; its DC side has
// to be asked for. Done BEFORE anything is registered on the batch: change_algorithm
// clears the object.
template<class Batch>
void select_batch_family(Batch & batch, const std::string & algo_name, bool dc)
{
    if(!dc) return;
    if(algo_name == "SparseLU") batch.change_algorithm(AlgorithmType::DC_SparseLU);
    else batch.change_algorithm(AlgorithmType::DC_KLU);
}

// The chronics: every load's P and Q moved ~2% around its base value, differently
// on every row -- the same deterministic pattern profile_cached_pf applies between
// two of its solves, laid out as the matrices TimeSeries takes.
void fill_time_series(const LSGrid & grid, int nb_rows,
                      RealMat & gen_p, RealMat & load_p, RealMat & load_q)
{
    const RealVect p0 = grid.get_loads().get_target_p();
    const RealVect q0 = grid.get_loads().get_target_q();
    const RealVect g0 = grid.get_generators().get_target_p();
    const int nb_load = static_cast<int>(p0.size());
    gen_p = RealMat(nb_rows, g0.size());
    gen_p.rowwise() = g0.transpose();
    load_p = RealMat(nb_rows, nb_load);
    load_q = RealMat(nb_rows, nb_load);
    for(int step = 0; step < nb_rows; ++step){
        for(int load_id = 0; load_id < nb_load; ++load_id){
            const real_type phase = static_cast<real_type>((load_id * 7919 + step * 104729) % 1000) / 1000.;
            const real_type factor = 1. + 0.02 * (phase - 0.5) * 2.;
            load_p(step, load_id) = p0(load_id) * factor;
            load_q(step, load_id) = q0(load_id) * factor;
        }
    }
}

// the first `nb` branches (lines first, then trafos), one N-1 contingency each
std::vector<int> first_n1(const LSGrid & grid, int nb)
{
    const int nb_branch = static_cast<int>(grid.nb_powerline() + grid.nb_trafo());
    std::vector<int> res;
    for(int i = 0; i < nb && i < nb_branch; ++i) res.push_back(i);
    return res;
}

struct Trace {
    std::vector<int> flags;              // one per step: the converged flag
    std::vector<std::vector<std::pair<real_type, real_type> > > rows;
};

void trace_voltages(Trace & trace, const CplxMat & V, const std::vector<char> & converged)
{
    for(Eigen::Index i = 0; i < V.rows(); ++i){
        trace.flags.push_back(static_cast<size_t>(i) < converged.size() ? converged[i] : 0);
        std::vector<std::pair<real_type, real_type> > row;
        row.reserve(V.cols());
        for(Eigen::Index j = 0; j < V.cols(); ++j) row.emplace_back(V(i, j).real(), V(i, j).imag());
        trace.rows.push_back(std::move(row));
    }
}

void trace_flows(Trace & trace, const RealMat & amps, const RealMat & mw, const std::vector<char> & converged)
{
    for(Eigen::Index i = 0; i < amps.rows(); ++i){
        trace.flags.push_back(static_cast<size_t>(i) < converged.size() ? converged[i] : 0);
        std::vector<std::pair<real_type, real_type> > row;
        row.reserve(amps.cols());
        for(Eigen::Index j = 0; j < amps.cols(); ++j) row.emplace_back(amps(i, j), mw(i, j));
        trace.rows.push_back(std::move(row));
    }
}

bool write_trace(const std::string & path, const Trace & trace)
{
    std::ofstream out(path);
    if(!out) return false;
    out << std::setprecision(17) << std::scientific;
    for(std::size_t step = 0; step < trace.rows.size(); ++step){
        out << "step " << step << " iter " << trace.flags[step] << "\n";
        for(const auto & v : trace.rows[step]) out << v.first << " " << v.second << "\n";
    }
    return true;
}

long peak_rss_kb()
{
    struct rusage usage;
    if(getrusage(RUSAGE_SELF, &usage) != 0) return -1;
    return usage.ru_maxrss;
}

}  // namespace


int main(int argc, char ** argv)
{
    if(argc < 3){
        std::cerr << "usage: " << argv[0]
                  << " <grid.lsb> <ts_ac|ts_dc|ts_flows|ca_ac|ca_dc|ca_ac_mask|ca_dc_mask|ca_flows|ca_construct>"
                     " [nb_rows] [KLU|SparseLU|<registry name>] [always] [trace_file]\n";
        return 2;
    }
    const std::string path = argv[1];
    const std::string phase = argv[2];
    const int nb_rows = (argc > 3) ? std::atoi(argv[3]) : 10;
    const std::string algo_name = (argc > 4) ? argv[4] : "KLU";
    // argv[5]: the refactorization policy of profile_cached_pf, accepted for the
    // scripts' sake and ignored
    const std::string trace_path = (argc > 6) ? argv[6] : "";

    const bool is_ts = phase.rfind("ts_", 0) == 0;
    const bool is_ca = phase.rfind("ca_", 0) == 0;
    const bool dc = phase == "ts_dc" || phase == "ca_dc" || phase == "ca_dc_mask";
    const bool mask = phase == "ca_ac_mask" || phase == "ca_dc_mask";
    const bool flows = phase == "ts_flows" || phase == "ca_flows";
    if(!is_ts && !is_ca){
        std::cerr << "unknown phase '" << phase << "'\n";
        return 2;
    }

    LSGrid grid = LSGrid::load_binary(path);
    try {
        select_algo(grid, algo_name);
    } catch (const std::exception & exc){
        std::cerr << "cannot select " << algo_name << ": " << exc.what() << "\n";
        return 3;
    }

    // The seed every batch starts from is a solved grid: that is what the python
    // wrappers hand over (backend.V), and what a grid2op loop holds at every step.
    const int nb_bus = static_cast<int>(grid.total_bus());
    const CplxVect V0 = CplxVect::Constant(nb_bus, grid.get_init_vm_pu());
    const CplxVect Vsolved = grid.ac_pf(V0, MAX_ITER, TOL);
    if(Vsolved.size() == 0){ std::cerr << "the seed solve diverged\n"; return 4; }

    Trace trace;
    long total_refactor = 0;
    int nb_solved = 0;
    int nb_converged = 0;
    int nb_measured = 0;   // rows (or constructions) the collected region covers
    double secs = 0.;

    if(phase == "ca_construct"){
        // ---- one fresh object per "row", nothing warmed up on purpose ----------
        const std::vector<int> conts = first_n1(grid, NB_CONT_CONSTRUCT);
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_ZERO_STATS;
        for(int step = 0; step < nb_rows; ++step){
            const auto t0 = std::chrono::steady_clock::now();
            CALLGRIND_TOGGLE_COLLECT;
            ContingencyAnalysis ca(grid);
            ca.add_multiple_n1(conts);
            ca.compute(Vsolved, MAX_ITER, TOL);
            CALLGRIND_TOGGLE_COLLECT;
            secs += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            nb_solved += ca.nb_solved();
            nb_converged += ca.nb_converged();
            total_refactor += static_cast<long>(ca.get_linear_solver_stats().nb_refactorize);
            ++nb_measured;
            if(!trace_path.empty() && step == nb_rows - 1){
                trace_voltages(trace, ca.get_voltages(), ca.converged_mask());
            }
        }
        CALLGRIND_STOP_INSTRUMENTATION;
    } else if(is_ts){
        TimeSeries ts(grid);
        select_batch_family(ts, algo_name, dc);
        RealMat gen_p, load_p, load_q;
        fill_time_series(grid, nb_rows, gen_p, load_p, load_q);
        ts.modify_gen_p(gen_p);
        ts.modify_load_p(load_p);
        ts.modify_load_q(load_q);

        // the wall clock brackets exactly what callgrind collects
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_ZERO_STATS;
        if(flows){
            ts.compute(Vsolved, MAX_ITER, TOL);
            const auto t0 = std::chrono::steady_clock::now();
            CALLGRIND_TOGGLE_COLLECT;
            ts.compute_flows();
            ts.compute_power_flows();
            CALLGRIND_TOGGLE_COLLECT;
            secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        } else {
            const auto t0 = std::chrono::steady_clock::now();
            CALLGRIND_TOGGLE_COLLECT;
            ts.compute(Vsolved, MAX_ITER, TOL);
            CALLGRIND_TOGGLE_COLLECT;
            secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        CALLGRIND_STOP_INSTRUMENTATION;

        nb_solved = ts.nb_solved();
        nb_converged = ts.nb_converged();
        total_refactor = static_cast<long>(ts.get_linear_solver_stats().nb_refactorize);
        nb_measured = nb_rows;
        if(ts.get_status() != 1){ std::cerr << phase << ": the time series did not converge\n"; return 4; }
        if(!trace_path.empty()){
            if(flows) trace_flows(trace, ts.get_flows(), ts.get_power_flows(), ts.converged_mask());
            else trace_voltages(trace, ts.get_voltages(), ts.converged_mask());
        }
    } else {
        ContingencyAnalysis ca(grid);
        select_batch_family(ca, algo_name, dc);
        if(mask) ca.set_handle_disconnected_grid(true);
        ca.add_multiple_n1(first_n1(grid, nb_rows));

        // the wall clock brackets exactly what callgrind collects
        CALLGRIND_START_INSTRUMENTATION;
        CALLGRIND_ZERO_STATS;
        if(flows){
            ca.compute(Vsolved, MAX_ITER, TOL);
            const auto t0 = std::chrono::steady_clock::now();
            CALLGRIND_TOGGLE_COLLECT;
            ca.compute_flows();
            ca.compute_power_flows();
            CALLGRIND_TOGGLE_COLLECT;
            secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        } else {
            const auto t0 = std::chrono::steady_clock::now();
            CALLGRIND_TOGGLE_COLLECT;
            ca.compute(Vsolved, MAX_ITER, TOL);
            CALLGRIND_TOGGLE_COLLECT;
            secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        }
        CALLGRIND_STOP_INSTRUMENTATION;

        nb_solved = ca.nb_solved();
        nb_converged = ca.nb_converged();
        total_refactor = static_cast<long>(ca.get_linear_solver_stats().nb_refactorize);
        nb_measured = static_cast<int>(ca.my_defaults().size());
        if(!trace_path.empty()){
            if(flows) trace_flows(trace, ca.get_flows(), ca.get_power_flows(), ca.converged_mask());
            else trace_voltages(trace, ca.get_voltages(), ca.converged_mask());
        }
    }

    if(!trace_path.empty() && !write_trace(trace_path, trace)){
        std::cerr << "cannot write the trace to '" << trace_path << "'\n";
        return 6;
    }

    // One refactorization per Newton iteration under the default policy, so this
    // is the iteration count a batch does not otherwise report (the "n" solve
    // included). The same line format as profile_cached_pf: ab_wallclock.sh reads
    // the ms/solve figure off it.
    std::cout << phase << ": " << nb_measured << " solves, "
              << nb_solved << " attempted / " << nb_converged << " converged, "
              << (nb_measured > 0 ? static_cast<double>(total_refactor) / nb_measured : 0.)
              << " refactorizations/solve, "
              << (nb_measured > 0 ? secs / nb_measured * 1e3 : 0.)
              << " ms/solve (wall, meaningless under valgrind), peak rss "
              << peak_rss_kb() << " kB\n";
    return 0;
}
