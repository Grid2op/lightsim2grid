// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// Step 1b of the C++-side powerflow profile: an ordinary ".lsb" -> the same grid
// with every NRSystem feature switched on at once.
//
//     benchmark_matpower_to_binary.py   case9241pegase.m  -> case9241pegase.lsb
//     make_exotic_grid                  case9241pegase.lsb -> case9241_exotic.lsb
//     benchmark_binary_powerflow.py     runs the powerflows
//
// Why it exists: NONE of the MATPOWER benchmark cases carries an hvdc line, an
// SVC or a distributed slack --
//
//     case118 / case1354pegase / case2869pegase / case9241pegase:
//         hvdc = 0, svc = 0, slack generators = 1
//
// -- so `NRSystem`'s extensions (MultiSlack, VoltageControl, Hvdc) contribute
// nothing at all on them: the feature half of the Jacobian is empty and any change
// to it measures as exactly zero. Anything touching that code needs a grid that
// actually exercises it, at a size where the cost is visible.
//
// What it adds, on top of whatever grid it is handed:
//   - a distributed slack over every ON generator with a non-zero target p, the
//     weight being |target_p| (LSGrid::update_slack_weights' own rule);
//   - 5 groups of up to 4 electrically close generators (BFS over lines + trafos)
//     jointly regulating one bus -- remote voltage control with reactive sharing;
//   - 5 voltage-mode SVCs, half of them with a slope;
//   - 4 hvdc lines covering both station types (VSC and LCC) in all three
//     combinations, two of them with angle droop.
//
// Two deliberate stress cases: the two droop lines SHARE an end bus, so their
// dP/dtheta Jacobian entries collide on one (P row, theta col) position; and the
// voltage-regulating VSC stations sit on the generator groups' controlled buses,
// so a generator group and an hvdc station co-regulate one bus.
//
// Two constraints the grid has to respect, both enforced by
// LSGrid::fill_voltage_control_solver_data, and both worth knowing before editing
// this file:
//   - every controller of one bus must ask for the SAME voltage, so the setpoint
//     each generator group chooses is propagated to whoever else lands on its bus;
//   - an SVC must be ALONE in its control group ("not supported in v1"), so the
//     SVCs cannot share a regulated bus with the generator groups, nor with each
//     other. They are placed on generator-free buses found by walking out from a
//     group. An SVC co-regulating with generators is simply not expressible today.
//
// Build (needs a built core library):
//     g++ -O2 -DNDEBUG -std=c++17 -I eigen -I src/core -I src/core/powerflow_algorithm \
//         benchmarks/make_exotic_grid.cpp -o make_exotic_grid \
//         -L <build> -llightsim2grid_core -Wl,-rpath,<build>
//
// Run:
//     ./make_exotic_grid case9241pegase.lsb case9241_exotic.lsb

#include <cstdio>
#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include "LSGrid.hpp"
using namespace ls2g;

int main(int argc, char** argv){
    if(argc < 3){
        std::printf("usage: %s <in.lsb> <out.lsb>\n", argv[0]);
        return 2;
    }
    LSGrid grid = LSGrid::load_binary(argv[1]);
    const int nb_bus = static_cast<int>(grid.total_bus());
    const auto & gens = grid.get_generators_as_data();
    const int nb_gen = static_cast<int>(gens.nb());
    std::printf("base: %d buses, %d gens, %d lines, %d trafos\n", nb_bus, nb_gen,
                (int)grid.get_powerlines_as_data().nb(), (int)grid.get_trafos_as_data().nb());

    // ---- adjacency over lines + trafos, for "electrically close" ----------------
    std::vector<std::vector<int> > adj(nb_bus);
    auto add_edges = [&](const auto & c){
        const auto & b1 = c.get_bus_id_side_1(); const auto & b2 = c.get_bus_id_side_2();
        const auto & s1 = c.get_status_side_1(); const auto & s2 = c.get_status_side_2();
        for(int k = 0; k < (int)c.nb(); ++k){
            if(!s1[k] || !s2[k]) continue;
            const int u = b1[k].cast_int(), v = b2[k].cast_int();
            if(u < 0 || v < 0) continue;
            adj[u].push_back(v); adj[v].push_back(u);
        }
    };
    add_edges(grid.get_powerlines_as_data());
    add_edges(grid.get_trafos_as_data());

    // gens per bus (only the ones that are on and voltage-regulating)
    std::vector<std::vector<int> > gen_at(nb_bus);
    std::vector<int> gen_bus(nb_gen, -1);
    for(int g = 0; g < nb_gen; ++g){
        const GenInfo gi(gens, g);
        if(!gi.connected || !gi.voltage_regulator_on) continue;
        gen_bus[g] = gi.bus_id;
        if(gi.bus_id >= 0 && gi.bus_id < nb_bus) gen_at[gi.bus_id].push_back(g);
    }

    // ---- the reference solution -------------------------------------------------
    // Solved BEFORE anything is added; every voltage setpoint below is read off it,
    // so a controlled bus is asked to hold the magnitude it already holds. The
    // setpoint is a property of THE BUS, not of whoever controls it, which is what
    // makes several controllers of one bus agree by construction -- disagreement is
    // the one thing fill_voltage_control_solver_data refuses outright.
    //
    // Not a detail: taking each group's target from its first generator instead (a
    // value belonging to a different bus) asks buses to hold voltages nothing was
    // holding them at, and the Newton-Raphson gives up -- SolverFactor after ~20
    // iterations, or 300 iterations without converging.
    const CplxVect v_flat = CplxVect::Constant((Eigen::Index)nb_bus, cplx_type(1., 0.));
    const CplxVect v_base = grid.ac_pf(v_flat, 30, 1e-8);
    if(v_base.size() == 0){ std::printf("the base grid does not solve\n"); return 1; }
    auto v_set = [&](int bus){ return std::abs(v_base(bus)); };

    // The bus each generator would control: one branch away. Computed up front so
    // the SVCs can avoid them (an SVC must be alone in its control group).
    std::vector<int> ctrl_of_gen(nb_gen, -1);
    std::vector<int> candidates;
    for(int g = 0; g < nb_gen; ++g){
        if(gen_bus[g] < 0 || adj[gen_bus[g]].empty()) continue;
        ctrl_of_gen[g] = adj[gen_bus[g]].front();
        candidates.push_back(g);
    }

    // ---- distributed slack ------------------------------------------------------
    // Every ON generator with a non-zero target p, weight |target_p| -- which is
    // exactly what LSGrid::update_slack_weights does with an all-true mask.
    Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> all_ok(nb_gen);
    all_ok.setConstant(true);
    grid.update_slack_weights(all_ok);
    std::printf("distributed slack: %d generators\n",
                (int)grid.get_generators_as_data().get_slack_bus_id().size());

    // ---- hvdc: a pan-European model carries plenty of cross-border links --------
    // 20 of them, chained across a spread of generator buses so consecutive lines
    // SHARE an end bus -- the case that makes two droop lines' dP/dtheta entries
    // collide on one Jacobian position. Angle droop is only available VSC - VSC
    // (HvdcLineContainer::init rejects anything else), so the type pattern and the
    // droop pattern are one decision.
    std::vector<int> seeds;
    {
        std::vector<int> gb;
        for(int g = 0; g < nb_gen; ++g) if(gen_bus[g] >= 0) gb.push_back(gen_bus[g]);
        std::sort(gb.begin(), gb.end()); gb.erase(std::unique(gb.begin(), gb.end()), gb.end());
        const int want = 12, step = std::max(1, (int)gb.size() / want);
        for(size_t i = 0; i < gb.size() && (int)seeds.size() < want; i += step) seeds.push_back(gb[i]);
    }
    if((int)seeds.size() < 6){ std::printf("not enough generator buses\n"); return 1; }
    const int n_h = 20;
    Eigen::VectorXi h1(n_h), h2(n_h);
    std::vector<int> t1(n_h), t2(n_h);
    std::vector<bool> vr1(n_h), vr2(n_h), droop(n_h);
    RealVect psp(n_h), p0(n_h), kdroop(n_h), vm1(n_h), vm2(n_h);
    for(int k = 0; k < n_h; ++k){
        h1(k) = seeds[k % seeds.size()];
        h2(k) = seeds[(k + 1) % seeds.size()];
        const int combo = k % 4;                // 0,1 = VSC/VSC  2 = VSC/LCC  3 = LCC/LCC
        t1[k] = (combo == 3) ? 1 : 0;           // 0 = VSC, 1 = LCC
        t2[k] = (combo >= 2) ? 1 : 0;
        vr1[k] = (t1[k] == 0); vr2[k] = (t2[k] == 0);   // only a VSC regulates voltage
        droop[k] = (combo < 2);
        psp(k) = 20. + 2. * (k % 7);
        p0(k) = droop[k] ? psp(k) : 0.;
        kdroop(k) = droop[k] ? (3. + (k % 4)) : 0.;
        vm1(k) = v_set(h1(k)); vm2(k) = v_set(h2(k));
    }
    std::vector<int> cmode(n_h, 0);
    RealVect lf = RealVect::Constant(n_h, 0.01);
    RealVect qz = RealVect::Zero(n_h);
    RealVect qmin = RealVect::Constant(n_h, -60.), qmax = RealVect::Constant(n_h, 60.);
    RealVect pf = RealVect::Constant(n_h, 0.95);
    RealVect r = RealVect::Constant(n_h, 0.1), vnom = RealVect::Constant(n_h, 320.);
    RealVect pmax = RealVect::Constant(n_h, 400.);
    grid.init_hvdc_lines(h1, h2, t1, t2, lf, lf, vr1, vr2, vm1, vm2, qz, qz,
                         qmin, qmax, qmin, qmax, pf, pf, cmode, psp, r, vnom,
                         droop, p0, kdroop, pmax, pmax);
    int nb_droop = 0; for(int k = 0; k < n_h; ++k) nb_droop += droop[k];
    std::printf("hvdc: %d lines (%d on droop) over %d buses\n", n_h, nb_droop, (int)seeds.size());

    // ---- 5 SVCs, each ALONE on the bus it regulates -----------------------------
    // NOT sharing a controlled bus with the generators, although that is the
    // interaction one would want: fill_voltage_control_solver_data refuses it --
    // "an SVC may only be ALONE in its control group ... not supported in v1". Two
    // SVCs on one bus is refused for the same reason.
    std::set<int> taken(seeds.begin(), seeds.end());
    for(int g = 0; g < nb_gen; ++g){
        if(gen_bus[g] >= 0) taken.insert(gen_bus[g]);
        if(ctrl_of_gen[g] >= 0) taken.insert(ctrl_of_gen[g]);
    }
    const int n_svc = 5;
    std::vector<int> svc_mode(n_svc, 1);        // SvcContainer::VOLTAGE
    RealVect svc_v(n_svc), svc_q(n_svc), svc_slope(n_svc), svc_bmin(n_svc), svc_bmax(n_svc);
    Eigen::VectorXi svc_reg(n_svc), svc_bus(n_svc);
    for(int i = 0; i < n_svc; ++i){
        int host = -1;
        std::set<int> seen_b; std::queue<int> qb;
        qb.push(seeds[i % seeds.size()]); seen_b.insert(seeds[i % seeds.size()]);
        while(!qb.empty() && host < 0){
            const int b = qb.front(); qb.pop();
            for(int cand : adj[b]){
                if(seen_b.count(cand)) continue;
                seen_b.insert(cand); qb.push(cand);
                if(!taken.count(cand)){ host = cand; break; }
            }
        }
        if(host < 0){ std::printf("no free bus for svc %d\n", i); return 1; }
        taken.insert(host);
        svc_bus(i) = host; svc_reg(i) = host;
        svc_v(i) = v_set(host); svc_q(i) = 0.;
        svc_slope(i) = (i % 2 == 0) ? 0.01 : 0.;
        svc_bmin(i) = -50.; svc_bmax(i) = 50.;
        std::printf("svc %d: bus %d @ %.4f pu, slope %.3f\n", i, host, svc_v(i), svc_slope(i));
    }
    grid.init_svcs(svc_mode, svc_v, svc_q, svc_slope, svc_bmin, svc_bmax, svc_reg, svc_bus);

    // ---- remote voltage control, LAST and verified ------------------------------
    // In a TSO snapshot essentially every generator regulates a bus that is not its
    // own -- typically the HV side of its step-up transformer -- and that shape is
    // what sizes the VoltageControl extension: one custom Q column per controller,
    // one custom voltage row per controlled bus, (units - 1) sharing rows per bus.
    //
    // Applied in chunks, each kept only if the grid STILL SOLVES with everything
    // else already in place. A minority of generator/bus pairs make a system the
    // Newton-Raphson cannot solve, and with 1445 generators some are always
    // included: taking 20 controllers at a time along the generator list, seven
    // windows out of eight converge in 6 iterations and one does not -- so it is the
    // individual pairs, not the count. Rather than guess a criterion, the tool
    // checks, which is also why this runs last: a chunk validated against a grid
    // without the hvdc, the SVCs and the distributed slack proves nothing about the
    // grid that gets saved.
    std::map<int, std::vector<int> > group_of_bus;
    const int CHUNK = 25;
    int kept = 0, dropped = 0;
    for(size_t i = 0; i < candidates.size(); i += CHUNK){
        const size_t hi = std::min(candidates.size(), i + CHUNK);
        std::vector<std::pair<int, real_type> > undo;
        for(size_t j = i; j < hi; ++j){
            const int g = candidates[j];
            undo.emplace_back(g, GenInfo(gens, g).target_vm_pu);
            grid.set_gen_regulated_bus(g, ctrl_of_gen[g]);
            grid.change_v_gen(g, v_set(ctrl_of_gen[g]));
        }
        if(grid.ac_pf(v_flat, 30, 1e-8).size() > 0){
            for(size_t j = i; j < hi; ++j)
                group_of_bus[ctrl_of_gen[candidates[j]]].push_back(candidates[j]);
            kept += (int)undo.size();
        } else {
            for(const auto & u : undo){
                grid.set_gen_regulated_bus(u.first, gen_bus[u.first]);
                grid.change_v_gen(u.first, u.second);
            }
            dropped += (int)undo.size();
        }
    }
    int biggest = 0;
    for(const auto & kv : group_of_bus) biggest = std::max(biggest, (int)kv.second.size());
    std::printf("remote voltage control: %d generators over %d controlled buses "
                "(largest group: %d units; %d dropped as unsolvable)\n",
                kept, (int)group_of_bus.size(), biggest, dropped);

    // ---- the grid that gets written must solve ----------------------------------
    const CplxVect v_final = grid.ac_pf(v_flat, 30, 1e-8);
    if(v_final.size() == 0){ std::printf("FINAL GRID DOES NOT SOLVE -- not saving\n"); return 1; }
    std::printf("final grid solves in %d iterations\n", grid.get_algo().get_nb_iter());

    grid.save_binary(argv[2]);
    std::printf("saved %s\n", argv[2]);
    return 0;
}
