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

    // ---- 5 groups of 3-4 close generators, all regulating the seed's bus -------
    std::vector<int> seeds;               // one bus per group: the controlled bus
    std::vector<std::vector<int> > groups;
    std::set<int> used_gen;
    const int stride = nb_gen / 7;        // spread the seeds across the grid
    for(int s = 0; s < nb_gen && (int)groups.size() < 5; s += stride){
        if(gen_bus[s] < 0 || used_gen.count(s)) continue;
        // BFS from the seed bus, collecting other unused generators within 3 hops
        std::vector<int> members{s};
        std::set<int> seen{gen_bus[s]};
        std::queue<std::pair<int,int> > q; q.push({gen_bus[s], 0});
        while(!q.empty() && (int)members.size() < 4){
            auto [b, d] = q.front(); q.pop();
            if(d >= 3) continue;
            for(int nb_ : adj[b]){
                if(seen.count(nb_)) continue; seen.insert(nb_);
                for(int g2 : gen_at[nb_]){
                    if(used_gen.count(g2) || g2 == s) continue;
                    members.push_back(g2);
                    if((int)members.size() >= 4) break;
                }
                q.push({nb_, d + 1});
                if((int)members.size() >= 4) break;
            }
        }
        if((int)members.size() < 3) continue;      // want 3-4, skip a lonely seed
        for(int g : members) used_gen.insert(g);
        seeds.push_back(gen_bus[s]);
        groups.push_back(members);
    }
    // point every member at the seed bus, with one common voltage target.
    // Every OTHER controller that ends up on one of these buses -- an SVC, a
    // voltage-regulating VSC station -- must agree on that setpoint, or
    // fill_voltage_control_solver_data rejects the grid: one controlled bus has one
    // voltage. So the chosen value is recorded per bus and reused below.
    std::map<int, real_type> v_of_bus;
    for(size_t g = 0; g < groups.size(); ++g){
        const real_type v_target = GenInfo(gens, groups[g][0]).target_vm_pu;
        v_of_bus[seeds[g]] = v_target;
        for(int gid : groups[g]){
            grid.set_gen_regulated_bus(gid, seeds[g]);
            grid.change_v_gen(gid, v_target);
        }
        std::printf("group %zu: %d gens -> bus %d (v=%.4f)\n", g, (int)groups[g].size(), seeds[g], v_target);
    }

    // ---- 5 SVCs, each ALONE on the bus it regulates ----------------------------
    // NOT sharing a controlled bus with the generator groups, although that is the
    // interaction one would want to exercise: LSGrid::fill_voltage_control_solver_data
    // refuses it outright -- "an SVC may only be ALONE in its control group ... not
    // supported in v1". Two SVCs on one bus is refused for the same reason. So each
    // SVC gets a bus of its own that carries no generator, picked next to a group so
    // it still sits in the electrically interesting part of the grid.
    //
    // The generator + HVDC-station case IS supported and is exercised below: the VSC
    // stations sit on the group buses and regulate them alongside the generators.
    const int n_svc = 5;
    std::vector<int> svc_mode(n_svc, 1);   // SvcContainer::VOLTAGE
    RealVect svc_v(n_svc), svc_q(n_svc), svc_slope(n_svc), svc_bmin(n_svc), svc_bmax(n_svc);
    Eigen::VectorXi svc_reg(n_svc), svc_bus(n_svc);
    std::set<int> taken(seeds.begin(), seeds.end());
    for(int i = 0; i < n_svc; ++i){
        // BFS out from a group bus until a bus turns up that no voltage controller
        // owns and that carries no generator of its own
        int host = -1;
        {
            std::set<int> seen_b; std::queue<int> qb;
            qb.push(seeds[i % seeds.size()]); seen_b.insert(seeds[i % seeds.size()]);
            while(!qb.empty() && host < 0){
                const int b = qb.front(); qb.pop();
                for(int cand : adj[b]){
                    if(seen_b.count(cand)) continue;
                    seen_b.insert(cand); qb.push(cand);
                    if(gen_at[cand].empty() && !taken.count(cand)){ host = cand; break; }
                }
            }
        }
        if(host < 0){ std::printf("no free bus for svc %d\n", i); return 1; }
        taken.insert(host);
        svc_bus(i) = host; svc_reg(i) = host;               // local control, alone
        svc_v(i) = 1.02; svc_q(i) = 0.;
        svc_slope(i) = (i % 2 == 0) ? 0.01 : 0.;            // some sloped, some not
        svc_bmin(i) = -50.; svc_bmax(i) = 50.;
        std::printf("svc %d: at bus %d, regulates bus %d @ %.4f pu, slope %.3f\n",
                    i, svc_bus(i), svc_reg(i), svc_v(i), svc_slope(i));
    }
    grid.init_svcs(svc_mode, svc_v, svc_q, svc_slope, svc_bmin, svc_bmax, svc_reg, svc_bus);

    // ---- 4 HVDC lines ----------------------------------------------------------
    const int n_h = 4;
    Eigen::VectorXi h1(n_h), h2(n_h);
    // lines 0 and 1 share end bus seeds[0]: two droop lines on one bus, so their
    // dP/dtheta feature entries collide on the same Jacobian position
    h1 << seeds[0], seeds[0], seeds[2], seeds[3];
    h2 << seeds[1], seeds[2], seeds[3], seeds[4];
    std::vector<int> t1{0, 0, 0, 1}, t2{0, 0, 1, 1};        // 0 = VSC, 1 = LCC
    std::vector<bool> vr1{true, true, true, false}, vr2{true, true, false, false};
    std::vector<bool> droop{true, true, false, false};
    std::vector<int> cmode(n_h, 0);
    RealVect lf1 = RealVect::Constant(n_h, 0.01), lf2 = RealVect::Constant(n_h, 0.01);
    RealVect vm1(n_h), vm2(n_h);   // same rule: match whoever else controls that bus
    for(int k = 0; k < n_h; ++k){
        vm1(k) = v_of_bus.count(h1(k)) ? v_of_bus[h1(k)] : 1.02;
        vm2(k) = v_of_bus.count(h2(k)) ? v_of_bus[h2(k)] : 1.02;
    }
    RealVect q1 = RealVect::Zero(n_h), q2 = RealVect::Zero(n_h);
    RealVect qmin = RealVect::Constant(n_h, -60.), qmax = RealVect::Constant(n_h, 60.);
    RealVect pf = RealVect::Constant(n_h, 0.95);
    RealVect psp(n_h); psp << 40., 30., 25., 20.;
    RealVect r = RealVect::Constant(n_h, 0.1), vnom = RealVect::Constant(n_h, 320.);
    RealVect p0(n_h); p0 << 40., 30., 0., 0.;
    RealVect kdroop(n_h); kdroop << 5., 4., 0., 0.;
    RealVect pmax = RealVect::Constant(n_h, 200.);
    grid.init_hvdc_lines(h1, h2, t1, t2, lf1, lf2, vr1, vr2, vm1, vm2, q1, q2,
                         qmin, qmax, qmin, qmax, pf, pf, cmode, psp, r, vnom,
                         droop, p0, kdroop, pmax, pmax);
    for(int k = 0; k < n_h; ++k)
        std::printf("hvdc %d: bus %d (type %d) -> bus %d (type %d), droop %d\n",
                    k, h1(k), t1[k], h2(k), t2[k], (int)droop[k]);

    // ---- distributed slack: every ON generator with p != 0, weight = |p| -------
    Eigen::Array<bool, Eigen::Dynamic, Eigen::RowMajor> all_ok(nb_gen);
    all_ok.setConstant(true);
    grid.update_slack_weights(all_ok);
    std::printf("slack gens: %d\n", (int)grid.get_generators_as_data().get_slack_bus_id().size());

    grid.save_binary(argv[2]);
    std::printf("saved %s\n", argv[2]);
    return 0;
}
