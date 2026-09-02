#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Step 2: reload grids from ".lsb" and run N powerflows -- a C++-side profile.

Everything this script spends time in is C++: ``LSGrid.load_binary`` deserialises
the grid and ``ac_pf`` runs the solve. Python only holds one numpy array of
initial voltages, so timing (or profiling) this script measures lightsim2grid's
C++ side and essentially nothing else. Build the ".lsb" inputs with
benchmark_matpower_to_binary.py first.

Two modes, because they stress completely different code:

  ``--rebuild``  (default) calls ``tell_solver_need_reset()`` before each solve,
                 so every powerflow re-derives the bus labelling, re-assembles
                 Ybus / Sbus, rebuilds the Jacobian sparsity and re-runs the
                 linear solver's symbolic analysis. This is what a grid2op step
                 that changed the topology pays.
  ``--reuse``    lets the solver keep its caches and its symbolic factorisation,
                 which is what a TimeSeries / repeated-solve workload gets.

Besides the wall clock it prints the solver's OWN timers (``get_timers_jacobian``
/ ``get_linear_solver_stats``), which are accumulated inside C++ and split a
solve into the algorithm-side work (dS assembly, Jacobian fill, mismatch,
voltage update) and the linear-solver work. Those are the numbers to look at
first: they say WHICH C++ phase costs what, without needing a profiler at all.

To get function-level attribution on top of that, run this script under
callgrind -- the extension is stripped in a normal install, but exported
template instantiations still resolve, which is enough to name the hot
functions::

    valgrind --tool=callgrind --callgrind-out-file=cg.out \\
             --toggle-collect='*ac_pf*' --cache-sim=no --branch-sim=no \\
             python benchmark_binary_powerflow.py grids/case9241pegase.lsb -n 3
    callgrind_annotate --threshold=90 cg.out

``--toggle-collect`` keeps python start-up, the numpy import and load_binary out
of the profile, so what is collected is the powerflows only.

Usage:  benchmark_binary_powerflow.py <grid.lsb> [more.lsb ...]
                                      [-n N] [--algo NAME] [--reuse|--rebuild]
"""

import argparse
import os
import time

import numpy as np

from lightsim2grid.network import LSGrid


def run_one(path, n_pf, algo, rebuild):
    t0 = time.perf_counter()
    grid = LSGrid.load_binary(path)
    t_load = time.perf_counter() - t0

    grid.change_algorithm(algo)
    v0 = np.full(grid.total_bus(), 1.0 + 0.0j, dtype=complex)

    # one solve outside the measurement: it pays the very first Ybus build and
    # symbolic factorisation, which is not what "N powerflows" is asking about
    if grid.ac_pf(v0, 30, 1e-8).shape[0] == 0:
        print(f"{os.path.basename(path)}: DIVERGED with {algo}")
        return None

    per_pf = []
    for _ in range(n_pf):
        if rebuild:
            grid.tell_solver_need_reset()
        t0 = time.perf_counter()
        v = grid.ac_pf(v0, 30, 1e-8)
        per_pf.append(time.perf_counter() - t0)
        assert v.shape[0] != 0, "diverged mid-benchmark"

    algo_obj = grid.get_algo()
    tj = algo_obj.get_timers_jacobian()
    name = os.path.basename(path).replace(".lsb", "")

    print(f"\n=== {name}: {grid.total_bus()} buses, {algo}, "
          f"{'rebuild every solve' if rebuild else 'solver cache reused'}, "
          f"{algo_obj.get_nb_iter()} NR iterations ===")
    print(f"  load_binary (once)   : {t_load * 1e3:9.2f} ms")
    print(f"  powerflow wall clock : {np.mean(per_pf) * 1e3:9.3f} ms mean, "
          f"{np.min(per_pf) * 1e3:.3f} ms best, over {n_pf} solves")

    algo_rows = [("dS assembly (fill_internal_variables)", tj.timer_dSbus),
                 ("Jacobian fill (fill_J)",                tj.timer_fillJ),
                 ("mismatch / residual",                   tj.timer_mismatch),
                 ("voltage update (apply_step)",           tj.timer_Va_Vm),
                 ("step scaling",                          tj.timer_scale),
                 ("convergence check",                     tj.timer_check),
                 ("pre-processing (Ybus/Sbus, J sparsity)", tj.timer_pre_proc)]
    ls_rows = [("linear solver: analyze + factorize", tj.timer_initialize + tj.timer_factor),
               ("linear solver: refactorize",         tj.timer_refactor),
               ("linear solver: solve",               tj.timer_solve)]
    total = tj.timer_total_nr
    pct = (lambda v: 100.0 * v / total if total > 0 else 0.0)

    print("  --- C++ solver timers, last solve (ms) ---")
    for label, val in algo_rows + ls_rows:
        if val >= 0:
            print(f"    {label:40s} {val * 1e3:9.4f}  {pct(val):5.1f}%")
    algo_side = sum(v for _, v in algo_rows if v > 0)
    ls_side = sum(v for _, v in ls_rows if v > 0)
    print(f"    {'-' * 40} {'-' * 9}")
    print(f"    {'ALGORITHM SIDE (lightsim2grid)':40s} {algo_side * 1e3:9.4f}  {pct(algo_side):5.1f}%")
    print(f"    {'LINEAR SOLVER':40s} {ls_side * 1e3:9.4f}  {pct(ls_side):5.1f}%")
    print(f"    {'total (timer_total_nr)':40s} {total * 1e3:9.4f}")
    return dict(name=name, n_bus=grid.total_bus(), mean=float(np.mean(per_pf)),
                algo_side=algo_side, ls_side=ls_side, total=total)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("grids", nargs="+", help='".lsb" files written by benchmark_matpower_to_binary.py')
    parser.add_argument("-n", "--n-powerflow", type=int, default=10)
    parser.add_argument("--algo", default="NR_KLU")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--rebuild", dest="rebuild", action="store_true", default=True)
    group.add_argument("--reuse", dest="rebuild", action="store_false")
    args = parser.parse_args()

    results = [r for r in (run_one(p, args.n_powerflow, args.algo, args.rebuild)
                           for p in args.grids) if r is not None]
    if len(results) > 1:
        print(f"\n=== summary ({args.algo}, "
              f"{'rebuild' if args.rebuild else 'cache reused'}) ===")
        print(f"  {'grid':18s} {'buses':>7s} {'ms/pf':>9s} {'algo side':>10s} {'lin solver':>11s}")
        for r in results:
            print(f"  {r['name']:18s} {r['n_bus']:7d} {r['mean'] * 1e3:9.3f} "
                  f"{100 * r['algo_side'] / r['total'] if r['total'] else 0:9.1f}% "
                  f"{100 * r['ls_side'] / r['total'] if r['total'] else 0:10.1f}%")


if __name__ == "__main__":
    main()
