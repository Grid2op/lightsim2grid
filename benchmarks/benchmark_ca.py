# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import warnings
import copy
import pandapower as pp
import numpy as np        
import hashlib
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from grid2op import make, Parameters
from grid2op.Chronics import FromNPY
from grid2op.Backend import PandaPowerBackend
from lightsim2grid import LightSimBackend, TimeSerie
try:
    from lightsim2grid import ContingencyAnalysis
except ImportError:
    from lightsim2grid import SecurityAnalysis as ContingencyAnalysis
    
from lightsim2grid.solver import SolverType

from tqdm import tqdm
import os
from utils_benchmark import print_configuration, get_env_name_displayed
from benchmark_solvers import solver_names

try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")
    TABULATE_AVAIL = False

VERBOSE = True

from benchmark_grid_size import (
    make_grid2op_env
)

case_names = [
    "case14.json",
              "case14.json",
              "case118.json",
              "case_illinois200.json",
              "case300.json",
              "case1354pegase.json",
              "case1888rte.json",
              # "GBnetwork.json",  # 2224 buses
              "case2848rte.json",
              "case2869pegase.json",
              "case3120sp.json",
              "case6495rte.json",
              "case6515rte.json",
              "case9241pegase.json"
              ]

MAX_CONT = 1000  # cap the number of n-1 contingencies simulated, for speed


def _set_benchmark_limits(grid):
    """Set voltage / current limits on the grid: loose enough that they are not expected
    to be violated (so this does not change convergence behaviour), but present, so that
    `compute_limit_violations=True` actually exercises its per-element checks instead of
    early-returning on empty limit vectors."""
    vn_kv = grid.get_bus_vn_kv()
    grid.set_bus_voltage_limits(0.5 * vn_kv, 1.5 * vn_kv)
    n_line = len(grid.get_lines())
    n_trafo = len(grid.get_trafos())
    grid.set_line_current_limit_side1(np.full(n_line, 1e4))
    grid.set_line_current_limit_side2(np.full(n_line, 1e4))
    grid.set_trafo_current_limit_side1(np.full(n_trafo, 1e4))
    grid.set_trafo_current_limit_side2(np.full(n_trafo, 1e4))


def _run_one_config(env_lightsim, *, init_from_n_powerflow=False,
                     handle_disconnected_grid=False, compute_limit_violations=False,
                     nb_thread=1):
    """Run all n-1 contingencies (capped at MAX_CONT) for one (env, configuration) pair and
    return (pf / s, nb_solved) so the caller can report timings.

    .. note::
        `sa.close()` also resets the underlying `computer`'s counters, so the values needed
        by the caller are extracted *before* closing, not returned as the (now-reset) computer
        object itself.
    """
    env_lightsim.reset()
    if compute_limit_violations:
        _set_benchmark_limits(env_lightsim.backend._grid)
    sa = ContingencyAnalysis(env_lightsim, compute_limit_violations=compute_limit_violations)
    for i in range(env_lightsim.n_line):
        sa.add_single_contingency(i)
        if i >= MAX_CONT:
            break
    sa.init_from_n_powerflow = init_from_n_powerflow
    sa.handle_disconnected_grid = handle_disconnected_grid
    sa.nb_thread = nb_thread
    sa.get_flows()
    computer_sa = sa.computer
    total_time = computer_sa.total_time() + computer_sa.amps_computation_time()
    nb_solved = computer_sa.nb_solved()
    pf_per_s = nb_solved / total_time if total_time > 0. else float("nan")
    sa.close()
    return pf_per_s, nb_solved


# the 4 configurations compared below, each toggling exactly one flag away from the
# defaults (`ContingencyAnalysis.get_flows()` with none of these options set)
FEATURE_CONFIGS = [
    ("default", dict()),
    ("init_from_n_powerflow", dict(init_from_n_powerflow=True)),
    ("handle_disconnected_grid", dict(handle_disconnected_grid=True)),
    ("compute_limit_violations", dict(compute_limit_violations=True)),
]


if __name__ == "__main__":
    tab_features = []  # one row per grid size, one column per entry in FEATURE_CONFIGS

    for case_name in tqdm(case_names):

        if not os.path.exists(case_name):
            import pandapower.networks as pn
            case = getattr(pn, os.path.splitext(case_name)[0])()
            pp.to_json(case, case_name)

        # load the case file
        case = pp.from_json(case_name)
        pp.runpp(case)  # for slack
        
        # extract reference data
        load_p_init = case.load["p_mw"].to_numpy().copy()
        load_q_init = case.load["q_mvar"].to_numpy().copy()
        gen_p_init = case.gen["p_mw"].to_numpy().copy()
        sgen_p_init = case.sgen["p_mw"].to_numpy().copy()
        
        nb_ts = 1
        # add slack !
        slack_gens =  np.zeros((nb_ts, case.ext_grid.shape[0]))
        if "res_ext_grid" in case:
            slack_gens += np.tile(case.res_ext_grid["p_mw"].to_numpy(), (nb_ts, 1))
        gen_p_g2op = np.concatenate((gen_p_init, slack_gens.reshape(-1)))  
        
        env_lightsim = make_grid2op_env(case,
                                case_name,
                                load_p_init.reshape((1,-1)),
                                load_q_init.reshape((1,-1)),
                                gen_p_g2op.reshape((1,-1)),
                                sgen_p_init.reshape((1,-1)))
        cases_by_threads = {}
        for nb_threads in [1, 2, 3, 4, 5, 6, 7, 8]:
            env_lightsim.reset()
            sa = ContingencyAnalysis(env_lightsim)
            for i in range(env_lightsim.n_line):
                sa.add_single_contingency(i)
                if i >= 1000:
                    break
            sa.init_from_n_powerflow = True
            sa.nb_thread = nb_threads
            p_or, a_or, voltages = sa.get_flows()
            computer_sa = sa.computer
            res_time = 1.
            res_unit = "s"
            
            if VERBOSE:
                # print detailed results if needed
                print("=====================================================")
                print(f"For environment: {case_name} ({env_lightsim.n_sub} substations) [{computer_sa.nb_solved()} powerflows], using {nb_threads} threads")
                # if nb_threads == 1:
                print(f"Total time spent in \"computer\" to solve everything: {res_time*computer_sa.total_time():.2f}{res_unit} "
                    f"({computer_sa.nb_solved() / computer_sa.total_time():.0f} pf / s), "
                    f"{1000.*computer_sa.total_time() / computer_sa.nb_solved():.2f} ms / pf)")

                total_time = computer_sa.total_time() + computer_sa.amps_computation_time()
                print(f"Compute time: {computer_sa.solve_time():.2e} / {total_time:.2e} "
                      f"({100. * computer_sa.solve_time() / total_time:.1f} % of total time)")
                print(f"Compute time (per thread): {computer_sa.solve_time() / nb_threads:.2e} / {total_time:.2e} "
                      f"({100. * computer_sa.solve_time() / nb_threads / total_time:.1f} % of total time per thread)")
                print(f"\tExtra threading time: {computer_sa.thread_init_time():.2e} / {total_time:.2e} "
                      f"({100. * computer_sa.thread_init_time() / total_time:.1f} % of total time)")
                print(f"\tExtra pre proc time: {computer_sa.preprocessing_time():.2e} / {total_time:.2e} "
                      f"({100. * computer_sa.preprocessing_time() / total_time:.1f} % of total time)")
                print(f"\tExtra post proc time: {computer_sa.amps_computation_time():.2e} / {total_time:.2e} "
                      f"({100. * computer_sa.amps_computation_time() / total_time:.1f} % of total time)")
                print("=====================================================\n")
                
            # sa_times.append(computer_sa.total_time() + computer_sa.amps_computation_time())
            # sa_speeds.append(computer_sa.nb_solved() / (computer_sa.total_time() + computer_sa.amps_computation_time()) )
            # sa_sizes.append(env_lightsim.n_sub)
            cases_by_threads[nb_threads] = computer_sa.total_time() + computer_sa.amps_computation_time()
            # close the env
            linear_solver_used_str = solver_names[env_lightsim.backend._grid.get_solver_type()]
        for k, v in cases_by_threads.items():
            print(f"{case_name}, nb_threads={k}: {cases_by_threads[1] / v:.1f}x speed-up vs 1 thread")

        # compare the "default", "init_from_n_powerflow", "handle_disconnected_grid" and
        # "compute_limit_violations" options against one another (single thread, so the
        # comparison isolates the cost of each feature from thread-scaling effects above)
        row = [get_env_name_displayed(case_name)]
        pf_per_s_default = None
        for config_name, kwargs in FEATURE_CONFIGS:
            pf_per_s, nb_solved = _run_one_config(env_lightsim, **kwargs)
            if config_name == "default":
                pf_per_s_default = pf_per_s
            if VERBOSE:
                print(f"[{case_name}] {config_name}: {nb_solved} pf solved, "
                      f"{pf_per_s:.0f} pf / s ({1000. / pf_per_s:.3f} ms / pf)")
            speed_ratio = pf_per_s / pf_per_s_default if pf_per_s_default else float("nan")
            row.append(f"{pf_per_s:.0f} pf/s ({speed_ratio:.2f}x default)")
        tab_features.append(row)

        env_lightsim.close()

    print("\nComparison of `default` / `init_from_n_powerflow` / `handle_disconnected_grid` / "
          "`compute_limit_violations` (single thread, up to "
          f"{MAX_CONT + 1} n-1 contingencies per grid)")
    if TABULATE_AVAIL:
        print(tabulate(tab_features,
                        headers=["grid"] + [name for name, _ in FEATURE_CONFIGS],
                        tablefmt="rst"))
    else:
        print(tab_features)