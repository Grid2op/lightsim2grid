# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

"""``ContingencyAnalysis.nb_thread`` scaling benchmark: how much does splitting an
n-1 security analysis across 1 to 8 CPU threads actually speed things up, on a
single, reasonably large real-topology grid (``case6515rte``, ~6515 buses)?

Results are reported in ``docs/security_analysis.rst`` (section
"Running the contingencies on multiple threads").
"""

import os
import pandapower as pp
import numpy as np
from tqdm import tqdm

from lightsim2grid import ContingencyAnalysis
from benchmark_grid_size import make_grid2op_env
from utils_benchmark import print_configuration, get_env_name_displayed

try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")
    TABULATE_AVAIL = False

CASE_NAME = "case6515rte.json"
NB_THREADS_RANGE = [1, 2, 3, 4, 5, 6, 7, 8]
MAX_CONT = 1000  # cap the number of n-1 contingencies simulated, for speed


if __name__ == "__main__":
    print_configuration()

    if not os.path.exists(CASE_NAME):
        import pandapower.networks as pn
        case = getattr(pn, os.path.splitext(CASE_NAME)[0])()
        pp.to_json(case, CASE_NAME)

    case = pp.from_json(CASE_NAME)
    pp.runpp(case)  # for slack

    load_p_init = case.load["p_mw"].to_numpy().copy()
    load_q_init = case.load["q_mvar"].to_numpy().copy()
    gen_p_init = case.gen["p_mw"].to_numpy().copy()
    sgen_p_init = case.sgen["p_mw"].to_numpy().copy()

    slack_gens = np.zeros((1, case.ext_grid.shape[0]))
    if "res_ext_grid" in case:
        slack_gens += np.tile(case.res_ext_grid["p_mw"].to_numpy(), (1, 1))
    gen_p_g2op = np.concatenate((gen_p_init, slack_gens.reshape(-1)))

    env_lightsim = make_grid2op_env(case,
                                     CASE_NAME,
                                     load_p_init.reshape((1, -1)),
                                     load_q_init.reshape((1, -1)),
                                     gen_p_g2op.reshape((1, -1)),
                                     sgen_p_init.reshape((1, -1)))

    rows = []
    total_time_1_thread = None
    for nb_threads in tqdm(NB_THREADS_RANGE):
        env_lightsim.reset()
        sa = ContingencyAnalysis(env_lightsim)
        for i in range(env_lightsim.n_line):
            sa.add_single_contingency(i)
            if i >= MAX_CONT:
                break
        sa.init_from_n_powerflow = True
        sa.nb_thread = nb_threads
        sa.get_flows()
        computer_sa = sa.computer
        total_time = computer_sa.total_time() + computer_sa.amps_computation_time()
        nb_solved = computer_sa.nb_solved()
        pf_per_s = nb_solved / total_time if total_time > 0. else float("nan")
        if nb_threads == 1:
            total_time_1_thread = total_time
        speed_up = total_time_1_thread / total_time if total_time > 0. else float("nan")
        rows.append([nb_threads, nb_solved, f"{1000. * total_time:.2f}", f"{pf_per_s:.0f}", f"{speed_up:.2f}x"])
        sa.close()

    env_lightsim.close()

    print(f"\nnb_thread scaling on {get_env_name_displayed(CASE_NAME)} "
          f"({env_lightsim.n_sub} substations, up to {MAX_CONT + 1} n-1 contingencies)")
    if TABULATE_AVAIL:
        print(tabulate(rows,
                        headers=["nb_thread", "nb solved", "time (ms)", "pf / s", "speed-up vs 1 thread"],
                        tablefmt="rst"))
    else:
        print(rows)
