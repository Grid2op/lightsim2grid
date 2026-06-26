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



if __name__ == "__main__":
    
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
        env_lightsim.close()
        break