# Copyright (c) 2024, RTE (https://www.rte-france.com)
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
from grid2op.Chronics import FromNPY, ChangeNothing
from grid2op.Backend import PandaPowerBackend
from grid2op.Exceptions import Grid2OpException
import lightsim2grid
from lightsim2grid import LightSimBackend
from benchmark_grid_size import (get_loads_gens,
                                 make_grid2op_env_pp,
                                 run_grid2op_env,
                                 make_grid2op_env)
from benchmark_solvers import solver_gs, solver_names, order_solver_print
    
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
    
VERBOSE = False
MAKE_PLOT = False
WITH_PP = False
DEBUG = False

case_names = [
              "case14.json",
              "case118.json",
              "case_illinois200.json",
              "case300.json",
              "case1354pegase.json",
              "case1888rte.json",
            # #   "GBnetwork.json",  # 2224 buses
            #   "case2848rte.json",
            #   "case2869pegase.json",
            #   "case3120sp.json",
            #   "case6495rte.json",
            #   "case6515rte.json",
            #   "case9241pegase.json"
              ]


if __name__ == "__main__":
    prng = np.random.default_rng(42)
    case_names_displayed = [get_env_name_displayed(el) for el in case_names]
    nb_iters = []
    ts_sizes = []
    errors = {}
    for case_name in tqdm(case_names):

        if not os.path.exists(case_name):
            import pandapower.networks as pn
            case = getattr(pn, os.path.splitext(case_name)[0])()
            pp.to_json(case, case_name)

        # load the case file
        case = pp.from_json(case_name)
        ts_sizes.append(case.bus.shape[0])
        pp.runpp(case)  # for slack
        
        # create the env
        param = Parameters.Parameters()
        param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
            
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env_pp = make("blank",
                        param=param, test=True,
                        backend=PandaPowerBackend(lightsim2grid=False),
                        chronics_class=ChangeNothing,
                        grid_path=case_name,
                        _add_to_name=f"{case_name}",
                        )
            env_ls = make("blank",
                        param=param, test=True,
                        backend=LightSimBackend(),
                        chronics_class=ChangeNothing,
                        grid_path=case_name,
                        _add_to_name=f"{case_name}",
                        )
        env_ls.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidel)
        all_iters = [1, 3, 10, 30, 100, 300, 1_000, 3_000,
                     10_000, 30_000,
                     100_000, 300_000
                    ]
        iters = []
        errors_p = []
        errors_q = []
        for max_iter in all_iters:
            env_ls.backend.set_solver_max_iter(max_iter)
            env_ls.backend._grid.tell_solver_need_reset()
            conv = True
            try:
                obs = env_ls.reset()
            except Grid2OpException as exc_:
                conv = False
            iters.append(env_ls.backend._grid.get_solver().get_nb_iter())
            v_tmp = env_ls.backend._grid.get_solver().get_V()
            res_tmp = env_ls.backend._grid.check_solution(v_tmp, False)
            error_p = 1. * np.abs(res_tmp.real).max()
            error_q = 1. * np.abs(res_tmp.imag).max()
            errors_p.append(error_p)
            errors_q.append(error_q)
            if conv:
                break
        if conv:
            nb_iters.append(iters[-1])
        else:
            nb_iters.append(None)
            
        errors[case.bus.shape[0]] = (errors_p, errors_q)
        
    print("Configuration:")
    print_configuration()
    print(f"Solver used for linear algebra: {lightsim2grid.SolverType.GaussSeidel}")
    print()
    hds = ["grid size (nb bus)", "gauss seidel max iter"]
    tab = []
    for sz, nb_it in zip(ts_sizes, nb_iters):
        tab.append([sz, nb_it])
        
    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab, headers=hds,  tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab)
        
    print(errors[118][0])
    print(errors[118][1])
    import pickle
    with open("res_gauss_seidel.pickle", "wb") as f:
        pickle.dump(errors, file=f)
    with open("res_gauss_seidel_nb_iters.pickle", "wb") as f:
        pickle.dump(nb_iters, file=f)
    print()
    print()


# total computation time : 1h27min16s
# Configuration:

# - date: 2024-12-02 18:46  CET
# - system: Linux 5.15.0-56-generic
# - OS: ubuntu 20.04
# - processor: Intel(R) Core(TM) i7-6820HQ CPU @ 2.70GHz
# - python version: 3.8.10.final.0 (64 bit)
# - numpy version: 1.24.3
# - pandas version: 2.0.3
# - pandapower version: 2.14.0
# - grid2op version: 1.11.0.dev2
# - lightsim2grid version: 0.9.2.post2
# - lightsim2grid extra information: 

# 	- klu_solver_available: True 
# 	- nicslu_solver_available: False 
# 	- cktso_solver_available: False 
# 	- compiled_march_native: False 
# 	- compiled_o3_optim: False 

# Solver used for linear algebra: SolverType.GaussSeidel

# ====================  =======================
#   grid size (nb bus)    gauss seidel max iter
# ====================  =======================
#                   14                      278
#                  118                     3274
#                  200                     8360
#                  300                    40783
#                 1354                   122169
#                 1888
# ====================  =======================
# [31.858705410410803, 13.801689961508492, 7.912199121114395, 6.387621207822959, 4.5494311573542525, 1.3539274305627065, 0.01652457790687702, 5.5928201247405206e-08, 9.957519963773673e-09]
# [111.7637849724719, 52.1105433668106, 6.3902552555152345, 1.1851759157023143, 0.8457897295792693, 0.25197455746676584, 0.0030761171444685202, 1.0415372012959338e-08, 1.8561325626140559e-09]
