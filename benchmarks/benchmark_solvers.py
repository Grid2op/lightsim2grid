# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import os
import warnings
import pandas as pd
from grid2op import make
from grid2op.Backend import PandaPowerBackend
from grid2op.Agent import DoNothingAgent
from grid2op.Chronics import ChangeNothing
import re
from packaging import version
import pandapower
if version.parse(pandapower.__version__) > version.parse("3.0.0"):
    PP_ORIG_FILE = "pandapower_v3"
else:
    PP_ORIG_FILE = "pandapower_v2"
    
try:
    from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance as GridStateFromFile
except ImportError:
    print("Be carefull: there might be maintenance")
    from grid2op.Chronics import GridStateFromFile

try:
    from pypowsybl2grid import PyPowSyBlBackend
    pypowbk_error = None
except ImportError as exc_:
    pypowbk_error = exc_
    print("Backend based on pypowsybl will not be benchmarked")
 
from grid2op.Parameters import Parameters
import lightsim2grid
from lightsim2grid.lightSimBackend import LightSimBackend
from utils_benchmark import run_env, str2bool, get_env_name_displayed, print_configuration
TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")

MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"
DONT_SAVE = "__DONT_SAVE"
NICSLU_LICENSE_AVAIL = os.path.exists("./nicslu.lic") and os.path.isfile("./nicslu.lic")

solver_names = {lightsim2grid.SolverType.GaussSeidel: "GS",
                lightsim2grid.SolverType.GaussSeidelSynch: "GS synch",
                lightsim2grid.SolverType.SparseLU: "NR (SLU)",
                lightsim2grid.SolverType.KLU: "NR (KLU)",
                lightsim2grid.SolverType.NICSLU: "NR (NICSLU *)",
                lightsim2grid.SolverType.CKTSO: "NR (CKTSO *)",
                lightsim2grid.SolverType.SparseLUSingleSlack: "NR single (SLU)",
                lightsim2grid.SolverType.KLUSingleSlack: "NR single (KLU)",
                lightsim2grid.SolverType.NICSLUSingleSlack: "NR single (NICSLU *)",
                lightsim2grid.SolverType.CKTSOSingleSlack: "NR single (CKTSO *)",
                lightsim2grid.SolverType.FDPF_XB_SparseLU: "FDPF XB (SLU)",
                lightsim2grid.SolverType.FDPF_BX_SparseLU: "FDPF BX (SLU)",
                lightsim2grid.SolverType.FDPF_XB_KLU: "FDPF XB (KLU)",
                lightsim2grid.SolverType.FDPF_BX_KLU: "FDPF BX (KLU)",
                lightsim2grid.SolverType.FDPF_XB_NICSLU: "FDPF XB (NICSLU *)",
                lightsim2grid.SolverType.FDPF_BX_NICSLU: "FDPF BX (NICSLU *)",
                lightsim2grid.SolverType.FDPF_XB_CKTSO: "FDPF XB (CKTSO *)",
                lightsim2grid.SolverType.FDPF_BX_CKTSO: "FDPF BX (CKTSO *)",
                # lightsim2grid.SolverType.DC: "LS+DC",
                # lightsim2grid.SolverType.KLUDC: "LS+SLU",
                # lightsim2grid.SolverType.NICSLUDC: "LS+SLU"
                }
solver_gs = {lightsim2grid.SolverType.GaussSeidelSynch, lightsim2grid.SolverType.GaussSeidel}
solver_fdpf = {lightsim2grid.SolverType.FDPF_XB_SparseLU, lightsim2grid.SolverType.FDPF_BX_SparseLU,
               lightsim2grid.SolverType.FDPF_XB_KLU, lightsim2grid.SolverType.FDPF_BX_KLU,
               lightsim2grid.SolverType.FDPF_XB_NICSLU, lightsim2grid.SolverType.FDPF_BX_NICSLU,
               lightsim2grid.SolverType.FDPF_XB_CKTSO, lightsim2grid.SolverType.FDPF_BX_CKTSO,
               }
res_times = {}

order_solver_print = [
    lightsim2grid.SolverType.GaussSeidel,
    lightsim2grid.SolverType.GaussSeidelSynch,
    lightsim2grid.SolverType.SparseLUSingleSlack,
    lightsim2grid.SolverType.SparseLU,
    lightsim2grid.SolverType.KLUSingleSlack,
    lightsim2grid.SolverType.KLU,
    lightsim2grid.SolverType.NICSLUSingleSlack,
    lightsim2grid.SolverType.NICSLU,
    lightsim2grid.SolverType.CKTSOSingleSlack,
    lightsim2grid.SolverType.CKTSO,
    lightsim2grid.SolverType.FDPF_XB_SparseLU,
    lightsim2grid.SolverType.FDPF_BX_SparseLU,
    lightsim2grid.SolverType.FDPF_XB_KLU,
    lightsim2grid.SolverType.FDPF_BX_KLU,
    lightsim2grid.SolverType.FDPF_XB_NICSLU,
    lightsim2grid.SolverType.FDPF_BX_NICSLU,
    lightsim2grid.SolverType.FDPF_XB_CKTSO,
    lightsim2grid.SolverType.FDPF_BX_CKTSO,
]


def main(max_ts,
         env_name_input,
         test=True,
         no_gs=False,
         no_gs_synch=False,
         no_pp=False,
         save_results=DONT_SAVE):
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
    aor_pp = None  # needed in case the user does not want to compute results for pandapower

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if re.match("^.*\\.json$", env_name_input) is None:
            # i provided an environment name
            env_pp = make(env_name_input, param=param, test=test,
                          backend=PandaPowerBackend(lightsim2grid=False, with_numba=True),
                          data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            env_pp_no_numba = make(env_name_input, param=param, test=test,
                                   backend=PandaPowerBackend(lightsim2grid=False, with_numba=False),
                                   data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            env_pp_ls_numba = make(env_name_input, param=param, test=test,
                                   backend=PandaPowerBackend(lightsim2grid=True, with_numba=True),
                                   data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            env_lightsim = make(env_name_input, backend=LightSimBackend(loader_kwargs={"pp_orig_file": PP_ORIG_FILE}), param=param, test=test,
                                data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            if pypowbk_error is None:
                env_pypow = make(env_name_input, param=param, test=test,
                                 backend=PyPowSyBlBackend(),
                                 data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
        else:
            # I provided an environment path
            env_pp = make("blank", param=param, test=True,
                          data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                          grid_path=env_name_input,
                          backend=PandaPowerBackend(lightsim2grid=False, with_numba=True)
                          )
            env_pp_no_numba = make("blank", param=param, test=True,
                                   data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                   grid_path=env_name_input,
                                   backend=PandaPowerBackend(lightsim2grid=False, with_numba=False)
                                   )
            env_pp_ls_numba = make("blank", param=param, test=True,
                                   data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                   grid_path=env_name_input,
                                   backend=PandaPowerBackend(lightsim2grid=True, with_numba=True)
                                   )
            if pypowbk_error is None:
                env_pypow = make("blank", param=param, test=True,
                                 data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                 grid_path=env_name_input,
                                 backend=PyPowSyBlBackend())
            env_lightsim = make("blank", param=param, test=True,
                                backend=LightSimBackend(loader_kwargs={"pp_orig_file": PP_ORIG_FILE}),
                                data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                grid_path=env_name_input)
            _, env_name_input = os.path.split(env_name_input)

    agent = DoNothingAgent(action_space=env_pp.action_space)
    if no_pp is False:
        print("Start using Pandapower")
        nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0, env_seed=0)
        pp_comp_time = env_pp.backend.comp_time
        pp_time_pf = env_pp._time_powerflow
        if hasattr(env_pp, "_time_step"):
            # for oldest grid2op version where this was not stored
            time_pp = env_pp._time_step
        
        tmp_no_numba = run_env(env_pp_no_numba, max_ts, agent, chron_id=0, env_seed=0)
        nb_ts_pp_no_numba, time_pp_no_numba, aor_pp_no_numba, gen_p_pp_no_numba, gen_q_pp_no_numba = tmp_no_numba
        pp_no_numba_comp_time = env_pp_no_numba.backend.comp_time
        pp_no_numba_time_pf = env_pp_no_numba._time_powerflow
        if hasattr(env_pp_no_numba, "_time_step"):
            # for oldest grid2op version where this was not stored
            time_pp_no_numba = env_pp_no_numba._time_step
        
        tmp_ls_numba = run_env(env_pp_ls_numba, max_ts, agent, chron_id=0, env_seed=0)
        nb_ts_pp_ls_numba, time_pp_ls_numba, aor_pp_ls_numba, gen_p_ls_numba, gen_q_ls_numba = tmp_ls_numba
        pp_ls_numba_comp_time = env_pp_ls_numba.backend.comp_time
        pp_ls_numba_time_pf = env_pp_ls_numba._time_powerflow
        if hasattr(env_pp_ls_numba, "_time_step"):
            # for oldest grid2op version where this was not stored
            time_pp_ls_numba = env_pp_ls_numba._time_step

    if pypowbk_error is None:
        # also benchmark pypowsybl backend
        nb_ts_pypow, time_pypow, aor_pypow, gen_p_pypow, gen_q_pypow = run_env(env_pypow, max_ts, agent, chron_id=0, env_seed=0)
        pypow_comp_time = env_pypow.backend.comp_time
        pypow_time_pf = env_pypow._time_powerflow
        if hasattr(env_pypow, "_time_step"):
            # for oldest grid2op version where this was not stored
            time_pypow = env_pypow._time_step
    
    wst = True  # print extra info in the run_env function
    solver_types = env_lightsim.backend.available_solvers
    for solver_type in solver_types:
        if solver_type not in solver_names:
            continue
        print(f"Start using {solver_type}")
        env_lightsim.backend.set_solver_type(solver_type)
        if solver_type in solver_gs:
            # gauss seidel sovler => more iterations
            env_lightsim.backend.set_solver_max_iter(10000)
            if lightsim2grid.SolverType.GaussSeidel == solver_type and no_gs:
                # I don't study the gauss seidel solver
                continue
            elif lightsim2grid.SolverType.GaussSeidelSynch  == solver_type and no_gs_synch:
                # I don't study the gauss seidel synch solver
                continue
        elif solver_type in solver_fdpf:
            # gauss seidel sovler => more iterations
            env_lightsim.backend.set_solver_max_iter(30)
        else:
            # NR based solver => less iterations
            env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                with_type_solver=wst, env_seed=0)
        gs_comp_time = env_lightsim.backend.comp_time
        gs_time_pf = env_lightsim._time_powerflow
        if hasattr(env_lightsim, "_time_step"):
            # for oldest grid2op version where this was not stored
            time_gs = env_lightsim._time_step
        res_times[solver_type] = (solver_names[solver_type],
                                  nb_ts_gs, time_gs, aor_gs, gen_p_gs,
                                  gen_q_gs, gs_comp_time, gs_time_pf)

    # NOW PRINT THE RESULTS
    print("Configuration:")
    config_str = print_configuration(pypowbk_error)
    if save_results != DONT_SAVE:
        with open(save_results+"config_info.txt", "w", encoding="utf-8") as f:
            f.write(config_str)
    # order on which the solvers will be 
    this_order =  [el for el in res_times.keys() if el not in order_solver_print] + order_solver_print

    env_name = get_env_name_displayed(env_name_input)
    hds = [f"{env_name}", "grid2op speed (it/s)", "grid2op 'backend.runpf' time (ms)", "time in 'algo' (ms / pf)"]
    tab = []
    if no_pp is False:
        tab.append(["PP", f"{nb_ts_pp/time_pp:.2e}",
                    f"{1000.*pp_time_pf/nb_ts_pp:.2e}",
                    f"{1000.*pp_comp_time/nb_ts_pp:.2e}"])
        tab.append(["PP (no numba)", f"{nb_ts_pp_no_numba/time_pp_no_numba:.2e}",
                    f"{1000.*pp_no_numba_time_pf/nb_ts_pp_no_numba:.2e}",
                    f"{1000.*pp_no_numba_comp_time/nb_ts_pp_no_numba:.2e}"])
        tab.append(["PP (with lightsim)", f"{nb_ts_pp_ls_numba/time_pp_ls_numba:.2e}",
                    f"{1000.*pp_ls_numba_time_pf/nb_ts_pp_ls_numba:.2e}",
                    f"{1000.*pp_ls_numba_comp_time/nb_ts_pp_ls_numba:.2e}"])
    if pypowbk_error is None:
        tab.append(["pypowsybl", f"{nb_ts_pypow/time_pypow:.2e}",
                    f"{1000.*pypow_time_pf/nb_ts_pypow:.2e}",
                    f"{1000.*pypow_comp_time/nb_ts_pypow:.2e}"])
        
    for key in this_order:
        if key not in res_times:
            continue
        solver_name, nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs, gs_comp_time, gs_time_pf = res_times[key]
        tab.append([solver_name,
                    f"{nb_ts_gs/time_gs:.2e}",
                    f"{1000.*gs_time_pf/nb_ts_gs:.2e}",
                    f"{1000.*gs_comp_time/nb_ts_gs:.2e}"])

    if TABULATE_AVAIL:
        res_use_with_grid2op_1 = tabulate(tab, headers=hds,  tablefmt="rst")
        print(res_use_with_grid2op_1)
    else:
        print(tab)
        
    if save_results != DONT_SAVE:
        dt = pd.DataFrame(tab, columns=hds)
        dt.to_csv(save_results+"speed.csv", index=False, header=True, sep=";")
    print()

    if TABULATE_AVAIL:
        res_github_readme = tabulate(tab, headers=hds,  tablefmt="github")
        print(res_github_readme)
    else:
        print(tab)
    print()

    if aor_pp is not None:
        nb_ts_this_table = res_times[solver_types[0]][1]
        hds = [f"{env_name} ({nb_ts_this_table} iter)", "Δ aor (amps)", "Δ gen_p (MW)", "Δ gen_q (MVAr)"]
        if no_pp is False:
            tab = [["PP (ref)", "0.00", "0.00", "0.00"]]
            
        for key in this_order:
            if key not in res_times:
                continue
            solver_name, nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs, gs_comp_time, gs_time_pf = res_times[key]
            tab.append([solver_name,
                        f"{np.max(np.abs(aor_gs - aor_pp)):.2e}",
                        f"{np.max(np.abs(gen_p_gs - gen_p_pp)):.2e}",
                        f"{np.max(np.abs(gen_q_gs - gen_q_pp)):.2e}"])

        if TABULATE_AVAIL:
            res_use_with_grid2op_2 = tabulate(tab, headers=hds,  tablefmt="rst")
            print(res_use_with_grid2op_2)
        else:
            print(tab)
            
        if save_results != DONT_SAVE:
            dt = pd.DataFrame(tab, columns=hds)
            dt.to_csv(save_results+"diff.csv", index=False, header=True, sep=";")
    print()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark of lightsim with a "do nothing" agent '
                                                 '(compare multiple lightsim solvers with default grid2op backend '
                                                 'PandaPower)')
    parser.add_argument('--env_name', default=ENV_NAME, type=str,
                        help='Environment name to be used for the benchmark.')
    parser.add_argument('--number', type=int, default=MAX_TS,
                        help='Maximum number of time steps for which the benchmark will be run.')
    parser.add_argument('--no_test', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not use \"test=True\" keyword argument when building the grid2op environments'
                             ' for the benchmark (default False: use \"test=True\"  environment)')
    parser.add_argument('--no_gs_synch', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not benchmark gauss seidel (synch) method (default: evaluate it)')
    parser.add_argument('--no_gs', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not benchmark gauss seidel (regular) method (default: evaluate it)')
    parser.add_argument('--no_pp', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not benchmark pandapower method (default: evaluate it)')
    parser.add_argument("--save_results", default=DONT_SAVE, type=str,
                        help='Name of the file in which you want to save the result table')
    args = parser.parse_args()

    max_ts = int(args.number)
    env_name = str(args.env_name)
    test_env = not args.no_test
    main(max_ts,
         env_name,
         test_env,
         no_gs=args.no_gs,
         no_gs_synch=args.no_gs_synch,
         no_pp=args.no_pp,
         save_results=args.save_results)
