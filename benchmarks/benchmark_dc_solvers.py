# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import time
import numpy as np
import os
import warnings
import pandas as pd
import re

from packaging import version
import pandapower
if version.parse(pandapower.__version__) > version.parse("3.0.0"):
    PP_ORIG_FILE = "pandapower_v3"
else:
    PP_ORIG_FILE = "pandapower_v2"

from grid2op import make
from grid2op.Exceptions import Grid2OpException
from grid2op.Backend import PandaPowerBackend
from grid2op.Agent import DoNothingAgent
from grid2op.Chronics import ChangeNothing
from grid2op.Environment import MultiMixEnvironment
try:
    from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance as GridStateFromFile
except ImportError:
    print("Be carefull: there might be maintenance")
    from grid2op.Chronics import GridStateFromFile
from grid2op.Parameters import Parameters
from grid2op.dtypes import dt_float

import lightsim2grid
from lightsim2grid import solver
from lightsim2grid import LightSimBackend, TimeSerie
try:
    from lightsim2grid import ContingencyAnalysis
except ImportError:
    from lightsim2grid import SecurityAnalysis as ContingencyAnalysis
    
from utils_benchmark import print_res, run_env, str2bool, get_env_name_displayed, print_configuration
TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")
    
try:
    from pypowsybl2grid import PyPowSyBlBackend
    PYPOW_ERROR = None
except ImportError as exc_:
    PYPOW_ERROR = exc_
    print("Backend based on pypowsybl will not be benchmarked")
    
MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"
DONT_SAVE = "__DONT_SAVE"
NICSLU_LICENSE_AVAIL = os.path.exists("./nicslu.lic") and os.path.isfile("./nicslu.lic")

solver_names = {lightsim2grid.SolverType.DC: "DC",
                lightsim2grid.SolverType.KLUDC: "DC (KLU)",
                lightsim2grid.SolverType.NICSLUDC: "DC (NICSLU *)",
                lightsim2grid.SolverType.CKTSODC: "DC (CKTSO *)"
                }
solver_gs = {}
solver_fdpf = {}
res_times = {}

order_solver_print = [
    lightsim2grid.SolverType.DC,
    lightsim2grid.SolverType.KLUDC,
    lightsim2grid.SolverType.NICSLUDC,
    lightsim2grid.SolverType.CKTSODC,
    
]


def main(max_ts,
         env_name_input,
         test=True,
         no_gs=False,
         no_gs_synch=False,
         no_pp=False,
         save_results=DONT_SAVE):
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True, "ENV_DC": True, "FORECAST_DC": True})
    pypow_error = PYPOW_ERROR
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if re.match("^.*\\.json$", env_name_input) is None:
            # i provided an environment name
            env_pp = make(env_name_input, param=param, test=test,
                          backend=PandaPowerBackend(lightsim2grid=False, with_numba=True),
                          data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            env_lightsim = make(env_name_input, backend=LightSimBackend(loader_kwargs={"pp_orig_file": PP_ORIG_FILE}), param=param, test=test,
                                data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            if pypow_error is None:
                try:
                    env_pypow = make(env_name_input, param=param, test=test,
                                    backend=PyPowSyBlBackend(),
                                    data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
                except Grid2OpException as exc_:
                    pypow_error = exc_
        else:
            # I provided an environment path
            env_pp = make("blank", param=param, test=True,
                          data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                          grid_path=env_name_input
                          )
            env_lightsim = make("blank", param=param, test=True,
                                backend=LightSimBackend(loader_kwargs={"pp_orig_file": PP_ORIG_FILE}),
                                data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                grid_path=env_name_input)
            if pypow_error is None:
                try:
                    env_pypow = make("blank", param=param, test=True,
                                    data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                    grid_path=env_name_input,
                                    backend=PyPowSyBlBackend())
                except Grid2OpException as exc_:
                    pypow_error = exc_
            _, env_name_input = os.path.split(env_name_input)

    agent = DoNothingAgent(action_space=env_pp.action_space)
    if no_pp is False:
        print("Start using Pandapower")
        nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0, env_seed=0)
        pp_comp_time = env_pp.backend.comp_time
        pp_time_pf = env_pp._time_powerflow
    
    if pypow_error is None:
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
        res_times[solver_type] = (solver_names[solver_type],
                                  nb_ts_gs, time_gs, aor_gs, gen_p_gs,
                                  gen_q_gs, gs_comp_time, gs_time_pf)

    env_name = get_env_name_displayed(env_name_input)

    real_env_ls = env_lightsim
    if isinstance(real_env_ls, MultiMixEnvironment):
        # get the first (in alphabetical order) env in case of multimix
        real_env_ls = real_env_ls[next(iter(sorted(real_env_ls.keys())))]
        
    # Perform the computation using TimeSerie
    load_p = 1. * real_env_ls.chronics_handler.real_data.data.load_p[:nb_ts_gs]
    load_q = 1. * real_env_ls.chronics_handler.real_data.data.load_q[:nb_ts_gs]
    prod_p = 1. * real_env_ls.chronics_handler.real_data.data.prod_p[:nb_ts_gs]
    time_serie = TimeSerie(real_env_ls)
    computer_ts = time_serie.computer
    computer_ts.change_solver(lightsim2grid.SolverType.KLUDC)
    v_init = real_env_ls.backend.V
    status = computer_ts.compute_Vs(prod_p,
                                    np.zeros((nb_ts_gs, 0), dtype=dt_float),
                                    load_p,
                                    load_q,
                                    v_init,
                                    real_env_ls.backend.max_it,
                                    real_env_ls.backend.tol)
    time_serie._TimeSerie__computed = True
    a_or = time_serie.compute_A()
    p_or = time_serie.compute_P()
    assert status, f"some powerflow diverge for Time Series for {env_name}: {computer_ts.nb_solved()} "
    ts_time = 1e3 * (computer_ts.total_time() + computer_ts.amps_computation_time()) / computer_ts.nb_solved()
    ts_algo_time = 1e3 * (computer_ts.solver_time()) / computer_ts.nb_solved()        
    
    # perform the computation using PTDF
    obs = real_env_ls.reset()
    load_bus = real_env_ls.local_bus_to_global(obs.load_bus, obs.load_to_subid)
    gen_bus = real_env_ls.local_bus_to_global(obs.gen_bus, obs.gen_to_subid)
    Sbus = np.zeros((nb_ts_gs,  real_env_ls.backend._grid.total_bus()), dtype=float)
    Sbus[:, load_bus] -= load_p
    Sbus[:, gen_bus] += prod_p
    T_Sbus = 1. * Sbus.T
    
    PTDF_ = 1.0 * real_env_ls.backend._grid.get_ptdf()
    beg_ = time.perf_counter()
    flows = np.dot(PTDF_, T_Sbus).T
    end_ = time.perf_counter()
    time_only_ptdf, *_ = real_env_ls.backend._grid.get_dc_solver().get_timers_ptdf_lodf()
    time_only_ptdf *= 1000. / Sbus.shape[0]
    time_flow_ptdf = 1000. * (end_ - beg_) / Sbus.shape[0]
    # ts_time = 1e3 * (computer_ts.total_time() + computer_ts.amps_computation_time()) / computer_ts.nb_solved()
    # ts_algo_time = 1e3 * (computer_ts.solver_time()) / computer_ts.nb_solved()  
    
    # Perform a securtiy analysis (up to 1000 contingencies)
    real_env_ls.reset()
    sa = ContingencyAnalysis(real_env_ls)
    computer_sa = sa.computer
    computer_sa.change_solver(lightsim2grid.SolverType.KLUDC)
    for i in range(real_env_ls.n_line):
        sa.add_single_contingency(i)
        if i >= 1000:
            break
    p_or_sa, a_or_sa, voltages = sa.get_flows()
    sa_time = 1e3 * (computer_sa.total_time() + computer_sa.amps_computation_time()) / computer_sa.nb_solved() 
    sa_algo_time = 1e3 * (computer_sa.solver_time()) / computer_sa.nb_solved()    
    
    # perform the computation using LODF
    init_powerflow = p_or[0]
    init_powerflow_diag = np.diag(init_powerflow)
    LODF_mat = 1.0 * real_env_ls.backend._grid.get_lodf()
    beg_ = time.perf_counter()
    por_lodf = init_powerflow + LODF_mat * init_powerflow_diag
    end_ = time.perf_counter()
    _, time_only_lodf, _ = real_env_ls.backend._grid.get_dc_solver().get_timers_ptdf_lodf()
    time_only_lodf *= 1000. / init_powerflow.shape[0]
    time_flow_lodf = 1000. * (end_ - beg_) / init_powerflow.shape[0]
    
    # NOW PRINT THE RESULTS
    print("Configuration:")
    config_str = print_configuration()
    if save_results != DONT_SAVE:
        with open(save_results+"config_info.txt", "w", encoding="utf-8") as f:
            f.write(config_str)
    # order on which the solvers will be 
    this_order =  [el for el in res_times.keys() if el not in order_solver_print] + order_solver_print

    hds = [f"{env_name}", f"grid2op speed (it/s)", f"grid2op 'backend.runpf' time (ms / pf)", f"time in 'algo' (ms / pf)"]
    tab = []
    if no_pp is False:
        tab.append(["PP DC", f"{nb_ts_pp/time_pp:.2e}",
                    f"{1000.*pp_time_pf/nb_ts_pp:.2e}",
                    f"{1000.*pp_comp_time/nb_ts_pp:.2e}"])
        
        
    if pypow_error is None:
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
    tab.append(("time serie **", None, ts_time, ts_algo_time))
    tab.append(("PTDF **", None, time_only_ptdf + time_flow_ptdf, time_flow_ptdf))
    tab.append(("contingency analysis ***", None, sa_time, sa_algo_time))
    tab.append(("LODF ***", None, time_only_lodf + time_flow_lodf, time_flow_lodf))

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

    if no_pp is False:
        hds = [f"{env_name} ({nb_ts_pp} iter)", f"Δ aor (amps)", f"Δ gen_p (MW)", f"Δ gen_q (MVAr)"]
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
