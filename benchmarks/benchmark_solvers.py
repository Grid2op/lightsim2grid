# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import os

from grid2op import make
from grid2op.Agent import DoNothingAgent
try:
    from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance as GridStateFromFile
except ImportError:
    print("Be carefull: there might be maintenance")
    from grid2op.Chronics import GridStateFromFile

from grid2op.Parameters import Parameters
import lightsim2grid
from lightsim2grid.LightSimBackend import LightSimBackend
from utils_benchmark import print_res, run_env, str2bool, get_env_name_displayed
TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabluate package is not installed. Some output might not work properly")

import pdb

MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"


def main(max_ts, ENV_NAME, test=True):
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

    env_pp = make(ENV_NAME, param=param, test=test,
                   data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
    agent = DoNothingAgent(action_space=env_pp.action_space)
    nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0, env_seed=0)
    pp_comp_time = env_pp.backend.comp_time
    pp_time_pf = env_pp._time_powerflow

    wst = True  # print extra info in the run_env function
    env_lightsim = make(ENV_NAME, backend=LightSimBackend(), param=param, test=test,
                        data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
    solver_types = env_lightsim.backend.available_solvers
    if lightsim2grid.SolverType.KLU in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.KLU)
        env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_klu, time_klu, aor_klu, gen_p_klu, gen_q_klu = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        klu_comp_time = env_lightsim.backend.comp_time
        klu_time_pf = env_lightsim._time_powerflow
    if lightsim2grid.SolverType.SparseLU in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.SparseLU)
        env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_slu, time_slu, aor_slu, gen_p_slu, gen_q_slu = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        slu_comp_time = env_lightsim.backend.comp_time
        slu_time_pf = env_lightsim._time_powerflow
    if lightsim2grid.SolverType.GaussSeidel in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidel)
        env_lightsim.backend.set_solver_max_iter(10000)
        nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                with_type_solver=wst, env_seed=0)
        gs_comp_time = env_lightsim.backend.comp_time
        gs_time_pf = env_lightsim._time_powerflow

    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidelSynch)
        env_lightsim.backend.set_solver_max_iter(10000)
        nb_ts_gsa, time_gsa, aor_gsa, gen_p_gsa, gen_q_gsa = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        gsa_comp_time = env_lightsim.backend.comp_time
        gsa_time_pf = env_lightsim._time_powerflow

    # NOW PRINT THE RESULTS
    env_name = get_env_name_displayed(ENV_NAME)
    hds = [f"{env_name}", f"grid2op speed (it/s)", f"grid2op powerflow time (ms)", f"solver powerflow time (ms)"]
    tab = [["PP", int(nb_ts_pp/time_pp),
            f"{1000.*pp_time_pf/nb_ts_pp:.2e}",
            f"{1000.*pp_comp_time/nb_ts_pp:.2e}"]]
    if lightsim2grid.SolverType.GaussSeidel in solver_types:
        tab.append(["LS+GS", int(nb_ts_gs/time_gs),
                    f"{1000.*gs_time_pf/nb_ts_gs:.2e}",
                    f"{1000.*gs_comp_time/nb_ts_gs:.2e}"])
    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types:
        tab.append(["LS+GS A", int(nb_ts_gsa/time_gsa),
                    f"{1000.*gsa_time_pf/nb_ts_gsa:.2e}",
                    f"{1000.*gsa_comp_time/nb_ts_gsa:.2e}"])
    if lightsim2grid.SolverType.SparseLU in solver_types:
        tab.append(["LS+SLU", int(nb_ts_slu/time_slu),
                    f"{1000.*slu_time_pf/nb_ts_slu:.2e}",
                    f"{1000.*slu_comp_time/nb_ts_slu:.2e}"])
    if lightsim2grid.SolverType.KLU in solver_types:
        tab.append(["LS+KLU", int(nb_ts_klu/time_klu),
                    f"{1000.*klu_time_pf/nb_ts_klu:.2e}",
                    f"{1000.*klu_comp_time/nb_ts_klu:.2e}"])
    res_use_with_grid2op_1 = tabulate(tab, headers=hds,  tablefmt="rst")
    print(res_use_with_grid2op_1)
    print()

    res_github_readme = tabulate(tab, headers=hds,  tablefmt="github")
    print(res_github_readme)
    print()

    hds = [f"{env_name} ({nb_ts_pp} iter)", f"Δ aor (amps)", f"Δ gen_p (MW)", f"Δ gen_q (MVAr)"]
    tab = [["PP", "0.00", "0.00", "0.00"]]
    if lightsim2grid.SolverType.GaussSeidel in solver_types:
        tab.append(["LS+GS",
                    f"{np.max(np.abs(aor_gs - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_gs - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_gs - gen_q_pp)):.2e}"])
    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types:
        tab.append(["LS+GS A",
                    f"{np.max(np.abs(aor_gsa - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_gsa - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_gsa - gen_q_pp)):.2e}"])
    if lightsim2grid.SolverType.SparseLU in solver_types:
        tab.append(["LS+SLU",
                    f"{np.max(np.abs(aor_slu - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_slu - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_slu - gen_q_pp)):.2e}"])
    if lightsim2grid.SolverType.KLU in solver_types:
        tab.append(["LS+KLU",
                    f"{np.max(np.abs(aor_klu - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_klu - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_klu - gen_q_pp)):.2e}"])

    res_use_with_grid2op_2 = tabulate(tab, headers=hds,  tablefmt="rst")
    print(res_use_with_grid2op_2)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark of lightsim with a "do nothing" agent '
                                                 '(compare multiple lightsim solvers with default grid2op backend '
                                                 'PandaPower)')
    parser.add_argument('--name', default=ENV_NAME, type=str,
                        help='Environment name to be used for the benchmark.')
    parser.add_argument('--number', type=int, default=MAX_TS,
                        help='Maximum number of time steps for which the benchmark will be run.')
    parser.add_argument('--no_test', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not use test environment for the benchmark (default False: use test environment)')

    args = parser.parse_args()

    max_ts = int(args.number)
    name = str(args.name)
    test_env = not args.no_test
    main(max_ts, name, test_env)
