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

from grid2op import make
from grid2op.Agent import DoNothingAgent
from grid2op.Chronics import ChangeNothing
import re
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
    print("The tabulate package is not installed. Some output might not work properly")

MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"

NICSLU_LICENSE_AVAIL = os.path.exists("./nicslu.lic") and os.path.isfile("./nicslu.lic")


def main(max_ts, env_name_input, test=True,
         no_gs=False, no_gs_synch=False):
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if re.match("^.*\\.json$", env_name_input) is None:
            # i provided an environment name
            env_pp = make(env_name_input, param=param, test=test,
                          data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
            env_lightsim = make(env_name_input, backend=LightSimBackend(), param=param, test=test,
                                data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
        else:
            # I provided an environment path
            env_pp = make("blank", param=param, test=True,
                          data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                          grid_path=env_name_input
                          )
            env_lightsim = make("blank", param=param, test=True,
                                backend=LightSimBackend(),
                                data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                grid_path=env_name_input)
            _, env_name_input = os.path.split(env_name_input)

    agent = DoNothingAgent(action_space=env_pp.action_space)
    nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0, env_seed=0)
    pp_comp_time = env_pp.backend.comp_time
    pp_time_pf = env_pp._time_powerflow

    wst = True  # print extra info in the run_env function
    solver_types = env_lightsim.backend.available_solvers
    if lightsim2grid.SolverType.KLU in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.KLU)
        env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_klu, time_klu, aor_klu, gen_p_klu, gen_q_klu = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        klu_comp_time = env_lightsim.backend.comp_time
        klu_time_pf = env_lightsim._time_powerflow

    if lightsim2grid.SolverType.NICSLU in solver_types and NICSLU_LICENSE_AVAIL:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.NICSLU)
        env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_nicslu, time_nicslu, aor_nicslu, gen_p_nicslu, gen_q_nicslu = run_env(env_lightsim,
                                                                                    max_ts,
                                                                                    agent, chron_id=0,
                                                                                    with_type_solver=wst,
                                                                                    env_seed=0)
        nicslu_comp_time = env_lightsim.backend.comp_time
        nicslu_time_pf = env_lightsim._time_powerflow

    if lightsim2grid.SolverType.SparseLU in solver_types:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.SparseLU)
        env_lightsim.backend.set_solver_max_iter(10)
        nb_ts_slu, time_slu, aor_slu, gen_p_slu, gen_q_slu = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        slu_comp_time = env_lightsim.backend.comp_time
        slu_time_pf = env_lightsim._time_powerflow

    if lightsim2grid.SolverType.GaussSeidel in solver_types and no_gs is False:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidel)
        env_lightsim.backend.set_solver_max_iter(10000)
        nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                with_type_solver=wst, env_seed=0)
        gs_comp_time = env_lightsim.backend.comp_time
        gs_time_pf = env_lightsim._time_powerflow

    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types and no_gs_synch is False:
        env_lightsim.backend.set_solver_type(lightsim2grid.SolverType.GaussSeidelSynch)
        env_lightsim.backend.set_solver_max_iter(10000)
        nb_ts_gsa, time_gsa, aor_gsa, gen_p_gsa, gen_q_gsa = run_env(env_lightsim, max_ts, agent, chron_id=0,
                                                                     with_type_solver=wst, env_seed=0)
        gsa_comp_time = env_lightsim.backend.comp_time
        gsa_time_pf = env_lightsim._time_powerflow

    # NOW PRINT THE RESULTS
    env_name = get_env_name_displayed(env_name_input)
    hds = [f"{env_name}", f"grid2op speed (it/s)", f"grid2op powerflow time (ms)", f"solver powerflow time (ms)"]
    tab = [["PP", f"{nb_ts_pp/time_pp:.2e}",
            f"{1000.*pp_time_pf/nb_ts_pp:.2e}",
            f"{1000.*pp_comp_time/nb_ts_pp:.2e}"]]
    if lightsim2grid.SolverType.GaussSeidel in solver_types and no_gs is False:
        tab.append(["LS+GS", f"{nb_ts_gs/time_gs:.2e}",
                    f"{1000.*gs_time_pf/nb_ts_gs:.2e}",
                    f"{1000.*gs_comp_time/nb_ts_gs:.2e}"])
    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types and no_gs_synch is False:
        tab.append(["LS+GS S", f"{nb_ts_gsa/time_gsa:.2e}",
                    f"{1000.*gsa_time_pf/nb_ts_gsa:.2e}",
                    f"{1000.*gsa_comp_time/nb_ts_gsa:.2e}"])
    if lightsim2grid.SolverType.SparseLU in solver_types:
        tab.append(["LS+SLU", f"{nb_ts_slu/time_slu:.2e}",
                    f"{1000.*slu_time_pf/nb_ts_slu:.2e}",
                    f"{1000.*slu_comp_time/nb_ts_slu:.2e}"])
    if lightsim2grid.SolverType.KLU in solver_types:
        tab.append(["LS+KLU", f"{nb_ts_klu/time_klu:.2e}",
                    f"{1000.*klu_time_pf/nb_ts_klu:.2e}",
                    f"{1000.*klu_comp_time/nb_ts_klu:.2e}"])
    if lightsim2grid.SolverType.NICSLU in solver_types:
        tab.append(["LS+NICSLU", f"{nb_ts_nicslu/time_nicslu:.2e}",
                    f"{1000.*nicslu_time_pf/nb_ts_nicslu:.2e}",
                    f"{1000.*nicslu_comp_time/nb_ts_nicslu:.2e}"])

    if TABULATE_AVAIL:
        res_use_with_grid2op_1 = tabulate(tab, headers=hds,  tablefmt="rst")
        print(res_use_with_grid2op_1)
    else:
        print(tab)
    print()

    if TABULATE_AVAIL:
        res_github_readme = tabulate(tab, headers=hds,  tablefmt="github")
        print(res_github_readme)
    else:
        print(tab)
    print()

    hds = [f"{env_name} ({nb_ts_pp} iter)", f"Δ aor (amps)", f"Δ gen_p (MW)", f"Δ gen_q (MVAr)"]
    tab = [["PP", "0.00", "0.00", "0.00"]]
    if lightsim2grid.SolverType.GaussSeidel in solver_types and no_gs is False:
        tab.append(["LS+GS",
                    f"{np.max(np.abs(aor_gs - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_gs - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_gs - gen_q_pp)):.2e}"])
    if lightsim2grid.SolverType.GaussSeidelSynch in solver_types and no_gs_synch is False:
        tab.append(["LS+GS S",
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
    if lightsim2grid.SolverType.NICSLU in solver_types:
        tab.append(["LS+NICSLU",
                    f"{np.max(np.abs(aor_nicslu - aor_pp)):.2e}",
                    f"{np.max(np.abs(gen_p_nicslu - gen_p_pp)):.2e}",
                    f"{np.max(np.abs(gen_q_nicslu - gen_q_pp)):.2e}"])

    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab, headers=hds,  tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab)
    print()


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
                        help='Do not use \"test=True\" keyword argument when building the grid2op environments'
                             ' for the benchmark (default False: use \"test=True\"  environment)')
    parser.add_argument('--no_gs_synch', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not benchmark gauss seidel (synch) method (default: evaluate it)')
    parser.add_argument('--no_gs', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not benchmark gauss seidel (regular) method (default: evaluate it)')

    args = parser.parse_args()

    max_ts = int(args.number)
    name = str(args.name)
    test_env = not args.no_test
    main(max_ts, name, test_env, no_gs =args.no_gs, no_gs_synch=args.no_gs_synch)
