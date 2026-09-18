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
from lightsim2grid import AlgorithmType
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

solver_names = {AlgorithmType.GaussSeidel: "GS",
                AlgorithmType.GaussSeidelSynch: "GS synch",
                AlgorithmType.NR_SparseLU: "NR (SLU)",
                AlgorithmType.NR_KLU: "NR (KLU)",
                AlgorithmType.NR_NICSLU: "NR (NICSLU\\*)",
                AlgorithmType.NR_CKTSO: "NR (CKTSO\\*)",
                AlgorithmType.NRSing_SparseLU: "NR single (SLU)",
                AlgorithmType.NRSing_KLU: "NR single (KLU)",
                AlgorithmType.NRSing_NICSLU: "NR single (NICSLU\\*)",
                AlgorithmType.NRSing_CKTSO: "NR single (CKTSO\\*)",
                AlgorithmType.FDPF_XB_SparseLU: "FDPF XB (SLU)",
                AlgorithmType.FDPF_BX_SparseLU: "FDPF BX (SLU)",
                AlgorithmType.FDPF_XB_KLU: "FDPF XB (KLU)",
                AlgorithmType.FDPF_BX_KLU: "FDPF BX (KLU)",
                AlgorithmType.FDPF_XB_NICSLU: "FDPF XB (NICSLU\\*)",
                AlgorithmType.FDPF_BX_NICSLU: "FDPF BX (NICSLU\\*)",
                AlgorithmType.FDPF_XB_CKTSO: "FDPF XB (CKTSO\\*)",
                AlgorithmType.FDPF_BX_CKTSO: "FDPF BX (CKTSO\\*)",
                }
solver_gs = {AlgorithmType.GaussSeidelSynch, AlgorithmType.GaussSeidel}
solver_fdpf = {AlgorithmType.FDPF_XB_SparseLU, AlgorithmType.FDPF_BX_SparseLU,
               AlgorithmType.FDPF_XB_KLU, AlgorithmType.FDPF_BX_KLU,
               AlgorithmType.FDPF_XB_NICSLU, AlgorithmType.FDPF_BX_NICSLU,
               AlgorithmType.FDPF_XB_CKTSO, AlgorithmType.FDPF_BX_CKTSO,
               }
solver_nr = set(solver_names.keys()) - solver_gs - solver_fdpf
res_times = {}
order_solver_print = [
    AlgorithmType.GaussSeidel,
    AlgorithmType.GaussSeidelSynch,
    AlgorithmType.NRSing_SparseLU,
    AlgorithmType.NR_SparseLU,
    AlgorithmType.NRSing_KLU,
    AlgorithmType.NR_KLU,
    AlgorithmType.NRSing_NICSLU,
    AlgorithmType.NR_NICSLU,
    AlgorithmType.NRSing_CKTSO,
    AlgorithmType.NR_CKTSO,
    AlgorithmType.FDPF_XB_SparseLU,
    AlgorithmType.FDPF_BX_SparseLU,
    AlgorithmType.FDPF_XB_KLU,
    AlgorithmType.FDPF_BX_KLU,
    AlgorithmType.FDPF_XB_NICSLU,
    AlgorithmType.FDPF_BX_NICSLU,
    AlgorithmType.FDPF_XB_CKTSO,
    AlgorithmType.FDPF_BX_CKTSO,
]


def _speed_metrics(nb_ts, elapsed_time, comp_time, time_pf):
    """returns (grid2op speed [it/s], grid2op 'backend.runpf' time [ms], time in 'algo' [ms]) for one backend"""
    speed = nb_ts / elapsed_time
    runpf_ms = 1000. * time_pf / nb_ts
    algo_ms = 1000. * comp_time / nb_ts
    return speed, runpf_ms, algo_ms


def generate_narrative(env_name,
                        res_times,
                        no_pp, nb_ts_pp=None, time_pp=None, pp_time_pf=None, pp_comp_time=None,
                        aor_pp=None, gen_p_pp=None, gen_q_pp=None,
                        pypowbk_error=True, nb_ts_pypow=None, time_pypow=None,
                        pypow_time_pf=None, pypow_comp_time=None):
    """Generates (from the numbers actually measured in this run) the descriptive text that goes
    under the "computation time" and "differences" tables in docs/benchmarks.rst.

    This avoids the manual, error prone, floating point formatting and copy / pasting that was
    previously needed to keep this text in sync with the tables above it.
    """
    ls_metrics = {}
    for algo_type, data in res_times.items():
        solver_name, nb_ts_gs, time_gs, aor_gs, gen_p_gs, gen_q_gs, gs_comp_time, gs_time_pf = data
        speed, runpf_ms, algo_ms = _speed_metrics(nb_ts_gs, time_gs, gs_comp_time, gs_time_pf)
        ls_metrics[algo_type] = {"name": solver_name, "speed": speed, "runpf_ms": runpf_ms, "algo_ms": algo_ms,
                                  "aor": aor_gs, "gen_p": gen_p_gs, "gen_q": gen_q_gs, "nb_ts": nb_ts_gs}
    if not ls_metrics:
        return ""

    def _get_speed(algo_type):
        return ls_metrics[algo_type]["speed"] if algo_type in ls_metrics else None

    best_type = max(ls_metrics, key=lambda k: ls_metrics[k]["speed"])
    best = ls_metrics[best_type]

    speed_paragraphs = []
    have_pp = (not no_pp) and nb_ts_pp

    if have_pp:
        pp_speed, pp_runpf_ms, pp_algo_ms = _speed_metrics(nb_ts_pp, time_pp, pp_comp_time, pp_time_pf)
        speedup = best["speed"] / pp_speed
        speed_paragraphs.append(
            f"From a grid2op perspective, lightsim2grid allows to compute up to ~{best['speed']:.0f} steps each "
            f"second (column `grid2op speed`, row `{best['name']}`) on the {env_name} and \"only\" ~{pp_speed:.0f} "
            f"for the default PandaPower Backend (column `grid2op speed`, row `PP`), leading to a speed up of "
            f"**~{speedup:.0f}** ({best['speed']:.0f} / {pp_speed:.0f}) in this case "
            f"(lightsim2grid Backend is ~{speedup:.0f} times faster than pandapower Backend when comparing "
            f"grid2op speed)."
        )

    if pypowbk_error is None and nb_ts_pypow:
        pypow_speed, pypow_runpf_ms, pypow_algo_ms = _speed_metrics(nb_ts_pypow, time_pypow, pypow_comp_time,
                                                                     pypow_time_pf)
        speedup_pypow = best["speed"] / pypow_speed
        speed_paragraphs.append(
            f"When compared to powsybl (with the pypowsybl backend), lightsim2grid (with newton raphson) is "
            f"around **~{speedup_pypow:.1f}** times faster ({pypow_speed:.0f} vs {best['speed']:.0f})."
        )

    slu = _get_speed(AlgorithmType.NR_SparseLU)
    klu = _get_speed(AlgorithmType.NR_KLU)
    if slu is not None and klu is not None:
        rel = abs(klu - slu) / max(klu, slu)
        qualifier = "no sensible" if rel < 0.1 else "a sensible"
        speed_paragraphs.append(
            f"For this environment there is {qualifier} difference in using `KLU` linear solver "
            f"(rows `NR single (KLU)` or `NR (KLU)`) compared to using the SparseLU solver of Eigen "
            f"(rows `NR single (SLU)` or `NR (SLU)`) ({slu:.0f} vs {klu:.0f} iterations on the reported runs, "
            f"might slightly vary across runs)."
        )

    nicslu = _get_speed(AlgorithmType.NR_NICSLU)
    cktso = _get_speed(AlgorithmType.NR_CKTSO)
    present = [v for v in (klu, nicslu, cktso) if v is not None]
    if len(present) >= 2:
        spread = (max(present) - min(present)) / max(present)
        if spread < 0.1:
            speed_paragraphs.append(
                "Linear solvers `KLU`, `NICSLU` and `CKTSO` achieve almost identical performances, "
                "at least we think the observed differences are within error margins."
            )
        else:
            speed_paragraphs.append(
                f"Linear solvers `KLU`, `NICSLU` and `CKTSO` show up to {100. * spread:.0f}% difference in "
                f"performance in this run."
            )

    single_klu = _get_speed(AlgorithmType.NRSing_KLU)
    if klu is not None and single_klu is not None:
        rel = abs(klu - single_klu) / max(klu, single_klu)
        qualifier = "very little" if rel < 0.1 else "some noticeable"
        speed_paragraphs.append(
            f"There are also {qualifier} differences between non distributed slack (`NR Single` rows) "
            f"and distributed slack (`NR` rows) for the tested linear solvers."
        )

    fdpf_candidates = [v for v in (_get_speed(AlgorithmType.FDPF_XB_KLU), _get_speed(AlgorithmType.FDPF_BX_KLU))
                       if v is not None]
    if fdpf_candidates and klu is not None:
        fdpf_avg = sum(fdpf_candidates) / len(fdpf_candidates)
        rel = (klu - fdpf_avg) / klu
        if abs(rel) < 0.1:
            qualifier = "equivalent"
        elif rel > 0:
            qualifier = "slightly less performant"
        else:
            qualifier = "slightly more performant"
        speed_paragraphs.append(
            f"Finally, the \"fast decoupled\" methods also lead to {qualifier} performances (compared to the "
            f"Newton Raphson one) for almost all linear solvers."
        )

    step_ms = 1000. / best["speed"]
    algo_ms = best["algo_ms"]
    step_us = step_ms * 1000.
    algo_us = algo_ms * 1000.
    ext_us = step_us - algo_us
    ext_pct = 100. * ext_us / step_us if step_us else 0.
    speed_paragraphs.append(
        "For this environment, for lightsim2grid backend (and if we don't take into account the \"agent time\"), "
        "the computation time is vastly dominated by factor external to the powerflow solver. Indeed, doing a "
        f"'env.step' (column `grid2op speed (it/s)`) takes {step_ms:.3g}ms (`1. / {best['speed']:.0f}. * 1000.`) "
        f"on average and on this {step_us:.0f} µs (or {step_ms:.3g}ms), only {algo_us:.0f} µs are spent "
        "in the backend (column `time in 'algo' (ms / pf)`). Meaning that "
        f"~{ext_us:.0f} µs are spent in the grid2op extra layer or in the backend implementation in this "
        f"case (`{ext_pct:.0f}%` of the computation time - `={ext_us:.0f} / {step_us:.0f}`- is external to the "
        "powerflow algorithm)"
    )

    diff_paragraphs = []
    if aor_pp is not None and gen_p_pp is not None and gen_q_pp is not None:
        all_daor = {t: np.max(np.abs(m["aor"] - aor_pp)) for t, m in ls_metrics.items()}
        all_dgenp = {t: np.max(np.abs(m["gen_p"] - gen_p_pp)) for t, m in ls_metrics.items()}
        all_dgenq = {t: np.max(np.abs(m["gen_q"] - gen_q_pp)) for t, m in ls_metrics.items()}
        max_daor = max(all_daor.values())
        nb_ts_ref = best["nb_ts"]
        diff_paragraphs.append(
            f"Here are the results for the {env_name} (max over {nb_ts_ref} powerflows). Flows are here measured "
            f"in amps (and not kA). The maximum difference of flows (across all tested solvers) is approximately "
            f"{max_daor * 1000.:.3g}mA or {max_daor:.1e} A. This difference is totally neglectible on power "
            f"transportation side where the current is usually around 1kA (1e3 A)."
        )

        nr_daor = [v for t, v in all_daor.items() if t in solver_nr]
        nr_dgenp = [v for t, v in all_dgenp.items() if t in solver_nr]
        nr_dgenq = [v for t, v in all_dgenq.items() if t in solver_nr]
        if nr_daor:
            max_nr = max(max(nr_dgenp), max(nr_dgenq))
            diff_paragraphs.append(
                "As we can see, the difference when using lightsim and pandapower is rather small, even when "
                "using a different algorithm to solve the powerflow. When using Newton Raphson solvers, the "
                "difference in absolute values when using lightsim2grid compared with using PandaPowerBackend "
                f"is {max_nr:.1e} for the generators active/reactive productions and {max(nr_daor):.1e} A for "
                "the powerline flows (`PP (ref)` row has, by construction, a difference of 0. in all columns)."
            )

    return "\n\n".join(speed_paragraphs), "\n\n".join(diff_paragraphs)


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
    solver_types = env_lightsim.backend.available_default_algorithms
    for solver_type in solver_types:
        if solver_type not in solver_names:
            continue
        print(f"Start using {solver_type}")
        env_lightsim.backend.set_algo_type(solver_type)
        if solver_type in solver_gs:
            # gauss seidel sovler => more iterations
            env_lightsim.backend.set_solver_max_iter(10000)
            if AlgorithmType.GaussSeidel == solver_type and no_gs:
                # I don't study the gauss seidel solver
                continue
            elif AlgorithmType.GaussSeidelSynch  == solver_type and no_gs_synch:
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

    # generate (from the numbers measured above) the descriptive text that goes under the tables
    # in docs/benchmarks.rst, so that it no longer needs to be re-derived and copy / pasted by hand.
    speed_text, diff_text = generate_narrative(
        env_name, res_times,
        no_pp=no_pp,
        nb_ts_pp=nb_ts_pp if no_pp is False else None,
        time_pp=time_pp if no_pp is False else None,
        pp_time_pf=pp_time_pf if no_pp is False else None,
        pp_comp_time=pp_comp_time if no_pp is False else None,
        aor_pp=aor_pp, gen_p_pp=gen_p_pp, gen_q_pp=gen_q_pp,
        pypowbk_error=pypowbk_error,
        nb_ts_pypow=nb_ts_pypow if pypowbk_error is None else None,
        time_pypow=time_pypow if pypowbk_error is None else None,
        pypow_time_pf=pypow_time_pf if pypowbk_error is None else None,
        pypow_comp_time=pypow_comp_time if pypowbk_error is None else None,
    )
    print("Description (computation time):")
    print(speed_text)
    print()
    print("Description (differences):")
    print(diff_text)
    print()
    if save_results != DONT_SAVE:
        with open(save_results+"description.rst", "w", encoding="utf-8") as f:
            f.write("Description (computation time):\n\n")
            f.write(speed_text)
            f.write("\n\nDescription (differences):\n\n")
            f.write(diff_text)
            f.write("\n")


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
