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
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid import LightSimBackend, TimeSerie, ContingencyAnalysis
from utils_benchmark import print_res, run_env, str2bool, get_env_name_displayed, print_configuration
TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Some output might not work properly")
    
try:
    from pypowsybl2grid import PyPowSyBlBackend
    PyPowSyBlBackend.shunts_data_available = False
    PYPOW_ERROR = None
except ImportError as exc_:
    PYPOW_ERROR = exc_
    print("Backend based on pypowsybl will not be benchmarked")
    
MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"
DONT_SAVE = "__DONT_SAVE"
NICSLU_LICENSE_AVAIL = os.path.exists("./nicslu.lic") and os.path.isfile("./nicslu.lic")

solver_names = {AlgorithmType.DC_SparseLU: "DC (SparseLU)",
                AlgorithmType.DC_KLU: "DC (KLU)",
                AlgorithmType.DC_NICSLU: "DC (NICSLU\\*)",
                AlgorithmType.DC_CKTSO: "DC (CKTSO\\*)"
                }
solver_gs = {}
solver_fdpf = {}
res_times = {}

order_solver_print = [
    AlgorithmType.DC_SparseLU,
    AlgorithmType.DC_KLU,
    AlgorithmType.DC_NICSLU,
    AlgorithmType.DC_CKTSO,
    
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
                        pypow_error=True, nb_ts_pypow=None, time_pypow=None,
                        pypow_time_pf=None, pypow_comp_time=None,
                        ts_time=None, ts_algo_time=None,
                        ptdf_time=None,
                        sa_time=None, sa_algo_time=None,
                        lodf_time=None):
    """Generates (from the numbers actually measured in this run) the descriptive text that should go
    under the tables of docs/benchmarks_dc.rst, the same way `generate_narrative` in
    benchmark_solvers.py does for docs/benchmarks.rst.

    This is what used to be summarized (and only partially, with numbers that do not even match the
    tables) by hand in the "TL;DR" section of that page.
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
            f"From a grid2op perspective, lightsim2grid allows to compute up to ~{best['speed']:.0f} DC steps "
            f"each second (column `grid2op speed`, row `{best['name']}`) on the {env_name} and \"only\" "
            f"~{pp_speed:.0f} for the default PandaPower Backend (column `grid2op speed`, row `PP DC`), "
            f"leading to a speed up of **~{speedup:.0f}** ({best['speed']:.0f} / {pp_speed:.0f}) in this case."
        )

    if pypow_error is None and nb_ts_pypow:
        pypow_speed, pypow_runpf_ms, pypow_algo_ms = _speed_metrics(nb_ts_pypow, time_pypow, pypow_comp_time,
                                                                     pypow_time_pf)
        speedup_pypow = best["speed"] / pypow_speed
        speed_paragraphs.append(
            f"When compared to powsybl (with the pypowsybl backend), lightsim2grid is around "
            f"**~{speedup_pypow:.1f}** times faster ({pypow_speed:.0f} vs {best['speed']:.0f})."
        )

    slu = _get_speed(AlgorithmType.DC_SparseLU)
    klu = _get_speed(AlgorithmType.DC_KLU)
    if slu is not None and klu is not None:
        rel = abs(klu - slu) / max(klu, slu)
        qualifier = "no sensible" if rel < 0.1 else "a sensible"
        speed_paragraphs.append(
            f"For this environment there is {qualifier} difference in using `KLU` linear solver "
            f"(row `DC (KLU)`) compared to using the SparseLU solver of Eigen (row `DC`) "
            f"({slu:.0f} vs {klu:.0f} iterations on the reported runs, might slightly vary across runs)."
        )

    nicslu = _get_speed(AlgorithmType.DC_NICSLU)
    cktso = _get_speed(AlgorithmType.DC_CKTSO)
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

    if ts_time:
        speedup_ts = best["runpf_ms"] / ts_time
        speed_paragraphs.append(
            f"The `TimeSerie` module performs one DC powerflow in {ts_time:.3g} ms on average (row `time serie`, "
            f"column `grid2op 'backend.runpf' time`), compared to {best['runpf_ms']:.3g} ms for the fastest "
            f"grid2op DC backend (`{best['name']}`), a **~{speedup_ts:.0f}x** speed up."
        )
    if sa_time:
        speedup_sa = best["runpf_ms"] / sa_time
        speed_paragraphs.append(
            f"Similarly, the `ContingencyAnalysis` module performs one DC contingency in {sa_time:.3g} ms on "
            f"average (row `contingency analysis`), a **~{speedup_sa:.0f}x** speed up compared to the fastest "
            f"grid2op DC backend."
        )
    if ptdf_time:
        speedup_ptdf = best["runpf_ms"] / ptdf_time
        speed_paragraphs.append(
            f"Using the PTDF matrix directly (row `PTDF`) is even faster: {ptdf_time:.3g} ms per powerflow, a "
            f"**~{speedup_ptdf:.0f}x** speed up compared to the fastest grid2op DC backend."
        )
    if lodf_time:
        speedup_lodf = best["runpf_ms"] / lodf_time
        speed_paragraphs.append(
            f"Likewise, using the LODF matrix (row `LODF`) to perform the contingency analysis takes "
            f"{lodf_time:.3g} ms per contingency, a **~{speedup_lodf:.0f}x** speed up compared to the fastest "
            f"grid2op DC backend."
        )

    diff_paragraphs = []
    if aor_pp is not None and gen_p_pp is not None and gen_q_pp is not None:
        all_daor = {t: np.max(np.abs(m["aor"] - aor_pp)) for t, m in ls_metrics.items()}
        all_dgenp = {t: np.max(np.abs(m["gen_p"] - gen_p_pp)) for t, m in ls_metrics.items()}
        all_dgenq = {t: np.max(np.abs(m["gen_q"] - gen_q_pp)) for t, m in ls_metrics.items()}
        max_daor = max(all_daor.values())
        max_dgenp = max(all_dgenp.values())
        max_dgenq = max(all_dgenq.values())
        nb_ts_ref = best["nb_ts"]
        diff_paragraphs.append(
            f"Here are the results for the {env_name} (max over {nb_ts_ref} DC powerflows). The maximum "
            f"difference against pandapower (PP DC, the reference), across all tested DC solvers, is "
            f"{max_daor:.2e} A for the line flows, {max_dgenp:.2e} MW for gen_p and {max_dgenq:.2e} MVAr "
            f"for gen_q."
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
                    bk = PyPowSyBlBackend()
                    env_pypow = make(env_name_input, param=param, test=test,
                                    backend=bk,
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
                    bk = PyPowSyBlBackend()
                    env_pypow = make("blank", param=param, test=True,
                                    data_feeding_kwargs={"gridvalueClass": ChangeNothing},
                                    grid_path=env_name_input,
                                    backend=bk)
                except Grid2OpException as exc_:
                    pypow_error = exc_
            _, env_name_input = os.path.split(env_name_input)
    agent = DoNothingAgent(action_space=env_pp.action_space)
    if no_pp is False:
        print("Start using Pandapower")
        nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0, env_seed=0,
                                                                is_dc=True)
        pp_comp_time = env_pp.backend.comp_time
        pp_time_pf = env_pp._time_powerflow
    
    if pypow_error is None:
        # also benchmark pypowsybl backend
        nb_ts_pypow, time_pypow, aor_pypow, gen_p_pypow, gen_q_pypow = run_env(env_pypow, max_ts, agent, chron_id=0, env_seed=0,
                                                                               is_dc=True)
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
                                                                with_type_solver=wst, env_seed=0, is_dc=True)
        gs_comp_time = env_lightsim.backend.comp_time
        gs_time_pf = env_lightsim._time_powerflow
        res_times[solver_type] = (solver_names[solver_type],
                                  nb_ts_gs, time_gs, aor_gs, gen_p_gs,
                                  gen_q_gs, gs_comp_time, gs_time_pf)
    # sys.exit(0)
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
    computer_ts.change_algorithm(AlgorithmType.DC_KLU)
    v_init = real_env_ls.backend.V
    status = computer_ts.compute_Vs(prod_p,
                                    np.zeros((nb_ts_gs, 0), dtype=dt_float),
                                    load_p,
                                    load_q,
                                    v_init,
                                    real_env_ls.backend.max_it,
                                    real_env_ls.backend.tol)
    time_serie._TimeSerie__computed = True
    a_or = time_serie.compute_A()  # noqa: F841
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
    flows = np.dot(PTDF_, T_Sbus).T  # noqa: F841
    end_ = time.perf_counter()
    time_only_ptdf, *_ = real_env_ls.backend._grid.get_dc_solver().get_timers_ptdf_lodf()
    time_only_ptdf *= 1000. / Sbus.shape[0]
    time_flow_ptdf = 1000. * (end_ - beg_) / Sbus.shape[0]
    # ts_time = 1e3 * (computer_ts.total_time() + computer_ts.amps_computation_time()) / computer_ts.nb_solved()
    # ts_algo_time = 1e3 * (computer_ts.solver_time()) / computer_ts.nb_solved()  
    
    # Perform a securtiy analysis (up to 1000 contingencies)
    real_env_ls.reset()
    sa = ContingencyAnalysis(real_env_ls)
    print("BEGINNING SA")
    computer_sa = sa.computer
    computer_sa.change_algorithm(AlgorithmType.DC_KLU)
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
    por_lodf = init_powerflow + LODF_mat * init_powerflow_diag  # noqa: F841
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

    hds = [f"{env_name}", "grid2op speed (it/s)", "grid2op 'backend.runpf' time (ms / pf)", "time in 'algo' (ms / pf)"]
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
    tab.append(("time serie \\*\\*", None, ts_time, ts_algo_time))
    tab.append(("PTDF \\*\\*", None, time_only_ptdf + time_flow_ptdf, time_flow_ptdf))
    tab.append(("contingency analysis \\*\\*\\*", None, sa_time, sa_algo_time))
    tab.append(("LODF \\*\\*\\*", None, time_only_lodf + time_flow_lodf, time_flow_lodf))

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
        hds = [f"{env_name} ({nb_ts_pp} iter)", "Δ aor (amps)", "Δ gen_p (MW)", "Δ gen_q (MVAr)"]
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

    # generate (from the numbers measured above) the descriptive text that should go under the
    # tables in docs/benchmarks_dc.rst, so that it no longer needs to be re-derived and copy / pasted
    # by hand (which is how the current "TL;DR" section on that page drifted out of sync with its
    # own tables).
    speed_text, diff_text = generate_narrative(
        env_name, res_times,
        no_pp=no_pp,
        nb_ts_pp=nb_ts_pp if no_pp is False else None,
        time_pp=time_pp if no_pp is False else None,
        pp_time_pf=pp_time_pf if no_pp is False else None,
        pp_comp_time=pp_comp_time if no_pp is False else None,
        aor_pp=aor_pp if no_pp is False else None,
        gen_p_pp=gen_p_pp if no_pp is False else None,
        gen_q_pp=gen_q_pp if no_pp is False else None,
        pypow_error=pypow_error,
        nb_ts_pypow=nb_ts_pypow if pypow_error is None else None,
        time_pypow=time_pypow if pypow_error is None else None,
        pypow_time_pf=pypow_time_pf if pypow_error is None else None,
        pypow_comp_time=pypow_comp_time if pypow_error is None else None,
        ts_time=ts_time, ts_algo_time=ts_algo_time,
        ptdf_time=time_only_ptdf + time_flow_ptdf,
        sa_time=sa_time, sa_algo_time=sa_algo_time,
        lodf_time=time_only_lodf + time_flow_lodf,
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
