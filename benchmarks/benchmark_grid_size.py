# Copyright (c) 2020, RTE (https://www.rte-france.com)
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
            #   "GBnetwork.json",  # 2224 buses
              "case2848rte.json",
              "case2869pegase.json",
              "case3120sp.json",
              "case6495rte.json",
              "case6515rte.json",
              "case9241pegase.json"
              ]

def make_grid2op_env(pp_case, case_name, load_p, load_q, gen_p, sgen_p):
    param = Parameters.Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env_lightsim = make("blank",
                            param=param, test=True,
                            backend=LightSimBackend(),
                            chronics_class=FromNPY,
                            data_feeding_kwargs={"load_p": load_p,
                                                 "load_q": load_q,
                                                 "prod_p": gen_p
                                                },
                            grid_path=case_name,
                            _add_to_name=f"{case_name}",
                            )
    return env_lightsim


def make_grid2op_env_pp(pp_case, case_name, load_p, load_q, gen_p, sgen_p):
    param = Parameters.Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env_pp = make("blank",
                      param=param, test=True,
                      backend=PandaPowerBackend(lightsim2grid=False),
                      chronics_class=FromNPY,
                      data_feeding_kwargs={"load_p": load_p,
                                           "load_q": load_q,
                                           "prod_p": gen_p
                                          },
                      grid_path=case_name,
                      _add_to_name=f"{case_name}",
                     )
    return env_pp


def get_loads_gens(load_p_init, load_q_init, gen_p_init, sgen_p_init, prng):
    # scale loads

    # use some French time series data for loads
    # see https://github.com/BDonnot/data_generation for where to find this file
    coeffs = {"sources": {
    "country": "France",
    "year": "2012",
    "web": "http://clients.rte-france.com/lang/fr/visiteurs/vie/vie_stats_conso_inst.jsp"
    },
    "month": {
    "jan": 1.21,
    "feb": 1.40,
    "mar": 1.05,
    "apr": 1.01,
    "may": 0.86,
    "jun": 0.84,
    "jul": 0.84,
    "aug": 0.79,
    "sep": 0.85,
    "oct": 0.94,
    "nov": 1.01,
    "dec": 1.20
    },
    "day": {
    "mon": 1.01,
    "tue": 1.05,
    "wed": 1.05,
    "thu": 1.05,
    "fri": 1.03,
    "sat": 0.93,
    "sun": 0.88
    },
    "hour": {
    "00:00": 1.00,
    "01:00": 0.93,
    "02:00": 0.91,
    "03:00": 0.86,
    "04:00": 0.84,
    "05:00": 0.85,
    "06:00": 0.90,
    "07:00": 0.97,
    "08:00": 1.03,
    "09:00": 1.06,
    "10:00": 1.08,
    "11:00": 1.09,
    "12:00": 1.09,
    "13:00": 1.09,
    "14:00": 1.06,
    "15:00": 1.03,
    "16:00": 1.00,
    "17:00": 1.00,
    "18:00": 1.04,
    "19:00": 1.09,
    "20:00": 1.05,
    "21:00": 1.01,
    "22:00": 0.99,
    "23:00": 1.03
    }
    }
    vals = list(coeffs["hour"].values())
    x_final = np.arange(12 * len(vals))

    # interpolate them at 5 minutes resolution (instead of 1h)
    vals.append(vals[0])
    vals = np.array(vals) * coeffs["month"]["oct"] * coeffs["day"]["mon"]
    x_interp = 12 * np.arange(len(vals))
    coeffs = interp1d(x=x_interp, y=vals, kind="cubic")
    all_vals = coeffs(x_final).reshape(-1, 1)
    if DEBUG:
        all_vals[:] = 1
    
    # compute the "smooth" loads matrix
    load_p_smooth = all_vals * load_p_init.reshape(1, -1)
    load_q_smooth = all_vals * load_q_init.reshape(1, -1)

    # add a bit of noise to it to get the "final" loads matrix
    load_p = load_p_smooth * prng.lognormal(mean=0., sigma=0.003, size=load_p_smooth.shape)
    load_q = load_q_smooth * prng.lognormal(mean=0., sigma=0.003, size=load_q_smooth.shape)
    if DEBUG:
        load_p[:] = load_p_smooth
        load_q[:] = load_q_smooth
    
    # scale generators accordingly
    gen_p = load_p.sum(axis=1).reshape(-1, 1) / load_p_init.sum() * gen_p_init.reshape(1, -1)
    sgen_p = load_p.sum(axis=1).reshape(-1, 1) / load_p_init.sum() * sgen_p_init.reshape(1, -1)
    return load_p, load_q, gen_p, sgen_p


def run_grid2op_env(env_lightsim, case, reset_solver,
                    solver_preproc_solver_time, 
                    g2op_speeds,
                    g2op_step_time,
                    ls_solver_time,
                    ls_gridmodel_time,
                    g2op_sizes,
                    sgen_p,
                    nb_ts
                    ):
    _ = env_lightsim.reset()
    done = False
    nb_step = 0
    changed_sgen = case.sgen["in_service"].values
    while not done:
        # hack for static gen...
        changed_sgen = copy.deepcopy(case.sgen["in_service"].values)
        this_sgen = sgen_p[nb_step, :].astype(np.float32)
        # this_sgen = sgen_p_init[changed_sgen].astype(np.float32)
        env_lightsim.backend._grid.update_sgens_p(changed_sgen, this_sgen)
        obs, reward, done, info = env_lightsim.step(env_lightsim.action_space())
        if reset_solver:
            env_lightsim.backend._grid.tell_solver_need_reset()
        nb_step += 1
        
    # NB lightsim2grid does not handle "static gen" because I cannot set "p" in gen in grid2op
    # so results will vary between TimeSeries and grid2op !
    # env_lightsim.backend._grid.tell_solver_need_reset()
    # env_lightsim.backend._grid.dc_pf(env_lightsim.backend.V, 1, 1e-7)
    # env_lightsim.backend._grid.get_bus_status()
    if nb_step != nb_ts:
        warnings.warn(f"only able to make {nb_step} (out of {nb_ts}) for {case_name} in grid2op. Results will not be availabe for grid2op step")
        solver_preproc_solver_time.append(None)
        g2op_speeds.append(None)
        g2op_step_time.append(None)
        ls_solver_time.append(None)
        ls_gridmodel_time.append(None)
    else:
        total_time = env_lightsim.backend._timer_preproc + env_lightsim.backend._timer_solver # + env_lightsim.backend._timer_postproc
        # total_time = env_lightsim._time_step
        solver_preproc_solver_time.append(total_time)
        g2op_speeds.append(1.0 * nb_step / total_time)
        g2op_step_time.append(1.0 * env_lightsim._time_step / nb_step)
        ls_solver_time.append(env_lightsim.backend.comp_time)
        ls_gridmodel_time.append(env_lightsim.backend.timer_gridmodel_xx_pf)
    g2op_sizes.append(env_lightsim.n_sub)
    return nb_step
        
        
if __name__ == "__main__":
    prng = np.random.default_rng(42)
    case_names_displayed = [get_env_name_displayed(el) for el in case_names]
    solver_preproc_solver_time = []
    g2op_speeds = []
    g2op_sizes = []
    g2op_step_time = []
    ls_solver_time = []
    ls_gridmodel_time = []
    
    solver_preproc_solver_time_reset = []
    g2op_speeds_reset = []
    g2op_sizes_reset = []
    g2op_step_time_reset = []
    ls_solver_time_reset = []
    ls_gridmodel_time_reset = []
    
    ts_times = []
    ts_speeds = []
    ts_sizes = []
    sa_times = []
    sa_speeds = []
    sa_sizes = []
    for case_name in tqdm(case_names):

        if not os.path.exists(case_name):
            import pandapower.networks as pn
            case = getattr(pn, os.path.splitext(case_name)[0])()
            pp.to_json(case, case_name)

        # load the case file
        case = pp.from_json(case_name)
        pp.runpp(case)  # for slack
        
        # extract reference data
        load_p_init = 1.0 * case.load["p_mw"].values
        load_q_init = 1.0 * case.load["q_mvar"].values
        gen_p_init = 1.0 * case.gen["p_mw"].values
        sgen_p_init = 1.0 * case.sgen["p_mw"].values

        res_time = 1.
        res_unit = "s"
        if len(load_p_init) <= 1000:
            # report results in ms if there are less than 1000 loads
            # only affects "verbose" printing
            res_time = 1e3
            res_unit = "ms"

        # simulate the data
        load_p, load_q, gen_p, sgen_p = get_loads_gens(load_p_init, load_q_init, gen_p_init, sgen_p_init, prng)
        if DEBUG:
            hash_fun = hashlib.blake2b(digest_size=16)
            hash_fun.update(load_p.tobytes())
            hash_fun.update(load_q.tobytes())
            hash_fun.update(gen_p.tobytes())
            hash_fun.update(sgen_p.tobytes())
            print(hash_fun.hexdigest())        
        # create the grid2op env
        nb_ts = gen_p.shape[0]
        # add slack !
        slack_gens =  np.zeros((nb_ts, case.ext_grid.shape[0]))
        if "res_ext_grid" in case:
            slack_gens += np.tile(case.res_ext_grid["p_mw"].values.reshape(1,-1), (nb_ts, 1))
        gen_p_g2op = np.concatenate((gen_p, slack_gens), axis=1)  
        # get the env        
        if WITH_PP:
            env_pp = make_grid2op_env_pp(case,
                                         case_name,
                                         load_p,
                                         load_q,
                                         gen_p_g2op,
                                         sgen_p)
            _ = env_pp.reset()
            done = False
            nb_step_pp = 0
            changed_sgen = case.sgen["in_service"].values
            while not done:
                # hack for static gen...
                env_pp.backend._grid.sgen["p_mw"] = sgen_p[nb_step_pp, :]
                obs, reward, done, info = env_pp.step(env_pp.action_space())
                nb_step_pp += 1
            if nb_step_pp != nb_ts:
                warnings.warn("Divergence even with pandapower !")
            print("Pandapower stops, lightsim starts")
            
        env_lightsim = make_grid2op_env(case,
                                        case_name,
                                        load_p,
                                        load_q,
                                        gen_p_g2op,
                                        sgen_p)
        # Perform the computation using grid2op
        reset_solver = True  # non default
        nb_step_reset = run_grid2op_env(env_lightsim, case, reset_solver,
                                        solver_preproc_solver_time_reset, 
                                        g2op_speeds_reset,
                                        g2op_step_time_reset,
                                        ls_solver_time_reset,
                                        ls_gridmodel_time_reset,
                                        g2op_sizes_reset, sgen_p, nb_ts
                                        )
        
        reset_solver = False  # default
        nb_step = run_grid2op_env(env_lightsim, case, reset_solver,
                                  solver_preproc_solver_time, 
                                  g2op_speeds,
                                  g2op_step_time,
                                  ls_solver_time,
                                  ls_gridmodel_time,
                                  g2op_sizes, sgen_p, nb_ts
                                  )
        
        # Perform the computation using TimeSerie
        env_lightsim.reset()
        time_serie = TimeSerie(env_lightsim)
        computer = time_serie.computer
        computer_ts = time_serie.computer
        v_init = env_lightsim.backend.V
        status = computer.compute_Vs(gen_p,
                                     sgen_p,
                                     load_p,
                                     load_q,
                                     v_init,
                                     env_lightsim.backend.max_it,
                                     env_lightsim.backend.tol)
        time_serie._TimeSerie__computed = True
        a_or = time_serie.compute_A()
        assert status, f"some powerflow diverge for Time Series for {case_name}: {computer.nb_solved()} "

        if VERBOSE:
            # print detailed results if needed
            print(f"For environment: {case_name} ({env_lightsim.n_sub} substations) [{computer.nb_solved()} powerflows]")
            print(f"Total time spent in \"computer\" to solve everything: {res_time*computer.total_time():.2f}{res_unit} "
                f"({computer.nb_solved() / computer.total_time():.0f} pf / s), "
                f"{1000.*computer.total_time() / computer.nb_solved():.2f} ms / pf)")
            print(f"\t - time to pre process the injections: {res_time * computer.preprocessing_time():.2f}{res_unit}")
            print(f"\t - time to perform powerflows: {res_time * computer.solver_time():.2f} {res_unit} "
                f"({computer.nb_solved() / computer.solver_time():.0f} pf / s, "
                f"{1000.*computer.solver_time() / computer.nb_solved():.2f} ms / pf)")
            print(f"In addition, it took {res_time * computer.amps_computation_time():.2f} {res_unit} to retrieve the current "
                f"from the complex voltages (in total "
                f"{computer.nb_solved() / ( computer.total_time() + computer.amps_computation_time()):.1f} "
                "pf /s, "
                f"{1000.*( computer.total_time() + computer.amps_computation_time()) / computer.nb_solved():.2f} ms / pf)")
        
        ts_times.append(computer.total_time() + computer.amps_computation_time())
        ts_speeds.append(computer.nb_solved() / (computer.total_time() + computer.amps_computation_time()) )
        ts_sizes.append(env_lightsim.n_sub)

        # Perform a securtiy analysis (up to 1000 contingencies)
        env_lightsim.reset()
        sa = ContingencyAnalysis(env_lightsim)
        for i in range(env_lightsim.n_line):
            sa.add_single_contingency(i)
            if i >= 1000:
                break
        p_or, a_or, voltages = sa.get_flows()
        computer_sa = sa.computer
        sa_times.append(computer_sa.total_time() + computer_sa.amps_computation_time())
        sa_speeds.append(computer_sa.nb_solved() / (computer_sa.total_time() + computer_sa.amps_computation_time()) )
        sa_sizes.append(env_lightsim.n_sub)
                
        # close the env
        linear_solver_used_str = solver_names[env_lightsim.backend._grid.get_solver_type()]
        env_lightsim.close()

    print("Configuration:")
    print_configuration()
    print(f"Solver used for linear algebra: {linear_solver_used_str}")
    print()
        
    print("Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (no recycling allowed, non default)")
    tab_g2op = []
    for i, nm_ in enumerate(case_names_displayed):
        tab_g2op.append((nm_,
                         ts_sizes[i],
                         1000. * g2op_step_time_reset[i] if g2op_step_time_reset[i] else None,
                         1000. / g2op_speeds_reset[i] if g2op_speeds_reset[i] else None,
                         g2op_speeds_reset[i],
                         1000. * ls_gridmodel_time_reset[i] / nb_step_reset if ls_gridmodel_time_reset[i] else None,
                         1000. * ls_solver_time_reset[i] / nb_step_reset if ls_solver_time_reset[i] else None,
                         ))
    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab_g2op,
                                          headers=["grid",
                                                   "size (nb bus)",
                                                   "avg step duration (ms)",
                                                   "time [DC + AC] (ms / pf)",
                                                   "speed (pf / s)",
                                                   "time in 'solver' (ms / pf)",
                                                   "time in 'algo' (ms / pf)",
                                                   ], 
                                          tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab_g2op)
    print()
    
    
    print("Results using grid2op.steps (288 consecutive steps, only measuring 'dc pf [init] + ac pf') (recyling allowed, default)")
    tab_g2op = []
    for i, nm_ in enumerate(case_names_displayed):
        tab_g2op.append((nm_,
                         ts_sizes[i],
                         1000. * g2op_step_time[i] if g2op_step_time[i] else None,
                         1000. / g2op_speeds[i] if g2op_speeds[i] else None,
                         g2op_speeds[i],
                         1000. * ls_gridmodel_time[i] / nb_step if ls_gridmodel_time[i] else None,
                         1000. * ls_solver_time[i] / nb_step if ls_solver_time[i] else None,
                         ))
    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab_g2op,
                                          headers=["grid",
                                                   "size (nb bus)",
                                                   "avg step duration (ms)",
                                                   "time [DC + AC] (ms / pf)",
                                                   "speed (pf / s)",
                                                   "time in 'gridmodel' (ms / pf)",
                                                   "time in 'pf algo' (ms / pf)",
                                                   ], 
                                          tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab_g2op)
    print()
        
    print("Results for TimeSeries (288 consecutive steps)")
    tab_ts = []
    for i, nm_ in enumerate(case_names_displayed):
        tab_ts.append((nm_, ts_sizes[i], 1000. / ts_speeds[i] if ts_speeds[i] else None, ts_speeds[i]))
    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab_ts,
                                          headers=["grid", "size (nb bus)", "time (ms / pf)", "speed (pf / s)"], 
                                          tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab_ts)
    print()

    print("Results for Contingency Analysis (up to 1000 contingencies - or all the powerlines)")
    tab_sa = []
    for i, nm_ in enumerate(case_names_displayed):
        tab_sa.append((nm_, sa_sizes[i], 1000. / sa_speeds[i] if sa_speeds[i] else None, sa_speeds[i]))
    if TABULATE_AVAIL:
        res_use_with_grid2op_2 = tabulate(tab_sa,
                                          headers=["grid", "size (nb bus)", "time (ms / cont.)", "speed (cont. / s)"], 
                                          tablefmt="rst")
        print(res_use_with_grid2op_2)
    else:
        print(tab_sa)
    print()
        
    if MAKE_PLOT:
        # make the plot summarizing all results
        plt.plot(g2op_sizes, solver_preproc_solver_time, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Time taken (s)")
        plt.title(f"Time to compute {g2op_sizes[0]} powerflows using Grid2Op.step (dc pf [init] + ac pf)")
        plt.show()

        plt.plot(g2op_sizes, g2op_speeds, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Speed (pf / s)")
        plt.title(f"Computation speed using Grid2Op.step (dc pf [init] + ac pf)")
        plt.yscale("log")
        plt.show()

        plt.plot(g2op_sizes, ls_solver_time, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Speed (solver time)")
        plt.title(f"Computation speed for solving the powerflow only")
        plt.yscale("log")
        plt.show()

        plt.plot(g2op_sizes, ls_gridmodel_time, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Speed (solver time)")
        plt.title(f"Computation speed for solving the powerflow only")
        plt.yscale("log")
        plt.show()
        
        # make the plot summarizing all results
        plt.plot(ts_sizes, ts_times, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Time taken (s)")
        plt.title(f"Time to compute {computer_ts.nb_solved()} powerflows using TimeSeries")
        plt.show()

        plt.plot(ts_sizes, ts_speeds, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Speed (pf / s)")
        plt.title(f"Computation speed for TimeSeries")
        plt.yscale("log")
        plt.show()
        
        # make the plot summarizing all results
        plt.plot(sa_sizes, [1000. / el for el in sa_speeds], linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Time taken per contingency (ms)")
        plt.title(f"Average time per contingencies")
        plt.show()

        plt.plot(sa_sizes, sa_speeds, linestyle='solid', marker='+', markersize=8)
        plt.xlabel("Size (number of substation)")
        plt.ylabel("Speed (contingency / s)")
        plt.title(f"Computation speed for Security Analysis")
        plt.yscale("log")
        plt.show()
