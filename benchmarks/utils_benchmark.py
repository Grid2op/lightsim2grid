# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import time
import re
import numpy as np
import argparse
import datetime
import importlib
try:
    from tqdm import tqdm
    from grid2op.Environment import MultiMixEnvironment # type: ignore
except ImportError:
    # above packages not mandatory for compare_lightsim2grid_pypowsybl.py
    pass

try:
    from lightsim2grid import LightSimBackend
except ImportError:
    LightSimBackend = None


def get_env_name_displayed(env_name):
    res = re.sub("^l2rpn_", "", env_name)
    res = re.sub("_small$", "", res)
    res = re.sub("_large$", "", res)
    res = re.sub("\\.json$", "", res)
    return res


def print_res(env_klu, env_pp,
              nb_ts_klu, nb_ts_pp,
              time_klu, time_pp,
              aor_klu, aor_pp,
              gen_p_klu, gen_p_pp,
              gen_q_klu, gen_q_pp):
    print("Overall speed-up of KLU vs pandapower (for grid2opbackend) {:.2f}\n".format(time_pp / time_klu))
    print("PyKLU Backend {} time steps in {}s ({:.2f} it/s)".format(nb_ts_klu, time_klu, nb_ts_klu/time_klu))
    print("\tTime apply act: {:.2f}ms".format(1000. * env_klu._time_apply_act / nb_ts_klu))
    print("\tTime powerflow: {:.2f}ms".format(1000. * env_klu._time_powerflow / nb_ts_klu))
    print("\tTime extract observation: {:.2f}ms".format(1000. * env_klu._time_extract_obs / nb_ts_klu))

    print("Pandapower Backend {} time steps in {}s ({:.2f} it/s)".format(nb_ts_pp, time_pp, nb_ts_pp/time_pp))
    print("\tTime apply act: {:.2f}ms".format(1000. * env_pp._time_apply_act / nb_ts_pp))
    print("\tTime powerflow: {:.2f}ms".format(1000. * env_pp._time_powerflow / nb_ts_pp))
    print("\tTime extract observation: {:.2f}ms".format(1000. * env_pp._time_extract_obs / nb_ts_pp))

    print("Absolute value of the difference (max) for aor: {}".format(np.max(np.abs(aor_klu - aor_pp))))
    print("Absolute value of the difference (max) for gen_p: {}".format(np.max(np.abs(gen_p_klu - gen_p_pp))))
    print("Absolute value of the difference (max) for gen_q: {}".format(np.max(np.abs(gen_q_klu - gen_q_pp))))


def get_rest(env_pp, env_KLU, env_SLU, env_GS):
    pass


def run_env(
    env,
    max_ts,
    agent,
    chron_id=None,
    keep_forecast=False,
    with_type_solver=False,
    env_seed=None,
    detailed_ls_output=False,
    is_dc=False):
    pbar_desc = None
    if with_type_solver:
        try:
            pbar_desc = f"{env.backend._grid.get_solver_type()}".split(".")[1]
        except Exception:
            # it just means I will not be able to print the fancy name on the progress bar...
            pass

    nb_rows = min(env.chronics_handler.max_timestep(), max_ts)
    if nb_rows == -1:
        # -1 indicated "infinite data"
        nb_rows = max_ts
    aor = np.zeros((nb_rows, env.n_line))
    gen_p = np.zeros((nb_rows, env.n_gen))
    gen_q = np.zeros((nb_rows, env.n_gen))
    need_reset = False
    if isinstance(env, MultiMixEnvironment):
        # get the first (in alphabetical order) env in case of multimix
        env = env[next(iter(sorted(env.keys())))]
    if env_seed is not None:
        env.seed(env_seed)
        need_reset = True

    if chron_id is not None:
        # reset the environment
        env.chronics_handler.tell_id(chron_id-1)
        # deactivate the forecast (not used here)
        if not keep_forecast:
            env.deactivate_forecast()
        need_reset = True
    
    if need_reset:
        obs = env.reset()
    else:
        obs = env.get_obs()

    # don't forget to reset the timers
    env.backend.comp_time = 0.
    # it's not 0. because a powerflow is run when the backend is reset, and this is one more powerflow than the
    # number of steps

    done = False
    reward = env.reward_range[0]
    nb_ts = 0
    prev_act = None
    time_ls = 0.
    ls_timer_Fx = 0.
    ls_timer_solve = 0.
    ls_timer_initialize = 0.
    ls_timer_check = 0.
    ls_timer_dSbus = 0.
    ls_timer_fillJ = 0.
    ls_timer_Va_Vm = 0.
    ls_timer_pre_proc = 0.
    ls_timer_total_nr = 0.
    ls_timer_factor = 0.
    ls_timer_refactor = 0.
    ls_timer_per_refactor = 0.
    ls_timer_per_solve = 0.
    ls_timer_per_scale = 0.
    ls_timer_per_mismatch = 0.
    if LightSimBackend is None or not isinstance(env.backend, LightSimBackend):
        detailed_ls_output = False
    beg_ = time.perf_counter()
    with tqdm(total=nb_rows, desc=pbar_desc) as pbar:
        while not done:
            act = agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            aor[nb_ts, :] = obs.a_or
            gen_p[nb_ts, :] = obs.prod_p
            gen_q[nb_ts, :] = obs.prod_q
            nb_ts += 1
            if detailed_ls_output:
                ls_grid = env.backend._grid
                if is_dc:
                    timers_jac = ls_grid.get_dc_solver().get_timers_jacobian()
                    nr_iter = 1
                else:
                    timers_jac = ls_grid.get_solver().get_timers_jacobian()
                    nr_iter =  ls_grid.get_solver().get_nb_iter()
                ls_timer_Fx += timers_jac.timer_Fx
                ls_timer_solve += timers_jac.timer_solve
                ls_timer_initialize += timers_jac.timer_initialize
                ls_timer_check += timers_jac.timer_check
                ls_timer_dSbus += timers_jac.timer_dSbus
                ls_timer_fillJ += timers_jac.timer_fillJ
                ls_timer_Va_Vm += timers_jac.timer_Va_Vm
                ls_timer_pre_proc += timers_jac.timer_pre_proc
                ls_timer_total_nr += timers_jac.timer_total_nr
                ls_timer_refactor += timers_jac.timer_refactor
                ls_timer_factor += timers_jac.timer_factor
                if abs(timers_jac.timer_factor) > 1e-8:
                    # factor is used
                    nb_denom = nr_iter - 1
                else:
                    # refactor is used accross all iterations
                    nb_denom = nr_iter
                ls_timer_per_refactor += timers_jac.timer_refactor / nb_denom 
                ls_timer_per_solve += timers_jac.timer_solve / nr_iter  # solve is called each iter
                ls_timer_per_scale += timers_jac.timer_scale
                ls_timer_per_mismatch += timers_jac.timer_mismatch
            pbar.update(1)
            if nb_ts >= max_ts or done:
                break
            # if np.sum(obs.line_status) < obs.n_line - 1 * (nb_ts % 2 == 1):
            #     print("There is a bug following action; {}".format(act))
            prev_act = act  # noqa: F841
    end_ = time.perf_counter()
    total_time = end_ - beg_
    if detailed_ls_output:
        time_ls = env.backend._timer_solver
        print("--------------------------------------")
        print("Detailed lightsim2grid timings: ")
        print(f"Total time {time_ls * 1000.:.2e}ms => {time_ls / nb_ts * 1e3:.2e} ms/pf | {nb_ts / time_ls:.0f} pf /s")
        print(f"Total time spent in the solver: {1e3 * ls_timer_total_nr:.2e} ms "
            f"({100. * ls_timer_total_nr / time_ls:.0f}% of total time spent in lightsim2grid)")
        print(f"Total time spent in linear solver {1e3 * (ls_timer_initialize + ls_timer_factor + ls_timer_refactor + ls_timer_solve):.2e} ms "
            f"({100. * (ls_timer_initialize + ls_timer_factor + ls_timer_refactor + ls_timer_solve) / ls_timer_total_nr:.0f}% of total time in the solver)")
        print(f"\t Time to pre process Ybus, Sbus etc.: {1e3 * ls_timer_pre_proc:.2e} ms ({100. * ls_timer_pre_proc / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to initialize linear solver {1e3 * ls_timer_initialize:.2e} ms ({100. * ls_timer_initialize / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to factor the linear solver {1e3 * ls_timer_factor:.2e} ms ({100. * ls_timer_factor / ls_timer_total_nr:.0f} % of time in solver)")
        if not is_dc:
            print(f"\t Time to fill the Jacobian {1e3 * ls_timer_fillJ:.2e} ms ({100. * ls_timer_fillJ / ls_timer_total_nr:.0f} % of time in solver)")
        print(f"\t Time to refactor the Jacobian linear system: {1e3 * ls_timer_refactor:.2e} ms ({100. * ls_timer_refactor / ls_timer_total_nr:.0f} % of time in solver)")
        if not is_dc:
            print(f"\t\t This equates to {1e3 * ls_timer_per_refactor:.2e} ms per refactor")
        print(f"\t Time to solve the Jacobian linear system: {1e3 * ls_timer_solve:.2e} ms ({100. * ls_timer_solve / ls_timer_total_nr:.0f} % of time in solver)")
        if not is_dc:
            print(f"\t\t This equates to {1e3 * ls_timer_per_solve:.2e} ms per solve")
            print(f"\t Time to evaluate p,q mismmatch at each bus {1e3*ls_timer_Fx:.2e} ms ({100. * ls_timer_Fx / ls_timer_total_nr:.0f} % of time in solver)")
            print(f"\t Time to evaluate cvg criteria {1e3*ls_timer_check:.2e} ms ({100. * ls_timer_check / ls_timer_total_nr:.0f} % of time in solver)")
            print(f"\t Time to build mismatch vector {1e3*ls_timer_per_mismatch:.2e} ms ({100. * ls_timer_per_mismatch / ls_timer_total_nr:.0f} % of time in solver)")
            print(f"\t Time to scale mismatch vector {1e3*ls_timer_per_scale:.2e} ms ({100. * ls_timer_per_scale / ls_timer_total_nr:.0f} % of time in solver)")
        else:
            # in DC timer_mismatch is used for the post processing
            print(f"\t Time in the post processing: {1e3*ls_timer_per_mismatch:.2e} ms ({100. * ls_timer_per_mismatch / ls_timer_total_nr:.0f} % of time in solver)")
        print("--------------------------------------\n")

    return nb_ts, total_time, aor, gen_p, gen_q


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def print_configuration(
    pypowbk_error=True,
    pypowsybl_error=True
    ):
    res = []
    print()
    tmp = f"- date: {datetime.datetime.now():%Y-%m-%d %H:%M %z} {time.localtime().tm_zone}"
    res.append(tmp)
    print(tmp)
    try:
        import platform
        tmp = f"- system: {platform.system()} {platform.release()}"
        res.append(tmp)
        print(tmp)
    except ImportError:
        tmp = "- system: please install the `platform` to have this information"
        res.append(tmp)
        print(tmp)

    try:
        import distro
        tmp = (f"- OS: {distro.linux_distribution(full_distribution_name=False)[0]} "
               f"{distro.linux_distribution(full_distribution_name=False)[1]}")
        res.append(tmp)
        print(tmp)
    except ImportError:
        tmp = ("- OS: please install the `distro` to have this information")
        res.append(tmp)
        print(tmp)

    try:
        import cpuinfo
        info_ = cpuinfo.get_cpu_info()
        tmp = (f"- processor: {info_['brand_raw']}")
        res.append(tmp)
        print(tmp)
        tmp = (f"- python version: {info_['python_version']}")
        res.append(tmp)
        print(tmp)

    except ImportError:
        tmp = ("- processor: please install the `py-cpuinfo` to have this information")
        res.append(tmp)
        print(tmp)
        tmp = ("- python version: please install the `py-cpuinfo` to have this information")
        res.append(tmp)
        print(tmp)

    import pandas as pd
    import lightsim2grid
    tmp = (f"- numpy version: {np.__version__}")
    res.append(tmp)
    print(tmp)
    tmp = (f"- pandas version: {pd.__version__}")
    res.append(tmp)
    print(tmp)
    try:
        import pandapower as pp
        tmp = (f"- pandapower version: {pp.__version__}")
        res.append(tmp)
        print(tmp)
    except ImportError:
        # pandapower not used for compare_lightsim2grid_pypowsybl
        pass
    if pypowbk_error is None:
        # print both pypowsybl and pypowsybl2grid info
        tmp = (f"- pypowsybl version: {importlib.metadata.version('pypowsybl')}")
        res.append(tmp)
        print(tmp)
        tmp = (f"- pypowsybl2grid version: {importlib.metadata.version('pypowsybl2grid')}")
        res.append(tmp)
        print(tmp)
    elif pypowsybl_error is None:
        # print only pypowsybl info
        tmp = (f"- pypowsybl version: {importlib.metadata.version('pypowsybl')}")
        res.append(tmp)
        print(tmp)
    try:
        import grid2op # type: ignore
        tmp = (f"- grid2op version: {grid2op.__version__}")
        res.append(tmp)
        print(tmp)
    except ImportError:
        # grid2op not used for compare_lightsim2grid_pypowsybl
        pass
    tmp = (f"- lightsim2grid version: {lightsim2grid.__version__}")
    res.append(tmp)
    print(tmp)
    try:
        from lightsim2grid import compilation_options
        tmp = ("- lightsim2grid extra information: ")
        res.append(tmp)
        print(tmp)
        print()
        tmp = (f"\t- klu_solver_available: {compilation_options.klu_solver_available} ")
        res.append(tmp)
        print(tmp)
        tmp = (f"\t- nicslu_solver_available: {compilation_options.nicslu_solver_available} ")
        res.append(tmp)
        print(tmp)
        tmp = (f"\t- cktso_solver_available: {compilation_options.cktso_solver_available} ")
        res.append(tmp)
        print(tmp)
        tmp = (f"\t- compiled_march_native: {compilation_options.compiled_march_native} ")
        res.append(tmp)
        print(tmp)
        tmp = (f"\t- compiled_o3_optim: {compilation_options.compiled_o3_optim} ")
        res.append(tmp)
        print(tmp)
    except ImportError:
        # before it was introduced
        pass
    print()
    return '\n'.join(res)
