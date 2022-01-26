# Copyright (c) 2020-2022, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import pdb
import pandapower as pp
import os
import numpy as np
import grid2op
from grid2op.Runner import Runner
from grid2op.Environment import MultiMixEnvironment, Environment
from lightsim2grid import LightSimBackend
from grid2op.MakeEnv.UpdateEnv import _hash_env
import warnings


def fix_pp_net(pp_net):
    ## sgen
    if "min_p_mw" in pp_net.sgen:
        min_p_mw = pp_net.sgen["min_p_mw"].values
        min_p_mw[~np.isfinite(min_p_mw)] = 0.
        pp_net.sgen["min_p_mw"][:] = min_p_mw

    if "max_p_mw" in pp_net.sgen:
        max_p_mw = pp_net.sgen["max_p_mw"].values
        max_p_mw[~np.isfinite(max_p_mw)] = 0.
        pp_net.sgen["max_p_mw"][:] = max_p_mw

    if "min_q_mvar" in pp_net.sgen:
        min_q_mvar = pp_net.sgen["min_q_mvar"].values
        min_q_mvar[~np.isfinite(min_q_mvar)] = 0.
        pp_net.sgen["min_q_mvar"][:] = min_q_mvar

    if "max_q_mvar" in pp_net.sgen:
        max_q_mvar = pp_net.sgen["max_q_mvar"].values
        max_q_mvar[~np.isfinite(max_q_mvar)] = 0.
        pp_net.sgen["max_q_mvar"][:] = max_q_mvar

    ## slack bus
    if np.any(pp_net.gen["slack"].values):
        # in this case the ext grid is not taken into account, i raise a warning if
        # there is one
        slack_bus_ids = pp_net.ext_grid["bus"].values
        if pp_net.ext_grid.shape[0] >= 1:
            del pp_net.ext_grid
        print("INFO: ext_grid deleted because pp_net.gen has some slacks")
    else:
        slack_bus_ids = pp_net.ext_grid["bus"].values
        if pp_net.ext_grid.shape[0] >= 2:
            raise RuntimeError("We found multiple slack buses in the ext_grid. We cannot modify "
                               "the grid automatically when this is the case")

        if np.all(np.isin(slack_bus_ids, pp_net.gen["bus"].values)):
            # all slack buses have a generator connected to them
            # so i assume it was just a computation artifact, and assign these generators as slack buses
            slack_gen_ids = np.isin(pp_net.gen["bus"].values, slack_bus_ids)  # id of generators connected to slack bus
            slack_gen_ids = np.where(slack_gen_ids)[0]  # keep only the id of the generators
            if "slack_weight" in pp_net.gen:
                slack_coeff = pp_net.gen["slack_weight"].values[slack_gen_ids]
            print("INFO: flag \"slack=True\" added in the original grid generators")
        else:
            print("INFO: new generators added to serve as slack !")
            nb_slack = len(slack_bus_ids)
            if "slack_weight" in pp_net.ext_grid:
                slack_coeff = 1.0 * pp_net.ext_grid["slack_weight"].values
            else:
                slack_coeff = np.ones(nb_slack)

            slack_coeff_norm = slack_coeff / slack_coeff.sum()
            slack_gen_ids = np.arange(nb_slack) + pp_net.gen.shape[0]
            slack_contrib = (np.sum(pp_net.gen["p_mw"]) - np.sum(pp_net.load["p_mw"]) ) * slack_coeff_norm
            nb_init_gen = pp_net.gen.shape[0]
            for gen_added_id in range(nb_slack):
                bus_id = slack_bus_ids[gen_added_id]
                vm_pu = pp_net.ext_grid["vm_pu"].values[gen_added_id]
                p_mw = slack_contrib[gen_added_id]
                min_q_mvar = -9999.
                max_q_mvar = +9999.
                pp.create_gen(net=pp_net,
                              name=f"gen_{bus_id}_{gen_added_id+nb_init_gen}",
                              bus=bus_id,
                              p_mw=p_mw,
                              vm_pu=vm_pu,
                              min_q_mvar=min_q_mvar,
                              max_q_mvar=max_q_mvar,
                              slack=True,
                              )
                warnings.warn("slack_weight not taken into account !")
        del pp_net.ext_grid

    ## trafo
    if pp_net.trafo.shape[0] > 0:
        tap_neutral = pp_net.trafo["tap_neutral"].values
        tap_neutral[~np.isfinite(tap_neutral)] = 0.
        pp_net.trafo["tap_neutral"][:] = tap_neutral

        tap_step_percent = pp_net.trafo["tap_step_percent"].values
        tap_step_percent[~np.isfinite(tap_step_percent)] = 0.
        pp_net.trafo["tap_step_percent"][:] = tap_step_percent

        tap_pos = pp_net.trafo["tap_pos"].values
        tap_pos[~np.isfinite(tap_pos)] = 0.
        pp_net.trafo["tap_pos"][:] = tap_pos

        shift_degree = pp_net.trafo["shift_degree"].values
        shift_degree[~np.isfinite(shift_degree)] = 0.
        pp_net.trafo["shift_degree"][:] = shift_degree

        is_tap_broken = ~np.array([el == "hv" or el == "lv" for el in pp_net.trafo["tap_side"].values])
        pp_net.trafo["tap_side"][is_tap_broken] = False

        tap_step_degree = pp_net.trafo["tap_step_degree"].values
        tap_step_degree[~np.isfinite(tap_step_degree)] = 0.
        pp_net.trafo["tap_step_degree"][:] = tap_step_degree


def check_env(env, new_grid_path):
    # check all warnings are removed
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        env_nowarn = grid2op.make(env_path,
                                  grid_path=new_grid_path,
                                  backend=LightSimBackend(),
                                  _add_to_name="_no_warn",)

    # basic check to "make sure" I did not mess with the file
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env_fix = grid2op.make(env_path,
                               grid_path=new_grid_path,
                               backend=LightSimBackend(),
                               _add_to_name="_new",
                               test=True)
    assert env_fix._init_grid_path == new_grid_path
    assert type(env_fix).__name__ != type(env).__name__

    # check the grid is loaded the same way
    assert env.n_gen == env_fix.n_gen
    assert len(env.name_gen) == len(env_fix.name_gen)
    assert np.all(env.name_gen == env_fix.name_gen)
    dict_ref = env.cls_to_dict()
    del dict_ref["env_name"]
    dict_new = env_fix.cls_to_dict()
    del dict_new["env_name"]
    assert dict_new == dict_ref

    # check it leads to the same results
    runner = Runner(**env.get_params_for_runner())
    runner_fix = Runner(**env_fix.get_params_for_runner())

    nb_episode = 10
    max_iter = 288*2
    res = runner.run(nb_episode=nb_episode,
                     pbar=True,
                     max_iter=max_iter,
                     env_seeds=[0] * nb_episode)
    res_fix = runner_fix.run(nb_episode=nb_episode,
                             pbar=True,
                             max_iter=max_iter,
                             env_seeds=[0] * nb_episode)

    for i, (path_, nm_, cum_reward, nb_step, total_step) in enumerate(res):
        path_fix, nm_fix, cum_reward_fix, nb_step_fix, total_step_fix = res_fix[i]
        assert cum_reward == cum_reward_fix, f"{cum_reward_fix = } vs {cum_reward = } for {nm_}"
        assert nb_step == nb_step_fix, f"{nb_step_fix = } vs {nb_step = }"


    assert res == res_fix


def fix_for_unit_env(env_path, init_grid_path, env, env_pp):
    pp_net = pp.from_json(init_grid_path)

    # check that i can perform the tests
    runner = Runner(**env.get_params_for_runner())
    runner_pp = Runner(**env_pp.get_params_for_runner())
    nb_episode = 4
    max_iter = 100
    res = runner.run(nb_episode=nb_episode, pbar=True, max_iter=max_iter, env_seeds=[0] * nb_episode)
    res_pp = runner_pp.run(nb_episode=nb_episode, pbar=True, max_iter=max_iter, env_seeds=[0] * nb_episode)
    assert res == res_pp, "pandapower and lightsim2grid does not give the same result. Stopping there"


    # back it up (if not done already)
    grid_backup = os.path.join(env_path, "grid.json.save")
    if not os.path.exists(grid_backup):
        pp.to_json(pp_net, grid_backup)

    # fix it
    fix_pp_net(pp_net)

    # check it works
    new_grid_path = os.path.abspath(f"./{env_name}_grid.json")
    pp.to_json(pp_net, new_grid_path)
    check_env(env, new_grid_path)

    # save the new grid
    pp.to_json(pp_net, init_grid_path)
    return new_grid_path


if __name__ == "__main__":
    env_name = "rte_case14_redisp"
    env_name = "l2rpn_case14_sandbox"
    env_name = "l2rpn_2019"
    env_name = "rte_case14_realistic"
    env_name = "l2rpn_wcci_2020"
    env_name = "l2rpn_neurips_2020_track1_small"
    env_name = "l2rpn_neurips_2020_track1_large"
    env_name = "l2rpn_icaps_2021_small"
    env_name = "l2rpn_icaps_2021_large"
    env_name = "l2rpn_neurips_2020_track2_small"
    env_name = "l2rpn_neurips_2020_track2_large"

    # env_name = "l2rpn_neurips_2020_track2_small"

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(env_name, backend=LightSimBackend())
        env_pp = grid2op.make(env_name)

    env_path = env.get_path_env()
    final_path = env_path = env.get_path_env()
    if isinstance(env, MultiMixEnvironment):
        key = next(iter(env.keys()))
        super_env = env
        super_env_pp = env_pp
        env = super_env[key]
        env_pp = super_env_pp[key]
        
        env_path = env.get_path_env()
        final_path = super_env.get_path_env()
    elif not isinstance(env, Environment):
        raise RuntimeError("Unuspported environment type !")
    
    init_grid_path = os.path.join(env_path, "grid.json")
    new_grid_path = fix_for_unit_env(env_path, init_grid_path, env, env_pp)

    # compute the hash:
    print(f"hash for {env_name} {_hash_env(final_path).hexdigest()}")
    print(f"new file for {env_name} {new_grid_path}")