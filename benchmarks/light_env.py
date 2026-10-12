# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

"""Compare the light environment (:class:`lightsim2grid.lightEnv.LightEnv`) with a grid2op
environment using :class:`lightsim2grid.LightSimBackend`, on the same grid and the same
chronics.

Three workloads are timed:

- ``do nothing``: an agent that never acts, ``env.step(do_nothing)`` vs ``light_env.step(0)``;
- ``topology``: every ``--act_every`` steps the agent plays a random unitary ``set_bus``
  action (the same sequence of action ids on both sides), do nothing otherwise. Such an agent
  quickly ends its episode, so ``--nb_episode_topo`` episodes are played (cycling through the
  chronics, each with its own random sequence);
- ``lookahead``: at every step ``--nb_candidates`` actions are evaluated one step ahead,
  ``obs.simulate(act)`` for grid2op vs ``light_env.copy().step(act_id)``, then the episode
  moves on with do nothing. Only the evaluation is timed.

Usage (from the ``benchmarks`` folder)::

    python light_env.py                    # l2rpn_case14_sandbox, the test chronics shipped with grid2op
    python light_env.py --no_test          # the full dataset (downloaded by grid2op)
"""

import argparse
import time
import warnings

import numpy as np
from tqdm import tqdm

import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.lightEnv import LightEnv, Protections

from utils_benchmark import print_configuration, str2bool

ENV_NAME = "l2rpn_case14_sandbox"


def make_g2op_env(env_name, test):
    param = Parameters()
    # the light environment only has the "soft" overflow protection (a line in overflow for
    # more than NB_TIMESTEP_OVERFLOW_ALLOWED steps is disconnected), not the instantaneous one
    param.HARD_OVERFLOW_THRESHOLD = 1e6
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(env_name, test=test, backend=LightSimBackend(), param=param)
    return env


def make_light_env(env, actions):
    """a light env on the grid of `env` (at the start of an episode) with the same thermal
    limits, protections and cooldowns"""
    light_env = LightEnv(env.backend._grid)
    param = env.parameters
    th_lim_ka = env.get_thermal_limit().astype(float) * 1e-3  # grid2op: A, light env: kA
    protections = Protections()
    protections.set_thermal_limit_or(th_lim_ka)
    protections.set_thermal_limit_ex(th_lim_ka * env.backend.lines_or_pu_to_kv / env.backend.lines_ex_pu_to_kv)
    protections.set_max_line_time_step_overflow(np.full(env.n_line, param.NB_TIMESTEP_OVERFLOW_ALLOWED,
                                                        dtype=np.int32))
    light_env.protections = protections
    light_env.nb_timestep_cooldown_sub = param.NB_TIMESTEP_COOLDOWN_SUB
    light_env.nb_timestep_cooldown_line = param.NB_TIMESTEP_COOLDOWN_LINE
    light_env.nb_timestep_reconnection = param.NB_TIMESTEP_RECONNECTION
    light_env.init_actions(actions)
    return light_env


def assign_chronic(light_env, env, max_ts):
    """give the light env the injections of the chronic `env` is currently playing"""
    data = env.chronics_handler.real_data.data
    nb_row = min(data.load_p.shape[0], max_ts + 1)
    light_env.assign_time_series(
        data.load_p[:nb_row].astype(float),
        data.load_q[:nb_row].astype(float),
        data.prod_p[:nb_row].astype(float),
        (data.prod_v[:nb_row] / env.backend.prod_pu_to_kv).astype(float),  # the light env wants pu
        np.full((nb_row, env.n_storage), np.nan),
        np.full((nb_row, env.n_shunt), np.nan),
        np.full((nb_row, env.n_shunt), np.nan),
        np.full((nb_row, 0), np.nan),
        np.full((nb_row, 0), np.nan),
    )


def action_sequence(rng, max_ts, nb_actions, act_every):
    """the action id played at each step: 0 (do nothing) except every `act_every` steps"""
    res = np.zeros(max_ts + 1, dtype=int)
    if act_every > 0 and nb_actions > 1:
        idx = np.arange(act_every, max_ts + 1, act_every)
        res[idx] = rng.integers(1, nb_actions, size=idx.shape[0])
    return res


class Timings:
    def __init__(self):
        self.nb_step = 0
        self.nb_episode = 0
        self.nb_game_over = 0
        self.total = 0.   # wall clock, seen from python
        self.inner = 0.   # time spent in the env itself (grid2op: env._time_step, light env: step_time)
        self.powerflow = 0.

    def ms_per_step(self, attr="total"):
        return 1e3 * getattr(self, attr) / max(self.nb_step, 1)

    def step_per_s(self):
        return self.nb_step / self.total if self.total > 0. else float("nan")


def run_g2op(env, g2op_actions, act_ids_per_chron, chron_ids, max_ts):
    res = Timings()
    for chron_id, act_ids in tqdm(zip(chron_ids, act_ids_per_chron), total=len(chron_ids), desc="grid2op"):
        env.set_id(chron_id)
        env.reset()
        res.nb_episode += 1
        done = False
        nb_ts = 0
        beg_ = time.perf_counter()
        while not done and nb_ts < max_ts:
            nb_ts += 1
            _, _, done, info = env.step(g2op_actions[act_ids[nb_ts]])
        res.total += time.perf_counter() - beg_
        res.inner += env._time_step
        res.powerflow += env._time_powerflow
        res.nb_step += nb_ts
        res.nb_game_over += int(done and len(info["exception"]) > 0)
    return res


def run_light(env, light_env, act_ids_per_chron, chron_ids, max_ts):
    res = Timings()
    for chron_id, act_ids in tqdm(zip(chron_ids, act_ids_per_chron), total=len(chron_ids), desc="light env"):
        env.set_id(chron_id)
        env.reset()  # only to read the chronic
        assign_chronic(light_env, env, max_ts)
        light_env.reset()
        res.nb_episode += 1
        done = False
        nb_ts = 0
        beg_ = time.perf_counter()
        while not done and nb_ts < max_ts:
            nb_ts += 1
            _, _, done, truncated, info = light_env.step(int(act_ids[nb_ts]))
        res.total += time.perf_counter() - beg_
        res.inner += light_env.step_time
        res.powerflow += light_env.protections.powerflow_time
        # the last step of a chronic returns done without playing anything
        res.nb_step += nb_ts - int(done and info.get("success") == "true")
        res.nb_game_over += int(done and info.get("failure") == "true")
    return res


def run_lookahead_g2op(env, g2op_actions, candidates, chron_ids, max_ts):
    res = Timings()
    do_nothing = g2op_actions[0]
    for chron_id in tqdm(chron_ids, desc="grid2op simulate"):
        env.set_id(chron_id)
        obs = env.reset()
        res.nb_episode += 1
        done = False
        nb_ts = 0
        nb_row = min(max_ts, env.chronics_handler.max_timestep())
        while not done and nb_ts < nb_row - 1:
            beg_ = time.perf_counter()
            for act_id in candidates:
                obs.simulate(g2op_actions[act_id])
            res.total += time.perf_counter() - beg_
            res.nb_step += len(candidates)
            obs, _, done, _ = env.step(do_nothing)
            nb_ts += 1
    return res


def run_lookahead_light(env, light_env, candidates, chron_ids, max_ts):
    res = Timings()
    for chron_id in tqdm(chron_ids, desc="light env copy + step"):
        env.set_id(chron_id)
        env.reset()
        assign_chronic(light_env, env, max_ts)
        light_env.reset()
        res.nb_episode += 1
        done = False
        nb_ts = 0
        # the step reaching max_step ends the episode without a powerflow: not a lookahead
        while not done and nb_ts < light_env.max_step - 2:
            beg_ = time.perf_counter()
            for act_id in candidates:
                light_env.copy().step(int(act_id))
            res.total += time.perf_counter() - beg_
            res.nb_step += len(candidates)
            _, _, done, _, _ = light_env.step(0)
            nb_ts += 1
    return res


def print_table(title, rows):
    print(f"\n{title}\n")
    header = ("", "steps", "game over", "steps / s", "ms / step", "ms / step (in env)", "ms / step (powerflow)")
    table = [header]
    for name, t in rows:
        table.append((name,
                      f"{t.nb_step}",
                      f"{t.nb_game_over}",
                      f"{t.step_per_s():.0f}",
                      f"{t.ms_per_step():.3f}",
                      f"{t.ms_per_step('inner'):.3f}" if t.inner > 0. else "",
                      f"{t.ms_per_step('powerflow'):.3f}" if t.powerflow > 0. else ""))
    widths = [max(len(r[i]) for r in table) for i in range(len(header))]
    sep = "+" + "+".join("-" * (w + 2) for w in widths) + "+"
    print(sep)
    for i, r in enumerate(table):
        print("| " + " | ".join(c.ljust(w) for c, w in zip(r, widths)) + " |")
        if i == 0:
            print(sep.replace("-", "="))
    print(sep)
    if len(rows) == 2 and rows[1][1].total > 0.:
        print(f"speed-up: x{rows[0][1].ms_per_step() / rows[1][1].ms_per_step():.1f}")


def main(env_name=ENV_NAME, test=True, max_ts=10_000, nb_chronics=3, act_every=10, nb_episode_topo=50,
         nb_candidates=10, seed=0):
    env = make_g2op_env(env_name, test)
    env.seed(seed)
    env.reset()  # the grid of the backend is the initial state of the episode

    # positions in the list of chronics (what env.set_id takes)
    chron_ids = list(range(min(nb_chronics, len(env.chronics_handler.real_data.available_chronics()))))

    # do nothing first, then the unitary topologies (only set_bus: what the light env models)
    g2op_actions = [env.action_space()] + list(env.action_space.get_all_unitary_topologies_set(env.action_space))
    light_env = make_light_env(env, g2op_actions)
    print(f"{env_name} ({'test' if test else 'full'} chronics): {len(chron_ids)} chronics, "
          f"{len(g2op_actions)} actions (do nothing + unitary set_bus)")

    rng = np.random.default_rng(seed)
    no_act = [np.zeros(max_ts + 1, dtype=int) for _ in chron_ids]
    topo_chron_ids = [chron_ids[i % len(chron_ids)] for i in range(nb_episode_topo)]
    topo_act = [action_sequence(rng, max_ts, len(g2op_actions), act_every) for _ in topo_chron_ids]
    candidates = rng.choice(np.arange(1, len(g2op_actions)), size=min(nb_candidates, len(g2op_actions) - 1),
                            replace=False)

    print_configuration()

    env.deactivate_forecast()  # not used by the first two workloads, grid2op would still read them
    g2op_dn = run_g2op(env, g2op_actions, no_act, chron_ids, max_ts)
    light_dn = run_light(env, light_env, no_act, chron_ids, max_ts)
    g2op_topo = run_g2op(env, g2op_actions, topo_act, topo_chron_ids, max_ts)
    light_topo = run_light(env, light_env, topo_act, topo_chron_ids, max_ts)
    env.reactivate_forecast()
    g2op_look = run_lookahead_g2op(env, g2op_actions, candidates, chron_ids, max_ts)
    light_look = run_lookahead_light(env, light_env, candidates, chron_ids, max_ts)
    env.close()

    print_table("Do nothing", [("grid2op + LightSimBackend", g2op_dn), ("LightEnv", light_dn)])
    print_table(f"Topology (a random unitary set_bus action every {act_every} steps, {nb_episode_topo} episodes)",
                [("grid2op + LightSimBackend", g2op_topo), ("LightEnv", light_topo)])
    print_table(f"Lookahead ({len(candidates)} candidate actions per step, one step ahead)",
                [("grid2op obs.simulate", g2op_look), ("LightEnv copy + step", light_look)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark of the light environment against grid2op")
    parser.add_argument("--env_name", type=str, default=ENV_NAME,
                        help="grid2op environment to use (default: l2rpn_case14_sandbox)")
    parser.add_argument("--no_test", type=str2bool, nargs="?", const=True, default=False,
                        help="use the full dataset (downloaded by grid2op) instead of its test chronics")
    parser.add_argument("--max_ts", type=int, default=10_000, help="maximum number of steps per chronic")
    parser.add_argument("--nb_chronics", type=int, default=3, help="number of chronics played")
    parser.add_argument("--act_every", type=int, default=10,
                        help="topology workload: a random action is played every `act_every` steps")
    parser.add_argument("--nb_episode_topo", type=int, default=50,
                        help="topology workload: number of episodes played")
    parser.add_argument("--nb_candidates", type=int, default=10,
                        help="lookahead workload: number of actions evaluated at every step")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    main(env_name=args.env_name, test=not args.no_test, max_ts=args.max_ts, nb_chronics=args.nb_chronics,
         act_every=args.act_every, nb_episode_topo=args.nb_episode_topo,
         nb_candidates=args.nb_candidates, seed=args.seed)
