# Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests of the light environment: registering grid2op actions, playing them, and the
cooldowns. The behaviour is compared step by step with a grid2op environment running the
same chronics on the same grid (``LightSimBackend``)."""

import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.lightEnv import (LightEnv,
                                    TopoAction,
                                    Protections,
                                    ElementType,
                                    topo_action_from_grid2op)

COOLDOWN_SUB = 3
COOLDOWN_LINE = 2
NB_TS_RECO = 5


class TestLightEnvActions(unittest.TestCase):
    def setUp(self) -> None:
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True  # the protections are tested on the light env alone
        param.NB_TIMESTEP_COOLDOWN_SUB = COOLDOWN_SUB
        param.NB_TIMESTEP_COOLDOWN_LINE = COOLDOWN_LINE
        param.NB_TIMESTEP_RECONNECTION = NB_TS_RECO
        param.MAX_SUB_CHANGED = 99  # the light env has no such limit
        param.MAX_LINE_STATUS_CHANGED = 99
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox",
                                    test=True,
                                    backend=LightSimBackend(),
                                    param=param)
        self.env.seed(0)
        self.env.set_id(0)
        self.obs = self.env.reset()
        self.cls = type(self.env)
        self.light_env = self.make_light_env()
        self.do_nothing = self.env.action_space({})

    def tearDown(self) -> None:
        self.env.close()

    def make_light_env(self, thermal_limit_ka=None, max_overflow=99, nb_ts=None) -> LightEnv:
        """A light env on the grid of the grid2op env, replaying its chronics (the first
        `nb_ts` steps only if given). By default the protections never trip."""
        env = self.env
        light_env = LightEnv(env.backend._grid)
        data = env.chronics_handler.real_data.data
        if nb_ts is None:
            nb_ts = data.load_p.shape[0]
        light_env.assign_time_series(
            data.load_p[:nb_ts].astype(float),
            data.load_q[:nb_ts].astype(float),
            data.prod_p[:nb_ts].astype(float),
            (data.prod_v[:nb_ts] / env.backend.prod_pu_to_kv).astype(float),  # the light env wants pu
            np.full((nb_ts, env.n_storage), np.nan),
            np.full((nb_ts, env.n_shunt), np.nan),
            np.full((nb_ts, env.n_shunt), np.nan),
            np.full((nb_ts, 0), np.nan),
            np.full((nb_ts, 0), np.nan),
        )
        if thermal_limit_ka is None:
            thermal_limit_ka = np.full(env.n_line, 1e9)
        protections = Protections()
        protections.set_thermal_limit_or(thermal_limit_ka.astype(float))
        protections.set_thermal_limit_ex((thermal_limit_ka *
                                          env.backend.lines_or_pu_to_kv /
                                          env.backend.lines_ex_pu_to_kv).astype(float))
        protections.set_max_line_time_step_overflow(np.full(env.n_line, max_overflow, dtype=np.int32))
        light_env.protections = protections
        light_env.nb_timestep_cooldown_sub = COOLDOWN_SUB
        light_env.nb_timestep_cooldown_line = COOLDOWN_LINE
        light_env.nb_timestep_reconnection = NB_TS_RECO
        return light_env

    def gridmodel_bus(self, sub_id, local_bus) -> int:
        """the bus id of the LSGrid for a grid2op (substation, local busbar)"""
        if local_bus == -1:
            return -1
        return int(sub_id + (local_bus - 1) * self.cls.n_sub)

    def assert_same_state(self, obs, msg=""):
        """the grid of the light env is in the same state as the grid2op observation"""
        grid = self.light_env.grid
        cls = self.cls
        line_status = np.concatenate([grid.get_lines_status(), grid.get_trafo_status()])
        np.testing.assert_array_equal(obs.line_status, line_status, err_msg=f"line status {msg}")
        for load_id in range(cls.n_load):
            self.assertEqual(grid.get_bus_load(load_id),
                             self.gridmodel_bus(cls.load_to_subid[load_id], obs.load_bus[load_id]),
                             f"bus of load {load_id} {msg}")
        for gen_id in range(cls.n_gen):
            self.assertEqual(grid.get_bus_gen(gen_id),
                             self.gridmodel_bus(cls.gen_to_subid[gen_id], obs.gen_bus[gen_id]),
                             f"bus of gen {gen_id} {msg}")
        # flows: the light env is in kA, grid2op in A
        a_or = 1e3 * np.concatenate([grid.get_line_res1()[3], grid.get_trafo_res1()[3]])
        np.testing.assert_allclose(a_or, obs.a_or, rtol=1e-4, atol=1e-2, err_msg=f"a_or {msg}")
        np.testing.assert_array_equal(obs.time_before_cooldown_sub,
                                      self.light_env.time_before_cooldown_sub,
                                      err_msg=f"cooldown sub {msg}")
        np.testing.assert_array_equal(obs.time_before_cooldown_line,
                                      self.light_env.time_before_cooldown_line,
                                      err_msg=f"cooldown line {msg}")

    def step_both(self, action, act_id):
        """one step in grid2op and in the light env, they should agree on everything"""
        obs, reward, done, info = self.env.step(action)
        lobs, lreward, ldone, ltrunc, linfo = self.light_env.step(act_id)
        msg = f"at step {self.light_env.current_step}"
        self.assertFalse(done, f"grid2op game over {msg}")
        self.assertFalse(ldone, f"light env game over {msg}")
        self.assertEqual(info["is_illegal"], linfo["is_illegal"] == "true", f"is_illegal differs {msg}")
        self.assert_same_state(obs, msg)
        return obs, info

    # --- init_actions ---

    def test_init_actions_valid(self):
        act_sub = self.env.action_space({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2)]}})
        act_disco = self.env.action_space({"set_line_status": [(3, -1)]})
        act_reco = TopoAction()
        act_reco.set_line_status(3, 1)
        self.light_env.init_actions([self.do_nothing, act_sub, act_disco, act_reco])
        self.assertEqual(self.light_env.nb_actions, 4)
        actions = self.light_env.get_actions()
        self.assertTrue(all(act.has_been_checked() for act in actions))
        self.assertTrue(actions[0].is_do_nothing())
        self.assertFalse(actions[1].is_do_nothing())
        self.assertEqual(actions[1].nb_set_bus, 2)
        self.assertEqual(actions[1].nb_set_line_status, 0)
        self.assertEqual(actions[2].nb_set_line_status, 1)
        self.assertEqual(actions[3].nb_set_line_status, 1)

    def test_topo_action_from_grid2op(self):
        act = self.env.action_space({"set_bus": {"loads_id": [(0, 2)], "lines_ex_id": [(0, 2)]},
                                     "set_line_status": [(3, -1)]})
        topo = topo_action_from_grid2op(act)
        self.assertFalse(topo.has_been_checked())
        topo.check_validity(self.env.backend._grid)
        self.assertTrue(topo.has_been_checked())
        self.assertEqual(topo.nb_set_bus, 2)
        self.assertEqual(topo.nb_set_line_status, 1)
        self.assertTrue(topo_action_from_grid2op(self.do_nothing).is_do_nothing())
        with self.assertRaises(TypeError):
            topo_action_from_grid2op("not an action")

    def test_init_actions_invalid(self):
        cls = self.cls
        valid = self.env.action_space({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2)]}})
        self.light_env.init_actions([valid])
        self.assertEqual(self.light_env.nb_actions, 1)

        # busbar -2 (grid2op refuses it on assignment, so it is forced in)
        bus_m2 = self.env.action_space({"set_bus": {"loads_id": [(0, 1)]}})
        bus_m2._set_topo_vect[cls.load_pos_topo_vect[0]] = -2
        # busbar 3 on a grid with 2 busbars per substation
        bus_3 = self.env.action_space({"set_bus": {"loads_id": [(0, 1)]}})
        bus_3._set_topo_vect[cls.load_pos_topo_vect[0]] = 3
        # load 18 on a grid with 11 loads
        load_18 = TopoAction()
        load_18.add_element(ElementType.load, 18, 1)
        # negative id
        load_m1 = TopoAction()
        load_m1.add_element(ElementType.load, -1, 1)
        # line 25 on a grid with 20 lines (set_bus and set_line_status)
        line_25_bus = TopoAction()
        line_25_bus.add_element(ElementType.line_or, 25, 1)
        line_25_status = TopoAction()
        line_25_status.set_line_status(25, -1)
        # transformer 5 on a grid with 5 transformers
        trafo_5 = TopoAction()
        trafo_5.add_element(ElementType.trafo_hv, 5, 1)
        # status 2
        status_2 = TopoAction()
        status_2.set_line_status(0, 2)
        # contradiction: line disconnected and its origin connected to a busbar
        ambiguous = TopoAction()
        ambiguous.set_line_status(2, -1)
        ambiguous.add_element(ElementType.line_or, 2, 1)
        # what the light env does not model
        change_bus = self.env.action_space({"change_bus": {"loads_id": [0]}})
        change_status = self.env.action_space({"change_line_status": [0]})
        redisp = self.env.action_space({"redispatch": [(0, 1.)]})

        invalid = {"bus -2": bus_m2, "bus 3": bus_3, "load 18": load_18, "load -1": load_m1,
                   "line 25 (set_bus)": line_25_bus, "line 25 (set_line_status)": line_25_status,
                   "trafo 5": trafo_5, "status 2": status_2, "ambiguous": ambiguous,
                   "change_bus": change_bus, "change_line_status": change_status,
                   "redispatch": redisp, "not an action": 1}
        for name, act in invalid.items():
            with self.subTest(name=name):
                with self.assertRaises(ValueError) as cm:
                    self.light_env.init_actions([valid, act])
                self.assertIn("action 1", str(cm.exception))
                # nothing registered, the previous actions are kept
                self.assertEqual(self.light_env.nb_actions, 1)

        # element types the light env does not support are refused when added
        for el_type in (ElementType.shunt, ElementType.static_gen, ElementType.dc_line_or, ElementType.dc_line_ex):
            with self.subTest(el_type=el_type):
                with self.assertRaises(ValueError):
                    TopoAction().add_element(el_type, 0, 1)

        # an unchecked action cannot be applied nor an out of range id checked directly
        with self.assertRaises(IndexError):
            load_18.check_validity(self.env.backend._grid)

    # --- step ---

    def test_step_without_actions(self):
        self.light_env.reset()
        self.assert_same_state(self.obs, "at reset")
        obs, reward, done, truncated, info = self.light_env.step(0)
        self.assertFalse(done)
        self.assertEqual(info["is_illegal"], "false")
        self.assertEqual(self.light_env.current_step, 1)
        with self.assertRaises(IndexError):
            self.light_env.step(1)
        self.assertEqual(self.light_env.current_step, 1)  # the env is untouched

    def test_step_unknown_action_id(self):
        self.light_env.init_actions([self.do_nothing, self.do_nothing])
        self.light_env.reset()
        for act_id in (2, -1, 42):
            with self.assertRaises(IndexError):
                self.light_env.step(act_id)
        self.assertEqual(self.light_env.current_step, 0)

    def test_set_bus_applied_and_sub_cooldown(self):
        cls = self.cls
        sub_id = 1  # load 0, gen 0, origin of lines 2, 3, 4 and extremity of line 0
        act_sub = self.env.action_space({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2)]}})
        # touches substation 2 without changing anything (load 1 stays on busbar 1):
        # grid2op still counts it as an impact on the substation
        act_other_sub = self.env.action_space({"set_bus": {"loads_id": [(1, 1)]}})
        self.assertEqual(cls.load_to_subid[1], 2)
        self.light_env.init_actions([self.do_nothing, act_sub, act_other_sub])
        self.light_env.reset()
        self.assert_same_state(self.obs, "at reset")

        obs, info = self.step_both(act_sub, 1)
        self.assertFalse(info["is_illegal"])
        self.assertEqual(obs.load_bus[0], 2)
        self.assertEqual(obs.line_or_bus[2], 2)
        self.assertEqual(self.light_env.grid.get_bus_load(0), self.gridmodel_bus(sub_id, 2))
        self.assertEqual(self.light_env.time_before_cooldown_sub[sub_id], COOLDOWN_SUB)
        self.assertTrue((self.light_env.time_before_cooldown_line == 0).all())

        # the substation is in cooldown: the same action is illegal for COOLDOWN_SUB steps
        for k in range(COOLDOWN_SUB):
            obs, info = self.step_both(act_sub, 1)
            self.assertTrue(info["is_illegal"], f"should be illegal, {k + 1} step(s) after the action")
            self.assertEqual(self.light_env.time_before_cooldown_sub[sub_id], COOLDOWN_SUB - 1 - k)
            self.assertEqual(obs.load_bus[0], 2)  # the topology is kept
        # another substation is free
        obs, info = self.step_both(act_other_sub, 2)
        self.assertFalse(info["is_illegal"])
        self.assertEqual(self.light_env.time_before_cooldown_sub[2], COOLDOWN_SUB)
        # and now the first one is legal again
        obs, info = self.step_both(act_sub, 1)
        self.assertFalse(info["is_illegal"])
        self.assertEqual(self.light_env.time_before_cooldown_sub[sub_id], COOLDOWN_SUB)

    def test_line_status_cooldown(self):
        line_id = 3
        disco = self.env.action_space({"set_line_status": [(line_id, -1)]})
        reco = self.env.action_space({"set_line_status": [(line_id, 1)]})
        self.light_env.init_actions([self.do_nothing, disco, reco])
        self.light_env.reset()

        obs, info = self.step_both(disco, 1)
        self.assertFalse(info["is_illegal"])
        self.assertFalse(obs.line_status[line_id])
        self.assertFalse(self.light_env.grid.get_lines_status()[line_id])
        self.assertEqual(self.light_env.time_before_cooldown_line[line_id], COOLDOWN_LINE)
        self.assertTrue((self.light_env.time_before_cooldown_sub == 0).all())  # a status change is not a topology change

        for k in range(COOLDOWN_LINE):
            obs, info = self.step_both(reco, 2)
            self.assertTrue(info["is_illegal"], f"should be illegal, {k + 1} step(s) after the action")
            self.assertFalse(obs.line_status[line_id])
        obs, info = self.step_both(reco, 2)
        self.assertFalse(info["is_illegal"])
        self.assertTrue(obs.line_status[line_id])
        self.assertTrue(self.light_env.grid.get_lines_status()[line_id])
        self.assertEqual(obs.line_or_bus[line_id], 1)  # back on the busbars it left
        self.assertEqual(self.light_env.time_before_cooldown_line[line_id], COOLDOWN_LINE)

    def test_line_status_through_set_bus(self):
        """disconnecting a line with set_bus = -1 on one end, reconnecting it with set_bus on
        both ends: a line status change, not a substation change"""
        line_id = 4
        disco = self.env.action_space({"set_bus": {"lines_ex_id": [(line_id, -1)]}})
        reco = self.env.action_space({"set_bus": {"lines_or_id": [(line_id, 1)], "lines_ex_id": [(line_id, 1)]}})
        self.light_env.init_actions([self.do_nothing, disco, reco])
        self.light_env.reset()

        obs, info = self.step_both(disco, 1)
        self.assertFalse(info["is_illegal"])
        self.assertFalse(obs.line_status[line_id])
        self.assertEqual(self.light_env.time_before_cooldown_line[line_id], COOLDOWN_LINE)
        self.assertTrue((self.light_env.time_before_cooldown_sub == 0).all())
        for _ in range(COOLDOWN_LINE):
            obs, info = self.step_both(reco, 2)
            self.assertTrue(info["is_illegal"])
        obs, info = self.step_both(reco, 2)
        self.assertFalse(info["is_illegal"])
        self.assertTrue(obs.line_status[line_id])
        self.assertTrue((self.light_env.time_before_cooldown_sub == 0).all())

    def test_reconnection_cooldown_after_protection(self):
        """a line disconnected by the protections cannot be reconnected for
        nb_timestep_reconnection steps (light env only: grid2op's protections differ)"""
        line_id = 0  # carries ~0.16 kA in this scenario
        thermal_limit = np.full(self.env.n_line, 1e9)
        thermal_limit[line_id] = 0.1
        light_env = self.make_light_env(thermal_limit_ka=thermal_limit, max_overflow=0)
        reco = self.env.action_space({"set_line_status": [(line_id, 1)]})
        light_env.init_actions([self.do_nothing, reco])
        light_env.reset()
        self.assertTrue(light_env.grid.get_lines_status()[line_id])  # nothing trips at reset

        obs, reward, done, truncated, info = light_env.step(0)
        self.assertFalse(done)
        self.assertFalse(light_env.grid.get_lines_status()[line_id])
        self.assertEqual(list(light_env.protections.lines_disconnected_this_step()), [line_id])
        self.assertEqual(light_env.time_before_cooldown_line[line_id], NB_TS_RECO)
        self.assertEqual(light_env.protections.line_time_step_overflow[line_id], 0)  # not on overflow any more

        for k in range(NB_TS_RECO):
            obs, reward, done, truncated, info = light_env.step(1)
            self.assertFalse(done)
            self.assertEqual(info["is_illegal"], "true", f"should be illegal, {k + 1} step(s) after the disconnection")
            self.assertFalse(light_env.grid.get_lines_status()[line_id])
            self.assertEqual(light_env.time_before_cooldown_line[line_id], NB_TS_RECO - 1 - k)
            self.assertEqual(list(light_env.protections.lines_disconnected_this_step()), [])
        # legal now: the line is reconnected, and tripped again by the protections
        obs, reward, done, truncated, info = light_env.step(1)
        self.assertFalse(done)
        self.assertEqual(info["is_illegal"], "false")
        self.assertEqual(list(light_env.protections.lines_disconnected_this_step()), [line_id])
        self.assertFalse(light_env.grid.get_lines_status()[line_id])
        self.assertEqual(light_env.time_before_cooldown_line[line_id], NB_TS_RECO)

    def test_reset_restores_topology_and_cooldowns(self):
        act_sub = self.env.action_space({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2)]}})
        disco = self.env.action_space({"set_line_status": [(3, -1)]})
        self.light_env.init_actions([self.do_nothing, act_sub, disco])
        obs_reset, info_reset = self.light_env.reset()
        self.light_env.step(1)
        self.light_env.step(2)
        grid = self.light_env.grid
        self.assertEqual(grid.get_bus_load(0), self.gridmodel_bus(1, 2))
        self.assertFalse(grid.get_lines_status()[3])
        self.assertEqual(self.light_env.current_step, 2)

        obs_reset2, info_reset2 = self.light_env.reset()
        grid = self.light_env.grid
        self.assertEqual(grid.get_bus_load(0), self.gridmodel_bus(1, 1))
        self.assertTrue(grid.get_lines_status()[3])
        self.assertEqual(self.light_env.current_step, 0)
        self.assertTrue((self.light_env.time_before_cooldown_sub == 0).all())
        self.assertTrue((self.light_env.time_before_cooldown_line == 0).all())
        np.testing.assert_allclose(obs_reset2, obs_reset)
        # the observation of a reset is the rho of its own powerflow, not a stale one
        np.testing.assert_allclose(obs_reset2, self.light_env.protections.rho)
        self.assertTrue((obs_reset2 > 0).any())
        self.assert_same_state(self.obs, "after the second reset")
        # the actions are still there, and the episode can be replayed
        self.assertEqual(self.light_env.nb_actions, 3)
        obs, info = self.step_both(act_sub, 1)
        self.assertFalse(info["is_illegal"])

    # --- reward and end of episode ---

    def test_reward_is_fraction_survived(self):
        nb_ts = 4
        light_env = self.make_light_env(nb_ts=nb_ts)
        light_env.reset()
        self.assertEqual(light_env.max_step, nb_ts)
        for step in range(1, nb_ts):
            obs, reward, done, truncated, info = light_env.step(0)
            self.assertFalse(done)
            self.assertAlmostEqual(reward, step / nb_ts)
        obs, reward, done, truncated, info = light_env.step(0)
        self.assertTrue(done)
        self.assertFalse(truncated)
        self.assertEqual(reward, 1.)
        self.assertEqual(info["success"], "true")
        self.assertAlmostEqual(float(info["survival_time"]), 1.)

    def test_survival_time_on_failure(self):
        nb_ts = 4
        light_env = self.make_light_env(nb_ts=nb_ts)
        # disconnecting every line leaves no grid: the powerflow of the first step fails
        light_env.init_actions([self.env.action_space({"set_line_status": [(l_id, -1) for l_id in range(self.env.n_line)]})])
        light_env.reset()
        obs, reward, done, truncated, info = light_env.step(0)
        self.assertTrue(done)
        self.assertEqual(reward, 0.)
        self.assertEqual(info["failure"], "true")
        self.assertAlmostEqual(float(info["survival_time"]), 1 / nb_ts)


if __name__ == "__main__":
    unittest.main()
