# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import grid2op
from lightsim2grid import LightSimBackend
try:
    from grid2op.Space import DEFAULT_ALLOW_DETACHMENT
    CAN_TEST_DETACHMENT = True
except ImportError:
    CAN_TEST_DETACHMENT = False
import pdb


class TestDistSlackBackend(unittest.TestCase):
    def setUp(self) -> None:
        if not CAN_TEST_DETACHMENT:
            self.skipTest("Detachement is not available on this grid2op version.")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env0 = grid2op.make("l2rpn_case14_sandbox",
                                       backend=LightSimBackend(automatically_disconnect=False),
                                       allow_detachment=True,
                                       test=True)
            self.env1 = grid2op.make("l2rpn_case14_sandbox",
                                       backend=LightSimBackend(automatically_disconnect=True),
                                       allow_detachment=True,
                                       test=True)
            
        for env in [self.env0, self.env1]:
            self._prepare_env(env)
        self.max_iter_real = 10
        return super().setUp()
    
    def _prepare_env(self, env):
        param = env.parameters
        param.NO_OVERFLOW_DISCONNECTION = True
        param.MAX_SUB_CHANGED = 99999
        param.MAX_LINE_STATUS_CHANGED = 99999
        env.change_parameters(param)
        env.seed(0)
        env.set_id(0)
        env.reset(options={"init ts": 122})
    
    def tearDown(self):
        for env in [self.env0, self.env1]:
            env.close()
        return super().tearDown()

    def test_gen_automatically_disco(self, tol=3e-5):
        act = self.env1.action_space({"set_line_status": [("6_7_18", -1)]})
        # obs1, _, done1, info1 = self.env1.step(self.env1.action_space())
        # assert abs(obs1.gen_p[4] - 10.7) <= 1e-5
            
        obs0, _, done0, info0 = self.env0.step(act)
        obs1, _, done1, info1 = self.env1.step(act)
        assert done0  # standard behaviour: it should be "done"
        assert not done1  # new behaviour: it should not fail
        assert not obs1.line_status[18]  # line 18 is disconnected (the action)
        assert obs1.topo_vect[obs1.gen_pos_topo_vect[4]] == -1  # generator is disconnected
        assert abs(obs1.gen_p_detached[4] - 10.7) <= tol, f"{obs1.gen_p_detached[4]} vs 10.7"
        assert abs(obs1.gen_p_delta[4] + 10.7) <= tol, f"{obs1.gen_p_delta[4]} vs 10.7"
        
        obs1, _, done1, info1 = self.env1.step(act)
        assert not done1  # new behaviour: it should not fail
        assert not obs1.line_status[18]  # line 18 is disconnected (the action)
        assert obs1.topo_vect[obs1.gen_pos_topo_vect[4]] == -1  # generator is disconnected
        assert abs(obs1.gen_p_detached[4] - 11.7) <= tol, f"{obs1.gen_p_detached[4]} vs 11.7"
        assert abs(obs1.gen_p_delta[4] + 0.0) <= tol, f"{obs1.gen_p_delta[4]} vs 0.0"
        
    def test_load_automatically_disco(self, tol=3e-5):
        load_id = 9
        th_val = 11.7
        act = self.env1.action_space({"set_line_status": [("5_12_9", -1),
                                                          ("11_12_13", -1),
                                                          ("12_13_14", -1)]})
        # obs1, _, done1, info1 = self.env1.step(self.env1.action_space())
        # assert abs(obs1.load_p[4] - th_val) <= 1e-5
            
        obs0, _, done0, info0 = self.env0.step(act)
        obs1, _, done1, info1 = self.env1.step(act)
        assert done0  # standard behaviour: it should be "done"
        assert not done1  # new behaviour: it should not fail
        assert not obs1.line_status[9]
        assert not obs1.line_status[13]
        assert not obs1.line_status[14]
        assert obs1.topo_vect[obs1.load_pos_topo_vect[load_id]] == -1  # generator is disconnected
        assert abs(obs1.load_p_detached[load_id] - th_val) <= tol, f"{obs1.load_p_detached[load_id]} vs {th_val}"
        
        obs1, _, done1, info1 = self.env1.step(act)
        assert not done1  # new behaviour: it should not fail
        assert not obs1.line_status[9]
        assert not obs1.line_status[13]
        assert not obs1.line_status[14]
        assert obs1.topo_vect[obs1.load_pos_topo_vect[load_id]] == -1  # generator is disconnected
        assert abs(obs1.load_p_detached[load_id] - th_val) <= tol, f"{obs1.load_p_detached[load_id]} vs {th_val}"


class Test_LSMainComp(unittest.TestCase):
    def setUp(self) -> None:
        if not CAN_TEST_DETACHMENT:
            self.skipTest("Detachement is not available on this grid2op version.")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox",
                                    backend=LightSimBackend(automatically_disconnect=True, dist_slack_non_renew=True),
                                    allow_detachment=True,
                                    test=True)
            
        for env in [self.env]:
            self._prepare_env(env)
        self.max_iter_real = 10
        self.grid = self.env.backend._grid
        return super().setUp()
    
    def _prepare_env(self, env):
        param = env.parameters
        param.NO_OVERFLOW_DISCONNECTION = True
        env.change_parameters(param)
        env.seed(0)
        env.set_id(0)
        env.reset(options={"init ts": 122})
    
    def tearDown(self):
        for env in [self.env]:
            env.close()
        return super().tearDown()
    
    def test_nothing(self):
        self.grid.consider_only_main_component()
        for el in self.grid.get_loads():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_generators():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
            
    def test_load_not_same_bus(self):
        self.grid.change_bus_load(0, 15)
        self.grid.consider_only_main_component()
        for el in self.grid.get_loads():
            if el.id == 0:
                assert not el.connected, f"error for {el.id}"
            else:
                assert el.connected, f"error for {el.id}"
        for el in self.grid.get_generators():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
        # test no crash
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        assert len(info["exception"]) == 0
        assert abs(obs.load_p[0]) <= 1e-5
            
    def test_gen_slack_not_same_bus(self):
        self.grid.change_bus_gen(0, 15)
        self.grid.consider_only_main_component()
        for el in self.grid.get_generators():
            if el.id == 0:
                assert not el.connected, f"error for {el.id}"
            else:
                assert el.connected, f"error for {el.id}"
        for el in self.grid.get_loads():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
        # test no crash
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        assert len(info["exception"]) == 0
        assert abs(obs.gen_p[0]) <= 1e-5
            
    def test_gen_no_slack_not_same_bus(self):
        self.grid.change_bus_gen(3, 15)
        self.grid.consider_only_main_component()
        for el in self.grid.get_generators():
            if el.id == 3:
                assert not el.connected, f"error for {el.id}"
            else:
                assert el.connected, f"error for {el.id}"
        for el in self.grid.get_loads():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
        # test no crash
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        assert len(info["exception"]) == 0
        assert abs(obs.gen_p[3]) <= 1e-5
        
    def test_load_alone(self):
        # load load_12_9
        load_id = 9
        
        # isolate it
        self.grid.change_bus_powerline_ex(9, 12+14)  # 5_12_9
        self.grid.change_bus_powerline_ex(13, 12+14)  # 11_12_13
        self.grid.change_bus_powerline_or(14, 12+14)  # 12_13_14
        
        self.grid.consider_only_main_component()
        for el in self.grid.get_loads():
            if el.id == load_id:
                assert not el.connected, f"error for {el.id}"
            else:
                assert el.connected, f"error for {el.id}"
        for el in self.grid.get_generators():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
        # test no crash
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        assert len(info["exception"]) == 0
        assert abs(obs.load_p[load_id]) <= 1e-5
        
    def test_gen_alone(self):
        # gen gen_2_1
        gen_id = 1
        # isolate it
        self.grid.change_bus_powerline_ex(2, 16)  # 1_2_2
        self.grid.change_bus_powerline_or(5, 16)  # 2_3_5
        self.grid.change_bus_load(1, 16)          # load_2_1
        self.grid.consider_only_main_component()
        for el in self.grid.get_generators():
            if el.id == gen_id:
                assert not el.connected, f"error for {el.id}: {el.connected}"
            else:
                assert el.connected, f"error for {el.id}"
        for el in self.grid.get_loads():
            assert el.connected, f"error for {el.id}"
        for el in self.grid.get_lines():
            assert el.connected_global, f"error for {el.id}"
        for el in self.grid.get_trafos():
            assert el.connected_global, f"error for {el.id}"
        # test no crash
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        assert len(info["exception"]) == 0
        assert abs(obs.gen_p[gen_id]) <= 1e-5

    