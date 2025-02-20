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
