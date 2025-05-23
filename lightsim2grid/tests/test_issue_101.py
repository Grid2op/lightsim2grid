# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import numpy as np
from lightsim2grid import LightSimBackend

import grid2op
from grid2op.Action import CompleteAction
try:
    from grid2op.Space import DEFAULT_ALLOW_DETACHMENT
    CAN_TEST_ISSUE_101 = True
except ImportError:
    CAN_TEST_ISSUE_101 = False


class Issue101Tester(unittest.TestCase):
    """issue is still not replicated and these tests pass"""
    def setUp(self) -> None:
        if not CAN_TEST_ISSUE_101:
            self.skipTest("grid2op does not allow detachement, please upgrade grid2op")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("educ_case14_storage", 
                                    test=True,
                                    allow_detachment=True,
                                    backend=LightSimBackend(),
                                    action_class=CompleteAction)
        params = self.env.parameters
        params.ENV_DOES_REDISPATCHING = False
        self.env.change_parameters(params)
        self.env.change_forecast_parameters(params)
        self.env.reset(seed=0, options={"time serie id": 0})
        obs = self.env.reset()
        return super().setUp()
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()
    
    def test_disco_gen(self, tol=1e-5):
        """test i can disconnect a gen"""
        gen_id = 0
        obs, reward, done, info = self.env.step(self.env.action_space(
            {
                "set_bus": {"generators_id": [(gen_id, -1)]}
            }
        ))
        assert not done, info["exception"]
        assert obs.gen_detached[gen_id]
        assert np.abs(obs.gen_p[gen_id]) <= tol
        assert np.abs(obs.gen_v[gen_id]) <= tol
        assert np.abs(obs.gen_q[gen_id]) <= tol
        assert np.abs(obs.gen_theta[gen_id]) <= tol
        assert abs(obs.gen_p_detached[gen_id] - 79.8) <= tol
        
        obs1, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs1.gen_detached[gen_id]
        assert np.abs(obs1.gen_p[gen_id]) <= tol
        assert np.abs(obs1.gen_v[gen_id]) <= tol
        assert np.abs(obs1.gen_q[gen_id]) <= tol
        assert np.abs(obs1.gen_theta[gen_id]) <= tol
        assert abs(obs1.gen_p_detached[gen_id] - 80.5) <= tol
        
        obs2, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs2.gen_detached[gen_id]
        assert np.abs(obs2.gen_p[gen_id]) <= tol
        assert np.abs(obs2.gen_v[gen_id]) <= tol
        assert np.abs(obs2.gen_q[gen_id]) <= tol
        assert np.abs(obs2.gen_theta[gen_id]) <= tol
        assert abs(obs2.gen_p_detached[gen_id] - 80.5) <= tol
        
    def test_disco_load(self, tol=1e-5):
        """test i can disconnect a load"""
        load_id = 0
        obs, reward, done, info = self.env.step(self.env.action_space(
            {
                "set_bus": {"loads_id": [(load_id, -1)]}
            }
        ))
        assert not done, info["exception"]
        assert obs.load_detached[load_id]
        assert np.abs(obs.load_p[load_id]) <= tol
        assert np.abs(obs.load_q[load_id]) <= tol
        assert np.abs(obs.load_v[load_id]) <= tol
        assert abs(obs.load_p_detached[load_id] - 21.9) <= tol
        
        obs1, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs1.load_detached[load_id]
        assert np.abs(obs1.load_p[load_id]) <= tol
        assert np.abs(obs1.load_q[load_id]) <= tol
        assert np.abs(obs1.load_v[load_id]) <= tol
        assert abs(obs1.load_p_detached[load_id] - 21.7) <= tol
        
        obs2, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs2.load_detached[load_id]
        assert np.abs(obs2.load_p[load_id]) <= tol
        assert np.abs(obs2.load_q[load_id]) <= tol
        assert np.abs(obs2.load_v[load_id]) <= tol
        assert abs(obs2.load_p_detached[load_id] - 21.6) <= tol
    
    def test_disco_storage(self, tol=1e-5):
        """test i can disconnect a storage unit"""
        sto_id = 0
        obs, reward, done, info = self.env.step(self.env.action_space(
            {
                "set_bus": {"storages_id": [(sto_id, -1)]},
                "set_storage": [(0, 1.)]
            }
        ))
        assert not done, info["exception"]
        assert obs.storage_detached[sto_id]
        assert np.abs(obs.storage_power[sto_id]) <= tol, f"{obs.storage_power[sto_id]} vs 0."
        assert abs(obs.storage_p_detached[sto_id] - 1.) <= tol, f"{obs.storage_p_detached[sto_id]} vs 1."
        
        obs1, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs1.storage_detached[sto_id]
        assert abs(obs1.storage_power[sto_id]) <= tol, f"{obs1.storage_power[sto_id]} vs 0."
        
        obs2, _, done, info = self.env.step(self.env.action_space({}))
        assert not done, info["exception"]
        assert obs2.storage_detached[sto_id]
        assert abs(obs2.storage_power[sto_id]) <= tol, f"{obs2.storage_power[sto_id]} vs 0."
        
        
if __name__ == "__main__":
    unittest.main()
        