# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import numpy as np
import pdb


import numpy as np
import grid2op
from grid2op.Runner import Runner
from lightsim2grid import LightSimBackend
import pdb


class TestTurnedOffNoPv(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env_pv = grid2op.make("l2rpn_wcci_2022",
                                       backend=LightSimBackend(turned_off_pv=True),
                                       test=True)
            self.env_npv = grid2op.make("l2rpn_wcci_2022",
                                       backend=LightSimBackend(turned_off_pv=False),
                                       test=True)
            
        for env in [self.env_pv, self.env_npv]:
            self._prepare_env(env)
        self.max_iter_real = 10
        return super().setUp()
    
    def _prepare_env(self, env):
        param = env.parameters
        param.NO_OVERFLOW_DISCONNECTION = True
        env.change_parameters(param)
        env.seed(0)
        env.set_id(0)
    
    def _run_env(self, env):
        obs = env.reset()
        done = False
        ts = 0
        aor = np.zeros((self.max_iter_real, env.n_line))
        gen_v = np.zeros((self.max_iter_real, env.n_gen))
        gen_p = np.zeros((self.max_iter_real, env.n_gen))
        while not done:
            obs, reward, done, info = env.step(env.action_space())
            aor[ts,:] = 1. * obs.a_or
            gen_v[ts,:] = 1. * obs.gen_v
            gen_p[ts,:] = 1. * obs.gen_p
            ts += 1
            if ts >= self.max_iter_real:
                break
        return ts, done, aor, gen_v, gen_p
        
    def test_different(self):
        self._aux_test_different(self.env_pv, self.env_npv)
        
    def _aux_test_different(self, env_pv, env_npv):
        ts_pv, done_pv, aor_pv, gen_v_pv, gen_p_pv = self._run_env(env_pv)
        ts_npv, done_npv, aor_npv, gen_v_npv, gen_p_npv = self._run_env(env_npv)
        assert ts_pv == ts_npv
        assert done_pv == done_npv
        # redispatchable gen are affected
        for ts in range(self.max_iter_real):
            are_zero = gen_p_pv[ts,:] == 0.         
            # non zero should not be modified
            assert np.all(gen_v_pv[ts, ~are_zero] == gen_v_npv[ts, ~are_zero]), f"error at iteration {ts}"
            # at least one p=0 should be modified...
            # all gen with p=0 without other gen such that p=0 at the same bus should be modified
            are_all_zero = ~np.isin(self.env_pv.gen_to_subid[are_zero], self.env_pv.gen_to_subid[~are_zero])
            assert np.all(gen_v_pv[ts, are_zero][are_all_zero] != gen_v_npv[ts, are_zero][are_all_zero]), f"error at iteration {ts}"
            # all gens with p=0 with another gen with p > 0. at the same bus should be the same
            are_some_zero = np.isin(self.env_pv.gen_to_subid[are_zero], self.env_pv.gen_to_subid[~are_zero])
            assert np.all(gen_v_pv[ts, are_zero][are_some_zero] == gen_v_npv[ts, are_zero][are_some_zero]), f"error at iteration {ts}"
    
    def test_copy(self):
        """test the properties are maintained after a copy"""
        env_npv = self.env_npv.copy()
        self._prepare_env(env_npv)
        self._aux_test_different(self.env_pv, env_npv)
        
    def test_after_reset(self):
        """test the "reset" does not lose the property"""
        _ = self.env_npv.reset()
        _ = self.env_npv.reset()
        self._prepare_env(self.env_npv)
        self._aux_test_different(self.env_pv, self.env_npv)
        
    def test_after_runner(self):
        """test I can use the runner"""
        runner_pv = Runner(**self.env_pv.get_params_for_runner())
        runner_npv = Runner(**self.env_npv.get_params_for_runner())
        res_pv = runner_pv.run(nb_episode=1, max_iter=self.max_iter_real, add_detailed_output=True, env_seeds=[0])
        res_npv = runner_npv.run(nb_episode=1, max_iter=self.max_iter_real, add_detailed_output=True, env_seeds=[0])
        assert res_pv[0][3] == res_npv[0][3]  # same number of steps survived
        assert res_pv[0][2] != res_npv[0][2]  # not the same reward
        ep_pv = res_pv[0][-1]
        ep_npv = res_npv[0][-1]
        
        # redispatchable gen are affected
        for ts in range(self.max_iter_real):
            gen_v_pv = ep_pv.observations[ts].gen_v
            gen_v_npv = ep_npv.observations[ts].gen_v
            gen_p_pv = ep_pv.observations[ts].gen_p
            
            are_zero = gen_p_pv == 0.
            # non zero should not be modified
            assert np.all(gen_v_pv[~are_zero] == gen_v_npv[~are_zero]), f"error at iteration {ts}"
            # at least one p=0 should be modified...
            # all gen with p=0 without other gen such that p=0 at the same bus should be modified
            are_all_zero = ~np.isin(ep_pv.observations[ts].gen_to_subid[are_zero], ep_pv.observations[ts].gen_to_subid[~are_zero])
            assert np.all(gen_v_pv[are_zero][are_all_zero] != gen_v_npv[are_zero][are_all_zero]), f"error at iteration {ts}"
            # all gens with p=0 with another gen with p > 0. at the same bus should be the same
            are_some_zero = np.isin(ep_pv.observations[ts].gen_to_subid[are_zero], ep_pv.observations[ts].gen_to_subid[~are_zero])
            assert np.all(gen_v_pv[are_zero][are_some_zero] == gen_v_npv[are_zero][are_some_zero]), f"error at iteration {ts}"
    
        
if __name__ == "__main__":
    unittest.main()
