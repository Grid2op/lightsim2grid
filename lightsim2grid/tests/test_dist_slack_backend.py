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


class TestDistSlackBackend(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env_ss = grid2op.make("l2rpn_wcci_2022",
                                       backend=LightSimBackend(dist_slack_non_renew=False),
                                       test=True)
            self.env_ds = grid2op.make("l2rpn_wcci_2022",
                                       backend=LightSimBackend(dist_slack_non_renew=True),
                                       test=True)
            
        for env in [self.env_ss, self.env_ds]:
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
        gen_p = np.zeros((self.max_iter_real, env.n_gen))
        while not done:
            obs, reward, done, info = env.step(env.action_space())
            aor[ts,:] = obs.a_or
            gen_p[ts,:] = obs.gen_p
            ts += 1
            if ts >= self.max_iter_real:
                break
        return ts, done, aor, gen_p
        
    def test_different(self):
        self._aux_test_different(self.env_ss, self.env_ds)
        
    def _aux_test_different(self, env_ss, env_ds):
        ts_ss, done_ss, aor_ss, gen_p_ss = self._run_env(env_ss)
        ts_ds, done_ds, aor_ds, gen_p_ds = self._run_env(env_ds)
        
        assert ts_ss == ts_ds
        assert done_ss == done_ds
        
        # non redispatchable gen are not affected
        assert np.all(gen_p_ss[:,~self.env_ds.gen_redispatchable] == gen_p_ds[:,~self.env_ds.gen_redispatchable])
        
        # redispatchable gen are affected (unless p=0)
        for ts in range(self.max_iter_real):
            different = gen_p_ss[ts,self.env_ds.gen_redispatchable] != gen_p_ds[ts,self.env_ds.gen_redispatchable]
            are_zero = gen_p_ss[ts,self.env_ds.gen_redispatchable] == 0.
            assert np.all(different | are_zero), f"error at iteration {ts}"
            # if p ==0. then it's not modified
            assert np.all(gen_p_ds[ts,self.env_ds.gen_redispatchable][are_zero] == 0.), f"error at iteration {ts}"
    
    def test_copy(self):
        """test the properties are maintained after a copy"""
        env_ds = self.env_ds.copy()
        self._prepare_env(env_ds)
        self._aux_test_different(self.env_ss, env_ds)
        
    def test_after_reset(self):
        """test the "reset" does not lose the property"""
        _ = self.env_ds.reset()
        _ = self.env_ds.reset()
        self._prepare_env(self.env_ds)
        self._aux_test_different(self.env_ss, self.env_ds)
    
    def _aux_get_kwargs_runner(self):
        return dict(nb_episode=1,
                    max_iter=self.max_iter_real,
                    add_detailed_output=True,
                    env_seeds=[0])
    
    def test_after_runner(self):
        """test I can use the runner"""
        runner_ss = Runner(**self.env_ss.get_params_for_runner())
        runner_ds = Runner(**self.env_ds.get_params_for_runner())
        res_ss = runner_ss.run(**self._aux_get_kwargs_runner())
        res_ds = runner_ds.run(**self._aux_get_kwargs_runner())
        if res_ss[0][3] != res_ds[0][3]: # same number of steps survived
            raise RuntimeError(f"{res_ss[0][3]} vs {res_ds[0][3]}: ")
        assert res_ss[0][2] != res_ds[0][2]  # not the same reward
        ep_ss = res_ss[0][-1]
        ep_ds = res_ds[0][-1]
        
        # redispatchable gen are affected
        for ts in range(self.max_iter_real):
            gen_p_ss = ep_ss.observations[ts].gen_p
            gen_p_ds = ep_ds.observations[ts].gen_p
            # renewables are not impacted
            assert np.all(gen_p_ss[~self.env_ds.gen_redispatchable] == gen_p_ds[~self.env_ds.gen_redispatchable]), f"error at iteration {ts}"
            # non renewable are impacted (unless p=0.)
            different = gen_p_ss[self.env_ds.gen_redispatchable] != gen_p_ds[self.env_ds.gen_redispatchable]
            are_zero = gen_p_ss[self.env_ds.gen_redispatchable] == 0.
            assert np.all(different | are_zero), f"error at iteration {ts}"
            # if p ==0. then it's not modified
            assert np.all(gen_p_ds[self.env_ds.gen_redispatchable][are_zero] == 0.), f"error at iteration {ts}"
    
        
if __name__ == "__main__":
    unittest.main()
