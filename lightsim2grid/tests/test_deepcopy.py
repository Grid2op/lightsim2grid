# Copyright (c) 2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import grid2op
from grid2op.gym_compat import GymEnv

import numpy as np
from lightsim2grid import LightSimBackend
from lightsim2grid.gridmodel import GridModel
from lightsim2grid.gridmodel.compare_gridmodel import compare_gridmodel_input

import unittest
import warnings
import copy

class TestDistSlackBackend(unittest.TestCase):
    """see issues:
    
    - https://github.com/Grid2op/lightsim2grid/issues/36
    - https://github.com/Grid2op/lightsim2grid/issues/97
    """
    
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, _add_to_name=type(self).__name__,
                                    backend=LightSimBackend())
        return super().setUp()
    def tearDown(self):
        return super().tearDown()
    
    def _aux_test_obs_dict(self, obs, obs_cpy):
        for k in obs_cpy:
            assert k in obs, f"error : key {k} (in the copy) not in orginal obs"
        for k in obs:
            assert k in obs_cpy, f"error : key {k} (in the original obs) not in the copy"
        
        for k in obs:
            assert np.array_equal(obs[k], obs_cpy[k]), f"error for {k}"
            
    def test_deepcopy_glop(self):
        """issue 36: https://github.com/Grid2op/lightsim2grid/issues/36"""
        obs = self.env.reset(seed=0, options={"time serie id": 0})
        env_cpy = copy.deepcopy(self.env)
        new_obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done
        new_obs_2, reward_2, done_2, info_2 = env_cpy.step(self.env.action_space())
        assert not done_2, f'powerflow diverge, error was: \n{info_2["exception"]}'
        assert new_obs == new_obs_2
        
    def test_deepcopy_gym(self):
        """issue 97: https://github.com/Grid2op/lightsim2grid/issues/97"""
        env_gym = GymEnv(self.env)

        # This line raised an error
        env_gym_cpy = copy.deepcopy(env_gym)
        # check the copy is done correctly
        assert isinstance(env_gym_cpy.init_env.backend, LightSimBackend)
        
        assert env_gym_cpy.init_env.backend._grid is not None
        assert isinstance(env_gym_cpy.init_env.backend._grid, GridModel)
        assert env_gym_cpy.init_env.backend._grid is not env_gym.init_env.backend._grid
        tmp = compare_gridmodel_input(
            env_gym_cpy.init_env.backend._grid,
            env_gym.init_env.backend._grid
        )
        assert len(tmp) == 0
        
        assert env_gym_cpy.init_env.backend._LightSimBackend__me_at_init is not None
        assert isinstance(env_gym_cpy.init_env.backend._LightSimBackend__me_at_init, GridModel)
        assert env_gym_cpy.init_env.backend._LightSimBackend__me_at_init is not env_gym.init_env.backend._LightSimBackend__me_at_init
        tmp = compare_gridmodel_input(
            env_gym_cpy.init_env.backend._LightSimBackend__me_at_init,
            env_gym.init_env.backend._LightSimBackend__me_at_init
        )
        assert len(tmp) == 0
        
        # now check they are equal
        obs_cpy, info_cpy = env_gym_cpy.reset(seed=0, options={"time serie id": 0})
        obs, info = env_gym.reset(seed=0, options={"time serie id": 0})
        self._aux_test_obs_dict(obs, obs_cpy)
        
        obs_cpy, reward_cpy, truncated_cpy, done_cpy, info_cpy = env_gym_cpy.step({})
        obs, reward, truncated, done, info = env_gym.step({})
        self._aux_test_obs_dict(obs, obs_cpy)
        
        # assert backends are different
        assert env_gym.init_env is not env_gym_cpy.init_env
        assert env_gym.init_env.backend is not env_gym_cpy.init_env.backend
        
    def test_close_copy_orig_works(self):
        env_gym = GymEnv(self.env)
        env_gym_cpy = copy.deepcopy(env_gym)
        
        obs_cpy, info_cpy = env_gym_cpy.reset()
        obs, info = env_gym.reset()
        self._aux_test_obs_dict(obs, obs_cpy)
        
        obs_cpy, reward_cpy, truncated_cpy, done_cpy, info_cpy = env_gym_cpy.step({})
        obs, reward, truncated, done, info = env_gym.step({})
        self._aux_test_obs_dict(obs, obs_cpy)
        
        # close the copy and make sure the original still works
        env_gym_cpy.close()
        obs, reward, truncated, done, info = env_gym.step({})
        obs, info = env_gym.reset()
        obs, reward, truncated, done, info = env_gym.step({})
        obs, reward, truncated, done, info = env_gym.step({})
        
    def test_close_orig_copy_works(self):
        env_gym = GymEnv(self.env)
        env_gym_cpy = copy.deepcopy(env_gym)
        
        obs_cpy, info_cpy = env_gym_cpy.reset()
        obs, info = env_gym.reset()
        self._aux_test_obs_dict(obs, obs_cpy)
        
        obs_cpy, reward_cpy, truncated_cpy, done_cpy, info_cpy = env_gym_cpy.step({})
        obs, reward, truncated, done, info = env_gym.step({})
        self._aux_test_obs_dict(obs, obs_cpy)
        
        # close the original and make sure the copy still works
        env_gym.close()
        obs, reward, truncated, done, info = env_gym_cpy.step({})
        obs, info = env_gym_cpy.reset()
        obs, reward, truncated, done, info = env_gym_cpy.step({})
        obs, reward, truncated, done, info = env_gym_cpy.step({})
        
        
        

    