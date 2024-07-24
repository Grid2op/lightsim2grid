# Copyright (c) 2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.
import unittest
import warnings
import copy
import numpy as np

import grid2op
from grid2op.Backend import PandaPowerBackend
from grid2op.Action import CompleteAction
from grid2op.Reward import EpisodeDurationReward

from lightsim2grid import LightSimBackend
from lightsim2grid.rewards import N1ContingencyReward
from lightsim2grid.solver import SolverType


TH_LIM_A_REF = np.array([
        541.0,
        450.0,
        375.0,
        636.0,
        175.0,
        285.0,
        335.0,
        657.0,
        496.0,
        827.0,
        442.0,
        641.0,
        840.0,
        156.0,
        664.0,
        235.0,
        119.0,
        179.0,
        1986.0,
        1572.0,
    ])
 
 
class TestN1ContingencyReward_Base(unittest.TestCase):
    def get_env_nm(self):
        return "educ_case14_storage"
    
    def init_backend(self):
        return LightSimBackend(solver_type=SolverType.SparseLUSingleSlack)
    
    def is_dc(self):
        return False
    
    def threshold_margin(self):
        return 1.
    
    def l_ids(self):
        return None
    
    def setUp(self) -> None:
        import sys
        if sys.platform == 'win32':
            self.skipTest("Not working, see issue https://github.com/BDonnot/lightsim2grid/issues/85")
        reward = N1ContingencyReward(dc=self.is_dc(),
                                     threshold_margin=self.threshold_margin(),
                                     l_ids=self.l_ids())
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(self.get_env_nm(), 
                                    test=True,
                                    backend=self.init_backend(),
                                    reward_class=reward,
                                    action_class=CompleteAction,
                                    _add_to_name=type(self).__name__)
        self.env_params = self.env.parameters
        self.env_params.MAX_LINE_STATUS_CHANGED = 999999
        self.env_params.MAX_SUB_CHANGED = 999999
        self.env_params.NB_TIMESTEP_COOLDOWN_LINE = 0
        self.env_params.NB_TIMESTEP_COOLDOWN_SUB = 0
        self.env_params.ACTIVATE_STORAGE_LOSS = False
        self.env.change_parameters(self.env_params)
        self.forecast_params = copy.deepcopy(self.env_params)
        if self.is_dc():
            self.forecast_params.ENV_DC = True
            self.forecast_params.FORECAST_DC = True
        self.env.change_forecast_parameters(self.forecast_params)
        assert (self.env.get_thermal_limit() == TH_LIM_A_REF).all()
        self.my_ids = self.l_ids()
        if self.my_ids is None:
            self.my_ids = list(range(type(self.env).n_line))
        self.env.observation_space.change_reward(EpisodeDurationReward)
        self.aux_reset_correctly()
        return super().setUp()
    
    def aux_reset_correctly(self):
        self.env.change_parameters(self.env_params)
        self.env.change_forecast_parameters(self.forecast_params)
        obs = self.env.reset(seed=0, options={"time serie id": 0})
        return obs
        
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()
    
    def _aux_test_reward(self, obs, reward):
        unsafe_cont = 0
        if self.is_dc():
            th_lim = obs._thermal_limit
            # A = sqrt(p**2 + q**2) / (sqrt3) * v)
            # A**2 = (p**2 + q**2) / (3. v**2)
            # p**2 = A**2 * 3. * v**2 - q**2
            # and now the units...
            # W**2 = A * 3. * V**2 - VAr**2
            # MW**2 = 1e12 * A**2 * 3. * v**2 - 1e12 * VAr**2
            # MW**2 = 1e6*A**2 * 3. * 1e6*(V)**2 *  - MVAr**2
            # MW**2 = kA**2 * 3. * kV**2 - MVAr**2
            p_square = 3. * (1e-3*th_lim)**2 * (obs.v_or)**2 - (obs.q_or)**2
            p_square[p_square <= 0.] = 0.
            th_lim_p = np.sqrt(p_square) * self.threshold_margin()
        
        # grid_cpy = self.env.backend._grid.copy() # TODO DEBUG WINDOWS
        # obs._obs_env.backend._grid = grid_cpy # TODO DEBUG WINDOWS
        # obs._obs_env.backend._reset_res_pointers() # TODO DEBUG WINDOWS
        # obs._obs_env.backend._grid.tell_solver_need_reset() # TODO DEBUG WINDOWS
        # Vinit = np.ones(grid_cpy.total_bus(), dtype=complex) # TODO DEBUG WINDOWS
        # V = grid_cpy.dc_pf(Vinit, 1, 1e-5) # TODO DEBUG WINDOWS
        # print(f"{V.shape[0] = }, {V[0]:.3f}, {V[1]:.3f}") # TODO DEBUG WINDOWS
        # print(f"Obs env {obs._obs_env.backend._grid.get_solver_type()}") # TODO DEBUG WINDOWS
        # print(f"Obs env {obs._obs_env.backend._grid.get_dc_solver_type()}") # TODO DEBUG WINDOWS
        # print("HERE HERE HERE") # TODO DEBUG WINDOWS
        # print(f"env.storage_loss = {self.env.storage_loss}")
        # print(f"_obs_env.storage_loss = {obs._obs_env.storage_loss}")
        # print(f"_amount_storage = {obs._obs_env._amount_storage}")
        # print(f"_amount_storage_prev = {obs._obs_env._amount_storage_prev}")
        # print(f"storage p (from obs): {obs.storage_power}") # TODO DEBUG WINDOWS
        # print(f"storage_power_target p (from obs): {obs.storage_power_target}") # TODO DEBUG WINDOWS
        # print(f"_obs_env._storage_power : {obs._obs_env._storage_power}") # TODO DEBUG WINDOWS
        # print(f"_obs_env.delta_time_seconds : {obs._obs_env.delta_time_seconds}") # TODO DEBUG WINDOWS
        # print(f"storage p seen from backend of obs_env {obs._obs_env.backend._storage_res[0]}")  # TODO DEBUG WINDOWS
        sim_obs, sim_r, sim_d, sim_i = obs.simulate(self.env.action_space(), time_step=0)
        # print(f"without contingency: {sim_d = }, {sim_i['exception']}") # TODO DEBUG WINDOWS
        # print(f"without contingency: {sim_d = }, {sim_i}")
        assert not sim_d, "base powerflow should not diverge"
        # return   # TODO DEBUG WINDOWS
        # print("test:")
        for l_id in self.my_ids:
            sim_obs, sim_r, sim_d, sim_i = obs.simulate(self.env.action_space({"set_line_status": [(l_id, -1)]}),
                                                        time_step=0)
            if not self.is_dc():
                # print(f"for {l_id}: {sim_d = }, {sim_i['exception']}, {(sim_obs.a_or / obs._thermal_limit).max()}")
                if np.any(sim_obs.a_or > obs._thermal_limit * self.threshold_margin()) or sim_d:
                    unsafe_cont += 1       
            else:
                if np.any(np.abs(sim_obs.p_or) > th_lim_p) or sim_d:
                    unsafe_cont += 1   
        # print(f"{unsafe_cont = }")  # TODO DEBUG WINDOWS
        assert reward == (len(self.my_ids) - unsafe_cont), f"{reward} vs {(len(self.my_ids) - unsafe_cont)}"
        
    def test_do_nothing(self):
        obs = self.aux_reset_correctly()
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done, f'done with {info["exception"]}'
        assert len(info["exception"]) == 0
        self._aux_test_reward(obs, reward)
        
        obs, reward, done, info = self.env.step(self.env.action_space())
        assert not done, f'done with {info["exception"]}'
        assert len(info["exception"]) == 0
        self._aux_test_reward(obs, reward)
        
    def test_disconnected_line(self):
        obs = self.aux_reset_correctly()
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_line_status": [(0, -1)]}))
        assert not done
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_line_status": [(0, +1)]}))
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_line_status": [(4, -1)]}))
        assert not done
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_line_status": [(4, +1)]}))
        self._aux_test_reward(obs, reward)
        
    def test_modif_topo(self):
        obs = self.aux_reset_correctly()
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_bus": {"substations_id": [(1, (1, 2, 1, 2, 1, 2))]}}))
        assert not done
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_bus": {"substations_id": [(3, (1, 2, 1, 2, 1, 2))]}}))
        assert not done
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_bus": {"substations_id": [(1, (1, 1, 1, 1, 1, 1))]}}))
        assert not done
        self._aux_test_reward(obs, reward)
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_bus": {"substations_id": [(3, (1, 1, 1, 1, 1, 1))]}}))
        assert not done
        self._aux_test_reward(obs, reward)
    
    def test_copy(self):
        obs = self.aux_reset_correctly()
        env_cpy = self.env.copy()
        obs_cpy = env_cpy.copy()
        assert self.env._reward_helper.template_reward._backend is not env_cpy._reward_helper.template_reward._backend
        obs, reward, done, info = self.env.step(self.env.action_space(
            {"set_bus": {"substations_id": [(3, (1, 2, 1, 2, 1, 2))]}}))
        obs_cpy, reward_cpy, done_cpy, info_cpy = env_cpy.step(self.env.action_space(
            {"set_line_status": [(0, -1)]}))
        self._aux_test_reward(obs_cpy, reward_cpy)
        self._aux_test_reward(obs, reward)
        assert reward != reward_cpy


class TestN1ContingencyReward_DC(TestN1ContingencyReward_Base):
    def is_dc(self):
        return True
    # def l_ids(self):
    #     return [18]
    
    
class TestN1ContingencyReward_LIDS(TestN1ContingencyReward_Base):
    def l_ids(self):
        return [0, 1, 2, 18, 16]
    
    
class TestN1ContingencyReward_Margins(TestN1ContingencyReward_Base):
    def threshold_margin(self):
        return 0.9
    
    
class TestN1ContingencyReward_PP(TestN1ContingencyReward_Base):
    def init_backend(self):
        return PandaPowerBackend(with_numba=False, lightsim2grid=False)


if __name__ == "__main__":
    unittest.main()
