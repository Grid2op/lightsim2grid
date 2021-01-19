# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
import unittest
import os
import numpy as np
import warnings
import tempfile
from grid2op.tests.test_Runner import PATH_ADN_CHRONICS_FOLDER, dt_float, PATH_DATA_TEST_PP
from grid2op.tests.test_Runner import Multifolder
from grid2op.tests.test_Runner import HelperTests, L2RPNReward, make
from grid2op.Runner import Runner
from grid2op.Agent import RecoPowerlineAgent

from lightsim2grid.LightSimBackend import LightSimBackend

import timeit

DETAILED_TIMER_INFO = False


def timer(verbose):
    def timer_inner(function):
        def new_function(*args, **kwargs):
            if verbose:
                print('"{name}"'.format(name=function.__name__))
            start_time = timeit.default_timer()
            function(*args, **kwargs)
            elapsed = timeit.default_timer() - start_time
            if verbose:
                print('\t\t"{name}", {time:.3f}s'.format(name=function.__name__, time=elapsed))
        return new_function
    return timer_inner


class TestRunner(HelperTests):
    """Test that i can use all functionalities of the grid2op runner (including execution in multi processing)"""
    def setUp(self):
        self.init_grid_path = os.path.join(PATH_DATA_TEST_PP, "test_case14.json")
        self.path_chron = PATH_ADN_CHRONICS_FOLDER
        self.parameters_path = None
        self.max_iter = 10
        self.real_reward = dt_float(179.99818)
        self.names_chronics_to_backend = {"loads": {"2_C-10.61": 'load_1_0', "3_C151.15": 'load_2_1',
                                                    "14_C63.6": 'load_13_2', "4_C-9.47": 'load_3_3',
                                                    "5_C201.84": 'load_4_4',
                                                    "6_C-6.27": 'load_5_5', "9_C130.49": 'load_8_6',
                                                    "10_C228.66": 'load_9_7',
                                                    "11_C-138.89": 'load_10_8', "12_C-27.88": 'load_11_9',
                                                    "13_C-13.33": 'load_12_10'},
                                          "lines": {'1_2_1': '0_1_0', '1_5_2': '0_4_1', '9_10_16': '8_9_2',
                                                    '9_14_17': '8_13_3',
                                                    '10_11_18': '9_10_4', '12_13_19': '11_12_5', '13_14_20': '12_13_6',
                                                    '2_3_3': '1_2_7', '2_4_4': '1_3_8', '2_5_5': '1_4_9',
                                                    '3_4_6': '2_3_10',
                                                    '4_5_7': '3_4_11', '6_11_11': '5_10_12', '6_12_12': '5_11_13',
                                                    '6_13_13': '5_12_14', '4_7_8': '3_6_15', '4_9_9': '3_8_16',
                                                    '5_6_10': '4_5_17',
                                                    '7_8_14': '6_7_18', '7_9_15': '6_8_19'},
                                          "prods": {"1_G137.1": 'gen_0_4', "3_G36.31": "gen_2_1", "6_G63.29": "gen_5_2",
                                                    "2_G-56.47": "gen_1_0", "8_G40.43": "gen_7_3"},
                                          }
        self.gridStateclass = Multifolder
        self.backendClass = LightSimBackend

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.runner = Runner(init_grid_path=self.init_grid_path,
                                 path_chron=self.path_chron,
                                 parameters_path=self.parameters_path,
                                 names_chronics_to_backend=self.names_chronics_to_backend,
                                 gridStateclass=self.gridStateclass,
                                 backendClass=self.backendClass,
                                 rewardClass=L2RPNReward,
                                 max_iter=self.max_iter,
                                 name_env="test_runner_env")

    @timer(DETAILED_TIMER_INFO)
    def test_one_episode(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            _, cum_reward, timestep, ep_data = self.runner.run_one_episode(max_iter=self.max_iter)
        assert int(timestep) == self.max_iter
        assert np.abs(cum_reward - self.real_reward) <= self.tol_one

    @timer(DETAILED_TIMER_INFO)
    def test_one_process_par(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            res = Runner._one_process_parrallel(self.runner, [0], 0, None, None, self.max_iter)
        assert len(res) == 1
        _, el1, el2, el3, el4 = res[0]
        assert el1 == "1"
        assert np.abs(el2 - self.real_reward) <= self.tol_one
        assert el3 == 10
        assert el4 == 10

    @timer(DETAILED_TIMER_INFO)
    def test_2episode(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            res = self.runner._run_sequential(nb_episode=2, max_iter=self.max_iter)
        assert len(res) == 2
        for i, _, cum_reward, timestep, total_ts in res:
            assert int(timestep) == self.max_iter
            assert np.abs(cum_reward - self.real_reward) <= self.tol_one

    @timer(DETAILED_TIMER_INFO)
    def test_2episode_2process(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            res = self.runner._run_parrallel(nb_episode=2, nb_process=2, max_iter=self.max_iter)
        assert len(res) == 2
        for i, _, cum_reward, timestep, total_ts in res:
            assert int(timestep) == self.max_iter
            assert np.abs(cum_reward - self.real_reward) <= self.tol_one

    @timer(DETAILED_TIMER_INFO)
    def test_complex_agent(self):
        nb_episode = 4
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with make("rte_case5_example", test=True) as env:
                f = tempfile.mkdtemp()
                runner_params = env.get_params_for_runner()
                runner = Runner(**runner_params, agentClass=RecoPowerlineAgent)
                res = runner.run(path_save=f,
                                 nb_episode=nb_episode,
                                 nb_process=2,
                                 max_iter=self.max_iter)
        test_ = set()
        for id_chron, name_chron, cum_reward, nb_time_step, max_ts in res:
            test_.add(name_chron)
        assert len(test_) == nb_episode


if __name__ == "__main__":
    unittest.main()
