# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
import unittest
import warnings

import grid2op
from grid2op.tests.test_Runner import TestRunner as TestRunner_glop
from grid2op.tests.test_RunnerFast import TestRunner as TestRunnerFast_glop
from grid2op.tests.test_Runner import HelperTests, L2RPNReward
from grid2op.Runner import Runner

from lightsim2grid.lightSimBackend import LightSimBackend

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


class TestRunner(TestRunner_glop):
    """Test that i can use all functionalities of the grid2op runner (including execution in multi processing)"""
    def setUp(self):
        super().setUp()
        self.backendClass = LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            try:
                # grid2op <= 1.6.5
                self.runner = Runner(init_grid_path=self.init_grid_path,
                                    path_chron=self.path_chron,
                                    parameters_path=self.parameters_path,
                                    names_chronics_to_backend=self.names_chronics_to_backend,
                                    gridStateclass=self.gridStateclass,
                                    backendClass=self.backendClass,
                                    rewardClass=L2RPNReward,
                                    max_iter=self.max_iter,
                                    name_env="test_runner_env")
            except TypeError:
                # grid2op >= 1.6.6
                self.runner = Runner(init_grid_path=self.init_grid_path,
                                    path_chron=self.path_chron,
                                    parameters_path=self.parameters_path,
                                    names_chronics_to_backend=self.names_chronics_to_backend,
                                    gridStateclass=self.gridStateclass,
                                    backendClass=self.backendClass,
                                    rewardClass=L2RPNReward,
                                    max_iter=self.max_iter,
                                    name_env="test_runner_env",
                                    init_env_path=self.init_grid_path)


class TestRunnerFast(TestRunnerFast_glop):
    """Test that i can use all functionalities of the grid2op runner (including execution in multi processing)"""
    def setUp(self):
        super().setUp()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())

        self.runner = Runner(**self.env.get_params_for_runner())


if __name__ == "__main__":
    unittest.main()
