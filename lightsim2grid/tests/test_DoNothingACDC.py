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
from grid2op.Parameters import Parameters
from lightsim2grid import LightSimBackend

class TestEnvDN(unittest.TestCase):
    """
    This class tests the new functionality introduced in version 0.5.5, where the Ybus matrix
    is not recomputed when the topology does not change. This is an optimization that can
    drastically improve the performances in such cases. 
    """
    def _aux_test_do_nothing(self, is_dc):
        class LightSimBackend_Test(LightSimBackend):
            def __init__(self, detailed_infos_for_cascading_failures=False):
                super().__init__(detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
                self.is_dc = False
            def runpf(self, is_dc=False):
                self.is_dc = is_dc
                return super().runpf(is_dc=is_dc)

        param = Parameters()
        param.ENV_DC = is_dc
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_case14_sandbox",
                               test=True, 
                               param=param,
                               backend=LightSimBackend_Test())
        obs = env.reset()
        # print("\n\n\n")
        # print("##################")
        # print("env.reset() called")
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("1st env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("2nd env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("3rd env.step called")
        assert env.backend.is_dc == is_dc
        obs, *_ = env.step(env.action_space())
        # print("\n\n\n")
        # print("###################")
        # print("4th env.step called")
        assert env.backend.is_dc == is_dc
    
    def test_dc(self):
        self._aux_test_do_nothing(is_dc=True)
    
    def test_ac(self):
        self._aux_test_do_nothing(is_dc=False)
