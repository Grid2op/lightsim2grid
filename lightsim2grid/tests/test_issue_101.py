# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings
from lightsim2grid import LightSimBackend

import grid2op
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
            self.env = grid2op.make("l2rpn_case14_sandbox", 
                                     test=True,
                                     allow_detachment=True,
                                     backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        obs = self.env.reset()
        return super().setUp()
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()
    
    def test_disco_gen(self):
        """test i can disconnect a load"""
        obs, reward, done, info = self.env.step(self.env.action_space(
            {
                "set_bus": {"generators_id": [(0, -1)]}
            }
        ))
        assert not done
        _, _, done, info = self.env.step(self.env.action_space({}))
        assert not done
        _, _, done, info = self.env.step(self.env.action_space({}))
        assert not done
    
        
        
if __name__ == "__main__":
    unittest.main()
        