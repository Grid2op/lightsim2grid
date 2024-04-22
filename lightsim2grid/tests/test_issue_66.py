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
from grid2op.Action import PlayableAction
from grid2op.Rules import AlwaysLegal
import grid2op

class Issue66Tester(unittest.TestCase):
    """issue is still not replicated and these tests pass"""
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("educ_case14_storage", test=True, backend=LightSimBackend(),
                                    action_class=PlayableAction,
                                    gamerules_class=AlwaysLegal)
        param = self.env.parameters
        param.NB_TIMESTEP_COOLDOWN_LINE = 0
        param.NB_TIMESTEP_COOLDOWN_SUB = 0
        self.env.change_parameters(param)
        self.env.change_forecast_parameters(param)
        return super().setUp()
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()
    
    def test_disco_load(self):
        """test i can disconnect a load"""
        obs = self.env.reset()
        act = self.env.action_space({"set_bus": {"loads_id": [(0, -1)]}})
        obs, reward, done, info = self.env.step(act)
        assert done
        # should not raise any RuntimeError
    
    def test_disco_gen(self):
        """test i can disconnect a load"""
        obs = self.env.reset()
        act = self.env.action_space({"set_bus": {"generators_id": [(0, -1)]}})
        obs, reward, done, info = self.env.step(act)
        assert done
        # should not raise any RuntimeError

    def test_change_bus_load(self):
        """test i can disconnect a load"""
        obs = self.env.reset()
        act = self.env.action_space({"set_bus": {"loads_id": [(9, 2)],
                                                 "lines_or_id": [(14, 2)]}})
        obs, reward, done, info = self.env.step(act)
        assert len(info['exception']) == 0
        assert not done
        # should not raise any RuntimeError
        
        # isolate the load
        act = self.env.action_space({"set_bus": {"lines_or_id": [(14, 1)]}})
        obs, reward, done, info = self.env.step(act)
        assert done

    def test_change_bus_gen(self):
        """test i can disconnect a gen"""
        obs = self.env.reset()
        act = self.env.action_space({"set_bus": {"generators_id": [(0, 2)],
                                                 "lines_ex_id": [(0, 2)]}})
        obs, reward, done, info = self.env.step(act)
        assert len(info['exception']) == 0
        assert not done
        # should not raise any RuntimeError
        
        # isolate the generator
        act = self.env.action_space({"set_bus": {"lines_ex_id": [(0, 1)]}})
        obs, reward, done, info = self.env.step(act)
        assert done
        
    def test_disco_storage(self):
        """test i can disconnect a storage unit"""
        obs = self.env.reset()
        act = self.env.action_space({"set_bus": {"storages_id": [(0, -1)]}})
        obs, reward, done, info = self.env.step(act)
        assert len(info['exception']) == 0
        assert not done
        # should not raise any RuntimeError
        
        act = self.env.action_space({"set_storage": [(0, -1)]})
        obs, reward, done, info = self.env.step(act)
        assert len(info['exception']) == 0
        assert not done
        # should not raise any RuntimeError
        
        
if __name__ == "__main__":
    unittest.main()
        