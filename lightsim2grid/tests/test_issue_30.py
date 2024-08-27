# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import grid2op
import unittest
import numpy as np
from lightsim2grid import LightSimBackend
from grid2op.Action import CompleteAction
import pdb


class TestBug(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox",
                                    test=True,
                                    action_class=CompleteAction,
                                    backend=LightSimBackend())

    def test_diverge(self):
        obs = self.env.reset()

        # prove that it should converged
        sim_o, sim_r, sim_d, sim_i = obs.simulate(self.env.action_space())
        assert not sim_d, "should have converged"

        # put the grid in a terrible spot
        act10 = self.env.action_space({"set_bus": {"lines_or_id": [(2, 1), (3, 2), (4, 1)],
                                                   "lines_ex_id": [(0, 2)],
                                                   "loads_id": [(0, 2)],
                                                   "generators_id": [(0, 1)]}})
        obs, reward, done, info = self.env.step(act10)
        assert not done, "should have converged here"

        # i set the self.V (used to anyway) to a super low value with a grid
        # that barely converges (in the bug it was with a set of actions carefully chosen)
        obs._obs_env.backend.V[:] = np.array([9.66548436e-01-3.55087542e-01j, 9.39831081e-01-4.20738137e-01j,
                                              1.01642368e+00+1.64881592e-01j, 9.92535287e-01+2.97637335e-02j,
                                              7.32930199e-01-3.22013532e-01j, 1.00021220e+00-4.57794282e-01j,
                                              6.07088708e-18-2.53076324e-18j, 1.01531177e+00-4.23251771e-01j,
                                              1.00220145e-01-4.02978083e-02j, 9.17213244e-01-3.85631892e-01j,
                                              9.53039556e-01-4.21909511e-01j, 9.06932626e-01-4.29640029e-01j,
                                              8.31803137e-01-3.93133296e-01j, 3.52206227e-01-1.62093322e-01j,
                                              1.05583333e+00+0.00000000e+00j, 8.67128154e-01-3.45062851e-01j,
                                              8.81047620e-01-5.81841146e-01j, 6.26448055e-01-2.78619470e-01j,
                                              1.05583333e+00+0.00000000e+00j, 1.05583333e+00+0.00000000e+00j,
                                              1.05583333e+00+0.00000000e+00j, 1.05583333e+00+0.00000000e+00j,
                                              1.05583333e+00+0.00000000e+00j, 1.05583333e+00+0.00000000e+00j,
                                              1.05583333e+00+0.00000000e+00j, 1.05583333e+00+0.00000000e+00j,
                                              1.05583333e+00+0.00000000e+00j, 1.05583333e+00+0.00000000e+00j])

        # and it used to make this diverge
        sim_o, sim_r, sim_d, sim_i = obs.simulate(self.env.action_space())
        assert not sim_d, "should have converged"
           
          
if __name__ == "__main__":
    unittest.main()
            