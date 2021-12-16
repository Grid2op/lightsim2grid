# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings
from lightsim2grid import TimeSerie
import grid2op
from lightsim2grid.lightSimBackend import LightSimBackend
import numpy as np

class TimeSerieTester(unittest.TestCase):
    def test_behaviour(self):

        env_name = "l2rpn_case14_sandbox"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make(env_name, backend=LightSimBackend(), test=True)

        time_series = TimeSerie(env)
        time_series.compute_V(scenario_id=0)
        As = time_series.compute_A()  # will contain the flows, in amps at each step (rows) for each powerline (column)
        Ps = time_series.compute_P()  # will contain the flows, in amps at each step (rows) for each powerline (column)
        
        env.set_id(0)
        env.reset()
        # I got the same voltages as a normal pf
        for it_num in range(100):
            obs, *_ = env.step(env.action_space())
            if np.max(np.abs(As[1 + it_num] - obs.a_or))  > 1e-3:
                raise RuntimeError(f"error at it {it_num} for A")
            if np.max(np.abs(Ps[1 + it_num] - obs.p_or))  > 1e-3:
                raise RuntimeError(f"error at it {it_num} for P")
