# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import grid2op
from grid2op.Parameters import Parameters
import warnings
import numpy as np
from lightsim2grid import LightSimBackend
from lightsim2grid_cpp import TimeSeriesCPP


class TestTimeSeriesCPP(unittest.TestCase):
    def test_basic(self):
        # print(f"{lightsim2grid_cpp.__file__}")
        env_name = "l2rpn_case14_sandbox"
        # env_name = "l2rpn_neurips_2020_track2_small"
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make(env_name, backend=LightSimBackend(), param=param, test=True)

        nb_bus = env.n_sub
        obs = env.reset()
        grid = env.backend._grid
        Vinit = env.backend.V
        prod_p = 1.0 * env.chronics_handler.real_data.data.prod_p
        load_p = 1.0 * env.chronics_handler.real_data.data.load_p
        load_q = 1.0 * env.chronics_handler.real_data.data.load_q

        # now perform the computation
        computer = TimeSeriesCPP(grid)
        # print("start the computation")
        status = computer.compute_Vs(prod_p,
                                    np.zeros((prod_p.shape[0], 0)),  # no static generators for now !
                                    load_p,
                                    load_q,
                                    Vinit,
                                    env.backend.max_it,
                                    env.backend.tol)
        if status != 1:
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {computer.nb_solved()} step(s)")

        assert computer.nb_solved() == 576, f"Error should have made 576 powerflows, but did {computer.nb_solved()}"
        Vs = computer.get_voltages()

        # check i can call the method to get the buses
        sbuses = computer.get_sbuses()

        # results should be different
        assert np.any(np.abs(np.min(np.abs(Vs), axis=0) - np.max(np.abs(Vs), axis=0)) > 0.01)
        
        # I got the same voltages as a normal pf
        for it_num in range(100):
            obs, *_ = env.step(env.action_space())
            if np.max(np.abs(env.backend.V[:nb_bus] - Vs[1 + it_num, :nb_bus]))  > 1e-6:
                raise RuntimeError(f"error at it {it_num}")

    def test_amps(self):
                # print(f"{lightsim2grid_cpp.__file__}")
        env_name = "l2rpn_case14_sandbox"
        # env_name = "l2rpn_neurips_2020_track2_small"
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make(env_name, backend=LightSimBackend(), param=param, test=True)

        nb_bus = env.n_sub
        obs = env.reset()
        grid = env.backend._grid
        Vinit = env.backend.V
        prod_p = 1.0 * env.chronics_handler.real_data.data.prod_p
        load_p = 1.0 * env.chronics_handler.real_data.data.load_p
        load_q = 1.0 * env.chronics_handler.real_data.data.load_q

        # now perform the computation
        computer = TimeSeriesCPP(grid)
        # print("start the computation")
        status = computer.compute_Vs(prod_p,
                                    np.zeros((prod_p.shape[0], 0)),  # no static generators for now !
                                    load_p,
                                    load_q,
                                    Vinit,
                                    env.backend.max_it,
                                    env.backend.tol)
        if status != 1:
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {computer.nb_solved()} step(s)")

        ampss = computer.compute_flows()

        # I got the same voltages as a normal pf
        for it_num in range(100):
            obs, *_ = env.step(env.action_space())
            if np.max(np.abs(ampss[1 + it_num] - obs.a_or * 1e-3))  > 1e-6:
                raise RuntimeError(f"error at it {it_num}")


if __name__ == "__main__":
    unittest.main()
