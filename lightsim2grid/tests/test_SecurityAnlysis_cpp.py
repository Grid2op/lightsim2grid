# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
from grid2op import Backend
import numpy as np
import grid2op

from lightsim2grid_cpp import SecurityAnalysis
from lightsim2grid import LightSimBackend
import warnings
import pdb

# TODO: not tested when multiple contingencies are run
# TODO: not tested when the topology is changed

class TestSecurityAnalysisCPP(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()

    def test_can_create(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        assert len(SA.my_defaults()) == 0

    def test_add_n1(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_n1(0)
        all_def = SA.my_defaults()
        assert len(all_def) == 1
        assert len(all_def[0]) == 1
        assert all_def[0][0] == 0

        SA.add_n1(1)
        all_def = SA.my_defaults()
        assert len(all_def) == 2
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 1
        assert all_def[0][0] == 0
        assert all_def[1][0] == 1

        # test that if I insert something twice it counts as one:
        SA.add_n1(0)
        all_def = SA.my_defaults()
        assert len(all_def) == 2
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 1
        assert all_def[0][0] == 0
        assert all_def[1][0] == 1

        # test i cannot insert "out of range"
        with self.assertRaises(RuntimeError):
            SA.add_n1(-1)
        with self.assertRaises(RuntimeError):
            SA.add_n1(20)

    def test_add_multiple_n1(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_multiple_n1([0])
        all_def = SA.my_defaults()
        assert len(all_def) == 1
        assert len(all_def[0]) == 1
        assert all_def[0][0] == 0

        SA.add_multiple_n1([1, 2])
        all_def = SA.my_defaults()
        assert len(all_def) == 3
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 1
        assert len(all_def[2]) == 1
        assert all_def[0][0] == 0
        assert all_def[1][0] == 1
        assert all_def[2][0] == 2

        # test that if I insert something twice it counts as one:
        SA.add_multiple_n1([1])
        all_def = SA.my_defaults()
        assert len(all_def) == 3
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 1
        assert len(all_def[2]) == 1
        assert all_def[0][0] == 0
        assert all_def[1][0] == 1
        assert all_def[2][0] == 2

        # test i cannot insert "out of range"
        with self.assertRaises(RuntimeError):
            SA.add_multiple_n1([-1])
        with self.assertRaises(RuntimeError):
            SA.add_multiple_n1([20])

    def test_add_nk(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_nk([0])
        all_def = SA.my_defaults()
        assert len(all_def) == 1
        assert len(all_def[0]) == 1
        assert all_def[0][0] == 0

        SA.add_nk([0, 1])
        all_def = SA.my_defaults()
        assert len(all_def) == 2
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 2
        assert all_def[0][0] == 0
        assert all_def[1][0] == 0
        assert all_def[1][1] == 1

        # test that if I insert something twice it counts as one:
        SA.add_nk([0])
        all_def = SA.my_defaults()
        assert len(all_def) == 2
        assert len(all_def[0]) == 1
        assert len(all_def[1]) == 2
        assert all_def[0][0] == 0
        assert all_def[1][0] == 0
        assert all_def[1][1] == 1

        # test i cannot insert "out of range"
        with self.assertRaises(RuntimeError):
            SA.add_nk([-1])
        with self.assertRaises(RuntimeError):
            SA.add_nk([20])

    def test_add_all_n1(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_all_n1()
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line

    def test_remove_n1(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_all_n1()
        assert SA.remove_n1(0)  # this should remove it and return true (because the removing is a success)
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line - 1 
        assert not SA.remove_n1(0)  # this should not remove another one and return false
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line - 1 

    def test_clear(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_all_n1()
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line
        SA.clear()
        all_def = SA.my_defaults()
        assert len(all_def) == 0

    def test_remove_multiple_n1(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_all_n1()
        nb_removed = SA.remove_multiple_n1([0, 1, 2])
        assert nb_removed == 3
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line - 3
        nb_removed = SA.remove_multiple_n1([2, 3, 4])
        assert nb_removed == 2
        all_def = SA.my_defaults()
        assert len(all_def) == self.env.n_line - 5

    def test_remove_nk(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        SA.add_nk([0, 1])
        SA.add_nk([0, 2])
        SA.add_nk([0, 3])
        all_def = SA.my_defaults()
        assert len(all_def) == 3
        assert SA.remove_nk([0, 1])  # this should remove it and return true (because the removing is a success)
        all_def = SA.my_defaults()
        assert len(all_def) == 2

        assert not SA.remove_nk([0, 1])  # this should not remove another one and return false
        all_def = SA.my_defaults()
        assert len(all_def) == 2

        assert not SA.remove_nk([1, 2])  # this should not remove (never has been in the list)
        all_def = SA.my_defaults()
        assert len(all_def) == 2

    def test_compute(self):
        SA = SecurityAnalysis(self.env.backend._grid)
        lid_cont = [0, 1, 2, 3]
        nb_sub = self.env.n_sub
        SA.add_multiple_n1(lid_cont)
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        res_SA = SA.get_voltages()
        res_flows = SA.compute_flows()
        assert res_SA.shape[0] == len(lid_cont)

        # now check the voltages are correctly computed
        obs = self.env.get_obs()
        for cont_id, l_id in enumerate(lid_cont):
            sim_obs, *_ = obs.simulate(self.env.action_space({"set_line_status": [(l_id, -1)]}),
                                       time_step=0)
            Vref = obs._obs_env.backend.V
            assert np.max(np.abs(Vref[:nb_sub] - res_SA[cont_id, :nb_sub])) <= 1e-6, f"error in V when disconnecting line {l_id} (contingency nb {cont_id})"
            assert np.max(np.abs(res_flows[cont_id] - sim_obs.a_or*1e-3)) <= 1e-6, f"error in flows when disconnecting line {l_id} (contingency nb {cont_id})"
