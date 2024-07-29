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

from lightsim2grid import ContingencyAnalysis
from lightsim2grid import LightSimBackend
import warnings
import pdb


class TestSecurityAnalysis(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.obs = self.env.reset(seed=0,options={"time serie id": 0})
        
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()

    def test_can_create(self):
        sa = ContingencyAnalysis(self.env)
    
    def test_clear(self):
        sa = ContingencyAnalysis(self.env)

        # add simple contingencies
        sa.add_multiple_contingencies(0, 1, 2, 3)
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 4
        assert len(sa._contingency_order) == 4

        # test that clear is working
        sa.clear()
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 0
        assert len(sa._contingency_order) == 0

    def test_add_single_contingency(self):
        sa = ContingencyAnalysis(self.env)
        with self.assertRaises(RuntimeError):
            sa.add_single_contingency("toto")
        with self.assertRaises(RuntimeError):
            sa.add_single_contingency(-1)
            sa.add_single_contingency("toto")
        with self.assertRaises(RuntimeError):
            sa.add_single_contingency(self.env.n_line)

        sa.add_single_contingency(self.env.name_line[2])
        sa.add_single_contingency(0)
        sa.add_single_contingency(0, 1)
        sa.add_single_contingency(self.env.name_line[2], 1)
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 4
        assert len(sa._contingency_order) == 4

    def test_add_multiple_contingencies(self):
        sa = ContingencyAnalysis(self.env)
        # add simple contingencies
        sa.add_multiple_contingencies(0, 1, 2, 3)
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 4
        assert len(sa._contingency_order) == 4

        # clear everything
        sa.clear()
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 0
        assert len(sa._contingency_order) == 0

        # add n-2 contingencies
        sa.add_multiple_contingencies([0, 4], [5, 7])
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 2
        assert len(sa._contingency_order) == 2
    
        # add mixed of everything
        sa.clear()
        sa.add_multiple_contingencies(0, [0, 4], [5, 7], self.env.name_line[2])
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == 4
        assert len(sa._contingency_order) == 4

    def test_add_all_n1_contingencies(self):
        sa = ContingencyAnalysis(self.env)
        sa.add_all_n1_contingencies()
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == self.env.n_line
        assert len(sa._contingency_order) == self.env.n_line
    
    def test_get_flows_simple(self):
        """test the get_flows method in the most simplest way: ask for all contingencies,
        contingencies are given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 1, 2)
        res_p, res_a, res_v = sa.get_flows()
        assert res_a.shape == (3, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,1]) <= 1e-6
        assert abs(res_a[2,2]) <= 1e-6

    def test_get_flows_1(self):
        """test the get_flows method: ask for all contingencies , 
        contingencies are NOT given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 2, 1)
        res_p, res_a, res_v = sa.get_flows()
        assert res_a.shape == (3, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,2]) <= 1e-6
        assert abs(res_a[2,1]) <= 1e-6

    def test_get_flows_2(self):
        """test the get_flows method: don't ask for all contingencies (same order as given), 
        contingencies are given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 1, 2)
        res_p, res_a, res_v = sa.get_flows(0, 1)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,1]) <= 1e-6

    def test_get_flows_3(self):
        """test the get_flows method in the most simplest way: not all contingencies (not same order as given), 
        contingencies are given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 1, 2)
        res_p, res_a, res_v = sa.get_flows(0, 2)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,2]) <= 1e-6

    def test_get_flows_4(self):
        """test the get_flows method: don't ask for all contingencies (same order as given), 
        contingencies are NOT given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 2, 1)
        res_p, res_a, res_v = sa.get_flows(0, 2)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,2]) <= 1e-6

    def test_get_flows_5(self):
        """test the get_flows method in the most simplest way: not all contingencies (not same order as given), 
        contingencies are NOT given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, 2, 1)
        res_p, res_a, res_v = sa.get_flows(0, 1)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,1]) <= 1e-6

    def test_get_flows_multiple(self):
        """test the get_flows function when multiple contingencies"""
        sa = ContingencyAnalysis(self.env)
        sa.add_multiple_contingencies(0, [0, 4], [5, 7], 4)

        # everything
        res_p, res_a, res_v = sa.get_flows()
        assert res_a.shape == (4, self.env.n_line)
        assert abs(res_a[0, 0]) <= 1e-6
        assert abs(res_a[1, 0]) + abs(res_a[1, 4]) <= 1e-6
        assert abs(res_a[2, 5]) + abs(res_a[2, 7]) <= 1e-6
        assert abs(res_a[3, 4]) <= 1e-6

        # only some
        res_p, res_a, res_v = sa.get_flows([0, 4], [5, 7], 4, 0)
        assert res_a.shape == (4, self.env.n_line)
        assert abs(res_a[3, 0]) <= 1e-6
        assert abs(res_a[0, 0]) + abs(res_a[0, 4]) <= 1e-6
        assert abs(res_a[1, 5]) + abs(res_a[1, 7]) <= 1e-6
        assert abs(res_a[2, 4]) <= 1e-6

        # only some, but not all
        res_p, res_a, res_v = sa.get_flows([5, 7], [0, 4], 4)
        assert res_a.shape == (3, self.env.n_line)
        assert abs(res_a[1, 0]) + abs(res_a[1, 4]) <= 1e-6
        assert abs(res_a[0, 5]) + abs(res_a[0, 7]) <= 1e-6
        assert abs(res_a[2, 4]) <= 1e-6
        
        assert abs(res_p[1, 0]) + abs(res_p[1, 4]) <= 1e-6
        assert abs(res_p[0, 5]) + abs(res_p[0, 7]) <= 1e-6
        assert abs(res_p[2, 4]) <= 1e-6

        # ask for a non simulated contingencies
        with self.assertRaises(RuntimeError):
            res_p, res_a, res_v = sa.get_flows(5)

    def test_change_injection(self):
        """test the capacity of the things to handle different steps"""
        sa1 = ContingencyAnalysis(self.env)
        conts = [0, [0, 4], [5, 7], 4]
        sa1.add_multiple_contingencies(*conts)
        obs = self.env.reset()
        sa2 = ContingencyAnalysis(self.env)
        sa2.add_multiple_contingencies(*conts)

        res_p1, res_a1, res_v1 = sa1.get_flows()
        res_p2, res_a2, res_v2 = sa2.get_flows()
        
        # basic check that the right flows are 0.
        assert abs(res_a1[0, 0]) <= 1e-6
        assert abs(res_a1[1, 0]) + abs(res_a1[1, 4]) <= 1e-6
        assert abs(res_a1[2, 5]) + abs(res_a1[2, 7]) <= 1e-6
        assert abs(res_a1[3, 4]) <= 1e-6
        assert abs(res_a2[0, 0]) <= 1e-6
        assert abs(res_a2[1, 0]) + abs(res_a2[1, 4]) <= 1e-6
        assert abs(res_a2[2, 5]) + abs(res_a2[2, 7]) <= 1e-6
        assert abs(res_a2[3, 4]) <= 1e-6
        assert abs(res_a1[0, 0]) <= 1e-6
        
        assert abs(res_p1[1, 0]) + abs(res_p1[1, 4]) <= 1e-6
        assert abs(res_p1[2, 5]) + abs(res_p1[2, 7]) <= 1e-6
        assert abs(res_p1[3, 4]) <= 1e-6
        assert abs(res_p2[0, 0]) <= 1e-6
        assert abs(res_p2[1, 0]) + abs(res_p2[1, 4]) <= 1e-6
        assert abs(res_p2[2, 5]) + abs(res_p2[2, 7]) <= 1e-6
        assert abs(res_p2[3, 4]) <= 1e-6

        # check that indeed the matrix are different
        assert np.max(np.abs(res_a1 - res_a2)) > 1.
        assert np.max(np.abs(res_p1 - res_p2)) > 1.
        
        params = self.env.parameters
        params.MAX_LINE_STATUS_CHANGED = 2
        params.NO_OVERFLOW_DISCONNECTION = True
        self.obs.change_forecast_parameters(params)
        
        for cont_id, cont in enumerate(conts):
            if isinstance(cont, (list, tuple, np.ndarray)):
                act_dict = {"set_line_status": [(l_id, -1) for l_id in cont]}
            else:
                act_dict = {"set_line_status": [(cont, -1)]}
            sim_obs1 = self.obs.simulate(self.env.action_space(act_dict), time_step=0)[0]
            sim_obs2 = obs.simulate(self.env.action_space(act_dict), time_step=0)[0]
            assert (np.abs(res_p1[cont_id, :] - sim_obs1.p_or) <= 1e-5).all(), f"{res_p1[cont_id, :] - sim_obs1.p_or}"
            assert (np.abs(res_p2[cont_id, :] - sim_obs2.p_or) <= 1e-5).all(), f"{res_p2[cont_id, :] - sim_obs2.p_or}"
            assert (np.abs(res_a1[cont_id, :] - sim_obs1.a_or) <= 1e-4).all(), f"{res_a1[cont_id, :] - sim_obs1.a_or}"
            assert (np.abs(res_a2[cont_id, :] - sim_obs2.a_or) <= 1e-4).all(), f"{res_a2[cont_id, :] - sim_obs2.a_or}"


if __name__ == "__main__":
    unittest.main()
