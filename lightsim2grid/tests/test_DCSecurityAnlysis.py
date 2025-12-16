# Copyright (c) 2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import numpy as np
import scipy
from pandapower.pypower.makeLODF import update_LODF_diag

import grid2op

from lightsim2grid.solver import SolverType
from lightsim2grid import ContingencyAnalysis, LightSimBackend
import warnings
import pdb


class TestDCSecurityAnalysis(unittest.TestCase):
    def _aux_do_reset(self):
        return self.env.reset(seed=0,options={"time serie id": 0})
    
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.obs = self._aux_do_reset()
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()

    def test_can_create(self):
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
    
    def test_clear(self):
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)

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
        sa.change_solver(SolverType.DC)
        
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
        sa.change_solver(SolverType.DC)
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
        sa.change_solver(SolverType.DC)
        sa.add_all_n1_contingencies()
        all_conts = sa.computer.my_defaults()
        assert len(all_conts) == self.env.n_line
        assert len(sa._contingency_order) == self.env.n_line
    
    def test_get_flows_simple(self):
        """test the get_flows method in the most simplest way: ask for all contingencies,
        contingencies are given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
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
        sa.change_solver(SolverType.DC)
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
        sa.change_solver(SolverType.DC)
        sa.add_multiple_contingencies(0, 1, 2)
        res_p, res_a, res_v = sa.get_flows(0, 1)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,1]) <= 1e-6

    def test_get_flows_3(self):
        """test the get_flows method in the most simplest way: not all contingencies (not same order as given), 
        contingencies are given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
        sa.add_multiple_contingencies(0, 1, 2)
        res_p, res_a, res_v = sa.get_flows(0, 2)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,2]) <= 1e-6

    def test_get_flows_4(self):
        """test the get_flows method: don't ask for all contingencies (same order as given), 
        contingencies are NOT given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
        sa.add_multiple_contingencies(0, 2, 1)
        res_p, res_a, res_v = sa.get_flows(0, 2)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,2]) <= 1e-6

    def test_get_flows_5(self):
        """test the get_flows method in the most simplest way: not all contingencies (not same order as given), 
        contingencies are NOT given in the right order"""
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
        sa.add_multiple_contingencies(0, 2, 1)
        res_p, res_a, res_v = sa.get_flows(0, 1)
        assert res_a.shape == (2, self.env.n_line)
        assert abs(res_a[0,0]) <= 1e-6
        assert abs(res_a[1,1]) <= 1e-6

    def test_get_flows_multiple(self):
        """test the get_flows function when multiple contingencies"""
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
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
        sa1.change_solver(SolverType.DC)
        conts = [0, [0, 4], [5, 7], 4]
        sa1.add_multiple_contingencies(*conts)
        obs = self.env.reset()
        sa2 = ContingencyAnalysis(self.env)
        sa2.change_solver(SolverType.DC)
        sa2.add_multiple_contingencies(*conts)
        
        res_p1, res_a1, res_v1 = sa1.get_flows()
        res_p2, res_a2, res_v2 = sa2.get_flows()
        
        # basic check that the right flows are 0.
        assert abs(res_a1[0, 0]) <= 1e-6, f"{abs(res_a1[0, 0])} vs 0."
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
        
        # check that the flows on a line that is never disconnected changed
        assert np.max(res_a1[:,2]) - np.min(res_a1[:,2]) >= 1.
        assert np.max(res_a2[:,2]) - np.min(res_a2[:,2]) >= 1.
        assert np.max(res_p1[:,2]) - np.min(res_p1[:,2]) >= 1.
        assert np.max(res_p2[:,2]) - np.min(res_p2[:,2]) >= 1.
            
        # check that indeed the matrix are different
        assert np.max(np.abs(res_a1 - res_a2)) > 1.
        assert np.max(np.abs(res_p1 - res_p2)) > 1.      
        
        
        params = self.env.parameters
        params.ENV_DC = True
        params.MAX_LINE_STATUS_CHANGED = 2
        params.NO_OVERFLOW_DISCONNECTION = True
        self.obs.change_forecast_parameters(params)  
        for cont_id, cont in enumerate(conts):
            if isinstance(cont, (list, tuple, np.ndarray)):
                act_dict = {"set_line_status": [(l_id, -1) for l_id in cont]}
            else:
                act_dict = {"set_line_status": [(cont, -1)]}
            sim_obs1, reward1, done1, info1 = self.obs.simulate(self.env.action_space(act_dict), time_step=0)
            assert not done1
            assert not info1["is_illegal"]
            sim_obs2, reward2, done2, info2 = obs.simulate(self.env.action_space(act_dict), time_step=0)
            assert not done2
            assert not info2["is_illegal"]
            assert np.abs(sim_obs1.q_or).max() <= 1e-6, "In DC there should not be any reactive !"
            assert np.abs(sim_obs2.q_or).max() <= 1e-6, "In DC there should not be any reactive !"
            assert np.allclose(np.abs(sim_obs1.p_or[cont]), 0.), "there should not be any power flow on disconnected lines"
            assert np.allclose(np.abs(sim_obs2.p_or[cont]), 0.), "there should not be any power flow on disconnected lines"
            assert (np.abs(res_p1[cont_id, :] - sim_obs1.p_or) <= 1e-5).all(), f"{res_p1[cont_id, :] - sim_obs1.p_or}"
            assert (np.abs(res_p2[cont_id, :] - sim_obs2.p_or) <= 1e-5).all(), f"{res_p2[cont_id, :] - sim_obs2.p_or}"
            # below: does not pass because the voltages are not the same. Vm are init with results of AC PF for 
            # the contingency analysis and with flat voltage for obs.simulate
            # assert (np.abs(res_a1[cont_id, :] - sim_obs1.a_or) <= 1e-4).all(), f"{res_a1[cont_id, :] - sim_obs1.a_or}"
            # assert (np.abs(res_a2[cont_id, :] - sim_obs2.a_or) <= 1e-4).all(), f"{res_a2[cont_id, :] - sim_obs2.a_or}"

    def test_compare_lodf(self, act=None):
        obs = self._aux_do_reset()
        if act is not None:
            obs, reward, done, info = self.env.step(act)
            assert not done, "Unable to do the action, powerflow diverges"
            assert not info["is_ambiguous"]
            assert not info["is_illegal"]
            assert not info["is_illegal_reco"]
        
        sa = ContingencyAnalysis(self.env)
        sa.change_solver(SolverType.DC)
        sa.add_all_n1_contingencies()
        
        # compute with security analysis
        res_p1, res_a1, res_v1 = sa.get_flows()
        
        # compute with LODF
        gridmodel = self.env.backend._grid.copy()
        gridmodel.change_solver(SolverType.DC)
        res = gridmodel.dc_pf(1. * self.env.backend._debug_Vdc, 10, 1e-7)   
        lor_p, *_ = gridmodel.get_lineor_res()
        tor_p, *_ = gridmodel.get_trafohv_res()
        res_powerflow = np.concatenate((lor_p, tor_p))
        LODF_mat = 1. * gridmodel.get_lodf()
        mat_flow = np.tile(res_powerflow, LODF_mat.shape[0]).reshape(LODF_mat.shape)
        por_lodf = mat_flow + LODF_mat.T * mat_flow.T
        
        # debug ptdf
        PTDF_ = 1.0 * gridmodel.get_ptdf()
        dcSbus = 1.0 * gridmodel.get_dcSbus().real
        PTDF = 1. * PTDF_
        res_ptdf = np.dot(PTDF, dcSbus * gridmodel.get_sn_mva())
        assert np.abs(res_ptdf  - res_powerflow).max() <= 1e-6, "error for ptdf"

        # debug lodf
        nb = gridmodel.total_bus()  # number of buses
        nbr = len(gridmodel.get_lineor_res()[0]) + len(gridmodel.get_trafohv_res()[0])
        f_ = np.concatenate((1 * gridmodel.get_lines().get_bus_id_side_1(), 1 * gridmodel.get_trafos().get_bus_id_side_1()))
        t_ = np.concatenate((1 * gridmodel.get_lines().get_bus_id_side_2(), 1 * gridmodel.get_trafos().get_bus_id_side_2()))
        Cft = scipy.sparse.csr_matrix((np.r_[np.ones(nbr), -np.ones(nbr)],
                                        (np.r_[f_, t_], np.r_[np.arange(nbr), np.arange(nbr)])),
                                        (nb, nbr))
        H = PTDF * Cft
        h = np.diag(H, 0)
        den = (np.ones((nbr, 1)) * h.T * -1 + 1.)
        with np.errstate(divide='ignore', invalid='ignore'):
            LODF_pypower = (H / den)
        update_LODF_diag(LODF_pypower)  
        isfinite = np.isfinite(LODF_pypower)
        assert np.abs(LODF_pypower[isfinite] - LODF_mat[isfinite]).max() <= 1e-6, (f"problem with lodf computation: "
                                                                                    f"{np.abs(LODF_pypower - LODF_mat).max():.2e}")

        # compute with the reference
        nb_real_line = len(gridmodel.get_lineor_res()[0])
        has_conv = np.any(res_v1 != 0., axis=1)
        for l_id in range(type(self.env).n_line):
            if not has_conv[l_id]:
                # nothing to check in this case
                continue
            gridmodel_tmp = gridmodel.copy()
            if l_id < nb_real_line:
                gridmodel_tmp.deactivate_powerline(l_id)
            else:
                t_id = l_id - nb_real_line
                gridmodel_tmp.deactivate_trafo(t_id)
            res_tmp = gridmodel_tmp.dc_pf(1. * self.env.backend._debug_Vdc, 10, 1e-7)   
            if res_tmp.shape[0] == 0:
                continue
            lor_tmp, *_ = gridmodel_tmp.get_lineor_res()
            tor_tmpp, *_ = gridmodel_tmp.get_trafohv_res()
            powerflow_tmp = np.concatenate((lor_tmp, tor_tmpp))
            assert np.abs(res_p1[l_id] - powerflow_tmp).max() <= 1e-6, f"error for line / trafo {l_id}: {np.abs(res_p1[l_id] - powerflow_tmp)}"
            assert np.abs(por_lodf[l_id] - powerflow_tmp).max() <= 1e-6, f"error for line / trafo {l_id}: {np.abs(por_lodf[l_id] - powerflow_tmp)}"

        assert np.abs(por_lodf[has_conv] - res_p1[has_conv]).max() <= 1e-6
        
    def test_compare_lodf_topo(self):
        self.test_compare_lodf(act=self.env.action_space({"set_bus": {"substations_id": [(1, (1, 2, 1, 2, 1, 2))]}}))
    
    def test_compare_lodf_deact_line(self):
        self.test_compare_lodf(act=self.env.action_space({"set_line_status": [(1, -1)]}))


if __name__ == "__main__":
    unittest.main()
