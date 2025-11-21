# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import pdb


import numpy as np
import grid2op
from lightsim2grid import LightSimBackend, SolverType
from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
import pdb


class TestSADC_14(unittest.TestCase):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
            self.sa = ContingencyAnalysis(self.env)
            
    def test_dc(self):
        self.sa.add_all_n1_contingencies()
        res_p, res_a, res_v  = self.sa.get_flows()
        
        self.sa.clear()
        self.sa.computer.change_solver(SolverType.DC)
        self.sa.add_all_n1_contingencies()
        res_p_dc, res_a_dc, res_v_dc  = self.sa.get_flows()
        
        assert np.any(res_p != res_p_dc), "DC and AC solver leads to same results"
        assert np.any(res_a != res_a_dc), "DC and AC solver leads to same results"
        assert np.any(res_v != res_v_dc), "DC and AC solver leads to same results"
        assert self.sa.computer.get_solver_type() == SolverType.DC
        
        nb_bus = self.env.n_sub
        nb_powerline = len(self.env.backend._grid.get_lines())
        # now check with the DC computation
        for l_id in range(type(self.env).n_line):
            grid_model = self.env.backend._grid.copy()
            if l_id < nb_powerline:
                grid_model.deactivate_powerline(l_id)
            else:
                grid_model.deactivate_trafo(l_id - nb_powerline)
            grid_model.tell_solver_need_reset()
            V = 1.0 * self.env.backend.V  # np.ones(2 * self.env.n_sub, dtype=complex)
            res = grid_model.dc_pf(V, 10, 1e-8)
            if len(res):
                # model has converged, I check the results are the same
                # check voltages
                assert np.allclose(res_v_dc[l_id, :nb_bus], res[:nb_bus]), f"error for contingency {l_id}: {np.abs(res_v_dc[l_id, :nb_bus]-res[:nb_bus]).max():.2e}"
                # now check the flows
                pl_dc, ql_dc, vl_dc, al_dc = grid_model.get_lineor_res()
                pt_dc, qt_dc, vt_dc, at_dc = grid_model.get_trafohv_res()
                # check active power
                p_dc_ref = np.concatenate((pl_dc, pt_dc))
                assert np.allclose(res_p_dc[l_id], p_dc_ref), f"error for contingency {l_id}"
                # check amps
                a_dc_ref = np.concatenate((al_dc, at_dc)) * 1000.
                assert np.allclose(res_a_dc[l_id], a_dc_ref), f"error for contingency {l_id}"
            else:
                # model has diverged, I check that it has diverge the same way in the security_analysis
                assert np.all(np.abs(res_v_dc[l_id]) == 0.), f"error for contingency {l_id}"
                assert np.all(np.abs(res_p_dc[l_id]) == 0.), f"error for contingency {l_id}"
                assert np.all(np.isnan(res_a_dc[l_id])), f"error for contingency {l_id}"


class TestSADC_36(TestSADC_14):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_neurips_2020_track1", test=True, backend=LightSimBackend())
            self.sa = ContingencyAnalysis(self.env)


class TestSADC_118(TestSADC_14):
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_wcci_2022", test=True, backend=LightSimBackend())
            self.sa = ContingencyAnalysis(self.env)
           
          
if __name__ == "__main__":
    unittest.main()
            