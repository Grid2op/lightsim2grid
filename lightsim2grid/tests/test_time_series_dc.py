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
from lightsim2grid.timeSerie import TimeSerie
import pdb


class TestTSDC_14(unittest.TestCase):
    def make_env(self):
        return grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
    
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = self.make_env()
            self.ts = TimeSerie(self.env)
            
        self.scenario_id = 1
        self.seed = 0
        
    def tearDown(self) -> None:
        self.env.close()
        self.ts.close()
        return super().tearDown()
    
    def test_dc(self):
        V_init = 1.0 * np.abs(self.env.backend.V) + 0j
        
        self.ts.computer.change_solver(SolverType.KLUSingleSlack)
        self.ts.clear()
        res_p, res_a, res_v = self.ts.get_flows(scenario_id=self.scenario_id,
                                                seed=self.seed,
                                                v_init=V_init)
        assert self.ts.computer.get_solver_type() == SolverType.KLUSingleSlack
        self.ts.clear()
        self.ts.computer.change_solver(SolverType.KLUDC)
        res_p_dc, res_a_dc, res_v_dc  = self.ts.get_flows(scenario_id=self.scenario_id,
                                                          seed=self.seed,
                                                          v_init=V_init)
        assert self.ts.computer.get_solver_type() == SolverType.KLUDC
        assert np.any(res_p != res_p_dc)
        assert np.any(res_a != res_a_dc)
        assert np.any(res_v != res_v_dc)
            
        nb_bus = self.env.n_sub
        
        # now check with the DC computation
        all_load = np.ones(self.env.n_load).astype(bool)
        all_gen = np.ones(self.env.n_gen).astype(bool)
        V = 1.0 * V_init
        for ts in range(self.ts.load_p.shape[0]):
            grid_model = self.env.backend._grid.copy()
            grid_model.update_loads_p(all_load, self.ts.load_p[ts, :].astype(np.float32))
            grid_model.update_loads_q(all_load, self.ts.load_q[ts, :].astype(np.float32))
            grid_model.update_gens_p(all_gen, self.ts.prod_p[ts, :].astype(np.float32))
            res = grid_model.dc_pf(V, 10, 1e-8)
            if len(res):
                # model has converged, I check the results are the same
                # check voltages
                assert np.allclose(res_v_dc[ts, :nb_bus], res[:nb_bus]), f"error for step {ts}"
                # now check the flows
                pl_dc, ql_dc, vl_dc, al_dc = grid_model.get_lineor_res()
                pt_dc, qt_dc, vt_dc, at_dc = grid_model.get_trafohv_res()
                # check active power
                p_dc_ref = np.concatenate((pl_dc, pt_dc))
                assert np.allclose(res_p_dc[ts], p_dc_ref), f"error for step {ts}"
                # check amps
                a_dc_ref = np.concatenate((al_dc, at_dc)) * 1000.
                assert np.allclose(res_a_dc[ts], a_dc_ref), f"error for step {ts}"
            else:
                # model has diverged, I check that it has diverge the same way in the security_analysis
                raise AssertionError("The DC computation should have reached the end.")
            V[:] = res


class TestSADC_36(TestTSDC_14):
    def make_env(self):
        return grid2op.make("l2rpn_neurips_2020_track1", test=True, backend=LightSimBackend())


class TestSADC_118(TestTSDC_14):
    def make_env(self):
        return grid2op.make("l2rpn_wcci_2022", test=True, backend=LightSimBackend())
           
          
if __name__ == "__main__":
    unittest.main()

            