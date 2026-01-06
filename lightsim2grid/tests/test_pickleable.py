# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import tempfile
import pickle
import os
import unittest
import warnings

import numpy as np

from grid2op import make
from lightsim2grid.lightSimBackend import LightSimBackend
from lightsim2grid.gridmodel.compare_gridmodel import compare_gridmodel_input
import pdb


class TestPickle(unittest.TestCase):
    def _aux_test_2sides(self, grid1, grid2, method_name, test_results=False):
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in ["id", "name", "sub1_id", "sub2_id",
                            "pos1_topo_vect", "pos2_topo_vect",
                            "connected_global", "connected1", "connected2",
                            "bus1_id", "bus2_id", "has_res",
                            ]:
                assert getattr(line1, attr_nm) == getattr(line2, attr_nm), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            for attr_nm in ["r_pu", "x_pu", "h1_pu", "h2_pu",
                            "yac_11", "yac_12", "yac_21", "yac_22",
                            "ydc_11", "ydc_12", "ydc_21", "ydc_22",
                            ]:
                # TODO other attributes
                assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            # TODO add the "results" "dataframe"
            if test_results:
                for attr_nm in ["res_p1_mw", "res_q1_mvar", "res_theta1_deg", "res_v1_kv", "res_a1_ka",
                                "res_p2_mw", "res_q2_mvar", "res_theta2_deg", "res_v2_kv", "res_a2_ka"]:
                    # TODO other attributes
                    assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
                
            
    def _aux_test_1side(self, grid1, grid2, method_name, test_results=False,
                        add_attr_int=None,
                        add_attr_float=None):
        li_attr_to_test_int = [
            "id", "name", "sub_id",
            "pos_topo_vect", 
            "connected",
            "bus_id",
        ]
        if add_attr_int is not None:
            li_attr_to_test_int += add_attr_int
            
        li_attr_to_test_float = []
        if add_attr_float is not None:
            li_attr_to_test_float += add_attr_float
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in li_attr_to_test_int:
                assert getattr(line1, attr_nm) == getattr(line2, attr_nm), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            for attr_nm in li_attr_to_test_float:
                # TODO other attributes
                assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            # TODO add the "results" "dataframe"
            if not test_results:
                continue
            li_attr_to_test = ["res_p_mw", "res_q_mvar", "res_theta_deg", "res_v_kv"]
            for attr_nm in li_attr_to_test:
                # TODO other attributes
                if not np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)):
                    diff_ = np.abs(getattr(line1, attr_nm) - getattr(line2, attr_nm))
                    
                    raise AssertionError(f"{method_name} error for {attr_nm} for gen id {line1.id} ({line1.name}): {getattr(line1, attr_nm)} "
                                            f"vs {getattr(line2, attr_nm)} -> {diff_}")
                
        
    def aux_test_2sides(self, grid1, grid2, test_results=False):
        self._aux_test_2sides(grid1, grid2, "get_lines", test_results)
        self._aux_test_2sides(grid1, grid2, "get_trafos", test_results)
        
    def aux_test_1side(self, grid1, grid2, test_results=False):
        self._aux_test_1side(grid1, grid2, "get_loads", test_results)
        self._aux_test_1side(grid1, grid2, "get_generators", test_results,
                             add_attr_int=["is_slack", "voltage_regulator_on"],
                             add_attr_float=["slack_weight", "target_vm_pu", "min_q_mvar", "max_q_mvar"])
        self._aux_test_1side(grid1, grid2, "get_storages", test_results)
        self._aux_test_1side(grid1, grid2, "get_shunts", test_results)
                
    def test_save_load(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self.env = make("l2rpn_idf_2023", test=True, backend=LightSimBackend())
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(os.path.join(tmpdir, "test_pickle.pickle"), "wb") as f:
                    pickle.dump(self.env.backend, f)
                with open(os.path.join(tmpdir, "test_pickle.pickle"), "rb") as f:
                    backend_1 = pickle.load(f)
                assert backend_1._grid.get_solver_type() ==  self.env.backend._grid.get_solver_type()
                assert backend_1._grid.get_dc_solver_type() ==  self.env.backend._grid.get_dc_solver_type()
                
                self.aux_test_2sides(self.env.backend._grid, backend_1._grid)
                self.aux_test_1side(self.env.backend._grid, backend_1._grid)
                tmp = compare_gridmodel_input(self.env.backend._grid, backend_1._grid)
                assert len(tmp) == 0
                
                nb_bus_total = self.env.n_sub * 2
                max_it = 10
                tol = 1e-8
                # TODO test in case the pickle file is corrupted...

                # test dc_pf
                V_0 = np.ones(nb_bus_total, dtype=complex)
                V_1 = V_0.copy()
                V_0 = self.env.backend._grid.dc_pf(V_0, max_it, tol)
                V_1 = backend_1._grid.dc_pf(V_1, max_it, tol)

                assert np.all(np.abs(V_0 - V_1) <= 1e-7), "dc pf does not lead to same results"
                self.aux_test_2sides(self.env.backend._grid, backend_1._grid, True)
                self.aux_test_1side(self.env.backend._grid, backend_1._grid, True)

                # test ac_pf
                V_0 = self.env.backend._grid.ac_pf(V_0, max_it, tol)
                V_1 = backend_1._grid.ac_pf(V_1, max_it, tol)
                assert np.all(np.abs(V_0 - V_1) <= 1e-7), "ac pf does not lead to same results"
                self.aux_test_2sides(self.env.backend._grid, backend_1._grid, True)
                self.aux_test_1side(self.env.backend._grid, backend_1._grid, True)


if __name__ == "__main__":
    unittest.main()
