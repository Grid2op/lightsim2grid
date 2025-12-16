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
import pdb


class TestPickle(unittest.TestCase):
    def _aux_test_2sides(self, grid1, grid2, method_name, test_results=False):
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in ["id", "name", "sub_1_id", "sub_2_id",
                            "pos_1_topo_vect", "pos_2_topo_vect",
                            "connected_global", "connected_1", "connected_2",
                            "bus_1_id", "bus_2_id", "has_res",
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
                
            
    def _aux_test_1side(self, grid1, grid2, method_name, test_results=False):
        assert len(getattr(grid1, method_name)()) == len(getattr(grid2, method_name)())
        for line1, line2 in zip(getattr(grid1, method_name)(), getattr(grid2, method_name)()):
            for attr_nm in ["id", "name", "sub_id",
                            "pos_topo_vect", 
                            "connected",
                            "bus_id",
                            ]:
                assert getattr(line1, attr_nm) == getattr(line2, attr_nm), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            for attr_nm in []:
                # TODO other attributes
                assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
            
            # TODO add the "results" "dataframe"
            if test_results:
                for attr_nm in ["res_p_mw", "res_q_mvar", "res_theta_deg", "res_v_kv"]:
                    # TODO other attributes
                    assert np.allclose(getattr(line1, attr_nm), getattr(line2, attr_nm)), f"{method_name} error for {attr_nm}: {getattr(line1, attr_nm)} vs {getattr(line2, attr_nm)}"
                
        
    def aux_test_2sides(self, grid1, grid2, test_results=False):
        self._aux_test_2sides(grid1, grid2, "get_lines", test_results)
        self._aux_test_2sides(grid1, grid2, "get_trafos", test_results)
        
    def aux_test_1side(self, grid1, grid2, test_results=False):
        self._aux_test_1side(grid1, grid2, "get_loads", test_results)
        self._aux_test_1side(grid1, grid2, "get_generators", test_results)
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
                self.aux_test_2sides(self.env.backend._grid, backend_1._grid)
                self.aux_test_1side(self.env.backend._grid, backend_1._grid)
                
                nb_bus_total = self.env.n_sub * 2
                max_it = 10
                tol = 1e-8
                # TODO test in case the pickle file is corrupted...

                # test dc_pf
                V_0 = np.ones(nb_bus_total, dtype=complex)
                V_0 = self.env.backend._grid.dc_pf(V_0, max_it, tol)

                V_1 = np.ones(nb_bus_total, dtype=complex)
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
