# Copyright (c) 2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform

import copy
import unittest
import numpy as np
from packaging import version

import pandapower as pp
from pandapower_contingency import run_contingency, run_contingency_ls2g

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.gridmodel.from_pandapower._my_const import _MIN_PP_VERSION_ADV_GRID_MODEL


def run_for_from_bus_loading(net, **kwargs):
    # copied from https://github.com/pawellytaev/pandapower/blob/5abbea2573a85e7fd3e2d21f3ea6c1aee5ffa3cf/pandapower/test/contingency/test_contingency.py#L331
    pp.runpp(net, **kwargs)
    net.res_line["loading_percent"] = net.res_line.i_from_ka / net.line.max_i_ka * 100.
    # 8, 37, 67, 70, 85
    # print(net.res_bus.iloc[37]["vm_pu"])
    if len(net.trafo) > 0:
        max_i_ka_limit = net.trafo.sn_mva.values / (net.trafo.vn_hv_kv.values * np.sqrt(3))
        net.res_trafo["loading_percent"] = net.res_trafo.i_hv_ka / max_i_ka_limit * 100.
        

class GridModelInit(unittest.TestCase):
    def setUp(self) -> None:
        if version.parse(pp.__version__) < _MIN_PP_VERSION_ADV_GRID_MODEL:
            self.skipTest(f"Pandapower too old for this test: requires >={_MIN_PP_VERSION_ADV_GRID_MODEL}, found {pp.__version__}")
        return super().setUp()
    
    def _aux_test_init(self, case=None):
        if case is None:
            case = pp.networks.case118()
            
        gridmodel = init_from_pandapower(case)
        
        # run powerflow
        pp.runpp(case)
        resV = gridmodel.ac_pf(np.ones(case.bus.shape[0], dtype=np.complex128), 10, 1e-8)
        
        # check Ybus
        Ybus_pp = case["_ppc"]["internal"]["Ybus"]
        Ybus_ls = gridmodel.get_Ybus()
        assert Ybus_pp.shape == Ybus_ls.shape
        assert Ybus_ls.size == Ybus_pp.size
        assert (Ybus_ls.nonzero()[0] == Ybus_pp.nonzero()[0]).all()
        assert (Ybus_ls.nonzero()[1] == Ybus_pp.nonzero()[1]).all()
        if np.abs(Ybus_ls - Ybus_pp).max() >= 1e-8:
            tol = 1e-1
            while True:
                row_, col_ = (np.abs(Ybus_ls - Ybus_pp).toarray() >= tol).nonzero()
                if len(row_) >= 1:
                    raise AssertionError("Error for some Ybus coeffs. "
                                         f"There are {(np.abs(Ybus_ls - Ybus_pp) >= 1e-8).sum()} difference (>= 1e-8). "
                                         f"At tolerance {tol:.1e}, issues are with {row_}, {col_}")
                tol /= 10
        # check Sbus
        # TODO
        
        # check resulting voltages
        # TODO
        
        # check flows
        # TODO
        
    def test_case118(self):
        net = pp.networks.case118()   
        self._aux_test_init(net)    
        

class PandapowerContingencyIntegration(unittest.TestCase):
    def setUp(self) -> None:
        if version.parse(pp.__version__) < _MIN_PP_VERSION_ADV_GRID_MODEL:
            self.skipTest(f"Pandapower too old for this test: requires >={_MIN_PP_VERSION_ADV_GRID_MODEL}, found {pp.__version__}")
        return super().setUp()
    
    def _aux_test_case(self, case=None, nminus1_cases=None):
        """test inspired from the test provided in the github issue https://github.com/Grid2Op/lightsim2grid/issues/88#issuecomment-2265150641
        linked https://github.com/pawellytaev/pandapower/blob/da4b5ae6acd42e75dd37ef20d4dcbd823fef48d3/pandapower/test/contingency/test_contingency.py#L156"""
        if case is None:
            case = pp.networks.case118()
        net = copy.deepcopy(case)
        net2 = copy.deepcopy(net)
        if nminus1_cases is None:
            nminus1_cases = {"line": {"index": net.line.index.values},
                             "trafo": {"index": net.trafo.index.values}}

        res = run_contingency(net2, nminus1_cases, contingency_evaluation_function=run_for_from_bus_loading)

        run_contingency_ls2g(net, nminus1_cases)
        for s in ("min", "max"):
            # import pdb
            # pdb.set_trace()
            # (np.abs(res["bus"][f"{s}_vm_pu"] - net.res_bus[f"{s}_vm_pu"].values) >= 1e-5).nonzero()
            assert np.allclose(res["bus"][f"{s}_vm_pu"], net.res_bus[f"{s}_vm_pu"].values, atol=1e-8, rtol=0), f'for bus voltages for {s}, max error {np.abs(res["bus"][f"{s}_vm_pu"] - net.res_bus[f"{s}_vm_pu"].values).max():.2e}'
            assert np.allclose(np.nan_to_num(res["line"][f"{s}_loading_percent"]),
                            net.res_line[f"{s}_loading_percent"].values, atol=1e-6, rtol=0), s
            if len(net.trafo) > 0:
                assert np.allclose(np.nan_to_num(res["trafo"][f"{s}_loading_percent"]),
                                net.res_trafo[f"{s}_loading_percent"].values, atol=1e-6, rtol=0), s
            
    def test_case118(self):
        net = pp.networks.case118()
        cont_lines = net.line.index.values
        cont_trafo = net.trafo.index.values
        
        # not the same behaviour in pp and ls for contingency of these lines / trafo
        lines_diff_behaviour = [6, 7, 103, 121, 163, 164, 170]
        trafos_diff_behaviour = [11, 12]
        #################
        # check that indeed the behaviour is different (ie that at least one bus is disconnected in pandapower)
        for el in lines_diff_behaviour:
            test_net = copy.deepcopy(net)
            test_net.line.loc[el, "in_service"] = False
            pp.runpp(test_net)
            assert (~np.isfinite(test_net.res_bus["vm_pu"])).any()
        for el in trafos_diff_behaviour:
            test_net = copy.deepcopy(net)
            test_net.trafo.loc[el, "in_service"] = False
            pp.runpp(test_net)
            assert (~np.isfinite(test_net.res_bus["vm_pu"])).any()
        # end of the "test of the integration test"
        #############
        # now the contingency analysis
        cont_lines = np.delete(cont_lines, lines_diff_behaviour)
        cont_trafo = np.delete(cont_trafo, trafos_diff_behaviour)
        nminus1_cases = {"line": {"index": cont_lines},
                         "trafo": {"index": cont_trafo}
                        }
        self._aux_test_case(net, nminus1_cases)
                  
    def test_case14(self):
        net = pp.networks.case14()
        nminus1_cases = {"line": {"index": net.line.index.values},
                         "trafo": {"index": np.delete(net.trafo.index.values, -2) # not the same behaviour in pp and ls for contingency of penultimate trafo
                                  }}
        self._aux_test_case(net, nminus1_cases)
  
            
if __name__ == "__main__":
    unittest.main()
            