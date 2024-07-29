# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.


import lightsim2grid
from lightsim2grid.gridmodel import init_from_pandapower
import pandapower as pp
import pandapower.networks as pn
import unittest
import numpy as np
import warnings
import pdb


class TestDCLine(unittest.TestCase):
    def setUp(self) -> None:
        self.net = pn.case14()
        pp.create_dcline(self.net, from_bus=12, to_bus=13, p_mw=1e2, loss_percent=1.2, loss_mw=25, vm_from_pu=1.01, vm_to_pu=1.02)
        return super().setUp()
    
    def test_init(self):
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pandapower(self.net)
        # different convention in pandapower and lightsim for now
        assert model.get_dclines()[0].target_p_or_mw == -100.

    def _aux_test_dc(self):
        # init ls
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pandapower(self.net)
        assert model.get_dclines()[0].target_p_or_mw == -self.net.dcline["p_mw"].values
        
        # run the dc powerflow for the reference
        pp.rundcpp(self.net, lightsim2grid=False, init="flat")
        theta_ref = 1.0 * self.net.res_bus["va_degree"].values
        p_ref = 1.0 * self.net.res_bus["p_mw"].values
        
        # run the dc powerflow for lightsim2grid
        Vinit = np.ones(self.net.bus.shape[0], dtype=np.complex128) * model.get_init_vm_pu()
        V_ls = model.dc_pf(Vinit, 1, 1.)
        theta_ls = np.rad2deg(np.angle(V_ls))
        p_ls = 1.0 * model.get_Sbus_solver() * model.get_sn_mva()
        
        # different convention in pandapower and lightsim for now
        assert np.allclose(model.get_dclines()[0].res_p_or_mw, -self.net.res_dcline["p_from_mw"].values)
        assert np.allclose(model.get_dclines()[0].res_p_ex_mw, -self.net.res_dcline["p_to_mw"].values)
        # check results are the same
        assert np.allclose(theta_ref, theta_ls)
        # assert np.allclose(p_ref, -1.0 * np.real(p_ls))  # pp does not use that
        
        BDC_ref = self.net._ppc["internal"]["Bbus"].todense()
        BDC_ls = model.get_dcYbus_solver().todense()
        assert np.abs(BDC_ls - BDC_ref).max() <= 1e-6

    def _aux_test_ac(self, max_iteration=10, tolerance_mva=1e-08):
        # init ls
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pandapower(self.net)
        assert model.get_dclines()[0].target_p_or_mw == -self.net.dcline["p_mw"].values
        
        # run the dc powerflow for the reference
        pp.runpp(self.net, lightsim2grid=False, init="flat", v_debug=True,
                 max_iteration=max_iteration, tolerance_mva=tolerance_mva)
        theta_ref = 1.0 * self.net.res_bus["va_degree"].values
        p_ref = 1.0 * self.net.res_bus["p_mw"].values
        
        # run the dc powerflow for lightsim2grid (init from V or pp)
        Vinit_pp = np.ones(self.net.bus.shape[0], dtype=np.complex128)
        Vinit_pp *= self.net._ppc["internal"]["Vm_it"][:,0]
        Vinit_pp *= np.exp(1J * self.net._ppc["internal"]["Va_it"][:,0])
        V_ls = model.ac_pf(Vinit_pp, max_iteration, tolerance_mva)
        theta_ls = np.rad2deg(np.angle(V_ls))
        p_ls = 1.0 * model.get_Sbus_solver() * model.get_sn_mva()
        
        # different convention in pandapower and lightsim for now
        assert np.allclose(model.get_dclines()[0].res_p_or_mw, -self.net.res_dcline["p_from_mw"].values)
        assert np.allclose(model.get_dclines()[0].res_p_ex_mw, -self.net.res_dcline["p_to_mw"].values)
        # check results are the same
        assert np.allclose(theta_ref, theta_ls)
        # assert np.allclose(p_ref, -1.0 * np.real(p_ls))  # pp does not use that
        
        Y_ref = self.net._ppc["internal"]["Ybus"].todense()
        Y_ls = model.get_Ybus_solver().todense()
        assert np.abs(Y_ls - Y_ref).max() <= 1e-6
        
    def test_dc_powerflow_without_loss(self):
        self.net.dcline["loss_mw"] = 0.
        self.net.dcline["loss_percent"] = 0.
        self._aux_test_dc()

    def test_dc_powerflow_without_loss_mw(self):
        self.net.dcline["loss_mw"] = 0.
        self._aux_test_dc()

    def test_dc_powerflow_without_loss_pct(self):
        self.net.dcline["loss_percent"] = 0.
        self._aux_test_dc()

    def test_dc_powerflow(self):
        self._aux_test_dc()

    def test_ac_powerflow(self):
        self._aux_test_ac()
        self.net.dcline["p_mw"] = 10.
        self._aux_test_ac()
        self.net.dcline["p_mw"] = 30.
        self._aux_test_ac()
        self.net.dcline["p_mw"] = 300.
        self._aux_test_ac()
        
        
if __name__ == "__main__":
    unittest.main()
