# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import pypowsybl as pp
import pypowsybl.loadflow as lf
import numpy as np
import unittest
import tempfile
import os

from lightsim2grid.gridmodel.from_pypowsybl import init

import pandapower.networks as pn
import pandapower.converter as pc
from lightsim2grid.gridmodel import init as init_from_pp


class TestInitFromPyPowSyBl(unittest.TestCase):    
    def setUp(self) -> None:
        self.network = pp.network.create_ieee14()
        # self.network = pp.network.create_ieee118()
        
        with tempfile.TemporaryDirectory() as f:
            tmp_f = os.path.join(f, 'cgmes_ieee118.zip')
            self.network.dump(tmp_f, format='CGMES')
            # if tmp is not None:
            pp_grid = pc.from_cim.from_cim(tmp_f, sn_mva=100.)
            self.pp_gridmo = init_from_pp(pp_grid)
            self.can_pp = True
            
        # init lightsim2grid model
        self.gridmodel = init(self.network)
        
        # use some data
        self.nb_bus_total = self.network.get_buses().shape[0]
        self.V_init = 1.0 * self.network.get_buses()["v_mag"].values
        self.tol = 1e-7
        self.tol_eq = 1e-6
        return super().setUp()
    
    def test_basic(self):
        """check that all elements are ok"""
        assert len(self.gridmodel.get_lines()) == self.network.get_lines().shape[0]    
        assert len(self.gridmodel.get_loads()) == self.network.get_loads().shape[0]    
        assert len(self.gridmodel.get_generators()) == self.network.get_generators().shape[0]    
        assert len(self.gridmodel.get_trafos()) == self.network.get_2_windings_transformers().shape[0]    
        assert len(self.gridmodel.get_shunts()) == self.network.get_shunt_compensators().shape[0]    
    
    def test_compare_pp(self):
        if not self.can_pp:
            self.skipTest("no equivalent pandapower grid")
            
        for l_id, (l_pypo, l_pp) in enumerate(zip(self.gridmodel.get_lines(), self.pp_gridmo.get_lines())):
            assert abs(l_pypo.r_pu - l_pp.r_pu) <= self.tol_eq, f"error for powerline {l_id} r: {l_pypo.r_pu} vs {l_pp.r_pu}"
            assert abs(l_pypo.x_pu - l_pp.x_pu) <= self.tol_eq, f"error for powerline {l_id} x: {l_pypo.x_pu} vs {l_pp.x_pu}"
            assert abs(l_pypo.h_pu - l_pp.h_pu) <= self.tol_eq, f"error for powerline {l_id} h: {l_pypo.h_pu} vs {l_pp.h_pu}"
            # print(f"error for line {l_id} r: {l_pypo.r_pu} vs {l_pp.r_pu}")
            # print(f"error for line {l_id} x: {l_pypo.x_pu} vs {l_pp.x_pu}")
            # print(f"error for line {l_id} h: {l_pypo.h_pu} vs {l_pp.h_pu}")
            
        # for l_id, (l_pypo, l_pp) in enumerate(zip(self.gridmodel.get_trafos(), self.pp_gridmo.get_trafos())):
            # assert abs(l_pypo.r_pu - l_pp.r_pu) <= self.tol_eq, f"error for trafo {l_id} r: {l_pypo.r_pu} vs {l_pp.r_pu}"
            # assert abs(l_pypo.x_pu - l_pp.x_pu) <= self.tol_eq, f"error for trafo {l_id} x: {l_pypo.x_pu} vs {l_pp.x_pu}"
            # assert abs(l_pypo.h_pu - l_pp.h_pu) <= self.tol_eq, f"error for trafo {l_id} h: {l_pypo.h_pu} vs {l_pp.h_pu}"
            # print(f"error for trafo {l_id} r: {l_pypo.r_pu} vs {l_pp.r_pu}")
            # print(f"error for trafo {l_id} x: {l_pypo.x_pu} vs {l_pp.x_pu}")
            # print(f"error for trafo {l_id} h: {l_pypo.h_pu} vs {l_pp.h_pu}")

    def test_dc_pf(self):
        v_ls = self.gridmodel.dc_pf(self.V_init, 2, self.tol)
        lf.run_dc(self.network)
        p_or_ls = [el.res_p_or_mw for el in self.gridmodel.get_lines()]
        p_or_pypo = self.network.get_lines()["p1"].values
        
if __name__ == "__main__":
    unittest.main()
    