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
import warnings

from lightsim2grid.gridmodel.from_pypowsybl import init

import pandapower.networks as pn
import pandapower as pdp
from lightsim2grid.gridmodel import init as init_from_pp


class AuxInitFromPyPowSyBl:    
    def get_pypo_grid(self):
        return pp.network.create_ieee14()
         
    def get_gen_slack_id(self):
        return 0
    
    def get_slackbus_id(self):
        return 0
    
    def get_equiv_pdp_grid(self):
        return pn.case14()
    
    def get_tol_eq(self):
        return 1e-6
    
    def compare_pp(self):
        """will this test suite compare pypowsybl and pandapower (cannot be used for ieee57 or ieee118)"""
        return True
    
    def setUp(self) -> None:
        self.network_ref = self.get_pypo_grid()
        
        # init equivalent pandapower grid (if any)
        tmp = self.get_equiv_pdp_grid()
        if tmp is not None:
            self.pp_samecase = tmp
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self.ref_samecase = init_from_pp(self.pp_samecase)
            self.can_pp = True
        else:
            # TODO will crash later if no equiv grid
            self.can_pp = False
            self.pp_samecase = None
            self.ref_samecase = None
            
        # init lightsim2grid model
        self.gridmodel = init(self.network_ref, gen_slack_id=self.get_gen_slack_id())
        
        # use some data
        self.nb_bus_total = self.network_ref.get_buses().shape[0]
        self.V_init_dc = np.ones(self.nb_bus_total, dtype=np.complex_)
        self.V_init_ac = 1.04 * self.V_init_dc
        self.tol = 1e-7  # for the solver
        self.tol_eq = self.get_tol_eq()
        return super().setUp()
    
    def test_basic(self):
        """check that all elements are ok"""
        assert len(self.gridmodel.get_lines()) == self.network_ref.get_lines().shape[0]    
        assert len(self.gridmodel.get_loads()) == self.network_ref.get_loads().shape[0]    
        assert len(self.gridmodel.get_generators()) == self.network_ref.get_generators().shape[0]    
        assert len(self.gridmodel.get_trafos()) == self.network_ref.get_2_windings_transformers().shape[0]    
        assert len(self.gridmodel.get_shunts()) == self.network_ref.get_shunt_compensators().shape[0]    
    
    def test_compare_pp(self):
        """compare from the reference case14"""
        if not self.can_pp:
            self.skipTest("no equivalent pandapower grid")
        if not self.compare_pp():
            self.skipTest("It is expected that pandapower and pypowsybl gives different results") 
            
        # check the SBus and Ybus in DC (i need to run powerflow for that)      
        v_ls = self.gridmodel.dc_pf(1.0 * self.V_init_dc, 2, self.tol)
        v_ls_ref = self.ref_samecase.dc_pf(1.0 * self.V_init_dc, 2, self.tol)
        slack_id = self.get_slackbus_id()
        # (array([ 30,  31, 212, 218, 265]), array([265,  31, 212, 218,  30]))
        
        # for case 118
        # bus_or = [64]
        # bus_ex = [67]
        # lines = [el.id for el in self.gridmodel.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # lines_ref = [el.id for el in self.ref_samecase.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # if not lines_ref:
        #     lines_ref = [el.id for el in self.ref_samecase.get_trafos() if (el.bus_hv_id in bus_or and el.bus_lv_id in bus_ex) or (el.bus_hv_id in bus_ex and el.bus_lv_id in bus_or)]
        # import pdb
        # pdb.set_trace()
        # # self.pp_samecase["_ppc"]["internal"]["Bbus"]
        # self.pp_samecase["_ppc"]["internal"]["Ybus"][64,67]
        # self.pp_samecase["_ppc"]["internal"]["bus"]

        # for case 300
        # np.where(tmp_.todense() >= 2000)
        # (array([ 30,  31, 212, 218, 265]), array([265,  31, 212, 218,  30]))

        # bus_or = [64]
        # bus_ex = [67]
        # lines = [el.id for el in self.gridmodel.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # lines_ref = [el.id for el in self.ref_samecase.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # if not lines_ref:
        #     lines_ref = [el.id for el in self.ref_samecase.get_trafos() if (el.bus_hv_id in bus_or and el.bus_lv_id in bus_ex) or (el.bus_hv_id in bus_ex and el.bus_lv_id in bus_or)]
        # # self.pp_samecase["_ppc"]["internal"]["Bbus"]
        # self.pp_samecase["_ppc"]["internal"]["Ybus"][64,67]
        # self.pp_samecase["_ppc"]["internal"]["bus"]
                
        assert np.abs(v_ls - v_ls_ref).max() <= self.tol_eq, "error for vresults for dc"
        tmp_ = self.gridmodel.get_dcYbus() - self.ref_samecase.get_dcYbus()
        assert np.abs(tmp_).max() <= self.tol_eq, "error for dcYbus"
        # check Sbus without slack
        if slack_id > 0:
            assert np.abs(self.gridmodel.get_Sbus()[:slack_id] - self.ref_samecase.get_Sbus()[:slack_id]).max() <= self.tol_eq, "error for dc Sbus"
        if slack_id != self.gridmodel.get_Sbus().shape[0] - 1:
            assert np.abs(self.gridmodel.get_Sbus()[(slack_id+1):] - self.ref_samecase.get_Sbus()[(slack_id+1):]).max() <= self.tol_eq, "error for dc Sbus"

        # same in AC
        v_ls = self.gridmodel.ac_pf(self.V_init_ac, 2, self.tol)
        v_ls_ref = self.ref_samecase.ac_pf(self.V_init_ac, 2, self.tol)
        assert np.abs(self.gridmodel.get_Ybus() - self.ref_samecase.get_Ybus()).max() <= self.tol_eq, "error for Ybus"
        # check Sbus without slack
        if slack_id > 0:
            assert np.abs(self.gridmodel.get_Sbus()[:slack_id] - self.ref_samecase.get_Sbus()[:slack_id]).max() <= self.tol_eq, "error for dc Sbus"
        if slack_id != self.gridmodel.get_Sbus().shape[0] - 1:
            assert np.abs(self.gridmodel.get_Sbus()[(slack_id+1):] - self.ref_samecase.get_Sbus()[(slack_id+1):]).max() <= self.tol_eq, "error for dc Sbus"

    def test_dc_pf(self):
        """test I get the same results as pandapower in dc"""
        v_ls = self.gridmodel.dc_pf(self.V_init_dc, 2, self.tol)
        if self.compare_pp():
            v_ls_ref = self.ref_samecase.dc_pf(self.V_init_dc, 2, self.tol)
            assert np.abs(v_ls - v_ls_ref).max() <= self.tol_eq, "error for vresults for dc"
        lf.run_dc(self.network_ref)
        if self.compare_pp():
            pdp.rundcpp(self.pp_samecase)
        
        # for case 300
        # np.where(np.abs(v_ang_ls - v_ang_pypo) > 3)  # => 173, 177
        # bus_or = [173]
        # bus_ex = [177]
        # # lines = [el.id for el in self.gridmodel.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # lines = [el.id for el in self.gridmodel.get_lines() if (el.bus_or_id in bus_or or el.bus_or_id in bus_ex or el.bus_ex_id in bus_or or el.bus_ex_id in bus_ex)]
        # trafos = [el.id for el in self.gridmodel.get_trafos() if (el.bus_hv_id in bus_or or el.bus_hv_id in bus_ex or el.bus_lv_id in bus_or or el.bus_lv_id in bus_ex)]
        # shunts = [el.id for el in  self.gridmodel.get_shunts() if el.bus_id in bus_or or el.bus_id in bus_ex]
        # trafo 77 => 204    2040    in cdf file
        # trafo 85 => 
        
        # v_mag not really relevant in dc so i study only va
        v_ang_pypo = self.network_ref.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_ls))
        if self.compare_pp():
            v_ang_pp = self.pp_samecase.res_bus["va_degree"].values
            assert np.abs(v_ang_ls - v_ang_pp).max() <= self.tol_eq, "error for va results for dc"
        assert np.abs(v_ang_ls - v_ang_pypo).max() <= self.tol_eq
        
    def test_ac_pf(self):
        # run the powerflows
        v_ls = self.gridmodel.ac_pf(1.0 * self.V_init_ac, 10, self.tol)
        if self.compare_pp():
            v_ls_ref = self.ref_samecase.ac_pf(1.0 * self.V_init_ac, 10, self.tol)   
            assert np.abs(v_ls - v_ls_ref).max() <= self.tol_eq, "error for vresults for ac"
        
        param = lf.Parameters(voltage_init_mode=pp._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
                              transformer_voltage_control_on=False,
                              no_generator_reactive_limits=True,
                              phase_shifter_regulation_on=False,
                              simul_shunt=False,
                              distributed_slack=False,
                              provider_parameters={"slackBusSelectionMode": "NAME",
                                                   "slackBusesIds": self.network_ref.get_buses().iloc[self.get_slackbus_id()].name}
                              ) 
        res_pypow = lf.run_ac(self.network_ref, parameters=param)
        if self.compare_pp():
            pdp.runpp(self.pp_samecase, init="flat")
            
        # check voltage angles
        v_ang_pypo = self.network_ref.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_ls))
        if self.compare_pp():
            v_ang_pp = self.pp_samecase.res_bus["va_degree"].values - self.pp_samecase.ext_grid["va_degree"].values
            assert np.abs(v_ang_ls - v_ang_pp).max() <= self.tol_eq, "error for va results for ac"
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            assert np.abs(v_ang_ls - v_ang_pypo).max() <= self.tol_eq
        
        # check voltage magnitudes
        bus_ref_kv = self.network_ref.get_voltage_levels().loc[self.network_ref.get_buses()["voltage_level_id"].values]["nominal_v"].values
        v_mag_pypo = self.network_ref.get_buses()["v_mag"].values / bus_ref_kv
        v_mag_ls = np.abs(v_ls)
        if self.compare_pp():
            v_mag_pp = self.pp_samecase.res_bus["vm_pu"].values
            assert np.abs(v_mag_ls - v_mag_pp).max() <= self.tol_eq, "error for va results for dc"
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            assert np.abs(v_mag_ls - v_mag_pypo).max() <= self.tol_eq
        
        # check that pypow solution is "feasible" as seen by lightsim2grid
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            v_pypow = v_mag_pypo * np.exp(1j * np.deg2rad(v_ang_pypo))
            v_pypow_ls = self.gridmodel.check_solution(v_pypow, False)
            assert np.abs(v_pypow_ls).max() <= 10. * self.tol_eq
     
        
class TestCase14FromPypo(AuxInitFromPyPowSyBl, unittest.TestCase):
    pass
        
        
class TestCase30FromPypo(AuxInitFromPyPowSyBl, unittest.TestCase):
    # unittest.TestCase does not work because of https://github.com/powsybl/pypowsybl/issues/644
    def get_pypo_grid(self):
        return pp.network.create_ieee30()
    
    def get_equiv_pdp_grid(self):
        return pn.case_ieee30()
        
        
class TestCase57FromPypo(AuxInitFromPyPowSyBl, unittest.TestCase):
    # does not appear to be the same grid !
    def get_pypo_grid(self):
        res = pp.network.create_ieee57()
        # df = res.get_2_windings_transformers()[["rated_u1"]]
        # df["rated_u1"] = 1.0
        # res.update_2_windings_transformers(df)
        
        # df = res.get_lines()[["b1", "b2", "r"]]
        # df["b1"] = 0.
        # df["b2"] = 0.
        # df["r"] = 0.
        # res.update_lines(df)
        return res
    
    def get_equiv_pdp_grid(self):
        return pn.case57()
    
    def compare_pp(self):
        return False
    
    def get_tol_eq(self):
        return 1e-4  # otherwise vangle from pypowsybl and pandapower does not match
                
        
class TestCase118FromPypo(AuxInitFromPyPowSyBl, unittest.TestCase):
    # unittest.TestCase does not work because of https://github.com/powsybl/pypowsybl/issues/644
    def get_pypo_grid(self):
        return pp.network.create_ieee118()
    
    def get_equiv_pdp_grid(self):
        return pn.case118()
     
    def get_gen_slack_id(self):
        return 29
    
    def get_slackbus_id(self):
        return 68
    
    def get_tol_eq(self):
        # return 3e-3  # otherwise vangle from pypowsybl and pandapower does not match
        return 3e-5  # otherwise vangle from pypowsybl and pandapower does not match
    
    def compare_pp(self):
        return False
                
        
class TestCase300FromPypo(AuxInitFromPyPowSyBl):
    # does not work, probably grid not ordered the same
    # need further investigation
    def get_pypo_grid(self):
        return pp.network.create_ieee300()
    
    def get_equiv_pdp_grid(self):
        return pn.case300()
    
    def get_tol_eq(self):
        # return 3e-3  # otherwise vangle from pypowsybl and pandapower does not match
        return 3e-5  # otherwise vangle from pypowsybl and pandapower does not match
    
    def compare_pp(self):
        return False
     
    def get_gen_slack_id(self):
        # does not work with PP, probably bus not ordered the same
        return 55
    
    def get_slackbus_id(self):
        # does not work with PP, probably bus not ordered the same
        return 257


if __name__ == "__main__":
    unittest.main()
    