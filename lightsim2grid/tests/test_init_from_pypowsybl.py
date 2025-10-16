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

from lightsim2grid.gridmodel import init_from_pypowsybl

try:
    import pandapower.networks as pn
    import pandapower as pdp
    from lightsim2grid.gridmodel import init_from_pandapower
    PDP_AVAIL = True
except ImportError:
    # pandapower not available, eg if testing with numpy 2
    PDP_AVAIL = False

from global_var_tests import MAX_PP_DATAREADER_NOT_BROKEN, CURRENT_PP_VERSION


if CURRENT_PP_VERSION > MAX_PP_DATAREADER_NOT_BROKEN:
    # deactivate compat with recent pandapower version, for now
    PDP_AVAIL = False


class AuxInitFromPyPowSyBlBusesForSub:    
    def use_buses_for_sub(self):
        return True
    
    def get_pypo_grid(self):
        return pp.network.create_ieee14()
    
    def get_slackbus_id(self):
        # id in pandapower which is the same than the ID in the pypowsybl network when loaded from disk
        # but might not be the same as the lightsim2grid gridmodel (if sort_index is True, which is the default)
        return 0
    
    def get_equiv_pdp_grid(self):
        return pn.case14()
    
    def get_tol_eq(self):
        return 1e-6
    
    def compare_pp(self):
        """will this test suite compare pypowsybl and pandapower (cannot be used for ieee57 or ieee118)"""
        return PDP_AVAIL
    
    def setUp(self) -> None:
        self.network_ref = self.get_pypo_grid()
        
        # init equivalent pandapower grid (if any)
        if PDP_AVAIL:
            tmp = self.get_equiv_pdp_grid()
        else:
            tmp = None
            
        if tmp is not None:
            self.pp_samecase = tmp
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                self.ref_samecase = init_from_pandapower(self.pp_samecase)
            self.can_pp = True
        else:
            # TODO will crash later if no equiv grid
            self.can_pp = False
            self.pp_samecase = None
            self.ref_samecase = None
            
        # init lightsim2grid model
        self.gridmodel, self.el_ids = init_from_pypowsybl(self.network_ref,
                                                          slack_bus_id=self.get_slackbus_id(),
                                                          sort_index=False,
                                                          return_sub_id=True,
                                                          buses_for_sub=self.use_buses_for_sub())
        
        # use some data
        if self.use_buses_for_sub():
            self.nb_bus_total = self.network_ref.get_buses().shape[0]
        else:
            self.nb_bus_total = self.gridmodel.get_bus_vn_kv().shape[0]
        self.V_init_dc = np.ones(self.nb_bus_total, dtype=np.complex128)
        self.V_init_ac = 1.06 * self.V_init_dc
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
        v_ls_ref = None
        if self.ref_samecase is not None:
            v_ls_ref = self.ref_samecase.dc_pf(1.0 * self.V_init_dc, 2, self.tol)
        slack_id = self.get_slackbus_id()
        reorder = self.gridmodel._orig_to_ls.reshape(1, -1)
        
        # for case 118
        # reorder_flat = reorder.reshape(-1)
        # bus_or = [64]  # pandapower
        # bus_ex = [67]  # pandapower
        # lines = [el.id for el in self.gridmodel.get_lines() if (reorder_flat[el.bus_or_id] in bus_or and reorder_flat[el.bus_ex_id] in bus_ex) or (reorder_flat[el.bus_or_id] in bus_ex and reorder_flat[el.bus_ex_id] in bus_or)]
        # lines = [el.id for el in self.gridmodel.get_lines() if (el.bus_or_id in reorder_flat[bus_or] and el.bus_ex_id in reorder_flat[bus_ex]) or (el.bus_or_id in reorder_flat[bus_ex] and el.bus_ex_id in reorder_flat[bus_or])]
        # tmp_or = [reorder_flat[el.bus_or_id]  for el in self.gridmodel.get_lines()]
        # tmp_ex = [reorder_flat[el.bus_ex_id]  for el in self.gridmodel.get_lines()]
        # lines = [el.id for el in self.gridmodel.get_trafos() if (reorder_flat[el.bus_hv_id] in bus_or and reorder_flat[el.bus_lv_id] in bus_ex) or (reorder_flat[el.bus_hv_id] in bus_ex and reorder_flat[el.bus_lv_id] in bus_or)]
        # tmp_or = [reorder_flat[el.bus_hv_id]  for el in self.gridmodel.get_trafos()]
        # tmp_ex = [reorder_flat[el.bus_lv_id]  for el in self.gridmodel.get_trafos()]
        # lines_ref = [el.id for el in self.ref_samecase.get_lines() if (el.bus_or_id in bus_or and el.bus_ex_id in bus_ex) or (el.bus_or_id in bus_ex and el.bus_ex_id in bus_or)]
        # if not lines_ref:
        #     lines_ref = [el.id for el in self.ref_samecase.get_trafos() if (el.bus_hv_id in bus_or and el.bus_lv_id in bus_ex) or (el.bus_hv_id in bus_ex and el.bus_lv_id in bus_or)]
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
        # bus_id = self.el_ids[0]
        # reorder = np.argsort([int(el.lstrip("VL").rstrip("0").rstrip("_")) for el in bus_id.index]).reshape(1, -1)
        if v_ls_ref is not None:
            max_ = np.abs(v_ls[reorder] - v_ls_ref).max()
            assert max_ <= self.tol_eq, f"error for vresults for dc: {max_:.2e}"
        tmp_ = self.gridmodel.get_dcYbus_solver()[reorder.T, reorder].todense() - self.ref_samecase.get_dcYbus_solver().todense()
        max_ = np.abs(tmp_).max()
        mat_ls = self.gridmodel.get_dcYbus_solver()[reorder.T, reorder].todense()
        mat_pp = self.ref_samecase.get_dcYbus_solver().todense()
        assert max_ <= self.tol_eq, f"error for dcYbus: {max_:.2e}"
        # check Sbus without slack
        Sbus_ordered = self.gridmodel.get_dcSbus_solver()[reorder].reshape(-1)
        if slack_id > 0:
            max_ = np.abs(Sbus_ordered[:slack_id] - self.ref_samecase.get_dcSbus_solver()[:slack_id]).max()
            assert max_ <= self.tol_eq, f"error for dc Sbus: {max_:.2e}"
        if slack_id != self.gridmodel.get_dcSbus_solver().shape[0] - 1:
            max_ = np.abs(Sbus_ordered[(slack_id+1):] - self.ref_samecase.get_dcSbus_solver()[(slack_id+1):]).max()
            assert max_ <= self.tol_eq, f"error for dc Sbus: {max_:.2e}"

        # same in AC
        v_ls = self.gridmodel.ac_pf(self.V_init_ac, 2, self.tol)
        v_ls_ref = self.ref_samecase.ac_pf(self.V_init_ac, 2, self.tol)
        max_ = np.abs(self.gridmodel.get_Ybus_solver()[reorder.T, reorder] - self.ref_samecase.get_Ybus_solver()).max()
        assert max_ <= self.tol_eq, f"error for Ybus: {max_:.2e}"
        # check Sbus without slack
        Sbus_ordered = self.gridmodel.get_Sbus_solver()[reorder].reshape(-1)
        if slack_id > 0:
            max_ = np.abs(Sbus_ordered[:slack_id] - self.ref_samecase.get_Sbus_solver()[:slack_id]).max() 
            assert max_ <= self.tol_eq, f"error for dc Sbus: {max_:.2e}"
        if slack_id != self.gridmodel.get_Sbus_solver().shape[0] - 1:
            max_ = np.abs(Sbus_ordered[(slack_id+1):] - self.ref_samecase.get_Sbus_solver()[(slack_id+1):]).max()
            assert max_ <= self.tol_eq, f"error for dc Sbus : {max_:.2e}"

    def test_dc_pf(self):
        """test I get the same results as pandapower in dc"""
        v_ls = self.gridmodel.dc_pf(self.V_init_dc, 2, self.tol)
        reorder = self.gridmodel._orig_to_ls.reshape(1, -1)
        if self.compare_pp():
            v_ls_ref = self.ref_samecase.dc_pf(self.V_init_dc, 2, self.tol)
            assert np.abs(v_ls[reorder] - v_ls_ref).max() <= self.tol_eq, f"error for vresults for dc: {np.abs(v_ls[reorder] - v_ls_ref).max():.2e}"
        lf.run_dc(self.network_ref)
        if self.compare_pp():
            pdp.rundcpp(self.pp_samecase)
        
        # v_mag not really relevant in dc so i study only va
        v_ang_pypo = self.network_ref.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_ls))
        # np.where(np.abs(v_ang_ls[reorder] - v_ang_pypo) >= 9.)
        if self.compare_pp():
            v_ang_pp = self.pp_samecase.res_bus["va_degree"].values
            assert np.abs(v_ang_ls[reorder] - v_ang_pp).max() <= self.tol_eq, f"error for va results for dc: {np.abs(v_ang_ls[reorder] - v_ang_pp).max():.2e}"
        assert np.abs(v_ang_ls[reorder] - v_ang_pypo).max() <= self.tol_eq, f"error for va results for dc: {np.abs(v_ang_ls[reorder] - v_ang_pypo).max():.2e}"
        
    def test_ac_pf(self):
        # run the powerflows
        v_ls = self.gridmodel.ac_pf(1.0 * self.V_init_ac, 10, self.tol)
        reorder = self.gridmodel._orig_to_ls.reshape(1, -1)
        if self.compare_pp():
            v_ls_ref = self.ref_samecase.ac_pf(1.0 * self.V_init_ac, 10, self.tol)   
            assert np.abs(v_ls[reorder] - v_ls_ref).max() <= self.tol_eq, f"error for vresults for ac: {np.abs(v_ls[reorder] - v_ls_ref).max():.2e}"
        
        try:
            param = lf.Parameters(voltage_init_mode=pp._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
                                  transformer_voltage_control_on=False,
                                  use_reactive_limits=False,
                                  shunt_compensator_voltage_control_on=False,
                                  phase_shifter_regulation_on=False,
                                  distributed_slack=False,
                                  provider_parameters={"slackBusSelectionMode": "NAME",
                                                       "slackBusesIds": self.network_ref.get_buses().iloc[self.get_slackbus_id()].name}
                                  ) 
        except TypeError:
            param = lf.Parameters(voltage_init_mode=pp._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
                                  transformer_voltage_control_on=False,
                                  no_generator_reactive_limits=True,  # documented in the doc but apparently fails
                                  phase_shifter_regulation_on=False,
                                  simul_shunt=False,  # documented in the doc but apparently fails
                                  distributed_slack=False,
                                  provider_parameters={"slackBusSelectionMode": "NAME",
                                                       "slackBusesIds": self.network_ref.get_buses().iloc[self.get_slackbus_id()].name}
                                  ) 

        res_pypow = lf.run_ac(self.network_ref, parameters=param)
        bus_ref_kv = self.network_ref.get_voltage_levels().loc[self.network_ref.get_buses()["voltage_level_id"].values]["nominal_v"].values
        v_mag_pypo = self.network_ref.get_buses()["v_mag"].values / bus_ref_kv
        v_ang_pypo = self.network_ref.get_buses()["v_angle"].values
        if self.compare_pp():
            pdp.runpp(self.pp_samecase, init="flat")

        # check that pypow solution is "feasible" as seen by lightsim2grid
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            v_pypow = v_mag_pypo * np.exp(1j * np.deg2rad(v_ang_pypo))
            v_pypow_for_ls = self.V_init_dc.copy()
            v_pypow_for_ls[self.gridmodel._orig_to_ls] = v_pypow
            v_pypow_ls = self.gridmodel.check_solution(v_pypow_for_ls, False)
            assert np.abs(v_pypow_ls).max() <= 10. * self.tol_eq, f"error when checking results of pypowsybl in lightsim2grid: {np.abs(v_pypow_ls).max():.2e}"
                 
        # check voltage angles
        v_ang_ls = np.rad2deg(np.angle(v_ls))
        if self.compare_pp():
            v_ang_pp = self.pp_samecase.res_bus["va_degree"].values - self.pp_samecase.ext_grid["va_degree"].values
            assert np.abs(v_ang_ls[reorder] - v_ang_pp).max() <= self.tol_eq, f"error for va results for ac: {np.abs(v_ang_ls[reorder] - v_ang_pp).max():.2e}"
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            assert np.abs(v_ang_ls[reorder] - v_ang_pypo).max() <= self.tol_eq, f"error for va results for ac: {np.abs(v_ang_ls[reorder] - v_ang_pypo).max():.2e}"
        
        # check voltage magnitudes
        v_mag_ls = np.abs(v_ls)
        if self.compare_pp():
            v_mag_pp = self.pp_samecase.res_bus["vm_pu"].values
            assert np.abs(v_mag_ls[reorder] - v_mag_pp).max() <= self.tol_eq, f"error for va results for dc: {np.abs(v_mag_ls[reorder] - v_mag_pp).max():.2e}"
        if res_pypow[0].status == pp._pypowsybl.LoadFlowComponentStatus.CONVERGED:
            assert np.abs(v_mag_ls[reorder] - v_mag_pypo).max() <= self.tol_eq, f"error for va results for dc: {np.abs(v_mag_ls[reorder] - v_mag_pypo).max():.2e}"
        

        
class TestCase14FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    pass


class TestCase14FromPypo(TestCase14FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
        
        
class TestCase30FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 30"""
    # unittest.TestCase does not work because of https://github.com/powsybl/pypowsybl/issues/644
    def get_pypo_grid(self):
        return pp.network.create_ieee30()
    
    def get_equiv_pdp_grid(self):
        return pn.case_ieee30()
 
 
class TestCase30FromPypo(TestCase30FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
    
           
        
class TestCase57FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 57"""
    # does not appear to be the same grid !
    def get_pypo_grid(self):
        res = pp.network.create_ieee57()
        return res
    
    def get_equiv_pdp_grid(self):
        return pn.case57()
    
    def compare_pp(self):
        return super().compare_pp() and False
    
    def get_tol_eq(self):
        return 1e-4  # otherwise vangle from pypowsybl and pandapower does not match


class TestCase57FromPypo(TestCase57FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
                    
        
class TestCase118FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 118: does not work because of https://github.com/e2nIEE/pandapower/issues/2131"""
    # does not work because of https://github.com/e2nIEE/pandapower/issues/2131
    def get_pypo_grid(self):
        return pp.network.create_ieee118()
    
    def get_equiv_pdp_grid(self):
        return pn.case118()
    
    def get_slackbus_id(self):
        return 68
    
    def get_tol_eq(self):
        # return 3e-3  # otherwise vangle from pypowsybl and pandapower does not match
        return 3e-5  # otherwise vangle from pypowsybl and pandapower does not match
    
    def compare_pp(self):
        return super().compare_pp() and False


class TestCase118FromPypo(TestCase118FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
                
        
class TestCase300FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub):
    """compare from the ieee 300"""
    # does not work because of phase tap changer
    def get_pypo_grid(self):
        res = pp.network.create_ieee300()
        df = res.get_shunt_compensators()[["connected"]] # "b", "g", 
        df["connected"] = False
        res.update_shunt_compensators(df)
        return res
    
    def get_equiv_pdp_grid(self):
        return pn.case300()
    
    def get_tol_eq(self):
        # return 3e-3  # otherwise vangle from pypowsybl and pandapower does not match
        return 3e-5  # otherwise vangle from pypowsybl and pandapower does not match
    
    def compare_pp(self):
        return super().compare_pp() and False
    
    def get_slackbus_id(self):
        # does not work with PP, probably bus not ordered the same
        return 257


class TestCase300FromPypo(TestCase300FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False


if __name__ == "__main__":
    unittest.main()
    