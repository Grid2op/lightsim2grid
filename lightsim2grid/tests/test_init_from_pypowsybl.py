# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import pdb
from packaging import version
import pypowsybl as pp
import pypowsybl.loadflow as lf
import numpy as np
import unittest
import warnings

from lightsim2grid.gridmodel import init_from_pypowsybl, GridModel

try:
    import pandapower.networks as pn
    import pandapower as pdp
    from lightsim2grid.gridmodel import init_from_pandapower
    PDP_AVAIL = True
except ImportError:
    # pandapower not available, eg if testing with numpy 2
    PDP_AVAIL = False

from global_var_tests import (
    MAX_PP2_DATAREADER,
    CURRENT_PP_VERSION,
    VERSION_PHASESHIFT_OK_PYPOW,
    CURRENT_PYPOW_VERSION)
from test_match_with_pypowsybl.utils_for_slack import (
    get_pypowsybl_parameters,
    get_same_slack
)

if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
    # deactivate compat with recent pandapower version, for now
    PDP_AVAIL = False


class AuxInitFromPyPowSyBlBusesForSub:    
    def use_buses_for_sub(self):
        return True
    
    def get_pypo_grid_name(self):
        return "ieee14"
    
    def get_equiv_pdp_grid(self):
        return pn.case14()
    
    def get_tol_eq(self):
        return 1e-6
    
    def get_tol_eq_kcl(self):
        return 1e-5
    
    def compare_pp(self):
        """will this test suite compare pypowsybl and pandapower (cannot be used for ieee57 or ieee118)"""
        return PDP_AVAIL
    
    def setUp(self) -> None:
        self.pypo_grid_name = self.get_pypo_grid_name()
        self.pypo_slack_name, self.ls_slack_bus_id = get_same_slack(self.pypo_grid_name)
        self.network_ref = getattr(pp.network, f"create_{self.pypo_grid_name}")()
        
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
                                                          slack_bus_id=self.ls_slack_bus_id,
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
        self.tol_eq_kcl = self.get_tol_eq_kcl()
        return super().setUp()
    
    def test_basic(self):
        """check that all elements are ok"""
        assert len(self.gridmodel.get_lines()) == self.network_ref.get_lines().shape[0]    
        assert len(self.gridmodel.get_loads()) == self.network_ref.get_loads().shape[0]    
        assert len(self.gridmodel.get_generators()) == self.network_ref.get_generators().shape[0]    
        assert len(self.gridmodel.get_trafos()) == self.network_ref.get_2_windings_transformers().shape[0]    
        assert len(self.gridmodel.get_shunts()) == self.network_ref.get_shunt_compensators().shape[0]    
    
    def test_compare_pp(self):
        """test I get the same results than pypowsybl, in AC"""
        if not self.can_pp:
            self.skipTest("no equivalent pandapower grid")
        if not self.compare_pp():
            self.skipTest("It is expected that pandapower and pypowsybl gives different results") 
            
        # check the SBus and Ybus in DC (i need to run powerflow for that)      
        v_ls = self.gridmodel.dc_pf(1.0 * self.V_init_dc, 2, self.tol)
        v_ls_ref = None
        if self.ref_samecase is not None:
            v_ls_ref = self.ref_samecase.dc_pf(1.0 * self.V_init_dc, 2, self.tol)
        slack_id = self.ls_slack_bus_id
        reorder = self.gridmodel._orig_to_ls.reshape(1, -1)

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
        """test I get the same results as pypowsybl in dc"""
        # if CURRENT_PYPOW_VERSION >= MIN_PYPO_DC_NOT_WORKING:
            # self.skipTest("Test not correct: pypowsybl change DC approx, see https://github.com/powsybl/pypowsybl/issues/1127")
        v_ls = self.gridmodel.dc_pf(self.V_init_dc, 2, self.tol)
        reorder = self.gridmodel._orig_to_ls.reshape(1, -1)
        if self.compare_pp():
            v_ls_ref = self.ref_samecase.dc_pf(self.V_init_dc, 2, self.tol)
            assert np.abs(v_ls[reorder] - v_ls_ref).max() <= self.tol_eq, f"error for vresults for dc: {np.abs(v_ls[reorder] - v_ls_ref).max():.2e}"
        # see https://github.com/powsybl/pypowsybl/issues/1127#issuecomment-3581713875
        params = get_pypowsybl_parameters(self.pypo_slack_name)
        lf.run_dc(self.network_ref, parameters=params)
        if self.compare_pp():
            pdp.rundcpp(self.pp_samecase)
        
        # v_mag not really relevant in dc so i study only va
        v_ang_pypo = self.network_ref.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_ls))
        # for case300
        # self.gridmodel.get_2_windings_transformers()[85]
        # shift of 0.19896753472735357 radian
        # 161 (bus1_id) and 424 (bus2_id)
        # DC
        # -9.859724243884386
        # 9.859724243884614
        # AC
        # 83.57031489368106
        # -83.56243958793598
        
        # self.network_ref.get_2_windings_transformers(all_attributes=True).iloc[85]
        # alpha = 11.4 deg
        # bus1_id VL196_0, bus2_id VL196_1
        # self.network_ref.get_buses().loc[["VL196_0", "VL196_1"]]
        # DC
        # p1                         88.74429
        # p2                        -88.74429
        # AC
        # p1                          83.570175
        # p2                           -83.5623
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
        
        param = get_pypowsybl_parameters(self.pypo_slack_name)
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
            assert np.abs(v_pypow_ls).max() <= self.tol_eq_kcl, f"error when checking results of pypowsybl in lightsim2grid: {np.abs(v_pypow_ls).max():.2e}"
                 
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
    def get_pypo_grid_name(self):
        return "ieee30"
    
    def get_equiv_pdp_grid(self):
        return pn.case_ieee30()
 
 
class TestCase30FromPypo(TestCase30FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
    
           
        
class TestCase57FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 57"""
    # does not appear to be the same grid !
    def get_pypo_grid_name(self):
        return "ieee57"
    
    def get_equiv_pdp_grid(self):
        return pn.case57()
    
    def compare_pp(self):
        return False
    
    def get_tol_eq_kcl(self):
        return 1.53e-04
    
    def get_tol_eq(self):
        return 5.47e-05
                    
    
class TestCase57FromPypo(TestCase57FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
        
class TestCase118FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 118: does not work because of https://github.com/e2nIEE/pandapower/issues/2131"""
    # does not work because of https://github.com/e2nIEE/pandapower/issues/2131
    def get_pypo_grid_name(self):
        return "ieee118"
    
    def get_equiv_pdp_grid(self):
        return pn.case118()
    
    def get_tol_eq_kcl(self):
        return 2.72e-4
    
    def get_tol_eq(self):
        return 1.46e-5
    
    def compare_pp(self):
        return super().compare_pp() and False


class TestCase118FromPypo(TestCase118FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False
                
        
class TestCase300FromPypoBusesForSub(AuxInitFromPyPowSyBlBusesForSub, unittest.TestCase):
    """compare from the ieee 300"""
    def get_pypo_grid_name(self):
        if CURRENT_PYPOW_VERSION < VERSION_PHASESHIFT_OK_PYPOW:
            self.skipTest("phase shifter are not supported for "
                          "this version of pypowsybl")
        return "ieee300"
    
    def get_equiv_pdp_grid(self):
        return pn.case300()
    
    def get_tol_eq(self):
        return 2.2e-5
    
    def get_tol_eq_kcl(self):
        return 3.7e-4
    
    def compare_pp(self):
        return super().compare_pp() and False


class TestCase300FromPypo(TestCase300FromPypoBusesForSub):
    def use_buses_for_sub(self):
        return False


class TestBusesForSub_dosort(unittest.TestCase):
    def do_i_sort(self):
        return True
    
    def get_align_vector(self):
        return np.asarray(
            [ 5,  0,  1,  2,  3,  4,  6,  7,  8,  9, 10, 11, 12, 13],
            dtype=int)
    
    def get_ls_to_orig(self):
        return self.get_align_vector()
    
    def get_gen_slack_bus_id(self):
        return 5
    
    def get_sub_names(self):
        return np.asarray(
            ['VL1', 'VL10', 'VL11', 'VL12', 'VL13', 'VL14', 'VL2', 
             'VL3', 'VL4', 'VL5', 'VL6', 'VL7', 'VL8', 'VL9'],
            dtype=str
        )
        
    def get_sub_names_b4s(self):
        return np.asarray(
            ['VL10_0', 'VL11_0', 'VL12_0', 'VL13_0', 'VL14_0', 'VL1_0', 
            'VL2_0', 'VL3_0', 'VL4_0', 'VL5_0', 'VL6_0', 'VL7_0', 'VL8_0', 'VL9_0'],
            dtype=str
        )
        
    def setUp(self):
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # self.path = os.path.join(dir_path, "case_14_iidm")
        # self.file_name = "grid.xiidm"
        # self.pypow_grid = pp.network.load(os.path.join(self.path, self.file_name))
        self.pypo_slack_name, _ = get_same_slack("ieee14")
        self.pypow_grid = pp.network.create_ieee14()
        self.ls_grid_b4s : GridModel = init_from_pypowsybl(
            self.pypow_grid,
            buses_for_sub=True,
            sort_index=self.do_i_sort(), 
            gen_slack_id=0)
        self.ls_grid : GridModel = init_from_pypowsybl(
            self.pypow_grid,
            buses_for_sub=False,
            sort_index=self.do_i_sort(), 
            gen_slack_id=0)
        self.align_vect = self.get_align_vector()
        return super().setUp()
    
    def test_some_differences(self):
        # correct generator is slack
        assert self.ls_grid.get_generators()[0].is_slack
        assert self.ls_grid_b4s.get_generators()[0].is_slack
        # no other generators are slack
        assert abs(self.ls_grid.get_generators()[0].slack_weight - 1.) < 1e-6
        assert abs(self.ls_grid_b4s.get_generators()[0].slack_weight - 1.) < 1e-6
        # generator is at the correct bus
        assert self.ls_grid.get_generators()[0].bus_id == 0
        assert self.ls_grid_b4s.get_generators()[0].bus_id == self.get_gen_slack_bus_id()
        # 'correct' ls_to_orig vector
        assert (
            self.ls_grid._ls_to_orig ==
            self.get_ls_to_orig()
        ).all()
        assert (
            self.ls_grid_b4s._ls_to_orig ==
            [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13]
        ).all()      
          
        # correct substation names
        assert (
            self.ls_grid.get_substation_names() ==
            self.get_sub_names()
        ).all()
        assert (
            self.ls_grid_b4s.get_substation_names() ==
            self.get_sub_names_b4s()
        ).all()
        
    def test_same_res_dc(self):
        v_init = np.ones(self.ls_grid.get_bus_vn_kv().shape[0], dtype=complex)
        v_dc = self.ls_grid.dc_pf(v_init, 10, 1e-6)
        v_dc_b4s = self.ls_grid_b4s.dc_pf(v_init, 10, 1e-6)
        assert v_dc.shape[0] > 0
        assert v_dc_b4s.shape[0] > 0
        
        # test same Sbus
        sbus = self.ls_grid.get_dcSbus_solver()
        sbus_b4s = self.ls_grid_b4s.get_dcSbus_solver()
        assert (np.abs(sbus_b4s[self.align_vect] - sbus) < 1e-6).all()
        
        # test same Ybus
        Ybus = self.ls_grid.get_dcYbus_solver()
        Ybus_b4s = self.ls_grid_b4s.get_dcYbus_solver()
        assert (np.abs(Ybus_b4s[self.align_vect.reshape(-1,1), self.align_vect] - Ybus) > 1e-6).size == 0
        
        # test resulting vectors are the same
        assert (np.abs(v_dc_b4s[self.align_vect] - v_dc) < 1e-6).all()
        
        # now test they both match pypowsybl results
        param = get_pypowsybl_parameters(self.pypo_slack_name)
        res_pypow = lf.run_dc(self.pypow_grid, parameters=param)
        vl_pypow = self.pypow_grid.get_voltage_levels()
        vl_pypow["order_pypo"] = np.arange(14)
        
        # voltage angle
        v_ang_pypo = self.pypow_grid.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_dc))
        v_pypo_aligned = v_ang_pypo[vl_pypow.loc[self.ls_grid.get_substation_names(), "order_pypo"]]
        assert (np.abs(v_pypo_aligned - v_ang_ls) < 1e-6).all()
        
    def test_same_res_ac(self):
        v_init = np.ones(self.ls_grid.get_bus_vn_kv().shape[0], dtype=complex)
        v_ac = self.ls_grid.ac_pf(v_init, 10, 1e-6)
        v_ac_b4s = self.ls_grid_b4s.ac_pf(v_init, 10, 1e-6)
        assert v_ac.shape[0] > 0
        assert v_ac_b4s.shape[0] > 0
        
        # test same Sbus
        sbus = self.ls_grid.get_Sbus_solver()
        sbus_b4s = self.ls_grid_b4s.get_Sbus_solver()
        assert (np.abs(sbus_b4s[self.align_vect] - sbus) < 1e-6).all()
        
        # test same Ybus
        Ybus = self.ls_grid.get_Ybus_solver()
        Ybus_b4s = self.ls_grid_b4s.get_Ybus_solver()
        assert (np.abs(Ybus_b4s[self.align_vect.reshape(-1,1), self.align_vect] - Ybus) > 1e-6).size == 0
        
        # test resulting vectors are the same
        assert (np.abs(v_ac_b4s[self.align_vect] - v_ac) < 1e-6).all()
        
        # now test they both match pypowsybl results
        param = get_pypowsybl_parameters(self.pypo_slack_name)
        res_pypow = lf.run_ac(self.pypow_grid, parameters=param)
        vl_pypow = self.pypow_grid.get_voltage_levels()
        vl_pypow["order_pypo"] = np.arange(14)
        
        v_ang_pypo = self.pypow_grid.get_buses()["v_angle"].values
        v_ang_ls = np.rad2deg(np.angle(v_ac))
        v_pypo_aligned = v_ang_pypo[vl_pypow.loc[self.ls_grid.get_substation_names(), "order_pypo"]]
        assert (np.abs(v_pypo_aligned - v_ang_ls) < 1e-6).all()
        
        # voltage angle
        bus_df = self.pypow_grid.get_buses()
        v_ang_pypo = bus_df["v_mag"].values / self.pypow_grid.get_voltage_levels().loc[bus_df["voltage_level_id"], "nominal_v"].values
        v_ang_ls = np.abs(v_ac)
        v_pypo_aligned = v_ang_pypo[vl_pypow.loc[self.ls_grid.get_substation_names(), "order_pypo"]]
        assert (np.abs(v_pypo_aligned - v_ang_ls) < 1e-6).all()
        
        
class TestBusesForSub_nosort(TestBusesForSub_dosort):   
    def do_i_sort(self):
        return False     
    
    def get_gen_slack_bus_id(self):
        return 0
    
    def get_ls_to_orig(self):
        return np.arange(14)
    
    def get_align_vector(self):
        return np.arange(14)
    
    def get_sub_names(self):
        return np.asarray(
            ['VL1', 'VL2', 'VL3', 'VL4', 'VL5', 'VL6', 'VL7', 'VL8', 
             'VL9', 'VL10', 'VL11', 'VL12', 'VL13', 'VL14'],
            dtype=str
        )
        
    def get_sub_names_b4s(self):
        return np.asarray(
            ['VL1_0', 'VL2_0', 'VL3_0', 'VL4_0', 'VL5_0', 'VL6_0', 
             'VL7_0', 'VL8_0', 
             'VL9_0', 'VL10_0', 'VL11_0', 'VL12_0', 'VL13_0', 'VL14_0'],
            dtype=str
        )
        
    
if __name__ == "__main__":
    unittest.main()
    