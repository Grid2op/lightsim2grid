# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import numpy as np
import pandas as pd
from packaging import version as version_packaging 

import pypowsybl.network as pp_network
import pypowsybl as pp
import pypowsybl.loadflow as pp_lf

from lightsim2grid.gridmodel import init_from_pypowsybl

from global_var_tests import CURRENT_PYPOW_VERSION, VERSION_PHASESHIFT_OK_PYPOW
from test_match_with_pypowsybl.utils_for_slack import (
    get_pypowsybl_parameters,
    get_same_slack
)
import warnings
import pdb


class BaseTests:    
    def setUp(self, grid_nm="ieee118"):
        fun_create = getattr(pp_network, f"create_{grid_nm}")
        self.net_ref = fun_create()
        self.net_datamodel = fun_create()

        # initialize constant stuff
        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 3e-4  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)
        self.tol_V = 1e-6
        self.slack_vl_id_pypow, self.ls_slack = get_same_slack(grid_nm)
        self.gen_slack_id = 29
        # initialize and use converters
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.model = init_from_pypowsybl(self.net_datamodel,
                                             slack_bus_id=self.ls_slack,
                                             buses_for_sub=True,
                                             sort_index=False)

    def assert_equal(self, tmp, ref, error="", tol=None):
        assert np.all(tmp.shape == ref.shape), f"vector does not have the same shape for {error}"
        if tol is None:
            tol = self.tol_test
        assert np.max(np.abs(tmp - ref)) <= tol, f"{error}: {np.max(np.abs(tmp - ref))} > {tol}"
        assert np.mean(np.abs(tmp - ref)) <= tol, f"{error}: {np.mean(np.abs(tmp - ref))} > {tol}"

    def check_res(self, Vfinal, net):
        assert Vfinal.shape[0] > 0, "powerflow diverged !"

        # final voltage vector
        va_rad_ls = np.angle(Vfinal)
        vm_pu_ls = np.abs(Vfinal)
        va_rad_pypow = np.deg2rad(net.get_buses()["v_angle"].values)
        vm_pu_pypow = (
            net.get_buses()["v_mag"].values / 
            net.get_voltage_levels().loc[net.get_buses()["voltage_level_id"], "nominal_v"].values
            )
        self.assert_equal(va_rad_ls, va_rad_pypow, "error for voltage angle", self.tol_V)
        self.assert_equal(vm_pu_ls, vm_pu_pypow, "error for voltage magnitude", self.tol_V)
    
        # check lines
        l_is = self.net_ref.get_lines()["connected1"]
        por, qor, vor, aor = self.model.get_line_res1()
        self.assert_equal(por[l_is], net.get_lines()["p1"].values[l_is], "error for p_from")
        self.assert_equal(qor[l_is], net.get_lines()["q1"].values[l_is], "error for q_from")
        self.assert_equal(aor[l_is], net.get_lines()["i1"].values[l_is] * 1e-3, "error for i_from_ka")
        vor_pp = net.get_buses().loc[net.get_lines().loc[l_is, "bus1_id"].values]["v_mag"].values
        self.assert_equal(vor[l_is], vor_pp, "error for vor_pp")

        # check trafo
        f_is = self.net_ref.get_2_windings_transformers()["connected2"]
        plv, qlv, vlv, alv = self.model.get_trafo_res2()
        self.assert_equal(plv[f_is], net.get_2_windings_transformers()["p2"].values[f_is], "error for p_lv_mw")
        self.assert_equal(qlv[f_is], net.get_2_windings_transformers()["q2"].values[f_is], "error for q_lv_mvar")
        self.assert_equal(alv[f_is], net.get_2_windings_transformers()["i2"].values[f_is] * 1e-3, "error for i_lv_ka")
        vlv_pp = net.get_buses().loc[net.get_2_windings_transformers().loc[f_is, "bus2_id"].values]["v_mag"].values
        self.assert_equal(vlv[f_is], vlv_pp, "error for vlv_pp")

        # check loads
        l_is = self.net_ref.get_loads()["connected"].values
        load_p, load_q, load_v = self.model.get_loads_res()
        self.assert_equal(load_p[l_is], net.get_loads()["p"].values[l_is], "error for load p_mw")
        self.assert_equal(load_q[l_is], net.get_loads()["q"].values[l_is], "error for load q_mvar")

        # check generators
        g_is = self.net_ref.get_generators()["connected"]
        prod_p, prod_q, prod_v = self.model.get_gen_res()
        
        # # test the slack bus is properly modeled as a generator
        # assert np.abs(np.sum(net._ppc["gen"][:, 1]) - np.sum(prod_p)) <= self.tol_test
        # if len(prod_p) != g_is.shape[0]:
        #     # it means a generator has been added for the slack bus
        #     prod_p = prod_p[:-1]
        #     prod_q = prod_q[:-1]
        #     prod_v = prod_v[:-1]
        self.assert_equal(prod_p[g_is], -net.get_generators()["p"].values[g_is], "error for gen p_mw")
        self.assert_equal(prod_q[g_is], -net.get_generators()["q"].values[g_is], "error for gen q_mvar")
        v_gen_pp = net.get_buses().loc[net.get_generators().loc[g_is, "bus_id"].values]["v_mag"].values
        self.assert_equal(prod_v[g_is], v_gen_pp, "error for prod_v")
        
        # check shunts
        s_is = self.net_ref.get_shunt_compensators()["connected"]
        shunt_p, shunt_q, shunt_v = self.model.get_shunts_res()
        self.assert_equal(shunt_p[s_is], net.get_shunt_compensators()["p"].values[s_is], "error for shunt p_mw")
        self.assert_equal(shunt_q[s_is], net.get_shunt_compensators()["q"].values[s_is], "error for shunt q_mvar")

    def make_v0(self, net):
        V0 = np.full(net.get_buses().shape[0],
                     fill_value=1.0,
                     dtype=complex)
        return V0

    def run_me_pf(self, V0):
        res = self.model.ac_pf(V0, self.max_it, self.tol)
        return res

    def run_ref_pf(self, net):
        params = get_pypowsybl_parameters(self.slack_vl_id_pypow)
        res_pypow = pp_lf.run_ac(net, params)
        # overrun in test_AC or test_DC below !

    def do_i_skip(self, func_name):
        # self.skipTest("dev")
        pass

    def _run_both_pf(self, net):
        V0 = self.make_v0(net)
        self.run_ref_pf(net)
        Vfinal = self.run_me_pf(V0)
        return Vfinal

    def test_function_works(self):
        self.do_i_skip("test_function_works")
        self.model.get_loads_status()
        self.model.get_shunts_status()
        self.model.get_gen_status()
        self.model.get_lines_status()
        self.model.get_trafo_status()

    def test_deactivate_index_out_of_bound(self):
        self.do_i_skip("test_deactivate_index_out_of_bound")
        n_loads = self.net_datamodel.get_loads().shape[0]
        with self.assertRaises(IndexError):
            self.model.deactivate_load(n_loads)
        n_gen = self.net_datamodel.get_generators().shape[0]
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.deactivate_gen(n_gen+1)
        n_trafo = self.net_datamodel.get_2_windings_transformers().shape[0]
        with self.assertRaises(IndexError):
            self.model.deactivate_trafo(n_trafo)
        n_line = self.net_datamodel.get_lines().shape[0]
        with self.assertRaises(IndexError):
            self.model.deactivate_powerline(n_line)
        n_shunt = self.net_datamodel.get_shunt_compensators().shape[0]
        with self.assertRaises(IndexError):
            self.model.deactivate_shunt(n_shunt)

    def test_changebus_index_out_of_bound(self):
        self.do_i_skip("test_changebus_index_out_of_bound")
        n_loads = self.net_datamodel.get_loads().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus_load(n_loads, 1)
        n_gen = self.net_datamodel.get_generators().shape[0]
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_bus_gen(n_gen+1, 1)
        n_shunt = self.net_datamodel.get_shunt_compensators().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus_shunt(n_shunt, 1)
        n_line = self.net_datamodel.get_lines().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus1_powerline(n_line, 1)
        with self.assertRaises(IndexError):
            self.model.change_bus2_powerline(n_line, 1)
        n_trafo = self.net_datamodel.get_2_windings_transformers().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus1_trafo(n_trafo, 1)
        with self.assertRaises(IndexError):
            self.model.change_bus2_trafo(n_trafo, 1)

    def test_changebus_newbus_out_of_bound(self):
        self.do_i_skip("test_changebus_newbus_out_of_bound")
        newbusid = self.net_datamodel.get_buses().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_bus_load(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_gen(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus_shunt(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus1_powerline(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus2_powerline(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus1_trafo(0, newbusid)
        with self.assertRaises(IndexError):
            self.model.change_bus2_trafo(0, newbusid)

    def test_changesetpoint_out_of_bound(self):
        self.do_i_skip("test_changesetpoint_out_of_bound")
        n_loads = self.net_datamodel.get_loads().shape[0]
        n_gens = self.net_datamodel.get_generators().shape[0]
        n_shunt = self.net_datamodel.get_shunt_compensators().shape[0]
        with self.assertRaises(IndexError):
            self.model.change_p_load(n_loads, 1)
        with self.assertRaises(IndexError):
            self.model.change_q_load(n_loads, 1)
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_p_gen(n_gens + 1, 1)
        with self.assertRaises(IndexError):
            # +1 is added for gen because in some cases, gen is added in gridmodel for the slack bus
            self.model.change_v_gen(n_gens + 1, 1)
        with self.assertRaises(IndexError):
            self.model.change_p_shunt(n_shunt, 1)
        with self.assertRaises(IndexError):
            self.model.change_p_shunt(n_shunt, 1)

    def test_pf(self):
        """
        Reference without modifying anything
        """
        self.do_i_skip("test_pf")
        # compute a powerflow on a net without anything
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)
        
        # check no error when retrieving these              
        self.model.get_V_solver()          
        self.model.get_Va_solver()            
        self.model.get_Vm_solver()
        if hasattr(self, "dc") and not self.dc:
            # does not make sense in dc powerflow              
            self.model.get_J_solver()
            self.model.get_V()     
            self.model.get_Va()  
            self.model.get_Vm()            

    def test_pf_disco_gen(self):
        self.do_i_skip("test_pf_disco_gen")
        # self.net_ref.gen["in_service"][0] = False
        self.net_ref.update_generators(
            id=self.net_ref.get_generators().index[0],
            connected=False
        )
        self.model.deactivate_gen(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_load(self):
        self.do_i_skip("test_pf_disco_load")
        # self.net_ref.load["in_service"][0] = False
        self.net_ref.update_loads(
            id=self.net_ref.get_loads().index[0],
            connected=False
        )
        self.model.deactivate_load(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_line(self):
        self.do_i_skip("test_pf_disco_line")
        # self.net_ref.line["in_service"][0] = False
        self.net_ref.update_lines(
            id=self.net_ref.get_lines().index[0],
            connected1=False,
            connected2=False
        )
        self.model.deactivate_powerline(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_shunt(self):
        self.do_i_skip("test_pf_disco_shunt")
        # self.net_ref.shunt["in_service"][0] = False
        self.net_ref.update_shunt_compensators(
            id=self.net_ref.get_shunt_compensators().index[0],
            connected=False,
        )
        self.model.deactivate_shunt(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_disco_trafo(self):
        self.do_i_skip("test_pf_disco_trafo")
        # self.net_ref.trafo["in_service"][0] = False
        self.net_ref.update_2_windings_transformers(
            id=self.net_ref.get_2_windings_transformers().index[0],
            connected1=False,
            connected2=False
        )
        self.model.deactivate_trafo(0)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_reactivate(self):
        # i deactivate everything, run a powerflow, and check that reactivating everything and supposes that the results
        # is the same
        self.do_i_skip("test_reactivate")
        self.run_ref_pf(self.net_ref)
        V0 = self.make_v0(self.net_ref)

        # i disconnect a load, the reconnect it
        self.model.deactivate_load(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_load(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_gen(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_gen(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_powerline(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_powerline(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

        self.model.deactivate_trafo(0)
        Vfinal = self.run_me_pf(V0)
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        self.model.reactivate_trafo(0)
        Vfinal = self.run_me_pf(V0)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_gen(self):
        self.do_i_skip("test_pf_changebus_gen")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        # self.net_ref.gen["bus"][0] = 2
        # self.net_ref.update_generators(
        #     id=self.net_ref.get_generators().index[0],
        #     bus_id="VL3_0"
        # )
        df_gen = self.net_ref.get_generators()
        # TODO test create a generator !
        self.net_ref.create_generators(
            id=f"{df_gen.index[0]}_cpy",
            voltage_level_id="VL3",
            bus_id="VL3_0",
            connectable_bus_id="VL3_0",
            target_p=df_gen.iloc[0]["target_p"],
            target_v=df_gen.iloc[0]["target_v"],
            max_p=df_gen.iloc[0]["max_p"],
            min_p=df_gen.iloc[0]["min_p"],
            voltage_regulator_on=df_gen.iloc[0]["voltage_regulator_on"],
        )
        self.model.change_bus_gen(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_load(self):
        self.do_i_skip("test_pf_changebus_load")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.load["bus"][0] = 2
        self.model.change_bus_load(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_shunt(self):
        self.do_i_skip("test_pf_changebus_shunt")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.shunt["bus"][0] = 2
        self.model.change_bus_shunt(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_lineor(self):
        self.do_i_skip("test_pf_changebus_lineor")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.line["from_bus"][0] = 2
        self.model.change_bus1_powerline(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_lineex(self):
        self.do_i_skip("test_pf_changebus_lineex")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.line["to_bus"][0] = 2
        self.model.change_bus2_powerline(0, 2)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_trafolv(self):
        self.do_i_skip("test_pf_changebus_trafolv")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.trafo["lv_bus"][0] = 5  # was 4 initially, and 4 is connected to 5
        self.model.change_bus2_trafo(0, 5)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changebus_trafohv(self):
        self.do_i_skip("test_pf_changebus_trafohv")
        self.skipTest("changebus does not work with pypowsybl at the moment")
        self.net_ref.trafo["hv_bus"][0] = 29  # was 7 initially, and 7 is connected to 29
        self.model.change_bus1_trafo(0, 29)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeloadp(self):
        self.do_i_skip("test_pf_changeloadp")
        # self.net_ref.load["p_mw"][0] = 50
        self.net_ref.update_loads(
            id=self.net_ref.get_loads().index[0],
            p0=50
        )
        self.model.change_p_load(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeloadq(self):
        self.do_i_skip("test_pf_changeloadq")
        # self.net_ref.load["q_mvar"][0] = 50
        self.net_ref.update_loads(
            id=self.net_ref.get_loads().index[0],
            q0=50
        )
        self.model.change_q_load(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeprodp(self):
        self.do_i_skip("test_pf_changeprodp")
        # self.net_ref.gen["p_mw"][0] = 50
        self.net_ref.update_generators(
            id=self.net_ref.get_generators().index[0],
            target_p=50
        )
        self.model.change_p_gen(0, 50)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeprodv(self):
        self.do_i_skip("test_pf_changeprodv")
        # self.net_ref.gen["vm_pu"][0] = 1.06
        df_vl = self.net_ref.get_voltage_levels()
        self.net_ref.update_generators(
            id=self.net_ref.get_generators().index[0],
            target_v=1.06*df_vl.loc[self.net_ref.get_generators().iloc[0]["voltage_level_id"]]["nominal_v"]
        )
        self.model.change_v_gen(0, 1.06)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeshuntp(self):
        self.do_i_skip("test_pf_changeshuntp")
        p_mw = 10.
        df_shunt = self.net_ref.get_shunt_compensators()
        shunt_0 = df_shunt.iloc[0]
        section_shunt_0 = self.net_ref.get_linear_shunt_compensator_sections().loc[shunt_0.name]
        shunt_kv = self.net_ref.get_voltage_levels().loc[shunt_0.voltage_level_id, "nominal_v"]
        self.net_ref.update_linear_shunt_compensator_sections(
            id=section_shunt_0.name,
            g_per_section=p_mw / shunt_kv**2,
        )
        self.model.change_p_shunt(0, p_mw)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)

    def test_pf_changeshuntq(self):
        self.do_i_skip("test_pf_changeshuntq")
        # self.net_ref.shunt["q_mvar"][0] = 10
        df_shunt = self.net_ref.get_shunt_compensators()
        shunt_0 = df_shunt.iloc[0]
        section_shunt_0 = self.net_ref.get_linear_shunt_compensator_sections().loc[shunt_0.name]
        self.net_ref.update_linear_shunt_compensator_sections(
            id=section_shunt_0.name,
            max_section_count=2,
        )
        self.net_ref.update_shunt_compensators(
            id=shunt_0.name,
            section_count=2  # double the reactive power
        )
        self.model.change_q_shunt(0, 80)  # it was 40 put it to 80
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)
        
    def test_change_trafo_ratio(self):
        self.do_i_skip("test_change_trafo_ratio")
        # update gridmodel
        old_val = self.model.get_trafos()[0].ratio
        nm_trafo = self.model.get_trafos()[0].name
        assert abs(old_val -1.) > 1e-6
        new_val = 1. + 2. * (old_val - 1.)  # multiply actual ratio by 2
        self.model.change_ratio_trafo(0, new_val)
        assert np.allclose(self.model.get_trafos()[0].ratio, new_val)
        
        # check it has an impact on the load flow
        Vfinal = self._run_both_pf(self.net_ref)
        with self.assertRaises(AssertionError):
            self.check_res(Vfinal, self.net_ref)
            
        # update pypowsybl
        pp_tr = self.net_ref.get_2_windings_transformers().loc[nm_trafo]
        ref_vnom1 = self.net_ref.get_voltage_levels().loc[pp_tr["voltage_level1_id"], "nominal_v"]
        
        self.net_ref.update_2_windings_transformers(
            id=nm_trafo,
            rated_u1=ref_vnom1 / new_val)
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)
        
    def test_change_trafo_shift(self):
        self.do_i_skip("test_change_trafo_shift")
        if CURRENT_PYPOW_VERSION < VERSION_PHASESHIFT_OK_PYPOW:
            self.skipTest("phase shifter are not supported for "
                          "this version of pypowsybl")
        self.setUp("ieee300")  # the smaller grid does not have phase shifters
        
        # update gridmodel
        trafo_id = 85
        old_val = self.model.get_trafos()[trafo_id].shift_rad
        nm_trafo = self.model.get_trafos()[trafo_id].name
        assert abs(old_val) > 1e-6
        new_val = 2. * old_val  # multiply actual phase shift by 2
        self.model.change_ratio_trafo(trafo_id, new_val)
        assert np.allclose(self.model.get_trafos()[trafo_id].ratio, new_val)
        
        # check it has an impact on the load flow
        Vfinal = self._run_both_pf(self.net_ref)
        with self.assertRaises(AssertionError):
            self.check_res(Vfinal, self.net_ref)
            
        self.skipTest("Hard to change the alpha of phase tap changer in pypowsbl")
        
        # update pypowsybl
        df = self.net_ref.get_phase_tap_changer_steps().loc[[(nm_trafo, 0)]]
        df["alpha"] *= 2
        df = self.net_ref.get_phase_tap_changer_steps()
        self.net_ref.update_phase_tap_changer_steps(df)  # .reset_index().set_index("id"))
        Vfinal = self._run_both_pf(self.net_ref)
        self.check_res(Vfinal, self.net_ref)


class MakeDCTests(BaseTests, unittest.TestCase):
    def run_me_pf(self, V0):
        self.dc = True
        return self.model.dc_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        params = get_pypowsybl_parameters(self.slack_vl_id_pypow)
        pp_lf.run_dc(net, params)

    def do_i_skip(self, test_nm):
        pass

    def check_res(self, Vfinal, net):
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        va_deg = net.get_buses()["v_angle"].values
        self.assert_equal(np.angle(Vfinal), np.deg2rad(va_deg))

    def test_pf_changeshuntp(self):
        self.skipTest("pypowsybl ignores shunt in dc, lightsim2grid does not")
        
        
class MakeACTests(BaseTests, unittest.TestCase):
    def run_me_pf(self, V0):
        self.dc = False
        return self.model.ac_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        params = get_pypowsybl_parameters(self.slack_vl_id_pypow)
        res_pypow = pp_lf.run_ac(net, params)
        
        # and update the slack...
        try:
            slack_abs = res_pypow[0].slack_bus_results[0].active_power_mismatch
        except AttributeError:
            # legacy pypowsybl version
            slack_abs = res_pypow[0].slack_bus_active_power_mismatch
            
        net.update_generators(
            id=net.get_generators().index[self.gen_slack_id],
            p=net.get_generators().iloc[self.gen_slack_id]["p"]-slack_abs)

    def do_i_skip(self, test_nm):
        pass

    def test_pf_disco_trafo(self):
        self.tol_V = 3e-5
        self.tol_test = 1e-2
        super().test_pf_disco_trafo()

    def test_pf_changeshuntq(self):
        # self.tol_V = 4e-4
        # self.tol_test = 1e-2
        super().test_pf_changeshuntq()
        
if __name__ == "__main__":
    unittest.main()
