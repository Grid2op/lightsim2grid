# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import pickle
import tempfile
import unittest
import warnings
import numpy as np

import pypowsybl.network as pp_network
import pypowsybl.loadflow as pp_lf

import grid2op
from grid2op.Chronics import ChangeNothing

from lightsim2grid import LightSimBackend
from lightsim2grid.network import init_from_pypowsybl


from test_match_with_pypowsybl.utils_for_slack import (
    get_pypowsybl_parameters,
    get_same_slack
)

class BaseDiscoOneSide:      
    """Test the information are correctly computed on the gridmodel side in all 4 cases"""  
    def synch_status_both_side(self):
        return True
    def ignore_status_global(self):
        return False
    
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.model = self.env.backend._grid
        self.model.set_ignore_status_global(self.ignore_status_global())
        self.model.set_synch_status_both_side(self.synch_status_both_side())
    
    def test_pickle(self):
        """test pickle does not lose the info"""
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "test_pickle.pickle"), "wb") as f:
                pickle.dump(self.model, f)
            with open(os.path.join(tmpdir, "test_pickle.pickle"), "rb") as f:
                backend_1 = pickle.load(f)
        assert self.model.get_ignore_status_global() == backend_1.get_ignore_status_global()
        assert self.model.get_synch_status_both_side() == backend_1.get_synch_status_both_side()
    
    def test_copy(self):
        """test copy preserve the flags"""
        backend_1 = self.model.copy()
        assert self.model.get_ignore_status_global() == backend_1.get_ignore_status_global()
        assert self.model.get_synch_status_both_side() == backend_1.get_synch_status_both_side()
        
    def test_gridmodel_line_global(self):
        """test disconnecting globally the line has the correct effect"""
        el_id = 0
        self.model.deactivate_powerline(el_id)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        else:
            assert not self.model.get_lines()[el_id].connected_global
        assert not self.model.get_lines()[el_id].connected1
        assert not self.model.get_lines()[el_id].connected2
        
        self.model.reactivate_powerline(el_id)
        assert self.model.get_lines()[el_id].connected_global
        assert self.model.get_lines()[el_id].connected1
        assert self.model.get_lines()[el_id].connected2
        
    def test_gridmodel_line_side1_topo(self):
        """test disconnecting the line side1 when updated from topology"""
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id]
        change[el_tp] = True
        new_values[el_tp] = -1
        
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert not self.model.get_lines()[el_id].connected1

        if self.synch_status_both_side():
            if not self.ignore_status_global():
                assert not self.model.get_lines()[el_id].connected_global
            assert not self.model.get_lines()[el_id].connected2
        else:
            assert self.model.get_lines()[el_id].connected2

        line = self.model.get_lines()[el_id]
        if self.synch_status_both_side():
            # both sides disconnected: all effective coeffs must be zero
            assert np.isclose(line.yac_eff_11, 0.)
            assert np.isclose(line.yac_eff_12, 0.)
            assert np.isclose(line.yac_eff_21, 0.)
            assert np.isclose(line.yac_eff_22, 0.)
        else:
            # side 1 open, side 2 connected: Kron-reduce the open end out
            assert np.isclose(line.yac_eff_11, 0.)
            assert np.isclose(line.yac_eff_12, 0.)
            assert np.isclose(line.yac_eff_21, 0.)
            assert np.isclose(line.yac_eff_22,
                              line.yac_22 - line.yac_21 * line.yac_12 / line.yac_11)

        # now reconnects it
        new_values[el_tp] = 1
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert self.model.get_lines()[el_id].connected1

        if self.synch_status_both_side():
            assert self.model.get_lines()[el_id].connected_global
            assert self.model.get_lines()[el_id].connected2
        else:
            assert self.model.get_lines()[el_id].connected2

        line = self.model.get_lines()[el_id]
        # both sides connected: effective coeffs must equal raw coeffs
        assert np.isclose(line.yac_eff_11, line.yac_11)
        assert np.isclose(line.yac_eff_12, line.yac_12)
        assert np.isclose(line.yac_eff_21, line.yac_21)
        assert np.isclose(line.yac_eff_22, line.yac_22)

    def test_gridmodel_line_side2_topo(self):
        """test disconnecting the line side2 when updated from topology"""
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_ex_pos_topo_vect[el_id]
        change[el_tp] = True
        new_values[el_tp] = -1
        
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert not self.model.get_lines()[el_id].connected2

        if self.synch_status_both_side():
            if not self.ignore_status_global():
                assert not self.model.get_lines()[el_id].connected_global
            assert not self.model.get_lines()[el_id].connected1
        else:
            assert self.model.get_lines()[el_id].connected1

        line = self.model.get_lines()[el_id]
        if self.synch_status_both_side():
            # both sides disconnected: all effective coeffs must be zero
            assert np.isclose(line.yac_eff_11, 0.)
            assert np.isclose(line.yac_eff_12, 0.)
            assert np.isclose(line.yac_eff_21, 0.)
            assert np.isclose(line.yac_eff_22, 0.)
        else:
            # side 1 connected, side 2 open: Kron-reduce the open end out
            assert np.isclose(line.yac_eff_11,
                              line.yac_11 - line.yac_21 * line.yac_12 / line.yac_22)
            assert np.isclose(line.yac_eff_12, 0.)
            assert np.isclose(line.yac_eff_21, 0.)
            assert np.isclose(line.yac_eff_22, 0.)

        # now reconnects it
        new_values[el_tp] = 1
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert self.model.get_lines()[el_id].connected2

        if self.synch_status_both_side():
            assert self.model.get_lines()[el_id].connected_global
            assert self.model.get_lines()[el_id].connected1
        else:
            assert self.model.get_lines()[el_id].connected1

        line = self.model.get_lines()[el_id]
        # both sides connected: effective coeffs must equal raw coeffs
        assert np.isclose(line.yac_eff_11, line.yac_11)
        assert np.isclose(line.yac_eff_12, line.yac_12)
        assert np.isclose(line.yac_eff_21, line.yac_21)
        assert np.isclose(line.yac_eff_22, line.yac_22)
            
    # def test_gridmodel_line_side1_changebus(self):
    #     """test disconnecting the line side1, when updated from change_bus"""
    #     # TODO lightsim2grid forbid this atm
    #     el_id = 0        
    #     self.model.change_bus_powerline_or(el_id, -1)
    #     if self.ignore_status_global():
    #         assert self.model.get_lines()[el_id].connected_global
    #     assert not self.model.get_lines()[el_id].connected1
        
    #     if self.synch_status_both_side():
    #         if not self.ignore_status_global():
    #             assert not self.model.get_lines()[el_id].connected_global
    #         assert not self.model.get_lines()[el_id].connected2
    #     else:
    #         assert self.model.get_lines()[el_id].connected2
        
    #     # now reconnects it
    #     self.model.change_bus_powerline_or(el_id, self.model.get_lines()[el_id].sub1_id)
    #     if self.ignore_status_global():
    #         assert self.model.get_lines()[el_id].connected_global
    #     assert self.model.get_lines()[el_id].connected1
        
    #     if self.synch_status_both_side():
    #         assert self.model.get_lines()[el_id].connected_global
    #         assert self.model.get_lines()[el_id].connected2
    #     else:
    #         assert self.model.get_lines()[el_id].connected2
            
    def test_gridmodel_line_both_sides_topo(self):
        """test disconnecting both sides of the line"""
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp1 = self.env.line_or_pos_topo_vect[el_id]
        el_tp2 = self.env.line_ex_pos_topo_vect[el_id]
        change[el_tp1] = True
        new_values[el_tp1] = -1
        change[el_tp2] = True
        new_values[el_tp2] = -1
        
        self.model.update_topo(change, new_values)
        assert not self.model.get_lines()[el_id].connected1
        assert not self.model.get_lines()[el_id].connected2
        if self.ignore_status_global():
            # flag not modified
            assert self.model.get_lines()[el_id].connected_global
        else:
            # automatic disconnection: both sides are disconnected
            assert not self.model.get_lines()[el_id].connected_global

        line = self.model.get_lines()[el_id]
        # both sides disconnected: all effective coeffs must be zero
        assert np.isclose(line.yac_eff_11, 0.)
        assert np.isclose(line.yac_eff_12, 0.)
        assert np.isclose(line.yac_eff_21, 0.)
        assert np.isclose(line.yac_eff_22, 0.)

        # now reconnects it
        new_values[el_tp1] = 1
        new_values[el_tp2] = 1
        self.model.update_topo(change, new_values)
        assert self.model.get_lines()[el_id].connected1
        assert self.model.get_lines()[el_id].connected2
        # both sides are reconnected, so this should reconnect this automatically
        # if it was disconnected
        assert self.model.get_lines()[el_id].connected_global

        line = self.model.get_lines()[el_id]
        # both sides connected: effective coeffs must equal raw coeffs
        assert np.isclose(line.yac_eff_11, line.yac_11)
        assert np.isclose(line.yac_eff_12, line.yac_12)
        assert np.isclose(line.yac_eff_21, line.yac_21)
        assert np.isclose(line.yac_eff_22, line.yac_22)
        
    
class GridModelDiscoOneSideTF(BaseDiscoOneSide, unittest.TestCase):
    pass


class GridModelDiscoOneSideTT(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return True
    def ignore_status_global(self):
        return True
    
    
class GridModelDiscoOneSideFT(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return False
    def ignore_status_global(self):
        return True
    
    
class GridModelDiscoOneSideFF(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return False
    def ignore_status_global(self):
        return False
    
    
class TestPFOk(unittest.TestCase):    
    """test powerflows are correctly working when powerlines are disconnected"""
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
        if grid_nm == "ieee118":
            self.gen_slack_id = 29
        elif grid_nm == "ieee300":
            self.gen_slack_id = "B7049-G"
            self.gen_slack_id = 55
        else:
            self.gen_slack_id = 0
            
        # not really used as the grid is loaded elsewhere
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.path = os.path.join(dir_path, "case_14_iidm")
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(
                "blank", 
                test=True,
                grid_path=self.path,
                n_busbar=1,
                _add_to_name=f"{type(self).__name__}_{grid_nm}",
                chronics_class=ChangeNothing,
                backend = LightSimBackend(
                    loader_method="pypowsybl",
                    gen_slack_id=self.gen_slack_id,
                    loader_kwargs={
                        "grid": self.net_datamodel,
                        "use_buses_for_sub": True,
                        "sort_index": False,
                        "use_grid2op_default_names": False,
                    }))
        self.model = self.env.backend._grid
        self.model.set_ignore_status_global(True)
        self.model.set_synch_status_both_side(False)
    
    def assert_equal(self, tmp, ref, error="", tol=None):
        assert np.all(tmp.shape == ref.shape), f"vector does not have the same shape for {error}"
        if tol is None:
            tol = self.tol_test
        assert np.max(np.abs(tmp - ref)) <= tol, f"{error}: {np.max(np.abs(tmp - ref))} > {tol}"
        assert np.mean(np.abs(tmp - ref)) <= tol, f"{error}: {np.mean(np.abs(tmp - ref))} > {tol}"
        
    def check_res(self, Vfinal, net, is_dc):
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
        if not is_dc:
            self.assert_equal(vm_pu_ls, vm_pu_pypow, "error for voltage magnitude", self.tol_V)
    
        # check lines
        l_is = self.net_ref.get_lines()["connected1"]
        por, qor, vor, aor = self.model.get_line_res1()
        self.assert_equal(por[l_is], net.get_lines()["p1"].values[l_is], "error for p_from")
        if not is_dc:
            self.assert_equal(qor[l_is], net.get_lines()["q1"].values[l_is], "error for q_from")
            self.assert_equal(aor[l_is], net.get_lines()["i1"].values[l_is] * 1e-3, "error for i_from_ka")
            vor_pp = net.get_buses().loc[net.get_lines().loc[l_is, "bus1_id"].values]["v_mag"].values
            self.assert_equal(vor[l_is], vor_pp, "error for vor_pp")

        # check trafo
        f_is = self.net_ref.get_2_windings_transformers()["connected2"]
        plv, qlv, vlv, alv = self.model.get_trafo_res2()
        self.assert_equal(plv[f_is], net.get_2_windings_transformers()["p2"].values[f_is], "error for p_lv_mw")
        if not is_dc:
            self.assert_equal(qlv[f_is], net.get_2_windings_transformers()["q2"].values[f_is], "error for q_lv_mvar")
            self.assert_equal(alv[f_is], net.get_2_windings_transformers()["i2"].values[f_is] * 1e-3, "error for i_lv_ka")
            vlv_pp = net.get_buses().loc[net.get_2_windings_transformers().loc[f_is, "bus2_id"].values]["v_mag"].values
            self.assert_equal(vlv[f_is], vlv_pp, "error for vlv_pp")

        # check loads
        l_is = self.net_ref.get_loads()["connected"].values
        load_p, load_q, load_v = self.model.get_loads_res()
        self.assert_equal(load_p[l_is], net.get_loads()["p"].values[l_is], "error for load p_mw")
        if not is_dc:
            self.assert_equal(load_q[l_is], net.get_loads()["q"].values[l_is], "error for load q_mvar")

        # check generators
        g_is = self.net_ref.get_generators()["connected"]
        prod_p, prod_q, prod_v = self.model.get_gen_res()
        
        self.assert_equal(prod_p[g_is], -net.get_generators()["p"].values[g_is], "error for gen p_mw")     
        if not is_dc:
            self.assert_equal(prod_q[g_is], -net.get_generators()["q"].values[g_is], "error for gen q_mvar")
            v_gen_pp = net.get_buses().loc[net.get_generators().loc[g_is, "bus_id"].values]["v_mag"].values
            self.assert_equal(prod_v[g_is], v_gen_pp, "error for prod_v")
        
        # check shunts
        s_is = self.net_ref.get_shunt_compensators()["connected"]
        shunt_p, shunt_q, shunt_v = self.model.get_shunts_res()
        self.assert_equal(shunt_p[s_is], net.get_shunt_compensators()["p"].values[s_is], "error for shunt p_mw")
        if not is_dc:
            self.assert_equal(shunt_q[s_is], net.get_shunt_compensators()["q"].values[s_is], "error for shunt q_mvar")

    def run_ref_pf(self, net, is_dc=False):
        if is_dc:
            params = get_pypowsybl_parameters(self.slack_vl_id_pypow)
            res_pypow = pp_lf.run_dc(net, params)
        else: 
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
    
    def test_gridmodel_line_side1_ac(self, is_dc=False):
        """test disconnecting the line side1"""
        # update lightsim2grid
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id]
        change[el_tp] = True
        new_values[el_tp] = -1
        self.model.update_topo(change, new_values)
        
        # sanity check
        # by setting the side 1 alone
        # on busbar 2
        model2 = init_from_pypowsybl(
            self.net_datamodel,
            slack_bus_id=self.ls_slack,
            buses_for_sub=True,
            sort_index=False,
            n_busbar_per_sub=2)
        model2.change_bus1_powerline(el_id, model2.get_lines()[el_id].sub1_id + len(model2.get_substations()))
        V0 = np.full(model2.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        if is_dc:
            Vres = model2.dc_pf(V0, 1, self.tol)
        else:
            Vres = model2.ac_pf(V0, self.max_it, self.tol)
            
        # update pypowsybl
        self.net_ref.update_lines(
            id=self.model.get_lines()[el_id].name,
            connected1=False,
        )
        
        # compute powerflows
        self.run_ref_pf(self.net_ref, is_dc)
        V0 = np.full(self.model.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        if is_dc:
            Vfinal = self.model.dc_pf(V0, 1, self.tol)
        else:
            Vfinal = self.model.ac_pf(V0, self.max_it, self.tol)
        
        
        # sanity check (for Vres, should match)
        model_init = self.model
        self.model = model2
        self.check_res(Vres[:118], self.net_ref, is_dc)
        
        # compare Vres[:118] and Vfinal
        self.model = model_init
        assert np.abs(Vres[:self.model.total_bus()] - Vfinal).max() <= self.tol, f"{np.abs(Vres[:self.model.total_bus()] - Vfinal).max()}"
        self.check_res(Vfinal, self.net_ref, is_dc)

    def test_gridmodel_line_side1_dc(self):
        self.test_gridmodel_line_side1_ac(is_dc=True)
        
    def test_gridmodel_trafo_side1_ac(self, is_dc=False):
        """test disconnecting the trafo side1, with ratio =/= 1"""
        self.tol_V = 1.8e-5  # TODO a bit large...
        self.tol_test = 9.9e-3  # TODO a bit large...
         
        # update lightsim2grid
        el_id = 0
        assert abs(self.model.get_trafos()[el_id].ratio - 1.) > 0.01
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id + len(self.model.get_lines())]
        change[el_tp] = True
        new_values[el_tp] = -1
        self.model.update_topo(change, new_values)
        
        # sanity check
        # by setting the side 1 alone
        # on busbar 2
        model2 = init_from_pypowsybl(
            self.net_datamodel,
            slack_bus_id=self.ls_slack,
            buses_for_sub=True,
            sort_index=False,
            n_busbar_per_sub=2)
        model2.change_bus1_trafo(el_id, model2.get_trafos()[el_id].sub1_id + len(model2.get_substations()))
        V0 = np.full(model2.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        if is_dc:
            Vres = model2.dc_pf(V0, 1, self.tol)
        else:
            Vres = model2.ac_pf(V0, self.max_it, self.tol)
            
        # update pypowsybl
        self.net_ref.update_2_windings_transformers(
            id=self.model.get_trafos()[el_id].name,
            connected1=False,
        )
        
        # compute powerflows
        self.run_ref_pf(self.net_ref, is_dc)
        V0 = np.full(self.model.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        if is_dc:
            Vfinal = self.model.dc_pf(V0, 1, self.tol)
        else:
            Vfinal = self.model.ac_pf(V0, self.max_it, self.tol)
        
        
        # sanity check (for Vres, should match)
        model_init = self.model
        self.model = model2
        self.check_res(Vres[:118], self.net_ref, is_dc)
        
        # compare Vres[:118] and Vfinal
        self.model = model_init
        assert np.abs(Vres[:self.model.total_bus()] - Vfinal).max() <= self.tol, f"{np.abs(Vres[:self.model.total_bus()] - Vfinal).max()}"
        self.check_res(Vfinal, self.net_ref, is_dc)
    
    def test_gridmodel_trafo_side1_dc(self):
        self.test_gridmodel_trafo_side1_ac(is_dc=True)
        
    def test_gridmodel_trafo_side1_alpha_ac(self, is_dc=False):
        """test disconnecting the trafo side1, with phase shift != 0."""
        self.setUp("ieee300")  # the smaller grid does not have phase shifters
        
        # self.tol_V = 1.8e-5  # TODO a bit large...
        # self.tol_test = 9.9e-3  # TODO a bit large...
         
        # update lightsim2grid
        el_id = 85
        assert abs(self.model.get_trafos()[el_id].shift_rad) > 0.01
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id + len(self.model.get_lines())]
        change[el_tp] = True
        new_values[el_tp] = -1
        self.model.update_topo(change, new_values)
        
        # sanity check
        # by setting the side 1 alone
        # on busbar 2
        # model2 = init_from_pypowsybl(
        #     self.net_datamodel,
        #     slack_bus_id=self.ls_slack,
        #     buses_for_sub=True,
        #     sort_index=False,
        #     n_busbar_per_sub=2)
        # model2.change_bus1_trafo(
        #     el_id,
        #     model2.get_trafos()[el_id].sub1_id + len(model2.get_substations()))
        # V0 = np.full(model2.total_bus(),
        #              fill_value=1.0,
        #              dtype=complex)
        # if is_dc:
        #     Vres = model2.dc_pf(V0, 1, self.tol)
        # else:
        #     Vres = model2.ac_pf(V0, self.max_it, self.tol)
            
        # update pypowsybl
        self.net_ref.update_2_windings_transformers(
            id=self.model.get_trafos()[el_id].name,
            connected1=False,
        )
        
        # compute powerflows
        self.run_ref_pf(self.net_ref, is_dc)
        V0 = np.full(self.model.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        if is_dc:
            Vfinal = self.model.dc_pf(V0, 1, self.tol)
        else:
            Vfinal = self.model.ac_pf(V0, self.max_it, self.tol)
        
        # sanity check (for Vres, should match)
        # model_init = self.model
        # self.model = model2
        # self.check_res(Vres[:len(model2.get_substations())], self.net_ref, is_dc)
        
        # compare Vres[:118] and Vfinal
        # self.model = model_init
        # assert np.abs(Vres[:self.model.total_bus()] - Vfinal).max() <= self.tol, f"{np.abs(Vres[:self.model.total_bus()] - Vfinal).max()}"
        self.check_res(Vfinal, self.net_ref, is_dc)
    
    def test_gridmodel_trafo_side1_alpha_dc(self):
        self.test_gridmodel_trafo_side1_alpha_ac(is_dc=True)


class TestImportHalfOpen(unittest.TestCase):
    """`init_from_pypowsybl(keep_half_open_lines=True)` models a branch connected on a
    single terminal as half-open (energized side kept, open end Kron-reduced), instead
    of fully deactivating it. Covers the static import path (the `TestPFOk` class above
    only exercises the grid2op `update_topo` path)."""

    LID = "L6-7-1"   # a line: side 1 opened
    TID = "T8-5-1"   # a transformer: side 2 opened
    GEN_SLACK_ID = 29

    def _open_sides(self, n):
        n.update_lines(id=self.LID, connected1=False)
        n.update_2_windings_transformers(id=self.TID, connected2=False)

    def _line(self, model):
        li = list(pp_network.create_ieee118().get_lines().index).index(self.LID)
        return model.get_lines()[li]

    def _trafo(self, model):
        ti = list(pp_network.create_ieee118().get_2_windings_transformers().index).index(self.TID)
        return model.get_trafos()[ti]

    def test_import_creates_half_open(self):
        n = pp_network.create_ieee118()
        self._open_sides(n)
        model = init_from_pypowsybl(n, gen_slack_id=self.GEN_SLACK_ID,
                                    sort_index=False, keep_half_open_lines=True)
        # line: side 1 open, side 2 energized, branch still globally connected
        line = self._line(model)
        assert not line.connected1
        assert line.connected2
        assert line.connected_global
        # open end Kron-reduced out; only the side-2 self admittance survives
        assert np.isclose(line.yac_eff_11, 0.)
        assert np.isclose(line.yac_eff_12, 0.)
        assert np.isclose(line.yac_eff_21, 0.)
        assert np.isclose(line.yac_eff_22,
                          line.yac_22 - line.yac_21 * line.yac_12 / line.yac_11)
        # transformer: side 2 open, side 1 energized
        tr = self._trafo(model)
        assert tr.connected1
        assert not tr.connected2
        assert tr.connected_global
        assert np.isclose(tr.yac_eff_22, 0.)
        assert np.isclose(tr.yac_eff_11,
                          tr.yac_11 - tr.yac_21 * tr.yac_12 / tr.yac_22)

    def test_import_keep_false_deactivates(self):
        """Default (keep_half_open_lines=False): a one-side-open branch is fully off."""
        n = pp_network.create_ieee118()
        self._open_sides(n)
        model = init_from_pypowsybl(n, gen_slack_id=self.GEN_SLACK_ID, sort_index=False)
        line = self._line(model)
        assert not line.connected1
        assert not line.connected2
        assert not line.connected_global
        tr = self._trafo(model)
        assert not tr.connected_global

    def _check_vs_olf(self, is_dc):
        slack_vl, _ = get_same_slack("ieee118")
        # reference: pypowsybl with the same one-sided outages, slack pinned by name
        nref = pp_network.create_ieee118()
        self._open_sides(nref)
        params = get_pypowsybl_parameters(slack_vl)
        res = pp_lf.run_dc(nref, params) if is_dc else pp_lf.run_ac(nref, params)
        slack_abs = res[0].slack_bus_results[0].active_power_mismatch
        gname = nref.get_generators().index[self.GEN_SLACK_ID]
        nref.update_generators(
            id=gname, p=nref.get_generators().iloc[self.GEN_SLACK_ID]["p"] - slack_abs)

        n = pp_network.create_ieee118()
        self._open_sides(n)
        model = init_from_pypowsybl(n, gen_slack_id=self.GEN_SLACK_ID,
                                    sort_index=False, keep_half_open_lines=True)
        nb = model.total_bus()
        V0 = np.full(nb, 1.0, dtype=complex)
        V = model.dc_pf(V0, 1, 1e-8) if is_dc else model.ac_pf(V0, 10, 1e-8)
        assert V.shape[0] > 0, "powerflow diverged"

        buses = nref.get_buses()
        vls = nref.get_voltage_levels()
        vm_ref = (buses["v_mag"] /
                  vls.loc[buses["voltage_level_id"], "nominal_v"].values).to_numpy()
        va_ref = np.deg2rad(buses["v_angle"].values)
        # tolerances follow TestPFOk (the trafo half-open is the looser one)
        assert np.abs(np.angle(V[:len(va_ref)]) - va_ref).max() <= 3e-5
        if not is_dc:
            assert np.abs(np.abs(V[:len(vm_ref)]) - vm_ref).max() <= 3e-5

    def test_import_half_open_ac_matches_olf(self):
        self._check_vs_olf(is_dc=False)

    def test_import_half_open_dc_matches_olf(self):
        self._check_vs_olf(is_dc=True)

# TODO trafo with alpha (phase shift)
# TODO FDPF powerflow too
