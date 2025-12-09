# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import warnings
import pandapower as pp
import numpy as np
import pandapower.networks as pn
import grid2op
import unittest
from lightsim2grid import LightSimBackend
from lightsim2grid.gridmodel import init_from_pandapower
from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance as GridStateFromFile


from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


VAR_GEN = ["bus", "p_mw", "vm_pu", "sn_mva", "name", "index", "max_q_mvar", "min_q_mvar", "min_p_mw",
            "max_p_mw", "scaling", "type", "slack", "controllable", "vn_kv", "xdss_pu", "rdss_pu",
            "cos_phi", "in_service"]
            
# adding a slack bus
def make_grid_multiple_slack(case):
    """create an equivalent grid (for case 14) with generator as slack buses"""
    pp.runpp(case)  # forced to do it to retrieve the power of the slack bus
    var_pp = ["bus", 'in_service', "name", "max_p_mw", "min_p_mw", "max_q_mvar", "min_q_mvar", "slack_weight"]
    if "slack_weight" not in case.ext_grid:
        var_pp = var_pp[:-1]

    slack_bus_gen_id_ppc = 0  # checked manually
    pp.create_gen(case,
                **case.ext_grid[var_pp].iloc[0],
                p_mw=case._ppc['gen'][slack_bus_gen_id_ppc, 1],
                slack=True)
    # "deactivate" the "ext_grid" 
    case.ext_grid.loc[0, "in_service"] = False
    pp.runpp(case)

    # now create a copy of it, by removing the ext_grid completely (to be sure)
    net = pp.create_empty_network("case118_custom", sn_mva=1.0 * case.sn_mva, f_hz= 1.0 * case.f_hz)
    # create bus
    for i in range(case.bus.shape[0]):
        pp.create_bus(net, **case.bus.iloc[i])

    # create lines
    var_line = ["from_bus", "to_bus", "length_km", "r_ohm_per_km", "x_ohm_per_km", "c_nf_per_km", "max_i_ka", "name",
                "index", "type", "geodata", "in_service", "df", "parallel", "g_us_per_km",
                "max_loading_percent", "alpha", "temperature_degree_celsius"]
    var_line = [el for el in var_line if el in case.line]
    for i in range(case.line.shape[0]):
        pp.create_line_from_parameters(net, **case.line[var_line].iloc[i])

    # create trafos
    if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
        var_trafo = ["hv_bus", "lv_bus", "sn_mva", "vn_hv_kv", "vn_lv_kv", "vkr_percent", "vk_percent", "pfe_kw", "i0_percent", "shift_degree",
                    "tap_side", "tap_neutral", "tap_max", "tap_min", "tap_step_percent", "tap_step_degree", 
                    "tap_pos", "tap_phase_shifter", "in_service", "name", "index", "max_loading_percent",
                    "parallel", "df"]
    else:
        var_trafo = ['name', 'std_type', 'hv_bus', 'lv_bus', 'sn_mva', 'vn_hv_kv',
        'vn_lv_kv', 'vk_percent', 'vkr_percent', 'pfe_kw', 'i0_percent',
        'shift_degree', 'tap_side', 'tap_neutral', 'tap_min', 'tap_max',
        'tap_step_percent', 'tap_step_degree', 'tap_pos', 'tap_changer_type',
        'id_characteristic_table', 'tap_dependency_table', 'parallel', 'df',
        'in_service', 'max_loading_percent', 'oltc']
    var_trafo = [el for el in var_trafo if el in case.trafo]
    for i in range(case.trafo.shape[0]):
        pp.create_transformer_from_parameters(net, **case.trafo[var_trafo].iloc[i])

    # create shunts
    var_shunt = ["bus", "q_mvar", "p_mw", "vn_kv", "step", "max_step", "name",
                "in_service", "index"]
    var_shunt = [el for el in var_shunt if el in case.shunt]
    for i in range(case.shunt.shape[0]):
        pp.create_shunt(net, **case.shunt[var_shunt].iloc[i])

    # create loads
    var_load = ["bus", "p_mw", "q_mvar", "const_z_percent", "const_i_percent", "sn_mva", "name",
                "scaling", "index", "in_service", "type", "max_p_mw", "min_p_mw", "max_q_mvar",
                "min_q_mvar", "controllable"]
    var_load = [el for el in var_load if el in case.load]
    for i in range(case.load.shape[0]):
        pp.create_load(net, **case.load[var_load].iloc[i])

    # create gens
    var_gen = [el for el in VAR_GEN if el in case.gen]
    for i in range(case.gen.shape[0]):
        pp.create_gen(net, **case.gen[var_gen].iloc[i])

    id_ref_slack = net.gen.shape[0]-1  # initial generator added as the slack bus added
    net.gen.loc[[id_ref_slack], "min_p_mw"] = 0.
    net.gen.loc[[id_ref_slack], "max_p_mw"] = 300.
    assert (net.gen["min_p_mw"] >= 0.).all()
    assert (net.gen["min_p_mw"] <= 400.).all()

    # we can test that is the same grid if we want to
    # pp.runpp(case, init="flat")
    # pp.runpp(net, init="flat")

    # assert ((case.res_gen - net.res_gen) <= 1e-6).all().all()
    # assert ((case.res_load - net.res_load) <= 1e-6).all().all()
    # assert ((case.res_shunt - net.res_shunt) <= 1e-6).all().all()
    # assert ((case.res_line - net.res_line) <= 1e-6).all().all()
    # assert ((case.res_trafo - net.res_trafo) <= 1e-6).all().all()
    # assert ((case.res_bus.iloc[:net.bus.shape[0]] - net.res_bus) <= 1e-6).all().all()

    return net

class TestMultipleL2RPN(unittest.TestCase):
    def setUp(self) -> None:
        # if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
        #     self.skipTest("Test not correct: pp changed the way it computed trafo params")
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 118

    def test_neurips_track2(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make("l2rpn_neurips_2020_track2",
                               test=True,
                               data_feeding_kwargs={"gridvalueClass": GridStateFromFile})  

        li_envs = list(env.keys())

        for el in li_envs:
            self.pp_net = env[el].backend._grid
            self.idx_slack = np.where(self.pp_net.gen["slack"].values)[0]
            if "slack_weight" not in self.pp_net.gen:
                warnings.warn("Unable to peform required tests because install pandapower version does not "
                              "allow multiple slack.")
                return
            self.pp_net.gen.loc[1, "slack_weight"] = 1.
            pp.rundcpp(self.pp_net)  # to forget the result
            pp.runpp(self.pp_net, distributed_slack=True, init_vm_pu="flat", init_va_degree="flat")
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                ls_grid = init_from_pandapower(self.pp_net, pp_orig_file="pandapower_v3")
            ls_grid.tell_solver_need_reset()
            V = np.ones(2 * self.nb_bus_total, dtype=complex)
            V = ls_grid.ac_pf(V, self.max_it, self.tol)
            self.check_results(V[:self.nb_bus_total], ls_grid, self.pp_net)

    def check_results(self, V_ls, ls_grid, pp_net):
        assert len(V_ls), "lightsim diverged !"
        my_ref = np.where(np.angle(V_ls) == 0.)[0][0]
        V_pp = pp_net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  pp_net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        V_pp = V_pp[:self.nb_bus_total]
        
        Ybus_pp = pp_net._ppc["internal"]["Ybus"]
        Ybus_ls = ls_grid.get_Ybus_solver()
        assert np.abs((Ybus_pp - Ybus_ls).toarray()).max() <= 1e-6, f"wrong Ybus {np.abs((Ybus_pp - Ybus_ls).toarray()).max()}"
        assert np.abs(V_pp - V_ls).max() <= 1e-6, f"wrong voltages: {np.abs(V_pp - V_ls).max()}"
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line["p_from_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_a_or_ka for el in ls_grid.get_lines()] - pp_net.res_line["i_from_ka"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values) <= 1e-6)
        

class Test118LightsimBackend(unittest.TestCase):
    def setUp(self) -> None:
        # if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
            # self.skipTest("Test not correct: pp changed the way it computed trafo params")
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 118

    def test_make_and_pf(self):
        env_name_input = "l2rpn_neurips_2020_track2"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env_ls = grid2op.make(env_name_input, backend=LightSimBackend(), test=True,
                                       data_feeding_kwargs={"gridvalueClass": GridStateFromFile})

class TestMultipleSlack118(unittest.TestCase):
    def setUp(self) -> None:
        # if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
            # self.skipTest("Test not correct: pp changed the way it computed trafo params")
        # LF PARAMETERS
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 118

        # retrieve the case14 and remove the "ext_grid" => put a generator as slack bus instead
        self.case = pn.case118()
        self.net = make_grid_multiple_slack(self.case)
        if "slack_weight" in self.net.gen:
            id_ref_slack = self.net.gen.shape[0]-1  # initial generator added as the slack bus added
            self.net.gen.loc[id_ref_slack, "slack_weight"] = 0.5

        self.unable_to_run_due_to_pp = False
        if "slack_weight" not in self.net.gen:
            msg_ = "Unable to peform required tests because install pandapower version does not " \
                   "allow multiple slack."
            warnings.warn(msg_)
            self.unable_to_run_due_to_pp = True
            self.skipTest(msg_)

    def test_single_slack(self):
        """check pandapower and lightsim get the same results when there is only one
           slack bus
        """
        pp.runpp(self.net,
                 init_vm_pu="flat",
                 init_va_degree="flat")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid_single = init_from_pandapower(self.net, pp_orig_file="pandapower_v3")
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid_single.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid_single, self.net)

        # now run with the option "dist_slack=True"
        pp.runpp(self.net,
                 distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat")
        self.check_results(V, ls_grid_single, self.net)

    def test_two_slacks_diff_bus(self):
        """test the results when there are two slacks, in most simple setting"""
        # now activate more slack bus
        self.net.gen.loc[1, "slack"] = True
        id_ref_slack = self.net.gen.shape[0] - 1
        if "slack_weight" in self.net.gen:
            self.net.gen.loc[[1, id_ref_slack], "slack_weight"] = 0.5
        
        # just to make sure pp forgot previous results
        pp.runpp(self.net,
                 distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat")  
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init_from_pandapower(self.net, pp_orig_file="pandapower_v3")
            
        self.net._ppc["internal"]["Ybus"]
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid, self.net)

    def check_results(self, V_ls, ls_grid, pp_net):
        assert len(V_ls), "lightsim diverged !"
        
        Ybus_pp = pp_net._ppc["internal"]["Ybus"]
        Ybus_ls = ls_grid.get_Ybus_solver()
        # (np.abs((Ybus_pp - Ybus_ls).toarray())> 3.).nonzero()
        assert np.abs((Ybus_pp - Ybus_ls).toarray()).max() <= 1e-6, f"wrong Ybus {np.abs((Ybus_pp - Ybus_ls).toarray()).max()}"
        
        my_ref = np.where(np.angle(V_ls) == 0.)[0][0]
        V_pp = pp_net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  pp_net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        assert np.abs(V_pp - V_ls).max() <= 1e-6, f"wrong voltages: {np.abs(V_pp - V_ls).max()}"
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line["p_from_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_a_or_ka for el in ls_grid.get_lines()] - pp_net.res_line["i_from_ka"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values) <= 1e-6)
        

if __name__ == "__main__":
    unittest.main()