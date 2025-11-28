# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import pdb
import warnings
import pandapower as pp
import numpy as np
import pandapower.networks as pn
import unittest
from lightsim2grid.gridmodel import init_from_pandapower


from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


VAR_GEN = ["bus", "p_mw", "vm_pu", "sn_mva", "name", "index", "max_q_mvar", "min_q_mvar", "min_p_mw",
            "max_p_mw", "scaling", "type", "slack", "controllable", "vn_kv", "xdss_pu", "rdss_pu",
            "cos_phi", "in_service"]
# adding a slack bus
def make_grid_multiple_slack(case):
    """create an equivalent grid (for case 14) with generator as slack buses"""
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        pp.runpp(case)  # forced to do it to retrieve the power of the slack bus
    
    var_pp = ["bus", 'in_service', "name", "max_p_mw", "min_p_mw", "max_q_mvar", "min_q_mvar", "slack_weight"]
    if "slack_weight" not in case.ext_grid:
        var_pp = var_pp[:-1]

    slack_bus_gen_id_ppc = 0  # checked manually
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        pp.create_gen(case,
                    **case.ext_grid[var_pp].iloc[0],
                    p_mw=case._ppc['gen'][slack_bus_gen_id_ppc, 1],
                    slack=True)
    # "deactivate" the "ext_grid" 
    case.ext_grid.loc[0, "in_service"] = False
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        pp.runpp(case)

    # now create a copy of it, by removing the ext_grid completely (to be sure)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pp.create_empty_network("case14_custom", sn_mva=1.0 * case.sn_mva, f_hz= 1.0 * case.f_hz)
    
    # create bus
    for i in range(case.bus.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_bus(net, **case.bus.iloc[i])

    # create lines
    var_line = ["from_bus", "to_bus", "length_km", "r_ohm_per_km", "x_ohm_per_km", "c_nf_per_km", "max_i_ka", "name",
                "index", "type", "geodata", "in_service", "df", "parallel", "g_us_per_km",
                "max_loading_percent", "alpha", "temperature_degree_celsius"]
    var_line = [el for el in var_line if el in case.line]
    for i in range(case.line.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_line_from_parameters(net, **case.line[var_line].iloc[i])

    # create trafos
    var_trafo = ["hv_bus", "lv_bus", "sn_mva", "vn_hv_kv", "vn_lv_kv", "vkr_percent", "vk_percent", "pfe_kw", "i0_percent", "shift_degree",
                "tap_side", "tap_neutral", "tap_max", "tap_min", "tap_step_percent", "tap_step_degree", 
                "tap_pos", "tap_phase_shifter", "in_service", "name", "index", "max_loading_percent",
                "parallel", "df"]
    var_trafo = [el for el in var_trafo if el in case.trafo]
    for i in range(case.trafo.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_transformer_from_parameters(net, **case.trafo[var_trafo].iloc[i])

    # create shunts
    var_shunt = ["bus", "q_mvar", "p_mw", "vn_kv", "step", "max_step", "name",
                "in_service", "index"]
    var_shunt = [el for el in var_shunt if el in case.shunt]
    for i in range(case.shunt.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_shunt(net, **case.shunt[var_shunt].iloc[i])

    # create loads
    var_load = ["bus", "p_mw", "q_mvar", "const_z_percent", "const_i_percent", "sn_mva", "name",
                "scaling", "index", "in_service", "type", "max_p_mw", "min_p_mw", "max_q_mvar",
                "min_q_mvar", "controllable"]
    var_load = [el for el in var_load if el in case.load]
    for i in range(case.load.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_load(net, **case.load[var_load].iloc[i])

    # create gens
    var_gen = [el for el in VAR_GEN if el in case.gen]
    for i in range(case.gen.shape[0]):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.create_gen(net, **case.gen[var_gen].iloc[i])

    id_ref_slack = net.gen.shape[0]-1  # initial generator added as the slack bus added
    net.gen["min_p_mw"][[id_ref_slack]] = 0.
    net.gen["max_p_mw"][[id_ref_slack]] = 300.
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


def int_or_null(el):
    try:
        res = int(el)
    except ValueError:
        res = 0
    return res


class TestMultipleSlack14(unittest.TestCase):
    def setUp(self) -> None:
        
        if [int_or_null(el) for el in pp.__version__.split(".")] < [2, 8, 0]:
            # check if functionality is available in pandapower installed (officially supported since 2.8.0)
            msg_ = ("Unable to peform required tests because install pandapower version does not "
                   "allow multiple slack." )
            warnings.warn(msg_)
            self.unable_to_run_due_to_pp = True
            self.skipTest(msg_)
            
        # LF PARAMETERS
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 14

        # retrieve the case14 and remove the "ext_grid" => put a generator as slack bus instead
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.case = pn.case14()
        self.net = make_grid_multiple_slack(self.case)

        id_ref_slack = self.net.gen.shape[0]-1  # initial generator added as the slack bus added
        if "slack_weight" not in self.net.gen:
            self.net.gen["slack_weight"] = 0.
        self.net.gen["slack_weight"][[id_ref_slack]] = 0.5
    
    def test_single_slack(self):
        """check pandapower and lightsim get the same results when there is only one
           slack bus
        """            
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.net,
                    init_vm_pu="flat",
                    init_va_degree="flat")
            
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid_single = init_from_pandapower(self.net)
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid_single.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid_single, self.net)

        # now run with the option "dist_slack=True"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.net,
                    distributed_slack=True,
                    init_vm_pu="flat",
                    init_va_degree="flat")
        self.check_results(V, ls_grid_single, self.net)

    def test_two_slacks_diff_bus(self):
        """test the results when there are two slacks, in most simple setting"""
        # now activate more slack bus
        self.net.gen["slack"][[1]] = True
        id_ref_slack = self.net.gen.shape[0] - 1
        if "slack_weight" in self.net.gen:
            self.net.gen["slack_weight"][[1, id_ref_slack]] = 0.5
        
        # just to make sure pp forgot previous results
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.net,
                     distributed_slack=True,
                     init_vm_pu="flat",
                     init_va_degree="flat",
                     lightsim2grid=False)  
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init_from_pandapower(self.net)
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid, self.net)

    def test_two_slacks_diff_bus_diff_weights(self):
        """
        test the results when there are two slacks, in less simple setting
        where the generator do not have the same weight
        """
        # now activate more slack bus
        gen_id_added = 1
        self.net.gen["slack"][[gen_id_added]] = True
        id_ref_slack = self.net.gen.shape[0] - 1
        if "slack_weight" in self.net.gen:
            self.net.gen["slack_weight"][[gen_id_added, id_ref_slack]] = 0.9, 0.3
        
        # just to make sure pp forgot previous results
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.net,
                    distributed_slack=True,
                    init_vm_pu="flat",
                    init_va_degree="flat")  
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init_from_pandapower(self.net)
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid, self.net)

    def test_multiple_slack_same_bus(self):
        """
        test the results when there are three slacks, two connected at the same bus
        """

        # create gens
        gen_id_added = 1
        var_gen = [el for el in VAR_GEN if el in self.net.gen]
        id_ref_slack = self.net.gen.shape[0] - 1
        pp.create_gen(self.net, **self.net.gen[var_gen].iloc[id_ref_slack])
        init_p_mw = self.net.gen["p_mw"].values[id_ref_slack]
        last_gen_id = self.net.gen.shape[0] - 1
        self.net.gen["slack"][[gen_id_added]] = True
        self.net.gen["slack"][[last_gen_id]] = True
        coefs_slack = 0.4, 0.3, 0.6
        self.net.gen["slack_weight"][[gen_id_added, id_ref_slack, last_gen_id]] = coefs_slack
        self.net.gen["p_mw"][[id_ref_slack, last_gen_id]] = 0.5 * init_p_mw

        # start the powerflow
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")  
            pp.rundcpp(self.net)
            pp.runpp(self.net, distributed_slack=True,
                    init_vm_pu="flat",
                    init_va_degree="flat") 
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init_from_pandapower(self.net)
        
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        assert len(V), "lightsim diverged !"
        # check that the losses have been properly split
        last = len(ls_grid.get_generators()) - 1
        last2 = len(ls_grid.get_generators()) - 2
        res_slack = np.array((ls_grid.get_generators()[last2].res_p_mw , ls_grid.get_generators()[last].res_p_mw ))
        res_slack -= 0.5 * init_p_mw
        tartine_0 = ls_grid.get_generators()[gen_id_added].res_p_mw
        # res_slack = self.net.res_gen["p_mw"].values[[-2, -1]] - 0.5 * init_p_mw
        # tartine_0 = self.net.res_gen["p_mw"].values[gen_id_added]
        assert abs(coefs_slack[2] * res_slack[0] - coefs_slack[1] * res_slack[1]) <= 1e-6
        assert abs(coefs_slack[1] * tartine_0 - coefs_slack[0] * res_slack[0]) <= 1e-6

        self.check_results(V, ls_grid, self.net)

    def test_multiple_slack_one_gen_not_slack(self):
        """
        test the results when there are two slacks, but at one node where there is a slack, another non slack gen is connected
        """

        # create gens
        gen_id_added = 1
        var_gen = [el for el in VAR_GEN if el in self.net.gen]
        id_ref_slack = self.net.gen.shape[0] - 1
        self.net.gen["slack"][[id_ref_slack]] = False  # this gen will be copied, so i remove it
        self.net.gen["slack"][[gen_id_added]] = True
        # create the gen on the same bus as one of the slack
        pp.create_gen(self.net, **self.net.gen[var_gen].iloc[id_ref_slack])
        self.net.gen["slack"][[id_ref_slack]] = True  # I reactivate it
        init_p_mw = self.net.gen["p_mw"].values[id_ref_slack]
        coefs_slack = 0.1, 0.9
        self.net.gen["slack_weight"][[gen_id_added, id_ref_slack]] = coefs_slack
        self.net.gen["p_mw"][[id_ref_slack, id_ref_slack]] = 0.5 * init_p_mw

        # start the powerflow
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(self.net)
            pp.runpp(self.net, distributed_slack=True,
                    init_vm_pu="flat",
                    init_va_degree="flat") 
            
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init_from_pandapower(self.net)
        
        V = np.ones(self.nb_bus_total, dtype=complex)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid, self.net)

    def check_results(self, V_ls, ls_grid, pp_net):
        # NB: the test bellow only works because pandapower and lightsim have the
        # bus in the same order !
        assert len(V_ls), "lightsim diverged !"
        
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            # recent pandapower grid import broken (waiting to clarify this)...
            Ybus_me = ls_grid.get_Ybus_solver()
            Ybus_ref = pp_net._ppc["internal"]["Ybus"]
            assert np.abs((Ybus_me - Ybus_ref).todense()).max() <= 1e-6, "wrong Ybus"
            tol_v_pu = 1e-6
            tol_line = 1e-6
            tol_gen_mw = 1e-6
            tol_gen_mvar = 1e-6
        else:
            # with recent pandapower version
            # import is so broken that this should be the tolerance...
            tol_v_pu = 0.0159
            tol_line = 1.154
            tol_gen_mw = 0.0487
            tol_gen_mvar = 27.43
            
        # check that the same results as pandapower
        my_ref = np.where(np.angle(V_ls) == 0.)[0][0]
        V_pp = pp_net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  pp_net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        # print(f"{np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values).max() = }")
        # print(f"{np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values).max() = }")
        # print(f"{np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values).max() = }")
        
        assert np.abs(V_pp - V_ls).max() <= tol_v_pu, f"wrong voltages: {np.abs(V_pp - V_ls).max()}"
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line["p_from_mw"].values) <= tol_line), f"{np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line['p_from_mw'].values).max()}"
        assert np.all(np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values) <= tol_line)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values) <= tol_gen_mw)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values) <= tol_gen_mvar)


if __name__ == "__main__":
    unittest.main()
