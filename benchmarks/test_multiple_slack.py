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
import pdb
import copy
import tempfile
import unittest
from lightsim2grid.initGridModel import init

# perf
# net = pn.case300()
# ls_grid_single = init(net)
# ls_grid_single.deactivate_result_computation()
# V = np.ones(net.bus.shape[0], dtype=np.complex_)
# # Vdc = ls_grid_single.dc_pf(copy.deepcopy(V), max_it, tol)
# ls_grid_single.reactivate_result_computation()
# V = ls_grid_single.ac_pf(V, 10, 1e-8)
# import sys
# sys.exit()

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
    case.ext_grid["in_service"][0] = False
    pp.runpp(case)

    # now create a copy of it, by removing the ext_grid completely (to be sure)
    net = pp.create_empty_network("case14_custom", sn_mva=1.0 * case.sn_mva, f_hz= 1.0 * case.f_hz)
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
    var_trafo = ["hv_bus", "lv_bus", "sn_mva", "vn_hv_kv", "vn_lv_kv", "vkr_percent", "vk_percent", "pfe_kw", "i0_percent", "shift_degree",
                "tap_side", "tap_neutral", "tap_max", "tap_min", "tap_step_percent", "tap_step_degree", 
                "tap_pos", "tap_phase_shifter", "in_service", "name", "index", "max_loading_percent",
                "parallel", "df"]
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

class TestMultipleSlack(unittest.TestCase):
    def setUp(self) -> None:
        # LF PARAMETERS
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 14

        # retrieve the case14 and remove the "ext_grid" => put a generator as slack bus instead
        self.case = pn.case14()
        self.net = make_grid_multiple_slack(self.case)
        if "slack_weight" in self.net.gen:
            id_ref_slack = self.net.gen.shape[0]-1  # initial generator added as the slack bus added
            self.net.gen["slack_weight"][[id_ref_slack]] = 0.5
    
    def test_single_slack(self):
        """check pandapower and lightsim get the same results when there is only one
           slack bus
        """
        # TODO SLACK uncomment
        pp.runpp(self.net,
                 init_vm_pu="flat",
                 init_va_degree="flat")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid_single = init(self.net)
        V = np.ones(self.nb_bus_total, dtype=np.complex_)
        V = ls_grid_single.ac_pf(V, self.max_it, self.tol)
        J_me = ls_grid_single.get_J()
        J_pp = self.net._ppc["internal"]["J"]
        self.check_results(V, ls_grid_single, self.net)  # TODO SLACK uncomment

        # now run with the option "dist_slack=True"
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
        pp.runpp(self.net,
                 distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat")  
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init(self.net)
        V = np.ones(self.nb_bus_total, dtype=np.complex_)
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
        pp.runpp(self.net,
                 distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat")  
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init(self.net)
        V = np.ones(self.nb_bus_total, dtype=np.complex_)
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
        pp.rundcpp(self.net)
        pp.runpp(self.net, distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat") 
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init(self.net)
        
        V = np.ones(self.nb_bus_total, dtype=np.complex_)
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
        pp.rundcpp(self.net)
        pp.runpp(self.net, distributed_slack=True,
                 init_vm_pu="flat",
                 init_va_degree="flat") 
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            ls_grid = init(self.net)
        
        V = np.ones(self.nb_bus_total, dtype=np.complex_)
        V = ls_grid.ac_pf(V, self.max_it, self.tol)
        self.check_results(V, ls_grid, self.net)

    def check_results(self, V_ls, ls_grid, pp_net):
        # NB: the test bellow only works because pandapower and lightsim have the
        # bus in the same order !
        assert len(V_ls), "lightsim diverged !"
        Ybus_me = ls_grid.get_Ybus()
        Ybus_ref = pp_net._ppc["internal"]["Ybus"]
        assert np.abs((Ybus_me - Ybus_ref).todense()).max() <= 1e-6, "wrong Ybus"

        my_ref = np.where(np.angle(V_ls) == 0.)[0][0]
        V_pp = pp_net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  pp_net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        assert np.abs(V_pp - V_ls).max() <= 1e-6, "wrong voltages"
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line["p_from_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values) <= 1e-6)

if __name__ == "__main__":
    unittest.main()
    if False:
        if False:  # test pass now
            print()
            print()
            print()
            print()
            print("With only one slack (case 2) ")
            pp.runpp(case, init="flat")
            print("\n")
            print("Real testcase")
            pp.runpp(net, init="flat")

            assert ((case.res_gen - net.res_gen) <= 1e-6).all().all()
            assert ((case.res_load - net.res_load) <= 1e-6).all().all()
            assert ((case.res_shunt - net.res_shunt) <= 1e-6).all().all()
            assert ((case.res_line - net.res_line) <= 1e-6).all().all()
            assert ((case.res_trafo - net.res_trafo) <= 1e-6).all().all()
            assert ((case.res_bus.iloc[:net.bus.shape[0]] - net.res_bus) <= 1e-6).all().all()

            gen_without_single_slack = copy.deepcopy(net.res_gen)
            print(gen_without_single_slack["p_mw"])

            ls_grid_single = init(net)
            nb_bus_total = 14
            ls_grid_single.deactivate_result_computation()
            V = np.ones(nb_bus_total, dtype=np.complex_)
            # Vdc = ls_grid_single.dc_pf(copy.deepcopy(V), max_it, tol)
            ls_grid_single.reactivate_result_computation()
            V = ls_grid_single.ac_pf(V, max_it, tol)
            real_J = np.load("J_dist_slack_one_slack_firstIter.npy")  # from pandapower
            real_J_ls = np.load(file="J_ref_single_slack_ls_firstIter.npy")  # from ls single slack
            my_J = ls_grid_single.get_J()
            my_first_row = my_J[0].todense()
            ref_first_row = real_J[0]
            ref_first_row_ls = real_J_ls[0]
            Ybus_me = ls_grid_single.get_Ybus()
            Ybus_ref = net._ppc["internal"]["Ybus"]
            assert np.abs((Ybus_me - Ybus_ref).todense()).max() <= 1e-6, "wrong Ybus"
            # n_pv = 5
            # n_pq = 8
            pdb.set_trace()
            assert np.all(np.abs([el.res_p_or_mw for el in ls_grid_single.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)

        print()
        print()
        print()
        print()
        print("Run PF with distributed slack (but still one slack)")
        pp.runpp(net, init="flat", distributed_slack=True)
        # now activate more slack bus
        net.gen["slack"][[1]] = True
        if "slack_weight" in net.gen:
            net.gen["slack_weight"][[1, id_ref_slack]] = 0.5
        pp.rundcpp(net)

        if False:
            # run the powerflow without distributed_slack=True
            print("Multiple slack, but nothing change in runpp (case 3)")
            pp.runpp(net, init="flat")
            gen_without_dist_slack = copy.deepcopy(net.res_gen)
            print(gen_without_dist_slack["p_mw"])

        # run the powerflow with distributed_slack=True (same weights)
        print("Multiple slack, and distributed_slack=True (case 4)")
        pp.runpp(net, init="flat", distributed_slack=True)
        print(net.res_gen["p_mw"])

        ls_grid = init(net)
        nb_bus_total = 14
        # ls_grid.deactivate_result_computation()
        V = np.ones(nb_bus_total, dtype=np.complex_)  # * ls_grid.get_init_vm_pu()
        # Vdc = ls_grid.dc_pf(copy.deepcopy(V), max_it, tol)
        # ls_grid.reactivate_result_computation()
        V = ls_grid.ac_pf(V, max_it, tol)
        my_ref = np.where(np.angle(V) == 0.)[0][0]
        V_pp = net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        assert np.abs(V_pp - V).max() <= 1e-6
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - net.res_gen["q_mvar"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - net.res_gen["p_mw"].values) <= 1e-6)

        load_p = ls_grid.get_loads_res()[0]
        gen_p = ls_grid.get_gen_res()[0]
        p_or = ls_grid.get_lineor_res()[0]
        p_ex = ls_grid.get_lineex_res()[0]
        p_hv = ls_grid.get_trafohv_res()[0]
        p_lv = ls_grid.get_trafolv_res()[0]
        Sbus = ls_grid.get_Sbus()
        Ybus = ls_grid.get_Ybus()
        mis = V * (Ybus * V).conjugate() - Sbus
        mis_pp = V_pp * (Ybus * V_pp).conjugate() - Sbus

        # bus id 2
        # gen 1
        # load 1
        # l2 ex
        # l5 or
        assert abs(-gen_p[1] + load_p[1] + p_ex[2] + p_or[5]) <= 1e-6

        # bus id 0
        # l0 or
        # l1 or
        # gen 4
        assert abs(-gen_p[4] + p_or[0] + p_or[1]) <= 1e-6

        # bus id 3
        # l3 => to
        # l5 => to
        # l6 => from
        # load 2
        # t0 hv
        # t1 hv

        assert abs(load_p[2] + p_ex[3] + p_ex[5] + p_or[6] + p_hv[0] + p_hv[1]) <= 1e-6

        # bus 5
        # l7 or
        # l8 or
        # l9 or
        # g2
        # l4
        # t2 lv
        assert abs(p_or[7] + p_or[8] + p_or[9] - gen_p[2] + load_p[4] + p_lv[2]) <= 1e-6

        mis = ls_grid.check_solution(V, False)
        assert np.abs(mis).max() <= 1e-6, "error for lighsim2grid"
        mis_pp = ls_grid.check_solution(V_pp, False)
        assert np.abs(mis_pp).max() <= 1e-6, "error for pandapower"


        ### re run the powerflow but with different weights
        net.gen["slack"][[1]] = True
        if "slack_weight" in net.gen:
            net.gen["slack_weight"][[1, id_ref_slack]] = 0.7, 0.3
        pp.rundcpp(net)

        pp.runpp(net, init="flat", distributed_slack=True)
        ls_grid2 = init(net)
        V = np.ones(nb_bus_total, dtype=np.complex_) 
        V = ls_grid2.ac_pf(V, max_it, tol)
        my_ref = np.where(np.angle(V) == 0.)[0][0]
        V_pp = net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        assert np.abs(V_pp - V).max() <= 1e-6
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid2.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)

        assert np.all(np.abs([el.res_q_mvar for el in ls_grid2.get_generators()] - net.res_gen["q_mvar"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid2.get_generators()] - net.res_gen["p_mw"].values) <= 1e-6)

        # TODO SLACK: hardest test 1: multiple slack gen at the same bus
        # TODO SLACK: hardest test 2: at a bus, have both slack gen and non slack gen
        # TODO SLACK: measure the performance of all this