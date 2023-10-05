# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import warnings

import pandapower.networks
import pandapower as pp

from pandapower.pypower.makeBdc import makeBdc
from pandapower.pypower.makeSbus import makeSbus
from pandapower.pypower.idx_bus import GS
from pandapower.pf.ppci_variables import _get_pf_variables_from_ppci
from pandapower.pd2ppc import _pd2ppc
from pandapower.auxiliary import _init_runpp_options

import unittest

from lightsim2grid.gridmodel import init
from lightsim2grid.solver import SolverType


class BaseMVOberrheinTester(unittest.TestCase):
    def get_network(self):
        res = pandapower.networks.mv_oberrhein()
        res.switch["closed"] = True
        res.trafo["shift_degree"] = 0.  # ignored in pandapower since `calculate_voltage_angles=False` 
        # (otherwise powerflow diverge)
        return res

    def _aux_get_init_pp_data(self):
        # init pp options
        _init_runpp_options(self.net, **self.get_pp_options())
        
        # clear lookups
        self.net._pd2ppc_lookups = {"bus": np.array([], dtype=np.int64), "ext_grid": np.array([], dtype=np.int64),
                                    "gen": np.array([], dtype=np.int64), "branch": np.array([], dtype=np.int64)}

        # convert pandapower net to ppc
        ppc, ppci = _pd2ppc(self.net)
        self.net["_ppc"] = ppc
        baseMVA, bus, gen, branch, svc, tcsc, ref, pv, pq, *_, gbus, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
        return pv, pq, V0
        
    def get_pp_options(self):
        return dict(algorithm="nr",
                    calculate_voltage_angles=False,
                    init="flat",
                    max_iteration="auto",
                    # max_iteration="auto",
                    tolerance_mva=self.tol_solver,
                    trafo_model="t",
                    trafo_loading="current",
                    enforce_q_lims=False,
                    check_connectivity=True,
                    voltage_depend_loads=True, 
                    consider_line_temperature=False,
                    run_control=False,
                    distributed_slack=False,
                    tdpf=False,
                    tdpf_delay_s=None)
        
    def setUp(self) -> None:
        self.net = self.get_network()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.gridmodel = init(self.net)
        self.tol = 1e-7
        self.tol_solver = 1e-8
        
    def test_Ybus(self):
        # Ybus seen by pandapower
        pp.runpp(self.net, **self.get_pp_options())
        Ybus_pp_tmp = self.net["_ppc"]["internal"]["Ybus"]
        id_bus_pp = self.net._pd2ppc_lookups["bus"][self.gridmodel._ls_to_pp]
        assert (id_bus_pp == np.arange(id_bus_pp.size)).all()  # assert bus are sorted correctly and in the same order
        # Ybus_pp = Ybus_pp_tmp[id_bus_pp, id_bus_pp]  # does not work in don't know why... anyway
        Ybus_pp = Ybus_pp_tmp
        *_, V0 = self._aux_get_init_pp_data()
        # Ybus from lightsim2grid
        V_ls = self.gridmodel.ac_pf(1. * V0, 30, self.tol_solver)
        Ybus_ls = self.gridmodel.get_Ybus()
        # np.where(np.abs(Ybus_pp - Ybus_ls).todense() >=100.)
        assert np.abs((Ybus_pp - Ybus_ls)).max() <= self.tol, f"error in Ybus: max {np.abs((Ybus_pp - Ybus_ls).todense()).max():.2e}"
        
    def test_Sbus(self):
        # Ybus seen by pandapower
        pp.runpp(self.net, **self.get_pp_options())
        bus_lookup = self.net._pd2ppc_lookups["bus"]  #[self.gridmodel._ls_to_pp]
        Sbus_pp = self.net["_ppc"]["internal"]["Sbus"]
        *_, V0 = self._aux_get_init_pp_data()
        # Ybus from lightsim2grid
        V_ls = self.gridmodel.ac_pf(1. * V0, 30, self.tol_solver)
        Sbus_ls = self.gridmodel.get_Sbus()
        slack_id = bus_lookup[self.net.ext_grid["bus"].values]
        is_not_slack = np.ones(self.net.bus.shape[0], dtype=bool)
        is_not_slack[slack_id] = False
        # np.where((np.abs(Sbus_pp - Sbus_ls) >= 1e-3) & is_not_slack)
        assert np.abs(Sbus_pp[is_not_slack] - Sbus_ls[is_not_slack]).max() <= self.tol, f"error in Sbus: max {np.abs((Sbus_pp - Sbus_ls)).max():.2e}"

    def test_ac_solver(self):
        *_, V0 = self._aux_get_init_pp_data()
        self.gridmodel.change_solver(SolverType.KLUSingleSlack)
        V_ls = self.gridmodel.ac_pf(1. * V0, 30, self.tol_solver)
        pp.runpp(self.net, **self.get_pp_options()) 
        assert V_ls.size > 0, f"lightsim2grid powerflow has diverged {self.gridmodel.get_solver().get_error()}"
        nb_iter_ls = self.gridmodel.get_solver().get_nb_iter()
        nb_iter_pp = self.net._ppc['iterations']
        assert nb_iter_pp == nb_iter_ls, f"mismatch nb_iter {nb_iter_pp} vs {nb_iter_ls}"
        V_pp = self.net.res_bus["vm_pu"].values * np.exp(1j * np.deg2rad(self.net.res_bus["va_degree"].values))
        assert np.abs(V_ls - V_pp).max() <= self.tol, f"mismatch in V: {np.abs(V_ls - V_pp).max():.2e}"

    def test_dc_Ybus(self):
        # Ybus seen by pandapower
        pp.rundcpp(self.net)
        Ybus_pp = self.net["_ppc"]["internal"]["Bbus"]
        *_, V0 = self._aux_get_init_pp_data()
        # Ybus from lightsim2grid
        V_ls = self.gridmodel.dc_pf(1. * V0, 30, self.tol_solver)
        Ybus_ls = self.gridmodel.get_dcYbus()
        # np.where(np.abs( (Ybus_pp - Ybus_ls).todense()) >= 1e-4)
        assert np.abs( (Ybus_pp - Ybus_ls).todense()).max() <= self.tol, f"error in Ybus (dc): max {np.abs((Ybus_pp - Ybus_ls).todense()).max():.2e}"

    def test_dc_Sbus(self):
        self.skipTest("I don't really know what pandapower is doing")
        # Ybus seen by pandapower
        pp.rundcpp(self.net)
        ppci = self.net._ppc
        baseMVA = ppci['baseMVA']
        bus = ppci['bus']
        gen = ppci['gen']
        branch = ppci['branch']
        B, Bf, Pbusinj, Pfinj, Cft = makeBdc(bus, branch)
        Sbus_pp = makeSbus(baseMVA, bus, gen) - Pbusinj - bus[:, GS] / baseMVA  # named Pbus in pandapower
        self.net["_ppc"]["internal"]["Bbus"] * np.deg2rad(self.net.res_bus["va_degree"].values) - Sbus_pp
        
        *_, V0 = self._aux_get_init_pp_data()
        # Ybus from lightsim2grid
        V_ls = self.gridmodel.dc_pf(1. * V0, 30, self.tol_solver)
        Sbus_ls = self.gridmodel.get_dcSbus()
        slack_id = self.net.ext_grid["bus"].values
        is_not_slack = np.ones(self.net.bus.shape[0], dtype=bool)
        is_not_slack[slack_id] = False
        # np.where((np.abs(Sbus_pp - Sbus_ls) >= 1e-3) & is_not_slack)
        assert np.abs(Sbus_pp[is_not_slack] - Sbus_ls[is_not_slack]).max() <= self.tol, f"error in Sbus: max {np.abs((Sbus_pp - Sbus_ls)).max():.2e}"

    def test_dc_solver(self):
        self.skipTest("I don't really know what pandapower is doing")
        *_, V0 = self._aux_get_init_pp_data()
        V_ls = self.gridmodel.dc_pf(1. * V0, 30, self.tol_solver)
        pp.rundcpp(self.net, tol=self.tol_solver, init="flat")
        assert V_ls.size > 0, f"lightsim2grid powerflow has diverged {self.gridmodel.get_dc_solver().get_error()}"
        V_pp = self.net.res_bus["vm_pu"].values * np.exp(1j * np.deg2rad(self.net.res_bus["va_degree"].values))
        import pdb
        pdb.set_trace()
        assert np.abs(V_ls - V_pp).max() <= self.tol, f"mismatch in V (for DC): {np.abs(V_ls - V_pp).max():.2e}"
        
        
if __name__ == "__main__":
    unittest.main()
