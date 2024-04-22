# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
from scipy.sparse import csc_matrix

from lightsim2grid_cpp import FDPFMethod
from lightsim2grid.gridmodel import init
from lightsim2grid.solver import SolverType 

import pandapower.networks as pn
import pandapower as pp
from pandapower.pypower.makeB import makeB
from pandapower.pf.ppci_variables import _get_pf_variables_from_ppci
from pandapower.pd2ppc import _pd2ppc
from pandapower.auxiliary import _init_runpp_options

import unittest



class BaseFDPFTester:
    def get_solving_method(self):
        return FDPFMethod.XB
    
    def get_algo_pp(self):
        return "fdxb" if self.fdpf_meth == FDPFMethod.XB else "fdbx"  # in theory
        # return "fdbx" if self.fdpf_meth == FDPFMethod.XB else "fdxb"  # but.... https://github.com/e2nIEE/pandapower/issues/2142
    
    def get_network(self):
        return pn.case14()
        # self.net = pn.case118()
        # self.net = pn.case300()
    
    def get_pp_options(self):
        return dict(algorithm=self.get_algo_pp(),
                    calculate_voltage_angles="auto",
                    init="flat",
                    max_iteration="auto",
                    tolerance_mva=1e-8,
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
        
        # XB
        self.fdpf_meth = self.get_solving_method()
        self.alg = 2 if self.fdpf_meth == FDPFMethod.XB else 3
        
        # alg = 3  # BX  # TODO later
        # self.fdpf_meth = FDPFMethod.BX  # TODO later
        return super().setUp()
    
    def _aux_get_Bp_Bpp(self):
        # init pp options
        _init_runpp_options(self.net, **self.get_pp_options())
        
        # clear lookups
        self.net._pd2ppc_lookups = {"bus": np.array([], dtype=np.int64), "ext_grid": np.array([], dtype=np.int64),
                                    "gen": np.array([], dtype=np.int64), "branch": np.array([], dtype=np.int64)}

        # convert pandapower net to ppc
        ppc, ppci = _pd2ppc(self.net)
        self.net["_ppc"] = ppc
        baseMVA, bus, gen, branch, svc, tcsc, ssc, ref, pv, pq, on, gbus, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
        pp_Bp, pp_Bpp = makeB(baseMVA, bus, np.real(branch), self.alg)
        return pp_Bp, pp_Bpp, pv, pq, V0
        
    def test_bp_bpp(self):
        """test that Bp and Bpp are correct (from the grid)"""
        pp_Bp, pp_Bpp, *_= self._aux_get_Bp_Bpp()
        
        # now get the matrices 
        self.gridmodel.ac_pf(1.04 * np.ones(self.net.bus.shape[0], dtype=complex), 10, 1e-7)  # need to init the "model bus id" id_grid_to_solver
        ls_Bp = self.gridmodel.debug_get_Bp_python(self.fdpf_meth)
        ls_Bpp = self.gridmodel.debug_get_Bpp_python(self.fdpf_meth)
        
        assert np.abs(pp_Bp - ls_Bp).max() <= self.tol, f"error in Bp: max {np.abs(pp_Bp - ls_Bp).max():.2e}"
        assert np.abs(pp_Bpp - ls_Bpp).max() <= self.tol, f"error in Bpp: max {np.abs(pp_Bpp - ls_Bpp).max():.2e}"

    def test_Bp_Bpp_solver(self):
        """test that Bp and Bpp are correct (from the solver: where only some indices are used)"""
        # retrieve the matrix in the solver
        self.gridmodel.change_solver(SolverType.FDPF_XB_SparseLU if self.fdpf_meth == FDPFMethod.XB else SolverType.FDPF_BX_SparseLU)
        V_ls = self.gridmodel.ac_pf(1.04 * np.ones(self.net.bus.shape[0], dtype=complex), 30, 1e0)  # to ensure it can "converge" to be able to retrieve Bp and Bpp
        used_solver = self.gridmodel.get_solver().get_fdpf_xb_lu() if self.fdpf_meth == FDPFMethod.XB else self.gridmodel.get_solver().get_fdpf_bx_lu() 
        ls_Bp = used_solver.debug_get_Bp_python()
        ls_Bpp = used_solver.debug_get_Bpp_python()
        # recreate the matrices from pandapower
        grid_Bp, grid_Bpp, pv, pq, *_ = self._aux_get_Bp_Bpp()
        pvpq = np.r_[pv, pq]
        pp_Bp = grid_Bp[np.array([pvpq]).T, pvpq].tocsc()
        pp_Bpp = grid_Bpp[np.array([pq]).T, pq].tocsc()
        
        # check they match
        # NB: if pandapower change signature of `_get_pf_variables_from_ppci`
        # you might get weird error here
        assert np.abs(pp_Bp - ls_Bp).max() <= self.tol, f"error in Bp: max {np.abs(pp_Bp - ls_Bp).max():.2f}"
        assert np.abs(pp_Bpp - ls_Bpp).max() <= self.tol, f"error in Bpp: max {np.abs(pp_Bpp - ls_Bpp).max():.2f}"
        
    def test_solver(self):
        self.gridmodel.change_solver(SolverType.FDPF_XB_SparseLU if self.fdpf_meth == FDPFMethod.XB else SolverType.FDPF_BX_SparseLU)
        *_, V0 = self._aux_get_Bp_Bpp()
        V_ls = self.gridmodel.ac_pf(1. * V0, 30, self.tol_solver * self.gridmodel.get_sn_mva())  # no division by sn_mva in pypower implementation TODO issue !
        pp.runpp(self.net, algorithm=self.get_algo_pp(), tol=self.tol_solver, init="flat") 
        assert V_ls.size > 0, f"lightsim2grid powerflow has diverged {self.gridmodel.get_solver().get_error()}"
        nb_iter_ls = self.gridmodel.get_solver().get_nb_iter()
        nb_iter_pp = self.net._ppc['iterations']
        assert nb_iter_pp == nb_iter_ls, f"mismatch nb_iter {nb_iter_pp} vs {nb_iter_ls} maybe https://github.com/e2nIEE/pandapower/issues/2142 has been fixed ?"
        V_pp = self.net.res_bus["vm_pu"].values * np.exp(1j * np.deg2rad(self.net.res_bus["va_degree"].values))
        assert np.abs(V_ls - V_pp).max() <= self.tol, f"mismatch in V: {np.abs(V_ls - V_pp).max():.2e}"


class FDPFTester_Case14_XB(BaseFDPFTester, unittest.TestCase):
    pass


class FDPFTester_Case14_BX(BaseFDPFTester, unittest.TestCase):
    def get_solving_method(self):
        return FDPFMethod.BX


class FDPFTester_Case118_XB(BaseFDPFTester, unittest.TestCase):
    def get_network(self):
        return pn.case118()


class FDPFTester_Case118_BX(BaseFDPFTester, unittest.TestCase):
    def get_network(self):
        return pn.case118()
    def get_solving_method(self):
        return FDPFMethod.BX


class FDPFTester_Case300_XB(BaseFDPFTester, unittest.TestCase):
    def get_network(self):
        return pn.case300()


class FDPFTester_Case300_BX(BaseFDPFTester, unittest.TestCase):
    def get_network(self):
        return pn.case300()
    def get_solving_method(self):
        return FDPFMethod.BX


if __name__ == "__main__":
    unittest.main()
