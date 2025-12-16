# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import pandapower as pp
import pandapower.networks as pn
import numpy as np
import scipy
import warnings

from pandapower.pypower.makeLODF import update_LODF_diag
from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.solver import SolverType

import pdb

class TestLODFCase14SLU(unittest.TestCase):
    def make_grid(self):
        case14 = pn.case14()
        return case14

    def get_solver_type(self):
        return SolverType.DC
    
    def setUp(self) -> None:
        self.case = self.make_grid()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.gridmodel = init_from_pandapower(self.case)
        self.V_init = 1. * self.gridmodel.get_bus_vn_kv()
        solver_type = self.get_solver_type()
        if solver_type not in self.gridmodel.available_solvers():
            self.skipTest("Solver type not supported on this platform")
        self.gridmodel.change_solver(solver_type)
        V = self.gridmodel.dc_pf(1. * self.V_init, 1, 1e-8)
        assert len(V), f"dc pf has diverged with error {self.gridmodel.get_dc_solver().get_error()}"
        self.dcYbus = 1.0 * self.gridmodel.get_dcYbus_solver()
        self.dcSbus = 1.0 * self.gridmodel.get_dcSbus_solver().real
        self.Bbus = 1.0 * self.dcYbus.real
        self.res_powerflow = 1.0 * np.concatenate((self.gridmodel.get_lineor_res()[0], self.gridmodel.get_trafohv_res()[0]))
        self.nb = self.case.bus.shape[0]
        self.nbr = self.case.line.shape[0] + self.case.trafo.shape[0]
        self.slack_bus = self.case.ext_grid.iloc[0]["bus"]
        self.noref = np.arange(1, self.nb)      ## use bus 1 for voltage angle reference
        self.noslack = np.flatnonzero(np.arange(self.nb) != self.slack_bus)
        self.tol = 1e-6
        return super().setUp()

    def test_from_PTDF(self):
        """compare the ls implementation pypower (in pandapowerer)
        implementation and see if they match"""
        PTDF = 1.0 * self.gridmodel.get_ptdf_solver()
        
        # pypower implementation
        f_ = np.concatenate((1 * self.gridmodel.get_lines().get_bus_id_side_1(), 1 * self.gridmodel.get_trafos().get_bus_id_side_1()))
        t_ = np.concatenate((1 * self.gridmodel.get_lines().get_bus_id_side_2(), 1 * self.gridmodel.get_trafos().get_bus_id_side_2()))
        Cft = scipy.sparse.csr_matrix((np.r_[np.ones(self.nbr), -np.ones(self.nbr)],
                                       (np.r_[f_, t_], np.r_[np.arange(self.nbr), np.arange(self.nbr)])),
                                      (self.nb, self.nbr))

        H = PTDF * Cft
        h = np.diag(H, 0)
        den = (np.ones((self.nbr, 1)) * h.T * -1 + 1.)
        with np.errstate(divide='ignore', invalid='ignore'):
            LODF_pypower = (H / den)
        update_LODF_diag(LODF_pypower)
        # end pypower implementation
        
        LODF_ls = 1.0 * self.gridmodel.get_lodf()
        
        for l_id in range(LODF_ls.shape[0]):
            pypow = 1. * LODF_pypower[:,l_id]
            ls = 1. * LODF_ls[:,l_id]
            if np.all(np.isfinite(pypow)):
                # if everything works on pyower, it should work on lightsim2grid too
                assert np.all(np.isfinite(ls)), f"error for line id {l_id}"
            if np.all(np.isfinite(ls)):
                # and the opposite too
                assert np.all(np.isfinite(pypow)), f"error for line id {l_id}"
            if np.any(~np.isfinite(pypow)):
                # if it does not work for a given line, then 
                # every other flow are nan
                assert np.all(~np.isfinite(pypow)), f"error for line id {l_id}"
                # every flow for ls are nans too
                assert np.all(~np.isfinite(ls)), f"error for line id {l_id}"
            if np.any(~np.isfinite(ls)):
                # if it does not work for a given line, then 
                # every other flow are nan
                assert np.all(~np.isfinite(ls)), f"error for line id {l_id}"  
                # every flow for pypower are nans too
                assert np.all(~np.isfinite(pypow)), f"error for line id {l_id}" 
            if np.all(np.isfinite(pypow)):
                assert np.abs(pypow - ls).max() <= 1e-6, f"error for line id {l_id}"

    def test_with_powerflow(self):
        """compute powerflow with LODF and with "standard" powerflow
        and see if it matches
        """
        LODF_ls = 1.0 * self.gridmodel.get_lodf()
        nb_powerlines = len(self.gridmodel.get_lines())
        for l_id in range(LODF_ls.shape[0]):
            if l_id < nb_powerlines:
                self.gridmodel.deactivate_powerline(l_id)
            else:
                self.gridmodel.deactivate_trafo(l_id - nb_powerlines)
                
            V = self.gridmodel.dc_pf(1. * self.V_init, 1, 1e-8)
            if V.shape[0] > 0:
                # it has converged
                por_pow = np.concatenate((self.gridmodel.get_lineor_res()[0],
                                          self.gridmodel.get_trafohv_res()[0]))
                por_lodf = self.res_powerflow + LODF_ls[:, l_id] * self.res_powerflow[l_id]
                assert np.abs(por_lodf - por_pow).max() <= 1e-6, f"error for line id {l_id}"
            else:
                # it has diverged
                assert np.all(~np.isfinite(LODF_ls[:,l_id])), f"error for line id {l_id}"
                
            if l_id < nb_powerlines:
                self.gridmodel.reactivate_powerline(l_id)
            else:
                self.gridmodel.reactivate_trafo(l_id - nb_powerlines)

class TestLODFCase30SLU(TestLODFCase14SLU):
    def make_grid(self):
        res = pn.case30()
        return res


class TestLODFCase118SLU(TestLODFCase14SLU):
    def make_grid(self):
        res = pn.case118()
        return res
    
    
class TestLODFCase14KLU(TestLODFCase14SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
class TestLODFCase30KLU(TestLODFCase30SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
class TestLODFCase118KLU(TestLODFCase118SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
if __name__ == "__main__":
    unittest.main()
