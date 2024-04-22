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
import warnings
from scipy.sparse.linalg import spsolve

from lightsim2grid.gridmodel import init
from lightsim2grid.solver import SolverType

import pdb

class TestCase14SLU(unittest.TestCase):
    def make_grid(self):
        case14 = pn.case14()
        return case14

    def get_solver_type(self):
        return SolverType.DC
    
    def setUp(self) -> None:
        self.case = self.make_grid()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.gridmodel = init(self.case)
        self.V_init = 1. * self.gridmodel.get_bus_vn_kv()
        solver_type = self.get_solver_type()
        if solver_type not in self.gridmodel.available_solvers():
            self.skipTest("Solver type not supported on this platform")
        self.gridmodel.change_solver(solver_type)
        V = self.gridmodel.dc_pf(self.V_init, 1, 1e-8)
        assert len(V), f"dc pf has diverged with error {self.gridmodel.get_dc_solver().get_error()}"
        self.dcYbus = 1.0 * self.gridmodel.get_dcYbus()
        self.dcSbus = 1.0 * self.gridmodel.get_dcSbus().real
        self.Bbus = 1.0 * self.dcYbus.real
        self.res_powerflow = 1.0 * np.concatenate((self.gridmodel.get_lineor_res()[0], self.gridmodel.get_trafohv_res()[0]))
        self.nb = self.case.bus.shape[0]
        self.nbr = self.case.line.shape[0] + self.case.trafo.shape[0]
        self.slack_bus = self.case.ext_grid.iloc[0]["bus"]
        self.noref = np.arange(1, self.nb)      ## use bus 1 for voltage angle reference
        self.noslack = np.flatnonzero(np.arange(self.nb) != self.slack_bus)
        self.tol = 1e-6
        return super().setUp()

    def test_from_pp(self):
        # test the right computation of Bf matrix (PTDF derived from python, might be slower)
        Bf = 1.0 * self.gridmodel.get_Bf()
        PTDF = np.zeros((self.nbr, self.nb))
        PTDF[:, self.noslack] = spsolve(self.Bbus[np.ix_(self.noslack, self.noref)].T, Bf[:, self.noref].toarray().T).T
        # test the solver works correctly
        tmp_mat = self.Bbus[np.ix_(self.noslack, self.noref)].T.todense()
        for line_id in range(self.nbr):
            solve_error = np.abs(np.dot(tmp_mat, PTDF[line_id, self.noslack]).T - Bf[line_id, self.noref].toarray().T).max()
            assert solve_error <= self.tol, f"error for line {line_id}: {solve_error:.2e}MW"
        # test powerflow are correct
        res_ptdf = np.dot(PTDF, self.dcSbus * self.gridmodel.get_sn_mva())
        assert np.abs(res_ptdf - self.res_powerflow).max() <= self.tol, f"max error for powerflow: {np.abs(res_ptdf - self.res_powerflow).max():.2e}MW"
    
    def test_ptdf_from_ls(self):
        # now test the right computation of the PTDF
        Bf = 1.0 * self.gridmodel.get_Bf()
        PTDF2 = 1.0 * self.gridmodel.get_ptdf()
        # test the solver works correctly
        tmp_mat = self.Bbus[np.ix_(self.noslack, self.noref)].T.todense()
        for line_id in range(self.nbr):
            solve_error = np.abs(np.dot(tmp_mat, PTDF2[line_id, self.noslack]).T - Bf[line_id, self.noref].toarray().T).max()
            assert solve_error <= self.tol, f"error for line {line_id}: {solve_error:.2e}MW"
        # test powerflow are correct
        res_ptdf2 = np.dot(PTDF2, self.dcSbus * self.gridmodel.get_sn_mva())
        np.where(np.abs(res_ptdf2 - self.res_powerflow) >= 1e3)
        assert np.abs(res_ptdf2 - self.res_powerflow).max() <= self.tol, f"max error for powerflow: {np.abs(res_ptdf2 - self.res_powerflow).max():.2e}MW"


class TestCase30SLU(TestCase14SLU):
    def make_grid(self):
        res = pn.case30()
        return res


class TestCase118SLU(TestCase14SLU):
    def make_grid(self):
        res = pn.case118()
        return res
    
    
class TestCase14KLU(TestCase14SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
class TestCase30KLU(TestCase30SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
class TestCase118KLU(TestCase118SLU):
    def get_solver_type(self):
        return SolverType.KLUDC
    
    
if __name__ == "__main__":
    unittest.main()
