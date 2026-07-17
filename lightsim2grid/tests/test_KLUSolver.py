# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import unittest
import numpy as np
import pdb
import zipfile
from scipy import sparse
KLU_AVAILBLE = False
try:
    from lightsim2grid.lightsim2grid_cpp import NR_KLU
    KLU_AVAILBLE = True
except ImportError:
    # KLU solver is not available, these tests cannot be carried out
    pass


class MakeTests(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        unittest.TestCase.__init__(self, methodName=methodName)
        self.methodName = methodName

        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 1e-4  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)
        if not KLU_AVAILBLE:
            return
        self.solver = NR_KLU()

        self.path = None
        self.V_init = None
        self.pq = None
        self.pv = None
        self.Sbus = None
        self.Ybus = None

    def load_array(self, myzip, nm):
        arr = myzip.extract(nm)
        res = np.load(arr)
        os.remove(arr)
        return res

    def load_path(self, path):
        try:
            self.path = path
            with zipfile.ZipFile(self.path) as myzip:
                self.V_init = self.load_array(myzip, "V0.npy")
                self.pq = self.load_array(myzip, "pq.npy")
                self.pv = self.load_array(myzip, "pv.npy")
                self.Sbus = self.load_array(myzip, "Sbus.npy")
                self.Ybus = self.load_array(myzip, "Ybus.npy")
                self.Ybus = sparse.csc_matrix(self.Ybus)
            return True
        except:
            return False

    def _load_vect(self, myzip, iter_max, nm):
        try:
            res = self.load_array(myzip, "{}_{}.npy".format(nm, iter_max))
        except Exception as exc_:
            raise RuntimeError("{} is not defined for iteration {} for case {}".format(nm, iter_max, self.path))
        return res

    def load_res(self, iter_max):
        with zipfile.ZipFile(self.path) as myzip:
            J = self._load_vect(myzip, iter_max, "J")
            V = self._load_vect(myzip, iter_max, "V")
        return J, V

    def compare_sparse_mat(self, J, J_pp, pv, pq):
        """
        Test that the matrices J and J_pp are almost equal
        :param J:
        :param J_pp:
        :param pv:
        :param pq:
        :return:
        """
        pvpq = np.r_[pv, pq]

        comp_val = np.abs(J - J_pp)  # J already remapped to pandapower [pvpq | pq] order in solver_aux
        comp_val = comp_val
        assert np.sum(np.abs(comp_val[:len(pvpq), :len(pvpq)])) <= self.tol_test, "J11 (dS_dVa_r) are not equal"
        assert np.sum(np.abs(comp_val[len(pvpq):, :len(pvpq)])) <= self.tol_test, "J21 (dS_dVa_i) are not equal"
        assert np.sum(np.abs(comp_val[:len(pvpq), len(pvpq):])) <= self.tol_test, "J12 (dS_dVm_r) are not equal"
        assert np.sum(np.abs(comp_val[len(pvpq):, len(pvpq):])) <= self.tol_test, "J22 (dS_dVm_i) are not equal"

    def solver_aux(self):
        self.solver.reset()
        # start = time.perf_counter()
        ref = set(np.arange(self.Sbus.shape[0])) - set(self.pv) - set(self.pq)
        ref = np.array(list(ref))
        # build the slack weights
        slack_weights = np.zeros(self.Sbus.shape[0])
        slack_weights[ref] = 1.0 / ref.shape[0]

        has_conv = self.solver.compute_pf(self.Ybus, self.V_init, self.Sbus,  ref, slack_weights, 
                                          self.pv, self.pq, self.max_it, self.tol)
        assert has_conv, "the load flow has diverged for {}".format(self.path)
        J = self.solver.get_J()
        # new NR layout puts the slack at the end and orders the bus unknowns via the
        # bus_id -> J-column converters; remap to pandapower's [pvpq | pq] order.
        theta = np.asarray(self.solver.get_theta_to_J_col())
        vm = np.asarray(self.solver.get_vm_to_J_col())
        pvpq = np.r_[self.pv, self.pq]
        col_perm = np.concatenate([theta[pvpq], vm[self.pq]])
        J = J[col_perm.reshape(-1, 1), col_perm.reshape(1, -1)]
        Va = self.solver.get_Va()
        Vm = self.solver.get_Vm()
        J_pp, V_pp = self.load_res(iter_max=self.solver.get_nb_iter())
        self.compare_sparse_mat(J, J_pp, self.pv, self.pq)
        Va_pp = np.angle(V_pp)
        Vm_pp = np.abs(V_pp)
        assert np.sum(np.abs(Va - Va_pp)) <= self.tol_test, "voltages angles are not the same"
        assert np.sum(np.abs(Vm - Vm_pp)) <= self.tol_test, "voltages magnitude are not the same"

    def test_dir(self):
        if not KLU_AVAILBLE:
            self.skipTest("KLU is not installed")
        nb_tested = 0
        for path in os.listdir("."):
            _, ext = os.path.splitext(path)
            if ext == ".zip":
                path_ok = self.load_path(path)
                if path_ok:
                    self.solver_aux()
                    nb_tested += 1
        assert nb_tested == 5, "incorrect number of test cases found, found {} while there should be 5".format(nb_tested)

    def _find_and_load_one_case(self):
        for path in os.listdir("."):
            _, ext = os.path.splitext(path)
            if ext == ".zip" and self.load_path(path):
                return
        raise RuntimeError("no test case zip found in the current directory")

    def test_linear_solver_stats(self):
        """get_linear_solver_stats() (LinearSolverPolicy wrapper) reports non-trivial,
        consistent counters after a successful powerflow."""
        if not KLU_AVAILBLE:
            self.skipTest("KLU is not installed")
        self._find_and_load_one_case()
        self.solver_aux()
        stats = self.solver.get_linear_solver_stats()
        assert stats.nb_analyze >= 1, "analyze() should have been called at least once"
        assert stats.nb_factorize >= 1, "factorize() should have been called at least once"
        assert stats.nb_solve >= 1, "solve() should have been called at least once"
        assert stats.nb_refactorize_failed == 0, "no refactor should fail on a well-behaved test case"
        assert stats.nb_fallback_factorize == 0, "no fallback factorize expected on plain NR_KLU"
        assert stats.timer_factor >= 0.0
        assert stats.timer_solve >= 0.0

    def test_refactor_retry_solver(self):
        """NRRefactorRetry_KLU (RefactorRetryLinearSolver) converges the same as NR_KLU
        on a well-behaved case, with zero fallback fires (sanity, not a regression)."""
        if not KLU_AVAILBLE:
            self.skipTest("KLU is not installed")
        try:
            from lightsim2grid.lightsim2grid_cpp import NRRefactorRetry_KLU
        except ImportError:
            self.skipTest("NRRefactorRetry_KLU is not available")
        self._find_and_load_one_case()
        solver = NRRefactorRetry_KLU()
        ref = set(np.arange(self.Sbus.shape[0])) - set(self.pv) - set(self.pq)
        ref = np.array(list(ref))
        slack_weights = np.zeros(self.Sbus.shape[0])
        slack_weights[ref] = 1.0 / ref.shape[0]
        has_conv = solver.compute_pf(self.Ybus, self.V_init, self.Sbus, ref, slack_weights,
                                      self.pv, self.pq, self.max_it, self.tol)
        assert has_conv, "NRRefactorRetry_KLU failed to converge on a case where NR_KLU succeeds"
        stats = solver.get_linear_solver_stats()
        assert stats.nb_factorize >= 1
        assert stats.nb_refactorize_failed == 0, "unexpected refactor failure on a well-behaved test case"
        assert stats.nb_fallback_factorize == 0, "unexpected fallback factorize on a well-behaved test case"


if __name__ == "__main__":
    unittest.main()
