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
import zipfile
from scipy import sparse

SparseLUSolver_AVAILBLE = False
try:
    from lightsim2grid.lightsim2grid_cpp import NR_SparseLU
    SparseLUSolver_AVAILBLE = True
except ImportError:
    # SparseLU solver is not available, these tests cannot be carried out
    pass


class TestJColConverters(unittest.TestCase):
    """Check the bus_id -> Jacobian-column converters exposed by the NR solver.

    For each bus the solver returns the Jacobian column of that bus' theta / vm /
    q unknown, or -1 if the bus owns no unknown of that type. The Jacobian columns
    are laid out as: theta of pvpq buses in [0, n_pvpq), then vm of pq buses in
    [n_pvpq, n_pvpq + n_pq). q is a placeholder (always -1 for now).

    These invariants are structural (derived from pv/pq) and therefore independent
    of the numerical Jacobian values.
    """

    def __init__(self, methodName='runTest'):
        unittest.TestCase.__init__(self, methodName=methodName)
        self.max_it = 10
        self.tol = 1e-8
        if not SparseLUSolver_AVAILBLE:
            return
        self.solver = NR_SparseLU()
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
        except Exception:
            return False

    def check_converters(self):
        self.solver.reset()
        n_bus = self.Sbus.shape[0]
        ref = np.array(list(set(np.arange(n_bus)) - set(self.pv) - set(self.pq)))
        slack_weights = np.zeros(n_bus)
        slack_weights[ref] = 1.0 / ref.shape[0]

        has_conv = self.solver.compute_pf(self.Ybus, self.V_init, self.Sbus, ref,
                                          slack_weights, self.pv, self.pq,
                                          self.max_it, self.tol)
        assert has_conv, "the load flow has diverged for {}".format(self.path)

        J = self.solver.get_J()
        theta = np.asarray(self.solver.get_theta_to_J_col())
        vm = np.asarray(self.solver.get_vm_to_J_col())
        q = np.asarray(self.solver.get_q_to_J_col())

        # all three maps have one entry per bus
        assert theta.shape[0] == n_bus, "theta_to_J_col has wrong size"
        assert vm.shape[0] == n_bus, "vm_to_J_col has wrong size"
        assert q.shape[0] == n_bus, "q_to_J_col has wrong size"

        # q is a placeholder for now: no bus owns a reactive unknown
        assert np.all(q == -1), "q_to_J_col must be -1 everywhere for now"

        # expected base-block layout: theta of pvpq in [0, n_pvpq), vm of pq next
        pvpq = np.r_[self.pv, self.pq]
        n_pvpq = pvpq.shape[0]
        expected_theta = -np.ones(n_bus, dtype=int)
        expected_theta[pvpq] = np.arange(n_pvpq)
        expected_vm = -np.ones(n_bus, dtype=int)
        expected_vm[self.pq] = n_pvpq + np.arange(self.pq.shape[0])

        assert np.array_equal(theta, expected_theta), \
            "theta_to_J_col does not match the expected pvpq column layout"
        assert np.array_equal(vm, expected_vm), \
            "vm_to_J_col does not match the expected pq column layout"

        # the reference (slack) buses own no theta and no vm unknown
        assert np.all(theta[ref] == -1), "slack bus must not own a theta column"
        assert np.all(vm[ref] == -1), "slack bus must not own a vm column"

        # every pq bus owns both a theta and a vm unknown
        assert np.all(theta[self.pq] >= 0), "every pq bus must own a theta column"
        assert np.all(vm[self.pq] >= 0), "every pq bus must own a vm column"

        # every reported column is a valid column of J
        ncols = J.shape[1]
        assert np.all(theta[theta >= 0] < ncols), "theta column index out of range"
        assert np.all(vm[vm >= 0] < ncols), "vm column index out of range"

    def test_dir(self):
        if not SparseLUSolver_AVAILBLE:
            self.skipTest("NR_SparseLU is not installed")
        nb_tested = 0
        for path in os.listdir("."):
            _, ext = os.path.splitext(path)
            if ext == ".zip":
                if self.load_path(path):
                    self.check_converters()
                    nb_tested += 1
        assert nb_tested == 5, \
            "incorrect number of test cases found, found {} while there should be 5".format(nb_tested)


if __name__ == "__main__":
    unittest.main()
