# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import sys
import unittest

import numpy as np
from lightsim2grid.algorithm import (
    NR_KLU,
    NRSing_KLU
)

# the fixtures are stored as numpy-version-independent .npz archives (no pickle),
# see jacobian_case14_sandbox/jacobian_io.py
sys.path.append(os.path.join(os.path.dirname(__file__), "jacobian_case14_sandbox"))
from jacobian_io import load_res


class JacobianMultiSlackTester(unittest.TestCase):
    @classmethod
    def get_name(cls):
        return os.path.join("jacobian_case14_sandbox","saved_jacobian_multi.npz")
    
    def get_algo(self):
        return NR_KLU()
    
    @classmethod
    def setUpClass(cls) -> None:
        cls.res = load_res(cls.get_name())

        return super().setUpClass()
    
    def setUp(self) -> None:
        self.nr_algo = self.get_algo()
        return super().setUp()    
    
    def new_to_old_indexes(self):
        """
        in lightsim2grid 1.0.0, the Jacobian rows / columns are laid out by
        increasing bus index within each block: base block first (theta / P of
        sorted pvpq, then vm / Q of sorted pq), then the distributed-slack
        extension (sorted non-ref slacks, slack_absorbed / ref equation last).

        the stored reference uses the historical layout:
        slack eq, P at slack_ids[1:] (input order), P at pv then pq (input order), Q at pq

        this should be used J_new_indx[self.new_to_old_indexes().T, self.new_to_old_indexes()] == J_old_index
        """
        init = type(self).res["init_state"]
        pv = np.asarray(init["pv"]).reshape(-1)
        pq = np.asarray(init["pq"]).reshape(-1)
        slack_ids = np.asarray(init["slack_ids"]).reshape(-1)
        pv_slack = slack_ids[1:]  # non-ref slack buses (input order)

        pvpq_sorted = np.sort(np.r_[pv, pq])
        pq_sorted = np.sort(pq)
        n_pvpq = pvpq_sorted.shape[0]
        base_size = n_pvpq + pq_sorted.shape[0]
        return np.concatenate(
            (
                [base_size + pv_slack.shape[0]],  # slack (ref eq / slack_absorbed col)
                base_size + np.searchsorted(np.sort(pv_slack), pv_slack),  # pv_slack
                np.searchsorted(pvpq_sorted, np.r_[pv, pq]),  # P / theta at pv then pq
                n_pvpq + np.searchsorted(pq_sorted, pq),  # Q / vm at pq
                )
            ).reshape(-1,1)
    
    def _aux_test_iter(self, iter):
        cls = type(self)
        ref_J = cls.res[f"{iter}"]["J"].copy()
        
        self.nr_algo.compute_pf(
            cls.res["init_state"]["Ybus"],
            cls.res["init_state"]["v_init"],
            cls.res["init_state"]["Sbus"],
            cls.res["init_state"]["slack_ids"],
            cls.res["init_state"]["slack_weights"],
            cls.res["init_state"]["pv"],
            cls.res["init_state"]["pq"],
            iter,
            cls.res["init_state"]["tol"])
        J_wrong_order = self.nr_algo.get_J()
        if ref_J.shape[0] > 0:
            assert J_wrong_order.shape == ref_J.shape, f"error for iter {iter}: J.shape = {J_wrong_order.shape} != {ref_J.shape}"
            assert J_wrong_order.nnz == ref_J.nnz, f"error for iter {iter}: J.nnz = {J_wrong_order.nnz} != {ref_J.nnz}"
            
            ref_J_data = ref_J.data.copy()
            J_data = J_wrong_order.data.copy()
            ref_J_data.sort()
            J_data.sort()
            assert np.abs(ref_J_data - J_data).max()<= 1e-6, f"wrong sorted data between J and ref_J, {(np.abs(ref_J_data - J_data) >= 1e-6).nonzero()[0]}"
            
            # J is csc !
            J_for_copy = J_wrong_order.copy()
            J_for_copy.data[np.abs(J_for_copy.data) < 1e-15] = 0.  # remove too small element, might differ in ref_J
            J = J_for_copy[self.new_to_old_indexes().T, self.new_to_old_indexes()]
            # I do the following because apparently this removes the 0. coeff from original matrices
            # and I i remove this from J_wrong_order, I need to remove them from J too.
            ref_J.data[np.abs(ref_J.data) < 1e-15] = 0. # ok if this is not exactly the same at this precision...
            init_index = np.arange(ref_J.shape[0]).reshape(-1,1)
            ref_J = ref_J[init_index.T, init_index].copy()
            assert (J.indptr == ref_J.indptr).all(), f"error for iter {iter}: J.indptr = {J.indptr} != {ref_J.indptr}"
            assert (J.indices == ref_J.indices).all(), f"error for iter {iter}: J.indices = {(J.indices != ref_J.indices).nonzero()[0]}"
            if J.shape[0] > 0:
                assert np.abs(J - ref_J).max() <= 1e-6, f"error for iter {iter}: {np.abs(J - ref_J).max()}"
            
        V = self.nr_algo.get_V()
        assert V.shape == cls.res[f"{iter}"]["V"].shape, f"error for iter {iter}: V.shape = {V.shape} != {cls.res[f'{iter}']['V'].shape}"
        if  V.shape[0] > 0:
            assert np.abs(V - cls.res[f"{iter}"]["V"]).max() <= 1e-6, f"error for iter {iter}: {np.abs(V - cls.res[f'{iter}']['V']).max()}"
        
    def test_O_iter(self):
        self._aux_test_iter(0)
        
    def test_1_iter(self):
        self._aux_test_iter(1)
        
    def test_2_iter(self):
        self._aux_test_iter(2)
        
    def test_3_iter(self):
        self._aux_test_iter(3)
        
    def test_4_iter(self):
        self._aux_test_iter(4)

    def test_col_converters(self):
        """Check the bus_id -> Jacobian-column converters, including the
        distributed-slack (pv_slack) extension columns.

        Column layout (sorted by bus index within each block, slack at the end):
          - theta of sorted base pvpq buses : [0, n_pvpq)
          - vm of sorted pq buses           : [n_pvpq, n_pvpq + n_pq)  (== [n_pvpq, base_size))
          - theta of sorted pv_slack buses  : [base_size, base_size + n_pv_slack)
          - slack_absorbed (no bus)         : base_size + n_pv_slack   (last column)
        """
        cls = type(self)
        init = cls.res["init_state"]
        self.nr_algo.compute_pf(init["Ybus"], init["v_init"], init["Sbus"],
                                init["slack_ids"], init["slack_weights"],
                                init["pv"], init["pq"], 10, init["tol"])

        J = self.nr_algo.get_J()
        theta = np.asarray(self.nr_algo.get_theta_to_J_col())
        vm = np.asarray(self.nr_algo.get_vm_to_J_col())
        q = np.asarray(self.nr_algo.get_q_to_J_col())

        n_bus = init["Ybus"].shape[0]
        pv = np.asarray(init["pv"]).reshape(-1)
        pq = np.asarray(init["pq"]).reshape(-1)
        slack_ids = np.asarray(init["slack_ids"]).reshape(-1)
        ref_slack = slack_ids[0]
        pv_slack = slack_ids[1:]  # distributed-slack pv buses (empty for single slack)

        pvpq = np.r_[pv, pq]
        n_pvpq = pvpq.shape[0]
        n_pq = pq.shape[0]
        base_size = n_pvpq + n_pq

        # one entry per bus
        assert theta.shape[0] == n_bus, "theta_to_J_col has wrong size"
        assert vm.shape[0] == n_bus, "vm_to_J_col has wrong size"
        assert q.shape[0] == n_bus, "q_to_J_col has wrong size"
        # q is a placeholder for now
        assert np.all(q == -1), "q_to_J_col must be -1 everywhere for now"

        # expected full layout (sorted by bus index within each block):
        # base block + pv_slack extension columns
        expected_theta = -np.ones(n_bus, dtype=int)
        expected_theta[np.sort(pvpq)] = np.arange(n_pvpq)
        expected_theta[np.sort(pv_slack)] = base_size + np.arange(pv_slack.shape[0])
        expected_vm = -np.ones(n_bus, dtype=int)
        expected_vm[np.sort(pq)] = n_pvpq + np.arange(n_pq)

        assert np.array_equal(theta, expected_theta), \
            f"theta_to_J_col mismatch: {theta} != {expected_theta}"
        assert np.array_equal(vm, expected_vm), \
            f"vm_to_J_col mismatch: {vm} != {expected_vm}"

        # the reference slack owns neither a theta nor a vm unknown
        assert theta[ref_slack] == -1, "ref slack must not own a theta column"
        assert vm[ref_slack] == -1, "ref slack must not own a vm column"

        # pv_slack buses own a theta column in the extension region, but no vm column
        if pv_slack.shape[0] > 0:
            assert np.all(theta[pv_slack] >= base_size), \
                "pv_slack theta columns must be in the extension region"
            assert np.all(vm[pv_slack] == -1), "pv_slack buses must not own a vm column"

        # the converters must reference exactly the bus-unknown columns
        # [0, base_size + n_pv_slack); any remaining J columns are the
        # slack_absorbed column(s) (present only for the multi-slack extension),
        # which are not per-bus unknowns and so are not referenced by any map.
        used = np.sort(np.concatenate([theta[theta >= 0], vm[vm >= 0]]))
        assert np.array_equal(used, np.arange(base_size + pv_slack.shape[0])), \
            "converter columns must be exactly [0, base_size + n_pv_slack)"
        assert J.shape[1] >= base_size + pv_slack.shape[0], \
            "Jacobian has fewer columns than bus unknowns"
        assert theta.max() < J.shape[1], "theta column index out of range"
        assert vm.max() < J.shape[1], "vm column index out of range"


class JacobianSingleSlackTester(JacobianMultiSlackTester):
    @classmethod
    def get_name(cls):
        return os.path.join("jacobian_case14_sandbox","saved_jacobian_single.npz")
    
    def get_algo(self):
        return NRSing_KLU()
    
    def new_to_old_indexes(self):
        """
        single slack: no extension block; the stored reference uses pv-then-pq
        blocks while the new layout sorts each block by bus index.
        """
        init = type(self).res["init_state"]
        pv = np.asarray(init["pv"]).reshape(-1)
        pq = np.asarray(init["pq"]).reshape(-1)
        pvpq_sorted = np.sort(np.r_[pv, pq])
        pq_sorted = np.sort(pq)
        n_pvpq = pvpq_sorted.shape[0]
        return np.concatenate(
            (
                np.searchsorted(pvpq_sorted, np.r_[pv, pq]),  # P / theta at pv then pq
                n_pvpq + np.searchsorted(pq_sorted, pq),  # Q / vm at pq
                )
            ).reshape(-1,1)
