# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import unittest
import pickle

import numpy as np
from lightsim2grid.algorithm import (
    NR_KLU,
    NRSing_KLU
)


class JacobianMultiSlackTester(unittest.TestCase):
    @classmethod
    def get_name(cls):
        return os.path.join("jacobian_case14_sandbox","saved_jacobian_multi.pkl")
    
    def get_algo(self):
        return NR_KLU()
    
    @classmethod
    def setUpClass(cls) -> None:
        with open(cls.get_name(), "rb") as f:
            cls.res = pickle.load(f)
            
        return super().setUpClass()
    
    def setUp(self) -> None:
        self.nr_algo = self.get_algo()
        return super().setUp()    
    
    def new_to_old_indexes(self):
        """
        in lightsim2grid 0.14, order of the Jac change and slack is put at the end (instead of the beginning)
        
        this should be used J_new_indx[self.reorder_jac_row().T, self.reorder_jac_row()] == J_old_index
        
        new order:
        pv, pq, pq, pv_slack, slack
        
        old order:
        slack, pv_(slack + base), pq, pq
        
        """
        return np.concatenate(
            (
                [22],  # slack
                [20, 21],  # pv_slack
                np.arange(20)  # rest (pv_base, pq, pq)
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


class JacobianSingleSlackTester(JacobianMultiSlackTester):
    @classmethod
    def get_name(cls):
        return os.path.join("jacobian_case14_sandbox","saved_jacobian_single.pkl")
    
    def get_algo(self):
        return NRSing_KLU()
    
    def new_to_old_indexes(self):
        """
        No changes for single slack atm
        """
        return np.arange(22).reshape(-1,1)
