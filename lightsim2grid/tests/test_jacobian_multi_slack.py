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
from lightsim2grid.algorithm import NR_KLU


class JacobianTester(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with open(os.path.join("jacobian_case14_sandbox","saved_jacobian.pkl"), "rb") as f:
            cls.res = pickle.load(f)
            
        return super().setUpClass()
    
    def setUp(self) -> None:
        self.solver = NR_KLU()
        return super().setUp()    
    
    def new_to_old_indexes(self):
        """
        in lightsim2grid 0.14, order of the Jac change and slack is put at the end (instead of the beginning)
        
        this should be used J_new_indx[self.reorder_jac_row().T, self.reorder_jac_row()] == J_old_index
        """
        return np.concatenate(([22], np.arange(22))).reshape(-1,1)
    
    def _aux_test_iter(self, iter):
        cls = type(self)
        self.solver.compute_pf(
            cls.res["init_state"]["Ybus"],
            cls.res["init_state"]["v_init"],
            cls.res["init_state"]["Sbus"],
            cls.res["init_state"]["slack_ids"],
            cls.res["init_state"]["slack_weights"],
            cls.res["init_state"]["pv"],
            cls.res["init_state"]["pq"],
            iter,
            cls.res["init_state"]["tol"])
        
        J = self.solver.get_J()[self.new_to_old_indexes().T, self.new_to_old_indexes()]
        V = self.solver.get_V()
        assert J.shape == cls.res[f"{iter}"]["J"].shape, f"error for iter {iter}: J.shape = {J.shape} != {cls.res[f'{iter}']['J'].shape}"
        assert J.nnz == cls.res[f"{iter}"]["J"].nnz, f"error for iter {iter}: J.nnz = {J.nnz} != {cls.res[f'{iter}']['J'].nnz}"
        if J.shape[0] > 0:
            assert np.abs(J - cls.res[f"{iter}"]["J"]).max() <= 1e-6, f"error for iter {iter}: {np.abs(J - cls.res[f'{iter}']['J']).max()}"
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
