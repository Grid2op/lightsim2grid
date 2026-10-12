# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
How the sparse matrices reach python: none of them is copied on the way.

A matrix a getter builds for the call (get_Ybus, get_dcYbus, get_Bf, get_Bf_solver, and
the per-row get_Ybus / get_dcYbus of the batch classes) is handed over to python: its
arrays are the C++ ones, owned by a capsule that frees them with the last array using
them. Each call gives a matrix of its own, writeable, independent of the grid.
"""
import gc
import unittest
import warnings

import numpy as np
import scipy.sparse as sp
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.algorithm import AlgorithmType

from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


def _make_grid():
    if MAX_PP2_DATAREADER < CURRENT_PP_VERSION:
        loader_kwargs = {"pp_orig_file": "pandapower_v3"}
    else:
        loader_kwargs = {"pp_orig_file": "pandapower_v2"}
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case118(), **loader_kwargs)
    V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 10, 1e-8)
    assert V.shape[0] > 0, "AC powerflow did not converge"
    V = grid.dc_pf(np.ones(grid.total_bus(), dtype=complex), 1, 1e-8)
    assert V.shape[0] > 0, "DC powerflow did not converge"
    return grid


def _arrays(mat):
    return mat.data, mat.indices, mat.indptr


def _root_owner(arr):
    """What actually holds the memory of `arr`, through the views scipy may take of it
    (it slices `data` and `indices` to nnz). None if one of the arrays owns its data,
    which means it is a copy."""
    while isinstance(arr, np.ndarray):
        if arr.flags.owndata:
            return None
        arr = arr.base
    return arr


class TestSparseGettersHandedOver(unittest.TestCase):
    def setUp(self):
        self.grid = _make_grid()

    def _assert_handed_over(self, mat):
        """`mat` is a scipy matrix on C++ buffers, owned by one capsule, writeable."""
        self.assertIsInstance(mat, sp.csc_matrix)
        self.assertTrue(mat.has_canonical_format)
        # a copy (by pybind11, or by scipy) would own its data
        owner = _root_owner(mat.data)
        self.assertEqual(type(owner).__name__, "PyCapsule")
        for arr in _arrays(mat):
            self.assertIs(_root_owner(arr), owner)
            self.assertTrue(arr.flags.writeable)

    def _assert_same(self, mat, ref):
        self.assertEqual(mat.shape, ref.shape)
        for arr, arr_ref in zip(_arrays(mat), _arrays(ref)):
            np.testing.assert_array_equal(arr, arr_ref)

    def _getters(self):
        return {"get_Ybus": self.grid.get_Ybus,
                "get_dcYbus": self.grid.get_dcYbus,
                "get_Bf": self.grid.get_Bf,
                "get_Bf_solver": self.grid.get_Bf_solver}

    def test_not_copied(self):
        for nm, getter in self._getters().items():
            with self.subTest(getter=nm):
                self._assert_handed_over(getter())

    def test_values(self):
        """What reaches python is the relabelled solver matrix, entry for entry."""
        id_ac = np.array(self.grid.id_ac_solver_to_me())
        Ybus = self.grid.get_Ybus()
        self._assert_same(Ybus[id_ac.reshape(-1, 1), id_ac.reshape(1, -1)].tocsc(),
                          self.grid.get_Ybus_solver())
        id_dc = np.array(self.grid.id_dc_solver_to_me())
        dcYbus = self.grid.get_dcYbus()
        self._assert_same(dcYbus[id_dc.reshape(-1, 1), id_dc.reshape(1, -1)].tocsc(),
                          self.grid.get_dcYbus_solver())
        Bf = self.grid.get_Bf()
        self._assert_same(Bf[:, id_dc].tocsc(), self.grid.get_Bf_solver())

    def test_calls_are_independent(self):
        for nm, getter in self._getters().items():
            with self.subTest(getter=nm):
                first = getter()
                ref = first.copy()
                second = getter()
                for arr, arr2 in zip(_arrays(first), _arrays(second)):
                    self.assertFalse(np.shares_memory(arr, arr2))
                # writing into one changes neither the other nor what the grid gives next
                first.data[:] = 0.
                self._assert_same(second, ref)
                self._assert_same(getter(), ref)

    def test_outlives_the_grid(self):
        mats = {nm: getter() for nm, getter in self._getters().items()}
        refs = {nm: mat.copy() for nm, mat in mats.items()}
        del self.grid
        gc.collect()
        for nm, mat in mats.items():
            with self.subTest(getter=nm):
                self._assert_same(mat, refs[nm])

    def test_outlives_its_scipy_matrix(self):
        """The arrays keep the buffers alive on their own, without the scipy object."""
        data, indices, indptr = _arrays(self.grid.get_Ybus())
        ref = self.grid.get_Ybus().copy()
        gc.collect()
        np.testing.assert_array_equal(data, ref.data)
        np.testing.assert_array_equal(indices, ref.indices)
        np.testing.assert_array_equal(indptr, ref.indptr)


class TestBatchYbusHandedOver(TestSparseGettersHandedOver):
    """The per-row matrices of ContingencyAnalysisCPP (ScenarioSweepCPP shares the binding)."""

    def setUp(self):
        self.grid = _make_grid()
        self.batch = {}
        for ac, algo in ((True, AlgorithmType.NR_SparseLU), (False, AlgorithmType.DC_SparseLU)):
            ca = ContingencyAnalysisCPP(self.grid)
            ca.change_algorithm(algo)
            ca.add_n1(3)
            ca.add_n1(7)
            ca.compute(np.ones(self.grid.total_bus(), dtype=complex), 30, 1e-10)
            self.batch[ac] = ca

    def _getters(self):
        res = {}
        for row in range(2):
            for solver_numbering in (False, True):
                key = f"row {row}, solver_numbering={solver_numbering}"
                res[f"get_Ybus {key}"] = (lambda r=row, s=solver_numbering:
                                          self.batch[True].get_Ybus(r, solver_numbering=s))
                res[f"get_dcYbus {key}"] = (lambda r=row, s=solver_numbering:
                                            self.batch[False].get_dcYbus(r, solver_numbering=s))
        return res

    def test_values(self):
        """Each row removes its own branch: two rows give two different matrices."""
        for ac in (True, False):
            getter = self.batch[True].get_Ybus if ac else self.batch[False].get_dcYbus
            self.assertNotEqual((getter(0) != getter(1)).nnz, 0)

    def test_outlives_the_grid(self):
        mats = {nm: getter() for nm, getter in self._getters().items()}
        refs = {nm: mat.copy() for nm, mat in mats.items()}
        del self.batch
        del self.grid
        gc.collect()
        for nm, mat in mats.items():
            with self.subTest(getter=nm):
                self._assert_same(mat, refs[nm])


if __name__ == "__main__":
    unittest.main()
