# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Tests that _relabel_matrix (and _relabel_vector) correctly convert Ybus / Sbus
from "solver bus id" labelling to "gridmodel bus id" labelling, accessible via:
  - get_Ybus()       / get_Ybus_solver()       (AC)
  - get_dcYbus()     / get_dcYbus_solver()      (DC)
  - get_Sbus()       / get_Sbus_solver()        (AC)
  - get_dcSbus()     / get_dcSbus_solver()      (DC)

The mapping between the two labellings is given by id_ac_solver_to_me() and
id_dc_solver_to_me().
"""
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower

from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


def _make_gridmodel(pp_net):
    """Load a pandapower network into a lightsim2grid gridmodel and return it."""
    if MAX_PP2_DATAREADER < CURRENT_PP_VERSION:
        loader_kwargs = {"pp_orig_file": "pandapower_v3"}
    else:
        loader_kwargs = {"pp_orig_file": "pandapower_v2"}
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pp_net, **loader_kwargs)
    return grid


class TestRelabelMatrix(unittest.TestCase):
    """
    Verifies that get_Ybus() / get_dcYbus() are consistent with the corresponding
    solver-labelled matrices (get_Ybus_solver() / get_dcYbus_solver()) after
    applying the id_ac_solver_to_me / id_dc_solver_to_me re-mapping.
    """

    def _check_ac(self, gridmodel):
        id_ac = np.array(gridmodel.id_ac_solver_to_me())

        Ybus_solver = gridmodel.get_Ybus_solver()
        Ybus = gridmodel.get_Ybus()

        nb_conn = gridmodel.nb_connected_bus()
        total = gridmodel.total_bus()

        self.assertEqual(Ybus_solver.shape, (nb_conn, nb_conn))
        self.assertEqual(Ybus.shape, (total, total))
        self.assertEqual(Ybus.nnz, Ybus_solver.nnz,
                         "nnz mismatch between get_Ybus and get_Ybus_solver")

        diff = Ybus[id_ac.reshape(-1, 1), id_ac.reshape(1, -1)] - Ybus_solver
        self.assertEqual(diff.nnz, 0,
                         f"get_Ybus and get_Ybus_solver differ after relabelling: "
                         f"max abs diff = {np.abs(diff).max():.3e}")

    def _check_ac_sbus(self, gridmodel):
        id_ac = np.array(gridmodel.id_ac_solver_to_me())

        Sbus_solver = gridmodel.get_Sbus_solver()
        Sbus = gridmodel.get_Sbus()

        self.assertEqual(Sbus_solver.shape, (gridmodel.nb_connected_bus(),))
        self.assertEqual(Sbus.shape, (gridmodel.total_bus(),))
        np.testing.assert_array_equal(Sbus[id_ac], Sbus_solver)

    def _check_dc(self, gridmodel):
        id_dc = np.array(gridmodel.id_dc_solver_to_me())

        dcYbus_solver = gridmodel.get_dcYbus_solver()
        dcYbus = gridmodel.get_dcYbus()

        nb_conn = gridmodel.nb_connected_bus()
        total = gridmodel.total_bus()

        self.assertEqual(dcYbus_solver.shape, (nb_conn, nb_conn))
        self.assertEqual(dcYbus.shape, (total, total))
        self.assertEqual(dcYbus.nnz, dcYbus_solver.nnz,
                         "nnz mismatch between get_dcYbus and get_dcYbus_solver")

        diff = dcYbus[id_dc.reshape(-1, 1), id_dc.reshape(1, -1)] - dcYbus_solver
        self.assertEqual(diff.nnz, 0,
                         f"get_dcYbus and get_dcYbus_solver differ after relabelling: "
                         f"max abs diff = {np.abs(diff).max():.3e}")

    def _check_dc_sbus(self, gridmodel):
        id_dc = np.array(gridmodel.id_dc_solver_to_me())

        dcSbus_solver = gridmodel.get_dcSbus_solver()
        dcSbus = gridmodel.get_dcSbus()

        self.assertEqual(dcSbus_solver.shape, (gridmodel.nb_connected_bus(),))
        self.assertEqual(dcSbus.shape, (gridmodel.total_bus(),))
        np.testing.assert_array_equal(dcSbus[id_dc], dcSbus_solver)

    # ------------------------------------------------------------------
    # AC powerflow tests
    # ------------------------------------------------------------------

    def _run_ac_checks(self, pp_net):
        gridmodel = _make_gridmodel(pp_net)
        V = gridmodel.ac_pf(np.ones(gridmodel.get_bus_vn_kv().shape[0], dtype=complex), 10, 1e-7)
        self.assertTrue(V.shape[0] > 0, f"AC powerflow did not converge: {gridmodel.get_solver().get_error()}")
        self._check_ac(gridmodel)
        self._check_ac_sbus(gridmodel)

    def test_ac_case14(self):
        self._run_ac_checks(pn.case14())

    def test_ac_case118(self):
        self._run_ac_checks(pn.case118())

    # ------------------------------------------------------------------
    # DC powerflow tests
    # ------------------------------------------------------------------

    def _run_dc_checks(self, pp_net):
        gridmodel = _make_gridmodel(pp_net)
        V = gridmodel.dc_pf(np.ones(gridmodel.get_bus_vn_kv().shape[0], dtype=complex), 1, 1e-5)
        self.assertTrue(V.shape[0] > 0, f"DC powerflow did not converge: {gridmodel.get_dc_solver().get_error()}")
        self._check_dc(gridmodel)
        self._check_dc_sbus(gridmodel)

    def test_dc_case14(self):
        self._run_dc_checks(pn.case14())

    def test_dc_case118(self):
        self._run_dc_checks(pn.case118())

    # ------------------------------------------------------------------
    # After AC-only / DC-only run: cross-accessor behaviour
    # ------------------------------------------------------------------

    def test_ac_only_raises_for_dc_ybus(self):
        """After a pure AC powerflow, get_dcYbus should raise RuntimeError."""
        gridmodel = _make_gridmodel(pn.case14())
        V = gridmodel.ac_pf(np.ones(gridmodel.get_bus_vn_kv().shape[0], dtype=complex), 10, 1e-7)
        self.assertTrue(V.shape[0] > 0, f"AC powerflow did not converge: {gridmodel.get_solver().get_error()}")

        dcYbus_solver = gridmodel.get_dcYbus_solver()
        self.assertEqual(dcYbus_solver.shape, (0, 0))
        with self.assertRaises(RuntimeError):
            _ = gridmodel.get_dcYbus()

    def test_dc_only_raises_for_ac_ybus(self):
        """After a pure DC powerflow, get_Ybus should raise RuntimeError."""
        gridmodel = _make_gridmodel(pn.case14())
        V = gridmodel.dc_pf(np.ones(gridmodel.get_bus_vn_kv().shape[0], dtype=complex), 1, 1e-5)
        self.assertTrue(V.shape[0] > 0, f"DC powerflow did not converge: {gridmodel.get_dc_solver().get_error()}")

        Ybus_solver = gridmodel.get_Ybus_solver()
        self.assertEqual(Ybus_solver.shape, (0, 0))
        with self.assertRaises(RuntimeError):
            _ = gridmodel.get_Ybus()


if __name__ == "__main__":
    unittest.main()
