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
import scipy.sparse as sp
import pandapower as pp
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


def _make_gridmodel_with_gaps():
    """case14 with three disconnected buses, numbered first, in the middle and last.

    The solver then covers a strict subset of the grid buses, with a gap at each end and
    one in between, so `id_solver_to_me` is neither the identity nor a contiguous range.
    """
    pp_net = pn.case14()
    nb_bus = pp_net.bus.shape[0]
    for _ in range(3):
        pp.create_bus(pp_net, vn_kv=pp_net.bus["vn_kv"].iloc[0], in_service=False)
    old_ids = list(pp_net.bus.index)
    new_order = [nb_bus] + old_ids[:7] + [nb_bus + 1] + old_ids[7:nb_bus] + [nb_bus + 2]
    pp.toolbox.reindex_buses(pp_net, {old: new for new, old in enumerate(new_order)})
    pp_net.bus.sort_index(inplace=True)
    gridmodel = _make_gridmodel(pp_net)
    for bus_id in np.where(~pp_net.bus["in_service"].values)[0]:
        gridmodel.deactivate_bus(int(bus_id))
    return gridmodel


class TestRelabelMatrixDisconnectedBuses(unittest.TestCase):
    """
    get_Ybus / get_dcYbus / get_Bf when some grid buses are not in the solver: the result
    must be exactly what a relabelling of the solver matrix gives, empty rows and columns
    for the disconnected buses included, and in canonical CSC form (sorted, no duplicate).
    """

    def _relabel(self, mat_solver, id_solver_to_me, total_bus, relabel_row=True):
        coo = mat_solver.tocoo()
        rows = id_solver_to_me[coo.row] if relabel_row else coo.row
        shape = (total_bus if relabel_row else mat_solver.shape[0], total_bus)
        return sp.csc_matrix((coo.data, (rows, id_solver_to_me[coo.col])), shape=shape)

    def _assert_same_csc(self, res, expected):
        self.assertEqual(res.shape, expected.shape)
        self.assertTrue(res.has_canonical_format)
        expected = expected.copy()
        expected.sum_duplicates()
        # array equality, not a numerical tolerance: relabelling moves the values unchanged
        np.testing.assert_array_equal(res.indptr, expected.indptr)
        np.testing.assert_array_equal(res.indices, expected.indices)
        np.testing.assert_array_equal(res.data, expected.data)

    def test_mapping_has_gaps(self):
        gridmodel = _make_gridmodel_with_gaps()
        V = gridmodel.ac_pf(np.ones(gridmodel.total_bus(), dtype=complex), 10, 1e-8)
        self.assertTrue(V.shape[0] > 0)
        id_ac = np.array(gridmodel.id_ac_solver_to_me())
        self.assertEqual(gridmodel.total_bus(), 17)
        self.assertEqual(id_ac.shape[0], 14)
        self.assertNotIn(0, id_ac)
        self.assertNotIn(8, id_ac)
        self.assertNotIn(16, id_ac)

    def test_ac_ybus(self):
        gridmodel = _make_gridmodel_with_gaps()
        V = gridmodel.ac_pf(np.ones(gridmodel.total_bus(), dtype=complex), 10, 1e-8)
        self.assertTrue(V.shape[0] > 0)
        id_ac = np.array(gridmodel.id_ac_solver_to_me())
        expected = self._relabel(gridmodel.get_Ybus_solver(), id_ac, gridmodel.total_bus())
        self._assert_same_csc(gridmodel.get_Ybus(), expected)

    def test_dc_ybus_and_bf(self):
        gridmodel = _make_gridmodel_with_gaps()
        V = gridmodel.dc_pf(np.ones(gridmodel.total_bus(), dtype=complex), 1, 1e-8)
        self.assertTrue(V.shape[0] > 0)
        id_dc = np.array(gridmodel.id_dc_solver_to_me())
        expected = self._relabel(gridmodel.get_dcYbus_solver(), id_dc, gridmodel.total_bus())
        self._assert_same_csc(gridmodel.get_dcYbus(), expected)
        # Bf: one row per branch, only its columns are relabelled
        expected = self._relabel(gridmodel.get_Bf_solver(), id_dc, gridmodel.total_bus(), relabel_row=False)
        self._assert_same_csc(gridmodel.get_Bf(), expected)

    def test_after_topology_change(self):
        """A second bus disconnected after a first solve: the new mapping is used."""
        gridmodel = _make_gridmodel_with_gaps()
        V = gridmodel.ac_pf(np.ones(gridmodel.total_bus(), dtype=complex), 10, 1e-8)
        self.assertTrue(V.shape[0] > 0)
        # grid bus 9 is bus 7 of case14: a generator at the end of a transformer, with
        # nothing else on it, so it can be disconnected without splitting the grid
        gen_id = [gen.id for gen in gridmodel.get_generators() if gen.bus_id == 9]
        self.assertEqual(len(gen_id), 1)
        gridmodel.deactivate_gen(gen_id[0])
        trafo_id = [tr.id for tr in gridmodel.get_trafos() if 9 in (tr.bus1_id, tr.bus2_id)]
        for tr_id in trafo_id:
            gridmodel.deactivate_trafo(tr_id)
        line_id = [li.id for li in gridmodel.get_lines() if 9 in (li.bus1_id, li.bus2_id)]
        for li_id in line_id:
            gridmodel.deactivate_powerline(li_id)
        gridmodel.deactivate_bus(9)
        V = gridmodel.ac_pf(np.ones(gridmodel.total_bus(), dtype=complex), 10, 1e-8)
        self.assertTrue(V.shape[0] > 0)
        id_ac = np.array(gridmodel.id_ac_solver_to_me())
        self.assertEqual(id_ac.shape[0], 13)
        self.assertNotIn(9, id_ac)
        expected = self._relabel(gridmodel.get_Ybus_solver(), id_ac, gridmodel.total_bus())
        self._assert_same_csc(gridmodel.get_Ybus(), expected)


if __name__ == "__main__":
    unittest.main()
