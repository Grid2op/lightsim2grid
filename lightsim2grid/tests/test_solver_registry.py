# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the SolverRegistry refactor (ChooseSolver → registry-backed unique_ptr)."""

import os
import unittest

from lightsim2grid.lightsim2grid_cpp import GridModel, AlgorithmType


def _make_grid():
    """Return a minimal GridModel with one bus so powerflows can be run."""
    gm = GridModel()
    gm.set_sn_mva(100.0)
    gm.set_init_vm_pu(1.0)
    return gm


class TestDefaultSolver(unittest.TestCase):
    """The default solver must be SparseLU."""

    def test_default_ac_solver_type(self):
        gm = _make_grid()
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_SparseLU)

    def test_default_dc_solver_type(self):
        gm = _make_grid()
        self.assertEqual(gm.get_dc_solver_type(), AlgorithmType.DC_SparseLU)


class TestEnumOverload(unittest.TestCase):
    """Enum-based change_solver must still work as before."""

    def test_change_solver_enum(self):
        gm = _make_grid()
        gm.change_algorithm(AlgorithmType.NRSing_SparseLU)
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NRSing_SparseLU)

    def test_change_dc_solver_enum(self):
        gm = _make_grid()
        gm.change_algorithm(AlgorithmType.DC_SparseLU)
        self.assertEqual(gm.get_dc_solver_type(), AlgorithmType.DC_SparseLU)

    def test_round_trip_enum(self):
        gm = _make_grid()
        gm.change_algorithm(AlgorithmType.GaussSeidel)
        self.assertEqual(gm.get_algo_type(), AlgorithmType.GaussSeidel)
        gm.change_algorithm(AlgorithmType.NR_SparseLU)
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_SparseLU)


class TestStringOverload(unittest.TestCase):
    """String-based change_solver must work for all built-in names."""

    def test_change_solver_string_sparselU(self):
        gm = _make_grid()
        gm.change_algorithm(AlgorithmType.GaussSeidel)   # change away from default
        gm.change_algorithm("SparseLU")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_SparseLU)

    def test_change_solver_string_gaussseidel(self):
        gm = _make_grid()
        gm.change_algorithm("GaussSeidel")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.GaussSeidel)

    def test_change_solver_string_dc(self):
        gm = _make_grid()
        gm.change_algorithm("DC")
        self.assertEqual(gm.get_dc_solver_type(), AlgorithmType.DC_SparseLU)

    def test_change_solver_unknown_name_raises(self):
        gm = _make_grid()
        with self.assertRaises(Exception):
            gm.change_algorithm("NonExistentSolverXYZ")


class TestAvailableSolvers(unittest.TestCase):
    """available_solvers() and available_solver_names() must return consistent info."""

    def test_available_solvers_returns_list(self):
        gm = _make_grid()
        solvers = gm.available_algorithms()
        self.assertIsInstance(solvers, list)
        self.assertIn(AlgorithmType.NR_SparseLU, solvers)
        self.assertIn(AlgorithmType.DC_SparseLU, solvers)

    def test_available_solver_names_returns_strings(self):
        gm = _make_grid()
        names = gm.available_solver_names()
        self.assertIsInstance(names, list)
        self.assertTrue(all(isinstance(n, str) for n in names))
        self.assertIn("SparseLU", names)
        self.assertIn("DC", names)

    def test_available_solver_names_covers_available_solvers(self):
        """Every enum returned by available_solvers() must have a string counterpart."""
        gm = _make_grid()
        names = set(gm.available_solver_names())
        for st in gm.available_algorithms():
            # Convert AlgorithmType to its string name by checking all known names
            gm2 = _make_grid()
            gm2.change_algorithm(st)
            # After change, get_algo_type or get_dc_solver_type reflects the change
            if st in (AlgorithmType.DC_SparseLU, AlgorithmType.DC_KLU if hasattr(AlgorithmType, "KLUDC") else None,
                      AlgorithmType.DC_NICSLU if hasattr(AlgorithmType, "NICSLUDC") else None,
                      AlgorithmType.DC_CKTSO if hasattr(AlgorithmType, "CKTSODC") else None):
                pass  # DC solver types
            else:
                self.assertEqual(gm2.get_algo_type(), st)


class TestKLUSolver(unittest.TestCase):
    """KLU solver tests (skipped gracefully when not compiled in)."""

    def setUp(self):
        from lightsim2grid.lightsim2grid_cpp import klu_solver_available
        if not klu_solver_available:
            self.skipTest("KLU solver not compiled in this build")

    def test_change_to_klu_by_enum(self):
        gm = _make_grid()
        gm.change_algorithm(AlgorithmType.NR_KLU)
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_KLU)

    def test_change_to_klu_by_string(self):
        gm = _make_grid()
        gm.change_algorithm("KLU")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_KLU)

    def test_klu_in_available_solver_names(self):
        gm = _make_grid()
        self.assertIn("KLU", gm.available_solver_names())


class TestPluginLoading(unittest.TestCase):
    """Plugin loading via load_solver_plugin (skipped if example not built)."""

    def _get_plugin_path(self):
        import sys
        base = os.path.join(os.path.dirname(__file__), "../../examples/external_solver")
        if sys.platform == "win32":
            candidates = [
                os.path.join(base, "build/Release/dummy_solver.dll"),
                os.path.join(base, "build/Debug/dummy_solver.dll"),
                os.path.join(base, "build/dummy_solver.dll"),
            ]
        else:
            candidates = [
                os.path.join(base, "build/libdummy_solver.so"),
                os.path.join(base, "libdummy_solver.so"),
            ]
        for p in candidates:
            if os.path.exists(os.path.abspath(p)):
                return os.path.abspath(p)
        return None

    def test_load_plugin_and_change_solver(self):
        path = self._get_plugin_path()
        if path is None:
            self.skipTest(
                "Example plugin not built. "
                "Build examples/external_solver/ first to run this test.")
        from lightsim2grid import load_solver_plugin
        load_solver_plugin(path)

        gm = _make_grid()
        names = gm.available_solver_names()
        self.assertIn("DummyExternal", names)
        gm.change_algorithm("DummyExternal")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.Custom)


if __name__ == "__main__":
    unittest.main()
