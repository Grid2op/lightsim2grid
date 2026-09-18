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

from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType


def _make_grid():
    """Return a minimal LSGrid with one bus so powerflows can be run."""
    gm = LSGrid()
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
        gm.change_algorithm("NR_SparseLU")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_SparseLU)

    def test_change_solver_string_gaussseidel(self):
        gm = _make_grid()
        gm.change_algorithm("GaussSeidel")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.GaussSeidel)

    def test_change_solver_string_dc(self):
        gm = _make_grid()
        gm.change_algorithm("DC_SparseLU")
        self.assertEqual(gm.get_dc_solver_type(), AlgorithmType.DC_SparseLU)

    def test_change_solver_unknown_name_raises(self):
        gm = _make_grid()
        with self.assertRaises(Exception):
            gm.change_algorithm("NonExistentSolverXYZ")


class TestAvailableSolvers(unittest.TestCase):
    """available_solvers() and available_solver_names() must return consistent info."""

    def test_available_solvers_returns_list(self):
        gm = _make_grid()
        solvers = gm.available_solvers()
        self.assertIsInstance(solvers, list)
        self.assertIn(AlgorithmType.NR_SparseLU, solvers)
        self.assertIn(AlgorithmType.DC_SparseLU, solvers)

    def test_available_solver_names_returns_strings(self):
        gm = _make_grid()
        names = gm.available_solver_names()
        self.assertIsInstance(names, list)
        self.assertTrue(all(isinstance(n, str) for n in names))
        self.assertIn("NR_SparseLU", names)
        self.assertIn("DC_SparseLU", names)

    def test_available_solver_names_covers_available_solvers(self):
        """Every enum returned by available_solvers() must have a string counterpart."""
        gm = _make_grid()
        names = set(gm.available_solver_names())
        for st in gm.available_solvers():
            self.assertIn(st.name, names)


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
        gm.change_algorithm("NR_KLU")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.NR_KLU)

    def test_klu_in_available_solver_names(self):
        gm = _make_grid()
        self.assertIn("NR_KLU", gm.available_solver_names())


class TestPluginLoading(unittest.TestCase):
    """Plugin loading via load_algorithm_plugin (skipped if example not built).

    The DummyExternal solver lives in examples/external_algorithm/. Point the
    test at the built plugin with either ``LS2G_TEST_PLUGIN`` (full path to the
    .so/.dll) or ``LS2G_TEST_PLUGIN_DIR`` (its build directory) -- needed when
    running against the *installed* package, whose __file__ is in site-packages
    and cannot reach the source-tree examples/ by relative path. Without either,
    a source-tree relative path is tried, and the test skips if nothing is found.
    """

    @staticmethod
    def _candidates_in(build_dir):
        import sys
        if sys.platform == "win32":
            names = ["Release/dummy_solver.dll", "Debug/dummy_solver.dll", "dummy_solver.dll"]
        else:
            names = ["libdummy_solver.so"]  # .so on macOS too (see the plugin CMakeLists)
        return [os.path.join(build_dir, n) for n in names]

    def _get_plugin_path(self):
        env_file = os.environ.get("LS2G_TEST_PLUGIN")
        if env_file and os.path.exists(env_file):
            return os.path.abspath(env_file)

        build_dirs = []
        env_dir = os.environ.get("LS2G_TEST_PLUGIN_DIR")
        if env_dir:
            build_dirs.append(env_dir)
        # source-tree fallback (in-repo dev run, not the installed package)
        build_dirs.append(os.path.join(
            os.path.dirname(__file__), "..", "..", "examples", "external_algorithm", "build"))

        for build_dir in build_dirs:
            for p in self._candidates_in(build_dir):
                if os.path.exists(os.path.abspath(p)):
                    return os.path.abspath(p)
        return None

    def test_load_plugin_and_change_solver(self):
        path = self._get_plugin_path()
        if path is None:
            self.skipTest(
                "Example plugin not built. Build examples/external_algorithm/ (or set "
                "LS2G_TEST_PLUGIN / LS2G_TEST_PLUGIN_DIR) to run this test.")
        from lightsim2grid import load_algorithm_plugin

        # Idempotent: another test module in the same process may already have
        # loaded it, and loading a plugin whose name is registered now raises.
        if "DummyExternal" not in _make_grid().available_algorithm_names():
            load_algorithm_plugin(path)

        gm = _make_grid()
        names = gm.available_algorithm_names()
        self.assertIn("DummyExternal", names)
        gm.change_algorithm("DummyExternal")
        self.assertEqual(gm.get_algo_type(), AlgorithmType.Custom)


if __name__ == "__main__":
    unittest.main()
