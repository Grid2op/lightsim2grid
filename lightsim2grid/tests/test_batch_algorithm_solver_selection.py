# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for name-based (registry) algorithm selection on the batch containers
(TimeSeriesCPP / ContingencyAnalysisCPP), mirroring what LSGrid already supports
(see test_solver_registry.py)."""

import unittest

from lightsim2grid.lightsim2grid_cpp import (LSGrid,
                                              AlgorithmType,
                                              TimeSeriesCPP,
                                              ContingencyAnalysisCPP,
                                              klu_solver_available)


def _make_grid():
    """Return a minimal LSGrid with one bus so it can back a batch container."""
    gm = LSGrid()
    gm.set_sn_mva(100.0)
    gm.set_init_vm_pu(1.0)
    return gm


class TestTimeSeriesSolverSelection(unittest.TestCase):
    def test_default_algo_name_matches_grid(self):
        ts = TimeSeriesCPP(_make_grid())
        self.assertEqual(ts.get_algo_type(), AlgorithmType.NR_SparseLU)
        self.assertEqual(ts.get_algo_name(), "NR_SparseLU")

    def test_change_algorithm_enum(self):
        ts = TimeSeriesCPP(_make_grid())
        ts.change_algorithm(AlgorithmType.GaussSeidel)
        self.assertEqual(ts.get_algo_type(), AlgorithmType.GaussSeidel)
        self.assertEqual(ts.get_algo_name(), "GaussSeidel")

    def test_change_algorithm_string(self):
        ts = TimeSeriesCPP(_make_grid())
        ts.change_algorithm("GaussSeidel")
        self.assertEqual(ts.get_algo_type(), AlgorithmType.GaussSeidel)
        self.assertEqual(ts.get_algo_name(), "GaussSeidel")
        ts.change_algorithm("NR_SparseLU")
        self.assertEqual(ts.get_algo_type(), AlgorithmType.NR_SparseLU)
        self.assertEqual(ts.get_algo_name(), "NR_SparseLU")

    def test_change_algorithm_unknown_name_raises(self):
        ts = TimeSeriesCPP(_make_grid())
        with self.assertRaises(Exception):
            ts.change_algorithm("NonExistentSolverXYZ")

    def test_available_algorithm_names(self):
        ts = TimeSeriesCPP(_make_grid())
        names = ts.available_algorithm_names()
        self.assertIsInstance(names, list)
        self.assertIn("NR_SparseLU", names)
        self.assertIn("DC_SparseLU", names)

    @unittest.skipUnless(klu_solver_available, "KLU solver not compiled in this build")
    def test_inherits_no_enum_builtin_solver_by_name(self):
        """A grid running a built-in with no dedicated AlgorithmType member (eg
        NRRefactorRetry_KLU, which the enum collapses onto AlgorithmType.Custom)
        must still be usable to construct a TimeSeries: before the fix this
        raised because the constructor propagated the algo via the (lossy)
        enum instead of by registry name."""
        gm = _make_grid()
        gm.change_algorithm("NRRefactorRetry_KLU")
        ts = TimeSeriesCPP(gm)
        self.assertEqual(ts.get_algo_name(), "NRRefactorRetry_KLU")
        self.assertEqual(ts.get_algo_type(), AlgorithmType.Custom)


class TestContingencyAnalysisSolverSelection(unittest.TestCase):
    def test_default_algo_name_matches_grid(self):
        ca = ContingencyAnalysisCPP(_make_grid(), False)
        self.assertEqual(ca.get_algo_type(), AlgorithmType.NR_SparseLU)
        self.assertEqual(ca.get_algo_name(), "NR_SparseLU")

    def test_change_algorithm_enum(self):
        ca = ContingencyAnalysisCPP(_make_grid(), False)
        ca.change_algorithm(AlgorithmType.GaussSeidel)
        self.assertEqual(ca.get_algo_type(), AlgorithmType.GaussSeidel)
        self.assertEqual(ca.get_algo_name(), "GaussSeidel")

    def test_change_algorithm_string(self):
        ca = ContingencyAnalysisCPP(_make_grid(), False)
        ca.change_algorithm("GaussSeidel")
        self.assertEqual(ca.get_algo_type(), AlgorithmType.GaussSeidel)
        self.assertEqual(ca.get_algo_name(), "GaussSeidel")
        ca.change_algorithm("NR_SparseLU")
        self.assertEqual(ca.get_algo_type(), AlgorithmType.NR_SparseLU)
        self.assertEqual(ca.get_algo_name(), "NR_SparseLU")

    def test_change_algorithm_unknown_name_raises(self):
        ca = ContingencyAnalysisCPP(_make_grid(), False)
        with self.assertRaises(Exception):
            ca.change_algorithm("NonExistentSolverXYZ")

    def test_available_algorithm_names(self):
        ca = ContingencyAnalysisCPP(_make_grid(), False)
        names = ca.available_algorithm_names()
        self.assertIsInstance(names, list)
        self.assertIn("NR_SparseLU", names)
        self.assertIn("DC_SparseLU", names)

    @unittest.skipUnless(klu_solver_available, "KLU solver not compiled in this build")
    def test_inherits_no_enum_builtin_solver_by_name(self):
        """Same as TestTimeSeriesSolverSelection's counterpart: also exercises
        the multi-threaded per-thread solver rebuild (nb_thread > 1), which
        used to propagate the algorithm via the lossy enum too."""
        gm = _make_grid()
        gm.change_algorithm("NRRefactorRetry_KLU")
        ca = ContingencyAnalysisCPP(gm, False)
        self.assertEqual(ca.get_algo_name(), "NRRefactorRetry_KLU")
        self.assertEqual(ca.get_algo_type(), AlgorithmType.Custom)
        ca.nb_thread = 2


if __name__ == "__main__":
    unittest.main()
