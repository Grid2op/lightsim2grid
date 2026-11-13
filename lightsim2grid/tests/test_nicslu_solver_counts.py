# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""LinearSolverStats of the *_NICSLU algorithms: NICSLU chooses between a factorization and a
refactorization itself, and factorizes when a refactorization fails. The stats count what it
ran: a refactorization that ended in a factorization is a failed one and a fallback, as with
KLU, where they used to show as plain refactorizations."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import test_outer_batch as tob
from lightsim2grid import compilation_options
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP


@unittest.skipUnless(compilation_options.nicslu_solver_available and compilation_options.klu_solver_available,
                     "NICSLU or KLU is not available")
class TestNICSLUStats(unittest.TestCase):
    def setUp(self):
        try:
            import pypowsybl  # noqa: F401
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        self.addCleanup(setattr, tob, "OUTER", tob.OUTER)

    def _run(self, algo):
        # the transformer voltage control grid of test_constrained_ordering: a row switching
        # the tap's voltage row costs NROuter_KLU a fallback factorization
        tob.OUTER = algo
        grid = tob._tap_grid()
        batch = ContingencyAnalysisCPP(grid)
        for line_id in range(tob.N_ROWS):
            batch.add_n1(line_id)
        batch.compute(tob._v0(grid), tob.MAX_IT, tob.TOL)
        return batch

    def test_repivoting_shows_as_fallbacks(self):
        klu = self._run("NROuter_KLU")
        nicslu = self._run("NROuter_NICSLU")
        np.testing.assert_array_equal(nicslu.converged_mask(), klu.converged_mask())
        np.testing.assert_allclose(nicslu.get_voltages(), klu.get_voltages(), rtol=0., atol=1e-10)
        self.assertGreaterEqual(klu.get_linear_solver_stats().nb_fallback_factorize, 1)

        stats = nicslu.get_linear_solver_stats()
        self.assertGreaterEqual(stats.nb_fallback_factorize, 1)
        self.assertEqual(stats.nb_fallback_factorize_failed, 0)
        self.assertEqual(stats.nb_refactorize_failed, stats.nb_fallback_factorize)
        # the first factorization, plus one per fallback
        self.assertEqual(stats.nb_factorize, stats.nb_analyze + stats.nb_fallback_factorize)


if __name__ == "__main__":
    unittest.main()
