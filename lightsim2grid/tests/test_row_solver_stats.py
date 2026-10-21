# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""What a batch tells of each row's solve: the linear solver's work (get_row_linear_solver_stats),
the size of the Jacobian (get_row_J_dim / get_row_J_nnz) and the decisions the outer loops acted
on (OuterRowResult.decisions_taken); and the pivot-growth threshold of the refactor fallback,
part of the algorithm's config, so that a batch uses the grid's."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import test_outer_batch as tob
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

KLU = tob.OUTER == "NROuter_KLU"
# the first line of case14 cut off, these generators cannot hold the voltage: the rows that converge
LINES = range(1, tob.N_ROWS + 1)
COUNTERS = ("nb_reset", "nb_analyze", "nb_factorize", "nb_refactorize", "nb_refactorize_failed",
            "nb_fallback_factorize", "nb_fallback_factorize_failed", "nb_solve")


def _set_rgrowth_threshold(grid, val):
    cfg = grid.get_ac_algo_config()
    rp = list(cfg.real_params)
    rp[6] = val  # after the Newton's 6
    cfg.real_params = rp
    grid.set_ac_algo_config(cfg)
    return grid


class TestRowSolverStats(unittest.TestCase):
    def _batch(self, grid=None, nb_thread=1):
        grid = tob._grid() if grid is None else grid
        batch = ContingencyAnalysisCPP(grid)
        for line_id in LINES:
            batch.add_n1(line_id)
        batch.nb_thread = nb_thread
        batch.compute(tob._v0(grid), tob.MAX_IT, tob.TOL)
        self.assertEqual(batch.nb_converged(), tob.N_ROWS)
        return batch

    def _n_solve_stats(self):
        """the linear solver's work for the "n" solve the batch starts with"""
        grid = tob._grid()
        grid.ac_pf(tob._v0(grid), tob.MAX_IT, tob.TOL)
        return grid.get_algo().get_linear_solver_stats(), grid.get_algo().get_J()

    def _check_rows_add_up(self, batch, nb_worker=0):
        """the rows' counters and the "n" solve's are the whole compute()'s, but for the reset
        each worker's algorithm gets when it is made, before any row"""
        n_stats, _ = self._n_solve_stats()
        total = batch.get_linear_solver_stats()
        rows = batch.get_row_linear_solver_stats()
        self.assertEqual(len(rows), tob.N_ROWS)
        for name in COUNTERS:
            outside = nb_worker if name == "nb_reset" else 0
            self.assertEqual(sum(getattr(r, name) for r in rows) + getattr(n_stats, name) + outside,
                             getattr(total, name), name)

    def test_the_counters_are_each_row_own(self):
        batch = self._batch()
        self._check_rows_add_up(batch)
        for row, outer in zip(batch.get_row_linear_solver_stats(), batch.get_outer_results()):
            # the pattern is the batch's: no row analyzes
            self.assertEqual(row.nb_analyze, 0)
            # one linear solve per Newton iteration, over every inner solve of the row
            self.assertEqual(row.nb_solve, outer.nb_nr_iterations)
            self.assertGreater(row.timer_solve, 0.)

    def test_two_threads(self):
        batch = self._batch(nb_thread=2)
        self._check_rows_add_up(batch, nb_worker=2)
        # each worker analyzes and factorizes once, in the first row it solves
        rows = batch.get_row_linear_solver_stats()
        self.assertEqual(sum(r.nb_analyze for r in rows), 2)
        self.assertEqual(sum(r.nb_factorize for r in rows), 2)

    def test_jacobian_size(self):
        _, J = self._n_solve_stats()
        batch = self._batch()
        np.testing.assert_array_equal(batch.get_row_J_dim(), np.full(tob.N_ROWS, J.shape[0]))
        np.testing.assert_array_equal(batch.get_row_J_nnz(), np.full(tob.N_ROWS, J.nnz))

    def test_decisions_taken(self):
        # the generators are short of reactive power: every row switches some of them PQ
        for outer in self._batch().get_outer_results():
            taken = dict(outer.decisions_taken)
            self.assertEqual(list(taken), sorted(taken))
            self.assertGreater(taken.get("ReactiveLimits/PV_TO_PQ", 0), 0)
            self.assertTrue(all(nb > 0 for nb in taken.values()))

    @unittest.skipUnless(KLU, "KLU is not available")
    def test_min_rgrowth(self):
        # the outer loops keep the refactor fallback on: KLU measures the growth it checks
        for row in self._batch().get_row_linear_solver_stats():
            self.assertTrue(np.isfinite(row.min_rgrowth))
            self.assertGreater(row.min_rgrowth, 1e-10)

    @unittest.skipUnless(KLU, "KLU is not available")
    def test_the_batch_uses_the_grid_threshold(self):
        ref = self._batch()
        # above any growth: every refactorization falls back, for the same result
        batch = self._batch(_set_rgrowth_threshold(tob._grid(), 1e300))
        for row in batch.get_row_linear_solver_stats():
            self.assertGreater(row.nb_refactorize, 0)
            self.assertEqual(row.nb_refactorize_failed, row.nb_refactorize)
            self.assertEqual(row.nb_fallback_factorize, row.nb_refactorize)
        np.testing.assert_allclose(batch.get_voltages(), ref.get_voltages(), atol=1e-10, rtol=0.)
        # 0: nothing checked, so nothing measured
        batch = self._batch(_set_rgrowth_threshold(tob._grid(), 0.))
        for row in batch.get_row_linear_solver_stats():
            self.assertTrue(np.isnan(row.min_rgrowth))


if __name__ == "__main__":
    unittest.main()
