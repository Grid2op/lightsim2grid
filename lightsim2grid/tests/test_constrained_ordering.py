# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""constrained_ordering: the KLU-based linear solvers order the rows the outer loops switch by
value after the Vm columns that keep their pivot valid. On by default for *_KLUStaticPivot, off
for *_KLU; the 5th int parameter of the Newton's AlgoConfig, so that a grid copy and the batch
classes keep it. On the transformer voltage control grid of test_outer_batch, a contingency
analysis needs one fallback factorization without it and none with it, for the same results."""

import copy
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import test_outer_batch as tob
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

try:
    from lightsim2grid.algorithm import NR_KLU, NR_KLUStaticPivot  # noqa: F401
    HAS_KLU = True
except ImportError:
    HAS_KLU = False

ORDERING = 4  # its place among the Newton's int parameters


def _with_ordering(grid, flag):
    cfg = grid.get_ac_algo_config()
    ip = list(cfg.int_params)
    ip[ORDERING] = int(flag)
    cfg.int_params = ip
    grid.set_ac_algo_config(cfg)
    return grid


@unittest.skipUnless(HAS_KLU, "KLU is not available")
class TestConstrainedOrdering(unittest.TestCase):
    def setUp(self):
        try:
            import pypowsybl  # noqa: F401
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        self.addCleanup(setattr, tob, "OUTER", tob.OUTER)

    def _grid(self, algo, flag=None):
        tob.OUTER = algo
        grid = tob._tap_grid()
        return grid if flag is None else _with_ordering(grid, flag)

    def _run(self, grid):
        batch = ContingencyAnalysisCPP(grid)
        for line_id in range(tob.N_ROWS):
            batch.add_n1(line_id)
        batch.compute(tob._v0(grid), tob.MAX_IT, tob.TOL)
        # the tap moves somewhere, or the rows switch nothing
        self.assertTrue(any(r.ratio_tap for r in batch.get_outer_results()))
        return batch

    def test_defaults(self):
        self.assertEqual(self._grid("NROuter_KLU").get_ac_algo_config().int_params[ORDERING], 0)
        self.assertEqual(self._grid("NROuter_KLUStaticPivot").get_ac_algo_config().int_params[ORDERING], 1)
        self.assertFalse(NR_KLU().constrained_ordering)
        self.assertTrue(NR_KLUStaticPivot().constrained_ordering)

    def test_klu_opt_in(self):
        default = self._run(self._grid("NROuter_KLU"))
        ordered = self._run(self._grid("NROuter_KLU", True))
        self.assertEqual(default.get_linear_solver_stats().nb_fallback_factorize, 1)
        self.assertEqual(ordered.get_linear_solver_stats().nb_fallback_factorize, 0)
        self.assertEqual(ordered.get_linear_solver_stats().nb_factorize, 1)
        np.testing.assert_array_equal(default.converged_mask(), ordered.converged_mask())
        np.testing.assert_allclose(ordered.get_voltages(), default.get_voltages(), rtol=0., atol=1e-10)
        self.assertEqual([r.nb_outer_iterations for r in ordered.get_outer_results()],
                         [r.nb_outer_iterations for r in default.get_outer_results()])

    def test_static_pivot_opt_out(self):
        self.assertEqual(self._run(self._grid("NROuter_KLUStaticPivot", False))
                         .get_linear_solver_stats().nb_fallback_factorize, 1)

    def test_kept_by_copies_and_the_resume_flag_untouched(self):
        grid = tob._resume_on(self._grid("NROuter_KLU", True))
        for other in (grid.copy(), copy.deepcopy(grid)):
            self.assertEqual(list(other.get_ac_algo_config().int_params[ORDERING:]), [1, 1])

    def test_a_solver_that_cannot(self):
        from lightsim2grid.algorithm import NR_SparseLU
        solver = NR_SparseLU()
        with self.assertRaises(RuntimeError):
            solver.constrained_ordering = True
        # a config that asks for it (another algorithm's) is ignored
        cfg = NR_KLUStaticPivot().get_config()
        solver.set_config(cfg)
        self.assertFalse(solver.constrained_ordering)
        self.assertEqual(solver.get_config().int_params[ORDERING], 0)

    def test_a_config_without_it_keeps_the_current_choice(self):
        solver = NR_KLU()
        solver.constrained_ordering = True
        cfg = solver.get_config()
        cfg.int_params = list(cfg.int_params)[:ORDERING]
        solver.set_config(cfg)
        self.assertTrue(solver.constrained_ordering)

    def test_standalone_switch(self):
        solver = NR_KLU()
        solver.constrained_ordering = True
        self.assertTrue(solver.constrained_ordering)
        self.assertEqual(solver.get_config().int_params[ORDERING], 1)
        solver.constrained_ordering = False
        self.assertEqual(solver.get_config().int_params[ORDERING], 0)


if __name__ == "__main__":
    unittest.main()
