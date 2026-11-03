# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""test_outer_batch's equivalences again, with NROuter_KLUStaticPivot: KLU's factorization, then
refactorizations on its frozen pivot order. Each row of a batch still equals the single solve of
its grid (to the tolerance those tests use: the static-pivot factors are not bit-identical to
KLU's), and the batch still pays no analysis and no factorization beyond the plain Newton's.
And the rows the outer loops switch by value need no fallback factorization: the solver orders
them after the Vm columns that keep their pivot usable in both forms."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import test_outer_batch as tob

try:
    from lightsim2grid.algorithm import NROuter_KLUStaticPivot  # noqa: F401
    STATIC_PIVOT = "NROuter_KLUStaticPivot"
except ImportError:
    STATIC_PIVOT = None


class _StaticPivot:
    """runs the test_outer_batch case with the static-pivot algorithm in place of its own"""

    def setUp(self):
        if STATIC_PIVOT is None:
            self.skipTest("KLU is not available")
        saved = tob.OUTER
        tob.OUTER = STATIC_PIVOT
        self.addCleanup(setattr, tob, "OUTER", saved)
        super().setUp()


class TestStaticPivotContingencyAnalysis(_StaticPivot, tob.TestOuterContingencyAnalysis):
    pass


class TestStaticPivotScenarioSweep(_StaticPivot, tob.TestOuterScenarioSweep):
    pass


class TestStaticPivotInjectionSweep(_StaticPivot, tob.TestOuterInjectionSweep):
    pass


class TestStaticPivotTimeSeries(_StaticPivot, tob.TestOuterTimeSeries):
    pass


class TestStaticPivotContingencyAnalysisFromN(_StaticPivot, tob.TestOuterContingencyAnalysisFromN):
    pass


class TestStaticPivotInjectionSweepFromN(_StaticPivot, tob.TestOuterInjectionSweepFromN):
    pass


class TestStaticPivotTapsInContingencyAnalysis(_StaticPivot, tob.TestOuterTapsInContingencyAnalysis):
    def test_the_switched_rows_need_no_fallback(self):
        # the rows that move the tap switch the ratio group's row from "Vm = target" to the
        # ratio held: on AMD's order alone, its frozen pivot sits on the entry that form zeroes
        batch = self._batch()
        batch.change_algorithm(STATIC_PIVOT)
        batch.compute(tob._v0(self.base), tob.MAX_IT, tob.TOL)
        self.assertTrue(any(r.ratio_tap for r in batch.get_outer_results()))
        stats = batch.get_linear_solver_stats()
        self.assertEqual(stats.nb_refactorize_failed, 0)
        self.assertEqual(stats.nb_fallback_factorize, 0)
        self.assertEqual(stats.nb_factorize, 1)


if __name__ == "__main__":
    unittest.main()
