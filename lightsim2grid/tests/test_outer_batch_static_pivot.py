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
KLU's), and the batch still pays no analysis and no factorization beyond the plain Newton's."""

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
    pass


if __name__ == "__main__":
    unittest.main()
