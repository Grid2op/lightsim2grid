# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`ContingencyAnalysis(CPP).violations_new_or_worse_only`: get_violations() keeps the violations
a contingency adds to the base case (absent from it, or present with a worse value), and the
default is unchanged (the complete state of the contingency)."""

import unittest
import warnings

import numpy as np

from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP, ViolationElementType, LimitViolationType
from test_ContingencyAnalysis_limit_violations import _set_tight_limits

REL = 1e-6  # NEW_OR_WORSE_REL_TOL of BaseBatchSweep.hpp


def _sig(rows):
    """(key, value) per violation; NaN (the grid-level entries) compares equal"""
    return [(_key(v), repr(v.value)) for v in rows]


def _key(v):
    return (v.element_type, v.element_id, v.side, v.violation_type)


def _is_new_or_worse(v, n_values):
    if v.element_type == ViolationElementType.GRID:
        return True
    base = n_values.get(_key(v))
    if base is None:
        return True
    margin = REL * abs(base)
    if v.violation_type == LimitViolationType.LOW_VOLTAGE:
        return v.value < base - margin
    return v.value > base + margin


class TestNewOrWorseViolations(unittest.TestCase):
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        self.backend = self.env.backend
        # every element violates in the base case, so what a contingency changes is only
        # how far above its limit an element sits
        _set_tight_limits(self.backend._grid)

    def tearDown(self):
        self.env.close()

    def _run(self, new_or_worse=None):
        ca = ContingencyAnalysisCPP(self.backend._grid, True)
        if new_or_worse is not None:
            ca.violations_new_or_worse_only = new_or_worse
        ca.add_all_n1()
        ca.compute(self.backend.V, self.backend.max_it, self.backend.tol)
        return ca

    def test_default_is_off_and_reports_the_complete_state(self):
        ca = self._run()
        self.assertFalse(ca.violations_new_or_worse_only)
        explicit = self._run(False)
        self.assertEqual([_sig(row) for row in ca.get_violations()],
                         [_sig(row) for row in explicit.get_violations()])
        n_keys = {_key(v) for v in ca.get_violations_n()}
        # tight limits: the complete state repeats base case violations
        self.assertTrue(any(_key(v) in n_keys for row in ca.get_violations() for v in row))

    def test_only_new_or_worse_are_kept(self):
        full = self._run(False)
        reduced = self._run(True)
        self.assertTrue(reduced.violations_new_or_worse_only)
        n_values = {_key(v): v.value for v in full.get_violations_n()}
        # the base case itself is reported in full
        self.assertEqual(_sig(reduced.get_violations_n()),
                         _sig(full.get_violations_n()))
        n_rows = len(full.get_violations())
        self.assertEqual(len(reduced.get_violations()), n_rows)
        nb_worse, nb_dropped = 0, 0
        for row_full, row_red in zip(full.get_violations(), reduced.get_violations()):
            expected = _sig([v for v in row_full if _is_new_or_worse(v, n_values)])
            self.assertEqual(_sig(row_red), expected)
            nb_dropped += len(row_full) - len(row_red)
            nb_worse += sum(1 for v in row_red if _key(v) in n_values)
        self.assertGreater(nb_dropped, 0)  # unchanged pre-existing violations are gone
        self.assertGreater(nb_worse, 0)  # ... and those that got worse are not: not an "absent from n" filter

    def test_divergence_is_kept(self):
        ca = self._run(True)
        full = self._run(False)
        for row_full, row_red in zip(full.get_violations(), ca.get_violations()):
            sentinels_full = [v for v in row_full if v.element_type == ViolationElementType.GRID]
            sentinels_red = [v for v in row_red if v.element_type == ViolationElementType.GRID]
            self.assertEqual(len(sentinels_full), len(sentinels_red))

    def test_toggle_discards_results_and_recomputes(self):
        ca = self._run(False)
        nb_full = sum(len(row) for row in ca.get_violations())
        ca.violations_new_or_worse_only = True
        ca.compute(self.backend.V, self.backend.max_it, self.backend.tol)
        self.assertLess(sum(len(row) for row in ca.get_violations()), nb_full)
        ca.violations_new_or_worse_only = False
        ca.compute(self.backend.V, self.backend.max_it, self.backend.tol)
        self.assertEqual(sum(len(row) for row in ca.get_violations()), nb_full)

    def test_python_wrapper(self):
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        self.assertFalse(sa.violations_new_or_worse_only)
        sa.add_single_contingency(4)
        full = sa.run()
        sa.violations_new_or_worse_only = True
        self.assertTrue(sa.computer.violations_new_or_worse_only)
        reduced = sa.run()
        self.assertEqual(len(reduced.pre_contingency_result.limit_violations),
                         len(full.pre_contingency_result.limit_violations))
        self.assertLess(len(reduced.post_contingency_results[0].limit_violations),
                        len(full.post_contingency_results[0].limit_violations))


if __name__ == "__main__":
    unittest.main()
