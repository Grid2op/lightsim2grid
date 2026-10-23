# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Detection mode is the loop's own trigger: after a plain solve, what get_physical_violations
reports for a loop is exactly what that loop decides in its first check of an outer-loop solve
started from the same state -- the same elements, the same values, the same limits."""

import os
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import ReactiveLimits
from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP, LimitViolationType,
                                             ViolationElementType)

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER, SING = "NROuter_KLU", "NRSing_KLU"
except ImportError:
    OUTER, SING = "NROuter_SparseLU", "NRSing_SparseLU"

Q_TYPES = (LimitViolationType.LOW_Q, LimitViolationType.HIGH_Q)


def _grid(max_q_mvar):
    """case14, its generators' reactive range narrowed, ReactiveLimits its only loop"""
    net = pn.case14()
    net.gen["min_q_mvar"] = -max_q_mvar
    net.gen["max_q_mvar"] = max_q_mvar
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(net)
    grid.clear_outer_loops()
    grid.add_outer_loop(ReactiveLimits())
    return grid


def _v0(grid):
    return np.full(grid.total_bus(), 1.04 + 0j)


def _detected(grid):
    """(bus, type, value, limit) of what detection reports on the last solve"""
    return sorted((v.element_id, v.violation_type, v.value, v.limit)
                  for v in grid.get_physical_violations(True, 0., 0.)
                  if v.element_type == ViolationElementType.BUS and v.violation_type in Q_TYPES)


def _first_check(grid):
    """(bus, type, value, limit) of the buses ReactiveLimits switches PV -> PQ in its first
    check -- the strongest one it keeps PV included: the trigger fired there too"""
    stats = grid.get_algo().get_outer_loop_stats()
    return sorted((d.element_id, d.reason, d.value, d.limit) for d in stats.decisions
                  if d.loop == "ReactiveLimits" and d.outer_iteration == 0
                  and d.action in ("PV_TO_PQ", "KEPT_PV_STRONGEST"))


class TestDetectionIsTheTrigger(unittest.TestCase):
    def _compare(self, max_q_mvar):
        plain = _grid(max_q_mvar)
        plain.change_algorithm(SING)
        self.assertGreater(plain.ac_pf(_v0(plain), 30, 1e-10).size, 0)
        outer = _grid(max_q_mvar)
        outer.change_algorithm(OUTER)
        outer.ac_pf(_v0(outer), 30, 1e-10)
        detected = _detected(plain)
        triggered = _first_check(outer)
        self.assertEqual([d[:2] for d in detected], [t[:2] for t in triggered])
        for d, t in zip(detected, triggered):
            self.assertAlmostEqual(d[2], t[2], delta=1e-6)
            self.assertAlmostEqual(d[3], t[3], delta=1e-9)
        return detected

    def test_tight_limits(self):
        self.assertGreater(len(self._compare(15.)), 0)

    def test_wide_limits_trigger_nothing(self):
        self.assertEqual(self._compare(1000.), [])

    def test_a_batch_row_detects_on_its_own_grid(self):
        # a contingency row reports what a plain solve of the row's grid reports
        base = _grid(15.)
        base.change_algorithm(SING)
        ca = ContingencyAnalysisCPP(base)
        ca.change_algorithm(SING)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        for line_id in range(6):
            ca.add_n1(line_id)
        ca.compute(_v0(base), 30, 1e-10)
        rows = ca.get_physical_violations()
        converged = np.asarray(ca.converged_mask(), dtype=bool)
        for line_id in range(6):
            grid = _grid(15.)
            grid.change_algorithm(SING)
            grid.deactivate_powerline(line_id)
            V = grid.ac_pf(_v0(grid), 30, 1e-10)
            with self.subTest(line=line_id):
                self.assertEqual(converged[line_id], V.size > 0)
                if not V.size:
                    continue
                ref = _detected(grid)
                got = sorted((v.element_id, v.violation_type, v.value, v.limit) for v in rows[line_id]
                             if v.element_type == ViolationElementType.BUS and v.violation_type in Q_TYPES)
                self.assertEqual([r[:2] for r in got], [r[:2] for r in ref])
                for g, r in zip(got, ref):
                    self.assertAlmostEqual(g[2], r[2], delta=1e-6)
                    self.assertAlmostEqual(g[3], r[3], delta=1e-9)


if __name__ == "__main__":
    unittest.main()
