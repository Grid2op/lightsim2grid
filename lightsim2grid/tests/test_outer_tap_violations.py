# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A transformer whose taps the outer loops moved: its current limit is checked at the taps
the solve settled on, the current its flows report -- in a batch row and in a single solve."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import pypowsybl  # noqa: F401
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

OUTER = "NROuter_KLU"
TOL = 1e-10


def _pst_grid():
    """IEEE 14, its phase shifter regulating 7 MW above its flow without control: the solve
    rounds it to tap 3; its side-1 current limit at half its current, always violated"""
    import test_phase_control as P
    from lightsim2grid.algorithm import PhaseControl
    free = P._grid(P._net("ACTIVE_POWER_CONTROL", 0., regulating=False))
    free.change_algorithm(OUTER)
    free.ac_pf(np.full(free.total_bus(), 1. + 0j), 30, TOL)
    pst = next(i for i, t in enumerate(free.get_trafos()) if t.name == P.PST)
    grid = P._grid(P._net("ACTIVE_POWER_CONTROL", free.get_trafos()[pst].res_p1_mw + 7.))
    grid.add_outer_loop(PhaseControl())
    grid.change_algorithm(OUTER)
    lim = np.full(len(grid.get_trafos()), np.nan)
    lim[pst] = 0.5 * free.get_trafos()[pst].res_a1_ka
    grid.set_trafo_current_limit_side1(lim)
    return grid, pst


def _current(violations, pst):
    return [v.value for v in violations if v.element_type.name == "TRAFO" and v.element_id == pst
            and v.violation_type.name == "CURRENT" and v.side == 1]


class TestOuterTapViolations(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def test_batch_row(self):
        from lightsim2grid.lightsim2grid_cpp import ScenarioSweepCPP
        grid, pst = _pst_grid()
        grid.ac_pf(np.full(grid.total_bus(), 1. + 0j), 30, TOL)
        batch = ScenarioSweepCPP(grid)
        batch.change_algorithm(OUTER)
        batch.compute_limit_violations = True
        batch.modify_load_p(np.asarray(grid.get_load_target_p())[None, :])
        batch.compute(np.asarray(grid.get_V(), dtype=complex), 30, TOL)
        self.assertTrue(batch.get_outer_results()[0].phase_tap)  # the tap moved
        batch.compute_flows()
        n_line = len(grid.get_lines())
        flow_ka = np.asarray(batch.get_flows())[0, n_line + pst]
        got = _current(batch.get_violations()[0], pst)
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got[0], flow_ka, places=10)

    def test_single_solve(self):
        grid, pst = _pst_grid()
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1. + 0j), 30, TOL).size, 0)
        t = grid.get_trafos()[pst]
        self.assertNotEqual(t.res_phase_tap_position, t.phase_tap_position)
        got = _current(grid.get_violations(), pst)
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got[0], t.res_a1_ka, places=10)


if __name__ == "__main__":
    unittest.main()
