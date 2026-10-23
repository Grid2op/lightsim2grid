# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A batch row that islands a continuous control (handle_disconnected_grid): the island is not
solved, so a shunt controller, a phase shifter or a ratio group inside it does nothing -- the
row is the one of a grid without that control. Each grid below has one control only, inside
the island the contingency cuts off."""

import os
import sys
import unittest
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import pypowsybl  # noqa: F401
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

OUTER = "NROuter_KLU"
TOL = 1e-10


def _branch_ids(grid, names):
    lines = [el.name for el in grid.get_lines()]
    trafos = [el.name for el in grid.get_trafos()]
    return [lines.index(n) if n in lines else len(lines) + trafos.index(n) for n in names]


def _row(make_grid, loop, branches):
    """the voltages of the contingency row, and whether it converged"""
    grid = make_grid()
    grid.clear_outer_loops()
    if loop is not None:
        grid.add_outer_loop(loop())
    grid.change_algorithm(OUTER)
    grid.ac_pf(np.full(grid.total_bus(), 1. + 0j), 30, TOL)
    ca = ContingencyAnalysisCPP(grid)
    ca.handle_disconnected_grid = True
    ca.add_nk(_branch_ids(grid, branches))
    ca.change_algorithm(OUTER)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ca.compute(np.full(grid.total_bus(), 1. + 0j), 30, TOL)
    return ca.get_voltages()[0], bool(ca.converged_mask()[0])


class TestOuterIslandedControls(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _check(self, make_grid, loop, branches):
        V_ref, ok_ref = _row(make_grid, None, branches)
        self.assertTrue(ok_ref)
        self.assertTrue(np.any(np.abs(V_ref) == 0.))  # an island, or this checks nothing
        V, ok = _row(make_grid, loop, branches)
        self.assertTrue(ok)
        np.testing.assert_allclose(V, V_ref, rtol=0., atol=1e-8)

    def test_shunt_controller_on_an_island(self):
        import test_shunt_voltage_control as S
        from lightsim2grid.algorithm import ShuntVoltageControl
        self._check(lambda: S._grid(S._net("B14", 12.8)), ShuntVoltageControl, ["L9-14-1", "L13-14-1"])

    def test_phase_shifter_on_an_island(self):
        import test_phase_control as P
        from lightsim2grid.algorithm import PhaseControl
        self._check(lambda: P._grid(P._net("ACTIVE_POWER_CONTROL", 30.)), PhaseControl,
                    ["L1-5-1", "L2-5-1", "L4-5-1", "L9-14-1", "L10-11-1"])

    def test_ratio_group_on_an_island(self):
        import test_transformer_voltage_control as T
        from lightsim2grid.algorithm import TransformerVoltageControl
        self._check(lambda: T._grid(T._net("T4-7-1", 12.6, 0.1)), TransformerVoltageControl,
                    ["L2-4-1", "L3-4-1", "L4-5-1", "T4-9-1", "L7-9-1"])


if __name__ == "__main__":
    unittest.main()
