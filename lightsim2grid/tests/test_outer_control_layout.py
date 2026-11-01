# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""get_outer_control_layout: where the continuous controls of the outer loops (the phase
shifters, the ratio groups, the shunt controllers) sit in the Newton's system -- what another
solver of the same equations reads to pose them."""

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


def _solved(grid):
    grid.change_algorithm(OUTER)
    V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
    assert V.size > 0
    return grid


class TestOuterControlLayout(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _check_owned(self, grid, rows, cols):
        # every row / column it reports is one of J's, each owned once
        dim = grid.get_J_solver().shape[0]
        for idx in rows + cols:
            self.assertTrue(0 <= idx < dim)
        self.assertEqual(len(set(rows)), len(rows))
        self.assertEqual(len(set(cols)), len(cols))

    @staticmethod
    def _bus_columns(grid):
        """each solver bus' angle and magnitude columns (-1: none). The bus block of J is
        square: a bus' P row has the index of its angle column, its Q row that of its Vm one"""
        algo = grid.get_algo()
        return np.asarray(algo.get_theta_to_J_col()), np.asarray(algo.get_vm_to_J_col())

    def _check_own_column(self, grid, col):
        # an unknown of its own: no bus' angle nor magnitude
        theta, vm = self._bus_columns(grid)
        self.assertNotIn(col, set(theta.tolist()) | set(vm.tolist()))

    def test_phase_shifter(self):
        from test_phase_control import _net, _grid, PST
        grid = _solved(_grid(_net("ACTIVE_POWER_CONTROL", 50.)))
        from lightsim2grid.algorithm import PhaseControl  # noqa: F401
        grid.clear_outer_loops()
        grid.add_outer_loop(PhaseControl())
        grid = _solved(grid)
        layout = grid.get_algo().get_outer_control_layout()
        pst = next(i for i, t in enumerate(grid.get_trafos()) if t.name == PST)
        self.assertEqual([b.trafo for b in layout.branches], [pst])
        b = layout.branches[0]
        self.assertTrue(b.phase and b.shift_column)
        self.assertFalse(b.ratio or b.ratio_column)
        self.assertEqual(b.ratio_col, -1)
        tr = grid.get_trafos()[pst]
        to_solver = np.asarray(grid.id_ac_solver_to_me())
        self.assertEqual((to_solver[b.bus1_solver], to_solver[b.bus2_solver]), (tr.bus1_id, tr.bus2_id))
        self._check_owned(grid, [b.shift_row], [b.shift_col])
        # the shift column is the shifter's: the flow it carries moves with it, out of one of
        # its buses and into the other (dP/dshift on their P rows, of opposite signs)
        self._check_own_column(grid, b.shift_col)
        theta, _ = self._bus_columns(grid)
        J = grid.get_J_solver().tocsc()
        self.assertGreaterEqual(min(theta[b.bus1_solver], theta[b.bus2_solver]), 0)
        self.assertLess(J[theta[b.bus1_solver], b.shift_col] * J[theta[b.bus2_solver], b.shift_col], 0.)
        self.assertEqual(len(layout.ratio_groups), 0)
        self.assertEqual(len(layout.shunts), 0)

    def test_ratio_group(self):
        from test_outer_batch import _tap_grid
        grid = _solved(_tap_grid())
        layout = grid.get_algo().get_outer_control_layout()
        trafo = next(i for i, t in enumerate(grid.get_trafos()) if t.name == "T4-9-1")
        self.assertEqual(len(layout.ratio_groups), 1)
        g = layout.ratio_groups[0]
        self.assertEqual([layout.branches[m].trafo for m in g.members], [trafo])
        b = layout.branches[g.members[0]]
        self.assertTrue(b.ratio and b.ratio_column)
        self.assertEqual(len(g.rows), len(g.members))
        self.assertAlmostEqual(g.target_vm, 12.6 / grid.get_bus_vn_kv()[np.asarray(grid.id_ac_solver_to_me())[g.bus_solver]])
        self._check_owned(grid, list(g.rows), [b.ratio_col, g.vm_col])
        # the ratio column is the transformer's, it moves the power at both of its buses; the
        # group's Vm column is the regulated bus' own
        self._check_own_column(grid, b.ratio_col)
        theta, vm = self._bus_columns(grid)
        self.assertEqual(g.vm_col, vm[g.bus_solver])
        J = grid.get_J_solver().tocsc()
        for bus in (b.bus1_solver, b.bus2_solver):
            self.assertGreaterEqual(theta[bus], 0)
            self.assertNotEqual(J[theta[bus], b.ratio_col], 0.)

    def test_shunt_group(self):
        from test_shunt_voltage_control import _net, _grid
        from lightsim2grid.algorithm import ShuntVoltageControl
        grid = _grid(_net("B9", 12.6))
        grid.clear_outer_loops()
        grid.add_outer_loop(ShuntVoltageControl())
        grid = _solved(grid)
        layout = grid.get_algo().get_outer_control_layout()
        self.assertEqual(len(layout.shunts), 1)  # the regulating shunts of a bus: one controller
        sh = layout.shunts[0]
        self.assertEqual(len(sh.shunts), 2)
        self.assertTrue(sh.column)
        self.assertEqual(len(layout.shunt_groups), 1)
        g = layout.shunt_groups[0]
        self.assertEqual(list(g.members), [0])
        self._check_owned(grid, list(g.rows) + [sh.q_row], [sh.col, g.vm_col])
        # the susceptance column is the controller's, read by its bus' Q row (dQ/dB); the
        # group's Vm column is the regulated bus' own
        self._check_own_column(grid, sh.col)
        _, vm = self._bus_columns(grid)
        self.assertEqual(sh.q_row, vm[sh.bus_solver])
        self.assertEqual(g.vm_col, vm[g.bus_solver])
        self.assertNotEqual(grid.get_J_solver().tocsc()[sh.q_row, sh.col], 0.)
        self.assertEqual(len(layout.branches), 0)

    def test_nothing_without_such_loops(self):
        from test_outer_batch import _grid
        grid = _grid()
        grid.ac_pf(np.full(grid.total_bus(), 1.04 + 0j), 30, 1e-8)
        layout = grid.get_algo().get_outer_control_layout()
        self.assertEqual((len(layout.branches), len(layout.ratio_groups), len(layout.shunts), len(layout.shunt_groups)),
                         (0, 0, 0, 0))


if __name__ == "__main__":
    unittest.main()
