# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``LSGrid.set_hold_frozen_regulators``: a generator an outer loop froze at a reactive limit
(``can_be_pv``, regulation off) that would regulate a REMOTE bus keeps its seat in that bus'
voltage-control group, held at its frozen reactive output. The system solved must be the one
without the option -- same voltages, same reactive outputs, same physical checks, in a single
solve and in the batch algorithms -- while the Jacobian carries the held machine's reactive
unknown and its pinned row.

The grid: a 138 kV ring 0-1-2-3 with leaves 4, 5, 6 and 7 behind short lines.
- gen 0 on bus 0: the slack;
- gen 1 on leaf 4, frozen, would regulate bus 1: a group of held machines only;
- gen 2 on leaf 5, regulating bus 2 remotely (active), and gen 3 on leaf 6, frozen, would
  regulate bus 2 at the same set-point: a held machine joining an active group;
- gen 4 on leaf 7, frozen, would regulate bus 2 at ANOTHER set-point: left out.
"""

import pickle
import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid, ContingencyAnalysisCPP

LINES = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 4), (2, 5), (2, 6), (2, 7)]
STEP_UP = {5, 6, 7, 8}
# (bus, p, vm, q, regulating, min_q, max_q, regulated bus)
GENS = [(0, 0., 1.03, 0., True, -500., 500., 0),
        (4, 30., 1.02, 25., False, -20., 25., 1),
        (5, 20., 1.03, 0., True, -40., 40., 2),
        (6, 15., 1.03, -10., False, -10., 30., 2),
        (7, 10., 1.05, 5., False, -5., 5., 2)]
CAN_BE_PV = [False, True, False, True, True]
MAX_IT, TOL = 30, 1e-11


def _grid(hold=False):
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(8, 1, np.full(8, 138.), 0, 0)
    x = np.array([0.03 if k in STEP_UP else 0.08 for k in range(len(LINES))])
    grid.init_powerlines(np.full(len(LINES), 0.01), x, np.full(len(LINES), 0.02j),
                         np.array([a for a, _ in LINES]), np.array([b for _, b in LINES]))
    grid.init_loads(np.array([60., 50., 40.]), np.array([25., 20., 15.]), np.array([1, 2, 3]))
    grid.init_generators_full(np.array([g[1] for g in GENS]), np.array([g[2] for g in GENS]),
                              np.array([g[3] for g in GENS]), [g[4] for g in GENS],
                              np.array([g[5] for g in GENS]), np.array([g[6] for g in GENS]),
                              np.array([g[0] for g in GENS]))
    for k, g in enumerate(GENS):
        if g[7] != g[0]:
            grid.set_gen_regulated_bus(k, g[7])
    grid.set_gen_can_be_pv(np.array(CAN_BE_PV))
    grid.add_gen_slackbus(0, 1.)
    grid.set_hold_frozen_regulators(hold)
    grid.tell_solver_need_reset()
    return grid


def _solve(grid):
    V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), MAX_IT, TOL)
    assert V.shape[0] > 0
    return V


class TestHoldFrozenRegulators(unittest.TestCase):
    def test_default_and_copy(self):
        grid = _grid()
        self.assertFalse(grid.get_hold_frozen_regulators())
        grid.set_hold_frozen_regulators(True)
        self.assertTrue(grid.get_hold_frozen_regulators())
        self.assertTrue(grid.copy().get_hold_frozen_regulators())

    def test_held_controllers(self):
        """gens 1 and 3 are held, gen 4 (another set-point than its group's) is left out; the
        held one of the mixed group comes after the active one."""
        grid = _grid(hold=True)
        _solve(grid)
        elem = np.asarray(grid.get_controller_elem_id_solver())
        held = np.asarray(grid.get_controller_held_solver())
        self.assertEqual(sorted(elem[held == 1].tolist()), [1, 3])
        self.assertEqual(sorted(elem[held == 0].tolist()), [2])
        self.assertNotIn(4, elem.tolist())
        self.assertLess(elem.tolist().index(2), elem.tolist().index(3))
        q = np.asarray(grid.get_controller_q_solver()) * grid.get_sn_mva()
        for gen_id in (1, 3):
            self.assertAlmostEqual(q[elem.tolist().index(gen_id)], GENS[gen_id][3], places=8)
        # and without the option, no held controller at all
        grid_off = _grid()
        _solve(grid_off)
        self.assertEqual(np.asarray(grid_off.get_controller_held_solver()).tolist(), [0])

    def test_same_solution(self):
        grid_off, grid_on = _grid(), _grid(hold=True)
        V_off, V_on = _solve(grid_off), _solve(grid_on)
        np.testing.assert_allclose(V_on, V_off, atol=1e-10, rtol=0)
        q_off = np.asarray(grid_off.get_gen_res()[1])
        q_on = np.asarray(grid_on.get_gen_res()[1])
        np.testing.assert_allclose(q_on, q_off, atol=1e-8, rtol=0)
        # the held machines' unknowns and rows are in the Jacobian
        self.assertGreater(grid_on.get_J_solver().shape[0], grid_off.get_J_solver().shape[0])

    def test_dc_unchanged(self):
        grid_off, grid_on = _grid(), _grid(hold=True)
        v0 = np.full(grid_off.total_bus(), 1.0 + 0j)
        np.testing.assert_array_equal(grid_on.dc_pf(1. * v0, 10, 1e-8), grid_off.dc_pf(1. * v0, 10, 1e-8))

    def test_can_be_pv_change_rebuilds(self):
        """with the option on, can_be_pv decides who is held: changing it is seen by the next
        solve, and the solution still does not move"""
        grid = _grid(hold=True)
        V_a = _solve(grid)
        grid.set_gen_can_be_pv(np.array([False, False, False, True, True]))
        V_b = _solve(grid)
        elem = np.asarray(grid.get_controller_elem_id_solver())
        held = np.asarray(grid.get_controller_held_solver())
        self.assertEqual(elem[held == 1].tolist(), [3])
        np.testing.assert_allclose(V_b, V_a, atol=1e-10, rtol=0)

    def test_switching_the_option_rebuilds(self):
        grid = _grid()
        V_a = _solve(grid)
        grid.set_hold_frozen_regulators(True)
        V_b = _solve(grid)
        self.assertEqual(int(np.asarray(grid.get_controller_held_solver()).sum()), 2)
        grid.set_hold_frozen_regulators(False)
        V_c = _solve(grid)
        self.assertEqual(int(np.asarray(grid.get_controller_held_solver()).sum()), 0)
        np.testing.assert_allclose(V_b, V_a, atol=1e-10, rtol=0)
        np.testing.assert_allclose(V_c, V_a, atol=1e-10, rtol=0)

    def test_physical_violations_unchanged(self):
        """the reactive-capability check and the release check report the same either way: a
        held machine holds no bus"""
        def records(grid):
            _solve(grid)
            return sorted((str(v.violation_type), int(v.element_id), round(float(v.value), 6))
                          for v in grid.get_physical_violations(True, 0., 0.))
        rec_off, rec_on = records(_grid()), records(_grid(hold=True))
        self.assertEqual(rec_on, rec_off)
        self.assertTrue(any("AT_M" in r[0] for r in rec_on), "the scenario reports a release")

    def test_contingency_analysis_unchanged(self):
        """every contingency (islanding the leaves too, handle_disconnected_grid) gives the
        same voltages and physical violations with the option on"""
        res = {}
        for hold in (False, True):
            grid = _grid(hold=hold)
            _solve(grid)
            ca = ContingencyAnalysisCPP(grid)
            ca.handle_disconnected_grid = True
            ca.compute_physical_violations = True
            ca.physical_violation_tol_mva = 0.
            ca.physical_violation_tol_vm_pu = 0.
            for branch in range(len(LINES)):
                ca.add_n1(branch)
            ca.compute(np.full(grid.total_bus(), 1.0 + 0j), MAX_IT, TOL)
            res[hold] = (np.asarray(ca.get_voltages()),
                         [sorted((str(v.violation_type), int(v.element_id), round(float(v.value), 6))
                                 for v in row) for row in ca.get_physical_violations()])
        np.testing.assert_allclose(res[True][0], res[False][0], atol=1e-9, rtol=0, equal_nan=True)
        self.assertEqual(res[True][1], res[False][1])

    def test_pickle_does_not_carry_it(self):
        """like set_keep_vinit_at_group_controlled_buses: a solver option, not grid data"""
        grid = _grid(hold=True)
        self.assertFalse(pickle.loads(pickle.dumps(grid)).get_hold_frozen_regulators())


if __name__ == "__main__":
    unittest.main()
