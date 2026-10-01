# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``LSGrid.cap_slack_at_active_limits``: OpenLoadFlow's distributed-slack rule on a solved
grid. A unit of the distributed slack the last ``ac_pf`` pushed past an active limit leaves the
slack at that limit (flagged "can participate"), and the grid is solved again until no unit is
pushed out.

The reference is the grid built by hand the way the routine leaves it: the capped units at
their limit, outside the slack.

The grid: a radial 138 kV feeder 0-1-2-3 with an 80 MW load on bus 3 and generators on buses 0,
1 and 2 sharing the slack (generator 0 at 0 MW: the slack has ~70 MW to share).
"""

import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid

MAX_IT, TOL = 30, 1e-11
GEN, STO = 5, 6
HIGH_P, LOW_P = 7, 8


def _grid(targets=(0., 10., 10.), slack=(1., 1., 1.), pmax=(np.nan, np.nan, np.nan), storage=None):
    """``storage``: (target_p_mw load convention, slack weight, min_p, max_p) of a battery on bus 2."""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
    grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                         np.array([0, 1, 2]), np.array([1, 2, 3]))
    grid.init_loads(np.array([80.]), np.array([20.]), np.array([3]))
    n = len(targets)
    grid.init_generators(np.array(targets, dtype=float), np.array([1.02] + [1.01] * (n - 1)),
                         np.full(n, -1e3), np.full(n, 1e3), np.arange(n))
    for k, w in enumerate(slack):
        if w > 0.:
            grid.add_gen_slackbus(k, w)
    grid.set_gen_p_limits(np.full(n, np.nan), np.array(pmax, dtype=float))
    if storage is not None:
        t, w, lo, hi = storage
        grid.init_storages_full(np.array([t]), np.array([0.]), [False], np.array([1.]),
                                np.array([-1e3]), np.array([1e3]), np.array([2], dtype=np.int32))
        grid.add_storage_slackbus(0, w)
        grid.set_storage_p_limits(np.array([lo]), np.array([hi]))
    grid.tell_solver_need_reset()
    V = grid.ac_pf(np.full(4, 1.0 + 0j), MAX_IT, TOL)
    assert V.shape[0] > 0
    return grid, V


def _slack_p(grid):
    return [(int(v.element_type), int(v.element_id), int(v.violation_type))
            for v in grid.get_physical_violations(True, 1e-6, 0.)
            if int(v.violation_type) in (HIGH_P, LOW_P)]


class TestCapSlackAtActiveLimits(unittest.TestCase):
    def test_caps_and_resolves(self):
        grid, _ = _grid(pmax=(np.nan, 20., np.nan))
        self.assertEqual(_slack_p(grid), [(GEN, 1, HIGH_P)])
        out = grid.cap_slack_at_active_limits(MAX_IT, TOL)
        self.assertEqual([(int(v.element_type), int(v.element_id), int(v.violation_type)) for v in out],
                         [(GEN, 1, HIGH_P)])
        self.assertGreater(out[0].value, 20.)
        self.assertEqual(out[0].limit, 20.)
        self.assertEqual(_slack_p(grid), [])
        g1 = grid.get_generators()[1]
        self.assertFalse(g1.is_slack)
        self.assertTrue(g1.can_participate_slack)
        self.assertEqual(g1.can_participate_slack_weight, 1.)
        self.assertEqual(g1.target_p_mw, 20.)
        self.assertAlmostEqual(g1.res_p_mw, 20., places=9)
        # the grid built that way by hand
        ref, V_ref = _grid(targets=(0., 20., 10.), slack=(1., 0., 1.), pmax=(np.nan, 20., np.nan))
        np.testing.assert_allclose(np.abs(grid.get_V()), np.abs(V_ref), atol=1e-10, rtol=0)
        np.testing.assert_allclose([g.res_p_mw for g in grid.get_generators()],
                                   [g.res_p_mw for g in ref.get_generators()], atol=1e-8, rtol=0)

    def test_second_round(self):
        """capping generator 1 hands its share to the others, which pushes generator 2 past
        a max_p it was below at first"""
        grid0, _ = _grid(pmax=(np.nan, 20., np.nan))
        p2 = grid0.get_generators()[2].res_p_mw
        grid, _ = _grid(pmax=(np.nan, 20., p2 + 1.))
        self.assertEqual(_slack_p(grid), [(GEN, 1, HIGH_P)])
        out = grid.cap_slack_at_active_limits(MAX_IT, TOL)
        self.assertEqual([int(v.element_id) for v in out], [1, 2])
        self.assertEqual(_slack_p(grid), [])
        self.assertAlmostEqual(grid.get_generators()[2].res_p_mw, p2 + 1., places=9)
        self.assertEqual([g.is_slack for g in grid.get_generators()], [True, False, False])

    def test_nothing_to_cap(self):
        grid, V = _grid(pmax=(np.nan, 1e3, 1e3))
        self.assertEqual(grid.cap_slack_at_active_limits(MAX_IT, TOL), [])
        np.testing.assert_array_equal(grid.get_V(), V)
        self.assertEqual([g.is_slack for g in grid.get_generators()], [True, True, True])

    def test_every_unit_saturated_keeps_them(self):
        """OpenLoadFlow keeps every unit in the slack when all of them saturate"""
        grid, _ = _grid(pmax=(10., 20., 20.))
        before = _slack_p(grid)
        self.assertEqual(len(before), 3)
        self.assertEqual(grid.cap_slack_at_active_limits(MAX_IT, TOL), [])
        self.assertEqual(_slack_p(grid), before)

    def test_storage_unit(self):
        """a battery in the slack pushed past its max_p (generator convention): its set-point,
        in the container's load convention, moves to -max_p"""
        grid, _ = _grid(targets=(0., 10.), slack=(1., 1.), pmax=(np.nan, np.nan),
                        storage=(-5., 1., -1e3, 15.))
        self.assertEqual(_slack_p(grid), [(STO, 0, HIGH_P)])
        out = grid.cap_slack_at_active_limits(MAX_IT, TOL)
        self.assertEqual([(int(v.element_type), int(v.element_id)) for v in out], [(STO, 0)])
        sto = grid.get_storages()[0]
        self.assertFalse(sto.is_slack)
        self.assertTrue(sto.can_participate_slack)
        self.assertEqual(sto.target_p_mw, -15.)
        self.assertEqual(_slack_p(grid), [])

    def test_needs_a_solved_grid(self):
        grid, _ = _grid()
        fresh = LSGrid()
        with self.assertRaises(RuntimeError):
            fresh.cap_slack_at_active_limits()
        # and works on a copy, leaving the original as it was
        grid2, _ = _grid(pmax=(np.nan, 20., np.nan))
        cp = grid2.copy()
        cp.ac_pf(grid2.get_V(), MAX_IT, TOL)
        self.assertEqual(len(cp.cap_slack_at_active_limits(MAX_IT, TOL)), 1)
        self.assertTrue(grid2.get_generators()[1].is_slack)
        self.assertEqual(_slack_p(grid2), [(GEN, 1, HIGH_P)])


if __name__ == "__main__":
    unittest.main()
