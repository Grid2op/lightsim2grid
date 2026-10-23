# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A batch row whose contingency strands a generator taking part in the slack.

With handle_disconnected_grid, the island is masked: it is not solved, so its unit takes no
share of the slack -- the row is the single solve of the grid without the island, whatever
weight the stranded unit carries. Without it, the row is not solved at all, and its
OuterRowResult says so (`converged`), although its loops never ran into anything (`status`)."""

import os
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.algorithm import DistributedSlack, OuterLoopStatus

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER = "NROuter_KLU"
except ImportError:
    OUTER = "NROuter_SparseLU"

MAX_IT = 30
TOL = 1e-8
# case14's generator on bus 8 (bus 7 here) hangs on one transformer only
LEAF_BUS = 7


def _grid(leaf_weight=300.):
    """case14, its only outer loop DistributedSlack, every generator sharing the slack in
    proportion to its maximum output; the one on the leaf bus has `leaf_weight`"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    gen_bus = [el.bus_id for el in grid.get_generators()]
    p_max = np.full(len(gen_bus), 300.)
    p_max[gen_bus.index(LEAF_BUS)] = leaf_weight
    grid.set_gen_p_limits(np.zeros(len(gen_bus)), p_max)
    grid.set_gen_can_participate_slack([True] * len(gen_bus), p_max)
    grid.clear_outer_loops()
    grid.add_outer_loop(DistributedSlack())
    grid.change_algorithm(OUTER)
    return grid


def _v0(grid):
    return np.full(grid.total_bus(), 1.04 + 0j)


def _leaf_trafo(grid):
    return [el.id for el in grid.get_trafos() if LEAF_BUS in (el.bus1_id, el.bus2_id)][0]


def _leaf_gen(grid):
    return [el.bus_id for el in grid.get_generators()].index(LEAF_BUS)


def _batch(grid, handle_disconnected_grid):
    """the contingency opening the leaf's transformer, and an ordinary line outage: the batch,
    the row of the first and the row of the second (the batch orders its contingencies)"""
    grid.ac_pf(_v0(grid), MAX_IT, TOL)
    ca = ContingencyAnalysisCPP(grid)
    ca.handle_disconnected_grid = handle_disconnected_grid
    leaf_branch = len(grid.get_lines()) + _leaf_trafo(grid)
    ca.add_nk([leaf_branch])
    ca.add_n1(0)
    ca.change_algorithm(OUTER)
    ca.compute_gen_results = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ca.compute(_v0(grid), MAX_IT, TOL)
    rows = [list(el) for el in ca.my_defaults()]
    return ca, rows.index([leaf_branch]), rows.index([0])


def _reference(leaf_weight):
    """the single solve of the grid without the island"""
    grid = _grid(leaf_weight)
    grid.deactivate_trafo(_leaf_trafo(grid))
    grid.consider_only_main_component()
    V = grid.ac_pf(_v0(grid), MAX_IT, TOL)
    return grid, V


class TestOuterIslandedSlackUnit(unittest.TestCase):
    def _check(self, leaf_weight):
        ca, row, _ = _batch(_grid(leaf_weight), handle_disconnected_grid=True)
        ref, V_ref = _reference(leaf_weight)
        self.assertGreater(V_ref.size, 0)
        self.assertTrue(ca.converged_mask()[row])
        res = ca.get_outer_results()[row]
        self.assertEqual(res.status, OuterLoopStatus.STABLE)
        self.assertEqual(res.nb_outer_iterations, ref.get_algo().get_outer_loop_stats().nb_outer_iterations)
        # the stranded unit's target is left alone
        leaf = _leaf_gen(_grid(leaf_weight))
        self.assertNotIn(leaf, [gen_id for gen_id, _ in res.gen_target_p])
        # the rest is the single solve's
        V = ca.get_voltages()[row]
        np.testing.assert_array_equal(np.flatnonzero(np.abs(V) == 0.), [LEAF_BUS])  # the island
        main = np.abs(V) > 0.
        np.testing.assert_allclose(V[main], V_ref[main], rtol=0., atol=1e-10)
        p, _, _ = ref.get_gen_res()
        np.testing.assert_allclose(ca.get_gen_results()[row, :, 0], p, rtol=0., atol=1e-6)

    def test_row_is_the_single_solve_without_the_island(self):
        self._check(leaf_weight=300.)

    def test_stranded_unit_with_most_of_the_weight(self):
        # sharing with it, the slack left on the slack bus shrinks by its weight's fraction per
        # outer iteration only: the row ran out of iterations
        self._check(leaf_weight=30000.)


class TestOuterRowResultConverged(unittest.TestCase):
    def test_converged_is_the_row_outcome(self):
        # the island is not handled: that row cannot be solved, its loops never ran
        ca, islanded, line = _batch(_grid(), handle_disconnected_grid=False)
        converged = np.asarray(ca.converged_mask(), dtype=bool)
        self.assertFalse(converged[islanded])
        self.assertTrue(converged[line])
        res = ca.get_outer_results()
        self.assertEqual([r.converged for r in res], list(converged))
        # the loops' own status says nothing of the row's outcome
        self.assertEqual(res[islanded].status, OuterLoopStatus.STABLE)


if __name__ == "__main__":
    unittest.main()
