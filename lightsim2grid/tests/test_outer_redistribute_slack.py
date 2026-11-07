# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`redistribute_slack` with outer loops: the DistributedSlack loop shares what a contingency
loses before the row's first solve, as OpenLoadFlow's security analysis does (its Default
ContingencyActivePowerLossDistribution), so that the first solve does not put the whole loss on
the slack bus.

What the loop shares afterwards (the change in the losses) is shared, as OpenLoadFlow does, as a
cumulative mismatch from the targets the row started from, not from those the first share left:
the result is that of the whole imbalance shared at once. It only matters where a unit reaches a
bound and the two have opposite signs. Here the units at 0 MW cannot take a share of a lost load,
and the losses the loop shares afterwards do not bring them back up."""

import os
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
from lightsim2grid.algorithm import DistributedSlack

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER = "NROuter_KLU"
except ImportError:
    OUTER = "NROuter_SparseLU"

MAX_IT = 30
TOL = 1e-8
# case14's bus 14 (13 here): a load fed by two lines only
LEAF_BUS = 13


def _grid(slack_mismatch_mw=1.):
    """case14, its only outer loop DistributedSlack (stopping below `slack_mismatch_mw`, and
    sharing down to it when it is below the loop's default residue), every
    generator sharing the slack with the same weight, within [0, 300] MW"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    nb_gen = len(grid.get_generators())
    p_max = np.full(nb_gen, 300.)
    grid.set_gen_p_limits(np.zeros(nb_gen), p_max)
    grid.set_gen_can_participate_slack([True] * nb_gen, p_max)
    grid.clear_outer_loops()
    grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=slack_mismatch_mw,
                                         p_residue_eps_mw=min(slack_mismatch_mw, 1e-3)))
    grid.change_algorithm(OUTER)
    return grid


# the loop shares down to a tiny mismatch: what it leaves on the slack bus below its threshold
# would otherwise depend on how the imbalance was split between its steps
TIGHT_MW = 1e-6


def _v0(grid):
    return np.full(grid.total_bus(), 1.04 + 0j)


def _leaf_lines(grid):
    return [el.id for el in grid.get_lines() if LEAF_BUS in (el.bus1_id, el.bus2_id)]


def _batch(redistribute_slack):
    """the contingency islanding the leaf bus, solved by the batch"""
    grid = _grid(TIGHT_MW)
    grid.ac_pf(_v0(grid), MAX_IT, TOL)
    ca = ContingencyAnalysisCPP(grid)
    ca.handle_disconnected_grid = True
    ca.redistribute_slack = redistribute_slack
    ca.add_nk(_leaf_lines(grid))
    ca.change_algorithm(OUTER)
    ca.compute_gen_results = True
    ca.compute(_v0(grid), MAX_IT, TOL)
    return ca


def _reference(pre_share):
    """the single solve of the grid without the leaf, the units' targets moved first by the
    share of the lost load when `pre_share` (what the loop would give if it shared what is left
    from those targets)"""
    grid = _grid(TIGHT_MW)
    lost_mw = -sum(el.target_p_mw for el in grid.get_loads() if el.bus_id == LEAF_BUS)
    if pre_share:
        # every unit has the same weight, and the ones at 0 MW cannot go lower: the two
        # others take half of it each
        movable = [el for el in grid.get_generators() if el.target_p_mw > 0.]
        assert len(movable) == 2
        for el in movable:
            grid.change_p_gen(el.id, el.target_p_mw + 0.5 * lost_mw)
            # change_p_gen puts a unit that can take part in the slack in the Newton's slack
            # set, which the batch's single-slack Newton does not have: out again, and back
            # among the DistributedSlack loop's units
            if not el.is_slack:
                grid.remove_gen_slackbus(el.id)
        nb_gen = len(grid.get_generators())
        grid.set_gen_can_participate_slack([True] * nb_gen, np.full(nb_gen, 300.))
    for line_id in _leaf_lines(grid):
        grid.deactivate_powerline(line_id)
    grid.consider_only_main_component(False)
    V = grid.ac_pf(_v0(grid), MAX_IT, TOL)
    return grid, V


class TestOuterRedistributeSlack(unittest.TestCase):
    def _check(self, ca, pre_share):
        ref, V_ref = _reference(pre_share)
        self.assertGreater(V_ref.size, 0)
        self.assertTrue(ca.converged_mask()[0])
        V = ca.get_voltages()[0]
        np.testing.assert_array_equal(np.flatnonzero(np.abs(V) == 0.), [LEAF_BUS])  # the island
        main = np.abs(V) > 0.
        np.testing.assert_allclose(V[main], V_ref[main], rtol=0., atol=1e-10)
        # get_gen_res is a view into `ref`, which does not keep it alive: `p` outlives it
        p = ref.get_gen_res()[0].copy()
        np.testing.assert_allclose(ca.get_gen_results()[0, :, 0], p, rtol=0., atol=1e-6)
        return p

    def test_loss_shared_before_the_first_solve(self):
        ca = _batch(redistribute_slack=True)
        self.assertIn(("DistributedSlack/SHARE_LOST_P", 1), ca.get_outer_results()[0].decisions_taken)
        # the same state as the whole imbalance shared at once (pre_share=False): the units at
        # 0 MW, that the lost load sent to their bound, stay there (the three condensers: with
        # none of them below 10 MW, the check below would check nothing)
        p = self._check(ca, pre_share=False)
        self.assertEqual(int((p < 10.).sum()), 3)
        self.assertTrue(np.all(p[p < 10.] == 0.))

    def test_without_the_option_the_first_solve_takes_it_all(self):
        p = self._check(_batch(redistribute_slack=False), pre_share=False)
        # the loop shared the lost load and the losses together: what is left goes down, the
        # units at 0 MW stay there
        self.assertEqual(int((p < 10.).sum()), 3)
        self.assertTrue(np.all(p[p < 10.] == 0.))


class TestOuterRedistributeSlackFromN(unittest.TestCase):
    """With init_from_n_powerflow, a row counts what it loses at the targets in force after the
    "n" solve, as OpenLoadFlow does: a generator whose grid target is 0 MW, but that the "n"
    solve's DistributedSlack moved up, takes that much generation with it when tripped -- even
    though the loss at the grid's targets is nothing."""

    TRIPPED = 3  # a synchronous condenser, 0 MW in the grid

    def test_a_unit_at_zero_in_the_grid_moved_by_the_n_solve(self):
        grid = _grid()
        self.assertEqual(grid.get_generators()[self.TRIPPED].target_p_mw, 0.)
        grid.ac_pf(_v0(grid), MAX_IT, TOL)
        # the non-slack units produce their target exactly: the "n" solve's moved targets
        p_n = grid.get_gen_res()[0].copy()
        is_slack = np.array([el.is_slack for el in grid.get_generators()])
        self.assertGreater(p_n[self.TRIPPED], 1.)

        ca = ContingencyAnalysisCPP(grid)
        ca.init_from_n_powerflow = True
        ca.redistribute_slack = True
        ca.add_contingency_gen(self.TRIPPED)
        ca.change_algorithm(OUTER)
        ca.compute(_v0(grid), MAX_IT, TOL)
        res = ca.get_outer_results()[0]
        self.assertTrue(res.converged)
        moved = dict(res.gen_target_p)
        # every remaining unit has the same weight and stays within its range: each takes a
        # quarter of what the tripped one produced, before the first solve; what the losses
        # change afterwards is below the loop's 1 MW threshold, so nothing moves after it
        for gen_id in range(len(p_n)):
            if gen_id == self.TRIPPED or is_slack[gen_id]:
                continue
            self.assertAlmostEqual(moved[gen_id], p_n[gen_id] + 0.25 * p_n[self.TRIPPED], places=6,
                                   msg=f"generator {gen_id}")



class TestOuterRedistributeSlackAtABound(unittest.TestCase):
    """A lost generator: the share before the first solve sends a unit to its max_p, and the
    losses then drop, so the loop takes back part of what it shared. Shared as a cumulative
    mismatch from the targets the row started from (the "n" solve's, with
    init_from_n_powerflow), the total is still an increase: the unit stays at its bound, as with
    the whole imbalance shared at once. Shared from the targets the first share left, it would
    come back down by its share of the correction."""

    TRIPPED = 0  # the generator at bus 2 (index 1)
    CAPPED = 2   # a synchronous condenser, moved up by the "n" solve's distributed slack
    HEADROOM_MW = 2.  # what the capped unit can still take: less than its share of the loss

    def _grid(self, max_p_capped):
        # the loop shares what the losses change too (below 1 MW here), the weights stay _grid's
        grid = _grid(TIGHT_MW)
        nb_gen = len(grid.get_generators())
        p_max = np.full(nb_gen, 300.)
        p_max[self.CAPPED] = max_p_capped
        grid.set_gen_p_limits(np.zeros(nb_gen), p_max)
        return grid

    def test_a_unit_at_its_bound_stays_there(self):
        grid = self._grid(300.)
        grid.ac_pf(_v0(grid), MAX_IT, TOL)
        cap = grid.get_gen_res()[0][self.CAPPED] + self.HEADROOM_MW
        grid = self._grid(cap)
        grid.ac_pf(_v0(grid), MAX_IT, TOL)
        p_n = grid.get_gen_res()[0].copy()
        self.assertLess(p_n[self.CAPPED], cap - 1.)  # not at its bound in the "n" state
        # every remaining unit has the same weight: the capped one cannot take its quarter
        self.assertGreater(0.25 * p_n[self.TRIPPED], self.HEADROOM_MW)

        ca = ContingencyAnalysisCPP(grid)
        ca.init_from_n_powerflow = True
        ca.redistribute_slack = True
        ca.add_contingency_gen(self.TRIPPED)
        ca.change_algorithm(OUTER)
        ca.compute(_v0(grid), MAX_IT, TOL)
        res = ca.get_outer_results()[0]
        self.assertTrue(res.converged)
        moved = dict(res.gen_target_p)
        shared = sum(moved[g] - p_n[g] for g in moved if g != self.TRIPPED)
        # the losses dropped: the loop shared less than what was lost, taking some back
        self.assertLess(shared, p_n[self.TRIPPED] - 0.1)
        self.assertGreater(shared, 0.)
        self.assertAlmostEqual(moved[self.CAPPED], cap, places=6)


if __name__ == "__main__":
    unittest.main()
