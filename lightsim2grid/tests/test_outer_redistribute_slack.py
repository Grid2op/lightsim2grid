# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`redistribute_slack` with outer loops: the DistributedSlack loop shares what a contingency
loses before the row's first solve, as OpenLoadFlow's security analysis does (its Default
ContingencyActivePowerLossDistribution), and starts from those targets.

Without it the first solve puts the whole loss on the slack bus, and the loop then shares the
loss and the change in the losses together: where a unit reaches a bound, the order matters.
Here the units at 0 MW cannot take a share of a lost load, but they do take one of the losses
the loop shares afterwards."""

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


def _grid():
    """case14, its only outer loop DistributedSlack, every generator sharing the slack with the
    same weight, within [0, 300] MW"""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    nb_gen = len(grid.get_generators())
    p_max = np.full(nb_gen, 300.)
    grid.set_gen_p_limits(np.zeros(nb_gen), p_max)
    grid.set_gen_can_participate_slack([True] * nb_gen, p_max)
    grid.clear_outer_loops()
    grid.add_outer_loop(DistributedSlack())
    grid.change_algorithm(OUTER)
    return grid


def _v0(grid):
    return np.full(grid.total_bus(), 1.04 + 0j)


def _leaf_lines(grid):
    return [el.id for el in grid.get_lines() if LEAF_BUS in (el.bus1_id, el.bus2_id)]


def _batch(redistribute_slack):
    """the contingency islanding the leaf bus, solved by the batch"""
    grid = _grid()
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
    share of the lost load OpenLoadFlow gives them when `pre_share`"""
    grid = _grid()
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
        # a copy: get_gen_res views the grid's own buffers, and `ref` dies with this call
        p = np.array(ref.get_gen_res()[0])
        np.testing.assert_allclose(ca.get_gen_results()[0, :, 0], p, rtol=0., atol=1e-6)
        return p

    def test_loss_shared_before_the_first_solve(self):
        p = self._check(_batch(redistribute_slack=True), pre_share=True)
        # the units at 0 MW took a share of the losses only
        self.assertTrue(np.all(p[p < 10.] > 0.))

    def test_without_the_option_the_first_solve_takes_it_all(self):
        p = self._check(_batch(redistribute_slack=False), pre_share=False)
        # the loop shared the lost load and the losses together: what is left goes down, the
        # units at 0 MW stay there
        self.assertTrue(np.all(p[p < 10.] == 0.))


if __name__ == "__main__":
    unittest.main()
