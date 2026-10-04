# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The "can participate in the slack" flag (``LSGrid.set_gen_can_participate_slack``): a unit an
outer loop left out of the distributed slack only because it sat at an active limit takes part in
the bounded redistribution pre-pass -- so it moves away from that limit, never across it -- and
never in the Newton solve's slack; ``bake_outer_loops(..., return_details=True)`` /
``init_from_pypowsybl(can_participate_slack=...)`` flag the units OpenLoadFlow capped."""

import copy
import os
import pickle
import tempfile
import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl import bake_outer_loops, BakeResult
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

VN_KV = 138.
CAPPED = 1          # gen 1: at its max_p, out of the slack
LEAF_LINE = 3       # line 1-4: taking it out islands bus 4


def _grid(flagged=True, leaf_gen_mw=0., overshoot_mw=0.):
    """buses 0-1-2-3 in a row plus a leaf bus 4 off bus 1 (line 3); 60 MW of load on bus 3 and
    20 MW on bus 4 (and a `leaf_gen_mw` generator there). Gen 0, on bus 0, is the slack; gen 1,
    on bus 2, is dispatched at its max_p and out of the slack -- flagged (or not) "can
    participate in the slack", with the same weight as gen 0."""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(5, 1, np.full(5, VN_KV), 0, 0)
    grid.init_powerlines(np.full(4, 0.01), np.full(4, 0.1), np.zeros(4, dtype=complex),
                         np.array([0, 1, 2, 1]), np.array([1, 2, 3, 4]))
    grid.init_loads(np.array([60., 20.]), np.array([10., 5.]), np.array([3, 4]))
    p = [30., 40.] + ([leaf_gen_mw] if leaf_gen_mw else [])
    bus = [0, 2] + ([4] if leaf_gen_mw else [])
    n = len(p)
    grid.init_generators_full(np.array(p), np.full(n, 1.02), np.zeros(n), [True] * n,
                              np.full(n, -1e3), np.full(n, 1e3), np.array(bus))
    grid.set_gen_p_limits(np.zeros(n), np.array([500., 40.] + ([100.] if leaf_gen_mw else [])))
    grid.add_gen_slackbus(0, 0.5)
    if flagged:
        grid.set_gen_can_participate_slack(np.array([False, True] + [False] * (n - 2)),
                                           np.array([0., 0.5] + [0.] * (n - 2)))
        if overshoot_mw:
            grid.set_gen_can_participate_slack_overshoot(np.array([0., overshoot_mw] + [0.] * (n - 2)))
    grid.tell_solver_need_reset()
    return grid


def _target_p(grid):
    return np.array([g.target_p_mw for g in grid.get_generators()])


class TestFlag(unittest.TestCase):

    def test_default_set_and_refused(self):
        grid = _grid(flagged=False)
        gen = grid.get_generators()[CAPPED]
        self.assertFalse(gen.can_participate_slack)
        self.assertEqual(gen.can_participate_slack_weight, 0.)
        grid.set_gen_can_participate_slack(np.array([False, True]), np.array([0., 0.5]))
        gen = grid.get_generators()[CAPPED]
        self.assertTrue(gen.can_participate_slack)
        self.assertEqual(gen.can_participate_slack_weight, 0.5)
        self.assertFalse(gen.is_slack)   # the Newton solve's slack is unchanged
        with self.assertRaises(RuntimeError):   # wrong size
            grid.set_gen_can_participate_slack(np.array([True]), np.array([1.]))
        with self.assertRaises(RuntimeError):   # flagged without a weight
            grid.set_gen_can_participate_slack(np.array([False, True]), np.array([0., 0.]))

    def test_kept_by_copy_pickle_and_binary(self):
        grid = _grid(overshoot_mw=15.)
        for other in (grid.copy(), copy.deepcopy(grid), pickle.loads(pickle.dumps(grid))):
            self.assertEqual(other.get_generators()[CAPPED].can_participate_slack_weight, 0.5)
            self.assertEqual(other.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 15.)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            gen = LSGrid.load_binary(path).get_generators()[CAPPED]
            self.assertEqual(gen.can_participate_slack_weight, 0.5)
            self.assertEqual(gen.can_participate_slack_overshoot_mw, 15.)

    def test_overshoot_default_and_refused(self):
        grid = _grid()
        self.assertEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 0.)
        for bad in ([0., np.nan], [0., np.inf], [0.]):
            with self.assertRaises(RuntimeError):
                grid.set_gen_can_participate_slack_overshoot(np.array(bad))
        # signed: < 0 for a unit capped at its lower limit
        grid.set_gen_can_participate_slack_overshoot(np.array([0., -1.]))
        self.assertEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, -1.)


class TestPrepass(unittest.TestCase):
    """consider_only_main_component(redistribute_slack=True) on an island of load (the
    generation must go down) and of generation (it must go up)."""

    def test_flagged_unit_takes_a_share_away_from_its_limit(self):
        grid = _grid(flagged=True)
        grid.deactivate_powerline(LEAF_LINE)   # 20 MW of load leave: the units inject less
        report = grid.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, -20., places=9)
        self.assertEqual(report.nb_participants, 2)
        # equal weights: 10 MW each, the capped unit moving down from its max_p
        np.testing.assert_allclose(_target_p(grid), [20., 30.], atol=1e-9)
        gens = grid.get_generators()
        self.assertTrue(gens[0].is_slack)
        self.assertFalse(gens[CAPPED].is_slack)   # still out of the Newton solve's slack

    def test_overshoot_used_up_first(self):
        # OLF's p = clamp(v + shift * w / W): the capped unit sits at 40 + 15 MW unbounded, so
        # it only leaves its max_p once the common shift took those 15 MW; 20 MW to take with
        # equal weights: shift -35 MW, the slack unit gives 17.5, the capped one 2.5
        grid = _grid(flagged=True, overshoot_mw=15.)
        grid.deactivate_powerline(LEAF_LINE)
        report = grid.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, -20., places=9)
        self.assertAlmostEqual(report.not_distributed_mw, 0., places=6)
        np.testing.assert_allclose(_target_p(grid), [12.5, 37.5], atol=1e-6)

    def test_overshoot_larger_than_the_shift_keeps_it_capped(self):
        # 25 MW beyond its max_p: the 20 MW the units must give all come from the slack one
        grid = _grid(flagged=True, overshoot_mw=25.)
        grid.deactivate_powerline(LEAF_LINE)
        grid.consider_only_main_component(True)
        np.testing.assert_allclose(_target_p(grid), [10., 40.], atol=1e-6)

    def test_overshoot_used_up_across_redistributions(self):
        # 25 MW beyond its max_p. A first -20 MW: shift -40, the slack unit gives 20, the capped
        # one stays, 5 MW of its overshoot left. A second -10 MW starts from those 5 MW, as
        # OLF sharing -30 MW at once would: shift -15, the slack unit gives 7.5, the capped 2.5
        grid = _grid(flagged=True, overshoot_mw=25.)
        grid.redistribute_active_power(-20.)
        np.testing.assert_allclose(_target_p(grid), [10., 40.], atol=1e-6)
        self.assertAlmostEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 5., places=6)
        grid.redistribute_active_power(-10.)
        np.testing.assert_allclose(_target_p(grid), [2.5, 37.5], atol=1e-6)
        self.assertEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 0.)
        # the same as one -30 MW redistribution from the start
        once = _grid(flagged=True, overshoot_mw=25.)
        once.redistribute_active_power(-30.)
        np.testing.assert_allclose(_target_p(once), _target_p(grid), atol=1e-6)

    def test_overshoot_does_not_outlive_the_cap(self):
        # moved off its limit, the unit is back in the slack: the reference solve's overshoot
        # is gone. Saturated again by a redistribution, it carries no overshoot (a
        # redistribution never makes one, as OLF's distribution step), so it leaves its max_p
        # at once on a mismatch of the other sign
        grid = _grid(flagged=True, overshoot_mw=25.)
        grid.change_p_gen(CAPPED, 30.)
        self.assertTrue(grid.get_generators()[CAPPED].is_slack)
        self.assertEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 0.)
        # +30 MW, equal weights: the unit at 30 is capped at 40, the slack one takes the rest
        grid.redistribute_active_power(30.)
        np.testing.assert_allclose(_target_p(grid), [50., 40.], atol=1e-6)
        self.assertFalse(grid.get_generators()[CAPPED].is_slack)
        self.assertEqual(grid.get_generators()[CAPPED].can_participate_slack_overshoot_mw, 0.)
        # -10 MW: 5 MW each, the old 25 MW overshoot would have kept it at 40
        grid.redistribute_active_power(-10.)
        np.testing.assert_allclose(_target_p(grid), [45., 35.], atol=1e-6)

    def test_overshoot_batch_matches_the_single_solve(self):
        ref = _grid(flagged=True, overshoot_mw=15.)
        ref.deactivate_powerline(LEAF_LINE)
        ref.consider_only_main_component(True)
        V_ref = ref.ac_pf(np.full(ref.total_bus(), 1.0 + 0j), 30, 1e-11)
        self.assertGreater(V_ref.shape[0], 0)
        grid = _grid(flagged=True, overshoot_mw=15.)
        ca = ContingencyAnalysisCPP(grid, True)
        ca.handle_disconnected_grid = True
        ca.redistribute_slack = True
        ca.add_n1(LEAF_LINE)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert list(ca.converged()) == [True]
        V_batch = np.asarray(ca.get_voltages())[0]
        np.testing.assert_allclose(V_batch[:4], V_ref[:4], atol=1e-8)

    def test_not_flagged_takes_nothing(self):
        grid = _grid(flagged=False)
        grid.deactivate_powerline(LEAF_LINE)
        report = grid.consider_only_main_component(True)
        self.assertEqual(report.nb_participants, 1)
        np.testing.assert_allclose(_target_p(grid), [10., 40.], atol=1e-9)

    def test_flagged_unit_never_crosses_its_limit(self):
        # the island produces more than it consumes: the units must inject more, and the
        # unit at its max_p cannot -- the slack one takes it all
        grid = _grid(flagged=True, leaf_gen_mw=30.)
        grid.deactivate_powerline(LEAF_LINE)
        report = grid.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, 10., places=9)
        tp = _target_p(grid)
        np.testing.assert_allclose(tp[:2], [40., 40.], atol=1e-9)
        self.assertTrue(grid.get_generators()[0].is_slack)   # did not saturate

    def _grid_only_flagged_has_room(self):
        # the island produces more than it consumes (+10 MW): the units must inject more.
        # The slack unit can only take 2 MW (max_p 12), the flagged one has room to spare.
        # 10 + 40 + 30 MW of generation for 80 MW of load: balanced but for the losses
        grid = _grid(flagged=True, leaf_gen_mw=30.)
        grid.change_p_gen(0, 10.)
        grid.set_gen_p_limits(np.zeros(3), np.array([12., 100., 100.]))
        return grid

    def test_slack_kept_when_only_a_flagged_unit_has_room(self):
        # the flagged unit takes the rest, but it is not in the slack the powerflow
        # distributes on: the saturated slack unit must stay in it, or there is none left
        grid = self._grid_only_flagged_has_room()
        grid.deactivate_powerline(LEAF_LINE)
        report = grid.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, 10., places=9)
        self.assertTrue(report.all_saturated)
        np.testing.assert_allclose(_target_p(grid)[:2], [12., 48.], atol=1e-9)
        gens = grid.get_generators()
        self.assertTrue(gens[0].is_slack)
        self.assertFalse(gens[CAPPED].is_slack)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        self.assertGreater(V.shape[0], 0)

    def test_batch_keeps_the_slack_when_only_a_flagged_unit_has_room(self):
        # same as above, row by row: the slack unit keeps its share of the losses, and the
        # p-limit check sees it past its max_p, as on the single solve
        ref = self._grid_only_flagged_has_room()
        ref.deactivate_powerline(LEAF_LINE)
        ref.consider_only_main_component(True)
        V_ref = ref.ac_pf(np.full(ref.total_bus(), 1.0 + 0j), 30, 1e-11)
        self.assertGreater(V_ref.shape[0], 0)
        viol_ref = {(v.element_type, v.element_id, v.violation_type) for v in ref.get_physical_violations()}

        grid = self._grid_only_flagged_has_room()
        ca = ContingencyAnalysisCPP(grid, True)
        ca.handle_disconnected_grid = True
        ca.redistribute_slack = True
        ca.compute_physical_violations = True
        ca.add_n1(LEAF_LINE)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert list(ca.converged()) == [True]
        V_batch = np.asarray(ca.get_voltages())[0]
        np.testing.assert_allclose(V_batch[:4], V_ref[:4], atol=1e-9)
        viol_batch = {(v.element_type, v.element_id, v.violation_type) for v in ca.get_physical_violations()[0]}
        self.assertEqual(viol_batch, viol_ref)
        self.assertTrue(any(el_id == 0 for _, el_id, _ in viol_batch), "the slack unit ends past its max_p")

    def test_the_batch_matches_the_single_solve(self):
        ref = _grid(flagged=True)
        ref.deactivate_powerline(LEAF_LINE)
        ref.consider_only_main_component(True)
        V_ref = ref.ac_pf(np.full(ref.total_bus(), 1.0 + 0j), 30, 1e-11)
        self.assertGreater(V_ref.shape[0], 0)

        grid = _grid(flagged=True)
        ca = ContingencyAnalysisCPP(grid, True)
        ca.handle_disconnected_grid = True
        ca.redistribute_slack = True
        ca.add_n1(LEAF_LINE)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert list(ca.converged()) == [True]
        V_batch = np.asarray(ca.get_voltages())[0]
        # the islanded bus 4 is out of the single solve, NaN / masked in the batch
        np.testing.assert_allclose(V_batch[:4], V_ref[:4], atol=1e-9)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestFromPypowsybl(unittest.TestCase):
    """GTH1 of the four substations network, dispatched at its max_p, is capped by OLF's
    distribution of the (positive) losses: the bake leaves it out of the slack, and returns it."""

    def _baked(self):
        n = pp.network.create_four_substations_node_breaker_network()
        lf.run_ac(n)
        return n, bake_outer_loops(n, return_details=True)

    def test_bake_returns_the_capped_unit(self):
        n, res = self._baked()
        self.assertIsInstance(res, BakeResult)
        self.assertEqual(list(res.can_participate_slack), ["GTH1"])
        self.assertFalse(n.get_extensions("activePowerControl").loc["GTH1", "participate"])
        # the plain return is unchanged
        n2 = pp.network.create_four_substations_node_breaker_network()
        lf.run_ac(n2)
        self.assertEqual(list(bake_outer_loops(n2)), list(res.can_be_pv))

    def test_bake_returns_the_overshoot_and_init_sets_it(self):
        n, res = self._baked()
        over = res.can_participate_slack_overshoot
        self.assertEqual(list(over.index), list(res.can_participate_slack))
        # capped at its max_p by a positive mismatch: beyond its upper limit, > 0
        self.assertTrue(np.all(np.isfinite(over.to_numpy())) and np.all(over.to_numpy() >= 0.))
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                   can_participate_slack=res.can_participate_slack,
                                   can_participate_slack_overshoot=over)
        gens = {g.name: g for g in grid.get_generators()}
        self.assertAlmostEqual(gens["GTH1"].can_participate_slack_overshoot_mw, over["GTH1"], places=9)

    def test_overshoot_without_the_flag_refused(self):
        # the overshoot only means something for a flagged unit: given alone (or with an
        # explicit slack) it is refused, as `can_participate_slack` itself is
        n, res = self._baked()
        over = res.can_participate_slack_overshoot
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                can_participate_slack_overshoot=over)
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, gen_slack_id="GTH2",
                                can_participate_slack_overshoot=over)

    def test_init_flags_it_with_olf_weight(self):
        n, res = self._baked()
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                   can_participate_slack=res.can_participate_slack)
        gens = {g.name: g for g in grid.get_generators()}
        self.assertFalse(gens["GTH1"].is_slack)
        self.assertTrue(gens["GTH1"].can_participate_slack)
        self.assertTrue(gens["GTH2"].is_slack)
        # max_p / droop, normalised with the slack weights: GTH1 100 MW, GTH2 400 MW
        self.assertAlmostEqual(gens["GTH1"].can_participate_slack_weight / gens["GTH2"].slack_weight,
                               100. / 400., places=12)
        # nothing flagged out of the slack by default (a slack unit carries the flag on its
        # own), an unknown id is refused, and an explicit slack too
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False)
        self.assertFalse(any(g.can_participate_slack and not g.is_slack for g in grid.get_generators()))
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                can_participate_slack=["NOT-A-UNIT"])
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, gen_slack_id="GTH2",
                                can_participate_slack=res.can_participate_slack)


if __name__ == "__main__":
    unittest.main()
