# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Distributed slack across the boundary of ``consider_only_main_component``.

When islanding strands a generator of the distributed slack, ``consider_only_main_component``
deactivates it. It used to leave it in the slack: its (now disconnected) bus was still
listed as a slack bus and the next ``ac_pf`` threw "One of the slack bus is disconnected".
It now removes such generators from the slack first, as OpenLoadFlow, which only
distributes the slack on the main component.

With ``redistribute_slack=True`` (the default) it also shares the power the islanding took
out on the remaining slack units as OpenLoadFlow's ``DistributedSlack`` outer loop does,
with their ``[min_p, max_p]`` bounds and without crossing 0 MW
(``redistribute_active_power``): the tests below check it against a plain Python
re-implementation of that loop.
"""

import unittest
import warnings
import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower

# bus 7 of case14 is a leaf (degree 1) with a generator (synchronous condenser) on it
_LEAF_BUS = 7
_LEAF_P_MW = 40.


def olf_distribute(injection, weight, min_p, max_p, mismatch, eps=1e-6):
    """OpenLoadFlow's GenerationActivePowerDistributionStep, in plain Python.
    Returns (new injections, saturated flags, remaining)."""
    new = np.array(injection, dtype=float)
    active = np.ones(new.shape[0], dtype=bool)
    sat = np.zeros(new.shape[0], dtype=bool)
    lo = np.where(np.isfinite(min_p), min_p, -np.inf)
    hi = np.where(np.isfinite(max_p), max_p, np.inf)
    # "we don't want to change the generation sign": 0 is a bound on the other side
    lo = np.where(new < 0., lo, np.maximum(lo, 0.))
    hi = np.where(new < 0., np.minimum(hi, 0.), hi)
    remaining = float(mismatch)
    while active.any() and abs(remaining) > eps:
        factor_sum = weight[active].sum()
        done = 0.
        for k in np.flatnonzero(active):
            cand = new[k] + remaining * weight[k] / factor_sum
            if remaining > 0. and cand >= hi[k]:
                cand = max(hi[k], new[k])
                active[k] = False
                sat[k] = True
            elif remaining < 0. and cand <= lo[k]:
                cand = min(lo[k], new[k])
                active[k] = False
                sat[k] = True
            done += cand - new[k]
            new[k] = cand
        remaining -= done
    if not active.any():
        sat[:] = False
    return new, sat, remaining


class TestMainComponentSlack(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.net = pn.case14()
            self.model = init_from_pandapower(self.net)
        # distributed slack on every generator
        for gen in self.model.get_generators():
            self.model.add_gen_slackbus(gen.id, 1.)
        self.leaf_gen = [gen.id for gen in self.model.get_generators() if gen.bus_id == _LEAF_BUS]
        assert len(self.leaf_gen) == 1
        self.n_gen = len(self.model.get_generators())

    def _isolate_leaf(self):
        for line in self.model.get_lines():
            if _LEAF_BUS in (line.bus1_id, line.bus2_id):
                self.model.deactivate_powerline(line.id)
        for trafo in self.model.get_trafos():
            if _LEAF_BUS in (trafo.bus1_id, trafo.bus2_id):
                self.model.deactivate_trafo(trafo.id)

    def _ac_pf(self):
        return self.model.ac_pf(np.ones(self.net.bus.shape[0], dtype=np.complex128), 30, 1e-10)

    def _targets(self):
        return np.array([gen.target_p_mw for gen in self.model.get_generators()])

    def _slack_flags(self):
        return np.array([gen.is_slack for gen in self.model.get_generators()])

    def _give_leaf_power(self):
        """the leaf generator is a synchronous condenser (P = 0): give it some power, so
        that islanding it loses something"""
        self.model.change_p_gen(self.leaf_gen[0], _LEAF_P_MW)

    def test_stranded_gen_leaves_slack(self):
        gens = self.model.get_generators()
        assert gens[self.leaf_gen[0]].is_slack
        nb_slack = sum(gen.is_slack for gen in gens)
        self._isolate_leaf()
        self.model.consider_only_main_component()
        gens = self.model.get_generators()
        self.assertFalse(gens[self.leaf_gen[0]].connected)
        self.assertFalse(gens[self.leaf_gen[0]].is_slack, "the stranded generator must leave the slack")
        self.assertEqual(sum(gen.is_slack for gen in gens), nb_slack - 1)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0, "ac_pf diverged after islanding a slack generator")

    def test_stranded_forced_reference_is_cleared(self):
        self.model.set_reference_slack_bus(_LEAF_BUS)
        self._isolate_leaf()
        self.model.consider_only_main_component()
        self.assertEqual(self.model.get_reference_slack_bus(), -1)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)

    def test_nothing_stranded_unchanged(self):
        self.model.set_reference_slack_bus(_LEAF_BUS)
        nb_slack = sum(gen.is_slack for gen in self.model.get_generators())
        targets = self._targets()
        report = self.model.consider_only_main_component()
        self.assertEqual(sum(gen.is_slack for gen in self.model.get_generators()), nb_slack)
        self.assertEqual(self.model.get_reference_slack_bus(), _LEAF_BUS)
        # nothing lost: nothing redistributed
        self.assertEqual(report.mismatch_mw, 0.)
        self.assertEqual(report.nb_rounds, 0)
        np.testing.assert_array_equal(self._targets(), targets)

    # ---- redistribution of the lost power ----------------------------------------
    def test_redistribute_sum(self):
        self._give_leaf_power()
        targets = self._targets()
        self._isolate_leaf()
        report = self.model.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, _LEAF_P_MW, places=9)
        self.assertEqual(report.nb_participants, self.n_gen - 1)
        self.assertEqual(report.nb_saturated, 0)
        self.assertEqual(report.nb_rounds, 1)
        self.assertFalse(report.all_saturated)
        new = self._targets()
        others = np.arange(self.n_gen) != self.leaf_gen[0]
        # the lost 40 MW are on the other generators, by weight (all equal here)
        self.assertAlmostEqual(new[others].sum(), targets[others].sum() + _LEAF_P_MW, places=9)
        np.testing.assert_allclose(new[others] - targets[others], _LEAF_P_MW / (self.n_gen - 1))
        # the stranded generator is out: its target is not touched
        self.assertEqual(new[self.leaf_gen[0]], targets[self.leaf_gen[0]])
        # ... and nobody left the slack
        self.assertEqual(self._slack_flags()[others].sum(), self.n_gen - 1)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)

    def test_redistribute_false_is_old_behaviour(self):
        self._give_leaf_power()
        targets = self._targets()
        self._isolate_leaf()
        report = self.model.consider_only_main_component(False)
        # the loss is still reported, nothing is moved
        self.assertAlmostEqual(report.mismatch_mw, _LEAF_P_MW, places=9)
        self.assertEqual(report.nb_rounds, 0)
        np.testing.assert_array_equal(self._targets(), targets)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)

    def test_saturation(self):
        self._give_leaf_power()
        targets = self._targets()
        others = np.flatnonzero(np.arange(self.n_gen) != self.leaf_gen[0])
        # two of the remaining units can only take 3 MW (instead of their 10 MW share)
        max_p = np.full(self.n_gen, np.inf)
        min_p = np.full(self.n_gen, -np.inf)
        clamped = others[:2]
        max_p[clamped] = targets[clamped] + 3.
        self.model.set_gen_p_limits(min_p, max_p)
        self._isolate_leaf()
        report = self.model.consider_only_main_component(True)
        self.assertEqual(report.nb_saturated, 2)
        self.assertFalse(report.all_saturated)
        self.assertGreaterEqual(report.nb_rounds, 2)
        self.assertAlmostEqual(report.not_distributed_mw, 0., places=9)
        new = self._targets()
        flags = self._slack_flags()
        for g in clamped:
            self.assertAlmostEqual(new[g], max_p[g], places=9)
            self.assertFalse(flags[g], "a saturated unit leaves the distributed slack")
        free = [g for g in others if g not in clamped]
        for g in free:
            self.assertAlmostEqual(new[g] - targets[g], (_LEAF_P_MW - 6.) / len(free), places=9)
            self.assertTrue(flags[g])
        self.assertAlmostEqual(new[others].sum(), targets[others].sum() + _LEAF_P_MW, places=9)
        # the plain Python OLF loop says the same
        exp, sat, remaining = olf_distribute(targets[others], np.ones(others.size),
                                             min_p[others], max_p[others], _LEAF_P_MW)
        np.testing.assert_allclose(new[others], exp, atol=1e-9)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)
        # a saturated unit produces its clamped target and takes no share of the losses
        res_p = np.array([gen.res_p_mw for gen in self.model.get_generators()])
        for g in clamped:
            self.assertAlmostEqual(res_p[g], max_p[g], places=6)

    def test_large_limits_equal_no_limits(self):
        self._give_leaf_power()
        self._isolate_leaf()
        # (on a copy: limits far away)
        with_limits = self.model.copy()
        with_limits.set_gen_p_limits(np.full(self.n_gen, -1e4), np.full(self.n_gen, 1e4))
        rep_a = self.model.consider_only_main_component(True)
        rep_b = with_limits.consider_only_main_component(True)
        self.assertEqual(rep_b.nb_saturated, 0)
        self.assertEqual(rep_a.nb_saturated, 0)
        np.testing.assert_allclose(np.array([g.target_p_mw for g in with_limits.get_generators()]),
                                   self._targets(), atol=1e-12)

    def test_all_saturated_keeps_slack(self):
        self._give_leaf_power()
        targets = self._targets()
        others = np.flatnonzero(np.arange(self.n_gen) != self.leaf_gen[0])
        max_p = targets + 1.  # nobody can take more than 1 MW
        self.model.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        self._isolate_leaf()
        report = self.model.consider_only_main_component(True)
        self.assertTrue(report.all_saturated)
        self.assertEqual(report.nb_saturated, others.size)
        self.assertAlmostEqual(report.not_distributed_mw, _LEAF_P_MW - others.size, places=9)
        new = self._targets()
        np.testing.assert_allclose(new[others], max_p[others], atol=1e-12)
        # every unit is at its bound: they all stay in the slack, the powerflow still has one
        self.assertEqual(self._slack_flags()[others].sum(), others.size)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)

    def test_redistribute_active_power_standalone(self):
        # no islanding at all: share 30 MW then -25 MW on the slack units, with limits
        targets = self._targets()
        rng = np.random.default_rng(0)
        weights = rng.uniform(0.5, 2., self.n_gen)
        for g in range(self.n_gen):
            self.model.remove_gen_slackbus(g)
            self.model.add_gen_slackbus(g, float(weights[g]))
        max_p = targets + rng.uniform(2., 12., self.n_gen)
        min_p = targets - rng.uniform(2., 12., self.n_gen)
        max_p[1] = np.nan  # unbounded above
        # 30 MW over weights in [0.5, 2]: a share is at least 2 MW, so these two saturate
        max_p[0] = targets[0] + 1.
        max_p[2] = targets[2] + 1.5
        self.model.set_gen_p_limits(min_p, max_p)

        report = self.model.redistribute_active_power(30.)
        exp, sat, remaining = olf_distribute(targets, weights, min_p, max_p, 30.)
        self.assertAlmostEqual(report.mismatch_mw, 30.)
        self.assertEqual(report.nb_participants, self.n_gen)
        self.assertEqual(report.nb_saturated, int(sat.sum()))
        self.assertAlmostEqual(report.not_distributed_mw, remaining, places=9)
        np.testing.assert_allclose(self._targets(), exp, atol=1e-9)
        np.testing.assert_array_equal(self._slack_flags(), ~sat)
        self.assertGreater(int(sat.sum()), 0, "this test wants at least one saturation")

        # second call, the other way: the saturated units are out of the slack now
        still_in = ~sat
        targets2 = self._targets()
        report2 = self.model.redistribute_active_power(-25.)
        exp2, sat2, remaining2 = olf_distribute(targets2[still_in], weights[still_in],
                                                min_p[still_in], max_p[still_in], -25.)
        self.assertEqual(report2.nb_participants, int(still_in.sum()))
        np.testing.assert_allclose(self._targets()[still_in], exp2, atol=1e-9)
        np.testing.assert_allclose(self._targets()[~still_in], targets2[~still_in], atol=1e-12)
        self.assertAlmostEqual(report2.not_distributed_mw, remaining2, places=9)
        V = self._ac_pf()
        self.assertGreater(V.shape[0], 0)

    def test_zero_crossing_generator(self):
        # a generator with min_p < 0 < max_p (a pumped-storage machine, say) injecting a little:
        # a negative mismatch stops it at 0 MW, not at its min_p
        targets = self._targets()
        small = int(np.argmin(np.where(targets > 0., targets, np.inf)))  # the smallest producer
        self.assertGreater(targets[small], 0.)
        min_p = np.full(self.n_gen, -np.inf)
        max_p = np.full(self.n_gen, np.inf)
        min_p[small] = -100.
        max_p[small] = 100.
        self.model.set_gen_p_limits(min_p, max_p)
        mismatch = -(targets[small] * self.n_gen + 10.)  # its equal share would take it below 0
        report = self.model.redistribute_active_power(mismatch)
        new = self._targets()
        self.assertAlmostEqual(new[small], 0., places=9)
        self.assertFalse(self._slack_flags()[small], "stopped at 0 MW, the unit leaves the slack")
        exp, sat, remaining = olf_distribute(targets, np.ones(self.n_gen), min_p, max_p, mismatch)
        self.assertTrue(sat[small])
        # (the synchronous condensers, at 0 MW already, take no share of a negative mismatch)
        self.assertEqual(report.nb_saturated, int(sat.sum()))
        np.testing.assert_allclose(new, exp, atol=1e-9)
        self.assertAlmostEqual(new.sum(), targets.sum() + mismatch, places=9)

    def test_zero_crossing_storage(self):
        # a storage unit whose range straddles 0: discharging, a negative mismatch stops it at
        # 0 MW; charging, a positive one does. Its [min_p, max_p] alone would let it through.
        for p_mw, mismatch in [(-5., -40.), (5., 40.)]:  # pandapower: p_mw > 0 is charging
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                net = pn.case14()
                pp.create_storage(net, bus=3, p_mw=p_mw, max_e_mwh=100., min_p_mw=-20., max_p_mw=20.)
                model = init_from_pandapower(net)
            if len(model.get_storages()) == 0:
                self.skipTest("this pandapower converter has no storage unit")
            for gen in model.get_generators():
                model.add_gen_slackbus(gen.id, 1.)
            model.add_storage_slackbus(0, 1.)
            n_gen = len(model.get_generators())
            model.set_gen_p_limits(np.full(n_gen, -1e4), np.full(n_gen, 1e4))
            model.set_storage_p_limits(np.array([-20.]), np.array([20.]))
            gen_targets = np.array([g.target_p_mw for g in model.get_generators()])
            inj_sto = -model.get_storages()[0].target_p_mw  # generator convention
            report = model.redistribute_active_power(mismatch)
            sto = model.get_storages()[0]
            self.assertAlmostEqual(sto.target_p_mw, 0., places=9, msg=f"storage at {p_mw} MW, mismatch {mismatch}")
            self.assertFalse(sto.is_slack)
            exp, sat, remaining = olf_distribute(np.concatenate((gen_targets, [inj_sto])), np.ones(n_gen + 1),
                                                 np.concatenate((np.full(n_gen, -1e4), [-20.])),
                                                 np.concatenate((np.full(n_gen, 1e4), [20.])), mismatch)
            self.assertTrue(sat[-1])
            self.assertEqual(report.nb_saturated, int(sat.sum()))
            np.testing.assert_allclose(np.array([g.target_p_mw for g in model.get_generators()]),
                                       exp[:n_gen], atol=1e-9)
            V = model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)
            self.assertGreater(V.shape[0], 0)

    def test_redistribute_active_power_with_storage(self):
        # a battery takes part in the slack like a generator, in generator convention:
        # its target_p is in load convention, its limits in generator convention
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            net = pn.case14()
            pp.create_storage(net, bus=3, p_mw=-5., max_e_mwh=100., min_p_mw=-20., max_p_mw=20.)
            model = init_from_pandapower(net)
        if len(model.get_storages()) == 0:
            self.skipTest("this pandapower converter has no storage unit")
        for gen in model.get_generators():
            model.add_gen_slackbus(gen.id, 1.)
        model.add_storage_slackbus(0, 1.)
        n_gen = len(model.get_generators())
        sto = model.get_storages()[0]
        inj_sto = -sto.target_p_mw  # 5 MW injected
        # generator limits: huge; the battery cannot inject more than 8 MW
        model.set_gen_p_limits(np.full(n_gen, -1e4), np.full(n_gen, 1e4))
        model.set_storage_p_limits(np.array([-20.]), np.array([8.]))
        gen_targets = np.array([g.target_p_mw for g in model.get_generators()])

        report = model.redistribute_active_power(50.)
        inj = np.concatenate((gen_targets, [inj_sto]))
        exp, sat, remaining = olf_distribute(inj, np.ones(n_gen + 1),
                                             np.concatenate((np.full(n_gen, -1e4), [-20.])),
                                             np.concatenate((np.full(n_gen, 1e4), [8.])), 50.)
        self.assertTrue(sat[-1], "the battery must saturate at 8 MW")
        self.assertEqual(report.nb_participants, n_gen + 1)
        self.assertEqual(report.nb_saturated, 1)
        np.testing.assert_allclose(np.array([g.target_p_mw for g in model.get_generators()]),
                                   exp[:n_gen], atol=1e-9)
        sto = model.get_storages()[0]
        self.assertAlmostEqual(-sto.target_p_mw, 8., places=9)  # load convention in the model
        self.assertFalse(sto.is_slack)
        V = model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)


if __name__ == "__main__":
    unittest.main()
