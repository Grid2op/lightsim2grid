# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``redistribute_slack`` on ContingencyAnalysis / ScenarioSweep.

With it, the active power a row loses (a generator contingency, an island cut off with
``handle_disconnected_grid``) is shared on the remaining slack units BEFORE the solve,
OpenLoadFlow-style, with their ``[min_p, max_p]`` bounds, the saturated units leaving
that row's distributed slack. The batch must then land exactly where the one-off path
does: a copy of the grid, the elements really removed,
``consider_only_main_component(True)`` / ``redistribute_active_power``, ``ac_pf``.
"""

import unittest
import warnings
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
    from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP, ScenarioSweepCPP
    from lightsim2grid.algorithm import AlgorithmType

# bus 7 of case14 is a leaf (degree 1) with a generator (synchronous condenser) on it
_LEAF_BUS = 7
_LEAF_P_MW = 40.
_MAX_IT = 30
_TOL = 1e-10


def _angles_rel(V, ref_bus=0):
    return np.angle(V) - np.angle(V[ref_bus])


class _Base(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.net = pn.case14()
            self.grid = init_from_pandapower(self.net)
        gens = self.grid.get_generators()
        self.n_gen = len(gens)
        self.n_line = len(self.grid.get_lines())
        self.leaf_gen = [g.id for g in gens if g.bus_id == _LEAF_BUS][0]
        # the (only) branch to the leaf bus, in ContingencyAnalysis numbering (lines then trafos)
        branch = [l.id for l in self.grid.get_lines() if _LEAF_BUS in (l.bus1_id, l.bus2_id)]
        branch += [self.n_line + t.id for t in self.grid.get_trafos() if _LEAF_BUS in (t.bus1_id, t.bus2_id)]
        assert len(branch) == 1
        self.leaf_branch = branch[0]
        # distributed slack on every generator, the leaf one producing something
        for g in gens:
            self.grid.add_gen_slackbus(g.id, 1.)
        self.grid.change_p_gen(self.leaf_gen, _LEAF_P_MW)
        self.targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        self.solved = np.array([b for b in range(self.V0.shape[0]) if b != _LEAF_BUS])

    def _limits_two_clamped(self, extra=(3., 3.)):
        """two of the units that stay can only take `extra` MW more"""
        others = [g for g in range(self.n_gen) if g != self.leaf_gen]
        max_p = np.full(self.n_gen, np.inf)
        min_p = np.full(self.n_gen, -np.inf)
        for g, e in zip(others[:2], extra):
            max_p[g] = self.targets[g] + e
        self.grid.set_gen_p_limits(min_p, max_p)
        return others[:2], max_p

    def _one_off_island(self, dc=False):
        ref = self.grid.copy()
        if self.leaf_branch < self.n_line:
            ref.deactivate_powerline(self.leaf_branch)
        else:
            ref.deactivate_trafo(self.leaf_branch - self.n_line)
        report = ref.consider_only_main_component(True)
        V = ref.dc_pf(1. * self.V0, _MAX_IT, _TOL) if dc else ref.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        assert V.shape[0] > 0, "the one-off reference diverged"
        return V, report, ref

    def _assert_same_state(self, V_batch, V_ref, buses, dc=False):
        if not dc:
            np.testing.assert_allclose(np.abs(V_batch[buses]), np.abs(V_ref[buses]), rtol=0., atol=1e-6)
        # the angle reference may differ between the two paths (the batch picks its own
        # reference slack bus; the one-off moves it when the first slack unit saturates)
        np.testing.assert_allclose(_angles_rel(V_batch)[buses], _angles_rel(V_ref)[buses], rtol=0., atol=1e-6)


class TestContingencyAnalysisRedistributeSlack(_Base):
    def _ca(self, dc=False, redistribute=True, physical=False):
        ca = ContingencyAnalysisCPP(self.grid)
        if dc:
            ca.change_algorithm(AlgorithmType.DC_SparseLU)
        ca.add_n1(self.leaf_branch)
        ca.handle_disconnected_grid = True
        ca.redistribute_slack = redistribute
        if physical:
            ca.compute_physical_violations = True
        ca.compute(1. * self.V0, _MAX_IT, _TOL)
        self.assertTrue(ca.converged_mask()[0], "the contingency did not converge")
        return ca

    def test_flag_default_off(self):
        ca = ContingencyAnalysisCPP(self.grid)
        self.assertFalse(ca.redistribute_slack)
        ca.redistribute_slack = True
        self.assertTrue(ca.redistribute_slack)

    def test_matches_one_off_no_limits(self):
        ca = self._ca()
        V_ref, report, _ = self._one_off_island()
        self.assertAlmostEqual(report.mismatch_mw, _LEAF_P_MW, places=9)
        self._assert_same_state(ca.get_voltages()[0], V_ref, self.solved)
        self.assertEqual(np.abs(ca.get_voltages()[0][_LEAF_BUS]), 0.)

    def test_matches_one_off_with_saturation(self):
        clamped, max_p = self._limits_two_clamped()
        ca = self._ca(physical=True)
        V_ref, report, ref = self._one_off_island()
        self.assertEqual(report.nb_saturated, 2)
        self._assert_same_state(ca.get_voltages()[0], V_ref, self.solved)
        # the saturated units sit at max_p and are not reported above it
        viol = ca.get_physical_violations()[0]
        self.assertEqual([v for v in viol if "HIGH_P" in str(v.violation_type)], [])
        self.assertEqual([v for v in viol if "LOW_P" in str(v.violation_type)], [])

    def test_matches_one_off_with_saturation_dc(self):
        self._limits_two_clamped()
        ca = self._ca(dc=True, physical=True)
        V_ref, report, _ = self._one_off_island(dc=True)
        self.assertEqual(report.nb_saturated, 2)
        self._assert_same_state(ca.get_voltages()[0], V_ref, self.solved, dc=True)
        viol = ca.get_physical_violations()[0]
        self.assertEqual([v for v in viol if "HIGH_P" in str(v.violation_type)], [])

    def test_reference_slack_saturates(self):
        # the batch's reference slack bus is the first generator's (bus 0, fewest
        # strandings, equal weights): make it the one that saturates
        max_p = np.full(self.n_gen, np.inf)
        max_p[0] = self.targets[0] + 1.
        self.grid.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        ca = self._ca()
        V_ref, report, _ = self._one_off_island()
        self.assertEqual(report.nb_saturated, 1)
        self._assert_same_state(ca.get_voltages()[0], V_ref, self.solved)

    def test_off_is_bit_identical_to_before(self):
        self._limits_two_clamped()
        ca_off = self._ca(redistribute=False)
        ca_plain = ContingencyAnalysisCPP(self.grid)
        ca_plain.add_n1(self.leaf_branch)
        ca_plain.handle_disconnected_grid = True
        ca_plain.compute(1. * self.V0, _MAX_IT, _TOL)
        np.testing.assert_array_equal(ca_off.get_voltages(), ca_plain.get_voltages())
        # ... and differs from the redistributed one (the limits bite)
        ca_on = self._ca(redistribute=True)
        self.assertGreater(np.max(np.abs(ca_on.get_voltages()[0] - ca_off.get_voltages()[0])), 1e-6)

    def test_nothing_lost_is_unchanged(self):
        # a contingency that islands nothing: the option changes nothing at all
        self._limits_two_clamped()
        other_branch = 0 if self.leaf_branch != 0 else 1
        res = []
        for redistribute in (False, True):
            ca = ContingencyAnalysisCPP(self.grid)
            ca.add_n1(other_branch)
            ca.handle_disconnected_grid = True
            ca.redistribute_slack = redistribute
            ca.compute(1. * self.V0, _MAX_IT, _TOL)
            res.append(ca.get_voltages())
        np.testing.assert_array_equal(res[0], res[1])


class TestContingencyAnalysisHvdcIsland(unittest.TestCase):
    """the island cut off holds an HVDC converter station: its setpoint is part of what
    the row loses (the one-off path counts it too, the two must agree)"""
    _PSP = 30.

    def setUp(self):
        import os
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _aux_make_hvdc import make_case14_hvdc
        # side 2 (the leaf bus) rectifies: it draws _PSP from the grid it is cut off from
        self.net, self.grid = make_case14_hvdc(3, _LEAF_BUS, converters_mode=1, p_setpoint=self._PSP)
        gens = self.grid.get_generators()
        self.n_gen = len(gens)
        n_line = len(self.grid.get_lines())
        branch = [l.id for l in self.grid.get_lines() if _LEAF_BUS in (l.bus1_id, l.bus2_id)]
        branch += [n_line + t.id for t in self.grid.get_trafos() if _LEAF_BUS in (t.bus1_id, t.bus2_id)]
        assert len(branch) == 1
        self.leaf_branch = branch[0]
        self.n_line = n_line
        for g in gens:
            self.grid.add_gen_slackbus(g.id, 1.)
        # a bound that the share of the lost consumption reaches on one unit
        targets = np.array([g.target_p_mw for g in gens])
        min_p = np.full(self.n_gen, -np.inf)
        producing = [g.id for g in gens if g.target_p_mw > 0. and g.bus_id != _LEAF_BUS]
        min_p[producing[0]] = targets[producing[0]] - 2.
        self.grid.set_gen_p_limits(min_p, np.full(self.n_gen, np.inf))
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        self.solved = np.array([b for b in range(self.V0.shape[0]) if b != _LEAF_BUS])

    def test_matches_one_off(self):
        ref = self.grid.copy()
        if self.leaf_branch < self.n_line:
            ref.deactivate_powerline(self.leaf_branch)
        else:
            ref.deactivate_trafo(self.leaf_branch - self.n_line)
        report = ref.consider_only_main_component(True)
        self.assertAlmostEqual(report.mismatch_mw, -self._PSP, places=9)
        # the clamped unit, plus the synchronous condensers (0 MW, they cannot go below)
        at_zero = sum(1 for g in self.grid.get_generators() if g.target_p_mw == 0. and g.bus_id != _LEAF_BUS)
        self.assertEqual(report.nb_saturated, 1 + at_zero)
        V_ref = ref.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V_ref.shape[0], 0, "the one-off reference diverged")

        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_n1(self.leaf_branch)
        ca.handle_disconnected_grid = True
        ca.redistribute_slack = True
        ca.compute(1. * self.V0, _MAX_IT, _TOL)
        self.assertTrue(ca.converged_mask()[0], "the contingency did not converge")
        V_batch = ca.get_voltages()[0]
        np.testing.assert_allclose(np.abs(V_batch[self.solved]), np.abs(V_ref[self.solved]), rtol=0., atol=1e-6)
        np.testing.assert_allclose(_angles_rel(V_batch)[self.solved], _angles_rel(V_ref)[self.solved], rtol=0., atol=1e-6)
        # ... and it is not what the unbounded slack gives
        ca_off = ContingencyAnalysisCPP(self.grid)
        ca_off.add_n1(self.leaf_branch)
        ca_off.handle_disconnected_grid = True
        ca_off.compute(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(np.max(np.abs(ca_off.get_voltages()[0] - V_batch)), 1e-6)


class TestScenarioSweepRedistributeSlack(_Base):
    def setUp(self):
        super().setUp()
        # a non-slack generator with some power too (gen on bus 2), for the "lost power
        # of a machine that took no part in the slack" case
        self.non_slack = [g.id for g in self.grid.get_generators() if g.bus_id == 2][0]
        self.grid.remove_gen_slackbus(self.non_slack)
        self.grid.change_p_gen(self.non_slack, 20.)
        self.targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        self.slack_gen = [g.id for g in self.grid.get_generators() if g.bus_id == 1][0]
        self.all_buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _reference(self, gen_off, row_p=None):
        ref = self.grid.copy()
        if row_p is not None:
            for g in range(self.n_gen):
                ref.change_p_gen(g, float(row_p[g]))
        lost = float(row_p[gen_off]) if row_p is not None else float(self.targets[gen_off])
        ref.deactivate_gen(gen_off)
        report = ref.redistribute_active_power(lost)
        V = ref.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        assert V.shape[0] > 0, "the one-off reference diverged"
        return V, report

    def _sweep(self, gen_mask, gen_p=None, redistribute=True):
        sweep = ScenarioSweepCPP(self.grid)
        sweep.set_contingency_gens(gen_mask)
        if gen_p is not None:
            sweep.modify_gen_p(gen_p)
        sweep.redistribute_slack = redistribute
        sweep.compute(1. * self.V0, _MAX_IT, _TOL)
        return sweep

    def test_slack_gen_off_with_saturation(self):
        others = [g for g in range(self.n_gen) if g not in (self.slack_gen, self.non_slack)]
        max_p = np.full(self.n_gen, np.inf)
        max_p[others[0]] = self.targets[others[0]] + 2.
        self.grid.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, self.slack_gen] = True
        sweep = self._sweep(gen_mask)
        self.assertTrue(sweep.converged_mask()[0])
        V_ref, report = self._reference(self.slack_gen)
        self.assertEqual(report.nb_saturated, 1)
        np.testing.assert_allclose(sweep.get_voltages()[0][self.all_buses], V_ref[self.all_buses],
                                   rtol=0., atol=1e-6)

    def test_non_slack_gen_off(self):
        others = [g for g in range(self.n_gen) if g != self.non_slack]
        max_p = np.full(self.n_gen, np.inf)
        max_p[others[0]] = self.targets[others[0]] + 1.
        self.grid.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, self.non_slack] = True
        sweep = self._sweep(gen_mask)
        self.assertTrue(sweep.converged_mask()[0])
        V_ref, report = self._reference(self.non_slack)
        self.assertAlmostEqual(report.mismatch_mw, 20.)
        self.assertEqual(report.nb_saturated, 1)
        np.testing.assert_allclose(sweep.get_voltages()[0][self.all_buses], V_ref[self.all_buses],
                                   rtol=0., atol=1e-6)

    def test_row_setpoints_are_the_rows_own(self):
        # modify_gen_p: the lost power (and the units' starting points) are the ROW's
        others = [g for g in range(self.n_gen) if g not in (self.slack_gen, self.non_slack)]
        max_p = np.full(self.n_gen, np.inf)
        max_p[others[0]] = self.targets[others[0]] + 2.
        self.grid.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        gen_p = np.vstack((self.targets, 1.1 * self.targets))
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)
        gen_mask[:, self.slack_gen] = True
        sweep = self._sweep(gen_mask, gen_p=gen_p)
        for row in range(2):
            with self.subTest(row=row):
                self.assertTrue(sweep.converged_mask()[row])
                V_ref, _ = self._reference(self.slack_gen, row_p=gen_p[row])
                np.testing.assert_allclose(sweep.get_voltages()[row][self.all_buses], V_ref[self.all_buses],
                                           rtol=0., atol=1e-6)

    def test_off_is_bit_identical_to_before(self):
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, self.slack_gen] = True
        sweep_off = self._sweep(gen_mask, redistribute=False)
        sweep_plain = ScenarioSweepCPP(self.grid)
        sweep_plain.set_contingency_gens(gen_mask)
        sweep_plain.compute(1. * self.V0, _MAX_IT, _TOL)
        np.testing.assert_array_equal(sweep_off.get_voltages(), sweep_plain.get_voltages())


if __name__ == "__main__":
    unittest.main()
