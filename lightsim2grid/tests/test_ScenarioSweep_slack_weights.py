# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
ScenarioSweep.set_gen_slack_weights / set_storage_slack_weights: each row gives some units
distributed-slack weights of their own.

The oracle is the grid itself: a row must give what the grid gives once it holds the row's
weights (LSGrid.add_gen_slackbus / remove_gen_slackbus). The voltages are compared up to a
common rotation: the batch keeps the angle reference it picked for the whole sweep, while a
grid taking its first slack bus out of the slack moves its reference elsewhere.
"""

import itertools
import unittest
import warnings

import numpy as np
import pandapower.networks as pn
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.scenarioSweep import ScenarioSweep, ScenarioSweepCPP, LimitViolationType
from lightsim2grid.lightsim2grid_cpp import DistributedSlack, LSGrid, ViolationElementType
from lightsim2grid.network import init_from_pandapower

try:
    from lightsim2grid.algorithm import NROuter_KLU, NROuterDist_KLU  # noqa: F401
    OUTER = "NROuter_KLU"
    OUTER_DIST = "NROuterDist_KLU"
except ImportError:
    OUTER = "NROuter_SparseLU"
    OUTER_DIST = "NROuterDist_SparseLU"


def _gap_up_to_rotation(a, b):
    """largest gap between two voltage vectors once their best common rotation is taken out"""
    rot = np.vdot(b, a)
    rot = rot / abs(rot) if abs(rot) > 0. else 1.
    return np.max(np.abs(a - b * rot))


class TestScenarioSweepSlackWeightsCPP(unittest.TestCase):
    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.Vinit = self.env.backend.V
        self.max_it = self.env.backend.max_it
        self.tol = 1e-10
        # every generator shares the grid's slack
        self.grid = self.env.backend._grid.copy()
        self.grid.change_algorithm(AlgorithmType.NR_KLU)
        self.n_gen = len(self.grid.get_generators())
        for gen_id in range(self.n_gen):
            self.grid.add_gen_slackbus(gen_id, 1. + gen_id)

        # 5 rows: row 0 changes nothing (all NaN), row 3 takes three generators out of the
        # slack (the reference bus's among them)
        self.nb_rows = 5
        self.gen_ids = np.array([0, 2, 4])
        self.weights = np.full((self.nb_rows, 3), np.nan)
        self.weights[1] = [10., 0., 1.]
        self.weights[2, 1] = 0.
        self.weights[3] = [0., 0., 0.]
        self.weights[4, 0] = 50.

    def tearDown(self):
        self.env.close()

    def _sweep(self, grid=None, weights=True, algo=None):
        sweep = ScenarioSweepCPP(self.grid if grid is None else grid)
        sweep.compute_gen_results = True
        if algo is not None:
            sweep.change_algorithm(algo)
        if weights:
            sweep.set_gen_slack_weights(self.gen_ids, self.weights)
        return sweep

    def _grid_of_row(self, row, grid=None):
        """the grid holding row `row`'s weights, given the way LSGrid takes them"""
        res = (self.grid if grid is None else grid).copy()
        for k, gen_id in enumerate(self.gen_ids):
            w = self.weights[row, k]
            if np.isnan(w):
                continue
            if w > 0.:
                res.add_gen_slackbus(int(gen_id), float(w))
            else:
                res.remove_gen_slackbus(int(gen_id))
        return res

    def _check_rows_against_grid(self, ac):
        algo = AlgorithmType.NR_KLU if ac else AlgorithmType.DC_KLU
        sweep = self._sweep(algo=algo)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.all(sweep.converged_mask()))
        Vs = sweep.get_voltages()
        P = sweep.get_gen_results()[:, :, 0]
        for row in range(self.nb_rows):
            ref = self._grid_of_row(row)
            V_ref = ref.ac_pf(self.Vinit, self.max_it, self.tol) if ac else ref.dc_pf(self.Vinit, self.max_it, self.tol)
            self.assertEqual(V_ref.size, Vs.shape[1])
            used = (np.abs(V_ref) > 0.) & (np.abs(Vs[row]) > 0.)
            self.assertLess(_gap_up_to_rotation(Vs[row][used], V_ref[used]), 1e-9, f"row {row}")
            p_ref = np.array([el.res_p_mw for el in ref.get_generators()])
            np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0, err_msg=f"row {row}")
            if row > 0:
                # the weights did act
                self.assertGreater(np.max(np.abs(P[row] - P[0])), 1e-3, f"row {row}")

    def test_rows_match_the_grid_given_their_weights_ac(self):
        self._check_rows_against_grid(ac=True)

    def test_rows_match_the_grid_given_their_weights_dc(self):
        self._check_rows_against_grid(ac=False)

    def test_one_analysis_and_threads(self):
        sweep = self._sweep()
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_linear_solver_stats().nb_analyze, 1)
        Vs = 1. * sweep.get_voltages()
        threaded = self._sweep()
        threaded.nb_thread = 2
        threaded.compute(self.Vinit, self.max_it, self.tol)
        for stats in threaded.get_linear_solver_stats_per_algo():
            self.assertLessEqual(stats.nb_analyze, 1)
        np.testing.assert_allclose(threaded.get_voltages(), Vs, atol=1e-12, rtol=0)

    def test_a_row_changing_nothing_is_the_plain_row(self):
        """a row whose cells are all NaN, and a sweep cleared of its weights, are the plain
        sweep, bit for bit"""
        gen_p = np.tile([el.target_p_mw for el in self.grid.get_generators()], (self.nb_rows, 1))
        plain = self._sweep()
        plain.set_gen_slack_weights(np.array([], dtype=np.int32), np.zeros((0, 0)))  # an empty call clears
        plain.modify_gen_p(gen_p)
        plain.compute(self.Vinit, self.max_it, self.tol)
        sweep = self._sweep()
        sweep.modify_gen_p(gen_p)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertGreater(np.max(np.abs(sweep.get_gen_results()[1] - plain.get_gen_results()[1])), 1e-3)
        np.testing.assert_array_equal(sweep.get_voltages()[0], plain.get_voltages()[0])
        np.testing.assert_array_equal(sweep.get_gen_results()[0], plain.get_gen_results()[0])

    def test_with_injections_and_generator_contingencies(self):
        """a row also moving the set-points and tripping a generator it gives a weight: the
        tripped one takes no share"""
        target_p = np.array([el.target_p_mw for el in self.grid.get_generators()])
        gen_p = np.tile(target_p, (self.nb_rows, 1))
        gen_p[1] *= 0.8
        gen_p[4, 1] += 10.
        mask = np.zeros((self.nb_rows, self.n_gen), dtype=bool)
        mask[1, 0] = True   # weighted 10 in that row
        mask[4, 2] = True   # weighted by the grid only
        sweep = self._sweep()
        sweep.modify_gen_p(gen_p)
        sweep.set_contingency_gens(mask)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.all(sweep.converged_mask()))
        P = sweep.get_gen_results()[:, :, 0]
        for row in range(self.nb_rows):
            ref = self._grid_of_row(row)
            for gen_id in range(self.n_gen):
                ref.change_p_gen(gen_id, float(gen_p[row, gen_id]))
            for gen_id in np.nonzero(mask[row])[0]:
                ref.deactivate_gen(int(gen_id))
            V_ref = ref.ac_pf(self.Vinit, self.max_it, self.tol)
            self.assertEqual(V_ref.size, sweep.get_voltages().shape[1])
            p_ref = np.array([el.res_p_mw for el in ref.get_generators()])
            np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0, err_msg=f"row {row}")
        self.assertEqual(P[1, 0], 0.)

    def test_gen_p_violations_read_the_row_weights(self):
        """the active-power limit check shares the slack with the row's weights: a limit
        between what generator 0 produces in two rows is broken in one only"""
        grid = self.grid.copy()
        max_p = np.full(self.n_gen, 1000.)
        weights = np.array([[np.nan], [100.]])
        probe = ScenarioSweepCPP(grid)
        probe.compute_gen_results = True
        probe.set_gen_slack_weights(np.array([0]), weights)
        probe.compute(self.Vinit, self.max_it, self.tol)
        p0 = probe.get_gen_results()[:, 0, 0]
        self.assertGreater(abs(p0[0] - p0[1]), 0.5)
        max_p[0] = 0.5 * (p0[0] + p0[1])
        grid.set_gen_p_limits(np.full(self.n_gen, -1000.), max_p)
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_gen_results = True
        sweep.compute_physical_violations = True
        sweep.set_gen_slack_weights(np.array([0]), weights)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        broken = [[v for v in viols if v.violation_type in (LimitViolationType.HIGH_P, LimitViolationType.LOW_P)
                   and v.element_id == 0] for viols in sweep.get_physical_violations()]
        high_row = int(np.argmax(p0))
        self.assertEqual(len(broken[high_row]), 1)
        self.assertEqual(len(broken[1 - high_row]), 0)
        self.assertAlmostEqual(broken[high_row][0].value, sweep.get_gen_results()[high_row, 0, 0], places=8)

    def test_gen_p_violations_ignore_a_weight_the_solve_ignores(self):
        """issue #317: a row weight below the threshold at which the solve takes a unit out of the
        slack (BaseConstants::_tol_equal_float) solves exactly as 0, but the active-power check
        still counted the unit as sharing the slack, and reported a generator whose target is
        already above its maximum. The 4-bus feeder of the GenPCheck tests (80 MW on bus 3),
        generators 0 and 1 sharing the slack, generator 2 (bus 0) out of it with its target above
        its maximum; the rows give it 0, 1e-8 and 1."""
        grid = LSGrid()
        grid.set_sn_mva(100.)
        grid.set_init_vm_pu(1.0)
        grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
        grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                             np.array([0, 1, 2]), np.array([1, 2, 3]))
        grid.init_loads(np.array([80.]), np.array([60.]), np.array([3]))
        grid.init_generators(np.array([0., 10., 5.]), np.array([1.02, 1.05, 1.02]),
                             np.full(3, -1e3), np.full(3, 1e3), np.array([0, 1, 0]))
        grid.add_gen_slackbus(0, 1.)
        grid.add_gen_slackbus(1, 1.)
        grid.set_gen_p_limits(np.zeros(3), np.array([30., 30., 3.]))
        grid.tell_solver_need_reset()
        V0 = grid.ac_pf(np.ones(4, dtype=complex), 20, 1e-10)
        self.assertEqual(V0.size, 4)
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_physical_violations = True
        sweep.physical_violation_tol_mva = 0.
        sweep.modify_load_p(np.full((3, 1), 80.))
        sweep.set_gen_slack_weights(np.array([2]), np.array([[0.], [1e-8], [1.]]))
        sweep.compute(V0, 20, 1e-10)
        self.assertTrue(np.all(sweep.converged_mask()))
        Vs = sweep.get_voltages()
        np.testing.assert_array_equal(Vs[1], Vs[0])  # the solve ignores the 1e-8
        gen_2 = [[v.violation_type for v in viols
                  if v.element_type == ViolationElementType.GENERATOR and v.element_id == 2]
                 for viols in sweep.get_physical_violations()]
        self.assertEqual(gen_2, [[], [], [LimitViolationType.HIGH_P]])

    def test_setters_check_their_input(self):
        sweep = self._sweep(weights=False)
        for ids, w in ((np.array([self.n_gen]), np.ones((2, 1))),        # out of range
                       (np.array([0, 0]), np.ones((2, 2))),             # named twice
                       (np.array([0]), np.ones((2, 2))),                # shape
                       (np.array([0]), np.full((2, 1), -1.)),           # negative
                       (np.array([0]), np.full((2, 1), np.inf))):       # infinite
            with self.assertRaises(ValueError):
                sweep.set_gen_slack_weights(ids, w)
        # nothing was registered: another count of rows is still free ...
        sweep.set_gen_slack_weights(np.array([0]), np.full((3, 1), np.nan))
        # ... until now
        with self.assertRaises(RuntimeError):
            sweep.set_gen_slack_weights(np.array([0]), np.full((2, 1), np.nan))

    def test_what_compute_refuses(self):
        # a positive weight to a generator on a bus outside the grid's slack
        grid = self.env.backend._grid.copy()   # only its slack generator shares the slack
        grid.change_algorithm(AlgorithmType.NR_KLU)
        sweep = ScenarioSweepCPP(grid)
        sweep.set_gen_slack_weights(np.array([0]), np.array([[np.nan], [1.]]))
        with self.assertRaises(ValueError):
            sweep.compute(self.Vinit, self.max_it, self.tol)
        # ... where 0 changes nothing to the set of buses
        sweep.set_gen_slack_weights(np.array([0]), np.array([[np.nan], [0.]]))
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.all(sweep.converged_mask()))

    def test_outer_loops_solve_the_row(self):
        """under outer loops the row is solved (TestScenarioSweepSlackWeightsOuter checks how)"""
        sweep = self._sweep(weights=False)
        sweep.compute_limit_violations = True
        sweep.change_algorithm(OUTER)
        sweep.set_gen_slack_weights(np.array([0]), np.array([[np.nan], [3.]]))
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertEqual(list(sweep.get_batch_unsupported_reasons()), ["", ""])
        self.assertEqual(list(sweep.converged_mask()), [True, True])


class TestScenarioSweepSlackWeightsRedistribute(unittest.TestCase):
    """with redistribute_slack: the bounded pre-pass shares what a row's injections take out
    of the balance with the row's weights. The reference is the one-off path of
    test_batch_redistribute_slack: a copy of the grid holding the row's weights and
    set-points, LSGrid.redistribute_active_power of what the row loses, then the powerflow."""

    MAX_IT = 30
    TOL = 1e-10

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(pn.case14())
        self.n_gen = len(self.grid.get_generators())
        for gen_id in range(self.n_gen):
            self.grid.add_gen_slackbus(gen_id, 1.)
        self.targets = np.array([el.target_p_mw for el in self.grid.get_generators()])
        self.load_targets = np.array([el.target_p_mw for el in self.grid.get_loads()])
        # generator 1 can only move 2 MW up: a heavy weight saturates it
        min_p = np.full(self.n_gen, -np.inf)
        max_p = np.full(self.n_gen, np.inf)
        max_p[1] = self.targets[1] + 2.
        self.grid.set_gen_p_limits(min_p, max_p)
        self.V0 = np.full(self.grid.total_bus(), 1.04 + 0j)
        # every row loses 30 MW (a heavier load); the weights decide who makes it up
        self.nb_rows = 5
        self.load_p = np.tile(self.load_targets, (self.nb_rows, 1))
        self.load_p[:, 2] += 30.
        self.gen_ids = np.array([1, 2, 3])
        self.weights = np.full((self.nb_rows, 3), np.nan)
        self.weights[1] = [10., np.nan, np.nan]   # generator 1 saturates
        self.weights[2] = [np.nan, 0., 5.]        # generator 2 out, 3 heavier
        self.weights[3] = [0., 0., 0.]
        self.weights[4] = [0.1, 3., np.nan]       # generator 1 light: stays inside

    def _sweep(self, dc):
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(AlgorithmType.DC_KLU if dc else AlgorithmType.NR_KLU)
        sweep.compute_gen_results = True
        sweep.redistribute_slack = True
        sweep.modify_load_p(self.load_p)
        sweep.set_gen_slack_weights(self.gen_ids, self.weights)
        sweep.compute(self.V0, self.MAX_IT, self.TOL)
        return sweep

    def _reference(self, row, dc):
        ref = self.grid.copy()
        ref.change_algorithm(AlgorithmType.DC_KLU if dc else AlgorithmType.NR_KLU)
        for k, gen_id in enumerate(self.gen_ids):
            w = self.weights[row, k]
            if np.isnan(w):
                continue
            if w > 0.:
                ref.add_gen_slackbus(int(gen_id), float(w))
            else:
                ref.remove_gen_slackbus(int(gen_id))
        for load_id, p in enumerate(self.load_p[row]):
            ref.change_p_load(load_id, float(p))
        lost = -float(np.sum(self.load_targets - self.load_p[row]))
        report = ref.redistribute_active_power(lost)
        V = ref.dc_pf(self.V0, self.MAX_IT, self.TOL) if dc else ref.ac_pf(self.V0, self.MAX_IT, self.TOL)
        self.assertGreater(V.shape[0], 0, "the one-off reference diverged")
        return V, np.array([el.res_p_mw for el in ref.get_generators()]), report

    def _check(self, dc):
        sweep = self._sweep(dc)
        self.assertTrue(np.all(sweep.converged_mask()))
        Vs = sweep.get_voltages()
        P = sweep.get_gen_results()[:, :, 0]
        saturated = []
        for row in range(self.nb_rows):
            with self.subTest(row=row):
                V_ref, p_ref, report = self._reference(row, dc)
                saturated.append(report.nb_saturated)
                self.assertLess(_gap_up_to_rotation(Vs[row], V_ref), 1e-8)
                np.testing.assert_allclose(P[row], p_ref, atol=1e-6, rtol=0)
        # the weights decide who saturates: generator 1 in row 1 only
        self.assertEqual(saturated[1], 1)
        self.assertEqual(saturated[4], 0)
        self.assertAlmostEqual(P[1, 1], self.targets[1] + 2., places=6)
        # ... and what the others take: not the grid's split
        self.assertGreater(np.max(np.abs(P[2] - P[0])), 1.)

    def test_matches_one_off_ac(self):
        self._check(dc=False)

    def test_matches_one_off_dc(self):
        self._check(dc=True)


class TestScenarioSweepSlackWeightsOuter(unittest.TestCase):
    """under outer loops the DistributedSlack loop shares the slack with the row's weights. The
    grid of test_outer_batch: case14, its generators short of reactive power, every one of them
    able to take part in the slack ("can participate"), generator 1 only 5 MW below its maximum.
    The reference is that grid holding the row's weights (set_gen_can_participate_slack, the
    weights the loop shares with) and loads, solved once by the same algorithm."""

    MAX_IT = 30
    TOL = 1e-8

    def _grid(self, algo, wide_q=False, tight=False):
        net = pn.case14()
        net.gen["min_q_mvar"] = -1000. if wide_q else -15.
        net.gen["max_q_mvar"] = 1000. if wide_q else 15.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(net)
        n_gen = len(grid.get_generators())
        targets = np.array([el.target_p_mw for el in grid.get_generators()])
        max_p = np.full(n_gen, 300.)
        max_p[1] = targets[1] + 5.
        grid.set_gen_p_limits(np.zeros(n_gen), max_p)
        grid.set_gen_can_participate_slack([True] * n_gen, np.full(n_gen, 300.))
        grid.change_algorithm(algo)
        if tight:
            # the loop stops within its tolerance of the slack bus' mismatch: tight, two solves
            # starting from different voltages land on the same state
            loops = list(grid.get_outer_loops())
            grid.clear_outer_loops()
            for loop in loops:
                grid.add_outer_loop(DistributedSlack(slack_bus_p_max_mismatch_mw=1e-6)
                                    if loop.name() == "DistributedSlack" else loop)
        return grid

    def setUp(self):
        grid = self._grid(OUTER)
        self.n_gen = len(grid.get_generators())
        self.load_targets = np.array([el.target_p_mw for el in grid.get_loads()])
        self.V0 = np.full(grid.total_bus(), 1.04 + 0j)
        # every row loses 40 MW (a heavier load), shared by the loop with the row's weights
        self.nb_rows = 5
        self.load_p = np.tile(self.load_targets, (self.nb_rows, 1))
        self.load_p[:, 3] += 40.
        self.gen_ids = np.array([1, 2, 3])
        self.weights = np.full((self.nb_rows, 3), np.nan)
        self.weights[1] = [3000., np.nan, np.nan]   # generator 1 stops at its maximum
        self.weights[2] = [np.nan, 0., 900.]
        self.weights[3] = [0., 0., 0.]
        self.weights[4] = [30., 600., np.nan]

    def _sweep(self, grid, init_from_n=False):
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_gen_results = True
        sweep.init_from_n_powerflow = init_from_n
        sweep.modify_load_p(self.load_p)
        sweep.set_gen_slack_weights(self.gen_ids, self.weights)
        sweep.compute(self.V0, self.MAX_IT, self.TOL)
        self.assertTrue(np.all(sweep.converged_mask()))
        self.assertEqual(set(sweep.get_batch_unsupported_reasons()), {""})
        return sweep

    def _reference(self, make_grid, row, V_start, uniform_start=True):
        ref = make_grid()
        ref.set_uniform_voltage_start(uniform_start)
        w = np.full(self.n_gen, 300.)
        for k, gen_id in enumerate(self.gen_ids):
            if not np.isnan(self.weights[row, k]):
                w[gen_id] = self.weights[row, k]
        ref.set_gen_can_participate_slack([bool(x > 0.) for x in w], np.where(w > 0., w, 0.))
        for load_id, p in enumerate(self.load_p[row]):
            ref.change_p_load(load_id, float(p))
        V = ref.ac_pf(V_start, self.MAX_IT, self.TOL)
        self.assertGreater(V.shape[0], 0, "the one-off reference diverged")
        return V, np.array([el.res_p_mw for el in ref.get_generators()])

    def _check(self, algo):
        sweep = self._sweep(self._grid(algo))
        Vs = sweep.get_voltages()
        P = sweep.get_gen_results()[:, :, 0]
        for row in range(self.nb_rows):
            with self.subTest(row=row):
                V_ref, p_ref = self._reference(lambda: self._grid(algo), row, self.V0)
                self.assertLess(_gap_up_to_rotation(Vs[row], V_ref), 1e-9)
                np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0)
        # generator 1 stops at its maximum however heavy its weight; the weights move the rest
        self.assertAlmostEqual(P[1, 1], P[0, 1], places=6)
        self.assertGreater(np.max(np.abs(P[2] - P[0])), 1.)
        self.assertGreater(np.max(np.abs(P[4] - P[0])), 1.)

    def test_rows_match_the_grid_given_their_weights(self):
        self._check(OUTER)

    def test_rows_match_the_grid_given_their_weights_dist(self):
        self._check(OUTER_DIST)

    def test_seeded_from_n_a_row_with_weights_starts_from_the_grid_targets(self):
        """with init_from_n_powerflow a row giving weights does not resume the "n" solve's
        loops: it is the grid holding its weights solved from the "n" voltages, its loops starting
        from the grid's targets. Wide reactive limits and a tight DistributedSlack tolerance, so
        that the voltages the solve starts from decide nothing else."""
        for algo in (OUTER, OUTER_DIST):
            with self.subTest(algo=algo):
                make = lambda: self._grid(algo, wide_q=True, tight=True)  # noqa: E731
                sweep = self._sweep(make(), init_from_n=True)
                V_n = make().ac_pf(self.V0, self.MAX_IT, self.TOL)
                P = sweep.get_gen_results()[:, :, 0]
                for row in range(1, self.nb_rows):
                    V_ref, p_ref = self._reference(make, row, V_n, uniform_start=False)
                    self.assertLess(_gap_up_to_rotation(sweep.get_voltages()[row], V_ref), 1e-9, f"row {row}")
                    np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0, err_msg=f"row {row}")

    def test_seeded_from_n_a_row_with_weights_starts_from_the_set_points(self):
        """issue #316: with init_from_n_powerflow a row giving weights does not resume the "n"
        solve's loops, so its generator buses start PV again -- but from the "n" solve's voltages,
        where the buses that solve froze at a reactive limit sit off their set-point, and a PV bus
        keeps the magnitude it starts from. With less load the row's machines stay within their
        limits, nothing freezes those buses again, and the row converged holding the frozen
        magnitude instead of the target. Here the narrow reactive limits make the "n" solve freeze
        buses, and the row only gives generator 0 the weight it already has. The reference is the
        grid at the row's load, solved from the "n" voltages; a row giving set-points of its own
        (modify_gen_v) keeps those. Both row loops: the plain one and the masked one
        (handle_disconnected_grid)."""
        for algo in (OUTER, OUTER_DIST):
            make = lambda: self._grid(algo, tight=True)  # noqa: E731
            V_n = make().ac_pf(self.V0, self.MAX_IT, self.TOL)
            targets = np.array([el.target_vm_pu for el in make().get_generators()])
            buses = np.array([el.bus_id for el in make().get_generators()])
            # the premise: the "n" solve left some generator bus off its set-point
            self.assertGreater(np.max(np.abs(np.abs(V_n[buses]) - targets)), 1e-3)
            load_p = np.tile(0.6 * self.load_targets, (2, 1))
            gen_v = np.tile(targets, (2, 1))
            gen_v[1, 2] += 0.01
            for with_gen_v, masked in itertools.product((False, True), (False, True)):
                with self.subTest(algo=algo, with_gen_v=with_gen_v, masked=masked):
                    sweep = ScenarioSweepCPP(make())
                    sweep.init_from_n_powerflow = True
                    sweep.handle_disconnected_grid = masked
                    sweep.modify_load_p(load_p)
                    if with_gen_v:
                        sweep.modify_gen_v(gen_v)
                    # row 0 resumes the "n" solve's loops, row 1 gives weights
                    sweep.set_gen_slack_weights(np.array([0]), np.array([[np.nan], [300.]]))
                    sweep.compute(self.V0, self.MAX_IT, self.TOL)
                    self.assertTrue(np.all(sweep.converged_mask()))
                    for row in range(2):
                        ref = make()
                        ref.set_uniform_voltage_start(False)
                        for load_id, p in enumerate(load_p[row]):
                            ref.change_p_load(load_id, float(p))
                        if with_gen_v:
                            for gen_id, v in enumerate(gen_v[row]):
                                ref.change_v_gen(gen_id, float(v))
                        V_ref = ref.ac_pf(V_n, self.MAX_IT, self.TOL)
                        self.assertGreater(V_ref.shape[0], 0, "the one-off reference diverged")
                        self.assertLess(_gap_up_to_rotation(sweep.get_voltages()[row], V_ref), 1e-9, f"row {row}")

    def test_a_unit_outside_the_grid_slack_is_refused(self):
        grid = self._grid(OUTER)
        # generator 2 can no longer take part in the grid's slack
        flags = [True] * self.n_gen
        flags[2] = False
        grid.set_gen_can_participate_slack(flags, np.where(flags, 300., 0.))
        sweep = ScenarioSweepCPP(grid)
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_gen_slack_weights(np.array([2]), np.array([[np.nan], [10.]]))
        with self.assertRaises(ValueError):
            sweep.compute(self.V0, self.MAX_IT, self.TOL)
        # ... and 0 keeps it out, which it already is
        sweep.set_gen_slack_weights(np.array([2]), np.array([[np.nan], [0.]]))
        sweep.compute(self.V0, self.MAX_IT, self.TOL)
        self.assertTrue(np.all(sweep.converged_mask()))


class TestScenarioSweepSlackWeights(unittest.TestCase):
    """the python wrapper: ids or names"""
    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_names(self):
        gens = self.env.backend._grid.get_generators()
        # the grid's slack generator re-weighted, generator 0 kept out of it
        slack_gen = [el.id for el in gens if el.is_slack][0]
        weights = np.array([[2., 0.], [np.nan, 0.], [5., np.nan]])
        by_id = ScenarioSweep(self.env)
        by_id.set_gen_slack_weights([slack_gen, 0], weights)
        by_id.compute()
        by_name = ScenarioSweep(self.env)
        by_name.set_gen_slack_weights([gens[slack_gen].name, gens[0].name], weights)
        by_name.compute()
        self.assertEqual(by_id.get_voltages().shape[0], 3)
        np.testing.assert_array_equal(by_id.get_voltages(), by_name.get_voltages())
        with self.assertRaises(RuntimeError):
            by_id.set_gen_slack_weights(["not a generator"], np.ones((3, 1)))
        with self.assertRaises(RuntimeError):
            by_id.set_gen_slack_weights([0], np.ones(3))  # not a matrix


if __name__ == "__main__":
    unittest.main()
