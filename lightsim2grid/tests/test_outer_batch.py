# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The outer loops in the batch classes: a batch whose algorithm is an NROuter_* one solves
each row with the grid's outer loops, and gets what a single solve of that row's grid gets --
with one symbolic analysis and one factorization for the whole batch."""

import os
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP, InjectionSweepCPP,
                                             ScenarioSweepCPP, TimeSeriesCPP)

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER = "NROuter_KLU"
except ImportError:
    OUTER = "NROuter_SparseLU"

MAX_IT = 30
TOL = 1e-8
N_ROWS = 8


def _grid(distributed=True):
    """case14, its generators short of reactive power (ReactiveLimits switches buses PQ) and,
    unless not `distributed`, sharing the slack in proportion to their maximum output
    (DistributedSlack re-shares whatever a row changes in the balance)"""
    net = pn.case14()
    net.gen["min_q_mvar"] = -15.
    net.gen["max_q_mvar"] = 15.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(net)
    nb_gen = len(grid.get_generators())
    grid.set_gen_p_limits(np.zeros(nb_gen), np.full(nb_gen, 300.))
    if distributed:
        grid.set_gen_can_participate_slack([True] * nb_gen, np.full(nb_gen, 300.))
    grid.change_algorithm(OUTER)
    return grid


def _v0(grid):
    return np.full(grid.total_bus(), 1.04 + 0j)


def _plain(algo):
    """the Newton-Raphson an outer-loop algorithm wraps, as a plain algorithm (its batch counts
    the analyses and factorizations the outer one must not exceed)"""
    return algo.replace("NROuterDist", "NR").replace("NROuter", "NR")


def _resume_on(grid):
    """a grid whose next solves start their loops from where the last one left them"""
    cfg = grid.get_ac_algo_config()
    ip = list(cfg.int_params)
    ip[5] = 1  # resume_outer_state, after the Newton's 5
    cfg.int_params = ip
    grid.set_ac_algo_config(cfg)
    return grid


def _from_n(make_grid, change_row):
    """a row's reference when the batch starts it from its "n" solve: the "n" grid solved, then
    the row's change, solved from the "n" voltages AND equations (the buses its loops froze,
    ...). Returns the grid, to solve, and the voltages to start from."""
    grid = _resume_on(make_grid())
    v_n = grid.ac_pf(_v0(grid), MAX_IT, TOL)
    change_row(grid)
    return grid, v_n


def _branch_results(grid):
    """P and Q at both ends of every branch after a single solve, as the batches publish
    them: shape (nb_branch, 4), the lines then the trafos"""
    res = []
    for side1, side2 in ((grid.get_line_res1(), grid.get_line_res2()),
                         (grid.get_trafo_res1(), grid.get_trafo_res2())):
        res.append(np.stack([side1[0], side1[1], side2[0], side2[1]], axis=1))
    return np.concatenate(res, axis=0)


class _Equivalence:
    """each row of the batch equals a single solve of the row's own grid"""

    # whether row i starts from row i - 1's voltages (TimeSeries) rather than from the seed
    CHAINED = False
    # the algorithm of the batch and of the single solves
    ALGO = OUTER

    def _grid(self, distributed=True):
        grid = _grid(distributed)
        grid.change_algorithm(self.ALGO)
        return grid

    def _check_rows(self, make_batch, row_grids, nb_thread=1, v_starts=None):
        """`v_starts`: the voltages each row's reference starts from (the batch's seed, or the
        previous row's for a chained one, by default)"""
        batch = make_batch()
        if nb_thread > 1:
            batch.nb_thread = nb_thread
        batch.change_algorithm(self.ALGO)
        batch.compute_gen_results = True
        batch.compute(_v0(self.base), MAX_IT, TOL)
        Vs = batch.get_voltages()
        branches = batch.compute_branch_results()
        gens = batch.get_gen_results()
        converged = np.asarray(batch.converged_mask(), dtype=bool)
        nb_ok = 0
        for i, grid in enumerate(row_grids):
            # the loops stop on thresholds (DistributedSlack leaves less than its threshold
            # on the slack), so where they stop depends on where the solve starts: the same
            # start as the row's
            if v_starts is not None:
                v_start = v_starts[i]
            else:
                v_start = Vs[i - 1] if self.CHAINED and i > 0 else _v0(grid)
            if self.CHAINED and i > 0:
                # a chained row also starts from the previous row's equations: the reference
                # resumes from the previous row's solve (see resume_outer_state)
                grid = self._resumed(row_grids[0], grid)
            V = grid.ac_pf(v_start, MAX_IT, TOL)
            V = V if V.size else None
            with self.subTest(row=i):
                self.assertEqual(converged[i], V is not None)
                if V is None:
                    continue
                nb_ok += 1
                np.testing.assert_allclose(Vs[i], V, rtol=0., atol=1e-10)
                # the flows (at the taps the loops settled on) and the generators' output
                # (at the targets DistributedSlack moved them to)
                np.testing.assert_allclose(branches[i], _branch_results(grid), rtol=0., atol=1e-6)
                p, q, _ = grid.get_gen_res()
                np.testing.assert_allclose(gens[i, :, 0], p, rtol=0., atol=1e-6)
                np.testing.assert_allclose(gens[i, :, 1], q, rtol=0., atol=1e-6)
        self.assertGreater(nb_ok, 0)
        # the outer loops add no symbolic analysis and no factorization to what the plain
        # Newton-Raphson batch does (one of each per thread, and one for the base solve),
        # but the fallback ones, which are counted
        plain = make_batch()
        if nb_thread > 1:
            plain.nb_thread = nb_thread
        plain.change_algorithm(_plain(self.ALGO))
        plain.compute(_v0(self.base), MAX_IT, TOL)
        stats = batch.get_linear_solver_stats()
        ref = plain.get_linear_solver_stats()
        self.assertEqual(stats.nb_analyze, ref.nb_analyze)
        self.assertEqual(stats.nb_factorize - stats.nb_fallback_factorize,
                         ref.nb_factorize - ref.nb_fallback_factorize)
        self.assertGreater(stats.nb_refactorize, 0)


class TestOuterContingencyAnalysis(_Equivalence, unittest.TestCase):
    def setUp(self):
        self.base = self._grid()

    def _rows(self):
        rows = []
        for line_id in range(N_ROWS):
            grid = self._grid()
            grid.deactivate_powerline(line_id)
            rows.append(grid)
        return rows

    def _batch(self):
        batch = ContingencyAnalysisCPP(self.base)
        for line_id in range(N_ROWS):
            batch.add_n1(line_id)
        return batch

    def test_each_row_is_the_single_solve_of_its_grid(self):
        self._check_rows(self._batch, self._rows())

    def test_two_threads(self):
        self._check_rows(self._batch, self._rows(), nb_thread=2)


class TestOuterScenarioSweep(_Equivalence, unittest.TestCase):
    def setUp(self):
        self.base = self._grid()
        gen_p = np.array([g.target_p_mw for g in self.base.get_generators()])
        load_p = np.asarray(self.base.get_load_target_p())
        scale = 1. + 0.05 * np.arange(N_ROWS)
        self.gen_p = np.outer(np.ones(N_ROWS), gen_p)
        self.load_p = np.outer(scale, load_p)
        self.line_mask = np.zeros((N_ROWS, len(self.base.get_lines())), dtype=bool)
        for i in range(N_ROWS):
            self.line_mask[i, i + 1] = True

    def _rows(self):
        rows = []
        for i in range(N_ROWS):
            grid = self._grid()
            for load_id, p in enumerate(self.load_p[i]):
                grid.change_p_load(load_id, float(p))
            grid.deactivate_powerline(i + 1)
            rows.append(grid)
        return rows

    def _batch(self):
        batch = ScenarioSweepCPP(self.base)
        batch.modify_gen_p(self.gen_p)
        batch.modify_load_p(self.load_p)
        batch.set_contingency_lines(self.line_mask)
        return batch

    def test_each_row_is_the_single_solve_of_its_grid(self):
        self._check_rows(self._batch, self._rows())

    def test_two_threads(self):
        self._check_rows(self._batch, self._rows(), nb_thread=2)


class _InjectionRows(_Equivalence):
    """a load increase per row"""

    def setUp(self):
        self.base = self._grid()
        load_p = np.asarray(self.base.get_load_target_p())
        self.load_p = np.outer(1. + 0.05 * np.arange(N_ROWS), load_p)
        self.gen_p = np.outer(np.ones(N_ROWS), [g.target_p_mw for g in self.base.get_generators()])

    def _rows(self):
        rows = []
        for i in range(N_ROWS):
            grid = self._grid()
            for load_id, p in enumerate(self.load_p[i]):
                grid.change_p_load(load_id, float(p))
            rows.append(grid)
        return rows

    def _fill(self, batch):
        batch.modify_gen_p(self.gen_p)
        batch.modify_load_p(self.load_p)
        return batch


class TestOuterInjectionSweep(_InjectionRows, unittest.TestCase):
    def test_each_row_is_the_single_solve_of_its_grid(self):
        self._check_rows(lambda: self._fill(InjectionSweepCPP(self.base)), self._rows())

    def test_two_threads(self):
        self._check_rows(lambda: self._fill(InjectionSweepCPP(self.base)), self._rows(), nb_thread=2)


class TestOuterTimeSeries(_InjectionRows, unittest.TestCase):
    CHAINED = True

    @staticmethod
    def _resumed(previous, grid):
        """`grid` solved by `previous`'s algorithm, with what its loops last left (a grid has
        one algorithm: `previous` takes the row's loads instead)"""
        for load_id, p in enumerate(np.asarray(grid.get_load_target_p())):
            previous.change_p_load(load_id, float(p))
        return previous

    def _rows(self):
        # one grid solved row after row, its loops resuming from the previous solve
        grid = _resume_on(self._grid())
        rows = []
        for i in range(N_ROWS):
            row = self._grid()
            for load_id, p in enumerate(self.load_p[i]):
                row.change_p_load(load_id, float(p))
            rows.append(row)
        rows[0] = grid
        for load_id, p in enumerate(self.load_p[0]):
            grid.change_p_load(load_id, float(p))
        return rows

    def test_each_row_is_the_single_solve_of_its_grid(self):
        # each row starts from the previous one's voltages AND its equations (the buses its
        # loops froze, ...): a sequence of solves resuming from one another
        self._check_rows(lambda: self._fill(TimeSeriesCPP(self.base)), self._rows())

    def test_rows_resume_from_the_previous_one(self):
        # resuming, a row needs fewer outer iterations than starting the loops afresh
        ts = self._fill(TimeSeriesCPP(self.base))
        ts.change_algorithm(self.ALGO)
        ts.compute(_v0(self.base), MAX_IT, TOL)
        resumed = [r.nb_outer_iterations for r in ts.get_outer_results()]
        fresh = self._fill(InjectionSweepCPP(self.base))
        fresh.change_algorithm(self.ALGO)
        fresh.compute(_v0(self.base), MAX_IT, TOL)
        restarted = [r.nb_outer_iterations for r in fresh.get_outer_results()]
        self.assertLess(sum(resumed[1:]), sum(restarted[1:]))


class _FromN(_Equivalence):
    """init_from_n_powerflow: every row starts from the "n" solve -- its voltages, and the
    equations its loops left (frozen buses, taps, ...) -- not from the grid's"""

    def _make_grid(self):
        return self._grid()

    def _change_row(self, i, grid):
        raise NotImplementedError

    def _make(self):
        raise NotImplementedError

    def _batch_from_n(self):
        batch = self._make()
        batch.init_from_n_powerflow = True
        return batch

    def _references(self):
        refs = [_from_n(self._make_grid, lambda g, i=i: self._change_row(i, g)) for i in range(N_ROWS)]
        return [g for g, _ in refs], [v for _, v in refs]

    def test_each_row_is_the_single_solve_resumed_from_n(self):
        grids, v_starts = self._references()
        self._check_rows(self._batch_from_n, grids, v_starts=v_starts)

    def test_each_row_is_the_single_solve_resumed_from_n_two_threads(self):
        grids, v_starts = self._references()
        self._check_rows(self._batch_from_n, grids, nb_thread=2, v_starts=v_starts)

    def test_rows_need_fewer_outer_iterations(self):
        # what the "n" solve's loops did is not redone: fewer outer iterations than each row
        # starting its loops afresh from the "n" voltages
        batch = self._batch_from_n()
        batch.change_algorithm(self.ALGO)
        batch.compute(_v0(self.base), MAX_IT, TOL)
        resumed = sum(r.nb_outer_iterations for r in batch.get_outer_results())
        afresh = 0
        for i in range(N_ROWS):
            grid = self._make_grid()
            v_n = grid.ac_pf(_v0(grid), MAX_IT, TOL)
            self._change_row(i, grid)
            grid.ac_pf(v_n, MAX_IT, TOL)
            afresh += grid.get_algo().get_outer_loop_stats().nb_outer_iterations
        self.assertLess(resumed, afresh)


class TestOuterContingencyAnalysisFromN(_FromN, unittest.TestCase):
    # no DistributedSlack: a row that keeps the grid's targets starts from where the "n" solve's
    # moved them, which no single solve does (see TestOuterRowThatChangesNothing)
    def _make_grid(self):
        return self._grid(distributed=False)

    def setUp(self):
        self.base = self._make_grid()

    def _change_row(self, i, grid):
        grid.deactivate_powerline(i)

    def _make(self):
        batch = ContingencyAnalysisCPP(self.base)
        for line_id in range(N_ROWS):
            batch.add_n1(line_id)
        return batch


class TestOuterScenarioSweepFromN(_FromN, unittest.TestCase):
    def setUp(self):
        TestOuterScenarioSweep.setUp(self)

    def _change_row(self, i, grid):
        for load_id, p in enumerate(self.load_p[i]):
            grid.change_p_load(load_id, float(p))
        grid.deactivate_powerline(i + 1)

    def _make(self):
        return TestOuterScenarioSweep._batch(self)


class TestOuterInjectionSweepFromN(_FromN, _InjectionRows, unittest.TestCase):
    def _change_row(self, i, grid):
        for load_id, p in enumerate(self.load_p[i]):
            grid.change_p_load(load_id, float(p))

    def _make(self):
        return self._fill(InjectionSweepCPP(self.base))


class TestOuterRowThatChangesNothing(unittest.TestCase):
    """a row seeded from the "n" solve that changes nothing IS the "n" solve: its loops have
    nothing to do -- DistributedSlack included, its units starting from the targets the "n"
    solve moved them to, where the row does not set its own"""

    def test_it_is_the_n_solve(self):
        base = _grid()
        grid = _grid()
        v_n = grid.ac_pf(_v0(grid), MAX_IT, TOL)
        self.assertGreater(grid.get_algo().get_outer_loop_stats().nb_outer_iterations, 0)
        p_n, q_n, _ = grid.get_gen_res()
        sweep = InjectionSweepCPP(base)
        # the loads only: the rows keep the grid's generator targets
        sweep.modify_load_p(np.outer(np.ones(N_ROWS), base.get_load_target_p()))
        sweep.init_from_n_powerflow = True
        sweep.change_algorithm(OUTER)
        sweep.compute_gen_results = True
        sweep.compute(_v0(base), MAX_IT, TOL)
        gens = sweep.get_gen_results()
        for i, (res, V) in enumerate(zip(sweep.get_outer_results(), sweep.get_voltages())):
            with self.subTest(row=i):
                self.assertEqual(res.nb_outer_iterations, 0)
                np.testing.assert_allclose(V, v_n, rtol=0., atol=1e-10)
                np.testing.assert_allclose(gens[i, :, 0], p_n, rtol=0., atol=1e-6)
                np.testing.assert_allclose(gens[i, :, 1], q_n, rtol=0., atol=1e-6)


def _tap_grid():
    """IEEE 14 (pypowsybl's) with a ratio tap changer regulating the low voltage side of a
    transformer, solved with ReactiveLimits and TransformerVoltageControl: rows that open a
    line move the tap, so their flows are those of the tap the row settled on"""
    from test_transformer_voltage_control import _net, _grid as _tvc_grid
    from lightsim2grid.algorithm import ReactiveLimits, TransformerVoltageControl
    grid = _tvc_grid(_net("T4-9-1", 12.6, 0.1), reactive_limits=True)
    grid.clear_outer_loops()
    grid.add_outer_loop(ReactiveLimits())
    # its generators frozen while it acts (the default leaves this grid's tap where it is)
    grid.add_outer_loop(TransformerVoltageControl(use_initial_tap_position=True, max_controlled_nominal_voltage=120.))
    grid.change_algorithm(OUTER)
    return grid


class TestOuterTapsInContingencyAnalysis(_Equivalence, unittest.TestCase):
    def setUp(self):
        try:
            import pypowsybl  # noqa: F401
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        self.base = _tap_grid()

    def _rows(self):
        rows = []
        for line_id in range(N_ROWS):
            grid = _tap_grid()
            grid.deactivate_powerline(line_id)
            rows.append(grid)
        return rows

    def _batch(self):
        batch = ContingencyAnalysisCPP(self.base)
        for line_id in range(N_ROWS):
            batch.add_n1(line_id)
        return batch

    def test_each_row_is_the_single_solve_of_its_grid(self):
        self._check_rows(self._batch, self._rows())
        # the tap moves somewhere, or this test checks nothing of the taps
        batch = self._batch()
        batch.change_algorithm(OUTER)
        batch.compute(_v0(self.base), MAX_IT, TOL)
        self.assertTrue(any(r.ratio_tap for r in batch.get_outer_results()))


class TestOuterTripsTheControlledTransformer(_Equivalence, unittest.TestCase):
    """a row that opens the transformer TransformerVoltageControl regulates with: the row has
    no such transformer, whatever tap the "n" solve left it at -- its Ybus must not keep any
    part of it, nor its result list a tap for it (the "n" solve moves that tap, and a row
    resuming from it carries the positions the "n" solve settled on)"""

    TRAFO = 1  # T4-9-1, the regulating one
    LINE = 1   # a line whose opening the grid survives

    def setUp(self):
        try:
            import pypowsybl  # noqa: F401
        except ImportError:
            self.skipTest("pypowsybl is not installed")
        self.base = _tap_grid()
        self.nb_line = len(self.base.get_lines())

    # the rows in the batch's order: a contingency analysis sorts its contingencies by branch
    # id, the lines then the trafos
    def _change_row(self, i, grid):
        if i == 0:
            grid.deactivate_powerline(self.LINE)
        else:
            grid.deactivate_trafo(self.TRAFO)

    def _batch(self):
        batch = ContingencyAnalysisCPP(self.base)
        batch.add_n1(self.LINE)
        batch.add_n1(self.nb_line + self.TRAFO)
        return batch

    def _rows(self):
        rows = []
        for i in range(2):
            grid = _tap_grid()
            self._change_row(i, grid)
            rows.append(grid)
        return rows

    def _batch_from_n(self):
        batch = self._batch()
        batch.init_from_n_powerflow = True
        return batch

    def test_the_n_solve_moves_the_tap(self):
        # or the resumed rows would carry the grid's own tap, and check nothing
        grid = _tap_grid()
        grid.ac_pf(_v0(grid), MAX_IT, TOL)
        tr = grid.get_trafos()[self.TRAFO]
        self.assertNotEqual(tr.res_ratio_tap_position, tr.ratio_tap_position)

    def test_each_row_is_the_single_solve_of_its_grid(self):
        self._check_rows(self._batch, self._rows())

    def test_each_row_is_the_single_solve_resumed_from_n(self):
        refs = [_from_n(_tap_grid, lambda g, i=i: self._change_row(i, g)) for i in range(2)]
        self._check_rows(self._batch_from_n, [g for g, _ in refs], v_starts=[v for _, v in refs])

    def test_a_resumed_single_solve_has_no_such_transformer(self):
        # the reference itself: a solve resuming from the "n" one, the transformer opened in
        # between, solves the grid without it -- its power balance is the one of its own Ybus
        grid, v_n = _from_n(_tap_grid, lambda g: self._change_row(1, g))
        V = grid.ac_pf(v_n, MAX_IT, TOL)
        self.assertGreater(V.size, 0)
        mis = V * np.conj(grid.get_Ybus_solver() @ V) - np.asarray(grid.get_Sbus_solver())
        pq = np.asarray(grid.get_pq_solver())
        self.assertLess(np.abs(mis[pq]).max(), 1e-6)

    def test_a_second_single_solve_has_no_such_transformer(self):
        # not resumed either: a grid solved, then the transformer opened and solved again
        grid = _tap_grid()
        v_n = grid.ac_pf(_v0(grid), MAX_IT, TOL)
        self._change_row(1, grid)
        V = grid.ac_pf(v_n, MAX_IT, TOL)
        self.assertGreater(V.size, 0)
        mis = V * np.conj(grid.get_Ybus_solver() @ V) - np.asarray(grid.get_Sbus_solver())
        pq = np.asarray(grid.get_pq_solver())
        self.assertLess(np.abs(mis[pq]).max(), 1e-6)

    def test_no_tap_for_the_opened_transformer(self):
        batch = self._batch_from_n()
        batch.change_algorithm(OUTER)
        batch.compute(_v0(self.base), MAX_IT, TOL)
        res = batch.get_outer_results()
        # the line row keeps the tap the "n" solve moved, the trafo row has no such trafo
        self.assertIn(self.TRAFO, [t for t, _ in res[0].ratio_tap])
        self.assertNotIn(self.TRAFO, [t for t, _ in res[1].ratio_tap])


class TestOuterThreadsWhoseFirstRowStrandsABus(unittest.TestCase):
    """A worker's outer loops reserve their controls (the buses ReactiveLimits may switch,
    ...) for the base case, as the member algorithm does with its "n" solve, not for
    whichever row the worker happens to solve first: the results do not depend on the
    number of threads."""

    # trafo 18 alone feeds bus 7, whose generator sits at its reactive limit: its row strands
    # them. Two threads split the four rows two and two, so the second worker starts with it.
    BRANCHES = (0, 1, 18, 19)
    STRANDED_GEN = 3

    def _run(self, nb_thread):
        base = _grid()
        ca = ContingencyAnalysisCPP(base)
        ca.handle_disconnected_grid = True
        for branch in self.BRANCHES:
            ca.add_n1(branch)
        ca.nb_thread = nb_thread
        ca.compute_gen_results = True
        ca.compute(_v0(base), MAX_IT, TOL)
        self.assertEqual([list(c) for c in ca.my_defaults()], [[b] for b in self.BRANCHES])
        return ca.get_voltages(), ca.get_gen_results(), np.asarray(ca.converged_mask(), dtype=bool)

    def test_same_results_with_two_threads(self):
        V1, gens1, conv1 = self._run(1)
        V2, gens2, conv2 = self._run(2)
        np.testing.assert_array_equal(conv1, conv2)
        self.assertTrue(conv1[3])
        # the last row has bus 7 back: its generator is held at its limit again
        self.assertLessEqual(abs(gens2[3, self.STRANDED_GEN, 1]), 15. + 1e-6)
        np.testing.assert_allclose(V2[conv1], V1[conv1], rtol=0., atol=1e-10)
        np.testing.assert_allclose(gens2[conv1], gens1[conv1], rtol=0., atol=1e-6)


class TestWrappersKeepTheOuterAlgorithm(unittest.TestCase):
    """the grid2op wrappers keep a backend's NROuter_* algorithm (any other one is replaced by
    NR_KLU, as before)"""

    def _env(self, algo_type=None):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            kwargs = {} if algo_type is None else {"algo_type": algo_type}
            return grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend(**kwargs))

    def test_contingency_analysis_and_time_series(self):
        from lightsim2grid import ContingencyAnalysis, TimeSerie
        env = self._env(OUTER)
        self.assertEqual(ContingencyAnalysis(env).computer.get_algo_name(), OUTER)
        self.assertEqual(TimeSerie(env).computer.get_algo_name(), OUTER)
        env.close()
        env = self._env()
        self.assertNotEqual(ContingencyAnalysis(env).computer.get_algo_name(), OUTER)
        env.close()


if __name__ == "__main__":
    unittest.main()
