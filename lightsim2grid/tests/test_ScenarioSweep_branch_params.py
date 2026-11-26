# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
ScenarioSweep.set_line_parameters / set_trafo_parameters: each row gives some branches
parameters of their own.

The oracle is the grid itself: a row must give what the grid gives once it holds the row's
values (LSGrid.update_powerlines_parameters / update_trafos_parameters_full /
change_ratio_trafo / change_shift_trafo). A row that also changes the injections, or trips a
branch, is compared with a ScenarioSweep of that modified grid playing the same row -- that
part of the row is not what is under test here.
"""

import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.scenarioSweep import ScenarioSweep, ScenarioSweepCPP
from test_ContingencyAnalysis_limit_violations import _set_tight_limits
from test_ScenarioSweep_exotic import (build_grid as build_exotic_grid, HALF_OPEN_LINE, HALF_OPEN_TRAFO,
                                       MAX_IT, TOL)

NAN_C = complex(np.nan, np.nan)


def _nan(shape):
    return np.full(shape, np.nan)


def _nan_c(shape):
    return np.full(shape, NAN_C, dtype=complex)


def _grid_holding_row(grid_base, row, line_ids, line_mats, trafo_ids, trafo_mats, trafo_ratio, trafo_shift):
    """a copy of `grid_base` holding row `row`'s values, given the way LSGrid takes them
    (`line_mats` / `trafo_mats`: the r, x, h1, h2 matrices; a NaN keeps the grid's value)"""
    grid = grid_base.copy()
    for el_ids, mats, elements, update in (
            (line_ids, line_mats, grid_base.get_lines(), grid.update_powerlines_parameters),
            (trafo_ids, trafo_mats, grid_base.get_trafos(), grid.update_trafos_parameters_full)):
        vects = [np.array([getattr(el, attr) for el in elements]) for attr in ("r_pu", "x_pu", "h1_pu", "h2_pu")]
        for k, el_id in enumerate(el_ids):
            for vect, mat in zip(vects, mats):
                if not np.isnan(mat[row, k]):
                    vect[el_id] = mat[row, k]
        update(*vects)
    for k, el_id in enumerate(trafo_ids):
        if not np.isnan(trafo_ratio[row, k]):
            grid.change_ratio_trafo(int(el_id), trafo_ratio[row, k])
        if not np.isnan(trafo_shift[row, k]):
            grid.change_shift_trafo(int(el_id), trafo_shift[row, k])
    return grid


class TestScenarioSweepBranchParamsCPP(unittest.TestCase):
    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.Vinit = self.env.backend.V
        self.max_it = self.env.backend.max_it
        self.tol = 1e-10
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())

        # 6 rows: some lines and transformers changed, row 2 changes nothing (all NaN)
        self.nb_rows = 6
        lines = self.grid.get_lines()
        trafos = self.grid.get_trafos()
        rng = np.random.default_rng(0)
        self.line_ids = np.array([0, 3, 7])
        base_x = np.array([lines[i].x_pu for i in self.line_ids])
        self.line_r = _nan((self.nb_rows, 3))
        self.line_x = base_x * rng.uniform(0.5, 2., (self.nb_rows, 3))
        self.line_h1 = _nan_c((self.nb_rows, 3))
        self.line_h2 = _nan_c((self.nb_rows, 3))
        self.line_r[1, :] = [2. * lines[i].r_pu for i in self.line_ids]
        self.line_h1[3, 1] = 0.01j
        self.line_h2[3, 1] = 0.02j
        self.line_x[2, :] = np.nan
        self.trafo_ids = np.array([0, 2])
        self.trafo_r = _nan((self.nb_rows, 2))
        self.trafo_x = _nan((self.nb_rows, 2))
        self.trafo_h1 = _nan_c((self.nb_rows, 2))
        self.trafo_h2 = _nan_c((self.nb_rows, 2))
        self.trafo_ratio = _nan((self.nb_rows, 2))
        self.trafo_shift = _nan((self.nb_rows, 2))
        self.trafo_x[:, 0] = trafos[0].x_pu * np.linspace(0.6, 1.6, self.nb_rows)
        self.trafo_x[2, 0] = np.nan
        self.trafo_ratio[3, 1] = 1.03
        self.trafo_shift[4, :] = 0.05
        self.trafo_r[5, 1] = 3. * trafos[2].r_pu

    def tearDown(self):
        self.env.close()

    def _sweep(self, grid=None, params=True, algo=None, violations=False):
        sweep = ScenarioSweepCPP(self.grid if grid is None else grid)
        # first: toggling it clears everything registered
        sweep.compute_limit_violations = violations
        if algo is not None:
            sweep.change_algorithm(algo)
        if params:
            sweep.set_line_parameters(self.line_ids, self.line_r, self.line_x, self.line_h1, self.line_h2)
            sweep.set_trafo_parameters(self.trafo_ids, self.trafo_r, self.trafo_x, self.trafo_h1, self.trafo_h2,
                                       self.trafo_ratio, self.trafo_shift)
        return sweep

    def _grid_of_row(self, row):
        """the grid holding row `row`'s values, given the way LSGrid takes them"""
        return _grid_holding_row(self.grid, row,
                                 self.line_ids, (self.line_r, self.line_x, self.line_h1, self.line_h2),
                                 self.trafo_ids, (self.trafo_r, self.trafo_x, self.trafo_h1, self.trafo_h2),
                                 self.trafo_ratio, self.trafo_shift)

    def _check_rows_against_grid(self, ac):
        algo = None if ac else AlgorithmType.DC_KLU
        sweep = self._sweep(algo=algo)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        Vs = 1. * sweep.get_voltages()
        P = 1. * sweep.compute_power_flows()
        A = 1. * sweep.compute_flows()
        for row in range(self.nb_rows):
            grid = self._grid_of_row(row)
            if ac:
                V = grid.ac_pf(self.Vinit, self.max_it, self.tol)
            else:
                V = grid.dc_pf(self.Vinit, self.max_it, self.tol)
            self.assertEqual(V.shape[0], Vs.shape[1], f"row {row}: the grid's solve diverged")
            # the busbars no element uses: 0 in the sweep, the starting point in a solve
            on = np.abs(Vs[row]) > 0.
            self.assertGreater(on.sum(), 0)
            p_ref = np.array([el.res_p1_mw for el in grid.get_lines()] + [el.res_p1_mw for el in grid.get_trafos()])
            a_ref = np.array([el.res_a1_ka for el in grid.get_lines()] + [el.res_a1_ka for el in grid.get_trafos()])
            np.testing.assert_allclose(Vs[row][on], V[on], atol=1e-9, rtol=0, err_msg=f"row {row}")
            np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0, err_msg=f"row {row}")
            if ac:
                np.testing.assert_allclose(A[row], a_ref, atol=1e-8, rtol=0, err_msg=f"row {row}")

    def test_rows_match_the_grid_given_their_values_ac(self):
        self._check_rows_against_grid(ac=True)

    def test_rows_match_the_grid_given_their_values_dc(self):
        """DC: x, the ratio and the shift matter (a shift is a pair of injections)"""
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
        """a row whose cells are all NaN, and a sweep cleared of its parameters, are the
        plain sweep, bit for bit"""
        plain = self._sweep(params=False)
        mask = np.zeros((self.nb_rows, self.n_line), dtype=bool)
        plain.set_contingency_lines(mask)
        plain.compute(self.Vinit, self.max_it, self.tol)
        sweep = self._sweep()
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.array_equal(sweep.get_voltages()[2], plain.get_voltages()[2]))
        sweep.set_line_parameters(np.array([], dtype=np.int32), _nan((0, 0)), _nan((0, 0)), _nan_c((0, 0)), _nan_c((0, 0)))
        sweep.set_trafo_parameters(np.array([], dtype=np.int32), _nan((0, 0)), _nan((0, 0)), _nan_c((0, 0)),
                                   _nan_c((0, 0)), _nan((0, 0)), _nan((0, 0)))
        sweep.compute(self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.array_equal(sweep.get_voltages(), plain.get_voltages()))

    def test_with_injections_and_contingencies(self):
        """a row's parameters together with its injections and its contingency: the
        modified grid playing the same row. A line the row trips ignores its parameters."""
        data = self.env.chronics_handler.real_data.data
        load_p = 1. * data.load_p[:self.nb_rows]
        mask = np.zeros((self.nb_rows, self.n_line), dtype=bool)
        mask[1, 3] = True   # tripped AND given parameters in that row: the trip wins
        mask[4, 10] = True  # tripped, another line given parameters
        sweep = self._sweep()
        sweep.modify_load_p(load_p)
        sweep.set_contingency_lines(mask)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        Vs = 1. * sweep.get_voltages()
        A = 1. * sweep.compute_flows()
        for row in range(self.nb_rows):
            ref = ScenarioSweepCPP(self._grid_of_row(row))
            ref.modify_load_p(load_p[row:row + 1])
            ref.set_contingency_lines(mask[row:row + 1])
            ref.compute(self.Vinit, self.max_it, self.tol)
            np.testing.assert_allclose(Vs[row], ref.get_voltages()[0], atol=1e-9, rtol=0, err_msg=f"row {row}")
            np.testing.assert_allclose(A[row], ref.compute_flows()[0], atol=1e-8, rtol=0, err_msg=f"row {row}")
        self.assertEqual(A[1, 3], 0.)

    def test_violations_and_branch_results(self):
        """the current violations and the branch results read the row's parameters"""
        grid = self.grid.copy()
        _set_tight_limits(grid, cur_ka=0.05)
        self.grid = grid
        sweep = self._sweep(violations=True)
        sweep.compute(self.Vinit, self.max_it, self.tol)
        res = 1. * sweep.compute_branch_results()
        for row in range(self.nb_rows):
            ref = ScenarioSweepCPP(self._grid_of_row(row))
            ref.compute_limit_violations = True
            ref.set_contingency_lines(np.zeros((1, self.n_line), dtype=bool))
            ref.compute(self.Vinit, self.max_it, self.tol)
            np.testing.assert_allclose(res[row], ref.compute_branch_results()[0], atol=1e-7, rtol=0, err_msg=f"row {row}")
            key = lambda v: (int(v.element_type), v.element_id, v.side, int(v.violation_type))  # noqa: E731
            got = sorted(sweep.get_violations()[row], key=key)
            exp = sorted(ref.get_violations()[0], key=key)
            self.assertEqual([key(v) for v in got], [key(v) for v in exp], f"row {row}")
            self.assertGreater(len(got), 0)
            for v_got, v_exp in zip(got, exp):
                self.assertAlmostEqual(v_got.value, v_exp.value, delta=1e-7)

    def test_row_ybus(self):
        sweep = self._sweep()
        sweep.compute(self.Vinit, self.max_it, self.tol)
        for row in range(self.nb_rows):
            grid = self._grid_of_row(row)
            grid.ac_pf(self.Vinit, self.max_it, self.tol)
            diff = (sweep.get_Ybus(row) - grid.get_Ybus()).toarray()
            self.assertLess(np.abs(diff).max(), 1e-10, f"row {row}")

    def test_setters_check_their_input(self):
        sweep = ScenarioSweepCPP(self.grid)
        two = _nan((2, 1))
        two_c = _nan_c((2, 1))
        with self.assertRaises(ValueError):
            sweep.set_line_parameters(np.array([self.n_line]), two, two, two_c, two_c)
        with self.assertRaises(ValueError):
            sweep.set_line_parameters(np.array([0]), two, _nan((3, 1)), two_c, two_c)
        with self.assertRaises(ValueError):
            sweep.set_line_parameters(np.array([0]), np.zeros((2, 1)), np.zeros((2, 1)), two_c, two_c)
        with self.assertRaises(ValueError):
            sweep.set_trafo_parameters(np.array([0]), two, two, two_c, two_c, -np.ones((2, 1)), two)
        # nothing was locked by the refused calls
        sweep.set_line_parameters(np.array([0]), _nan((3, 1)), _nan((3, 1)), _nan_c((3, 1)), _nan_c((3, 1)))
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_lines(np.zeros((2, self.n_line), dtype=bool))
        # a non finite admittance is found when the row is resolved
        x0 = _nan((3, 1))
        x0[1, 0] = 0.
        sweep.set_line_parameters(np.array([0]), np.zeros((3, 1)) + 0.01, x0, _nan_c((3, 1)), _nan_c((3, 1)))
        sweep.change_algorithm(AlgorithmType.DC_KLU)
        with self.assertRaises(ValueError):
            sweep.compute(self.Vinit, self.max_it, self.tol)


class TestScenarioSweepBranchParamsGrid2op(unittest.TestCase):
    """the python wrapper: names, omitted matrices, shapes"""
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), test=True)
        self.grid = self.env.backend._grid

    def tearDown(self):
        self.env.close()

    def test_names_and_omitted_matrices(self):
        lines = self.grid.get_lines()
        x = np.array([[lines[0].x_pu * 1.5, lines[3].x_pu * 0.8]] * 3)
        by_id = ScenarioSweep(self.env)
        by_id.set_line_parameters([0, 3], x=x)
        by_id.compute()
        by_name = ScenarioSweep(self.env)
        by_name.set_line_parameters([lines[0].name, lines[3].name], x=x)
        by_name.compute()
        np.testing.assert_array_equal(by_id.get_voltages(), by_name.get_voltages())
        # the same row through the C++ class, every other matrix NaN
        cpp = ScenarioSweepCPP(self.grid)
        cpp.change_algorithm(by_id.computer.get_algo_type())
        cpp.set_line_parameters(np.array([0, 3]), _nan((3, 2)), x, _nan_c((3, 2)), _nan_c((3, 2)))
        cpp.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        np.testing.assert_allclose(by_id.get_voltages(), cpp.get_voltages(), atol=1e-12, rtol=0)

    def test_errors(self):
        sweep = ScenarioSweep(self.env)
        with self.assertRaises(RuntimeError):
            sweep.set_line_parameters([0])  # nothing given
        with self.assertRaises(RuntimeError):
            sweep.set_line_parameters([0], x=np.ones(3))  # not a matrix
        with self.assertRaises(RuntimeError):
            sweep.set_line_parameters([0], x=np.ones((3, 1)), r=np.ones((2, 1)))  # row counts differ
        with self.assertRaises(RuntimeError):
            sweep.set_trafo_parameters(["not a transformer"], x=np.ones((3, 1)))
        trafos = self.grid.get_trafos()
        sweep.set_trafo_parameters([0], ratio=np.full((3, 1), 1.02), shift_rad=np.full((3, 1), 0.01))
        sweep.compute()
        self.assertEqual(sweep.get_voltages().shape[0], 3)
        self.assertTrue(np.all(sweep.computer.converged_mask()))
        self.assertEqual(len(trafos), len(self.grid.get_trafos()))


class TestScenarioSweepBranchParamsHalfOpen(unittest.TestCase):
    """Parameters given to branches the grid holds open at one end: what such a branch
    stamps is its block with the open end Kron-reduced out, in the Ybus edit, the flows and
    the branch results alike. In DC it is in no matrix, whatever its parameters.

    With no shunt admittance a branch open at one end stamps nothing, whatever its r, x or
    ratio, so the rows that change those also give the branch shunts: otherwise they would
    pass however the parameters were read."""

    SHUNTS = (0.04j, 0.08j, 0.01 + 0.03j, 0.05j)  # line h1, h2, trafo h1, h2

    def setUp(self):
        self.grid = build_exotic_grid()
        lines = self.grid.get_lines()
        trafos = self.grid.get_trafos()
        line = lines[HALF_OPEN_LINE]
        trafo = trafos[HALF_OPEN_TRAFO]
        self.assertTrue(line.connected_global and line.connected1 and not line.connected2)
        self.assertTrue(trafo.connected_global and not trafo.connected1 and trafo.connected2)
        self.V0 = np.ones(self.grid.total_bus(), dtype=complex)
        self.nb_rows = 7
        n = self.nb_rows
        # the half-open line and a line closed at both ends; the half-open trafo
        self.line_ids = np.array([HALF_OPEN_LINE, 3])
        self.trafo_ids = np.array([HALF_OPEN_TRAFO])
        self.line_r, self.line_x = _nan((n, 2)), _nan((n, 2))
        self.line_h1, self.line_h2 = _nan_c((n, 2)), _nan_c((n, 2))
        self.trafo_r, self.trafo_x, self.trafo_ratio, self.trafo_shift = (_nan((n, 1)) for _ in range(4))
        self.trafo_h1, self.trafo_h2 = _nan_c((n, 1)), _nan_c((n, 1))
        # row 0 changes nothing
        self.line_x[1, 0] = 2. * line.x_pu
        self.line_h1[2, 0] = 0.05j                       # the shunts only, the open end's too
        self.line_h2[2, 0] = 0.1j
        self.line_r[3, 0] = 3. * line.r_pu + 0.01
        self.trafo_x[4, 0] = 1.5 * trafo.x_pu
        self.trafo_ratio[4, 0] = 1.05
        self.trafo_shift[5, 0] = 0.1
        self.trafo_h2[5, 0] = 0.02j
        self.line_x[6, :] = [0.5 * line.x_pu, 1.3 * lines[3].x_pu]
        self.trafo_r[6, 0] = 2. * trafo.r_pu + 0.002
        # the rows changing r, x or the ratio: shunts on both branches, so that those matter
        self.rx_rows = (1, 3, 4, 6)
        for row in self.rx_rows:
            self.line_h1[row, 0], self.line_h2[row, 0] = self.SHUNTS[0], self.SHUNTS[1]
            self.trafo_h1[row, 0], self.trafo_h2[row, 0] = self.SHUNTS[2], self.SHUNTS[3]

    def _sweep(self, algo):
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(algo)
        sweep.set_line_parameters(self.line_ids, self.line_r, self.line_x, self.line_h1, self.line_h2)
        sweep.set_trafo_parameters(self.trafo_ids, self.trafo_r, self.trafo_x, self.trafo_h1, self.trafo_h2,
                                   self.trafo_ratio, self.trafo_shift)
        return sweep

    def _grid_of_row(self, row):
        return _grid_holding_row(self.grid, row,
                                 self.line_ids, (self.line_r, self.line_x, self.line_h1, self.line_h2),
                                 self.trafo_ids, (self.trafo_r, self.trafo_x, self.trafo_h1, self.trafo_h2),
                                 self.trafo_ratio, self.trafo_shift)

    def _shunts_only(self, algo):
        """the voltages with the rx rows' shunts and the grid's r, x, ratio"""
        one = lambda v: np.full((1, 1), v, dtype=complex)  # noqa: E731
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(algo)
        sweep.set_line_parameters(np.array([HALF_OPEN_LINE]), _nan((1, 1)), _nan((1, 1)),
                                  one(self.SHUNTS[0]), one(self.SHUNTS[1]))
        sweep.set_trafo_parameters(self.trafo_ids, _nan((1, 1)), _nan((1, 1)), one(self.SHUNTS[2]),
                                   one(self.SHUNTS[3]), _nan((1, 1)), _nan((1, 1)))
        sweep.compute(self.V0, MAX_IT, TOL)
        return 1. * sweep.get_voltages()[0]

    def _check(self, ac):
        algo = AlgorithmType.NR_KLU if ac else AlgorithmType.DC_KLU
        if algo not in ScenarioSweepCPP(self.grid).available_default_algorithms():
            self.skipTest(f"{algo} is not available")
        sweep = self._sweep(algo)
        sweep.compute(self.V0, MAX_IT, TOL)
        self.assertTrue(np.all(sweep.converged_mask()))
        for stats in sweep.get_linear_solver_stats_per_algo():
            self.assertLessEqual(stats.nb_analyze, 1)
        Vs = 1. * sweep.get_voltages()
        P = 1. * sweep.compute_power_flows()
        A = 1. * sweep.compute_flows()
        res = 1. * sweep.compute_branch_results() if ac else None
        for row in range(self.nb_rows):
            grid = self._grid_of_row(row)
            V = grid.ac_pf(self.V0, MAX_IT, TOL) if ac else grid.dc_pf(self.V0, MAX_IT, TOL)
            self.assertEqual(V.shape[0], Vs.shape[1], f"row {row}: the grid's solve diverged")
            on = np.abs(Vs[row]) > 0.
            np.testing.assert_allclose(Vs[row][on], V[on], atol=1e-9, rtol=0, err_msg=f"row {row}")
            p_ref = np.array([el.res_p1_mw for el in grid.get_lines()] + [el.res_p1_mw for el in grid.get_trafos()])
            np.testing.assert_allclose(P[row], p_ref, atol=1e-7, rtol=0, err_msg=f"row {row}")
            if not ac:
                continue
            a_ref = np.array([el.res_a1_ka for el in grid.get_lines()] + [el.res_a1_ka for el in grid.get_trafos()])
            np.testing.assert_allclose(A[row], a_ref, atol=1e-8, rtol=0, err_msg=f"row {row}")
            diff = (sweep.get_Ybus(row) - grid.get_Ybus()).toarray()
            self.assertLess(np.abs(diff).max(), 1e-10, f"row {row}")
            # both ends of every branch (the side-2 results of the half-open trafo included)
            ref = ScenarioSweepCPP(grid)
            ref.change_algorithm(algo)
            ref.set_contingency_lines(np.zeros((1, len(grid.get_lines())), dtype=bool))
            ref.compute(self.V0, MAX_IT, TOL)
            np.testing.assert_allclose(res[row], ref.compute_branch_results()[0], atol=1e-7, rtol=0,
                                       err_msg=f"row {row}")
        return Vs, algo

    def test_ac(self):
        Vs, algo = self._check(ac=True)
        # guard the guard: the r, x and ratio of those rows did change the solve
        V_shunts = self._shunts_only(algo)
        for row in self.rx_rows:
            self.assertGreater(np.abs(Vs[row] - V_shunts).max(), 1e-6, f"row {row}")

    def test_dc(self):
        Vs, _ = self._check(ac=False)
        # nothing of the half-open branches reaches the DC matrix: rows 1 to 5 are the grid's
        for row in range(1, 6):
            np.testing.assert_array_equal(Vs[row], Vs[0], err_msg=f"row {row}")


if __name__ == "__main__":
    unittest.main()
