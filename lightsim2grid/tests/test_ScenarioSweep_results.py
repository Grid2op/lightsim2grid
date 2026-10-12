# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The per-row results of a batch beyond the voltages, against a one-at-a-time solve.

* compute_branch_results: P and Q at both ends of every branch, (nb_rows, nb_branch, 4),
  equal to get_line_res1/2 and get_trafo_res1/2 of the same topology solved alone;
* get_Ybus / get_dcYbus of a row: the pattern and values LSGrid.get_Ybus gives for that
  topology;
* get_row_solve_times / get_row_nb_iter: one entry per row.
"""

import functools
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ScenarioSweepCPP, ContingencyAnalysisCPP,
                                             InjectionSweepCPP)
from lightsim2grid.algorithm import AlgorithmType

MAX_IT = 30
TOL = 1e-10


@functools.lru_cache(maxsize=None)
def _make_net():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case118()
        # a few phase shifters, of both signs: case118 has none of its own
        net.trafo.loc[[0, 3, 6], "shift_degree"] = [5., -8., 12.]
        return net


def _make_grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(_make_net())


def _algos(ac):
    grid = _make_grid()
    available = ScenarioSweepCPP(grid).available_default_algorithms()
    if ac:
        res = [AlgorithmType.NR_SparseLU]
        if AlgorithmType.NR_KLU in available:
            res.append(AlgorithmType.NR_KLU)
    else:
        res = [AlgorithmType.DC_SparseLU]
        if AlgorithmType.DC_KLU in available:
            res.append(AlgorithmType.DC_KLU)
    return res


@functools.lru_cache(maxsize=None)
def _reference(line_off, trafo_off, gen_off, ac):
    """one-at-a-time solve of the grid with these elements out: (branch results
    (nb_branch, 4), Ybus, solved voltages), or None if it does not converge"""
    ref = _make_grid()
    for l_id in line_off:
        ref.deactivate_powerline(int(l_id))
    for t_id in trafo_off:
        ref.deactivate_trafo(int(t_id))
    for g_id in gen_off:
        ref.deactivate_gen(int(g_id))
    V0 = np.ones(ref.total_bus(), dtype=complex)
    V = ref.ac_pf(V0, MAX_IT, TOL) if ac else ref.dc_pf(V0, MAX_IT, TOL)
    if V.shape[0] == 0:
        return None
    l1, l2 = ref.get_line_res1(), ref.get_line_res2()
    t1, t2 = ref.get_trafo_res1(), ref.get_trafo_res2()
    res = np.stack([np.concatenate([l1[0], t1[0]]), np.concatenate([l1[1], t1[1]]),
                    np.concatenate([l2[0], t2[0]]), np.concatenate([l2[1], t2[1]])], axis=1)
    ybus = ref.get_Ybus() if ac else ref.get_dcYbus()
    return res, ybus.tocsr(), np.asarray(V)


def _random_rows(grid, n_row, seed):
    n_line, n_trafo = len(grid.get_lines()), len(grid.get_trafos())
    n_gen = len(grid.get_generators())
    slack_gens = {g.id for g in grid.get_generators() if g.is_slack}
    gen_candidates = [g.id for g in grid.get_generators() if g.id not in slack_gens]
    rng = np.random.default_rng(seed)
    line_mask = np.zeros((n_row, n_line), dtype=bool)
    trafo_mask = np.zeros((n_row, n_trafo), dtype=bool)
    gen_mask = np.zeros((n_row, n_gen), dtype=bool)
    for row in range(n_row):
        kind = rng.integers(0, 4)
        if kind in (0, 3):
            line_mask[row, rng.choice(n_line, size=rng.integers(1, 3), replace=False)] = True
        if kind == 1:
            trafo_mask[row, rng.integers(0, n_trafo)] = True
        if kind in (2, 3):
            gen_mask[row, rng.choice(gen_candidates, size=rng.integers(1, 3), replace=False)] = True
    return line_mask, trafo_mask, gen_mask


def _row_key(line_mask, trafo_mask, gen_mask, row):
    return (tuple(np.flatnonzero(line_mask[row]).tolist()),
            tuple(np.flatnonzero(trafo_mask[row]).tolist()),
            tuple(np.flatnonzero(gen_mask[row]).tolist()))


def _islanding_line(grid):
    """a line whose disconnection leaves a bus with no branch at all"""
    degree = {}
    branches = list(grid.get_lines()) + list(grid.get_trafos())
    for br in branches:
        if not br.connected_global:
            continue
        for b in (br.bus1_id, br.bus2_id):
            degree[b] = degree.get(b, 0) + 1
    for line in grid.get_lines():
        if line.connected_global and (degree[line.bus1_id] == 1 or degree[line.bus2_id] == 1):
            return line.id
    raise RuntimeError("no radial line in this grid")


def _sweep(grid, algo, line_mask, trafo_mask, gen_mask):
    sweep = ScenarioSweepCPP(grid)
    sweep.change_algorithm(algo)
    sweep.set_contingency_lines(line_mask)
    sweep.set_contingency_trafos(trafo_mask)
    sweep.set_contingency_gens(gen_mask)
    sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
    return sweep


class TestBranchResults(unittest.TestCase):
    n_row = 24

    def _check(self, ac):
        grid = _make_grid()
        n_line, n_trafo = len(grid.get_lines()), len(grid.get_trafos())
        line_mask, trafo_mask, gen_mask = _random_rows(grid, self.n_row, seed=1)
        # one row that islands a bus: it does not converge, every branch reads 0
        line_mask[0] = False
        trafo_mask[0] = False
        gen_mask[0] = False
        line_mask[0, _islanding_line(grid)] = True
        for algo in _algos(ac):
            sweep = _sweep(grid, algo, line_mask, trafo_mask, gen_mask)
            res = sweep.compute_branch_results()
            assert res.shape == (self.n_row, n_line + n_trafo, 4)
            assert np.array_equal(sweep.get_branch_results(), res)
            conv = np.asarray(sweep.converged_mask())
            assert not conv[0]
            assert np.all(res[0] == 0.)
            sweep.compute_power_flows()
            p_or = np.asarray(sweep.get_power_flows())
            nb_checked = 0
            for row in range(self.n_row):
                tripped = np.concatenate([np.flatnonzero(line_mask[row]),
                                          n_line + np.flatnonzero(trafo_mask[row])])
                assert np.all(res[row, tripped] == 0.), f"{algo}: row {row}, a tripped branch has a flow"
                if not conv[row]:
                    continue
                ref = _reference(*_row_key(line_mask, trafo_mask, gen_mask, row), ac)
                assert ref is not None, f"{algo}: row {row} converged, the reference did not"
                ref_res = ref[0]
                err = np.abs(res[row] - ref_res).max()
                assert err <= 1e-6, f"{algo}: row {row}, max error {err} MW / MVAr"
                # the origin side is the one compute_power_flows already reports
                assert np.allclose(res[row, :, 0], p_or[row], atol=1e-9, rtol=0.)
                if not ac:
                    assert np.all(res[row, :, 1] == 0.) and np.all(res[row, :, 3] == 0.)
                    assert np.allclose(res[row, :, 2], -res[row, :, 0], atol=1e-9, rtol=0.)
                nb_checked += 1
            assert nb_checked >= self.n_row // 2

    def test_ac(self):
        self._check(ac=True)

    def test_dc(self):
        self._check(ac=False)

    def test_contingency_analysis(self):
        """same accessor on ContingencyAnalysis: rows follow the sorted contingencies"""
        grid = _make_grid()
        n_line = len(grid.get_lines())
        br_ids = [3, 10, n_line, n_line + 6]  # two lines, a trafo, a phase shifter
        for ac in (True, False):
            for algo in _algos(ac):
                ca = ContingencyAnalysisCPP(grid)
                ca.change_algorithm(algo)
                for br_id in br_ids:
                    ca.add_n1(br_id)
                ca.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
                res = ca.compute_branch_results()
                for row, br_id in enumerate(sorted(br_ids)):
                    key = ((br_id,), (), ()) if br_id < n_line else ((), (br_id - n_line,), ())
                    ref_res = _reference(*key, ac)[0]
                    assert np.all(res[row, br_id] == 0.)
                    err = np.abs(res[row] - ref_res).max()
                    assert err <= 1e-6, f"{algo}: branch {br_id}, max error {err}"

    def test_injection_sweep(self):
        """a batch with no contingency: every row is the grid's own topology"""
        grid = _make_grid()
        n_load = len(grid.get_loads())
        load_p = np.array([l.target_p_mw for l in grid.get_loads()])
        rows = np.stack([load_p * f for f in (0.95, 1.0, 1.05)])
        sweep = InjectionSweepCPP(grid)
        sweep.modify_load_p(rows)
        sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        res = sweep.compute_branch_results()
        assert res.shape[0] == 3 and res.shape[2] == 4
        assert n_load > 0
        ref_res = _reference((), (), (), True)[0]
        assert np.abs(res[1] - ref_res).max() <= 1e-6
        sweep.compute_power_flows()
        assert np.allclose(res[:, :, 0], np.asarray(sweep.get_power_flows()), atol=1e-9, rtol=0.)


class TestRowYbus(unittest.TestCase):
    n_row = 12

    def _check(self, ac):
        grid = _make_grid()
        line_mask, trafo_mask, gen_mask = _random_rows(grid, self.n_row, seed=2)
        for algo in _algos(ac)[:1]:
            sweep = _sweep(grid, algo, line_mask, trafo_mask, gen_mask)
            getter = sweep.get_Ybus if ac else sweep.get_dcYbus
            other = sweep.get_dcYbus if ac else sweep.get_Ybus
            with self.assertRaises(RuntimeError):
                other(0)
            with self.assertRaises(RuntimeError):
                getter(self.n_row)
            # results of a compute() whose inputs changed since are refused
            stale = _sweep(grid, algo, line_mask, trafo_mask, gen_mask)
            stale.set_contingency_lines(line_mask[::-1].copy())
            with self.assertRaises(RuntimeError):
                (stale.get_Ybus if ac else stale.get_dcYbus)(0)
            nb_checked = 0
            for row in range(self.n_row):
                ref = _reference(*_row_key(line_mask, trafo_mask, gen_mask, row), ac)
                if ref is None:
                    continue
                ybus_ref = ref[1]
                ybus = getter(row).tocsr()
                assert ybus.shape == ybus_ref.shape == (grid.total_bus(), grid.total_bus())
                ybus.sort_indices()
                ybus_ref.sort_indices()
                assert np.array_equal(ybus.indptr, ybus_ref.indptr), f"row {row}: pattern differs"
                assert np.array_equal(ybus.indices, ybus_ref.indices), f"row {row}: pattern differs"
                scale = np.abs(ybus_ref.data).max()
                err = np.abs(ybus.data - ybus_ref.data).max()
                assert err <= 1e-12 * scale, f"row {row}: value error {err}"
                nb_checked += 1
                # the solver numbering: the same matrix, on the connected buses only
                ybus_solver = getter(row, solver_numbering=True)
                assert ybus_solver.shape[0] <= grid.total_bus()
                assert ybus_solver.nnz == ybus.nnz
            assert nb_checked >= self.n_row // 2

    def test_ac(self):
        self._check(ac=True)

    def test_dc(self):
        self._check(ac=False)

    def test_contingency_analysis(self):
        grid = _make_grid()
        n_line = len(grid.get_lines())
        ca = ContingencyAnalysisCPP(grid)
        ca.add_n1(3)
        ca.add_n1(n_line + 6)  # a phase shifter
        ca.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        for row, key in enumerate([((3,), (), ()), ((), (6,), ())]):
            ybus_ref = _reference(*key, True)[1]
            ybus = ca.get_Ybus(row).tocsr()
            ybus.sort_indices()
            ybus_ref.sort_indices()
            assert np.array_equal(ybus.indices, ybus_ref.indices)
            assert np.abs(ybus.data - ybus_ref.data).max() <= 1e-12 * np.abs(ybus_ref.data).max()


class TestRowSolveStats(unittest.TestCase):
    def test_one_entry_per_row(self):
        grid = _make_grid()
        n_row = 6
        line_mask, trafo_mask, gen_mask = _random_rows(grid, n_row, seed=3)
        line_mask[0] = False
        trafo_mask[0] = False
        gen_mask[0] = False
        line_mask[0, _islanding_line(grid)] = True  # skipped before the solver
        sweep = _sweep(grid, AlgorithmType.NR_SparseLU, line_mask, trafo_mask, gen_mask)
        times = np.asarray(sweep.get_row_solve_times())
        iters = np.asarray(sweep.get_row_nb_iter())
        conv = np.asarray(sweep.converged_mask())
        assert times.shape == (n_row,) and iters.shape == (n_row,)
        assert times[0] == 0. and iters[0] == 0
        assert np.all(times[conv.astype(bool)] > 0.)
        assert np.all(iters[conv.astype(bool)] >= 1)
        # a second compute() starts them again
        sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        assert np.asarray(sweep.get_row_nb_iter()).shape == (n_row,)

    def test_threads(self):
        """worker threads write their own rows: same results as single-threaded"""
        grid = _make_grid()
        line_mask, trafo_mask, gen_mask = _random_rows(grid, 16, seed=5)
        out = []
        for nb_thread in (1, 3):
            sweep = ScenarioSweepCPP(grid)
            sweep.change_algorithm(AlgorithmType.NR_SparseLU)
            sweep.nb_thread = nb_thread
            sweep.set_contingency_lines(line_mask)
            sweep.set_contingency_trafos(trafo_mask)
            sweep.set_contingency_gens(gen_mask)
            sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
            conv = np.asarray(sweep.converged_mask()).astype(bool)
            times = np.asarray(sweep.get_row_solve_times())
            assert np.all(times[conv] > 0.)
            out.append((sweep.compute_branch_results(), np.asarray(sweep.get_row_nb_iter())))
        assert np.array_equal(out[0][0], out[1][0])
        assert np.array_equal(out[0][1], out[1][1])


def _make_multi_slack_grid():
    """case118 with a second generator on the bus of generator 4 (a different reactive
    range, so the two split their bus' Q unevenly) and the slack spread over several
    units, two of them on that bus"""
    import copy
    import pandapower as pp
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = copy.deepcopy(_make_net())
        gen4 = net.gen.iloc[4]
        pp.create_gen(net, bus=gen4.bus, p_mw=120., vm_pu=gen4.vm_pu,
                      min_q_mvar=-50., max_q_mvar=150.)
        grid = init_from_pandapower(net)
    n_gen = len(grid.get_generators())
    extra = n_gen - 1
    slack = np.zeros(n_gen, dtype=bool)
    slack[[4, 5, 10, 11, extra]] = True
    slack[[g.id for g in grid.get_generators() if g.is_slack]] = True
    grid.update_slack_weights(slack)
    return grid, extra


@functools.lru_cache(maxsize=None)
def _reference_gen(gen_p, line_off, gen_off, ac):
    """one-at-a-time solve of the multi-slack grid: get_gen_res as (nb_gen, 2), or None"""
    ref, _ = _make_multi_slack_grid()
    for g_id, p in enumerate(gen_p):
        ref.change_p_gen(g_id, float(p))
    for l_id in line_off:
        ref.deactivate_powerline(int(l_id))
    for g_id in gen_off:
        ref.deactivate_gen(int(g_id))
    V0 = np.ones(ref.total_bus(), dtype=complex)
    V = ref.ac_pf(V0, MAX_IT, TOL) if ac else ref.dc_pf(V0, MAX_IT, TOL)
    if V.shape[0] == 0:
        return None
    p, q = ref.get_gen_res()[:2]
    return np.stack([p, q], axis=1)


class TestGenResults(unittest.TestCase):
    """compute_gen_results: P and Q of every generator per row, against get_gen_res of a
    one-at-a-time solve -- distributed slack over several units, two of them sharing a
    bus, rows that take participants out (one of the pair, both, another one)"""

    def _rows(self, grid, extra):
        n_line, n_gen = len(grid.get_lines()), len(grid.get_generators())
        rows = [([], []), ([7], []), ([], [4]), ([], [extra]), ([], [4, extra]),
                ([20], [10]), ([], [5, 11]), ([33], [])]
        line_mask = np.zeros((len(rows), n_line), dtype=bool)
        gen_mask = np.zeros((len(rows), n_gen), dtype=bool)
        for row, (lines, gens) in enumerate(rows):
            line_mask[row, lines] = True
            gen_mask[row, gens] = True
        rng = np.random.default_rng(4)
        gen_p = np.array(grid.get_gen_target_p()) * rng.uniform(0.9, 1.1, (len(rows), n_gen))
        return line_mask, gen_mask, gen_p

    def _check(self, ac):
        grid, extra = _make_multi_slack_grid()
        line_mask, gen_mask, gen_p = self._rows(grid, extra)
        n_row = line_mask.shape[0]
        for algo in _algos(ac):
            sweep = ScenarioSweepCPP(grid)
            sweep.change_algorithm(algo)
            sweep.compute_gen_results = True
            sweep.modify_gen_p(gen_p)
            sweep.set_contingency_lines(line_mask)
            sweep.set_contingency_gens(gen_mask)
            sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
            res = sweep.get_gen_results()
            assert res.shape == (n_row, len(grid.get_generators()), 2)
            conv = np.asarray(sweep.converged_mask())
            nb_checked = 0
            for row in range(n_row):
                ref = _reference_gen(tuple(gen_p[row]), tuple(np.flatnonzero(line_mask[row]).tolist()),
                                     tuple(np.flatnonzero(gen_mask[row]).tolist()), ac)
                assert conv[row] == (ref is not None), f"{algo}: row {row}, convergence differs"
                if ref is None:
                    assert np.all(res[row] == 0.)
                    continue
                assert np.all(res[row, gen_mask[row]] == 0.), f"{algo}: row {row}, a disconnected generator produces"
                err = np.abs(res[row] - ref).max()
                assert err <= 1e-6, f"{algo}: row {row}, generator error {err} (MW / MVAr)"
                if not ac:
                    assert np.all(res[row, :, 1] == 0.)
                nb_checked += 1
            assert nb_checked >= n_row - 1

    def test_ac(self):
        self._check(ac=True)

    def test_dc(self):
        self._check(ac=False)

    def test_pair_shares_its_bus(self):
        """the guard of the test above: the two generators of one bus do split its slack
        share by weight and its Q by reactive range, so the comparison is not vacuous"""
        grid, extra = _make_multi_slack_grid()
        grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        p, q = grid.get_gen_res()[:2]
        target = np.array(grid.get_gen_target_p())
        assert abs(p[4] - target[4]) > 1. and abs(p[extra] - target[extra]) > 1.
        assert abs(q[4] - q[extra]) > 1.

    def test_off_by_default(self):
        grid, extra = _make_multi_slack_grid()
        line_mask, gen_mask, gen_p = self._rows(grid, extra)
        sweep = ScenarioSweepCPP(grid)
        assert not sweep.compute_gen_results
        sweep.set_contingency_lines(line_mask)
        sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        assert sweep.get_gen_results().size == 0

    def test_contingency_analysis(self):
        grid, extra = _make_multi_slack_grid()
        gen_p = tuple(np.array(grid.get_gen_target_p()))
        for ac in (True, False):
            ca = ContingencyAnalysisCPP(grid)
            ca.change_algorithm(_algos(ac)[0])
            ca.compute_gen_results = True
            for l_id in (3, 20):
                ca.add_n1(l_id)
            ca.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
            res = ca.get_gen_results()
            for row, l_id in enumerate((3, 20)):
                ref = _reference_gen(gen_p, (l_id,), (), ac)
                assert np.abs(res[row] - ref).max() <= 1e-6

    def test_slack_shares(self):
        """LSGrid.get_gen_slack_shares: P - target = share * (total slack absorbed)"""
        grid, extra = _make_multi_slack_grid()
        shares = np.asarray(grid.get_gen_slack_shares())
        # update_slack_weights weighs by target P: a unit flagged with no target takes no part
        assert abs(shares.sum() - 1.) <= 1e-12 and np.count_nonzero(shares) >= 4
        grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        delta = np.asarray(grid.get_gen_res()[0]) - np.array(grid.get_gen_target_p())
        assert np.abs(delta - shares * delta.sum()).max() <= 1e-6


if __name__ == "__main__":
    unittest.main()
