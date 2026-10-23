# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""DC batch rows against a fresh `dc_pf` of the same topology.

Two things are pinned here:

* every row of a DC ScenarioSweep is solved on its OWN matrix, whatever the rows before
  it did. The DC algorithm keeps its reduced matrix incrementally (a row removes its
  branches, solves, adds them back); a row that edits no admittance must still see the
  matrix with every branch in, not the factorization the previous row left behind.
* the DC active flow of a phase-shifting transformer, as the batch reports it (flows,
  and the current-limit check), follows the model `dc_pf` publishes:
  `P_from = (theta_f - theta_t - shift) / (x . tau)`.
"""

import functools
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ScenarioSweepCPP, ContingencyAnalysisCPP,
                                             LimitViolationType, ViolationElementType)
from lightsim2grid.algorithm import AlgorithmType


def _dc_algos():
    res = [AlgorithmType.DC_SparseLU]
    if AlgorithmType.DC_KLU in ScenarioSweepCPP(_make_grid()).available_default_algorithms():
        res.append(AlgorithmType.DC_KLU)
    return res


@functools.lru_cache(maxsize=None)
def _make_net(phase_shift):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        net = pn.case118()
        if phase_shift:
            # a few phase shifters, of both signs: case118 has none of its own
            net.trafo.loc[[0, 3, 6], "shift_degree"] = [5., -8., 12.]
        return net


def _make_grid(phase_shift=False):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(_make_net(phase_shift))


@functools.lru_cache(maxsize=None)
def _reference_cached(line_off, trafo_off, gen_off, phase_shift, max_it, tol):
    """see _reference: one dc_pf per topology, shared by every algorithm and mode tested"""
    ref = _make_grid(phase_shift)
    for l_id in line_off:
        ref.deactivate_powerline(int(l_id))
    for t_id in trafo_off:
        ref.deactivate_trafo(int(t_id))
    for g_id in gen_off:
        ref.deactivate_gen(int(g_id))
    V = ref.dc_pf(np.ones(ref.total_bus(), dtype=complex), max_it, tol)
    assert V.shape[0], "the reference dc_pf diverged"
    p_or = np.concatenate([ref.get_line_res1()[0], ref.get_trafo_res1()[0]])
    return np.asarray(V), p_or


def _reference(line_off, trafo_off, gen_off, phase_shift, max_it=10, tol=1e-8):
    """a fresh dc_pf of the grid with these elements out: (voltages, origin-side P)"""
    return _reference_cached(tuple(int(x) for x in line_off), tuple(int(x) for x in trafo_off),
                             tuple(int(x) for x in gen_off), phase_shift, max_it, tol)


class TestScenarioSweepDCRowIndependence(unittest.TestCase):
    """Every DC row equals a fresh dc_pf of its own topology, whatever the order of
    line / trafo / generator rows."""
    max_it = 10
    tol = 1e-8

    def _sweep(self, grid, algo, line_mask, trafo_mask, gen_mask, handle_disconnected_grid=False):
        sweep = ScenarioSweepCPP(grid)
        sweep.change_algorithm(algo)
        sweep.handle_disconnected_grid = handle_disconnected_grid
        sweep.set_contingency_lines(line_mask)
        sweep.set_contingency_trafos(trafo_mask)
        sweep.set_contingency_gens(gen_mask)
        sweep.compute(np.ones(grid.total_bus(), dtype=complex), self.max_it, self.tol)
        return sweep

    def _check_random_mixes(self, phase_shift):
        grid = _make_grid(phase_shift)
        n_line, n_trafo = len(grid.get_lines()), len(grid.get_trafos())
        n_gen = len(grid.get_generators())
        # non-slack generators that regulate their own bus: the ones a row may remove
        slack_gens = {g.id for g in grid.get_generators() if g.is_slack}
        gen_candidates = [g.id for g in grid.get_generators() if g.id not in slack_gens]
        rng = np.random.default_rng(0)
        n_row = 30
        line_mask = np.zeros((n_row, n_line), dtype=bool)
        trafo_mask = np.zeros((n_row, n_trafo), dtype=bool)
        gen_mask = np.zeros((n_row, n_gen), dtype=bool)
        for row in range(n_row):
            # cycle through "line only", "trafo only", "generator only" and mixes, so that
            # every kind of row follows every other kind somewhere in the batch
            kind = rng.integers(0, 4)
            if kind in (0, 3):
                line_mask[row, rng.choice(n_line, size=rng.integers(1, 3), replace=False)] = True
            if kind == 1:
                trafo_mask[row, rng.integers(0, n_trafo)] = True
            if kind in (2, 3):
                gen_mask[row, rng.choice(gen_candidates, size=rng.integers(1, 3), replace=False)] = True

        for algo, hdg in ((algo, hdg) for algo in _dc_algos() for hdg in (False, True)):
            sweep = self._sweep(grid, algo, line_mask, trafo_mask, gen_mask, hdg)
            Vs = np.asarray(sweep.get_voltages())
            conv = np.asarray(sweep.converged_mask())
            sweep.compute_power_flows()
            flows = np.asarray(sweep.get_power_flows())
            nb_checked = 0
            for row in range(n_row):
                if not conv[row] or np.any(Vs[row] == 0.):
                    # a row that islands the grid (skipped, or solved with its smaller
                    # part masked out): no single dc_pf to compare it with
                    continue
                V_ref, p_ref = _reference(np.flatnonzero(line_mask[row]),
                                               np.flatnonzero(trafo_mask[row]),
                                               np.flatnonzero(gen_mask[row]),
                                               phase_shift)
                nb_checked += 1
                err_deg = np.degrees(np.abs(np.angle(Vs[row]) - np.angle(V_ref))).max()
                assert err_deg <= 1e-10, (f"{algo}, handle_disconnected_grid={hdg}: row {row} angle "
                                          f"error {err_deg} deg")
                err_mw = np.abs(flows[row] - p_ref).max()
                assert err_mw <= 1e-8, (f"{algo}, handle_disconnected_grid={hdg}: row {row} flow "
                                        f"error {err_mw} MW")
            assert nb_checked >= n_row // 2, "too few rows converged to test anything"

    def test_generator_row_after_line_row(self):
        """the minimal case: a row that edits no admittance right after one that did"""
        grid = _make_grid()
        n_line, n_trafo = len(grid.get_lines()), len(grid.get_trafos())
        n_gen = len(grid.get_generators())
        slack_gens = {g.id for g in grid.get_generators() if g.is_slack}
        gen_id = next(g.id for g in grid.get_generators()
                      if g.id not in slack_gens and abs(g.target_p_mw) > 1.)
        line_mask = np.zeros((2, n_line), dtype=bool)
        line_mask[0, 5] = True
        trafo_mask = np.zeros((2, n_trafo), dtype=bool)
        gen_mask = np.zeros((2, n_gen), dtype=bool)
        gen_mask[1, gen_id] = True
        V_ref, _ = _reference([], [], [gen_id], False)
        for algo in _dc_algos():
            for flip in (False, True):
                lm, tm, gm = line_mask, trafo_mask, gen_mask
                if flip:
                    lm, tm, gm = lm[::-1].copy(), tm[::-1].copy(), gm[::-1].copy()
                Vs = np.asarray(self._sweep(grid, algo, lm, tm, gm).get_voltages())
                row = 0 if flip else 1
                err_deg = np.degrees(np.abs(np.angle(Vs[row]) - np.angle(V_ref))).max()
                assert err_deg <= 1e-10, f"{algo}, flip={flip}: angle error {err_deg} deg"

    def test_random_mixes(self):
        self._check_random_mixes(phase_shift=False)

    def test_random_mixes_phase_shifters(self):
        self._check_random_mixes(phase_shift=True)


class TestBatchDCPhaseShifterFlows(unittest.TestCase):
    """the flow through a phase shifter, as the batch reports it, is the one dc_pf
    publishes"""
    max_it = 10
    tol = 1e-8

    def setUp(self):
        self.grid = _make_grid(phase_shift=True)
        self.n_line = len(self.grid.get_lines())
        self.shifters = np.flatnonzero([abs(t.shift_rad) > 0. for t in self.grid.get_trafos()])
        assert self.shifters.size == 3
        ref = _make_grid(phase_shift=True)
        ref.dc_pf(np.ones(ref.total_bus(), dtype=complex), self.max_it, self.tol)
        self.p_or_ref = np.concatenate([ref.get_line_res1()[0], ref.get_trafo_res1()[0]])
        self.p_ex_ref = np.concatenate([ref.get_line_res2()[0], ref.get_trafo_res2()[0]])
        self.a_or_ref = np.concatenate([ref.get_line_res1()[3], ref.get_trafo_res1()[3]])
        self.a_ex_ref = np.concatenate([ref.get_line_res2()[3], ref.get_trafo_res2()[3]])
        # make sure the test actually looks at a non-trivial shift contribution
        p_shift = self.p_or_ref[self.n_line + self.shifters]
        assert np.all(np.abs(p_shift) > 1.)

    def test_power_flows(self):
        for algo in _dc_algos():
            sweep = ScenarioSweepCPP(self.grid)
            sweep.change_algorithm(algo)
            sweep.set_contingency_lines(np.zeros((1, self.n_line), dtype=bool))
            sweep.compute(np.ones(self.grid.total_bus(), dtype=complex), self.max_it, self.tol)
            sweep.compute_power_flows()
            p_or = np.asarray(sweep.get_power_flows())[0]
            err = np.abs(p_or - self.p_or_ref).max()
            assert err <= 1e-8, f"{algo}: max flow error {err} MW"
            sweep.compute_flows()
            a_or = np.asarray(sweep.get_flows())[0]
            err = np.abs(a_or - self.a_or_ref).max()
            assert err <= 1e-6, f"{algo}: max current error {err} A"

    def test_contingency_analysis_trips_a_phase_shifter(self):
        """a DC contingency that disconnects a phase shifter takes its shift injection out
        with it"""
        for algo in _dc_algos():
            for hdg in (False, True):
                ca = ContingencyAnalysisCPP(self.grid)
                ca.change_algorithm(algo)
                ca.handle_disconnected_grid = hdg
                br_ids = [self.n_line + int(t) for t in self.shifters] + [0]
                for br_id in br_ids:
                    ca.add_n1(br_id)
                ca.compute(np.ones(self.grid.total_bus(), dtype=complex), self.max_it, self.tol)
                Vs = np.asarray(ca.get_voltages())
                ca.compute_power_flows()
                flows = np.asarray(ca.get_power_flows())
                # rows follow the sorted contingency ids
                for row, br_id in enumerate(sorted(br_ids)):
                    lines = [br_id] if br_id < self.n_line else []
                    trafos = [br_id - self.n_line] if br_id >= self.n_line else []
                    V_ref, p_ref = _reference(lines, trafos, [], True)
                    err_deg = np.degrees(np.abs(np.angle(Vs[row]) - np.angle(V_ref))).max()
                    assert err_deg <= 1e-10, f"{algo}, hdg={hdg}, branch {br_id}: angle error {err_deg} deg"
                    err_mw = np.abs(flows[row] - p_ref).max()
                    assert err_mw <= 1e-8, f"{algo}, hdg={hdg}, branch {br_id}: flow error {err_mw} MW"

    def test_current_limit_check(self):
        """a limit between the right current and the one the old formula gave is reported
        on the right side of it: one just below the dc_pf current is a violation, one
        just above is not"""
        trafo_id = int(self.shifters[0])
        br_id = self.n_line + trafo_id
        a_or, a_ex = self.a_or_ref[br_id], self.a_ex_ref[br_id]
        assert a_or > 1e-3 and a_ex > 1e-3
        n_trafo = len(self.grid.get_trafos())
        for algo in _dc_algos():
            for scale, expect in ((0.999, True), (1.001, False)):
                grid = _make_grid(phase_shift=True)
                lim1 = np.full(n_trafo, np.nan)
                lim2 = np.full(n_trafo, np.nan)
                lim1[trafo_id] = scale * a_or
                lim2[trafo_id] = scale * a_ex
                grid.set_trafo_current_limit_side1(lim1)
                grid.set_trafo_current_limit_side2(lim2)
                ca = ContingencyAnalysisCPP(grid, True)
                ca.change_algorithm(algo)
                ca.add_n1(0)
                ca.compute(np.ones(grid.total_bus(), dtype=complex), self.max_it, self.tol)
                found = [v for v in ca.get_violations_n()
                         if v.element_type == ViolationElementType.TRAFO and v.element_id == trafo_id
                         and v.violation_type == LimitViolationType.CURRENT]
                assert bool(found) == expect, (f"{algo}, limit x{scale}: expected violation={expect}, "
                                               f"got {found}")
                for v in found:
                    ref = a_or if v.side == 1 else a_ex
                    assert abs(v.value - ref) <= 1e-6 * max(1., ref), (v.value, ref)

    def test_lsgrid_get_violations(self):
        """LSGrid.get_violations(ac=False) shares the check: same answer after a dc_pf"""
        trafo_id = int(self.shifters[0])
        br_id = self.n_line + trafo_id
        a_or, a_ex = self.a_or_ref[br_id], self.a_ex_ref[br_id]
        n_trafo = len(self.grid.get_trafos())
        for scale, expect in ((0.999, True), (1.001, False)):
            grid = _make_grid(phase_shift=True)
            lim1 = np.full(n_trafo, np.nan)
            lim2 = np.full(n_trafo, np.nan)
            lim1[trafo_id] = scale * a_or
            lim2[trafo_id] = scale * a_ex
            grid.set_trafo_current_limit_side1(lim1)
            grid.set_trafo_current_limit_side2(lim2)
            grid.dc_pf(np.ones(grid.total_bus(), dtype=complex), self.max_it, self.tol)
            found = [v for v in grid.get_violations(ac=False)
                     if v.element_type == ViolationElementType.TRAFO and v.element_id == trafo_id
                     and v.violation_type == LimitViolationType.CURRENT]
            assert bool(found) == expect, f"limit x{scale}: expected violation={expect}, got {found}"


if __name__ == "__main__":
    unittest.main()
