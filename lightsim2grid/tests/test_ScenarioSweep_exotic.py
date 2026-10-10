# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Batch rows on a grid carrying the less common elements, against the same
modifications applied to the grid and solved one at a time.

The grid (case14, built here with no pypowsybl) carries at once:

* a remote voltage control group: two generators holding bus 9 remotely, with a
  storage unit and an HVDC VSC station regulating that same bus locally;
* SVCs in every mode: VOLTAGE remote, VOLTAGE local with a slope, REACTIVE_POWER, OFF;
* HVDC lines: VSC-VSC with angle droop (one station in the group, the other at a fixed
  Q), VSC-VSC without droop regulating at both ends, LCC-LCC;
* a storage unit regulating its bus (in the group) and one that does not;
* a line open at side 2 and a transformer open at side 1 (half-open branches).

Each row's generator P and Q (compute_gen_results) are compared too.

Each row changes the injections (generator and load set-points) and disconnects
branches and/or generators. Every converged row must equal a fresh powerflow of the
grid with the same changes: complex voltages, P and Q at both ends of every branch, and
the admittance matrix.
"""

import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ScenarioSweepCPP, ContingencyAnalysisCPP
from lightsim2grid.algorithm import AlgorithmType

MAX_IT = 30
TOL = 1e-10

V_SET = 1.04          # the setpoint every member of the group holds bus 9 at
GROUP_BUS = 9
GROUP_GENS = (2, 3)   # on buses 5 and 7, both holding bus 9 remotely
LOCAL_GENS = (0, 1)   # PV generators holding their own bus
REMOTE_SVC_BUS = 13   # held by the SVC standing on bus 8
HALF_OPEN_LINE = 14   # open at side 2
HALF_OPEN_TRAFO = 4   # open at side 1


def build_grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    for gen_id in GROUP_GENS:
        grid.set_gen_regulated_bus(gen_id, GROUP_BUS)
        grid.change_v_gen(gen_id, V_SET)
    # one SVC per mode, each the only controller of the bus it holds (an SVC cannot
    # share a group): VOLTAGE remote, VOLTAGE local with a slope, REACTIVE_POWER, OFF
    grid.init_svcs([1, 1, 2, 0],
                   np.array([1.03, 1.02, 1.0, 1.0]), np.array([0., 0., 5., 0.]),
                   np.array([0., 0.02, 0., 0.]),
                   np.full(4, -0.5), np.full(4, 0.5),
                   np.array([REMOTE_SVC_BUS, 4, 12, 3], dtype=np.int32),
                   np.array([8, 4, 12, 3], dtype=np.int32))
    # a storage unit in the group, and one that does not regulate
    grid.init_storages_full(np.array([5., -3.]), np.array([0., 1.]), [True, False],
                            np.array([V_SET, 1.]), np.array([-50., -50.]), np.array([50., 50.]),
                            np.array([GROUP_BUS, 12], dtype=np.int32))
    # HVDC: VSC-VSC with droop (side 1 in the group, side 2 at a fixed Q), VSC-VSC
    # without droop regulating at both ends, LCC-LCC
    big = np.full(3, 1e3)
    grid.init_hvdc_lines(np.array([GROUP_BUS, 10, 12], dtype=np.int32), np.array([3, 11, 13], dtype=np.int32),
                         [0, 0, 1], [0, 0, 1],
                         np.ones(3), np.ones(3),
                         [True, True, False], [False, True, False],
                         np.array([V_SET, 1.03, 1.]), np.array([1., 1.03, 1.]),
                         np.zeros(3), np.array([3., 0., 0.]),
                         -big, big, -big, big,
                         np.array([1., 1., 0.9]), np.array([1., 1., 0.9]),
                         [0, 0, 0], np.array([10., 5., 3.]),
                         np.ones(3), np.full(3, 12.),
                         [True, False, False], np.array([5., 0., 0.]), np.array([2., 0., 0.]),
                         np.full(3, 100.), np.full(3, 100.))
    grid.deactivate_powerline_side2(HALF_OPEN_LINE)
    grid.deactivate_trafo_side1(HALF_OPEN_TRAFO)
    grid.tell_solver_need_reset()
    return grid


def _injections(grid, n_row, seed):
    rng = np.random.default_rng(seed)
    gen_p = np.array(grid.get_gen_target_p())
    load_p = np.array(grid.get_load_target_p())
    load_q = np.array([ld.target_q_mvar for ld in grid.get_loads()])
    return (gen_p * rng.uniform(0.9, 1.1, (n_row, gen_p.size)),
            load_p * rng.uniform(0.9, 1.1, (n_row, load_p.size)),
            load_q * rng.uniform(0.9, 1.1, (n_row, load_q.size)))


def _rows(grid, with_group_gens):
    """One row per kind of change, so that every kind follows every other one."""
    n_line, n_trafo, n_gen = len(grid.get_lines()), len(grid.get_trafos()), len(grid.get_generators())
    kinds = [
        ([], [], []),                                   # injections only
        ([3], [], []),                                  # a line
        ([10], [], []),                                 # a line into the group bus
        ([12], [], []),                                 # another one
        ([], [1], []),                                  # a trafo (trafo 0 would island buses 6-7)
        ([HALF_OPEN_LINE], [], []),                     # the half-open line, fully out
        ([], [HALF_OPEN_TRAFO], []),                    # the half-open trafo, fully out
        ([], [], [LOCAL_GENS[0]]),                      # a local PV generator
        ([5], [], [LOCAL_GENS[1]]),                     # a line and a local generator
        ([13], [], []),                                 # a line next to the remote-SVC bus
    ]
    if with_group_gens:
        kinds += [
            ([], [], [GROUP_GENS[0]]),                  # one generator of the group
            ([], [], list(GROUP_GENS)),                 # every generator of the group
            ([6], [], [GROUP_GENS[1]]),                 # a line and one of the group
        ]
    n_row = len(kinds)
    line_mask = np.zeros((n_row, n_line), dtype=bool)
    trafo_mask = np.zeros((n_row, n_trafo), dtype=bool)
    gen_mask = np.zeros((n_row, n_gen), dtype=bool)
    for row, (lines, trafos, gens) in enumerate(kinds):
        line_mask[row, lines] = True
        trafo_mask[row, trafos] = True
        gen_mask[row, gens] = True
    return line_mask, trafo_mask, gen_mask


def _reference(gen_p, load_p, load_q, lines, trafos, gens, ac):
    """the grid with the same changes, solved alone: (V, branch results, Ybus)"""
    grid = build_grid()
    for i, p in enumerate(gen_p):
        grid.change_p_gen(i, float(p))
    for i, p in enumerate(load_p):
        grid.change_p_load(i, float(p))
    for i, q in enumerate(load_q):
        grid.change_q_load(i, float(q))
    for l_id in lines:
        grid.deactivate_powerline(int(l_id))
    for t_id in trafos:
        grid.deactivate_trafo(int(t_id))
    for g_id in gens:
        grid.deactivate_gen(int(g_id))
    V0 = np.ones(grid.total_bus(), dtype=complex)
    V = grid.ac_pf(V0, MAX_IT, TOL) if ac else grid.dc_pf(V0, MAX_IT, TOL)
    if V.shape[0] == 0:
        return None
    l1, l2 = grid.get_line_res1(), grid.get_line_res2()
    t1, t2 = grid.get_trafo_res1(), grid.get_trafo_res2()
    res = np.stack([np.concatenate([l1[0], t1[0]]), np.concatenate([l1[1], t1[1]]),
                    np.concatenate([l2[0], t2[0]]), np.concatenate([l2[1], t2[1]])], axis=1)
    ybus = (grid.get_Ybus() if ac else grid.get_dcYbus()).tocsr()
    ybus.sort_indices()
    gen_p, gen_q = grid.get_gen_res()[:2]
    return np.asarray(V), res, ybus, np.stack([gen_p, gen_q], axis=1)


def _algos(ac):
    available = ScenarioSweepCPP(build_grid()).available_default_algorithms()
    if ac:
        return [a for a in (AlgorithmType.NR_SparseLU, AlgorithmType.NR_KLU) if a in available]
    return [a for a in (AlgorithmType.DC_SparseLU, AlgorithmType.DC_KLU) if a in available]


class TestGridIsWhatItClaims(unittest.TestCase):
    """Guard the guard: if nothing regulated, every comparison below would pass for
    reasons that have nothing to do with the features under test."""

    def test_controls_act(self):
        grid = build_grid()
        V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        assert V.shape[0] > 0
        assert abs(abs(V[GROUP_BUS]) - V_SET) <= 1e-8
        assert abs(abs(V[REMOTE_SVC_BUS]) - 1.03) <= 1e-8
        line = grid.get_lines()[HALF_OPEN_LINE]
        assert line.connected_global and line.connected1 and not line.connected2
        trafo = grid.get_trafos()[HALF_OPEN_TRAFO]
        assert trafo.connected_global and not trafo.connected1 and trafo.connected2
        q_group = np.asarray(grid.get_gen_res()[1])[list(GROUP_GENS)]
        assert np.allclose(q_group, q_group[0])  # the two generators share the bus' reactive power


class TestScenarioSweepExotic(unittest.TestCase):

    def _check(self, ac, handle_disconnected_grid=False):
        grid = build_grid()
        line_mask, trafo_mask, gen_mask = _rows(grid, with_group_gens=not ac)
        n_row = line_mask.shape[0]
        gen_p, load_p, load_q = _injections(grid, n_row, seed=0)
        for algo in _algos(ac):
            sweep = ScenarioSweepCPP(grid)
            sweep.change_algorithm(algo)
            sweep.handle_disconnected_grid = handle_disconnected_grid
            sweep.compute_gen_results = True
            sweep.modify_gen_p(gen_p)
            sweep.modify_load_p(load_p)
            sweep.modify_load_q(load_q)
            sweep.set_contingency_lines(line_mask)
            sweep.set_contingency_trafos(trafo_mask)
            sweep.set_contingency_gens(gen_mask)
            sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
            Vs = np.asarray(sweep.get_voltages())
            conv = np.asarray(sweep.converged_mask())
            res = sweep.compute_branch_results()
            gen_res = sweep.get_gen_results()
            get_ybus = sweep.get_Ybus if ac else sweep.get_dcYbus
            for row in range(n_row):
                msg = f"{algo}, handle_disconnected_grid={handle_disconnected_grid}, row {row}"
                ref = _reference(gen_p[row], load_p[row], load_q[row],
                                 np.flatnonzero(line_mask[row]), np.flatnonzero(trafo_mask[row]),
                                 np.flatnonzero(gen_mask[row]), ac)
                assert conv[row] == (ref is not None), f"{msg}: convergence differs"
                if ref is None:
                    continue
                V_ref, res_ref, ybus_ref, gen_ref = ref
                err_v = np.abs(Vs[row] - V_ref).max()
                assert err_v <= 1e-8, f"{msg}: voltage error {err_v}"
                err_br = np.abs(res[row] - res_ref).max()
                assert err_br <= 1e-6, f"{msg}: branch error {err_br}"
                # generator P and Q, the reactive power of the group's members included
                err_gen = np.abs(gen_res[row] - gen_ref).max()
                assert err_gen <= 1e-6, f"{msg}: generator error {err_gen}"
                ybus = get_ybus(row).tocsr()
                ybus.sort_indices()
                assert np.array_equal(ybus.indices, ybus_ref.indices), f"{msg}: Ybus pattern differs"
                assert np.abs(ybus.data - ybus_ref.data).max() <= 1e-12 * np.abs(ybus_ref.data).max(), msg

    def test_ac(self):
        self._check(ac=True)

    def test_ac_handle_disconnected_grid(self):
        self._check(ac=True, handle_disconnected_grid=True)

    def test_dc(self):
        """DC also takes the generators of the group out (it has no PV / PQ switch)"""
        self._check(ac=False)

    def test_dc_handle_disconnected_grid(self):
        self._check(ac=False, handle_disconnected_grid=True)

    def test_ac_group_generator_contingency_is_refused(self):
        """Not supported yet in AC: disconnecting a generator of a remote control group
        needs Jacobian slots the batch does not reserve. It must say so, for one
        generator of the group as for all of them, rather than solve something else."""
        grid = build_grid()
        n_line, n_gen = len(grid.get_lines()), len(grid.get_generators())
        for gens in ([GROUP_GENS[0]], list(GROUP_GENS)):
            gen_mask = np.zeros((1, n_gen), dtype=bool)
            gen_mask[0, gens] = True
            sweep = ScenarioSweepCPP(grid)
            sweep.set_contingency_lines(np.zeros((1, n_line), dtype=bool))
            sweep.set_contingency_gens(gen_mask)
            with self.assertRaises(RuntimeError) as ctx:
                sweep.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
            assert "remote bus" in str(ctx.exception)


class TestContingencyAnalysisExotic(unittest.TestCase):
    def test_rows_match(self):
        grid = build_grid()
        n_line = len(grid.get_lines())
        br_ids = [3, 10, 12, 13, HALF_OPEN_LINE, n_line + 1, n_line + HALF_OPEN_TRAFO]
        gen_p = np.array(grid.get_gen_target_p())
        load_p = np.array(grid.get_load_target_p())
        load_q = np.array([ld.target_q_mvar for ld in grid.get_loads()])
        for ac in (True, False):
            for algo in _algos(ac):
                ca = ContingencyAnalysisCPP(grid)
                ca.change_algorithm(algo)
                for br_id in br_ids:
                    ca.add_n1(br_id)
                ca.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
                Vs = np.asarray(ca.get_voltages())
                res = ca.compute_branch_results()
                for row, br_id in enumerate(sorted(br_ids)):
                    lines = [br_id] if br_id < n_line else []
                    trafos = [br_id - n_line] if br_id >= n_line else []
                    ref = _reference(gen_p, load_p, load_q, lines, trafos, [], ac)
                    assert ref is not None, f"branch {br_id}: the reference did not converge"
                    assert np.abs(Vs[row] - ref[0]).max() <= 1e-8, f"{algo}, branch {br_id}"
                    assert np.abs(res[row] - ref[1]).max() <= 1e-6, f"{algo}, branch {br_id}"


if __name__ == "__main__":
    unittest.main()
