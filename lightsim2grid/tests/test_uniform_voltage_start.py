# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Issue #284: ``LSGrid.set_uniform_voltage_start``.

From a flat start, the Newton-Raphson diverged when a PV bus is tied to another bus by a
very low impedance branch (with some resistance): the PV bus starts at its target, its
neighbour at the flat start, and the large flow this puts on the branch makes the first
step throw the rest of the grid too far. OpenLoadFlow, from its uniform start, converges.

With the option, the first step is taken from the magnitudes passed to ``ac_pf`` and
moves the held magnitudes to their targets within that step. The solution is the same,
only the path to it changes.
"""

import copy
import os
import pickle
import tempfile
import unittest

import numpy as np

from lightsim2grid.network import LSGrid
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP, TimeSeriesCPP

MAX_IT = 20
TOL = 1e-8

# MATPOWER case14 reduced to the elements needed to show the divergence (per unit on
# 100 MVA): (from, to, r, x, b), plus a bus 15 tied to bus 2 (a PV bus) by a branch of
# impedance of the order of 1e-5 pu
LINES = [(1, 5, 0.05403, 0.22304, 0.0492), (2, 5, 0.05695, 0.17388, 0.0346), (6, 13, 0.06615, 0.13027, 0.),
         (9, 14, 0.12711, 0.27038, 0.), (13, 14, 0.17093, 0.34802, 0.), (5, 6, 0., 0.25202, 0.)]
STIFF = (2, 15, 2.7e-6, 1.1e-5, 0.)
GENS = [(1, 232.4, 1.06), (2, 40., 1.045)]  # bus, p_mw, vm_pu (bus 1 is the slack)
SHUNTS = [(9, -19.)]  # bus, q_mvar
# a second line 1-5, so that a contingency can trip one without splitting the grid
PARALLEL = (1, 5, 0.05403, 0.22304, 0.0492)


def _grid(lines, algo=AlgorithmType.NR_KLU, uniform=False):
    buses = sorted({b for e in lines for b in e[:2]})
    idx = {b: i for i, b in enumerate(buses)}
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(len(buses), 1, np.full(len(buses), 135.), 0, 0)
    grid.init_powerlines(np.array([e[2] for e in lines]), np.array([e[3] for e in lines]),
                         np.array([1j * e[4] for e in lines]),
                         np.array([idx[e[0]] for e in lines], dtype=np.int32),
                         np.array([idx[e[1]] for e in lines], dtype=np.int32))
    grid.init_generators(np.array([g[1] for g in GENS]), np.array([g[2] for g in GENS]),
                         np.full(len(GENS), -1e30), np.full(len(GENS), 1e30),
                         np.array([idx[g[0]] for g in GENS], dtype=np.int32))
    grid.init_shunt(np.zeros(len(SHUNTS)), np.array([s[1] for s in SHUNTS]),
                    np.array([idx[s[0]] for s in SHUNTS], dtype=np.int32))
    grid.add_gen_slackbus(0, 1.)
    grid.change_algorithm(algo)
    grid.set_uniform_voltage_start(uniform)
    return grid, buses


def _flat(buses):
    return np.ones(len(buses), dtype=np.complex128)


class TestUniformVoltageStart(unittest.TestCase):
    def test_default_and_accessors(self):
        grid, _ = _grid(LINES)
        self.assertFalse(grid.get_uniform_voltage_start())
        grid.set_uniform_voltage_start(True)
        self.assertTrue(grid.get_uniform_voltage_start())
        grid.set_uniform_voltage_start(False)
        self.assertFalse(grid.get_uniform_voltage_start())

    def test_kept_by_copy_pickle_and_binary(self):
        grid, _ = _grid(LINES, uniform=True)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            from_file = LSGrid.load_binary(path)
        for other in (grid.copy(), copy.deepcopy(grid), pickle.loads(pickle.dumps(grid)), from_file):
            self.assertTrue(other.get_uniform_voltage_start())
        self.assertFalse(pickle.loads(pickle.dumps(_grid(LINES)[0])).get_uniform_voltage_start())

    def test_state_saved_before_the_option(self):
        # a state written by 1.3.0 has one setting fewer: read back with the option off
        grid, _ = _grid(LINES, uniform=True)
        state = list(grid.__getstate__())
        inner = list(state[3])
        self.assertEqual(len(inner[29]), 6)
        inner[29] = inner[29][:5]
        state[3] = tuple(inner)
        other = LSGrid()
        other.__setstate__(tuple(state))
        self.assertFalse(other.get_uniform_voltage_start())

    def test_low_impedance_branch_next_to_a_pv_bus(self):
        # issue #284: fails from a flat start without the option, converges with it, to the
        # solution reached when bus 15 is started at bus 2's target
        for algo in (AlgorithmType.NR_KLU, AlgorithmType.NR_SparseLU, AlgorithmType.NRSing_KLU,
                     "NROuter_SparseLU"):
            with self.subTest(algo=algo):
                grid, buses = _grid(LINES + [STIFF], algo)
                V = grid.ac_pf(_flat(buses), MAX_IT, TOL)
                self.assertEqual(V.shape[0], 0, "the flat start is expected to diverge without the option")

                grid, buses = _grid(LINES + [STIFF], algo)
                v0 = _flat(buses)
                v0[buses.index(15)] = GENS[1][2]
                V_ref = grid.ac_pf(v0, MAX_IT, TOL)
                self.assertGreater(V_ref.shape[0], 0)

                grid, buses = _grid(LINES + [STIFF], algo, uniform=True)
                V = grid.ac_pf(_flat(buses), MAX_IT, TOL)
                self.assertGreater(V.shape[0], 0, f"diverged with the option: {grid.get_algo().get_error()}")
                np.testing.assert_allclose(V, V_ref, atol=1e-8, rtol=0.)
                for bus, _, vm in GENS:
                    self.assertAlmostEqual(abs(V[buses.index(bus)]), vm, places=10)

    def test_same_solution_on_an_ordinary_grid(self):
        # without the stiff branch both starts converge, to the same voltages
        for algo in (AlgorithmType.NR_KLU, "NROuter_SparseLU"):
            with self.subTest(algo=algo):
                grid, buses = _grid(LINES, algo)
                V_ref = grid.ac_pf(_flat(buses), MAX_IT, TOL)
                grid, buses = _grid(LINES, algo, uniform=True)
                V = grid.ac_pf(_flat(buses), MAX_IT, TOL)
                self.assertGreater(V_ref.shape[0], 0)
                self.assertGreater(V.shape[0], 0)
                np.testing.assert_allclose(V, V_ref, atol=1e-8, rtol=0.)

    def test_warm_start_is_unchanged(self):
        # a converged voltage already holds the targets: the option changes nothing
        grid, buses = _grid(LINES)
        V0 = grid.ac_pf(_flat(buses), MAX_IT, TOL)
        res = []
        for uniform in (False, True):
            grid, _ = _grid(LINES, uniform=uniform)
            res.append(grid.ac_pf(V0, MAX_IT, TOL))
            self.assertEqual(grid.get_algo().get_nb_iter(), 0)
        np.testing.assert_array_equal(res[0], res[1])

    def test_other_algorithms_ignore_it(self):
        # only the Newton-Raphson takes the lifted step: the others solve exactly as before
        for algo in (AlgorithmType.FDPF_XB_KLU, AlgorithmType.GaussSeidel, "dc_pf"):
            with self.subTest(algo=algo):
                res = []
                for uniform in (False, True):
                    if algo == "dc_pf":
                        grid, buses = _grid(LINES, uniform=uniform)
                        res.append(grid.dc_pf(_flat(buses), MAX_IT, TOL))
                    else:
                        grid, buses = _grid(LINES, algo, uniform=uniform)
                        res.append(grid.ac_pf(_flat(buses), 1000, TOL))
                self.assertGreater(res[0].shape[0], 0)
                np.testing.assert_array_equal(res[0], res[1])


class TestUniformStartBatch(unittest.TestCase):
    """The batches inherit the option, and their "n" solve and every row starting from the
    voltage passed take the uniform start; a row starting from a solved state does not."""

    def test_inherited_and_own_property(self):
        grid, _ = _grid(LINES)
        self.assertFalse(ContingencyAnalysisCPP(grid).uniform_voltage_start)
        self.assertFalse(TimeSeriesCPP(grid).uniform_voltage_start)
        grid.set_uniform_voltage_start(True)
        self.assertTrue(ContingencyAnalysisCPP(grid).uniform_voltage_start)
        ts = TimeSeriesCPP(grid)
        self.assertTrue(ts.uniform_voltage_start)
        ts.uniform_voltage_start = False
        self.assertFalse(ts.uniform_voltage_start)
        self.assertTrue(grid.get_uniform_voltage_start())

    def test_contingency_analysis(self):
        # issue #284 in a batch: the "n" solve diverges without the option (the analysis
        # then refuses to run), and with it every row matches a single solve of its own
        lines = LINES + [PARALLEL, STIFF]
        trips = [len(LINES), 0]  # either of the two parallel lines
        grid, buses = _grid(lines)
        ca = ContingencyAnalysisCPP(grid)
        for line_id in trips:
            ca.add_n1(line_id)
        with self.assertRaises(RuntimeError):
            ca.compute(_flat(buses), MAX_IT, TOL)

        for init_from_n in (False, True):
            with self.subTest(init_from_n=init_from_n):
                grid, buses = _grid(lines, uniform=True)
                ca = ContingencyAnalysisCPP(grid)
                ca.init_from_n_powerflow = init_from_n
                for line_id in trips:
                    ca.add_n1(line_id)
                ca.compute(_flat(buses), MAX_IT, TOL)
                Vs = np.array(ca.get_voltages())
                for row, line_id in enumerate(trips):
                    single, _ = _grid(lines, uniform=True)
                    single.deactivate_powerline(line_id)
                    V = single.ac_pf(_flat(buses), MAX_IT, TOL)
                    self.assertGreater(V.shape[0], 0)
                    np.testing.assert_allclose(Vs[row], V, atol=1e-8, rtol=0.)

    def test_time_series(self):
        # the first row takes the uniform start, the next ones start from the previous
        # solution, exactly as without the option
        scales = np.array([1.0, 1.05, 0.95])
        gen_mat = np.outer(scales, [g[1] for g in GENS])
        empty = np.zeros((len(scales), 0))

        grid, buses = _grid(LINES + [STIFF], uniform=True)
        ts = TimeSeriesCPP(grid)
        self.assertEqual(ts.compute_Vs(gen_mat, empty, empty, empty, _flat(buses), MAX_IT, TOL), 1)
        Vs = np.array(ts.get_voltages())
        for row, scale in enumerate(scales):
            single, _ = _grid(LINES + [STIFF], uniform=True)
            for gen_id, gen in enumerate(GENS):
                single.change_p_gen(gen_id, scale * gen[1])
            V = single.ac_pf(_flat(buses), MAX_IT, TOL)
            self.assertGreater(V.shape[0], 0)
            np.testing.assert_allclose(Vs[row], V, atol=1e-8, rtol=0.)

        # on a grid where both starts converge, a chained row is solved the same way
        nb_iter = []
        for uniform in (False, True):
            grid, buses = _grid(LINES, uniform=uniform)
            ts = TimeSeriesCPP(grid)
            self.assertEqual(ts.compute_Vs(gen_mat, empty, empty, empty, _flat(buses), MAX_IT, TOL), 1)
            nb_iter.append(list(ts.get_row_nb_iter()))
        self.assertEqual(nb_iter[0][1:], nb_iter[1][1:])

    def test_changing_it_drops_the_base_case(self):
        grid, buses = _grid(LINES + [PARALLEL])
        ca = ContingencyAnalysisCPP(grid)
        ca.add_n1(len(LINES))
        ca.compute(_flat(buses), MAX_IT, TOL)
        ca.compute(_flat(buses), MAX_IT, TOL)
        self.assertTrue(ca.base_case_was_reused())
        ca.uniform_voltage_start = True
        ca.compute(_flat(buses), MAX_IT, TOL)
        self.assertFalse(ca.base_case_was_reused())
        ca.uniform_voltage_start = True
        ca.compute(_flat(buses), MAX_IT, TOL)
        self.assertTrue(ca.base_case_was_reused())


class TestUniformStartVoltageControlGroup(unittest.TestCase):
    """A bus regulated remotely, or by an SVC, has a Vm unknown and a voltage row: with
    the option it starts from the magnitude passed (like set_keep_vinit_at_group_controlled_buses)
    and the solution does not change."""

    def test_same_solution(self):
        from test_voltage_control_batch import _make, KINDS
        for kind in KINDS:
            with self.subTest(kind=kind):
                net, model = _make(kind)
                vinit = np.full(net.bus.shape[0], 1.0, dtype=np.complex128)
                V_ref = model.ac_pf(vinit, MAX_IT, TOL)
                self.assertGreater(V_ref.shape[0], 0)
                net, model = _make(kind)
                model.set_uniform_voltage_start(True)
                V = model.ac_pf(vinit, MAX_IT, TOL)
                self.assertGreater(V.shape[0], 0)
                np.testing.assert_allclose(V, V_ref, atol=1e-7, rtol=0.)


if __name__ == "__main__":
    unittest.main()
