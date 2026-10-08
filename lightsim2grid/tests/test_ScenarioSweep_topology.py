# Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
ScenarioSweep's topological axis: one action per row (``set_topo_actions``, see
batch_algorithm/BaseBatchSweep.hpp and TopoPlan.hpp).

The oracle throughout is a plain, one-off powerflow: for each row, a fresh copy of the
grid with that row's injections and with the action really applied
(``TopoAction.apply_to_gridmodel``), solved by ``ac_pf``. A sweep row must land on the
same voltages and the same flows -- and, as for the generator contingencies, doing so
must NOT cost a symbolic re-factorization per row (``test_single_symbolic_analysis``).

This first version plays disconnections only; what it refuses is pinned here too.
"""

import copy
import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Action import BaseAction
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.scenarioSweep import ScenarioSweep, ScenarioSweepCPP, LimitViolationType, ViolationElementType
from lightsim2grid.lightEnv import TopoAction, ElementType, topo_action_from_grid2op


class _TopoSweepBase(unittest.TestCase):
    """the grid, the chronics and the reference / comparison helpers"""
    env_name = "l2rpn_case14_sandbox"

    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(self.env_name, backend=LightSimBackend(),
                                    param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.Vinit = self.env.backend.V
        self.max_it = self.env.backend.max_it
        self.tol = self.env.backend.tol
        self.n_gen = self.env.n_gen
        self.n_load = self.env.n_load
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())
        self.bus_of_gen = [g.bus_id for g in self.grid.get_generators()]

        # real per-row injections, so that the rows do vary and the injection side of a
        # disconnection (a load taken out) is exercised against the reference
        data = self.env.chronics_handler.real_data.data
        self.nb_steps = 6
        self.gen_p = 1.0 * data.prod_p[:self.nb_steps]
        self.load_p = 1.0 * data.load_p[:self.nb_steps]
        self.load_q = 1.0 * data.load_q[:self.nb_steps]

    def tearDown(self):
        self.env.close()

    # ------------------------------------------------------------------ helpers
    def _act(self, dict_=None):
        return self.env.action_space(dict_ if dict_ is not None else {})

    def _reference(self, row, action, grid=None):
        """One-off powerflow on a copy of the grid with the row's injections and the
        action really applied. Returns (V, flows in kA, the solved grid)."""
        grid = copy.deepcopy(self.grid if grid is None else grid)
        for gen_id in range(self.n_gen):
            grid.change_p_gen(gen_id, float(self.gen_p[row, gen_id]))
        for load_id in range(self.n_load):
            grid.change_p_load(load_id, float(self.load_p[row, load_id]))
            grid.change_q_load(load_id, float(self.load_q[row, load_id]))
        topo = topo_action_from_grid2op(action) if isinstance(action, BaseAction) else action
        topo.check_validity(grid)
        topo.apply_to_gridmodel(grid)
        V = grid.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertGreater(V.shape[0], 0, f"the reference powerflow itself diverged for row {row}")
        flows = np.concatenate([grid.get_line_res1()[3], grid.get_trafo_res1()[3]])
        return V, flows, grid

    @staticmethod
    def _topo(actions):
        """the C++ class takes TopoAction objects; the grid2op wrapper (ScenarioSweep)
        does this conversion for its users"""
        return [topo_action_from_grid2op(act) if isinstance(act, BaseAction) else act for act in actions]

    def _sweep(self, actions, grid=None, nb_thread=1, **kwargs):
        sweep = ScenarioSweepCPP(self.grid if grid is None else grid)
        nb_rows = len(actions)
        sweep.modify_gen_p(self.gen_p[:nb_rows])
        sweep.modify_load_p(self.load_p[:nb_rows])
        sweep.modify_load_q(self.load_q[:nb_rows])
        for name, val in kwargs.items():
            getattr(sweep, name)(val)
        sweep.set_topo_actions(self._topo(actions))
        sweep.nb_thread = nb_thread
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        return sweep

    def _assert_row_matches(self, sweep, row, action, grid=None):
        """the row's voltages on the buses the reference solved, and every flow"""
        ref_V, ref_flows, ref_grid = self._reference(row, action, grid)
        buses = np.asarray(ref_grid.id_ac_solver_to_me(), dtype=int)
        self.assertGreater(buses.size, 0, "no solved bus to compare")
        self.assertTrue(sweep.converged_mask()[row], f"row {row} did not converge")
        got = sweep.get_voltages()[row]
        np.testing.assert_allclose(got[buses], ref_V[buses], rtol=1e-8, atol=1e-8,
                                   err_msg=f"row {row}: voltages differ from the one-off powerflow")
        got_flows = sweep.compute_flows()[row]
        np.testing.assert_allclose(got_flows, ref_flows, rtol=1e-6, atol=1e-8,
                                   err_msg=f"row {row}: flows differ from the one-off powerflow")



class TestScenarioSweepTopology(_TopoSweepBase):
    """stage 1: disconnections through actions"""
    def test_branch_disconnection_matches_reference(self):
        """a branch off by set_line_status or by set_bus -1 on either end, lines and trafos"""
        actions = [
            self._act(),
            self._act({"set_line_status": [(3, -1)]}),
            self._act({"set_bus": {"lines_or_id": [(5, -1)]}}),
            self._act({"set_bus": {"lines_ex_id": [(7, -1)]}}),
            self._act({"set_line_status": [(self.n_line + 0, -1)]}),
            self._act({"set_line_status": [(1, -1), (4, -1)]}),
        ]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)
        self.assertEqual(list(sweep.get_row_disconnected_branches(0)), [])
        self.assertEqual(list(sweep.get_row_disconnected_branches(1)), [3])
        self.assertEqual(list(sweep.get_row_disconnected_branches(4)), [self.n_line])
        self.assertEqual(list(sweep.get_row_disconnected_branches(5)), [1, 4])
        # a disconnected branch carries no flow
        amps = sweep.compute_flows()
        self.assertEqual(amps[1, 3], 0.)
        self.assertEqual(amps[4, self.n_line], 0.)
        self.assertNotEqual(amps[0, 3], 0.)

    def test_gen_disconnection_matches_reference(self):
        """a generator off: its bus stays PV while another controller remains, turns PQ
        when the last one goes (gens 2 and 3 share bus 5 on case14)"""
        shared = [g for g in range(self.n_gen) if self.bus_of_gen.count(self.bus_of_gen[g]) > 1]
        self.assertGreaterEqual(len(shared), 2, "this test needs a bus with two generators")
        g_a, g_b = shared[0], shared[1]
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack and g not in (g_a, g_b)]
        actions = [
            self._act({"set_bus": {"generators_id": [(g_a, -1)]}}),
            self._act({"set_bus": {"generators_id": [(g_a, -1), (g_b, -1)]}}),
            self._act({"set_bus": {"generators_id": [(non_slack[0], -1)]}}),
            self._act(),
        ]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)
        # row 1 really turned the shared bus PQ: its magnitude left the setpoint
        vm_set = self.grid.get_generators()[g_a].target_vm_pu
        self.assertAlmostEqual(abs(sweep.get_voltages()[0][self.bus_of_gen[g_a]]), vm_set, places=6)
        self.assertNotAlmostEqual(abs(sweep.get_voltages()[1][self.bus_of_gen[g_a]]), vm_set, places=4)

    def test_load_disconnection_matches_reference(self):
        actions = [
            self._act({"set_bus": {"loads_id": [(0, -1)]}}),
            self._act({"set_bus": {"loads_id": [(3, -1), (5, -1)]}}),
            # a load, a line and a generator in one row
            self._act({"set_bus": {"loads_id": [(2, -1)], "generators_id": [(1, -1)]},
                       "set_line_status": [(2, -1)]}),
        ]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)

    def test_mixed_with_masks(self):
        """an action and the masks on the same row, different elements: both apply"""
        actions = [self._act({"set_line_status": [(3, -1)]}),
                   self._act({"set_bus": {"generators_id": [(1, -1)]}})]
        line_mask = np.zeros((2, self.n_line), dtype=bool)
        line_mask[0, 5] = True
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)
        gen_mask[1, 2] = True
        sweep = self._sweep(actions, set_contingency_lines=line_mask, set_contingency_gens=gen_mask)
        self.assertEqual(sweep.get_status(), 1)
        both = [self._act({"set_line_status": [(3, -1), (5, -1)]}),
                self._act({"set_bus": {"generators_id": [(1, -1), (2, -1)]}})]
        for row in range(2):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, both[row])
        self.assertEqual(list(sweep.get_row_disconnected_branches(0)), [3, 5])

    def test_mask_and_action_on_the_same_element_refused(self):
        # generator 2 is on in every base grid these tests use (generator 1 is not)
        actions = [self._act({"set_line_status": [(3, -1)]}),
                   self._act({"set_bus": {"generators_id": [(2, -1)]}})]
        line_mask = np.zeros((2, self.n_line), dtype=bool)
        line_mask[0, 3] = True
        with self.assertRaises(RuntimeError) as cm:
            self._sweep(actions, set_contingency_lines=line_mask)
        self.assertIn("row 0", str(cm.exception))
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)
        gen_mask[1, 2] = True
        with self.assertRaises(RuntimeError) as cm:
            self._sweep(actions, set_contingency_gens=gen_mask)
        self.assertIn("row 1", str(cm.exception))

    def test_no_op_reconnection_is_a_plain_row(self):
        sweep = self._sweep([self._act({"set_line_status": [(3, 1)]})])
        self.assertEqual(sweep.get_status(), 1)
        self._assert_row_matches(sweep, 0, self._act())

    def test_half_open_branch_placement_refused(self):
        """a branch with one end open in the base grid, which a row puts on (moves one of
        its ends, or reconnects it): what becomes of the open end is not decided
        (TopoAction.apply_to_gridmodel closes it on the bus it was last on), so the row is
        refused rather than solved on a guess. Taking it out is well defined, and solved."""
        half_open = 5
        grid = copy.deepcopy(self.grid)
        grid.deactivate_powerline_side2(half_open)
        self.assertTrue(grid.get_lines_status()[half_open])
        refused = {
            "an end moved": self._act({"set_bus": {"lines_or_id": [(half_open, 2)], "loads_id": [(0, 2)]}}),
            "reconnected": self._act({"set_line_status": [(half_open, 1)]}),
        }
        for name, action in refused.items():
            with self.subTest(name=name):
                with self.assertRaises(RuntimeError) as cm:
                    self._sweep([self._act(), action], grid=grid)
                self.assertIn("row 1", str(cm.exception))
                self.assertIn(f"powerline {half_open}", str(cm.exception))
        actions = [self._act({"set_line_status": [(half_open, -1)]}),
                   self._act({"set_bus": {"lines_or_id": [(half_open, -1)]}}),
                   self._act({"set_line_status": [(3, -1)]})]
        sweep = self._sweep(actions, grid=grid)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action, grid=grid)

    def test_slack_generator_move_refused(self):
        slack = [g for g in range(self.n_gen) if self.grid.get_generators()[g].is_slack]
        self.assertTrue(slack)
        with self.assertRaises(RuntimeError) as cm:
            self._sweep([self._act(), self._act({"set_bus": {"generators_id": [(slack[0], 2)]}})])
        self.assertIn("slack", str(cm.exception))
        self.assertIn("row 1", str(cm.exception))

    def test_invalid_action_refused_when_set(self):
        cls = type(self.env)
        bus_m2 = self._act({"set_bus": {"loads_id": [(0, 1)]}})
        bus_m2._set_topo_vect[cls.load_pos_topo_vect[0]] = -2
        load_18 = TopoAction()
        load_18.add_element(ElementType.load, 18, -1)
        sweep = ScenarioSweepCPP(self.grid)
        for name, action in {"bus -2": topo_action_from_grid2op(bus_m2), "load 18": load_18}.items():
            with self.subTest(name=name):
                with self.assertRaises(ValueError) as cm:
                    sweep.set_topo_actions([topo_action_from_grid2op(self._act()), action])
                self.assertIn("action 1", str(cm.exception))
        # nothing was registered: the row count is still free
        sweep.modify_load_p(self.load_p[:3])

    def test_dc_refused(self):
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(AlgorithmType.DC_SparseLU)
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_topo_actions(self._topo([self._act(), self._act({"set_line_status": [(3, -1)]})]))
        with self.assertRaises(RuntimeError) as cm:
            sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertIn("DC", str(cm.exception))
        # a do-nothing action is a plain row, DC or not
        sweep.set_topo_actions(self._topo([self._act(), self._act()]))
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_status(), 1)

    def test_algorithm_without_bus_masking_refused(self):
        """Gauss-Seidel cannot mask a bus: the error is about the topological actions, not
        about handle_disconnected_grid, which nobody turned on"""
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(AlgorithmType.GaussSeidel)
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_topo_actions(self._topo([self._act(), self._act({"set_line_status": [(3, -1)]})]))
        self.assertFalse(sweep.handle_disconnected_grid)
        with self.assertRaises(RuntimeError) as cm:
            sweep.compute(1.0 * self.Vinit, 10000, self.tol)
        self.assertIn("set_topo_actions", str(cm.exception))
        self.assertNotIn("handle_disconnected_grid", str(cm.exception))

    def test_do_nothing_rows_without_bus_masking(self):
        """a do-nothing action is a plain row on an algorithm that cannot mask a bus too
        (as on the DC one, see test_dc_refused): nothing is ever masked"""
        sweep = ScenarioSweepCPP(self.grid)
        sweep.change_algorithm(AlgorithmType.GaussSeidel)
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_topo_actions(self._topo([self._act(), self._act()]))
        sweep.compute(1.0 * self.Vinit, 10000, self.tol)
        self.assertEqual(sweep.get_status(), 1)
        plain = ScenarioSweepCPP(self.grid)
        plain.change_algorithm(AlgorithmType.GaussSeidel)
        plain.modify_load_p(self.load_p[:2])
        plain.compute(1.0 * self.Vinit, 10000, self.tol)
        np.testing.assert_array_equal(sweep.get_voltages(), plain.get_voltages())

    def test_do_nothing_rows_are_bit_identical(self):
        actions = [self._act() for _ in range(self.nb_steps)]
        with_actions = self._sweep(actions)
        plain = ScenarioSweepCPP(self.grid)
        plain.modify_gen_p(self.gen_p)
        plain.modify_load_p(self.load_p)
        plain.modify_load_q(self.load_q)
        plain.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertTrue(np.array_equal(with_actions.get_voltages(), plain.get_voltages()))

    def test_single_symbolic_analysis(self):
        """the point of the feature: N rows, ONE symbolic analysis"""
        rng = np.random.default_rng(0)
        non_slack = [g for g in range(self.n_gen) if not self.grid.get_generators()[g].is_slack]
        safe_lines = [0, 1, 2, 3, 4, 5, 6, 7]
        actions = []
        for row in range(self.nb_steps):
            kind = row % 3
            if kind == 0:
                actions.append(self._act({"set_line_status": [(int(rng.choice(safe_lines)), -1)]}))
            elif kind == 1:
                actions.append(self._act({"set_bus": {"generators_id": [(int(rng.choice(non_slack)), -1)]}}))
            else:
                actions.append(self._act({"set_bus": {"loads_id": [(int(rng.integers(self.n_load)), -1)]}}))
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        stats = sweep.get_linear_solver_stats()
        self.assertEqual(stats.nb_analyze, 1,
                         f"expected a single symbolic analysis for the whole sweep, got "
                         f"{stats.nb_analyze} -- a row is changing the sparsity pattern")
        self.assertGreater(stats.nb_refactorize, self.nb_steps)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)

    def test_threads_agree(self):
        actions = [self._act({"set_line_status": [(3, -1)]}),
                   self._act({"set_bus": {"generators_id": [(1, -1)]}}),
                   self._act({"set_bus": {"loads_id": [(0, -1)]}}),
                   self._act(),
                   self._act({"set_line_status": [(self.n_line + 0, -1)]}),
                   self._act({"set_bus": {"loads_id": [(4, -1)], "generators_id": [(2, -1)]}})]
        one = self._sweep(actions, nb_thread=1)
        four = self._sweep(actions, nb_thread=4)
        np.testing.assert_allclose(four.get_voltages(), one.get_voltages(), rtol=1e-10, atol=1e-10)

    def test_row_splitting_the_grid(self):
        """trafo 3 out strands a bus with generator 4 alone on it: the row is NOT_SIMULATED
        by default, solved on the main component in the handle_disconnected_grid mode -- as
        the masks do it. Taking generator 4 out as well leaves that bus with no element at
        all: nothing is stranded then, and the row is solved in either mode. Both answers
        are the grid without the island (a one-off powerflow does not solve one)."""
        split = self._act({"set_line_status": [(self.n_line + 3, -1)]})
        emptied = self._act({"set_line_status": [(self.n_line + 3, -1)],
                             "set_bus": {"generators_id": [(4, -1)]}})
        ref_V, _, ref_grid = self._reference(0, emptied)
        live = np.asarray(ref_grid.id_ac_solver_to_me(), dtype=int)
        # the base grid's own solved buses, off a solve (a copy modified since its
        # last powerflow answers with a stale labelling)
        base = copy.deepcopy(self.grid)
        self.assertGreater(base.ac_pf(1.0 * self.Vinit, self.max_it, self.tol).shape[0], 0)
        stranded = sorted(set(np.asarray(base.id_ac_solver_to_me(), dtype=int)) - set(live))
        self.assertTrue(stranded, "this test needs a contingency that strands a bus")

        sweep = ScenarioSweepCPP(self.grid)
        sweep.compute_limit_violations = True
        # both rows with the injections of the reference
        sweep.modify_load_p(np.repeat(self.load_p[:1], 2, axis=0))
        sweep.modify_load_q(np.repeat(self.load_q[:1], 2, axis=0))
        sweep.modify_gen_p(np.repeat(self.gen_p[:1], 2, axis=0))
        sweep.set_topo_actions(self._topo([split, emptied]))
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertFalse(sweep.converged_mask()[0])
        self.assertEqual(sweep.get_violations()[0][0].violation_type, LimitViolationType.NOT_SIMULATED)
        self.assertTrue(sweep.converged_mask()[1])

        sweep.handle_disconnected_grid = True
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        for row in range(2):
            with self.subTest(row=row):
                self.assertTrue(sweep.converged_mask()[row])
                V = sweep.get_voltages()[row]
                for b in stranded:
                    self.assertEqual(abs(V[b]), 0., f"stranded bus {b} should read 0")
                np.testing.assert_allclose(V[live], ref_V[live], rtol=1e-6, atol=1e-6)


class TestScenarioSweepGenReactivation(TestScenarioSweepTopology):
    """stage 2: a generator disconnected in the base grid, reactivated by a row on the
    bus it was last on -- the bus turns PQ -> PV for that row, at constant sparsity"""
    def setUp(self):
        super().setUp()
        # generator 1 stands alone on its bus (sub 2): off, that bus is PQ in the base grid
        self.gen_off = 1
        self.assertEqual(self.bus_of_gen.count(self.bus_of_gen[self.gen_off]), 1)
        self.grid = copy.deepcopy(self.grid)
        self.grid.deactivate_gen(self.gen_off)
        self.bus_off = self.bus_of_gen[self.gen_off]
        self.vm_set = self.grid.get_generators()[self.gen_off].target_vm_pu

    def _reco(self):
        return self._act({"set_bus": {"generators_id": [(self.gen_off, 1)]}})

    # the stage 1 tests run again on this base grid (inherited); on top of them:
    def test_reactivation_matches_reference(self):
        shared = [g for g in range(self.n_gen) if self.bus_of_gen.count(self.bus_of_gen[g]) > 1]
        actions = [
            self._reco(),
            self._act(),
            # reactivated while another generator goes out (PQ -> PV and PV -> PQ in one row)
            self._act({"set_bus": {"generators_id": [(self.gen_off, 1), (shared[0], -1), (shared[1], -1)]}}),
            # ... with a line out too
            self._act({"set_bus": {"generators_id": [(self.gen_off, 1)]}, "set_line_status": [(3, -1)]}),
            self._reco(),
        ]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)
        Vs = sweep.get_voltages()
        # the bus is held at the set-point where the generator is back, solved elsewhere
        self.assertAlmostEqual(abs(Vs[0][self.bus_off]), self.vm_set, places=8)
        self.assertNotAlmostEqual(abs(Vs[1][self.bus_off]), self.vm_set, places=3)
        self.assertEqual(sweep.get_linear_solver_stats().nb_analyze, 1)

    def test_reactivation_with_per_row_set_point(self):
        gen_v = np.tile([g.target_vm_pu for g in self.grid.get_generators()], (3, 1))
        gen_v[1, self.gen_off] = self.vm_set + 0.02
        gen_v[2, self.gen_off] = self.vm_set - 0.02
        actions = [self._reco(), self._reco(), self._reco()]
        sweep = self._sweep(actions, modify_gen_v=gen_v)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                # the reference: the set-point the row asked for, then the reactivation
                grid = copy.deepcopy(self.grid)
                grid.change_v_gen(self.gen_off, float(gen_v[row, self.gen_off]))
                self._assert_row_matches(sweep, row, action, grid=grid)
                self.assertAlmostEqual(abs(sweep.get_voltages()[row][self.bus_off]),
                                       gen_v[row, self.gen_off], places=8)

    def test_reactivation_refusals(self):
        # a slack generator: gen 5 carries the slack on case14
        slack = [g for g in range(self.n_gen) if self.grid.get_generators()[g].is_slack]
        self.assertTrue(slack)
        grid = copy.deepcopy(self.env.backend._grid)
        # the base grid keeps its slack: a second slack participant is added, then taken out
        grid.add_gen_slackbus(self.gen_off, 0.5)
        grid.deactivate_gen(self.gen_off)
        with self.assertRaises(RuntimeError) as cm:
            self._sweep([self._reco()], grid=grid)
        self.assertIn("slack", str(cm.exception))
        # keep_jacobian is not wired for it (the physical checks are, see
        # TestScenarioSweepTopologyPhysical)
        for name in ("keep_jacobian",):
            with self.subTest(name=name):
                sweep = ScenarioSweepCPP(self.grid)
                setattr(sweep, name, True)
                sweep.modify_load_p(self.load_p[:1])
                sweep.set_topo_actions(self._topo([self._reco()]))
                with self.assertRaises(RuntimeError) as cm:
                    sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
                self.assertIn(name, str(cm.exception))


class TestScenarioSweepTopologyMoves(_TopoSweepBase):
    """stage 3: elements moved between busbars (a bus created or merged), branches
    reconnected -- every row still runs on the one symbolic analysis of the union
    layout, and matches a one-off powerflow on a grid really rewired"""
    def setUp(self):
        super().setUp()
        cls = type(self.env)
        self.assertEqual(cls.n_busbar_per_sub, 2)
        # substation 1 of case14: load 0, generator 0, the origin of lines 2, 3 and 4,
        # the extremity of line 0
        self.assertEqual(cls.load_to_subid[0], 1)
        self.assertEqual(cls.gen_to_subid[0], 1)
        self.assertEqual(list(cls.line_or_to_subid[[2, 3, 4]]), [1, 1, 1])
        self.assertEqual(cls.line_ex_to_subid[0], 1)

    def test_bus_split_and_merge_match_reference(self):
        # substation 3: load 2 and the extremity of line 3 (whose origin is at substation 1)
        at_sub3 = self.env.action_space.get_obj_connect_to(substation_id=3)
        self.assertIn(3, list(at_sub3["lines_ex_id"]))
        sub3_load = int(at_sub3["loads_id"][0])
        actions = [
            # the substation split in two live buses: busbar 2 takes the load and two lines
            self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]}}),
            # ... merged back: a plain row
            self._act(),
            # only a line end moves: busbar 2 is a bus with one branch and no injection
            self._act({"set_bus": {"lines_or_id": [(2, 2)]}}),
            # the generator moves with a line: busbar 2 turns PV, busbar 1 turns PQ
            self._act({"set_bus": {"generators_id": [(0, 2)], "lines_or_id": [(4, 2)]}}),
            # a split, a disconnection elsewhere and the row's injections together
            self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]},
                       "set_line_status": [(7, -1)]}),
            # two substations rewired in one row: substation 3 as well, its load and the
            # extremity of line 3 on busbar 2 -- a chain of two new buses
            self._act({"set_bus": {"loads_id": [(0, 2), (sub3_load, 2)], "lines_or_id": [(2, 2), (3, 2)],
                                   "lines_ex_id": [(3, 2)]}}),
        ]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action)
        self.assertEqual(sweep.get_linear_solver_stats().nb_analyze, 1)
        Vs = sweep.get_voltages()
        n_sub = type(self.env).n_sub
        # busbar 2 of substation 1 is a live bus on the rows that use it, 0 elsewhere
        self.assertNotEqual(abs(Vs[0][1 + n_sub]), 0.)
        self.assertEqual(abs(Vs[1][1 + n_sub]), 0.)
        # the generator holds the bus it moved to
        vm_set = self.grid.get_generators()[0].target_vm_pu
        self.assertAlmostEqual(abs(Vs[3][1 + n_sub]), vm_set, places=8)
        self.assertNotAlmostEqual(abs(Vs[3][1]), vm_set, places=3)

    def test_reconnection_matches_reference(self):
        grid = copy.deepcopy(self.grid)
        grid.deactivate_powerline(3)
        actions = [
            self._act({"set_line_status": [(3, 1)]}),                              # back where it was
            self._act({"set_bus": {"lines_or_id": [(3, 1)], "lines_ex_id": [(3, 1)]}}),
            self._act(),
            # back on a new busbar, with the load: a bus created by the reconnection
            self._act({"set_bus": {"lines_or_id": [(3, 2)], "loads_id": [(0, 2)]}}),
        ]
        sweep = self._sweep(actions, grid=grid)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action, grid=grid)
        self.assertEqual(sweep.get_linear_solver_stats().nb_analyze, 1)
        amps = sweep.compute_flows()
        self.assertNotEqual(amps[0, 3], 0.)
        self.assertEqual(amps[2, 3], 0.)

    def test_clear_forgets_the_placements(self):
        """clear() drops the actions and what they were resolved into: a batch registered
        afterwards without any action is a plain batch, not one still putting back the
        branch an earlier action reconnected"""
        grid = copy.deepcopy(self.grid)
        grid.deactivate_powerline(3)
        sweep = self._sweep([self._act({"set_line_status": [(3, 1)]}), self._act()], grid=grid)
        self.assertEqual(sweep.get_status(), 1)
        sweep.clear()
        plain = ScenarioSweepCPP(grid)
        for computer in (sweep, plain):
            computer.modify_gen_p(self.gen_p[:2])
            computer.modify_load_p(self.load_p[:2])
            computer.modify_load_q(self.load_q[:2])
            computer.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_status(), 1)
        np.testing.assert_allclose(sweep.get_voltages(), plain.get_voltages(), rtol=1e-10, atol=1e-10)
        self.assertEqual(sweep.compute_flows()[0, 3], 0.)

    def test_merge_needs_no_disconnected_grid_mode(self):
        """a row merging back a busbar the base grid splits leaves that busbar with no
        element at all: nothing is stranded, so the row is solved without
        handle_disconnected_grid, as the split itself is"""
        grid = copy.deepcopy(self.grid)
        split = self._topo([self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]}})])[0]
        split.check_validity(grid)
        split.apply_to_gridmodel(grid)
        actions = [self._act({"set_bus": {"loads_id": [(0, 1)], "lines_or_id": [(2, 1), (3, 1)]}}),
                   self._act()]
        sweep = self._sweep(actions, grid=grid)
        self.assertFalse(sweep.handle_disconnected_grid)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action, grid=grid)
        # the busbar emptied by the merge reads 0, the one the base grid uses does not
        n_sub = type(self.env).n_sub
        self.assertEqual(abs(sweep.get_voltages()[0][1 + n_sub]), 0.)
        self.assertNotEqual(abs(sweep.get_voltages()[1][1 + n_sub]), 0.)

    def test_moved_generator_leaves_its_gen_v_group(self):
        """modify_gen_v: a generator the row moves off a bus it shares no longer has to agree
        with the generator left there, nor writes its set-point there (generators 2 and 3
        share a bus on case14; 3 is the one set_vm visits last)"""
        cls = type(self.env)
        stays, moves = 2, 3
        self.assertEqual(self.bus_of_gen[stays], self.bus_of_gen[moves])
        at_sub = self.env.action_space.get_obj_connect_to(substation_id=int(cls.gen_to_subid[moves]))
        line_end = ("lines_or_id", int(at_sub["lines_or_id"][0])) if len(at_sub["lines_or_id"]) \
            else ("lines_ex_id", int(at_sub["lines_ex_id"][0]))
        moved = self._act({"set_bus": {"generators_id": [(moves, 2)], line_end[0]: [(line_end[1], 2)]}})
        actions = [moved, self._act()]
        target_vm = np.array([g.target_vm_pu for g in self.grid.get_generators()])
        gen_v = np.tile(target_vm, (len(actions), 1))
        gen_v[0, moves] += 0.02
        sweep = self._sweep(actions, modify_gen_v=gen_v)
        self.assertEqual(sweep.get_status(), 1)
        # the reference: the move first (a grid refuses two set-points on one bus), then
        # the moved generator's own
        ref_grid = copy.deepcopy(self.grid)
        move = self._topo([moved])[0]
        move.check_validity(ref_grid)
        move.apply_to_gridmodel(ref_grid)
        ref_grid.change_v_gen(moves, float(gen_v[0, moves]))
        self._assert_row_matches(sweep, 0, self._act(), grid=ref_grid)
        self._assert_row_matches(sweep, 1, actions[1])
        # each bus at the magnitude of the generator that stands on it in the row
        n_sub = cls.n_sub
        V = sweep.get_voltages()[0]
        self.assertAlmostEqual(abs(V[self.bus_of_gen[stays]]), target_vm[stays], places=8)
        self.assertAlmostEqual(abs(V[self.bus_of_gen[moves] + n_sub]), gen_v[0, moves], places=8)

    def test_extra_busbar_not_checked_in_the_n_case(self):
        """the extra busbars of the union layout do not exist in the base grid: masked in
        the "n" case, they are not voltage-checked there either -- as in a plain row"""
        from test_ContingencyAnalysis_limit_violations import _set_tight_limits
        grid = copy.deepcopy(self.grid)
        _set_tight_limits(grid)
        extra = 1 + type(self.env).n_sub
        # a limit the extra busbar breaks wherever it is checked
        vn = grid.get_bus_vn_kv()
        vmax = 1. * vn
        vmax[extra] = 0.5 * vn[extra]
        grid.set_bus_voltage_limits(np.full(vn.shape[0], np.nan), vmax)
        actions = [self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]}}),
                   self._act()]
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_topo_actions(self._topo(actions))
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_status(), 1)

        def buses(violations):
            return {v.element_id for v in violations
                    if v.element_type == ViolationElementType.BUS and v.violation_type == LimitViolationType.HIGH_VOLTAGE}
        self.assertIn(extra, buses(sweep.get_violations()[0]))   # live in the row that uses it
        self.assertNotIn(extra, buses(sweep.get_violations()[1]))
        self.assertNotIn(extra, buses(sweep.get_violations_n()))
        self.assertEqual(buses(sweep.get_violations_n()), buses(sweep.get_violations()[1]))

    def test_isolated_bus_is_masked(self):
        """a load alone on busbar 2 is an island of one bus: masked, the row is solved
        without it -- the same as the load disconnected"""
        actions = [self._act({"set_bus": {"loads_id": [(0, 2)]}})]
        sweep = self._sweep(actions)
        self.assertEqual(sweep.get_status(), 1)
        self.assertTrue(sweep.converged_mask()[0])
        n_sub = type(self.env).n_sub
        self.assertEqual(abs(sweep.get_voltages()[0][1 + n_sub]), 0.)
        self._assert_row_matches(sweep, 0, self._act({"set_bus": {"loads_id": [(0, -1)]}}))

    def test_threads_agree_and_plain_rows_match(self):
        actions = [self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]}}),
                   self._act(),
                   self._act({"set_bus": {"generators_id": [(0, 2)], "lines_or_id": [(4, 2)]}}),
                   self._act({"set_line_status": [(7, -1)]}),
                   self._act({"set_bus": {"lines_or_id": [(2, 2)]}}),
                   self._act()]
        one = self._sweep(actions, nb_thread=1)
        four = self._sweep(actions, nb_thread=4)
        np.testing.assert_allclose(four.get_voltages(), one.get_voltages(), rtol=1e-10, atol=1e-10)
        # a plain row of the union layout is the plain row of the plain layout
        plain = ScenarioSweepCPP(self.grid)
        plain.modify_gen_p(self.gen_p)
        plain.modify_load_p(self.load_p)
        plain.modify_load_q(self.load_q)
        plain.compute(1.0 * self.Vinit, self.max_it, self.tol)
        buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)
        np.testing.assert_allclose(one.get_voltages()[1][buses], plain.get_voltages()[1][buses], rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(one.get_voltages()[5][buses], plain.get_voltages()[5][buses], rtol=1e-9, atol=1e-9)

    def test_violations_follow_the_row(self):
        from test_ContingencyAnalysis_limit_violations import _set_tight_limits
        grid = copy.deepcopy(self.grid)
        _set_tight_limits(grid)
        actions = [self._act({"set_bus": {"loads_id": [(0, 2)], "lines_or_id": [(2, 2), (3, 2)]}}),
                   self._act()]
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(self.load_p[:2])
        sweep.set_topo_actions(self._topo(actions))
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_status(), 1)
        # the currents of the moved lines are read with the row's buses: compared with
        # the same limits on the reference grid
        ref_V, ref_flows, ref_grid = self._reference(0, actions[0], grid)
        got = {(v.element_type, v.element_id, v.side): v.value for v in sweep.get_violations()[0]
               if v.violation_type == LimitViolationType.CURRENT}
        for line_id in (2, 3):
            self.assertIn((ViolationElementType.LINE, line_id, 1), got)
            self.assertAlmostEqual(got[(ViolationElementType.LINE, line_id, 1)], ref_flows[line_id], places=6)


class _RedistributeSlackBase(_TopoSweepBase):
    """``redistribute_slack`` with topological actions: the slack pre-pass shares what the
    row really loses, read off the row's own placement -- a generator moved still
    injects, a load the action disconnects (or leaves alone on a busbar) is gone, the
    elements moved off a busbar the row leaves empty are not stranded, and a unit that
    takes a share takes it on the bus the row gives it. The oracle is a one-off grid with
    the action applied and that lost power redistributed."""
    def setUp(self):
        super().setUp()
        self.grid = copy.deepcopy(self.grid)
        # a distributed slack (the environment's default algorithm is single-slack)
        self.grid.change_algorithm(AlgorithmType.NR_KLU)

    def _two_slack_units(self, grid, tight):
        """generator 5 (the slack of case14) and `tight`, which can only move 2 MW: every
        row losing power saturates it"""
        grid.add_gen_slackbus(tight, 1.)
        target = grid.get_generators()[tight].target_p_mw
        min_p = np.full(self.n_gen, -np.inf)
        max_p = np.full(self.n_gen, np.inf)
        min_p[tight] = target - 2.
        max_p[tight] = target + 2.
        grid.set_gen_p_limits(min_p, max_p)

    def _sweep_redistribute(self, grid, actions, redistribute=True):
        sweep = ScenarioSweepCPP(grid)
        sweep.change_algorithm(AlgorithmType.NR_KLU)
        sweep.set_topo_actions(self._topo(actions))
        sweep.redistribute_slack = redistribute
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        return sweep

    def _reference_redistributed(self, grid, action, lost_mw):
        ref = copy.deepcopy(grid)
        topo = self._topo([action])[0]
        topo.check_validity(ref)
        topo.apply_to_gridmodel(ref)
        if lost_mw != 0.:
            report = ref.redistribute_active_power(lost_mw)
            self.assertGreater(report.nb_saturated, 0, "the tight unit should saturate")
        V = ref.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertGreater(V.shape[0], 0, "the one-off reference diverged")
        return V, np.asarray(ref.id_ac_solver_to_me(), dtype=int)

    def _assert_rows(self, grid, sweep, rows):
        for row, (ref_action, lost_mw) in enumerate(rows):
            with self.subTest(row=row):
                self.assertTrue(sweep.converged_mask()[row], f"row {row} did not converge")
                ref_V, buses = self._reference_redistributed(grid, ref_action, lost_mw)
                got = sweep.get_voltages()[row]
                np.testing.assert_allclose(np.abs(got[buses]), np.abs(ref_V[buses]), rtol=0., atol=1e-6)
                np.testing.assert_allclose(np.angle(got[buses]) - np.angle(got[0]),
                                           np.angle(ref_V[buses]) - np.angle(ref_V[0]), rtol=0., atol=1e-6)



class TestScenarioSweepTopologyRedistributeSlack(_RedistributeSlackBase):
    def test_rows_lose_what_their_action_takes_out(self):
        grid = self.grid
        self._two_slack_units(grid, tight=1)
        at_sub4 = self.env.action_space.get_obj_connect_to(substation_id=4)
        p_load = [l.target_p_mw for l in grid.get_loads()]
        p_gen = [g.target_p_mw for g in grid.get_generators()]
        empty_sub4 = {"loads_id": [(int(l), 2) for l in at_sub4["loads_id"]],
                      "lines_or_id": [(int(l), 2) for l in at_sub4["lines_or_id"]],
                      "lines_ex_id": [(int(l), 2) for l in at_sub4["lines_ex_id"]]}
        actions = [
            # generator 0 moved with a line: it still injects, nothing is lost
            self._act({"set_bus": {"generators_id": [(0, 2)], "lines_or_id": [(4, 2)]}}),
            # a load disconnected by the action: its consumption is gone
            self._act({"set_bus": {"loads_id": [(2, -1)]}}),
            # a load alone on busbar 2: an island of one bus, as good as disconnected
            self._act({"set_bus": {"loads_id": [(0, 2)]}}),
            # every element of substation 4 moved to busbar 2: busbar 1 is left empty
            # (masked), nothing stands on it any more and nothing is lost
            self._act({"set_bus": empty_sub4}),
            # a generator disconnected by the action
            self._act({"set_bus": {"generators_id": [(2, -1)]}}),
            self._act(),
        ]
        rows = [(actions[0], 0.),
                (actions[1], -p_load[2]),
                (self._act({"set_bus": {"loads_id": [(0, -1)]}}), -p_load[0]),
                (actions[3], 0.),
                (actions[4], p_gen[2]),
                (actions[5], 0.)]
        sweep = self._sweep_redistribute(grid, actions)
        self.assertEqual(sweep.get_status(), 1)
        self._assert_rows(grid, sweep, rows)
        # not a vacuous check: where power is lost, the redistribution changes the answer
        plain = self._sweep_redistribute(grid, actions, redistribute=False)
        self.assertGreater(np.max(np.abs(plain.get_voltages()[1] - sweep.get_voltages()[1])), 1e-6)

    def test_reactivated_generator_is_a_gain(self):
        grid = self.grid
        self._two_slack_units(grid, tight=0)
        grid.deactivate_gen(1)
        p_gen1 = grid.get_generators()[1].target_p_mw
        actions = [self._act({"set_bus": {"generators_id": [(1, 1)]}}), self._act()]
        sweep = self._sweep_redistribute(grid, actions)
        self.assertEqual(sweep.get_status(), 1)
        self._assert_rows(grid, sweep, [(actions[0], -p_gen1), (actions[1], 0.)])

    def test_moved_participant_takes_its_share_where_it_lands(self):
        # generator 0 only "can participate in the slack", with a tight range: moved with a
        # line in a row that loses a load, it takes its share on busbar 2
        grid = self.grid
        self._two_slack_units(grid, tight=1)
        target = [g.target_p_mw for g in grid.get_generators()]
        min_p = np.array([g.min_p_mw for g in grid.get_generators()])
        max_p = np.array([g.max_p_mw for g in grid.get_generators()])
        min_p[0] = target[0] - 1.
        max_p[0] = target[0] + 1.
        grid.set_gen_p_limits(min_p, max_p)
        grid.set_gen_can_participate_slack([g == 0 for g in range(self.n_gen)],
                                           np.array([1. if g == 0 else 0. for g in range(self.n_gen)]))
        p_load2 = grid.get_loads()[2].target_p_mw
        actions = [self._act({"set_bus": {"generators_id": [(0, 2)], "lines_or_id": [(4, 2)],
                                          "loads_id": [(2, -1)]}})]
        sweep = self._sweep_redistribute(grid, actions)
        self.assertEqual(sweep.get_status(), 1)
        self._assert_rows(grid, sweep, [(actions[0], -p_load2)])


class TestScenarioSweepTopologyRedistributeSlackStorage(_RedistributeSlackBase):
    """the same with storage units (educ_case14_storage: storage 0 on substation 5, with
    the origin of line 7). Its action space has no set_bus: the rows are TopoAction."""
    env_name = "educ_case14_storage"

    @staticmethod
    def _topo_act(*elements):
        act = TopoAction()
        for el_type, el_id, bus in elements:
            act.add_element(el_type, el_id, bus)
        return act

    def test_storage_taken_out_is_lost(self):
        grid = self.grid
        self._two_slack_units(grid, tight=1)
        grid.change_p_storage(0, 10.)  # consuming
        off = self._topo_act((ElementType.storage, 0, -1))
        actions = [
            off,
            # alone on busbar 2: an island of one bus
            self._topo_act((ElementType.storage, 0, 2)),
            # moved with a line: still consuming
            self._topo_act((ElementType.storage, 0, 2), (ElementType.line_or, 7, 2)),
        ]
        sweep = self._sweep_redistribute(grid, actions)
        self.assertEqual(sweep.get_status(), 1)
        self._assert_rows(grid, sweep, [(off, -10.), (off, -10.), (actions[2], 0.)])

    def test_moved_participant_takes_its_share_where_it_lands(self):
        # storage 0 only "can participate in the slack", 1 MW of room each way, moved with a
        # line in a row that loses a load
        grid = self.grid
        self._two_slack_units(grid, tight=1)
        n_storage = len(grid.get_storages())
        grid.set_storage_p_limits(np.array([-1.] + [np.nan] * (n_storage - 1)),
                                  np.array([1.] + [np.nan] * (n_storage - 1)))
        grid.set_storage_can_participate_slack([s == 0 for s in range(n_storage)],
                                               np.array([1. if s == 0 else 0. for s in range(n_storage)]))
        p_load2 = grid.get_loads()[2].target_p_mw
        actions = [self._topo_act((ElementType.storage, 0, 2), (ElementType.line_or, 7, 2),
                                  (ElementType.load, 2, -1))]
        sweep = self._sweep_redistribute(grid, actions)
        self.assertEqual(sweep.get_status(), 1)
        self._assert_rows(grid, sweep, [(actions[0], -p_load2)])


class TestScenarioSweepTopologyStorage(_TopoSweepBase):
    """storage units a row takes out that take part in the voltage control
    (educ_case14_storage, whose action space has no set_bus: the rows are TopoAction)"""
    env_name = "educ_case14_storage"

    @staticmethod
    def _topo_act(*elements):
        act = TopoAction()
        for el_type, el_id, bus in elements:
            act.add_element(el_type, el_id, bus)
        return act

    def _regulating_storages(self, min_q=-50., max_q=50.):
        """storage 0 the only regulator of a load bus (reactive range [min_q, max_q]),
        storage 1 regulating the bus of generators 2 and 3 alongside them, at their
        set-point"""
        grid = copy.deepcopy(self.grid)
        shared_bus = self.bus_of_gen[2]
        self.assertEqual(self.bus_of_gen[3], shared_bus)
        vm_shared = grid.get_generators()[2].target_vm_pu
        load_bus = grid.get_loads()[2].bus_id
        self.assertNotIn(load_bus, self.bus_of_gen)
        grid.init_storages_full(np.array([5., 0.]), np.array([0., 0.]), [True, True],
                                np.array([1.04, vm_shared]), np.array([min_q, -50.]), np.array([max_q, 50.]),
                                np.array([load_bus, shared_bus], dtype=np.int32))
        grid.set_storage_to_subid(np.array([load_bus, shared_bus], dtype=np.int32))
        grid.tell_solver_need_reset()
        return grid, load_bus, shared_bus, vm_shared

    def test_regulating_storage_taken_out(self):
        """a bus whose last regulator the row takes out turns PQ, a storage unit included;
        and a bus a storage unit still holds stays PV when its generators go"""
        grid, load_bus, shared_bus, vm_shared = self._regulating_storages()
        actions = [self._topo_act((ElementType.storage, 0, -1)),
                   self._topo_act((ElementType.gen, 2, -1), (ElementType.gen, 3, -1)),
                   TopoAction()]
        sweep = self._sweep(actions, grid=grid)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                self._assert_row_matches(sweep, row, action, grid=grid)
        Vs = sweep.get_voltages()
        self.assertNotAlmostEqual(abs(Vs[0][load_bus]), 1.04, places=4)
        self.assertAlmostEqual(abs(Vs[2][load_bus]), 1.04, places=8)
        self.assertAlmostEqual(abs(Vs[1][shared_bus]), vm_shared, places=8)

    def test_regulating_storage_taken_out_leaves_the_reactive_check(self):
        """compute_physical_violations: the reactive range of a storage unit the row takes
        out is not its bus's any more -- given one that a bus solved as PQ (about 0 MVAr)
        falls short of, a check still counting it would report that bus"""
        grid, load_bus, _, _ = self._regulating_storages(min_q=10., max_q=20.)
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_physical_violations = True
        sweep.modify_load_p(self.load_p[:1])
        sweep.set_topo_actions([self._topo_act((ElementType.storage, 0, -1))])
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(sweep.get_status(), 1)
        q_buses = [v.element_id for v in sweep.get_physical_violations()[0]
                   if v.element_type == ViolationElementType.BUS and
                   v.violation_type in (LimitViolationType.LOW_Q, LimitViolationType.HIGH_Q)]
        self.assertNotIn(load_bus, q_buses)


class TestScenarioSweepTopologyPhysical(unittest.TestCase):
    """compute_physical_violations follows a row's generator placements: the reactive
    power a bus asks of its machines is checked where the row puts them"""
    def setUp(self):
        import pandapower.networks as pn
        from lightsim2grid.network import init_from_pandapower
        net = pn.case14()
        net.gen["min_q_mvar"] = -5.
        net.gen["max_q_mvar"] = 5.
        if "min_q_mvar" in net.ext_grid:
            net.ext_grid["min_q_mvar"] = -5.
            net.ext_grid["max_q_mvar"] = 5.
        # the loader expects every busbar in the network (n_sub x n_busbar buses, busbar 2
        # of substation k being bus k + n_sub), as grid2op's own backend builds it
        import pandapower as pp
        self.n_sub = len(net.bus)
        for sub_id in range(self.n_sub):
            pp.create_bus(net, vn_kv=net.bus["vn_kv"].iloc[sub_id], index=sub_id + self.n_sub)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(net, n_sub=self.n_sub, n_busbar_per_sub=2)
        self.Vinit = np.full(self.grid.total_bus(), self.grid.get_init_vm_pu() + 0j)
        self.assertGreater(self.grid.ac_pf(1.0 * self.Vinit, 30, 1e-10).shape[0], 0)
        gens = list(self.grid.get_generators())
        # a regulating, non-slack generator alone on its bus, and a line of that bus
        self.gen_id = next(g_id for g_id, g in enumerate(gens)
                           if g.voltage_regulator_on and not g.is_slack
                           and sum(1 for o in gens if o.bus_id == g.bus_id) == 1)
        self.gen_bus = gens[self.gen_id].bus_id
        line_or_bus = np.asarray(self.grid.get_lines().get_bus_id_side_1(), dtype=int)
        self.line_id = int(np.nonzero(line_or_bus == self.gen_bus)[0][0])

    def _reference_bus_q(self, base_grid, action):
        grid = copy.deepcopy(base_grid)
        action.check_validity(grid)
        action.apply_to_gridmodel(grid)
        V = grid.ac_pf(1.0 * self.Vinit, 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        per_bus = {}
        for gen in grid.get_generators():
            if not gen.connected or not gen.voltage_regulator_on:
                continue
            per_bus[gen.bus_id] = per_bus.get(gen.bus_id, 0.) + gen.res_q_mvar
        slack_buses = {gen.bus_id for gen in grid.get_generators() if gen.connected and gen.is_slack}
        return per_bus, slack_buses

    def _check_rows(self, grid, actions):
        sweep = ScenarioSweepCPP(grid)
        sweep.compute_physical_violations = True
        sweep.physical_violation_tol_mva = 0.
        sweep.modify_gen_p(np.tile([g.target_p_mw for g in grid.get_generators()], (len(actions), 1)))
        sweep.set_topo_actions(actions)
        sweep.compute(1.0 * self.Vinit, 30, 1e-10)
        self.assertEqual(sweep.get_status(), 1)
        for row, action in enumerate(actions):
            with self.subTest(row=row):
                expected, slack_buses = self._reference_bus_q(grid, action)
                viols = [v for v in sweep.get_physical_violations()[row]
                         if v.violation_type in (LimitViolationType.LOW_Q, LimitViolationType.HIGH_Q)]
                self.assertGreater(len(viols), 0)
                reported = set()
                for v in viols:
                    self.assertIn(v.element_id, expected, f"bus {v.element_id} holds no regulating generator")
                    self.assertAlmostEqual(v.value, expected[v.element_id], places=4)
                    reported.add(v.element_id)
                # every bus asking more than its machines own is reported (but the slack
                # bus: the reactive check does not cover it, plain row or not)
                for bus, q in expected.items():
                    if bus in slack_buses:
                        continue
                    if abs(q) > 5. + 1e-6:
                        self.assertIn(bus, reported, f"bus {bus} asks {q:.2f} MVAr and is not reported")
        return sweep

    def test_generator_moved(self):
        move = TopoAction()
        move.add_element(ElementType.gen, self.gen_id, 2)
        move.add_element(ElementType.line_or, self.line_id, 2)
        sweep = self._check_rows(self.grid, [move, TopoAction()])
        new_bus = self.gen_bus + self.n_sub
        reported = {v.element_id for v in sweep.get_physical_violations()[0]}
        self.assertIn(new_bus, reported)
        self.assertNotIn(self.gen_bus, reported)

    def test_generator_reactivated(self):
        grid = copy.deepcopy(self.grid)
        grid.deactivate_gen(self.gen_id)
        reco = TopoAction()
        reco.add_element(ElementType.gen, self.gen_id, 1)
        sweep = self._check_rows(grid, [reco, TopoAction()])
        self.assertIn(self.gen_bus, {v.element_id for v in sweep.get_physical_violations()[0]})
        self.assertNotIn(self.gen_bus, {v.element_id for v in sweep.get_physical_violations()[1]})


class TestScenarioSweepTopologyGrid2op(unittest.TestCase):
    """the grid2op wrapper: grid2op actions in, element ids out of run()"""
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), test=True)
        self.env.set_id(0)
        self.obs = self.env.reset()

    def tearDown(self):
        self.env.close()

    def test_unsupported_grid2op_action_refused(self):
        sweep = ScenarioSweep(self.env)
        for name, dict_ in {"change_bus": {"change_bus": {"loads_id": [0]}},
                            "redispatch": {"redispatch": [(0, 1.)]},
                            "change_line_status": {"change_line_status": [0]}}.items():
            with self.subTest(name=name):
                with self.assertRaises(ValueError) as cm:
                    sweep.set_topo_actions([self.env.action_space({}), self.env.action_space(dict_)])
                self.assertIn("action 1", str(cm.exception))
        with self.assertRaises(ValueError):
            sweep.set_topo_actions([1])

    def test_run_reports_the_disconnected_elements(self):
        from test_ContingencyAnalysis_limit_violations import _set_tight_limits
        _set_tight_limits(self.env.backend._grid)
        n_sim = 3
        n_line = len(self.env.backend._grid.get_lines())
        line_mask = np.zeros((n_sim, n_line), dtype=bool)
        line_mask[2, 5] = True
        actions = [self.env.action_space({}),
                   self.env.action_space({"set_line_status": [(2, -1)]}),
                   self.env.action_space({"set_bus": {"lines_ex_id": [(3, -1)], "loads_id": [(0, -1)]}})]
        sweep = ScenarioSweep(self.env)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(np.tile(self.obs.load_p, (n_sim, 1)))
        sweep.set_contingency_lines(line_mask)
        sweep.set_topo_actions(actions)
        res = sweep.run()
        self.assertEqual(len(res.post_contingency_results), n_sim)
        self.assertEqual(res.post_contingency_results[0].element_ids, [])
        self.assertEqual(res.post_contingency_results[1].element_ids, [2])
        self.assertEqual(res.post_contingency_results[2].element_ids, [3, 5])
        self.assertEqual(res.post_contingency_results[2].element_names,
                         [str(self.env.name_line[3]), str(self.env.name_line[5])])
        for row in res.post_contingency_results:
            self.assertTrue(row.converged)

    def test_run_reports_a_diverging_base_case(self):
        """the "n" powerflow diverging is reported by run(), with what each row
        disconnects, topological actions or not"""
        actions = [self.env.action_space({}), self.env.action_space({"set_line_status": [(2, -1)]})]
        sweep = ScenarioSweep(self.env)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(np.tile(self.obs.load_p, (len(actions), 1)))
        sweep.set_topo_actions(actions)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # a tolerance no powerflow reaches: the "n" case diverges
            sweep.compute(tol=1e-30, ignore_errors=True)
        self.assertNotEqual(sweep.computer.get_status(), 1)
        res = sweep.run()
        self.assertFalse(res.pre_contingency_result.converged)
        self.assertEqual([row.element_ids for row in res.post_contingency_results], [[], [2]])
        for row in res.post_contingency_results:
            self.assertFalse(row.converged)


if __name__ == "__main__":
    unittest.main()
