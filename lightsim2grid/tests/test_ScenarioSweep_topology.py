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
    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(),
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
        """trafo 3 out strands a bus (with generator 4 alone on it, taken out too): the
        row is NOT_SIMULATED by default, solved on the main component in the
        handle_disconnected_grid mode -- as the masks do it"""
        action = self._act({"set_line_status": [(self.n_line + 3, -1)],
                            "set_bus": {"generators_id": [(4, -1)]}})
        ref_V, _, ref_grid = self._reference(0, action)
        live = np.asarray(ref_grid.id_ac_solver_to_me(), dtype=int)
        # the base grid's own solved buses, off a solve (a copy modified since its
        # last powerflow answers with a stale labelling)
        base = copy.deepcopy(self.grid)
        self.assertGreater(base.ac_pf(1.0 * self.Vinit, self.max_it, self.tol).shape[0], 0)
        stranded = sorted(set(np.asarray(base.id_ac_solver_to_me(), dtype=int)) - set(live))
        self.assertTrue(stranded, "this test needs a contingency that strands a bus")

        sweep = ScenarioSweepCPP(self.grid)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(self.load_p[:1])
        sweep.modify_load_q(self.load_q[:1])
        sweep.modify_gen_p(self.gen_p[:1])
        sweep.set_topo_actions(self._topo([action]))
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertFalse(sweep.converged_mask()[0])
        self.assertEqual(sweep.get_violations()[0][0].violation_type, LimitViolationType.NOT_SIMULATED)

        sweep.handle_disconnected_grid = True
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertTrue(sweep.converged_mask()[0])
        V = sweep.get_voltages()[0]
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
        # keep_jacobian and the physical checks are not wired for it
        for name in ("keep_jacobian", "compute_physical_violations"):
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


if __name__ == "__main__":
    unittest.main()
