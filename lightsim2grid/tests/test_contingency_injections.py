# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Contingencies tripping shunts, loads and generators, not only branches:
``ContingencyAnalysisCPP.add_contingency`` / ``add_contingency_shunt`` / ``_load`` /
``_gen``, the ``ContingencyAnalysis`` wrapper's counterparts, and ``ScenarioSweep``'s
``set_contingency_shunts`` / ``set_contingency_loads``.

The oracle is a one-off powerflow on a copy of the grid with the elements really
deactivated.
"""

import copy
import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend, ContingencyAnalysis
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP
from lightsim2grid.scenarioSweep import ScenarioSweep


class _Base(unittest.TestCase):
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
        self.n_shunt = len(self.grid.get_shunts())
        self.assertGreater(self.n_shunt, 0)
        # generator 0 does not share its bus: taking it out turns that bus PQ
        self.gen_alone = 0

    def tearDown(self):
        self.env.close()

    def _reference_V(self, branches=(), shunts=(), loads=(), gens=()):
        grid = copy.deepcopy(self.grid)
        n_line = len(grid.get_lines())
        for br in branches:
            if br < n_line:
                grid.deactivate_powerline(int(br))
            else:
                grid.deactivate_trafo(int(br) - n_line)
        for el in shunts:
            grid.deactivate_shunt(int(el))
        for el in loads:
            grid.deactivate_load(int(el))
        for el in gens:
            grid.deactivate_gen(int(el))
        grid.tell_solver_need_reset()
        V = grid.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertGreater(V.shape[0], 0)
        # a bus the solve did not see (an unused busbar) reads 0, as the batches report it
        solved = np.zeros(V.shape[0], dtype=bool)
        solved[np.asarray(grid.id_ac_solver_to_me())] = True
        V[~solved] = 0.
        return V


class TestContingencyAnalysisCPP(_Base):
    TRIPS = [
        ((), (0,), (), ()),            # a shunt
        ((), (), (3,), ()),            # a load
        ((), (), (), (0,)),            # the last generator regulating its bus: PV -> PQ
        ((5,), (0,), (2, 4), ()),      # a line, the shunt and two loads together
    ]

    def test_matches_a_one_off_powerflow(self):
        ca = ContingencyAnalysisCPP(self.grid)
        for br, sh, ld, gn in self.TRIPS:
            ca.add_contingency(list(br), list(sh), list(ld), list(gn))
        ca.compute(1.0 * self.Vinit, self.max_it, self.tol)
        V = np.asarray(ca.get_voltages())
        conts = ca.my_contingencies()
        self.assertEqual(len(conts), len(self.TRIPS))
        for row, (br, sh, ld, gn, br1, br2) in enumerate(conts):
            self.assertEqual((br1, br2), ([], []))
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference_V(br, sh, ld, gn), atol=1e-8, rtol=0.)

    def test_shortcuts_take_one_id_or_several(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency_shunt(0)
        ca.add_contingency_load([1, 2])
        ca.add_contingency_gen(self.gen_alone)
        ca.add_contingency_load(1)
        conts = ca.my_contingencies()
        self.assertEqual(len(conts), 4)
        self.assertIn(([], [0], [], [], [], []), [tuple(c) for c in conts])
        self.assertIn(([], [], [1, 2], [], [], []), [tuple(c) for c in conts])
        self.assertIn(([], [], [], [self.gen_alone], [], []), [tuple(c) for c in conts])
        # my_defaults lists the branches only
        self.assertEqual(ca.my_defaults(), [[], [], [], []])
        self.assertTrue(ca.remove_contingency(load_ids=[1, 2]))
        self.assertEqual(len(ca.my_contingencies()), 3)
        with self.assertRaises(RuntimeError):
            ca.add_contingency_shunt(self.n_shunt)
        with self.assertRaises(RuntimeError):
            ca.add_contingency()


class TestOneSideTrips(_Base):
    """``branch_ids_side1`` / ``branch_ids_side2``: a branch opened at one end only keeps its
    other end, with its charging, as LSGrid.deactivate_powerline_side1 (and co) leave it"""

    def _reference_one_side(self, branch, side):
        grid = copy.deepcopy(self.grid)
        grid.set_synch_status_both_side(False)
        n_line = len(grid.get_lines())
        if branch < n_line:
            getattr(grid, f"deactivate_powerline_side{side}")(int(branch))
        else:
            getattr(grid, f"deactivate_trafo_side{side}")(int(branch) - n_line)
        grid.tell_solver_need_reset()
        V = grid.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertGreater(V.shape[0], 0)
        solved = np.zeros(V.shape[0], dtype=bool)
        solved[np.asarray(grid.id_ac_solver_to_me())] = True
        V[~solved] = 0.
        return V

    def test_matches_a_half_open_branch(self):
        n_line = len(self.grid.get_lines())
        # lines 1 and 3 have charging (opened at side 1, it shows); line 5 at side 2; the
        # transformers have no shunt admittance, half-open or tripped is the same for them
        trips = [(1, 1), (3, 1), (5, 2), (n_line, 1), (n_line + 1, 2)]
        ca = ContingencyAnalysisCPP(self.grid)
        for br, side in trips:
            ca.add_contingency([], [], [], [], [br] if side == 1 else [], [br] if side == 2 else [])
        ca.compute(1.0 * self.Vinit, self.max_it, self.tol)
        V = np.asarray(ca.get_voltages())
        row_of = {(tuple(c[4]), tuple(c[5])): row for row, c in enumerate(ca.my_contingencies())}
        for br, side in trips:
            row = row_of[((br,), ()) if side == 1 else ((), (br,))]
            self.assertTrue(ca.converged_mask()[row])
            V_ref = self._reference_one_side(br, side)
            np.testing.assert_allclose(V[row], V_ref, atol=1e-8, rtol=0.)
            if br in (1, 3):
                # not the same as tripping it at both ends
                self.assertGreater(np.abs(V_ref - self._reference_V(branches=(br,))).max(), 1e-4)

    def test_both_ends_is_a_trip(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(branch_ids_side1=[3, 5], branch_ids_side2=[3])
        ca.add_contingency(branch_ids=[4], branch_ids_side2=[4])
        conts = sorted(tuple(c) for c in ca.my_contingencies())
        self.assertEqual(conts, [([3], [], [], [], [5], []), ([4], [], [], [], [], [])])
        self.assertTrue(ca.remove_contingency(branch_ids=[3], branch_ids_side1=[5]))
        self.assertEqual(len(ca.my_contingencies()), 1)
        with self.assertRaises(RuntimeError):
            ca.add_contingency(branch_ids_side1=[len(self.grid.get_lines()) + len(self.grid.get_trafos())])

    def test_flows_read_zero(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(branch_ids_side1=[3])
        ca.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertEqual(float(np.asarray(ca.compute_flows())[0, 3]), 0.)


class TestContingencyAnalysisWrapper(_Base):
    def test_run_reports_what_each_contingency_trips(self):
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_single_contingency(5)
        sa.add_contingency(lines=[5], loads=[self.env.name_load[2]], name="line and its load")
        sa.add_contingency_shunt(0)
        sa.add_contingency_gen(self.env.name_gen[self.gen_alone])
        res = sa.run()
        self.assertEqual(len(res.post_contingency_results), 4)
        first, mixed, shunt, gen = res.post_contingency_results   # in the order they were added
        self.assertEqual(first.element_ids, [5])
        self.assertEqual(first.load_ids, [])
        self.assertEqual(mixed.element_ids, [5])
        self.assertEqual(mixed.load_ids, [2])
        self.assertEqual(mixed.contingency_name, "line and its load")
        self.assertEqual(shunt.shunt_ids, [0])
        self.assertEqual(gen.gen_ids, [self.gen_alone])
        for cont in res.post_contingency_results:
            self.assertTrue(cont.converged)

        # get_flows' rows follow the same order
        _, _, res_v = sa.get_flows()
        self.assertEqual(res_v.shape[0], 4)
        V_line = self._reference_V(branches=(5,))
        V_mixed = self._reference_V(branches=(5,), loads=(2,))
        self.assertGreater(np.abs(V_line - V_mixed).max(), 1e-5)

    def test_a_single_line_by_name_or_numpy_id(self):
        # one line, given alone rather than in a list: by name, by python int, by numpy int
        for lines in (self.env.name_line[5], 5, np.int64(5), np.int32(5)):
            sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
            sa.add_contingency(lines=lines, loads=[2])
            res = sa.run()
            self.assertEqual(len(res.post_contingency_results), 1, msg=repr(lines))
            cont = res.post_contingency_results[0]
            self.assertEqual(cont.element_ids, [5], msg=repr(lines))
            self.assertEqual(cont.load_ids, [2], msg=repr(lines))
            sa.close()


class TestScenarioSweep(_Base):
    def test_shunt_and_load_masks_match_a_one_off_powerflow(self):
        n_load, n_line = self.env.n_load, len(self.grid.get_lines())
        shunt_off = np.zeros((3, self.n_shunt), dtype=bool)
        load_off = np.zeros((3, n_load), dtype=bool)
        line_off = np.zeros((3, n_line), dtype=bool)
        shunt_off[0, 0] = True
        load_off[1, 3] = True
        shunt_off[2, 0] = True
        load_off[2, [2, 4]] = True
        line_off[2, 5] = True
        sweep = ScenarioSweep(self.env)
        sweep.set_contingency_shunts(shunt_off)
        sweep.set_contingency_loads(load_off)
        sweep.set_contingency_lines(line_off)
        sweep.compute()
        V = sweep.get_voltages()
        refs = [self._reference_V(shunts=(0,)),
                self._reference_V(loads=(3,)),
                self._reference_V(branches=(5,), shunts=(0,), loads=(2, 4))]
        for row, ref in enumerate(refs):
            np.testing.assert_allclose(V[row], ref, atol=1e-8, rtol=0.)

    def test_a_load_in_the_mask_and_the_action_is_refused(self):
        load_off = np.zeros((1, self.env.n_load), dtype=bool)
        load_off[0, 2] = True
        act = self.env.action_space({"set_bus": {"loads_id": [(2, -1)]}})
        sweep = ScenarioSweep(self.env)
        sweep.set_contingency_loads(load_off)
        sweep.set_topo_actions([act])
        with self.assertRaises(RuntimeError):
            sweep.compute()


if __name__ == "__main__":
    unittest.main()
