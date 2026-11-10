# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A contingency the batch does not support is reported BATCH_UNSUPPORTED, the others still run (#280).

Tripping a static var compensator with a standby automaton is not supported by the batch yet (its
voltage control lives in Jacobian rows and columns the batch does not reserve per row). Such a
contingency does not make ``compute()`` raise for the whole batch: that row alone is reported
``LimitViolationType.BATCH_UNSUPPORTED``, with its reason in ``get_batch_unsupported_reasons()``. With
``batch_unsupported_fallback``, those rows are then solved one at a time by ``ac_pf`` on a copy of the
grid, and must match that one-off solve.

Fixture: pandapower case14 with a standby SVC on the leaf bus 7, whose generator is off.
"""

import copy
import unittest
import warnings

import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP, LimitViolationType,  # noqa: E402
                                             ViolationCategory, ViolationElementType)

MAX_IT = 30
TOL = 1e-11
GEN_ON_LEAF = 3   # on bus 7, off
TRAFO_BEHIND = 3  # buses 6-7: bus 7 is behind it, alone
LOCAL_GEN = 1     # a generator regulating its own bus
SVC = 0           # with a standby automaton, on bus 7: what the batch does not support


def _grid():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case14())
    grid.deactivate_gen(GEN_ON_LEAF)
    grid.init_svcs([2], np.ones(1), np.array([5.]), np.zeros(1), np.full(1, -0.5), np.full(1, 0.5),
                   np.array([7], dtype=np.int32), np.array([7], dtype=np.int32))
    grid.set_svc_standby([True], np.full(1, 0.9), np.full(1, 1.1))
    grid.tell_solver_need_reset()
    return grid


def _grid_with_half_open_hvdc():
    """the same fixture, with an HVDC line whose end 1 is open in the base case, on a bus nothing else
    touches: the state the main-component clean-up leaves a line in"""
    net = pn.case14()
    dead_bus = pp.create_bus(net, vn_kv=net.bus.vn_kv.iloc[3])
    pp.create_dcline(net, from_bus=dead_bus, to_bus=3, p_mw=0.0, loss_percent=0.0, loss_mw=0.0,
                     vm_from_pu=1.0, vm_to_pu=1.02)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(net)
    grid.deactivate_dcline_side1(0)
    grid.deactivate_gen(GEN_ON_LEAF)
    grid.init_svcs([2], np.ones(1), np.array([5.]), np.zeros(1), np.full(1, -0.5), np.full(1, 0.5),
                   np.array([7], dtype=np.int32), np.array([7], dtype=np.int32))
    grid.set_svc_standby([True], np.full(1, 0.9), np.full(1, 1.1))
    grid.tell_solver_need_reset()
    return grid


class TestBatchUnsupportedContingencies(unittest.TestCase):
    def setUp(self):
        self.grid = _grid()
        self.n_line = len(self.grid.get_lines())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)
        self.buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _reference(self, branches=(), gens=()):
        ref = copy.deepcopy(self.grid)
        for br in branches:
            if br < self.n_line:
                ref.deactivate_powerline(int(br))
            else:
                ref.deactivate_trafo(int(br - self.n_line))
        for g in gens:
            ref.deactivate_gen(int(g))
        return ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)

    def _analysis(self, nb_thread=1, handle_disconnected_grid=False):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.handle_disconnected_grid = handle_disconnected_grid
        ca.nb_thread = nb_thread
        for line_id in (0, 3, 7):
            ca.add_n1(line_id)
        ca.add_contingency(svc_ids=[SVC])
        ca.add_contingency_gen(LOCAL_GEN)
        ca.add_contingency(branch_ids=[5], svc_ids=[SVC])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        return ca

    @staticmethod
    def _row(ca, branches=(), gens=(), svcs=()):
        for row, cont in enumerate(ca.my_contingencies()):
            if sorted(cont[0]) == sorted(branches) and sorted(cont[3]) == sorted(gens) and sorted(cont[8]) == sorted(svcs):
                return row
        raise AssertionError(f"no row for branches {branches}, generators {gens}")

    def _check(self, ca):
        reasons = ca.get_batch_unsupported_reasons()
        conts = ca.my_contingencies()
        self.assertEqual(len(reasons), len(conts))
        unsupported = {self._row(ca, svcs=[SVC]), self._row(ca, branches=[5], svcs=[SVC])}
        for row, cont in enumerate(conts):
            with self.subTest(row=row, contingency=cont):
                if row in unsupported:
                    self.assertIn("standby", reasons[row])
                    self.assertIn(f"static var compensator {SVC}", reasons[row])
                    self.assertFalse(ca.converged()[row])
                    self.assertFalse(ca.converged_mask()[row])
                    viol = ca.get_violations()[row]
                    self.assertEqual(len(viol), 1)
                    self.assertEqual(viol[0].violation_type, LimitViolationType.BATCH_UNSUPPORTED)
                    self.assertEqual(viol[0].element_type, ViolationElementType.GRID)
                    self.assertEqual(viol[0].category, ViolationCategory.SOLVER)
                    self.assertTrue(np.all(ca.get_voltages()[row] == 0.))
                else:
                    self.assertEqual(reasons[row], "")
                    self.assertTrue(ca.converged()[row])
                    ref = self._reference(branches=cont[0], gens=cont[3])
                    self.assertGreater(ref.shape[0], 0)
                    np.testing.assert_allclose(ca.get_voltages()[row][self.buses], ref[self.buses],
                                               rtol=1e-8, atol=1e-8)

    def test_compute_does_not_raise(self):
        self._check(self._analysis())

    def test_multithreaded(self):
        self._check(self._analysis(nb_thread=3))

    def test_handle_disconnected_grid(self):
        self._check(self._analysis(handle_disconnected_grid=True))

    def test_flows_of_an_unsupported_row_are_nan(self):
        """as those of a row the connectivity check skips"""
        ca = self._analysis()
        amps = ca.compute_flows()
        row = self._row(ca, svcs=[SVC])
        self.assertTrue(np.all(np.isnan(amps[row])))
        other = self._row(ca, gens=[LOCAL_GEN])
        self.assertTrue(np.any(amps[other] != 0.))

    def test_a_split_grid_stays_not_simulated(self):
        """The connectivity check wins: a contingency that strands the SVC's own
        bus is NOT_SIMULATED, whatever else it trips -- no one-off solve can do better."""
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.add_contingency(branch_ids=[self.n_line + TRAFO_BEHIND], svc_ids=[SVC])
        ca.add_n1(0)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        row = self._row(ca, branches=[self.n_line + TRAFO_BEHIND], svcs=[SVC])
        self.assertEqual(ca.get_batch_unsupported_reasons()[row], "")
        self.assertEqual([v.violation_type for v in ca.get_violations()[row]],
                         [LimitViolationType.NOT_SIMULATED])

    def test_reasons_follow_the_contingency_list(self):
        """Registering again (L2 cleared) and recomputing settles the reasons anew."""
        ca = self._analysis()
        self.assertEqual(sum(bool(r) for r in ca.get_batch_unsupported_reasons()), 2)
        ca.remove_contingency(svc_ids=[SVC])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertEqual(sum(bool(r) for r in ca.get_batch_unsupported_reasons()), 1)
        self.assertEqual(len(ca.get_batch_unsupported_reasons()), len(ca.my_contingencies()))

    def test_without_unsupported_contingency_nothing_changes(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_n1(0)
        ca.add_contingency_gen(LOCAL_GEN)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertEqual(list(ca.get_batch_unsupported_reasons()), ["", ""])
        self.assertTrue(all(ca.converged_mask()))


class TestBatchUnsupportedFallback(unittest.TestCase):
    """``batch_unsupported_fallback``: the rows the batch does not support are solved after it, one at a
    time, by ``ac_pf`` on a copy of the grid with the contingency's elements out."""

    def setUp(self):
        self.grid = _grid()
        self.n_line = len(self.grid.get_lines())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)
        self.buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _reference(self, grid, cont, init_from_n=False):
        """the one-off solve of that contingency, as the fallback is documented to run it"""
        ref = copy.deepcopy(grid)
        v_start = 1.0 * self.v0
        outer = "NROuter" in ref.get_algo_name()
        if init_from_n:
            v_start = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
            self.assertGreater(v_start.shape[0], 0)
            if outer:
                ref.set_resume_outer_state(True)
        branches, shunts, loads, gens, side1, side2 = (list(c) for c in cont[:6])
        for br in branches:
            if br < self.n_line:
                ref.deactivate_powerline(int(br))
            else:
                ref.deactivate_trafo(int(br - self.n_line))
        for br in side1:
            if br < self.n_line:
                ref.deactivate_powerline_side1(int(br))
            else:
                ref.deactivate_trafo_side1(int(br - self.n_line))
        for g in gens:
            ref.deactivate_gen(int(g))
        for s in shunts:
            ref.deactivate_shunt(int(s))
        for s in cont[8]:
            ref.deactivate_svc(int(s))
        for ld in loads:
            ref.deactivate_load(int(ld))
        return ref, ref.ac_pf(v_start, MAX_IT, TOL)

    def _analysis(self, grid=None, nb_thread=1, init_from_n=False, fallback=True):
        grid = self.grid if grid is None else grid
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_limit_violations = True
        ca.batch_unsupported_fallback = fallback
        ca.init_from_n_powerflow = init_from_n
        ca.nb_thread = nb_thread
        for line_id in (0, 3):
            ca.add_n1(line_id)
        ca.add_contingency(svc_ids=[SVC])
        ca.add_contingency(branch_ids=[5], svc_ids=[SVC])
        ca.add_contingency(branch_ids_side1=[7], svc_ids=[SVC])
        ca.add_contingency(load_ids=[2], svc_ids=[SVC])
        ca.add_contingency_gen(LOCAL_GEN)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        return ca

    def _check(self, ca, grid=None, init_from_n=False):
        grid = self.grid if grid is None else grid
        reasons = ca.get_batch_unsupported_reasons()
        amps = ca.compute_flows()
        nb_fallback = 0
        for row, cont in enumerate(ca.my_contingencies()):
            with self.subTest(row=row, contingency=cont):
                self.assertEqual(bool(reasons[row]), SVC in cont[8])
                nb_fallback += bool(reasons[row])
                ref, V_ref = self._reference(grid, cont, init_from_n)
                if V_ref.shape[0] == 0 and not reasons[row]:
                    # a row of the batch the one-off solve does not solve either (the outer
                    # loops, at this tolerance): nothing to compare
                    self.assertFalse(ca.converged()[row])
                    continue
                self.assertGreater(V_ref.shape[0], 0)
                self.assertTrue(ca.converged()[row])
                self.assertFalse(any(v.element_type == ViolationElementType.GRID
                                     for v in ca.get_violations()[row]))
                np.testing.assert_allclose(ca.get_voltages()[row][self.buses], V_ref[self.buses],
                                           rtol=1e-7, atol=1e-7)
                if reasons[row]:
                    # the flows read off the row's voltages, as for any other row
                    a_or = np.concatenate([ref.get_line_res1()[3], ref.get_trafo_res1()[3]])
                    out = set(cont[0]) | set(cont[4]) | set(cont[5])
                    keep = [br for br in range(a_or.shape[0]) if br not in out]
                    np.testing.assert_allclose(amps[row][keep], a_or[keep], rtol=1e-6, atol=1e-6)
        self.assertEqual(nb_fallback, 4)

    def test_every_row_is_solved(self):
        self._check(self._analysis())

    def test_multithreaded(self):
        self._check(self._analysis(nb_thread=3))

    def test_init_from_n_powerflow(self):
        self._check(self._analysis(init_from_n=True), init_from_n=True)

    def test_off_by_default(self):
        ca = ContingencyAnalysisCPP(self.grid)
        self.assertFalse(ca.batch_unsupported_fallback)
        ca = self._analysis(fallback=False)
        self.assertEqual(sum(not c for c in ca.converged()), 4)

    def test_toggling_back_off(self):
        ca = self._analysis()
        self.assertTrue(all(ca.converged()))
        ca.batch_unsupported_fallback = False
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertEqual(sum(not c for c in ca.converged()), 4)
        row = TestBatchUnsupportedContingencies._row(ca, svcs=[SVC])
        self.assertEqual([v.violation_type for v in ca.get_violations()[row]], [LimitViolationType.BATCH_UNSUPPORTED])

    def test_source_grid_untouched(self):
        V_before = copy.deepcopy(self.grid).ac_pf(1.0 * self.v0, MAX_IT, TOL)
        self._analysis()
        statuses = [l.connected_global for l in self.grid.get_lines()]
        self.assertTrue(all(statuses))
        self.assertTrue(all(g.connected or g.id == GEN_ON_LEAF for g in self.grid.get_generators()))
        np.testing.assert_allclose(self.grid.ac_pf(1.0 * self.v0, MAX_IT, TOL), V_before, atol=1e-12)

    def test_with_outer_loops(self):
        for init_from_n in (False, True):
            with self.subTest(init_from_n=init_from_n):
                grid = _grid()
                grid.change_algorithm("NROuter_SparseLU")
                ca = self._analysis(grid=grid, init_from_n=init_from_n)
                self.assertIn("NROuter", ca.get_algo_name())
                self._check(ca, grid=grid, init_from_n=init_from_n)

    def test_a_split_grid_keeps_its_main_component(self):
        """handle_disconnected_grid: the batch solves a split grid on its main component, and so
        does the fallback (consider_only_main_component on its own copy of the grid). The batch's
        own split row is checked against the same reference, so the two conventions agree."""
        trafo = self.n_line + TRAFO_BEHIND  # the SVC's own bus is behind it
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.handle_disconnected_grid = True
        ca.batch_unsupported_fallback = True
        ca.add_contingency(branch_ids=[trafo], svc_ids=[SVC])
        ca.add_contingency(branch_ids=[trafo], gen_ids=[LOCAL_GEN])  # the batch's own
        ca.add_n1(0)
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        amps = ca.compute_flows()
        for kw in ({"svcs": [SVC]}, {"gens": [LOCAL_GEN]}):
            with self.subTest(**kw):
                row = TestBatchUnsupportedContingencies._row(ca, branches=[trafo], **kw)
                self.assertEqual(bool(ca.get_batch_unsupported_reasons()[row]), "svcs" in kw)
                self.assertTrue(ca.converged()[row])
                ref = copy.deepcopy(self.grid)
                ref.deactivate_trafo(TRAFO_BEHIND)
                if "svcs" in kw:
                    ref.deactivate_svc(SVC)
                else:
                    ref.deactivate_gen(LOCAL_GEN)
                ref.consider_only_main_component(False)
                V_ref = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
                self.assertGreater(V_ref.shape[0], 0)
                # the stranded bus reads 0, as a masked one (ac_pf hands back its input there)
                main = [b for b in self.buses if b != 7]
                np.testing.assert_allclose(ca.get_voltages()[row][main], V_ref[main], rtol=1e-7, atol=1e-7)
                self.assertEqual(ca.get_voltages()[row][7], 0.)
                a_or = np.concatenate([ref.get_line_res1()[3], ref.get_trafo_res1()[3]])
                live = [br for br in range(a_or.shape[0]) if br != trafo]
                np.testing.assert_allclose(np.nan_to_num(amps[row][live]), a_or[live], rtol=1e-6, atol=1e-6)
        # the source grid of the analysis kept its element: only the fallback's copy lost it
        self.assertTrue(self.grid.get_svcs()[SVC].connected)

    def test_redistribute_slack_without_outer_loops_is_refused(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.compute_limit_violations = True
        ca.redistribute_slack = True
        ca.batch_unsupported_fallback = True
        ca.add_contingency(svc_ids=[SVC])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertFalse(ca.converged()[0])
        self.assertIn("redistribute_slack", ca.get_batch_unsupported_reasons()[0])
        self.assertEqual([v.violation_type for v in ca.get_violations()[0]], [LimitViolationType.BATCH_UNSUPPORTED])

    def test_undoing_a_half_open_hvdc_line_keeps_its_open_end_open(self):
        """Issue #281, PR #292: the fallback undoes each row on a grid copy shared with the next rows.
        Putting back an HVDC line that was already half-open in the base case reconnected its open end
        too, on a bus nothing else touches: every later row on that copy had a singular Jacobian and
        diverged, although each of them converges alone."""
        grid = _grid_with_half_open_hvdc()
        v0 = np.ones(grid.total_bus(), dtype=complex)
        rows = [dict(svc_ids=[SVC], hvdc_ids=[0]),   # the row whose undo matters
                dict(svc_ids=[SVC], load_ids=[2]),
                dict(svc_ids=[SVC], branch_ids=[1])]

        def run(contingencies):
            ca = ContingencyAnalysisCPP(grid)
            ca.batch_unsupported_fallback = True
            for cont in contingencies:
                ca.add_contingency(**cont)
            ca.compute(1.0 * v0, MAX_IT, TOL)
            return ca, np.asarray(ca.converged_mask()).astype(bool), np.asarray(ca.get_voltages())

        ca, converged, voltages = run(rows)
        # the object sorts its rows: check that the row undoing the HVDC line is not the last one
        hvdc_rows = [i for i, cont in enumerate(ca.my_contingencies()) if len(cont[9])]
        self.assertEqual(hvdc_rows, [0])
        self.assertTrue(all(ca.get_batch_unsupported_reasons()))  # all three went through the fallback
        self.assertTrue(converged.all(), f"rows converged: {converged}")
        for row, cont in enumerate(ca.my_contingencies()):
            kept = dict(svc_ids=list(cont[8]), hvdc_ids=list(cont[9]), load_ids=list(cont[2]), branch_ids=list(cont[0]))
            _, alone_converged, alone = run([kept])
            self.assertTrue(alone_converged[0])
            np.testing.assert_allclose(voltages[row], alone[0], atol=1e-9, rtol=0.)
        # the line stays half-open on the source grid, too
        self.assertFalse(grid.get_dclines()[0].connected1)


if __name__ == "__main__":
    unittest.main()
