# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Contingencies of ``ContingencyAnalysisCPP`` tripping what is neither a branch, a shunt, a
load nor a generator: a batch row against the one-off powerflow of the grid with the element
really deactivated."""

import copy
import unittest
import warnings

import grid2op
import numpy as np
import pandapower as pp
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP, ScenarioSweepCPP
from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.lightEnv import TopoAction, ElementType
from lightsim2grid import LightSimBackend, ContingencyAnalysis
from lightsim2grid.scenarioSweep import ScenarioSweep

MAX_IT = 30
TOL = 1e-10
LEAF_BUS = 7       # case14: a leaf, behind one transformer
GEN_REMOTE = 3     # the generator on the leaf bus, made a remote regulator where a row must be refused
REG_BUS = 9


def _leaf_branch(grid):
    n_line = len(grid.get_lines())
    branch = [n_line + t.id for t in grid.get_trafos() if LEAF_BUS in (t.bus1_id, t.bus2_id)]
    assert len(branch) == 1
    return branch[0]


class _ContingencyMoreElements(unittest.TestCase):
    """What every kind of element has to do. A subclass says how to build the grid (`_net`
    gives the pandapower network, `_setup_grid` what can only be set on the grid) and how
    to trip an element of it: `ARG` names the list of `add_contingency`, `_deactivate` does it
    on a grid."""
    ARG = None
    SETTER = None

    @staticmethod
    def _container(grid):
        raise NotImplementedError

    @staticmethod
    def _net():
        return pn.case14()

    def _setup_grid(self, grid):
        pass

    def _deactivate(self, grid, el_id):
        raise NotImplementedError

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(self._net())
        self._setup_grid(self.grid)
        self.n_line = len(self.grid.get_lines())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)

    def _add_droop_hvdc(self):
        """an AC emulation HVDC line, what the batch leaves out"""
        big = np.full(1, 1e3)
        self.grid.init_hvdc_lines(np.array([4], dtype=np.int32), np.array([11], dtype=np.int32), [0], [0], np.ones(1),
                                  np.ones(1), [False], [False], np.ones(1), np.ones(1), np.zeros(1), np.zeros(1),
                                  -big, big, -big, big, np.ones(1), np.ones(1), [0], np.array([3.]), np.ones(1),
                                  np.full(1, 12.), [True], np.array([3.]), np.array([2.]), np.full(1, 100.), np.full(1, 100.))
        self.grid.tell_solver_need_reset()

    def _reference(self, el_ids=(), branches=(), gens=(), dc=False, grid=None, hvdcs=()):
        ref = copy.deepcopy(self.grid if grid is None else grid)
        for br in branches:
            if br < self.n_line:
                ref.deactivate_powerline(int(br))
            else:
                ref.deactivate_trafo(int(br) - self.n_line)
        for gen_id in gens:
            ref.deactivate_gen(int(gen_id))
        for hvdc_id in hvdcs:
            ref.deactivate_dcline(int(hvdc_id))
        for el_id in el_ids:
            self._deactivate(ref, int(el_id))
        ref.tell_solver_need_reset()
        V = ref.dc_pf(1.0 * self.v0, MAX_IT, TOL) if dc else ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        self.assertGreater(V.shape[0], 0)
        solved = np.zeros(V.shape[0], dtype=bool)
        solved[np.asarray(ref.id_dc_solver_to_me() if dc else ref.id_ac_solver_to_me())] = True
        V[~solved] = 0.
        return V

    def _analysis(self, conts, grid=None, dc=False, **attributes):
        ca = ContingencyAnalysisCPP(self.grid if grid is None else grid)
        if dc:
            ca.change_algorithm(AlgorithmType.DC_SparseLU)
        for name, value in attributes.items():
            setattr(ca, name, value)
        for el_ids, branches in conts:
            ca.add_contingency(branch_ids=list(branches), **{self.ARG: list(el_ids)})
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        return ca

    def _sweep(self, trips, grid=None, dc=False, **attributes):
        """One row per `(element ids, branches)` of `trips`, as masks (a branch of the grid is a line here)."""
        grid = self.grid if grid is None else grid
        sweep = ScenarioSweepCPP(grid)
        if dc:
            sweep.change_algorithm(AlgorithmType.DC_SparseLU)
        for name, value in attributes.items():
            setattr(sweep, name, value)
        nb_el = len(self._container(grid))
        mask = np.zeros((len(trips), nb_el), dtype=bool)
        line_mask = np.zeros((len(trips), self.n_line), dtype=bool)
        for row, (el_ids, branches) in enumerate(trips):
            mask[row, list(el_ids)] = True
            line_mask[row, list(branches)] = True
        getattr(sweep, self.SETTER)(mask)
        sweep.set_contingency_lines(line_mask)
        return sweep

    def _check_sweep(self, trips, sweep, dc=False, angles_only=False):
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        V = np.asarray(sweep.get_voltages())
        for row, (el_ids, branches) in enumerate(trips):
            self.assertTrue(sweep.converged_mask()[row], (el_ids, branches))
            ref = self._reference(el_ids, branches, dc=dc)
            if angles_only:
                np.testing.assert_allclose(np.angle(V[row]), np.angle(ref), atol=1e-8, rtol=0.)
            else:
                np.testing.assert_allclose(V[row], ref, atol=1e-8, rtol=0.)

    def _row(self, ca, el_ids, branches=()):
        pos = self._arg_position()
        for row, cont in enumerate(ca.my_contingencies()):
            if cont[pos] == sorted(el_ids) and cont[0] == sorted(branches):
                return row
        raise AssertionError(f"no row for {el_ids}, {branches}")

    def _arg_position(self):
        return ["branch_ids", "shunt_ids", "load_ids", "gen_ids", "branch_ids_side1", "branch_ids_side2",
                "sgen_ids", "storage_ids", "svc_ids", "hvdc_ids"].index(self.ARG)


class TestStaticGenerators(_ContingencyMoreElements):
    ARG = "sgen_ids"
    SETTER = "set_contingency_sgens"
    TRIPS = [([0], []), ([1], []), ([0, 1], []), ([0], [5])]

    @staticmethod
    def _container(grid):
        return grid.get_static_generators()

    def test_scenario_sweep(self):
        self._check_sweep(self.TRIPS, self._sweep(self.TRIPS))

    def test_scenario_sweep_dc(self):
        self._check_sweep(self.TRIPS, self._sweep(self.TRIPS, dc=True), dc=True)

    def test_scenario_sweep_takes_the_row_s_own_power_out(self):
        # the tripped unit is given another power in its row: nothing of it must be left
        trips = [([0], []), ([1], [])]
        sweep = self._sweep(trips)
        base = np.array([s.target_p_mw for s in self.grid.get_static_generators()])
        sweep.modify_sgen_p(np.tile(base, (2, 1)) * np.array([[3., 1.], [1., 2.]]))
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        V = np.asarray(sweep.get_voltages())
        for row, (el_ids, branches) in enumerate(trips):
            ref = copy.deepcopy(self.grid)
            kept = 1 - el_ids[0]
            ref.change_p_sgen(kept, float(base[kept] * [[3., 1.], [1., 2.]][row][kept]))
            ref.deactivate_sgen(el_ids[0])
            if branches:
                ref.deactivate_powerline(branches[0])
            ref.tell_solver_need_reset()
            np.testing.assert_allclose(V[row], ref.ac_pf(1.0 * self.v0, MAX_IT, TOL), atol=1e-8, rtol=0.)

    def test_scenario_sweep_redistribute_slack(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        sweep = self._sweep([([0], [])], redistribute_slack=True)
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        ref = copy.deepcopy(self.grid)
        ref.deactivate_sgen(0)
        ref.redistribute_active_power(8.)
        ref.tell_solver_need_reset()
        np.testing.assert_allclose(np.asarray(sweep.get_voltages())[0], ref.ac_pf(1.0 * self.v0, MAX_IT, TOL),
                                   atol=1e-6, rtol=0.)

    def test_scenario_sweep_checks_the_mask(self):
        sweep = ScenarioSweepCPP(self.grid)
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_sgens(np.zeros((1, 3), dtype=bool))

    @staticmethod
    def _net():
        net = pn.case14()
        for bus, p_mw, q_mvar in ((3, 8., 2.), (LEAF_BUS, 5., 0.)):
            pp.create_sgen(net, bus, p_mw=p_mw, q_mvar=q_mvar)
        return net

    def _deactivate(self, grid, el_id):
        grid.deactivate_sgen(el_id)

    def test_matches_a_one_off_powerflow(self):
        ca = self._analysis(self.TRIPS)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            row = self._row(ca, el_ids, branches)
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_dc(self):
        ca = self._analysis(self.TRIPS, dc=True)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            np.testing.assert_allclose(V[self._row(ca, el_ids, branches)],
                                       self._reference(el_ids, branches, dc=True), atol=1e-8, rtol=0.)

    def test_registration(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(sgen_ids=[1, 0])
        ca.add_contingency(sgen_ids=[0, 1])
        self.assertEqual([c[6] for c in ca.my_contingencies()], [[0, 1]])
        self.assertTrue(ca.remove_contingency(sgen_ids=[0, 1]))
        self.assertEqual(len(ca.my_contingencies()), 0)
        with self.assertRaises(RuntimeError):
            ca.add_contingency(sgen_ids=[len(self.grid.get_static_generators())])

    def test_redistribute_slack_shares_what_it_loses(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        ca = self._analysis([([0], [])], redistribute_slack=True)
        ref = copy.deepcopy(self.grid)
        ref.deactivate_sgen(0)
        ref.redistribute_active_power(8.)
        ref.tell_solver_need_reset()
        V_ref = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], V_ref, atol=1e-6, rtol=0.)

    def test_a_tripped_element_is_not_lost_twice_when_stranded(self):
        # sgen 1 stands on the leaf bus: stranded by the branch, tripped or not, it is lost once
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        branch = _leaf_branch(self.grid)
        ca = self._analysis([([1], [branch])], handle_disconnected_grid=True, redistribute_slack=True)
        ca_ref = ContingencyAnalysisCPP(self.grid)
        ca_ref.handle_disconnected_grid = True
        ca_ref.redistribute_slack = True
        ca_ref.add_n1(branch)
        ca_ref.compute(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], np.asarray(ca_ref.get_voltages())[0],
                                   atol=1e-8, rtol=0.)

    def test_fallback(self):
        # an AC emulation line takes the row out of the batch: the fallback trips the rest too
        self._add_droop_hvdc()
        ca = ContingencyAnalysisCPP(self.grid)
        ca.batch_unsupported_fallback = True
        ca.add_contingency(hvdc_ids=[0], sgen_ids=[0])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertTrue(ca.converged_mask()[0])
        self.assertIn("AC emulation", ca.get_batch_unsupported_reasons()[0])
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([0], hvdcs=[0]),
                                   atol=1e-8, rtol=0.)


class TestStorageUnits(_ContingencyMoreElements):
    ARG = "storage_ids"
    SETTER = "set_contingency_storages"
    REGULATING, PQ = 0, 1
    BUS_REGULATING, BUS_PQ = 12, 13
    TRIPS = [([REGULATING], []), ([PQ], []), ([REGULATING, PQ], []), ([PQ], [5])]

    def _setup_grid(self, grid):
        # one unit alone regulating its bus (nothing else does there), one that does not regulate
        grid.init_storages_full(np.array([5., -3.]), np.array([0., 1.]), [True, False],
                                np.array([1.03, 1.]), np.array([-50., -50.]), np.array([50., 50.]),
                                np.array([self.BUS_REGULATING, self.BUS_PQ], dtype=np.int32))
        grid.tell_solver_need_reset()

    def _deactivate(self, grid, el_id):
        grid.deactivate_storage(el_id)

    def test_matches_a_one_off_powerflow(self):
        ca = self._analysis(self.TRIPS)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            row = self._row(ca, el_ids, branches)
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_dc(self):
        # the angles: DC does not solve a magnitude (and a batch keeps its seed's at the unit's bus)
        ca = self._analysis(self.TRIPS, dc=True)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in self.TRIPS:
            np.testing.assert_allclose(np.angle(V[self._row(ca, el_ids, branches)]),
                                       np.angle(self._reference(el_ids, branches, dc=True)), atol=1e-8, rtol=0.)

    def test_the_bus_it_regulated_turns_pq(self):
        ca = self._analysis([([self.REGULATING], [])])
        V = np.asarray(ca.get_voltages())[0]
        self.assertGreater(abs(abs(V[self.BUS_REGULATING]) - 1.03), 1e-4)

    def test_a_unit_of_the_slack_leaves_it(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        self.grid.add_storage_slackbus(self.PQ, 1.)
        self.grid.tell_solver_need_reset()
        ca = self._analysis([([self.PQ], [])], compute_gen_results=True)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([self.PQ]), atol=1e-8, rtol=0.)

    def test_redistribute_slack_shares_what_it_loses(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        self.grid.tell_solver_need_reset()
        ca = self._analysis([([self.PQ], [])], redistribute_slack=True)
        ref = copy.deepcopy(self.grid)
        ref.deactivate_storage(self.PQ)
        ref.redistribute_active_power(3.)  # it was producing 3 MW (-3 in the load convention)
        ref.tell_solver_need_reset()
        V_ref = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], V_ref, atol=1e-6, rtol=0.)

    def test_registration(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(storage_ids=[1, 0])
        self.assertEqual([c[7] for c in ca.my_contingencies()], [[0, 1]])
        self.assertTrue(ca.remove_contingency(storage_ids=[0, 1]))
        with self.assertRaises(RuntimeError):
            ca.add_contingency(storage_ids=[len(self.grid.get_storages())])

    def test_fallback(self):
        self._add_droop_hvdc()
        ca = ContingencyAnalysisCPP(self.grid)
        ca.batch_unsupported_fallback = True
        ca.add_contingency(hvdc_ids=[0], storage_ids=[self.PQ])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertTrue(ca.converged_mask()[0])
        self.assertIn("AC emulation", ca.get_batch_unsupported_reasons()[0])
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0],
                                   self._reference([self.PQ], hvdcs=[0]), atol=1e-8, rtol=0.)

    @staticmethod
    def _container(grid):
        return grid.get_storages()

    def test_scenario_sweep(self):
        self._check_sweep(self.TRIPS, self._sweep(self.TRIPS))

    def test_scenario_sweep_dc(self):
        self._check_sweep(self.TRIPS, self._sweep(self.TRIPS, dc=True), dc=True, angles_only=True)

    def test_scenario_sweep_with_the_islands_handled(self):
        self._check_sweep(self.TRIPS, self._sweep(self.TRIPS, handle_disconnected_grid=True))

    def test_scenario_sweep_a_unit_of_the_slack_leaves_it(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        self.grid.add_storage_slackbus(self.PQ, 1.)
        self.grid.tell_solver_need_reset()
        trips = [([self.PQ], [])]
        self._check_sweep(trips, self._sweep(trips))

    def test_scenario_sweep_a_regulating_unit_on_a_group_bus_is_left_out(self):
        self.grid.set_gen_regulated_bus(GEN_REMOTE, self.BUS_REGULATING)
        self.grid.change_v_gen(GEN_REMOTE, 1.03)
        self.grid.tell_solver_need_reset()
        sweep = self._sweep([([self.REGULATING], []), ([self.PQ], [])])
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertIn("control group", sweep.get_batch_unsupported_reasons()[0])
        self.assertFalse(sweep.converged_mask()[0])
        self.assertTrue(sweep.converged_mask()[1])

    def test_scenario_sweep_refuses_a_unit_also_in_the_action(self):
        self.grid.set_storage_to_subid(np.array([self.BUS_REGULATING, self.BUS_PQ], dtype=np.int32))
        act = TopoAction()
        act.add_element(ElementType.storage, self.PQ, -1)
        sweep = self._sweep([([self.PQ], [])])
        sweep.set_topo_actions([act])
        with self.assertRaisesRegex(RuntimeError, "set_contingency_storages"):
            sweep.compute(1.0 * self.v0, MAX_IT, TOL)

    def test_a_regulating_unit_on_a_group_bus_is_left_out(self):
        # a remote generator holds the bus the unit regulates: their group owns that bus
        self.grid.set_gen_regulated_bus(GEN_REMOTE, self.BUS_REGULATING)
        self.grid.change_v_gen(GEN_REMOTE, 1.03)
        self.grid.tell_solver_need_reset()
        ca = self._analysis([([self.REGULATING], []), ([self.PQ], [])], compute_limit_violations=True)
        row = self._row(ca, [self.REGULATING])
        self.assertIn("control group", ca.get_batch_unsupported_reasons()[row])
        self.assertFalse(ca.converged_mask()[row])
        self.assertTrue(ca.converged_mask()[self._row(ca, [self.PQ])])


class TestStaticVarCompensators(_ContingencyMoreElements):
    ARG = "svc_ids"
    SETTER = "set_contingency_svcs"
    FIXED_Q, OFF_B0, VOLTAGE = 0, 1, 2   # REACTIVE_POWER, OFF with a b0, VOLTAGE

    def _setup_grid(self, grid):
        grid.init_svcs([2, 0, 1], np.array([1., 1., 1.03]), np.array([5., 0., 0.]),
                       np.zeros(3), np.full(3, -0.5), np.full(3, 0.5),
                       np.array([12, 13, 4], dtype=np.int32), np.array([12, 13, 4], dtype=np.int32))
        grid.set_svc_b0(np.array([0.02, 0.03, 0.]))
        grid.tell_solver_need_reset()

    def _deactivate(self, grid, el_id):
        grid.deactivate_svc(el_id)

    def test_matches_a_one_off_powerflow(self):
        trips = [([self.FIXED_Q], []), ([self.OFF_B0], []), ([self.FIXED_Q, self.OFF_B0], []), ([self.FIXED_Q], [5])]
        ca = self._analysis(trips)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in trips:
            row = self._row(ca, el_ids, branches)
            self.assertEqual(ca.get_batch_unsupported_reasons()[row], "")
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    @staticmethod
    def _container(grid):
        return grid.get_svcs()

    def test_scenario_sweep(self):
        trips = [([self.FIXED_Q], []), ([self.OFF_B0], []), ([self.FIXED_Q, self.OFF_B0], []), ([self.FIXED_Q], [5])]
        self._check_sweep(trips, self._sweep(trips))

    def test_scenario_sweep_a_voltage_mode_one(self):
        trips = [([self.VOLTAGE], []), ([self.VOLTAGE, self.FIXED_Q], [5])]
        self._check_sweep(trips, self._sweep(trips))

    def test_scenario_sweep_with_a_standby_one_is_left_out(self):
        self._set_standby()
        sweep = self._sweep([([self.FIXED_Q], []), ([self.OFF_B0], [])])
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertIn("standby", sweep.get_batch_unsupported_reasons()[0])
        self.assertFalse(sweep.converged_mask()[0])
        self.assertTrue(sweep.converged_mask()[1])

    def test_dc_trips_nothing(self):
        # no active power, no b0 in DC: the row is the base case
        ca = self._analysis([([self.VOLTAGE], []), ([self.FIXED_Q], [5])], dc=True)
        self.assertEqual(list(ca.get_batch_unsupported_reasons()), ["", ""])
        V = np.asarray(ca.get_voltages())
        np.testing.assert_allclose(V[self._row(ca, [self.VOLTAGE])], self._reference([], dc=True), atol=1e-8, rtol=0.)

    def _set_standby(self):
        """an automaton on the fixed-Q one, which the batch cannot take out of a row"""
        self.grid.set_svc_standby([True, False, False], np.full(3, 0.9), np.full(3, 1.1))
        self.grid.tell_solver_need_reset()

    def test_a_voltage_mode_one(self):
        trips = [([self.VOLTAGE], []), ([self.VOLTAGE, self.FIXED_Q], []), ([self.VOLTAGE], [5])]
        ca = self._analysis(trips)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in trips:
            row = self._row(ca, el_ids, branches)
            self.assertEqual(ca.get_batch_unsupported_reasons()[row], "")
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_the_regulated_bus_goes_free(self):
        ca = self._analysis([([self.VOLTAGE], [])])
        self.assertGreater(abs(abs(np.asarray(ca.get_voltages())[0][4]) - 1.03), 1e-4)

    def test_a_voltage_mode_one_with_a_slope_regulating_a_remote_bus(self):
        self.grid.init_svcs([1], np.array([1.03]), np.zeros(1), np.array([0.02]), np.full(1, -0.5), np.full(1, 0.5),
                            np.array([13], dtype=np.int32), np.array([8], dtype=np.int32))
        self.grid.tell_solver_need_reset()
        ca = self._analysis([([0], [])])
        self.assertEqual(ca.get_batch_unsupported_reasons()[0], "")
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([0]), atol=1e-8, rtol=0.)

    def test_with_outer_loops(self):
        self.grid.change_algorithm("NROuter_SparseLU")
        self.grid.tell_solver_need_reset()
        trips = [([self.VOLTAGE], []), ([self.FIXED_Q], [])]
        ca = ContingencyAnalysisCPP(self.grid)
        ca.change_algorithm("NROuter_SparseLU")
        for el_ids, branches in trips:
            ca.add_contingency(branch_ids=list(branches), svc_ids=list(el_ids))
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in trips:
            row = self._row(ca, el_ids, branches)
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_with_a_standby_automaton_is_left_out(self):
        self._set_standby()
        ca = self._analysis([([self.FIXED_Q], []), ([self.OFF_B0], [])], compute_limit_violations=True)
        row = self._row(ca, [self.FIXED_Q])
        self.assertIn("standby", ca.get_batch_unsupported_reasons()[row])
        self.assertIn("static var compensator 0", ca.get_batch_unsupported_reasons()[row])
        self.assertFalse(ca.converged_mask()[row])
        self.assertTrue(ca.converged_mask()[self._row(ca, [self.OFF_B0])])

    def test_the_fallback_solves_a_standby_one(self):
        self._set_standby()
        ca = ContingencyAnalysisCPP(self.grid)
        ca.batch_unsupported_fallback = True
        ca.add_contingency(svc_ids=[self.FIXED_Q])
        ca.compute(1.0 * self.v0, MAX_IT, TOL)
        self.assertTrue(ca.converged_mask()[0])
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([self.FIXED_Q]), atol=1e-8, rtol=0.)

    def test_registration(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(svc_ids=[1, 0])
        self.assertEqual([c[8] for c in ca.my_contingencies()], [[0, 1]])
        self.assertTrue(ca.remove_contingency(svc_ids=[0, 1]))
        with self.assertRaises(RuntimeError):
            ca.add_contingency(svc_ids=[len(self.grid.get_svcs())])


class TestHvdcLines(_ContingencyMoreElements):
    ARG = "hvdc_ids"
    SETTER = "set_contingency_hvdcs"
    DROOP, REGULATING, FIXED = 0, 1, 2   # AC emulation, stations regulating at both ends, LCC at fixed P and Q

    def _setup_grid(self, grid):
        big = np.full(3, 1e3)
        grid.init_hvdc_lines(np.array([9, 10, 12], dtype=np.int32), np.array([3, 11, 13], dtype=np.int32),
                             [0, 0, 1], [0, 0, 1],
                             np.ones(3), np.ones(3),
                             [True, True, False], [False, True, False],
                             np.array([1.03, 1.03, 1.]), np.array([1., 1.03, 1.]),
                             np.zeros(3), np.array([3., 0., 0.]),
                             -big, big, -big, big,
                             np.array([1., 1., 0.9]), np.array([1., 1., 0.9]),
                             [0, 0, 0], np.array([10., 5., 3.]),
                             np.ones(3), np.full(3, 12.),
                             [True, False, False], np.array([5., 0., 0.]), np.array([2., 0., 0.]),
                             np.full(3, 100.), np.full(3, 100.))
        grid.tell_solver_need_reset()

    def _deactivate(self, grid, el_id):
        grid.deactivate_dcline(el_id)

    def test_matches_a_one_off_powerflow(self):
        trips = [([self.FIXED], []), ([self.FIXED], [5])]
        ca = self._analysis(trips)
        V = np.asarray(ca.get_voltages())
        for el_ids, branches in trips:
            row = self._row(ca, el_ids, branches)
            self.assertEqual(ca.get_batch_unsupported_reasons()[row], "")
            self.assertTrue(ca.converged_mask()[row])
            np.testing.assert_allclose(V[row], self._reference(el_ids, branches), atol=1e-8, rtol=0.)

    def test_dc(self):
        ca = self._analysis([([self.FIXED], []), ([self.REGULATING], [])], dc=True)
        self.assertEqual(ca.get_batch_unsupported_reasons()[self._row(ca, [self.REGULATING])], "")
        V = np.asarray(ca.get_voltages())
        for el_ids in ([self.FIXED], [self.REGULATING]):
            np.testing.assert_allclose(np.angle(V[self._row(ca, el_ids)]), np.angle(self._reference(el_ids, dc=True)),
                                       atol=1e-8, rtol=0.)

    @staticmethod
    def _container(grid):
        return grid.get_dclines()

    def test_scenario_sweep(self):
        trips = [([self.FIXED], []), ([self.FIXED], [5])]
        self._check_sweep(trips, self._sweep(trips))

    def test_scenario_sweep_dc(self):
        trips = [([self.FIXED], []), ([self.REGULATING], [])]
        self._check_sweep(trips, self._sweep(trips, dc=True), dc=True, angles_only=True)

    def test_scenario_sweep_what_it_leaves_out(self):
        sweep = self._sweep([([self.DROOP], []), ([self.REGULATING], []), ([self.FIXED], [])])
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        reasons = sweep.get_batch_unsupported_reasons()
        self.assertIn("AC emulation", reasons[0])
        self.assertIn("regulates a voltage", reasons[1])
        self.assertEqual(reasons[2], "")
        self.assertTrue(sweep.converged_mask()[2])

    def test_scenario_sweep_redistribute_slack(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        self.grid.tell_solver_need_reset()
        sweep = self._sweep([([self.FIXED], [])], redistribute_slack=True)
        sweep.compute(1.0 * self.v0, MAX_IT, TOL)
        info = self.grid.get_dclines()[self.FIXED]
        ref = copy.deepcopy(self.grid)
        ref.deactivate_dcline(self.FIXED)
        ref.redistribute_active_power(info.target_p1_mw + info.p2_mw)
        ref.tell_solver_need_reset()
        np.testing.assert_allclose(np.asarray(sweep.get_voltages())[0], ref.ac_pf(1.0 * self.v0, MAX_IT, TOL),
                                   atol=1e-6, rtol=0.)

    def test_what_the_batch_leaves_out(self):
        ca = self._analysis([([self.DROOP], []), ([self.REGULATING], []), ([self.FIXED], [])],
                            compute_limit_violations=True)
        reasons = ca.get_batch_unsupported_reasons()
        self.assertIn("AC emulation", reasons[self._row(ca, [self.DROOP])])
        self.assertIn("regulates a voltage", reasons[self._row(ca, [self.REGULATING])])
        self.assertEqual(reasons[self._row(ca, [self.FIXED])], "")
        self.assertFalse(ca.converged_mask()[self._row(ca, [self.DROOP])])

    def test_the_fallback_solves_them(self):
        for line in (self.DROOP, self.REGULATING):
            with self.subTest(line=line):
                ca = ContingencyAnalysisCPP(self.grid)
                ca.batch_unsupported_fallback = True
                ca.add_contingency(hvdc_ids=[line])
                ca.compute(1.0 * self.v0, MAX_IT, TOL)
                self.assertTrue(ca.converged_mask()[0])
                np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], self._reference([line]),
                                           atol=1e-8, rtol=0.)

    def test_redistribute_slack_shares_what_it_loses(self):
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 1.)
        self.grid.tell_solver_need_reset()
        ca = self._analysis([([self.FIXED], [])], redistribute_slack=True)
        # the stations' injections (generator convention), which the slack makes up
        info = self.grid.get_dclines()[self.FIXED]
        ref = copy.deepcopy(self.grid)
        ref.deactivate_dcline(self.FIXED)
        ref.redistribute_active_power(info.target_p1_mw + info.p2_mw)
        ref.tell_solver_need_reset()
        V_ref = ref.ac_pf(1.0 * self.v0, MAX_IT, TOL)
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[0], V_ref, atol=1e-6, rtol=0.)

    def test_registration(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_contingency(hvdc_ids=[1, 0])
        self.assertEqual([c[9] for c in ca.my_contingencies()], [[0, 1]])
        self.assertTrue(ca.remove_contingency(hvdc_ids=[0, 1]))
        with self.assertRaises(RuntimeError):
            ca.add_contingency(hvdc_ids=[len(self.grid.get_dclines())])


class TestRemoteVoltageControllers(unittest.TestCase):
    """Generators of a voltage-control group, tripped by a row of the batch: gone from the
    group, as in a one-off solve of the grid without them."""
    REG_BUS = 9
    GROUP = (1, 2, 3)
    TRIPS = ([1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3])

    def _grid(self, group=GROUP, algo=None, q_limits=None, reg_bus=None):
        reg_bus = self.REG_BUS if reg_bus is None else reg_bus
        net = pn.case14()
        for gen_pos, (min_q, max_q) in (q_limits or {}).items():
            net.gen.loc[net.gen.index[gen_pos - 1], ["min_q_mvar", "max_q_mvar"]] = [min_q, max_q]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(net)
        for gen_id in group:
            grid.set_gen_regulated_bus(gen_id, reg_bus)
            grid.change_v_gen(gen_id, 1.04)
        if algo is not None:
            grid.change_algorithm(algo)
        grid.tell_solver_need_reset()
        return grid

    @staticmethod
    def _one_off(grid, gens, v0):
        ref = copy.deepcopy(grid)
        for gen_id in gens:
            ref.deactivate_gen(int(gen_id))
        ref.tell_solver_need_reset()
        return ref.ac_pf(1.0 * v0, MAX_IT, TOL)

    def _check(self, grid, trips, algo=None, **attributes):
        v0 = np.ones(grid.total_bus(), dtype=complex)
        ca = ContingencyAnalysisCPP(grid)
        if algo is not None:
            ca.change_algorithm(algo)
        for name, value in attributes.items():
            setattr(ca, name, value)
        for gens in trips:
            ca.add_contingency(gen_ids=list(gens))
        ca.compute(1.0 * v0, MAX_IT, TOL)
        rows = {tuple(c[3]): row for row, c in enumerate(ca.my_contingencies())}
        V = np.asarray(ca.get_voltages())
        for gens in trips:
            with self.subTest(gens=gens):
                row = rows[tuple(gens)]
                self.assertEqual(ca.get_batch_unsupported_reasons()[row], "")
                self.assertTrue(ca.converged_mask()[row])
                np.testing.assert_allclose(V[row], self._one_off(grid, gens, v0), atol=1e-8, rtol=0.)

    def test_matches_a_one_off_powerflow(self):
        self._check(self._grid(), self.TRIPS)

    def test_multithreaded(self):
        self._check(self._grid(), self.TRIPS, nb_thread=3)

    def test_with_the_islands_handled(self):
        self._check(self._grid(), self.TRIPS, handle_disconnected_grid=True)

    def test_a_group_of_one(self):
        self._check(self._grid(group=(3,)), [[3]])

    def test_with_outer_loops(self):
        self._check(self._grid(algo="NROuter_SparseLU"), self.TRIPS, algo="NROuter_SparseLU")

    def test_scenario_sweep(self):
        grid = self._grid()
        v0 = np.ones(grid.total_bus(), dtype=complex)
        mask = np.zeros((len(self.TRIPS), len(grid.get_generators())), dtype=bool)
        for row, gens in enumerate(self.TRIPS):
            mask[row, gens] = True
        sweep = ScenarioSweepCPP(grid)
        sweep.set_contingency_gens(mask)
        sweep.compute(1.0 * v0, MAX_IT, TOL)
        V = np.asarray(sweep.get_voltages())
        for row, gens in enumerate(self.TRIPS):
            self.assertEqual(sweep.get_batch_unsupported_reasons()[row], "")
            np.testing.assert_allclose(V[row], self._one_off(grid, gens, v0), atol=1e-8, rtol=0.)

    def test_the_regulated_bus_goes_free_when_nobody_is_left(self):
        grid = self._grid()
        ca = ContingencyAnalysisCPP(grid)
        ca.add_contingency(gen_ids=list(self.GROUP))
        ca.add_contingency(gen_ids=[1])
        ca.compute(np.ones(grid.total_bus(), dtype=complex), MAX_IT, TOL)
        V = {tuple(c[3]): v for c, v in zip(ca.my_contingencies(), np.asarray(ca.get_voltages()))}
        self.assertGreater(abs(abs(V[self.GROUP][self.REG_BUS]) - 1.04), 1e-3)
        self.assertAlmostEqual(abs(V[(1,)][self.REG_BUS]), 1.04, places=8)

    def test_a_local_generator_of_a_group(self):
        # generator 1 stands on the bus generator 3 regulates remotely: both are of the group
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bus_of_1 = init_from_pandapower(pn.case14()).get_generators()[1].bus_id
        grid = self._grid(group=(1, 3), reg_bus=bus_of_1)
        self._check(grid, [[1], [3], [1, 3]])

    def test_outer_loops_with_a_controller_at_its_reactive_limit(self):
        # generator 1 cannot give the group what it asks: ReactiveLimits holds it
        limits = {1: (-1., 1.)}
        self._check(self._grid(algo="NROuter_SparseLU", q_limits=limits), self.TRIPS, algo="NROuter_SparseLU")

    def test_outer_loops_a_tripped_controller_whose_range_excludes_0(self):
        # Q = 0 is outside what generator 1 may produce: it is gone, nothing checks it
        limits = {1: (5., 6.)}
        self._check(self._grid(algo="NROuter_SparseLU", q_limits=limits), [[1], [1, 2], [1, 2, 3]],
                    algo="NROuter_SparseLU")

    def test_outer_loops_multithreaded(self):
        self._check(self._grid(algo="NROuter_SparseLU", q_limits={1: (-1., 1.)}), self.TRIPS,
                    algo="NROuter_SparseLU", nb_thread=3)

    def test_a_group_regulating_a_slack_bus(self):
        # gen 1 regulates the slack bus with the slack generator: one of them may go, not both
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(pn.case14())
        slack = next(g for g in grid.get_generators() if g.is_slack)
        grid.set_gen_regulated_bus(1, slack.bus_id)
        grid.change_v_gen(1, slack.target_vm_pu)
        grid.tell_solver_need_reset()
        v0 = np.ones(grid.total_bus(), dtype=complex)
        ca = ContingencyAnalysisCPP(grid)
        ca.add_contingency(gen_ids=[1])
        ca.add_contingency(gen_ids=[1, slack.id])
        ca.compute(1.0 * v0, MAX_IT, TOL)
        rows = {tuple(c[3]): row for row, c in enumerate(ca.my_contingencies())}
        one = rows[(1,)]
        self.assertEqual(ca.get_batch_unsupported_reasons()[one], "")
        np.testing.assert_allclose(np.asarray(ca.get_voltages())[one], self._one_off(grid, [1], v0), atol=1e-8, rtol=0.)
        both = rows[tuple(sorted((1, slack.id)))]
        self.assertIn("slack bus", ca.get_batch_unsupported_reasons()[both])
        self.assertFalse(ca.converged_mask()[both])

    def test_consecutive_computes_without_a_trip(self):
        # a second compute() with no such row is back to the plain pattern
        grid = self._grid()
        v0 = np.ones(grid.total_bus(), dtype=complex)
        ca = ContingencyAnalysisCPP(grid)
        ca.add_contingency(gen_ids=[1])
        ca.compute(1.0 * v0, MAX_IT, TOL)
        ca.remove_contingency(gen_ids=[1])
        ca.add_n1(0)
        ca.compute(1.0 * v0, MAX_IT, TOL)
        self.assertTrue(ca.converged_mask()[0])


class TestWrappers(unittest.TestCase):
    """`ContingencyAnalysis.add_contingency` and the `ScenarioSweep` setters, on a grid2op
    environment: what they forward is what the c++ classes do."""
    ENV_NAME = "l2rpn_case14_sandbox"

    def _env(self, env_name=ENV_NAME):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            env = grid2op.make(env_name, test=True, backend=LightSimBackend())
        self.addCleanup(env.close)
        return env

    def _setup_sgens(self, env):
        env.backend._grid.init_sgens(np.array([8., 5.]), np.array([2., 0.]), np.zeros(2), np.full(2, 50.),
                                     np.full(2, -50.), np.full(2, 50.), np.array([3, 12], dtype=np.int32))
        env.backend._grid.tell_solver_need_reset()

    def _check(self, env, keyword, cpp_kw, field, el_id, setter, nb_el):
        ca = ContingencyAnalysis(env, compute_limit_violations=True)
        ca.add_contingency(name="cont", **{keyword: [el_id]})
        result = ca.run().post_contingency_results[0]
        self.assertEqual(getattr(result, field), [el_id])
        self.assertEqual(result.contingency_name, "cont")
        self.assertTrue(result.converged)
        backend = ca._ls_backend
        cpp = ContingencyAnalysisCPP(backend._grid)
        cpp.add_contingency(**{cpp_kw: [el_id]})
        cpp.compute(1.0 * backend.V, backend.max_it, backend.tol)
        np.testing.assert_allclose(ca.compute_V()[0], np.asarray(cpp.get_voltages())[0], atol=1e-8, rtol=0.)
        # registered once, named or not
        ca.add_contingency(**{keyword: el_id})
        self.assertEqual(len(ca.computer.my_contingencies()), 1)
        with self.assertRaises(RuntimeError):
            ca.add_contingency(**{keyword: [nb_el]})

        sweep = ScenarioSweep(env)
        mask = np.zeros((1, nb_el), dtype=bool)
        mask[0, el_id] = True
        getattr(sweep, setter)(mask)
        sweep.compute(1.0 * env.backend.V)
        self.assertTrue(sweep.computer.converged_mask()[0])
        with self.assertRaises(RuntimeError):
            getattr(ScenarioSweep(env), setter)(np.zeros((1, nb_el + 1), dtype=bool))

    def test_static_generators(self):
        env = self._env()
        self._setup_sgens(env)
        self._check(env, "sgens", "sgen_ids", "sgen_ids", 1, "set_contingency_sgens", 2)

    def test_storage_units(self):
        env = self._env("educ_case14_storage")
        self._check(env, "storages", "storage_ids", "storage_ids", 1, "set_contingency_storages", env.n_storage)

    def test_storage_units_by_name(self):
        env = self._env("educ_case14_storage")
        ca = ContingencyAnalysis(env)
        ca.add_contingency(storages=[env.name_storage[1]])
        self.assertEqual([c[7] for c in ca.computer.my_contingencies()], [[1]])
        with self.assertRaises(RuntimeError):
            ca.add_contingency(storages=["not a storage unit"])

    def test_static_var_compensators(self):
        env = self._env()
        grid = env.backend._grid
        grid.init_svcs([2, 2], np.array([1., 1.]), np.array([5., 3.]), np.zeros(2), np.full(2, -0.5), np.full(2, 0.5),
                       np.array([12, 13], dtype=np.int32), np.array([12, 13], dtype=np.int32))
        grid.tell_solver_need_reset()
        self._check(env, "svcs", "svc_ids", "svc_ids", 1, "set_contingency_svcs", 2)

    def test_hvdc_lines(self):
        env = self._env()
        grid = env.backend._grid
        big = np.full(2, 1e3)
        grid.init_hvdc_lines(np.array([10, 12], dtype=np.int32), np.array([11, 13], dtype=np.int32),
                             [1, 1], [1, 1], np.ones(2), np.ones(2),
                             [False, False], [False, False],
                             np.ones(2), np.ones(2), np.zeros(2), np.zeros(2),
                             -big, big, -big, big,
                             np.full(2, 0.9), np.full(2, 0.9),
                             [0, 0], np.array([5., 3.]),
                             np.ones(2), np.full(2, 12.),
                             [False, False], np.zeros(2), np.zeros(2),
                             np.full(2, 100.), np.full(2, 100.))
        grid.tell_solver_need_reset()
        self._check(env, "hvdcs", "hvdc_ids", "hvdc_ids", 1, "set_contingency_hvdcs", 2)

    def test_an_empty_contingency_is_refused(self):
        with self.assertRaises(RuntimeError):
            ContingencyAnalysis(self._env()).add_contingency()

    def test_the_shortcuts(self):
        env = self._env()
        self._setup_sgens(env)
        ca = ContingencyAnalysis(env)
        ca.add_contingency_sgen(0)
        self.assertEqual([c[6] for c in ca.computer.my_contingencies()], [[0]])


if __name__ == "__main__":
    unittest.main()
