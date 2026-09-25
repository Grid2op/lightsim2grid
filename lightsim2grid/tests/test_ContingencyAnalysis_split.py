# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the `handle_disconnected_grid` mode of the contingency analysis.

When this mode is ON, a contingency that splits the grid into several connected
components is no longer skipped: the largest component is solved while the buses
of the other component(s) are "masked" (their NR equations are forced to identity
and their voltage reported as 0), without any extra symbolic factorization.
"""

import unittest
import warnings
import numpy as np
import pandapower as pp

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP, LimitViolationType
from lightsim2grid.algorithm import AlgorithmType


def _line(net, f, t):
    pp.create_line_from_parameters(net, f, t, length_km=1.,
                                   r_ohm_per_km=0.1, x_ohm_per_km=0.3,
                                   c_nf_per_km=0., max_i_ka=1.)


def _build_radial(n_bus, with_island=True):
    """A simple radial grid 0-1-2(-3). Disconnecting the last line (id 2)
    isolates bus 3. `with_island=False` builds the reference grid (no bus 3)."""
    net = pp.create_empty_network(sn_mva=1.)
    nb = n_bus if with_island else 3
    for _ in range(nb):
        pp.create_bus(net, vn_kv=20.)
    pp.create_ext_grid(net, 0, vm_pu=1.0)
    _line(net, 0, 1)               # line 0
    _line(net, 1, 2)               # line 1
    pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
    pp.create_load(net, 2, p_mw=1.0, q_mvar=0.2)
    if with_island:
        _line(net, 2, 3)           # line 2 (splits the grid -> isolates bus 3)
        pp.create_load(net, 3, p_mw=1.0, q_mvar=0.2)
    return net


class TestStrandedRemoteControllerEveryAlgo(unittest.TestCase):
    """case14, generator 3 (on bus 7) regulating bus 9 remotely; tripping trafo 3
    (buses 6-7) strands the controller alone on bus 7. In `handle_disconnected_grid`
    mode the stranded controller's row is repurposed at constant sparsity, which
    moves a pivot: KLU's refactorize halted on it and the row came back diverged under
    NR_KLU / NRSing_KLU while SparseLU (re-pivoting) and NRRefactorRetry_KLU solved it.
    The batch now enables the numeric-factorize fallback on its algorithm."""
    def setUp(self):
        import pandapower.networks as pn
        self.max_it = 30
        self.tol = 1e-8
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(pn.case14())
            ref = init_from_pandapower(pn.case14())
        self.nb = self.grid.total_bus()
        self.nb_line = len(self.grid.get_lines())
        self.grid.set_gen_regulated_bus(3, 9)
        self.V0 = self.grid.ac_pf(np.ones(self.nb, dtype=complex), self.max_it, self.tol)
        assert self.V0.shape[0] == self.nb, "base case did not converge"
        # single-shot reference of the post-contingency grid: trafo out, and the
        # controller it stranded out with it
        ref.set_gen_regulated_bus(3, 9)
        ref.deactivate_trafo(3)
        ref.deactivate_gen(3)
        self.Vref = ref.ac_pf(np.ones(self.nb, dtype=complex), self.max_it, self.tol)
        assert self.Vref.shape[0] == self.nb, "reference did not converge"
        self.live = [b for b in range(self.nb) if b != 7]

    def test_every_nr_algorithm(self):
        names = [nm for nm in ContingencyAnalysisCPP(self.grid).available_algorithm_names()
                 if nm.startswith("NR")]
        assert "NR_SparseLU" in names
        for name in names:
            with self.subTest(name):
                SA = ContingencyAnalysisCPP(self.grid, True)
                SA.change_algorithm(name)
                SA.add_n1(self.nb_line + 3)
                SA.handle_disconnected_grid = True
                SA.compute(1.0 * self.V0, self.max_it, self.tol)
                assert list(SA.converged()) == [True], f"{name}: row reported diverged"
                V = SA.get_voltages()[0]
                assert abs(V[7]) == 0., "the stranded bus reports 0"
                np.testing.assert_allclose(V[self.live], self.Vref[self.live], rtol=0., atol=1e-6)


class TestContingencySplitMode(unittest.TestCase):
    def setUp(self):
        self.max_it = 30
        self.tol = 1e-8
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(_build_radial(4))
            self.ref = init_from_pandapower(_build_radial(4, with_island=False))
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        # reference voltages for the surviving component {0, 1, 2}
        self.Vref = self.ref.ac_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                   self.max_it, self.tol)
        assert self.Vref.size == 3, "reference powerflow did not converge"

    def _run(self, handle_disconnected):
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # disconnect line 2 -> isolates bus 3
        SA.handle_disconnected_grid = handle_disconnected
        SA.compute(1. * self.V0, self.max_it, self.tol)
        return SA.get_voltages().copy()

    def test_default_is_off(self):
        SA = ContingencyAnalysisCPP(self.grid)
        assert SA.handle_disconnected_grid is False

    def test_flag_off_skips_split(self):
        # legacy behaviour: a splitting contingency is not simulated -> all 0.
        v = self._run(False)
        assert np.max(np.abs(v[0])) == 0., f"split contingency should be skipped, got {v[0]}"

    def test_flag_on_solves_main_component(self):
        v = self._run(True)
        # the surviving component {0, 1, 2} matches the reference powerflow
        assert np.max(np.abs(v[0, :3] - self.Vref)) <= 1e-6, \
            f"main component mismatch: {np.max(np.abs(v[0, :3] - self.Vref)):.2e}"
        # the isolated bus 3 is masked -> reported as exactly 0
        assert v[0, 3] == 0., f"masked bus should be 0, got {v[0, 3]}"

    def test_masked_setter_validates_bool(self):
        SA = ContingencyAnalysisCPP(self.grid)
        # the C++ property only accepts a bool; a clean round-trip is enough here
        SA.handle_disconnected_grid = True
        assert SA.handle_disconnected_grid is True
        SA.handle_disconnected_grid = False
        assert SA.handle_disconnected_grid is False

    def test_error_on_non_nr_ac_algorithm(self):
        # the mode requires an NR algorithm: an AC non-NR solver must be rejected
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.GaussSeidel)
        SA.add_n1(2)
        SA.handle_disconnected_grid = True
        with self.assertRaises(RuntimeError):
            SA.compute(1. * self.V0, self.max_it, self.tol)

    def test_dc_flag_off_skips_split(self):
        # legacy DC behaviour: a split contingency diverges (no slack reference in the
        # stranded island) and is not stored -> all 0.
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.DC_SparseLU)
        SA.add_n1(2)
        SA.handle_disconnected_grid = False
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()
        assert np.max(np.abs(v[0])) == 0., f"DC split contingency should be skipped, got {v[0]}"

    def test_dc_flag_on_solves_main_component(self):
        # with the flag ON, DC solves the largest component (masking the island).
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.DC_SparseLU)
        SA.add_n1(2)
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()[0]
        # the surviving component {0,1,2} matches a DC powerflow on the reference grid
        Vref_dc = self.ref.dc_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                 self.max_it, self.tol)
        assert np.max(np.abs(np.angle(v[:3]) - np.angle(Vref_dc))) <= 1e-6, \
            f"DC main component mismatch: {np.max(np.abs(np.angle(v[:3]) - np.angle(Vref_dc))):.2e}"
        # the isolated bus 3 is masked -> reported as exactly 0
        assert v[3] == 0., f"masked bus should be 0, got {v[3]}"

    def test_is_grid_connected_no_segfault(self):
        # regression: is_grid_connected_after_contingency() used to segfault because
        # it relied on the (empty) grid-model Ybus. It must now work both standalone
        # (before compute) and after compute, and agree.
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # splits the grid (isolates bus 3)
        SA.add_n1(0)  # also splits (isolates {1, 2, 3} from the slack)
        standalone = np.asarray(SA.is_grid_connected_after_contingency())
        assert standalone.shape == (2,)
        assert np.all(standalone == 0), f"both radial cuts disconnect the grid, got {standalone}"
        SA.compute(1. * self.V0, self.max_it, self.tol)
        after = np.asarray(SA.is_grid_connected_after_contingency())
        assert np.array_equal(standalone, after)

    def test_is_grid_connected_dc(self):
        # in DC the connectivity is now reported off the real Bbus (it used to always
        # claim "connected"); both radial cuts disconnect the grid.
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.DC_SparseLU)
        SA.add_n1(2)
        SA.add_n1(0)
        standalone = np.asarray(SA.is_grid_connected_after_contingency())
        assert np.all(standalone == 0), f"both radial cuts disconnect the grid, got {standalone}"
        SA.compute(1. * self.V0, self.max_it, self.tol)
        after = np.asarray(SA.is_grid_connected_after_contingency())
        assert np.array_equal(standalone, after)


class TestContingencySplitMultiSlack(unittest.TestCase):
    """Two slacks at the two ends of a line; a contingency strands one of them.

    Topology: ext_grid@0 - 1 - 2 - 3 - ext_grid@4. Disconnecting line 2 (2-3)
    splits it into {0,1,2} (with slack 0) and {3,4} (with slack 4). The largest
    component {0,1,2} is solved; the stranded slack 4 has its weight zeroed.
    """
    def setUp(self):
        self.max_it = 30
        self.tol = 1e-8
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(self._net(with_island=True))
            self.ref = init_from_pandapower(self._net(with_island=False))
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        self.Vref = self.ref.ac_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                   self.max_it, self.tol)
        assert self.Vref.size == 3

    @staticmethod
    def _net(with_island):
        net = pp.create_empty_network(sn_mva=1.)
        nb = 5 if with_island else 3
        for _ in range(nb):
            pp.create_bus(net, vn_kv=20.)
        pp.create_ext_grid(net, 0, vm_pu=1.0)
        _line(net, 0, 1)   # line 0
        _line(net, 1, 2)   # line 1
        pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
        pp.create_load(net, 2, p_mw=1.0, q_mvar=0.2)
        if with_island:
            _line(net, 2, 3)   # line 2 (splits the grid)
            _line(net, 3, 4)   # line 3
            pp.create_ext_grid(net, 4, vm_pu=1.0)   # second slack, stranded by the split
            pp.create_load(net, 3, p_mw=1.0, q_mvar=0.2)
        return net

    def test_stranded_slack_is_masked(self):
        SA = ContingencyAnalysisCPP(self.grid)
        SA.add_n1(2)  # isolates {3, 4} (which holds the second slack)
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()
        # surviving component {0,1,2} solved with slack 0 absorbing everything;
        # the reference grid keeps only slack 0, so the voltages must match.
        assert np.max(np.abs(v[0, :3] - self.Vref)) <= 1e-6, \
            f"main component mismatch: {np.max(np.abs(v[0, :3] - self.Vref)):.2e}"
        # the stranded island {3, 4} is masked -> reported as 0
        assert np.max(np.abs(v[0, 3:5])) == 0., f"masked island should be 0, got {v[0, 3:5]}"

    def test_stranded_slack_is_masked_dc(self):
        # same in DC: the stranded second slack's weight is zeroed, the live slack 0
        # absorbs the imbalance and the island is reported as 0.
        SA = ContingencyAnalysisCPP(self.grid)
        SA.change_algorithm(AlgorithmType.DC_SparseLU)
        SA.add_n1(2)
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()
        Vref_dc = self.ref.dc_pf(np.ones(self.ref.get_bus_vn_kv().shape[0], dtype=complex),
                                 self.max_it, self.tol)
        assert np.max(np.abs(np.angle(v[0, :3]) - np.angle(Vref_dc))) <= 1e-6, \
            f"DC main component mismatch: {np.max(np.abs(np.angle(v[0, :3]) - np.angle(Vref_dc))):.2e}"
        assert np.max(np.abs(v[0, 3:5])) == 0., f"masked island should be 0, got {v[0, 3:5]}"


class TestContingencySplitCase14(unittest.TestCase):
    """Integration test on l2rpn_case14_sandbox: enabling the mode must leave the
    connected contingencies bit-identical and additionally solve at least one of
    the contingencies that split the grid (skipped when the mode is off)."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        self.nb_sub = self.env.n_sub

    def tearDown(self):
        self.env.close()

    def _voltages(self, handle_disconnected, algorithm=None):
        SA = ContingencyAnalysisCPP(self.env.backend._grid)
        if algorithm is not None:
            SA.change_algorithm(algorithm)
        SA.add_all_n1()
        SA.handle_disconnected_grid = handle_disconnected
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        return SA.get_voltages().copy()

    def test_regression_and_new_contingencies(self):
        voff = self._voltages(False)
        von = self._voltages(True)
        nb_sub = self.nb_sub

        # contingencies solved when the mode is OFF (connected, converged) must be
        # bit-identical when the mode is ON (their masked set is empty)
        off_solved = np.array([np.any(np.abs(voff[i, :nb_sub]) > 1e-9) for i in range(voff.shape[0])])
        assert off_solved.any(), "sanity: some contingencies should be solved with the mode off"
        assert np.max(np.abs(von[off_solved] - voff[off_solved])) == 0., \
            "connected contingencies must be unchanged when the mode is enabled"

        # at least one contingency skipped when OFF is now solved on its largest
        # connected component when ON
        newly_solved = [i for i in range(voff.shape[0])
                        if (not off_solved[i]) and np.any(np.abs(von[i, :nb_sub]) > 1e-9)]
        assert len(newly_solved) >= 1, "the mode should solve at least one split contingency"

        for i in newly_solved:
            vm = np.abs(von[i, :nb_sub])
            live = vm[vm > 1e-9]
            assert np.all(np.isfinite(von[i, :nb_sub])), f"cont {i}: non-finite voltage"
            assert live.min() > 0.5 and live.max() < 1.5, \
                f"cont {i}: unrealistic voltages on the live component [{live.min()}, {live.max()}]"
            assert np.sum(vm <= 1e-9) >= 1, f"cont {i}: expected at least one masked (0) bus"

    def test_regression_and_new_contingencies_dc(self):
        # same as above but with the DC solver
        voff = self._voltages(False, algorithm=AlgorithmType.DC_SparseLU)
        von = self._voltages(True, algorithm=AlgorithmType.DC_SparseLU)
        nb_sub = self.nb_sub

        # contingencies solved when the mode is OFF must be bit-identical when ON
        off_solved = np.array([np.any(np.abs(voff[i, :nb_sub]) > 1e-9) for i in range(voff.shape[0])])
        assert off_solved.any(), "sanity: some DC contingencies should be solved with the mode off"
        assert np.max(np.abs(von[off_solved] - voff[off_solved])) == 0., \
            "connected DC contingencies must be unchanged when the mode is enabled"

        # at least one split contingency (skipped when OFF) is now solved on its
        # largest connected component
        newly_solved = [i for i in range(voff.shape[0])
                        if (not off_solved[i]) and np.any(np.abs(von[i, :nb_sub]) > 1e-9)]
        assert len(newly_solved) >= 1, "the DC mode should solve at least one split contingency"
        for i in newly_solved:
            assert np.all(np.isfinite(von[i, :nb_sub])), f"DC cont {i}: non-finite voltage"
            assert np.sum(np.abs(von[i, :nb_sub]) <= 1e-9) >= 1, \
                f"DC cont {i}: expected at least one masked (0) bus"


class TestContingencySplitReferenceSlack(unittest.TestCase):
    """The reference slack of the "handle disconnected grid" mode, on a grid whose
    gridmodel and solver bus numberings differ.

    Topology: buses 0, 1, 2 out of service (so solver id = gridmodel id - 3), a meshed
    core 3..7 with the ext_grid on bus 5 (solver 2) and a generator on bus 3, and a leaf
    bus 8 (solver 5) with a generator, connected by line 7 only. Every generator takes
    part in the distributed slack.

    The reference used to be chosen mixing both numberings: gridmodel bus 5 (a slack
    bus) was checked against the weight and the strand status of SOLVER bus 5 (the
    leaf), so it was picked, and line 7 (which strands the leaf) was skipped. Its pick
    did not even reach the solver, which kept its own reference.
    """
    MAX_IT = 30
    TOL = 1e-10
    LEAF_BUS = 8
    BIG_SLACK_BUS = 3
    LEAF_LINE = 7
    MESH_LINE = 0  # 3-4: the grid stays in one piece

    @staticmethod
    def _net():
        net = pp.create_empty_network(sn_mva=100.)
        for bus_id in range(9):
            pp.create_bus(net, vn_kv=20., in_service=bus_id >= 3)
        for f, t in [(3, 4), (4, 5), (5, 6), (6, 3), (4, 6), (6, 7), (7, 3), (7, 8)]:
            _line(net, f, t)
        pp.create_ext_grid(net, 5, vm_pu=1.0)
        pp.create_gen(net, 3, p_mw=1.0, vm_pu=1.0)
        pp.create_gen(net, 8, p_mw=2.0, vm_pu=1.0)
        for bus_id in (4, 6, 7):
            pp.create_load(net, bus_id, p_mw=1.5, q_mvar=0.3)
        return net

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(self._net())
        # the generator of bus 3 has the largest weight: it is the automatic choice
        # (stranded by no contingency, then the largest weight)
        for gen in self.grid.get_generators():
            self.grid.add_gen_slackbus(gen.id, 2. if gen.bus_id == self.BIG_SLACK_BUS else 1.)
        self.nb_bus = self.grid.get_bus_vn_kv().shape[0]
        self.V0 = np.ones(self.nb_bus, dtype=complex)
        assert self.grid.ac_pf(1. * self.V0, self.MAX_IT, self.TOL).shape[0] > 0
        # the numbering this test is about: if it changes, the test must be rewritten
        me2s = np.asarray(self.grid.id_me_to_ac_solver())
        s2me = np.asarray(self.grid.id_ac_solver_to_me())
        assert not np.array_equal(s2me, np.arange(s2me.shape[0])), "numberings should differ"
        slack_me = set(np.asarray(self.grid.get_slack_ids()).tolist())
        assert slack_me == {3, 5, self.LEAF_BUS}, f"unexpected slack buses {slack_me}"
        assert me2s[self.LEAF_BUS] == 5, "the gridmodel id of a slack bus (5) should be the solver id of the leaf"
        self.v_angles = np.linspace(0.01, 0.09, self.nb_bus)  # a distinct angle per bus
        self.V_ang = np.exp(1j * self.v_angles)

    def _analysis(self, grid=None, conts=(LEAF_LINE, MESH_LINE), algorithm=None):
        SA = ContingencyAnalysisCPP(self.grid if grid is None else grid, True)
        if algorithm is not None:
            SA.change_algorithm(algorithm)
        SA.handle_disconnected_grid = True
        SA.redistribute_slack = True
        SA.add_multiple_n1(list(conts))
        return SA

    @staticmethod
    def _not_simulated(SA, row):
        return any(v.violation_type == LimitViolationType.NOT_SIMULATED for v in SA.get_violations()[row])

    def _row_of(self, SA, line_id):
        rows = [i for i, cont in enumerate(SA.my_defaults()) if list(cont) == [line_id]]
        assert len(rows) == 1
        return rows[0]

    def test_island_with_slack_gen_is_solved(self):
        SA = self._analysis()
        SA.compute(1. * self.V0, self.MAX_IT, self.TOL)
        row = self._row_of(SA, self.LEAF_LINE)
        assert SA.converged()[row], "the contingency islanding the leaf slack generator should be solved"
        assert not self._not_simulated(SA, row)
        v_batch = SA.get_voltages()[row]
        assert v_batch[self.LEAF_BUS] == 0., "the islanded leaf should be reported as 0"

        # same as one contingency at a time
        one = self.grid.copy()
        one.deactivate_powerline(self.LEAF_LINE)
        one.consider_only_main_component(True)
        v_one = one.ac_pf(1. * self.V0, self.MAX_IT, self.TOL)
        assert v_one.shape[0] > 0
        live = np.arange(3, self.LEAF_BUS)
        # the two may use another reference bus: compare the angles relative to bus 3
        rel_batch = v_batch[live] * np.conj(v_batch[3]) / np.abs(v_batch[3])
        rel_one = v_one[live] * np.conj(v_one[3]) / np.abs(v_one[3])
        assert np.max(np.abs(rel_batch - rel_one)) <= 1e-6, \
            f"batch / one at a time mismatch: {np.max(np.abs(rel_batch - rel_one)):.2e}"

    def test_island_with_slack_gen_is_solved_dc(self):
        SA = self._analysis(algorithm=AlgorithmType.DC_SparseLU)
        SA.compute(1. * self.V0, self.MAX_IT, self.TOL)
        row = self._row_of(SA, self.LEAF_LINE)
        assert SA.converged()[row], "DC: the contingency islanding the leaf slack generator should be solved"
        assert not self._not_simulated(SA, row)
        assert SA.get_voltages()[row][self.LEAF_BUS] == 0.

    def test_pick_reference_slack_is_slack_and_not_stranded(self):
        SA = self._analysis()
        ref = SA.pick_reference_slack()
        assert ref in np.asarray(self.grid.get_slack_ids()).tolist(), f"{ref} is not a slack bus"
        assert ref != self.LEAF_BUS, "the leaf is stranded by a contingency, it should not be picked"
        assert ref == self.BIG_SLACK_BUS, \
            f"expected bus {self.BIG_SLACK_BUS} (not stranded, largest weight), got {ref}"

    def test_reference_used_by_solver(self):
        # the NR keeps the angle of its reference bus at its starting value
        SA = self._analysis()
        SA.init_from_n_powerflow = False
        ref = SA.pick_reference_slack()
        SA.compute(1. * self.V_ang, self.MAX_IT, self.TOL)
        for row in range(len(SA.my_defaults())):
            assert SA.converged()[row]
            v = SA.get_voltages()[row]
            assert abs(np.angle(v[ref]) - self.v_angles[ref]) <= 1e-12, \
                f"row {row}: bus {ref} (pick_reference_slack) is not the reference of the solver"

    def test_forced_reference_is_respected(self):
        self.grid.set_reference_slack_bus(self.LEAF_BUS)
        SA = self._analysis()
        SA.init_from_n_powerflow = False
        # pick_reference_slack stays a suggestion: the automatic choice
        assert SA.pick_reference_slack() != self.LEAF_BUS
        SA.compute(1. * self.V_ang, self.MAX_IT, self.TOL)
        # the forced reference is kept for the whole batch: the contingency stranding it is skipped
        row_leaf = self._row_of(SA, self.LEAF_LINE)
        assert not SA.converged()[row_leaf]
        assert self._not_simulated(SA, row_leaf)
        # ... and it is the reference of the other ones
        row_mesh = self._row_of(SA, self.MESH_LINE)
        assert SA.converged()[row_mesh]
        v = SA.get_voltages()[row_mesh]
        assert abs(np.angle(v[self.LEAF_BUS]) - self.v_angles[self.LEAF_BUS]) <= 1e-12, \
            "the forced reference should be the reference of the solver"

    def test_copy_keeps_forced_reference(self):
        self.grid.set_reference_slack_bus(self.LEAF_BUS)
        assert self.grid.copy().get_reference_slack_bus() == self.LEAF_BUS


if __name__ == "__main__":
    unittest.main()
