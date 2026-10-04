# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The physical check of an idle SVC carrying a standby automaton (``LSGrid.set_svc_standby``):
a flagged, non-regulating SVC whose regulated bus leaves the automaton's voltage thresholds --
which OpenLoadFlow's MonitoringVoltageOuterLoop would switch to voltage control -- is reported
as ``LOW_VOLTAGE_SVC_STANDBY`` / ``HIGH_VOLTAGE_SVC_STANDBY`` on the SVC, in kV, by a single
solve and by every batch alike; and ``bake_outer_loops`` / ``init_from_pypowsybl(can_be_pv=...)``
flag the SVCs a bake left idle."""

import copy
import os
import pickle
import tempfile
import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import (
    LSGrid,
    LimitViolationType,
    ViolationCategory,
    ViolationElementType,
)
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP
from lightsim2grid.timeSerie import TimeSeriesCPP

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    import pypowsybl.report as rp
    from lightsim2grid.network import bake_outer_loops, init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

VN_KV = 138.
STANDBY_TYPES = (LimitViolationType.LOW_VOLTAGE_SVC_STANDBY,
                 LimitViolationType.HIGH_VOLTAGE_SVC_STANDBY)
# SvcContainer.RegulationMode
OFF_MODE, VOLTAGE_MODE, REACTIVE_POWER_MODE = 0, 1, 2


def _standby(viols):
    return [v for v in viols if v.violation_type in STANDBY_TYPES]


def _feeder(svc_bus=2, svc_reg_bus=None, svc_mode=REACTIVE_POWER_MODE, nb_extra_leaf=0):
    """the 4-bus radial feeder 0-1-2-3 (80 MW / 60 MVAr load on bus 3), gen 0 the PV slack on
    bus 0, and one SVC (fixed Q = 0 by default) on `svc_bus`. `nb_extra_leaf` adds a leaf bus 4
    hanging off bus 1 (line 3)."""
    nb_bus = 4 + nb_extra_leaf
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(nb_bus, 1, np.full(nb_bus, VN_KV), 0, 0)
    fr = [0, 1, 2] + ([1] if nb_extra_leaf else [])
    to = [1, 2, 3] + ([4] if nb_extra_leaf else [])
    nb_line = len(fr)
    grid.init_powerlines(np.full(nb_line, 0.01), np.full(nb_line, 0.1), np.zeros(nb_line, dtype=complex),
                         np.array(fr), np.array(to))
    grid.init_loads(np.array([80.]), np.array([60.]), np.array([3]))
    grid.init_generators_full(np.array([0.]), np.array([1.02]), np.array([0.]), [True],
                              np.array([-1e3]), np.array([1e3]), np.array([0]))
    grid.set_gen_names(["slack"])
    grid.add_gen_slackbus(0, 1.)
    reg_bus = svc_bus if svc_reg_bus is None else svc_reg_bus
    grid.init_svcs([svc_mode], np.array([1.0]), np.array([0.]), np.array([0.]),
                   np.array([-1.]), np.array([1.]), np.array([reg_bus], dtype=np.int32),
                   np.array([svc_bus], dtype=np.int32))
    grid.set_svc_names(["svc"])
    grid.tell_solver_need_reset()
    return grid


def _solve(grid):
    V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
    assert V.shape[0] > 0
    return V


class TestSvcStandbyFlag(unittest.TestCase):
    """`set_svc_standby` itself: validation, what `SvcInfo` shows, and that it survives a
    copy, a pickle and a binary save."""

    def test_default_and_set(self):
        grid = _feeder()
        svc = grid.get_svcs()[0]
        self.assertFalse(svc.standby)
        self.assertTrue(np.isnan(svc.standby_low_vm_pu))
        self.assertTrue(np.isnan(svc.standby_high_vm_pu))
        grid.set_svc_standby(np.array([True]), np.array([0.95]), np.array([1.05]))
        svc = grid.get_svcs()[0]
        self.assertTrue(svc.standby)
        self.assertEqual(svc.standby_low_vm_pu, 0.95)
        self.assertEqual(svc.standby_high_vm_pu, 1.05)
        # an unflagged SVC keeps no threshold
        grid.set_svc_standby(np.array([False]), np.array([0.95]), np.array([1.05]))
        self.assertTrue(np.isnan(grid.get_svcs()[0].standby_low_vm_pu))

    def test_refused(self):
        grid = _feeder()
        with self.assertRaises(RuntimeError):  # wrong size
            grid.set_svc_standby(np.array([True, True]), np.array([0.95, 0.95]), np.array([1.05, 1.05]))
        with self.assertRaises(RuntimeError):  # low >= high
            grid.set_svc_standby(np.array([True]), np.array([1.05]), np.array([0.95]))
        with self.assertRaises(RuntimeError):  # no threshold
            grid.set_svc_standby(np.array([True]), np.array([np.nan]), np.array([1.05]))
        self.assertFalse(grid.get_svcs()[0].standby)  # nothing half-applied

    def test_kept_by_copy_pickle_and_binary(self):
        grid = _feeder()
        grid.set_svc_standby(np.array([True]), np.array([0.95]), np.array([1.05]))
        for other in (grid.copy(), copy.deepcopy(grid), pickle.loads(pickle.dumps(grid))):
            svc = other.get_svcs()[0]
            self.assertTrue(svc.standby)
            self.assertEqual((svc.standby_low_vm_pu, svc.standby_high_vm_pu), (0.95, 1.05))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            svc = LSGrid.load_binary(path).get_svcs()[0]
        self.assertTrue(svc.standby)
        self.assertEqual((svc.standby_low_vm_pu, svc.standby_high_vm_pu), (0.95, 1.05))
        # a fresh init of the SVCs describes other SVCs: nothing flagged
        grid.init_svcs([REACTIVE_POWER_MODE], np.array([1.0]), np.array([0.]), np.array([0.]),
                       np.array([-1.]), np.array([1.]), np.array([2], dtype=np.int32),
                       np.array([2], dtype=np.int32))
        self.assertFalse(grid.get_svcs()[0].standby)


class TestSvcStandbyCheck(unittest.TestCase):
    """The check on a hand-made grid: reported in kV, on the SVC, only when flagged, idle and
    outside the thresholds -- and the same by a single solve and by the batches."""

    @staticmethod
    def _grid_outside(high=True, **kwargs):
        """the feeder, its SVC flagged with thresholds the solved voltage of bus 2 is outside
        (above the high one when `high`, below the low one otherwise)"""
        grid = _feeder(**kwargs)
        vm = abs(_solve(grid)[2])
        if high:
            low, hi = vm - 0.05, vm - 0.01
        else:
            low, hi = vm + 0.01, vm + 0.05
        grid.set_svc_standby(np.array([True]), np.array([low]), np.array([hi]))
        return grid, vm, (low, hi)

    def test_single_solve_high(self):
        grid, vm, (low, hi) = self._grid_outside(high=True)
        V = _solve(grid)
        viols = _standby(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        v = viols[0]
        self.assertEqual(v.element_type, ViolationElementType.SVC)
        self.assertEqual(v.element_id, 0)
        self.assertEqual(v.violation_type, LimitViolationType.HIGH_VOLTAGE_SVC_STANDBY)
        self.assertEqual(v.category, ViolationCategory.PHYSICAL)
        self.assertEqual(v.name, "svc")
        self.assertAlmostEqual(v.value, abs(V[2]) * VN_KV, places=6)
        self.assertAlmostEqual(v.limit, hi * VN_KV, places=9)
        self.assertGreater(v.value, v.limit)

    def test_single_solve_low(self):
        grid, vm, (low, hi) = self._grid_outside(high=False)
        _solve(grid)
        viols = _standby(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].violation_type, LimitViolationType.LOW_VOLTAGE_SVC_STANDBY)
        self.assertAlmostEqual(viols[0].limit, low * VN_KV, places=9)
        self.assertLess(viols[0].value, viols[0].limit)

    def test_nothing_reported(self):
        # inside the thresholds
        grid = _feeder()
        vm = abs(_solve(grid)[2])
        grid.set_svc_standby(np.array([True]), np.array([vm - 0.01]), np.array([vm + 0.01]))
        _solve(grid)
        self.assertEqual(_standby(grid.get_physical_violations(True, 0., 0.)), [])
        # outside, but not flagged
        grid = _feeder()
        _solve(grid)
        self.assertEqual(_standby(grid.get_physical_violations(True, 0., 0.)), [])
        # outside, but within the tolerance
        grid, vm, _ = self._grid_outside(high=True)
        _solve(grid)
        self.assertEqual(_standby(grid.get_physical_violations(True, 0., 0.02)), [])
        # outside, but disconnected
        grid, vm, _ = self._grid_outside(high=True)
        grid.deactivate_svc(0)
        _solve(grid)
        self.assertEqual(_standby(grid.get_physical_violations(True, 0., 0.)), [])

    def test_regulating_svc_is_already_on(self):
        grid = _feeder(svc_mode=VOLTAGE_MODE)
        grid.set_svc_standby(np.array([True]), np.array([1.1]), np.array([1.2]))
        _solve(grid)
        self.assertEqual(_standby(grid.get_physical_violations(True, 0., 0.)), [])

    def test_remote_svc_checks_its_regulated_bus(self):
        # the SVC sits on bus 2 and regulates bus 1: bus 1 is what the automaton monitors
        grid = _feeder(svc_bus=2, svc_reg_bus=1)
        V = _solve(grid)
        vm1, vm2 = abs(V[1]), abs(V[2])
        self.assertGreater(vm1, vm2)
        # its own bus (2) inside the thresholds, the regulated one (1) above the high one:
        # reported, with the voltage of bus 1
        mid = 0.5 * (vm1 + vm2)
        grid.set_svc_standby(np.array([True]), np.array([vm2 - 0.1]), np.array([mid]))
        V = _solve(grid)
        viols = _standby(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertAlmostEqual(viols[0].value, abs(V[1]) * VN_KV, places=6)

    def test_batches_match_the_single_solve(self):
        grid, vm, _ = self._grid_outside(high=True)
        _solve(grid)
        ref = _standby(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(ref), 1)

        ts = TimeSeriesCPP(grid)
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.physical_violation_tol_vm_pu = 0.
        ts.modify_gen_p(np.array([[g.target_p_mw for g in grid.get_generators()]]))
        ts.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert ts.converged_mask()[0]
        row = _standby(ts.get_physical_violations()[0])
        self.assertEqual(len(row), 1)
        self.assertEqual(row[0].element_type, ViolationElementType.SVC)
        self.assertAlmostEqual(row[0].value, ref[0].value, places=6)
        self.assertAlmostEqual(row[0].limit, ref[0].limit, places=9)

        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.add_n1(2)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        n_case = _standby(ca.get_physical_violations_n())
        self.assertEqual(len(n_case), 1)
        self.assertAlmostEqual(n_case[0].value, ref[0].value, places=6)

    def test_stranded_svc_reports_nothing(self):
        # the SVC on the leaf bus 4 (off bus 1), regulating bus 1: taking line 3 out strands
        # it while the bus it regulates stays in the main component
        grid = _feeder(svc_bus=4, svc_reg_bus=1, nb_extra_leaf=1)
        vm = abs(_solve(grid)[1])
        grid.set_svc_standby(np.array([True]), np.array([vm - 0.05]), np.array([vm - 0.01]))
        _solve(grid)
        self.assertEqual(len(_standby(grid.get_physical_violations(True, 0., 0.))), 1,
                         "sanity: the N state reports the switch")

        ca = ContingencyAnalysisCPP(grid, True)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.handle_disconnected_grid = True
        ca.add_n1(3)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert list(ca.converged()) == [True]
        self.assertEqual(len(_standby(ca.get_physical_violations_n())), 1)
        self.assertEqual(_standby(ca.get_physical_violations()[0]), [],
                         "a stranded SVC should not be reported as switched on")


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestSvcStandbyFromPypowsybl(unittest.TestCase):
    """`bake_outer_loops` returns the standby SVC it left idle, `init_from_pypowsybl(can_be_pv=...)`
    flags it with the automaton's thresholds, and a change of the grid that makes OLF switch it
    on is reported by lightsim2grid on the baked grid."""

    # thresholds bracketing the voltage of the SVC's bus when the SVC idles (b0 = 0)
    LOW_KV, HIGH_KV = 398., 401.

    def _network(self, regulation_mode="VOLTAGE"):
        n = pp.network.create_four_substations_node_breaker_network()
        if regulation_mode != "VOLTAGE":
            n.update_static_var_compensators(id="SVC", target_q=0.)
            n.update_static_var_compensators(id="SVC", regulation_mode=regulation_mode)
        n.create_extensions("standbyAutomaton", id="SVC", standby=True, b0=0.,
                            low_voltage_threshold=self.LOW_KV, low_voltage_setpoint=self.LOW_KV + 1.,
                            high_voltage_threshold=self.HIGH_KV, high_voltage_setpoint=self.HIGH_KV - 1.)
        return n

    @staticmethod
    def _svc_bus_vn(n):
        bus = n.get_static_var_compensators().loc["SVC", "voltage_level_id"]
        return float(n.get_voltage_levels().loc[bus, "nominal_v"])

    def _baked(self):
        n = self._network()
        lf.run_ac(n)
        self.assertEqual(n.get_static_var_compensators().loc["SVC", "q"], 0.,
                         "sanity: OLF keeps the SVC idle")
        return n, bake_outer_loops(n)

    def test_bake_returns_the_idle_svc_and_init_flags_it(self):
        n, pinned = self._baked()
        self.assertIn("SVC", set(pinned))
        self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], "REACTIVE_POWER")
        self.assertEqual(len(bake_outer_loops(n)), 0)  # already baked
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=pinned)
        svc = grid.get_svcs()[0]
        vn = self._svc_bus_vn(n)
        self.assertTrue(svc.standby)
        self.assertAlmostEqual(svc.standby_low_vm_pu, self.LOW_KV / vn, places=12)
        self.assertAlmostEqual(svc.standby_high_vm_pu, self.HIGH_KV / vn, places=12)
        # not flagged without `can_be_pv`
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False)
        self.assertFalse(grid.get_svcs()[0].standby)
        # an unknown id is still refused
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                can_be_pv=list(pinned) + ["NOT-AN-ELEMENT"])

    def test_svc_without_automaton_is_flagged_as_frozen(self):
        # an SVC id of `can_be_pv` with no standby automaton is one frozen at a limit
        n = pp.network.create_four_substations_node_breaker_network()
        svc = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False,
                                  can_be_pv=["SVC"]).get_svcs()[0]
        self.assertFalse(svc.standby)
        self.assertTrue(svc.can_be_pv)

    def test_non_regulating_svc_ignores_its_automaton(self):
        # OLF only arms the automaton of an SVC regulating voltage
        n = self._network(regulation_mode="REACTIVE_POWER")
        lf.run_ac(n)
        self.assertNotIn("SVC", set(bake_outer_loops(n)))

    def _reported(self, n, pinned):
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=pinned)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        return _standby(grid.get_physical_violations())

    def test_switch_reported_where_olf_switches_it_on(self):
        n, pinned = self._baked()
        self.assertEqual(self._reported(n, pinned), [])  # the reference state: nothing

        for q0, expected in ((-100., LimitViolationType.HIGH_VOLTAGE_SVC_STANDBY),
                             (100., LimitViolationType.LOW_VOLTAGE_SVC_STANDBY)):
            with self.subTest(q0=q0):
                # OLF, outer loops on, on the raw network: the automaton switches the SVC on
                raw = self._network()
                raw.update_loads(id="LD6", q0=q0)
                lf.run_ac(raw)
                self.assertNotEqual(raw.get_static_var_compensators().loc["SVC", "q"], 0.)
                # lightsim2grid on the baked network cannot, and reports it
                baked, pinned = self._baked()
                baked.update_loads(id="LD6", q0=q0)
                viols = self._reported(baked, pinned)
                self.assertEqual(len(viols), 1)
                self.assertEqual(viols[0].violation_type, expected)
                self.assertEqual(viols[0].name, "SVC")
                limit = self.HIGH_KV if expected == LimitViolationType.HIGH_VOLTAGE_SVC_STANDBY else self.LOW_KV
                self.assertAlmostEqual(viols[0].limit, limit, places=6)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestSvcStandbyB0Range(unittest.TestCase):
    """An SVC carrying a standby automaton, in standby or not: OpenLoadFlow models its b0 as
    a fixed susceptance apart from the SVC, whose own susceptance stays in [b_min, b_max] --
    so the total output (what pypowsybl reports and lightsim2grid models) ranges over
    [b_min + b0, b_max + b0]."""
    V_KV = 400.

    def _net(self, b0_mvar):
        # the four substations SVC holds its bus producing a bit more than 12.5 MVAr
        n = pp.network.create_four_substations_node_breaker_network()
        n.update_static_var_compensators(id="SVC", b_max=15. / self.V_KV ** 2, b_min=-15. / self.V_KV ** 2)
        n.create_extensions("standbyAutomaton", id="SVC", b0=b0_mvar / self.V_KV ** 2, standby=False,
                            low_voltage_threshold=380., low_voltage_setpoint=390.,
                            high_voltage_threshold=420., high_voltage_setpoint=410.)
        return n

    def _olf_switches(self, n):
        report = rp.ReportNode()
        lf.run_ac(n, lf.Parameters(), report_node=report)
        return any("Switch bus" in line for line in str(report).splitlines())

    def _ls_high_q(self, n):
        grid = init_from_pypowsybl(n, sort_index=True)
        assert grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10).shape[0] > 0
        return [el for el in grid.get_physical_violations(True, 0., 0.)
                if el.violation_type == LimitViolationType.HIGH_Q]

    def test_reported_where_olf_switches(self):
        # an inductive b0 asks the SVC part for more than its b_max: OLF switches it; a
        # capacitive one leaves it room
        for b0_mvar, switched in [(0., False), (-5., True), (5., False)]:
            self.assertEqual(self._olf_switches(self._net(b0_mvar)), switched, f"b0 {b0_mvar}")
            viols = self._ls_high_q(self._net(b0_mvar))
            self.assertEqual(len(viols) == 1, switched, f"b0 {b0_mvar}")
            if switched:
                self.assertAlmostEqual(viols[0].limit, 15. + b0_mvar, places=3)

    def test_b0_only_shifts_a_voltage_controller(self):
        # OLF reads b0 in LfStaticVarCompensatorImpl.setupVoltageControl only: an SVC in
        # REACTIVE_POWER mode, or not regulating, keeps its own [b_min, b_max]
        for mode, regulating, shifted in [("VOLTAGE", True, True),
                                          ("REACTIVE_POWER", True, False),
                                          ("VOLTAGE", False, False)]:
            n = self._net(-5.)
            n.update_static_var_compensators(id="SVC", target_q=0.)
            n.update_static_var_compensators(id="SVC", regulation_mode=mode, regulating=regulating)
            grid = init_from_pypowsybl(n, sort_index=True)
            svc = grid.get_svcs()[0]
            # b in S at 400 kV -> pu: b_max (S) * V^2 / sn_mva is the range in MVAr / sn_mva
            expected = (10. if shifted else 15.) / grid.get_sn_mva()
            self.assertAlmostEqual(svc.b_max, expected, places=6, msg=f"{mode} regulating={regulating}")

    def test_bake_freezes_the_svc_at_its_shifted_limit(self):
        # OLF switched it at the shifted limit: the bake must see it there and freeze it
        n = self._net(-5.)
        lf.run_ac(n, lf.Parameters())
        v_ref = n.get_buses()["v_mag"].copy()
        bake_outer_loops(n)
        self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], "REACTIVE_POWER")
        from lightsim2grid.network import remove_outer_loops
        res = lf.run_ac(n, remove_outer_loops(lf.Parameters()))
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        self.assertLess((n.get_buses()["v_mag"] - v_ref).abs().max(), 1e-2)


if __name__ == "__main__":
    unittest.main()
