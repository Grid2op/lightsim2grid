# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The PQ -> PV release of an SVC an outer loop froze at a reactive limit
(``LSGrid.set_svc_can_be_pv``, the generators' ``can_be_pv``): a fixed-Q SVC at the absorbing
(resp. producing) end of its susceptance range whose regulated bus sits below (resp. above) its
target is reported as ``LOW_VOLTAGE_AT_MIN_Q`` / ``HIGH_VOLTAGE_AT_MAX_Q`` on the SVC, in kV, by
a single solve and by every batch alike; and ``bake_outer_loops`` /
``init_from_pypowsybl(can_be_pv=...)`` flag the SVCs a bake froze at a limit."""

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
    from lightsim2grid.network import bake_outer_loops, init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

VN_KV = 138.
RELEASE_TYPES = (LimitViolationType.LOW_VOLTAGE_AT_MIN_Q, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)
# SvcContainer.RegulationMode
VOLTAGE_MODE, REACTIVE_POWER_MODE = 1, 2
# the SVC's susceptance range (pu, base 100 MVA): it absorbs up to 5 MVAr, produces up to 50
B_MIN, B_MAX = -0.05, 0.5


def _release(viols, svc_id=0):
    return [v for v in viols if v.element_type == ViolationElementType.SVC
            and v.element_id == svc_id and v.violation_type in RELEASE_TYPES]


def _frozen_svc_grid(target_vm=1.0, at_min=True, flagged=True, mode=REACTIVE_POWER_MODE,
                     svc_bus=2, svc_reg_bus=None, leaf=False):
    """the 4-bus radial feeder 0-1-2-3 (80 MW / 60 MVAr load on bus 3), gen 0 the PV slack on
    bus 0, and an SVC on `svc_bus` frozen at the absorbing (`at_min`) or producing end of its
    range at `target_vm` -- the output it had while holding that target. `leaf` adds a leaf
    bus 4 off bus 1 (line 3)."""
    nb_bus = 5 if leaf else 4
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(nb_bus, 1, np.full(nb_bus, VN_KV), 0, 0)
    fr, to = [0, 1, 2] + ([1] if leaf else []), [1, 2, 3] + ([4] if leaf else [])
    nl = len(fr)
    grid.init_powerlines(np.full(nl, 0.01), np.full(nl, 0.1), np.zeros(nl, dtype=complex),
                         np.array(fr), np.array(to))
    grid.init_loads(np.array([80.]), np.array([60.]), np.array([3]))
    grid.init_generators_full(np.array([0.]), np.array([1.02]), np.array([0.]), [True],
                              np.array([-1e3]), np.array([1e3]), np.array([0]))
    grid.set_gen_names(["slack"])
    grid.add_gen_slackbus(0, 1.)
    q_limit = (B_MIN if at_min else B_MAX) * target_vm ** 2 * 100.
    reg = svc_bus if svc_reg_bus is None else svc_reg_bus
    grid.init_svcs([mode], np.array([target_vm]), np.array([q_limit]), np.array([0.]),
                   np.array([B_MIN]), np.array([B_MAX]), np.array([reg], dtype=np.int32),
                   np.array([svc_bus], dtype=np.int32))
    grid.set_svc_names(["frozen"])
    if flagged:
        grid.set_svc_can_be_pv(np.array([True]))
    grid.tell_solver_need_reset()
    return grid


def _solve(grid):
    V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
    assert V.shape[0] > 0
    return V


class TestSvcCanBePvFlag(unittest.TestCase):

    def test_default_set_and_kept(self):
        grid = _frozen_svc_grid(flagged=False)
        self.assertFalse(grid.get_svcs()[0].can_be_pv)
        grid.set_svc_can_be_pv(np.array([True]))
        self.assertTrue(grid.get_svcs()[0].can_be_pv)
        with self.assertRaises(RuntimeError):
            grid.set_svc_can_be_pv(np.array([True, True]))
        for other in (grid.copy(), copy.deepcopy(grid), pickle.loads(pickle.dumps(grid))):
            self.assertTrue(other.get_svcs()[0].can_be_pv)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            self.assertTrue(LSGrid.load_binary(path).get_svcs()[0].can_be_pv)


class TestSvcRelease(unittest.TestCase):

    def test_at_min_below_target(self):
        grid = _frozen_svc_grid(1.0, at_min=True)
        V = _solve(grid)
        self.assertLess(abs(V[2]), 1.0)
        viols = _release(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        v = viols[0]
        self.assertEqual(v.violation_type, LimitViolationType.LOW_VOLTAGE_AT_MIN_Q)
        self.assertEqual(v.category, ViolationCategory.CONTROL)
        self.assertEqual(v.name, "frozen")
        self.assertAlmostEqual(v.value, abs(V[2]) * VN_KV, places=6)
        self.assertAlmostEqual(v.limit, 1.0 * VN_KV, places=9)

    def test_at_max_above_target(self):
        grid = _frozen_svc_grid(0.7, at_min=False)
        V = _solve(grid)
        self.assertGreater(abs(V[2]), 0.7)
        viols = _release(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].violation_type, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)
        self.assertAlmostEqual(viols[0].limit, 0.7 * VN_KV, places=9)

    def test_nothing_reported(self):
        for kwargs, tol in [(dict(flagged=False), 0.),                # not flagged
                            (dict(target_vm=0.5, at_min=True), 0.),    # voltage above: not a release
                            (dict(mode=VOLTAGE_MODE), 0.),             # regulating already
                            (dict(), 1.)]:                             # within the tolerance
            with self.subTest(**kwargs, tol=tol):
                grid = _frozen_svc_grid(**kwargs)
                _solve(grid)
                self.assertEqual(_release(grid.get_physical_violations(True, 0., tol)), [])

    def test_batches_match_the_single_solve(self):
        grid = _frozen_svc_grid(1.0, at_min=True)
        _solve(grid)
        ref = _release(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(ref), 1)

        ts = TimeSeriesCPP(grid)
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.physical_violation_tol_vm_pu = 0.
        ts.modify_gen_p(np.array([[g.target_p_mw for g in grid.get_generators()]]))
        ts.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert ts.converged_mask()[0]
        row = _release(ts.get_physical_violations()[0])
        self.assertEqual(len(row), 1)
        self.assertAlmostEqual(row[0].value, ref[0].value, places=6)
        self.assertAlmostEqual(row[0].limit, ref[0].limit, places=9)

        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.add_n1(2)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        n_case = _release(ca.get_physical_violations_n())
        self.assertEqual(len(n_case), 1)
        self.assertAlmostEqual(n_case[0].value, ref[0].value, places=6)

    def test_stranded_svc_releases_nothing(self):
        # the SVC on the leaf bus 4 (off bus 1), regulating bus 1: taking line 3 out strands
        # it while the bus it regulates stays in the main component
        grid = _frozen_svc_grid(1.0, at_min=True, svc_bus=4, svc_reg_bus=1, leaf=True)
        _solve(grid)
        self.assertEqual(len(_release(grid.get_physical_violations(True, 0., 0.))), 1,
                         "sanity: the N state reports the release")
        ca = ContingencyAnalysisCPP(grid, True)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.handle_disconnected_grid = True
        ca.add_n1(3)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert list(ca.converged()) == [True]
        self.assertEqual(_release(ca.get_physical_violations()[0]), [])


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestSvcReleaseFromPypowsybl(unittest.TestCase):
    """`bake_outer_loops` returns an SVC it froze at a limit, `init_from_pypowsybl(can_be_pv=...)`
    flags it (not as a standby SVC), and a change of the grid that makes OLF put it back in
    voltage control is reported by lightsim2grid on the baked grid."""

    # the SVC of the four substations network holds its 400 kV bus producing a bit more than
    # 12.5 MVAr: with at most 10 it saturates and its bus sags below the target
    B_MAX_MVAR = 10.

    def _network(self):
        n = pp.network.create_four_substations_node_breaker_network()
        n.update_static_var_compensators(id="SVC", b_max=self.B_MAX_MVAR / 400. ** 2)
        return n

    def _baked(self):
        n = self._network()
        lf.run_ac(n)
        return n, bake_outer_loops(n)

    def test_bake_returns_the_frozen_svc_and_init_flags_it(self):
        n, pinned = self._baked()
        self.assertIn("SVC", set(pinned))
        self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], "REACTIVE_POWER")
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=pinned)
        svc = grid.get_svcs()[0]
        self.assertTrue(svc.can_be_pv)
        self.assertFalse(svc.standby)
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False)
        self.assertFalse(grid.get_svcs()[0].can_be_pv)

    def test_release_reported_where_olf_regulates_again(self):
        n, pinned = self._baked()
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=pinned)
        _solve(grid)
        self.assertEqual(_release(grid.get_physical_violations()), [])  # the reference state

        # less reactive load on the SVC's bus: its voltage rises above the target
        raw = self._network()
        raw.update_loads(id="LD6", q0=-100.)
        lf.run_ac(raw)
        s = raw.get_static_var_compensators().loc["SVC"]
        v = raw.get_buses().loc[s["bus_id"], "v_mag"]
        self.assertAlmostEqual(v, s["target_v"], places=6)   # OLF: back in voltage control
        n.update_loads(id="LD6", q0=-100.)
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=pinned)
        _solve(grid)
        viols = _release(grid.get_physical_violations())
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].violation_type, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)
        self.assertEqual(viols[0].name, "SVC")
        self.assertAlmostEqual(viols[0].limit, 400., places=6)

    def test_switched_on_standby_svc_is_no_longer_standby(self):
        # thresholds its bus is already outside of: OLF switches it on, the bake says so
        n = pp.network.create_four_substations_node_breaker_network()
        n.create_extensions("standbyAutomaton", id="SVC", standby=True, b0=0.,
                            low_voltage_threshold=401., low_voltage_setpoint=402.,
                            high_voltage_threshold=410., high_voltage_setpoint=405.)
        lf.run_ac(n)
        pinned = bake_outer_loops(n)
        self.assertNotIn("SVC", set(pinned))
        self.assertFalse(n.get_extensions("standbyAutomaton").loc["SVC", "standby"])
        self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], "VOLTAGE")


if __name__ == "__main__":
    unittest.main()
