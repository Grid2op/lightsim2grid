# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The trigger of ReactiveLimits' robust remote voltage control, in detection mode: a remote
controller whose own bus leaves the loop's realistic range -- which OpenLoadFlow's robust
remote voltage control switches to PQ -- is
reported as ``LOW_VOLTAGE_REMOTE_CONTROL`` / ``HIGH_VOLTAGE_REMOTE_CONTROL`` on the generator,
in kV, by a single solve and by the batches alike, exactly where OpenLoadFlow switches it."""

import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import (
    LimitViolationType,
    ViolationCategory,
    ViolationElementType,
)
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    import pypowsybl.report as rp
    from lightsim2grid.network import init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl._olf_const import (
        _OLF_MIN_REALISTIC_VOLTAGE_PU, _OLF_MAX_REALISTIC_VOLTAGE_PU)
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

REMOTE_TYPES = (LimitViolationType.LOW_VOLTAGE_REMOTE_CONTROL,
                LimitViolationType.HIGH_VOLTAGE_REMOTE_CONTROL)


def _remote(viols):
    return [el for el in viols if el.violation_type in REMOTE_TYPES]


def _ieee14_remote(target_v_kv):
    """IEEE 14 with B6-G (12 kV bus) regulating the 12 kV bus of a load further away, with
    reactive limits wide enough that only its own voltage can stop it."""
    n = pp.network.create_ieee14()
    loads = n.get_loads(attributes=["voltage_level_id"])
    load = loads.index[loads["voltage_level_id"] == "VL12"][0]
    n.update_generators(id="B6-G", regulated_element_id=load, target_v=target_v_kv,
                        min_q=-9999., max_q=9999.)
    return n


def _olf_switches(n):
    """whether OpenLoadFlow's robust remote voltage control switched a controller to PQ.

    Its realistic voltage range is passed explicitly, the one the tests give lightsim2grid's
    ReactiveLimits loop: pypowsybl builds ship different defaults for it (the public one a range so
    wide the robust mode never triggers here)."""
    params = lf.Parameters()
    params.provider_parameters = {"voltageRemoteControlRobustMode": "true",
                                  "minRealisticVoltage": str(_OLF_MIN_REALISTIC_VOLTAGE_PU),
                                  "maxRealisticVoltage": str(_OLF_MAX_REALISTIC_VOLTAGE_PU)}
    report = rp.ReportNode()
    lf.run_ac(n, params, report_node=report)
    return any("remote voltage target is maintained" in line for line in str(report).splitlines())


def _ls(n, **kwargs):
    grid = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=True, **kwargs)
    grid.set_gen_names(list(n.get_generators().index))
    V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10)
    assert V.shape[0] > 0, "lightsim2grid diverged"
    return grid


def _with_range(grid, vmin, vmax):
    """the grid, its only loop a ReactiveLimits with that realistic range"""
    from lightsim2grid.algorithm import ReactiveLimits
    grid.clear_outer_loops()
    grid.add_outer_loop(ReactiveLimits(min_realistic_voltage=vmin, max_realistic_voltage=vmax))
    return grid


class TestRemoteVoltageControlTypes(unittest.TestCase):
    def test_violation_types(self):
        for vt in REMOTE_TYPES:
            self.assertEqual(int(vt) in (13, 14), True)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestRemoteVoltageControlCheck(unittest.TestCase):
    def test_reported_where_olf_switches(self):
        # around the high bound, and far below the low one: reported iff OpenLoadFlow switches
        for target_v_kv in [9.5, 13.85, 13.8975, 13.9, 14.2]:
            olf = _olf_switches(_ieee14_remote(target_v_kv))
            grid = _with_range(_ls(_ieee14_remote(target_v_kv)), _OLF_MIN_REALISTIC_VOLTAGE_PU,
                               _OLF_MAX_REALISTIC_VOLTAGE_PU)
            viols = _remote(grid.get_physical_violations(True, 0., 0.))
            self.assertEqual(len(viols) == 1, olf, f"target {target_v_kv} kV")
            for el in viols:
                self.assertEqual(el.element_type, ViolationElementType.GENERATOR)
                self.assertEqual(el.name, "B6-G")
                self.assertEqual(el.category, ViolationCategory.PHYSICAL)

    def test_value_and_limit(self):
        # an explicit range (OpenLoadFlow's default differs from one build to the next)
        grid = _with_range(_ls(_ieee14_remote(14.2)), 0.8, 1.2)
        viols = _remote(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].violation_type, LimitViolationType.HIGH_VOLTAGE_REMOTE_CONTROL)
        # the bound in kV of its own (12 kV) bus
        self.assertAlmostEqual(viols[0].limit, 12. * 1.2 / 1.02, places=9)
        self.assertGreater(viols[0].value, viols[0].limit)

    def test_switched_off(self):
        # the trigger is ReactiveLimits' robust mode: without it no loop would switch the
        # controller, and nothing is reported; with a range wide enough, neither
        from lightsim2grid.algorithm import ReactiveLimits
        grid = _ls(_ieee14_remote(14.2))
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits(robust_mode=False))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits(min_realistic_voltage=0.5, max_realistic_voltage=1.5))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])

    def test_range_and_mode_follow_the_driver(self):
        # one parameter each in OpenLoadFlow: a ReactiveLimits given none of its own reads the
        # driver's, in a single solve and in a batch alike; one given its own keeps it
        from lightsim2grid.algorithm import OuterLoopDriverParams, ReactiveLimits
        grid = _ls(_ieee14_remote(14.2))
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits())
        grid.set_outer_loop_driver_params(OuterLoopDriverParams(min_realistic_voltage=0.5, max_realistic_voltage=2.))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])
        grid.set_outer_loop_driver_params(OuterLoopDriverParams(min_realistic_voltage=0.8, max_realistic_voltage=1.2))
        viols = _remote(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertAlmostEqual(viols[0].limit, 12. * 1.2 / 1.02, places=9)
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.add_n1(0)
        ca.compute(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10)
        self.assertEqual(len(_remote(ca.get_physical_violations_n())), 1)
        grid.set_outer_loop_driver_params(OuterLoopDriverParams(voltage_remote_control_robust_mode=False,
                                                                min_realistic_voltage=0.8, max_realistic_voltage=1.2))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])
        # the loop's own values win over the driver's
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits(robust_mode=True))
        self.assertEqual(len(_remote(grid.get_physical_violations(True, 0., 0.))), 1)
        grid.clear_outer_loops()
        grid.add_outer_loop(ReactiveLimits(robust_mode=True, min_realistic_voltage=0.5, max_realistic_voltage=1.5))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])

    def test_local_controller_not_checked(self):
        # the same generator holding its own bus far above the range: not a remote controller
        n = pp.network.create_ieee14()
        n.update_generators(id="B6-G", target_v=14.5, min_q=-9999., max_q=9999.)
        self.assertEqual(_remote(_ls(n).get_physical_violations(True, 0., 0.)), [])

    def test_contingency_analysis_matches_the_single_solve(self):
        grid = _with_range(_ls(_ieee14_remote(14.2)), 0.8, 1.2)
        ref = _remote(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(ref), 1)
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.physical_violation_tol_vm_pu = 0.
        ca.add_n1(0)
        ca.compute(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10)
        n_case = _remote(ca.get_physical_violations_n())
        self.assertEqual(len(n_case), 1)
        self.assertAlmostEqual(n_case[0].value, ref[0].value, places=6)
        self.assertEqual(len(_remote(ca.get_physical_violations()[0])), 1)


if __name__ == "__main__":
    unittest.main()
