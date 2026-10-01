# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The physical check of the generators regulating a remote bus
(``LSGrid.set_remote_voltage_control_vm_range``): a remote controller whose own bus leaves the
realistic range -- which OpenLoadFlow's robust remote voltage control switches to PQ -- is
reported as ``LOW_VOLTAGE_REMOTE_CONTROL`` / ``HIGH_VOLTAGE_REMOTE_CONTROL`` on the generator,
in kV, by a single solve and by the batches alike, exactly where OpenLoadFlow switches it."""

import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import (
    LSGrid,
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
    """whether OpenLoadFlow's robust remote voltage control switched a controller to PQ"""
    report = rp.ReportNode()
    lf.run_ac(n, lf.Parameters(), report_node=report)
    return any("remote voltage target is maintained" in line for line in str(report).splitlines())


def _ls(n, **kwargs):
    grid = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=True, **kwargs)
    grid.set_gen_names(list(n.get_generators().index))
    V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10)
    assert V.shape[0] > 0, "lightsim2grid diverged"
    return grid


class TestRemoteVoltageControlRange(unittest.TestCase):
    def test_default_and_set(self):
        grid = LSGrid()
        self.assertTrue(np.isnan(grid.get_remote_voltage_control_min_vm_pu()))
        self.assertTrue(np.isnan(grid.get_remote_voltage_control_max_vm_pu()))
        grid.set_remote_voltage_control_vm_range(0.8, 1.2)
        self.assertEqual(grid.get_remote_voltage_control_min_vm_pu(), 0.8)
        self.assertEqual(grid.get_remote_voltage_control_max_vm_pu(), 1.2)
        grid.set_remote_voltage_control_vm_range(np.nan, 1.2)  # one side only
        self.assertTrue(np.isnan(grid.get_remote_voltage_control_min_vm_pu()))
        # copied with the grid (a batch built from it inherits it), like
        # keep_vinit_at_group_controlled_buses not part of the state a pickle carries
        self.assertEqual(grid.copy().get_remote_voltage_control_max_vm_pu(), 1.2)

    def test_refused(self):
        grid = LSGrid()
        for bad in [(1.2, 0.8), (1., 1.), (-0.1, 1.2), (0.8, np.inf)]:
            with self.assertRaises(RuntimeError):
                grid.set_remote_voltage_control_vm_range(*bad)

    def test_violation_types(self):
        for vt in REMOTE_TYPES:
            self.assertEqual(int(vt) in (13, 14), True)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestRemoteVoltageControlCheck(unittest.TestCase):
    def test_reported_where_olf_switches(self):
        # around the high bound, and far below the low one: reported iff OpenLoadFlow switches
        for target_v_kv in [9.5, 13.85, 13.8975, 13.9, 14.2]:
            olf = _olf_switches(_ieee14_remote(target_v_kv))
            viols = _remote(_ls(_ieee14_remote(target_v_kv)).get_physical_violations(True, 0., 0.))
            self.assertEqual(len(viols) == 1, olf, f"target {target_v_kv} kV")
            for el in viols:
                self.assertEqual(el.element_type, ViolationElementType.GENERATOR)
                self.assertEqual(el.name, "B6-G")
                self.assertEqual(el.category, ViolationCategory.PHYSICAL)

    def test_value_and_limit(self):
        grid = _ls(_ieee14_remote(14.2))
        viols = _remote(grid.get_physical_violations(True, 0., 0.))
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].violation_type, LimitViolationType.HIGH_VOLTAGE_REMOTE_CONTROL)
        # the bound in kV of its own (12 kV) bus, OpenLoadFlow's default maxRealisticVoltage
        # with its margin
        self.assertAlmostEqual(viols[0].limit, 12. * 1.2 / 1.02, places=9)
        self.assertGreater(viols[0].value, viols[0].limit)

    def test_switched_off(self):
        grid = _ls(_ieee14_remote(14.2), remote_voltage_control_vm_range=None)
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])
        grid = _ls(_ieee14_remote(14.2), remote_voltage_control_vm_range=(0.5, 1.5))
        self.assertEqual(_remote(grid.get_physical_violations(True, 0., 0.)), [])

    def test_local_controller_not_checked(self):
        # the same generator holding its own bus far above the range: not a remote controller
        n = pp.network.create_ieee14()
        n.update_generators(id="B6-G", target_v=14.5, min_q=-9999., max_q=9999.)
        self.assertEqual(_remote(_ls(n).get_physical_violations(True, 0., 0.)), [])

    def test_contingency_analysis_matches_the_single_solve(self):
        grid = _ls(_ieee14_remote(14.2))
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
