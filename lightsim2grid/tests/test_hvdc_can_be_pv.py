# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The PQ -> PV release of a VSC converter station an outer loop froze at a reactive limit
(``LSGrid.set_hvdc_can_be_pv``, the generators' ``can_be_pv``): reported as
``LOW_VOLTAGE_AT_MIN_Q`` / ``HIGH_VOLTAGE_AT_MAX_Q`` on the HVDC line, ``side`` the station's end;
``bake_outer_loops`` / ``init_from_pypowsybl(can_be_pv=...)`` flag the stations a bake froze."""

import copy
import os
import pickle
import tempfile
import unittest

import numpy as np

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationElementType, LSGrid
    from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP
    from lightsim2grid.network import bake_outer_loops, init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _release(viols):
    return [v for v in viols if v.element_type == ViolationElementType.HVDC
            and v.violation_type in (LimitViolationType.LOW_VOLTAGE_AT_MIN_Q,
                                     LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)]


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestHvdcRelease(unittest.TestCase):
    """VSC2 (side 2 of HVDC1 in the four substations network) regulating 412 kV with at most
    130 MVAr: it saturates, its bus sags below the target, and the bake freezes it at max_q."""

    TARGET_KV = 412.
    # a target below the voltage the frozen station leaves: at max_q and above its target, it
    # produces too much for it -- OLF would regulate it again
    NEW_TARGET_KV = 409.

    def _network(self):
        n = pp.network.create_four_substations_node_breaker_network()
        n.update_vsc_converter_stations(id="VSC2", target_v=self.TARGET_KV, max_q=130.)
        n.update_vsc_converter_stations(id="VSC2", voltage_regulator_on=True)
        return n

    def _baked(self):
        n = self._network()
        lf.run_ac(n)
        return n, bake_outer_loops(n)

    def _grid(self, n, can_be_pv):
        grid = init_from_pypowsybl(n, sort_index=False, buses_for_sub=False, can_be_pv=can_be_pv)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        return grid

    def test_bake_returns_the_frozen_station_and_init_flags_its_side(self):
        n, pinned = self._baked()
        self.assertIn("VSC2", set(pinned))
        self.assertFalse(n.get_vsc_converter_stations().loc["VSC2", "voltage_regulator_on"])
        line = list(self._grid(n, pinned).get_dclines())[0]
        self.assertEqual(line.name, "HVDC1")
        self.assertFalse(line.station1.can_be_pv)
        self.assertTrue(line.station2.can_be_pv)
        self.assertFalse(list(self._grid(n, None).get_dclines())[0].station2.can_be_pv)
        # the reference state: at max_q with its bus BELOW the target, nothing to release
        self.assertEqual(_release(self._grid(n, pinned).get_physical_violations()), [])

    def test_release_reported_where_olf_regulates_again(self):
        raw = self._network()
        raw.update_vsc_converter_stations(id="VSC2", target_v=self.NEW_TARGET_KV)
        lf.run_ac(raw)
        vsc = raw.get_vsc_converter_stations().loc["VSC2"]
        self.assertAlmostEqual(raw.get_buses().loc[vsc["bus_id"], "v_mag"], self.NEW_TARGET_KV, places=6)
        self.assertLess(-vsc["q"], 130. - 1.)   # OLF: back inside its range, holding the target

        n, pinned = self._baked()
        n.update_vsc_converter_stations(id="VSC2", target_v=self.NEW_TARGET_KV)
        grid = self._grid(n, pinned)
        viols = _release(grid.get_physical_violations())
        self.assertEqual(len(viols), 1)
        v = viols[0]
        self.assertEqual(v.violation_type, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)
        self.assertEqual(v.name, "HVDC1")
        self.assertEqual(v.side, 2)
        self.assertAlmostEqual(v.limit, self.NEW_TARGET_KV, places=6)
        self.assertGreater(v.value, v.limit)
        # not flagged: nothing
        self.assertEqual(_release(self._grid(n, None).get_physical_violations()), [])

        # a batch's base case says the same
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.add_n1(0)
        ca.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        n_case = _release(ca.get_physical_violations_n())
        self.assertEqual([(x.element_id, x.side, x.violation_type) for x in n_case],
                         [(v.element_id, 2, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)])

    def test_flag_is_kept_and_checked(self):
        n, pinned = self._baked()
        grid = self._grid(n, None)
        # two hvdc lines: HVDC1 (VSC1 / VSC2) and HVDC2 (LCC1 / LCC2)
        grid.set_hvdc_can_be_pv(np.array([False, False]), np.array([True, False]))
        self.assertTrue(list(grid.get_dclines())[0].station2.can_be_pv)
        for other in (grid.copy(), copy.deepcopy(grid), pickle.loads(pickle.dumps(grid))):
            self.assertTrue(list(other.get_dclines())[0].station2.can_be_pv)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            self.assertTrue(list(LSGrid.load_binary(path).get_dclines())[0].station2.can_be_pv)
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_can_be_pv(np.array([True, True, True]), np.array([True, False]))


if __name__ == "__main__":
    unittest.main()
