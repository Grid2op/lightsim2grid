# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Per-bus voltage limits and per-side branch thermal (current) limits, read from
pypowsybl's `low_voltage_limit`/`high_voltage_limit` (voltage levels) and
`get_operational_limits()` (lines / transformers) when a grid is built with
`init_from_pypowsybl`. Skipped when pypowsybl is unavailable."""

import os
import pickle
import tempfile
import unittest
import numpy as np

from lightsim2grid.network import LSGrid

try:
    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

try:
    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    HAS_PANDAPOWER = True
except ImportError:
    HAS_PANDAPOWER = False


def _build_net():
    """3-bus network: slack gen (B1) -- line L12 -- B2 -- trafo T23 -- B3 (load).

    B1's voltage level carries no voltage limit (pypowsybl default: NaN).
    L12 carries an asymmetric per-side current limit (side 1: permanent + a higher
    temporary duration, to check the min-over-durations reduction; side 2: a single
    permanent limit). T23 carries a current limit on side 1 only (side 2 stays unset).
    """
    n = pp.network.create_empty()
    # T23 needs both its voltage levels on the same substation (IIDM constraint)
    n.create_substations(id=["S1", "S2"])
    n.create_voltage_levels(id=["VL1", "VL2", "VL3"], substation_id=["S1", "S2", "S2"],
                            topology_kind=["BUS_BREAKER"] * 3, nominal_v=[400.0, 400.0, 90.0],
                            low_voltage_limit=[float("nan"), 380.0, 80.0],
                            high_voltage_limit=[float("nan"), 420.0, 100.0])
    n.create_buses(id=["B1", "B2", "B3"], voltage_level_id=["VL1", "VL2", "VL3"])
    n.create_lines(id="L12", voltage_level1_id="VL1", voltage_level2_id="VL2",
                   bus1_id="B1", bus2_id="B2", r=1.0, x=20.0, g1=0.0, b1=0.0, g2=0.0, b2=0.0)
    n.create_2_windings_transformers(id="T23", voltage_level1_id="VL2", voltage_level2_id="VL3",
                                     bus1_id="B2", bus2_id="B3", rated_u1=400.0, rated_u2=90.0,
                                     r=1.0, x=10.0, g=0.0, b=0.0, rated_s=100.0)
    n.create_generators(id="G", voltage_level_id="VL1", bus_id="B1",
                        target_p=100.0, target_q=0.0, target_v=400.0,
                        voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_loads(id="LD", voltage_level_id="VL3", bus_id="B3", p0=80.0, q0=20.0)

    # side 1: permanent=500A + a higher 10min=600A temporary limit -> min=500A=0.5kA
    n.create_operational_limits(element_id=["L12", "L12"],
                                side=["ONE", "ONE"],
                                name=["permanent", "10min"],
                                type=["CURRENT"] * 2,
                                value=[500.0, 600.0],
                                acceptable_duration=[-1, 600])
    # side 2: single permanent limit -> 450A = 0.45kA
    n.create_operational_limits(element_id=["L12"], side=["TWO"], name=["permanent"],
                                type=["CURRENT"], value=[450.0], acceptable_duration=[-1])
    # trafo side 1 only (side 2 stays unset -> NaN)
    n.create_operational_limits(element_id=["T23"], side=["ONE"], name=["permanent"],
                                type=["CURRENT"], value=[123.4], acceptable_duration=[-1])
    return n


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestLimitsPypowsybl(unittest.TestCase):
    def setUp(self):
        self.model = init_from_pypowsybl(_build_net(), gen_slack_id="G", sort_index=True)

    def test_voltage_limits(self):
        vmin = np.asarray(self.model.get_bus_vmin_kv())
        vmax = np.asarray(self.model.get_bus_vmax_kv())
        self.assertEqual(vmin.shape[0], 3)
        sub_names = list(self.model.get_substation_names())
        by_name = {nm: i for i, nm in enumerate(sub_names)}
        self.assertTrue(np.isnan(vmin[by_name["VL1"]]))
        self.assertTrue(np.isnan(vmax[by_name["VL1"]]))
        self.assertAlmostEqual(vmin[by_name["VL2"]], 380.0)
        self.assertAlmostEqual(vmax[by_name["VL2"]], 420.0)
        self.assertAlmostEqual(vmin[by_name["VL3"]], 80.0)
        self.assertAlmostEqual(vmax[by_name["VL3"]], 100.0)

    def test_line_current_limit(self):
        line = [el for el in self.model.get_lines() if el.name == "L12"][0]
        self.assertAlmostEqual(line.limit_a1_ka, 0.5, places=6)  # min(500, 600) A
        self.assertAlmostEqual(line.limit_a2_ka, 0.45, places=6)

    def test_trafo_current_limit(self):
        trafo = [el for el in self.model.get_trafos() if el.name == "T23"][0]
        self.assertAlmostEqual(trafo.limit_a1_ka, 0.1234, places=6)
        self.assertTrue(np.isnan(trafo.limit_a2_ka))

    def test_no_limits_grid_gives_nan_not_crash(self):
        n = pp.network.create_ieee14()
        model = init_from_pypowsybl(n)
        vmin = np.asarray(model.get_bus_vmin_kv())
        self.assertEqual(vmin.shape[0], model.get_bus_status().__len__())
        for el in model.get_lines():
            self.assertTrue(np.isnan(el.limit_a1_ka))
            self.assertTrue(np.isnan(el.limit_a2_ka))

    def test_pickle_roundtrip_pypowsybl_grid(self):
        # see also TestPickleFromPypo in test_init_from_pypowsybl.py for the general
        # (limits-independent) pickle/binary regression test
        model2 = pickle.loads(pickle.dumps(self.model))
        np.testing.assert_allclose(np.asarray(model2.get_bus_vmin_kv()), np.asarray(self.model.get_bus_vmin_kv()), equal_nan=True)
        np.testing.assert_allclose(np.asarray(model2.get_bus_vmax_kv()), np.asarray(self.model.get_bus_vmax_kv()), equal_nan=True)
        line = [el for el in model2.get_lines() if el.name == "L12"][0]
        self.assertAlmostEqual(line.limit_a1_ka, 0.5, places=6)
        self.assertAlmostEqual(line.limit_a2_ka, 0.45, places=6)


def _build_ls_grid_with_limits():
    """Plain LSGrid (no pypowsybl involved), exercising the new setters directly --
    a minimal/fast complement to `test_pickle_roundtrip_pypowsybl_grid` above."""
    g = LSGrid()
    g.init_bus(2, 1, np.array([20.0, 20.0]), 0, 0)
    g.set_bus_voltage_limits(np.array([18.0, np.nan]), np.array([22.0, np.nan]))
    g.init_powerlines(np.array([0.01]), np.array([0.05]), np.array([0.0 + 0.0j]),
                      np.array([0], dtype=np.int32), np.array([1], dtype=np.int32))
    g.set_line_current_limit_side1(np.array([1.234]))
    g.set_line_current_limit_side2(np.array([np.nan]))
    return g


class TestLimitsSerialization(unittest.TestCase):
    def test_pickle_roundtrip(self):
        model = _build_ls_grid_with_limits()
        model2 = pickle.loads(pickle.dumps(model))
        np.testing.assert_allclose(np.asarray(model2.get_bus_vmin_kv()), np.asarray(model.get_bus_vmin_kv()), equal_nan=True)
        np.testing.assert_allclose(np.asarray(model2.get_bus_vmax_kv()), np.asarray(model.get_bus_vmax_kv()), equal_nan=True)
        a1_before = [el.limit_a1_ka for el in model.get_lines()]
        a1_after = [el.limit_a1_ka for el in model2.get_lines()]
        np.testing.assert_allclose(a1_after, a1_before, equal_nan=True)

    def test_binary_roundtrip(self):
        model = _build_ls_grid_with_limits()
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "model.bin")
            model.save_binary(path)
            model2 = LSGrid.load_binary(path)
        np.testing.assert_allclose(np.asarray(model2.get_bus_vmin_kv()), np.asarray(model.get_bus_vmin_kv()), equal_nan=True)
        a1_before = [el.limit_a1_ka for el in model.get_lines()]
        a1_after = [el.limit_a1_ka for el in model2.get_lines()]
        np.testing.assert_allclose(a1_after, a1_before, equal_nan=True)

    def test_unset_fields_roundtrip_empty(self):
        g = LSGrid()
        g.init_bus(1, 1, np.array([20.0]), 0, 0)
        g.init_powerlines(np.array([0.01]), np.array([0.05]), np.array([0.0 + 0.0j]),
                          np.array([0], dtype=np.int32), np.array([0], dtype=np.int32))
        g2 = pickle.loads(pickle.dumps(g))
        self.assertEqual(np.asarray(g2.get_bus_vmin_kv()).shape[0], 0)
        self.assertTrue(np.isnan([el.limit_a1_ka for el in g2.get_lines()][0]))


@unittest.skipUnless(HAS_PANDAPOWER, "pandapower is not installed")
class TestLimitsPandapowerUnset(unittest.TestCase):
    """pandapower-origin grids never call the new setters: getters must return
    empty arrays / per-element NaN, not crash."""

    def test_unset_limits(self):
        model = init_from_pandapower(pn.case14())
        self.assertEqual(np.asarray(model.get_bus_vmin_kv()).shape[0], 0)
        self.assertEqual(np.asarray(model.get_bus_vmax_kv()).shape[0], 0)
        for el in model.get_lines():
            self.assertTrue(np.isnan(el.limit_a1_ka))
            self.assertTrue(np.isnan(el.limit_a2_ka))


if __name__ == "__main__":
    unittest.main()
