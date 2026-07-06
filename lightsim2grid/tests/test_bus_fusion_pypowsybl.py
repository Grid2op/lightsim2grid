# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`fuse_zero_impedance_branches` (from_pypowsybl `init()`).

A real RTE grid was found to have a line with `r_pu = x_pu = 0.0`. lightsim2grid
computed `B = 1/X = Inf`, corrupting the Ybus/dcYbus with literal `+/-Inf` entries
and making the sparse LU factorization fail outright. OpenLoadFlow avoids this via
its `lowImpedanceBranchMode` (default `REPLACE_BY_ZERO_IMPEDANCE_LINE`): the two
terminal buses of a near-zero-impedance branch are fused into a single electrical
node instead of computing `1/Z`. This is the same behaviour, opt-in via
`fuse_zero_impedance_branches` (see HVDC_OLF_FINDINGS.md and the docstring of
`lightsim2grid.network.from_pypowsybl.init`).
"""

import unittest
import warnings

import numpy as np

try:
    import pypowsybl.network as pn
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _dc_pf_ok(model):
    nb_bus = model.get_bus_vn_kv().shape[0]
    flat = np.full(nb_bus, 1.0, dtype=complex)
    Vdc = model.dc_pf(flat, 30, 1e-8)
    return Vdc, Vdc.shape[0] > 0


class TestBusFusionZeroImpedanceLine(unittest.TestCase):
    """A single r=x=0 line in an otherwise normal (ieee14-sized) grid."""

    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        self.net = pn.create_ieee14()
        # L6-11-1 connects two non-slack, non-reference buses (unlike L1-2-1,
        # which happens to touch the reference bus and would not reproduce the
        # crash: eliminating the reference row/column also eliminates the Inf
        # entries by coincidence).
        self.net.update_lines(id="L6-11-1", r=0.0, x=0.0)

    def test_flag_off_reproduces_the_bug(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(self.net, sort_index=False)
        _, ok = _dc_pf_ok(model)
        self.assertFalse(ok, "dc_pf should fail without fuse_zero_impedance_branches "
                             "(Inf entries from B=1/X=1/0 break the LU factorization)")

    def test_flag_on_fixes_it_and_fuses_buses(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(self.net, sort_index=False,
                                        fuse_zero_impedance_branches=True)
        line = next(l for l in model.get_lines() if l.name == "L6-11-1")
        # the fusing line is kept (for bookkeeping) but fully deactivated
        self.assertFalse(line.connected1)
        self.assertFalse(line.connected2)
        self.assertEqual(line.bus1_id, -1)
        self.assertEqual(line.bus2_id, -1)

        Vdc, ok = _dc_pf_ok(model)
        self.assertTrue(ok, "dc_pf should succeed with fuse_zero_impedance_branches=True")
        V = model.ac_pf(Vdc, 30, 1e-8)
        self.assertGreater(V.shape[0], 0, "ac_pf should also succeed")

        # find the two ORIGINAL terminal buses of L6-11-1 (before fusion) using a
        # model built WITHOUT the flag (where L6-11-1 still has real bus ids), via
        # two OTHER lines each touching exactly one of the two original buses
        # (bus "6" also carries L6-12-1/L6-13-1, bus "11" also carries L10-11-1).
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model_off = init_from_pypowsybl(self.net, sort_index=False)
        by_name_off = {l.name: l for l in model_off.get_lines()}
        fused = by_name_off["L6-11-1"]
        bus_side1, bus_side2 = fused.bus1_id, fused.bus2_id

        def _anchor(name, target_bus):
            l = by_name_off[name]
            return l.bus1_id if l.bus1_id == target_bus else l.bus2_id

        anchor1_id = _anchor("L6-12-1", bus_side1)
        anchor2_id = _anchor("L10-11-1", bus_side2)
        self.assertEqual(anchor1_id, bus_side1)
        self.assertEqual(anchor2_id, bus_side2)

        by_name_on = {l.name: l for l in model.get_lines()}

        def _same_side_bus(name, target_bus):
            l = by_name_on[name]
            return l.bus1_id if l.bus1_id == target_bus or l.bus2_id != target_bus else l.bus2_id

        # since fusion never renumbers an unaffected bus, whichever of the two
        # original ids is NOT the (smaller-id) representative disappears from
        # every OTHER element too, replaced by the representative -- so the two
        # anchor lines' relevant side should now agree.
        b1_on = by_name_on["L6-12-1"].bus1_id if by_name_off["L6-12-1"].bus1_id == bus_side1 else by_name_on["L6-12-1"].bus2_id
        b2_on = by_name_on["L10-11-1"].bus1_id if by_name_off["L10-11-1"].bus1_id == bus_side2 else by_name_on["L10-11-1"].bus2_id
        self.assertEqual(b1_on, b2_on,
                          "the two original terminal buses of the fused line "
                          "should now resolve to the same solver bus")
        self.assertAlmostEqual(V[b1_on], V[b2_on])


class TestBusFusionCrossVoltageLineRaises(unittest.TestCase):
    """A zero-impedance line spanning two different nominal voltages is
    inconsistent data and must raise, never be silently fused."""

    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        self.net = pn.create_empty()
        self.net.create_substations(id="S1")
        self.net.create_voltage_levels(substation_id="S1", id="VL1", nominal_v=225.0,
                                       topology_kind="BUS_BREAKER")
        self.net.create_voltage_levels(substation_id="S1", id="VL2", nominal_v=90.0,
                                       topology_kind="BUS_BREAKER")
        self.net.create_buses(voltage_level_id="VL1", id="B1")
        self.net.create_buses(voltage_level_id="VL2", id="B2")
        self.net.create_lines(id="L1", voltage_level1_id="VL1", bus1_id="B1",
                              voltage_level2_id="VL2", bus2_id="B2",
                              r=0.0, x=0.0, g1=0.0, b1=0.0, g2=0.0, b2=0.0)

    def test_raises(self):
        with self.assertRaises(RuntimeError):
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                init_from_pypowsybl(self.net, sort_index=False,
                                    fuse_zero_impedance_branches=True,
                                    only_main_component=False)


class TestBusFusionTransformer(unittest.TestCase):
    """A zero-impedance, neutral-tap transformer between two SAME-nominal-voltage
    buses is a genuine ideal wire and gets fused; the same transformer between two
    DIFFERENT nominal voltages must NOT be fused, even though pypowsybl's per-unit
    `rho` still reports ~1 (rho is the deviation from the transformer's own rated
    ratio, not its absolute turns ratio)."""

    def _build(self, vn2):
        net = pn.create_empty()
        net.create_substations(id="S1")
        net.create_voltage_levels(substation_id="S1", id="VL1", nominal_v=225.0,
                                  topology_kind="BUS_BREAKER")
        net.create_voltage_levels(substation_id="S1", id="VL2", nominal_v=225.0,
                                  topology_kind="BUS_BREAKER")
        net.create_voltage_levels(substation_id="S1", id="VL3", nominal_v=vn2,
                                  topology_kind="BUS_BREAKER")
        net.create_buses(voltage_level_id="VL1", id="B1")
        net.create_buses(voltage_level_id="VL2", id="B2")
        net.create_buses(voltage_level_id="VL3", id="B3")
        net.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1",
                              target_p=10.0, target_v=225.0,
                              min_p=0.0, max_p=100.0, voltage_regulator_on=True)
        net.create_loads(id="LD1", voltage_level_id="VL3", bus_id="B3", p0=5.0, q0=1.0)
        # the candidate zero-impedance transformer: B1 (225kV) <-> B2 (225kV)
        net.create_2_windings_transformers(id="T1", voltage_level1_id="VL1", bus1_id="B1",
                                           voltage_level2_id="VL2", bus2_id="B2",
                                           r=0.0, x=0.0, g=0.0, b=0.0,
                                           rated_u1=225.0, rated_u2=225.0)
        # a real, non-zero-impedance line so the system doesn't collapse to a
        # single trivial bus after fusion (a pre-existing, unrelated lightsim2grid
        # issue with 0-non-slack-bus networks, not something this feature should
        # need to work around by construction).
        net.create_lines(id="L23", voltage_level1_id="VL2", bus1_id="B2",
                         voltage_level2_id="VL3", bus2_id="B3",
                         r=0.01, x=0.1, g1=0.0, b1=0.0, g2=0.0, b2=0.0)
        return net

    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")

    def test_same_nominal_voltage_is_fused(self):
        net = self._build(vn2=225.0)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(net, sort_index=False,
                                        fuse_zero_impedance_branches=True)
        trafo = next(t for t in model.get_trafos() if t.name == "T1")
        self.assertFalse(trafo.connected1)
        self.assertFalse(trafo.connected2)

        gen_bus = next(g.bus_id for g in model.get_generators())
        line_bus1 = next(l.bus1_id for l in model.get_lines())
        self.assertEqual(gen_bus, line_bus1,
                          "B1 and B2 should be fused to the same solver bus")

        Vdc, ok = _dc_pf_ok(model)
        self.assertTrue(ok)
        V = model.ac_pf(Vdc, 30, 1e-8)
        self.assertGreater(V.shape[0], 0)

    def test_different_nominal_voltage_is_not_fused(self):
        # B1 is 225kV, "VL3"/B3 stays 225kV too here but T1 itself is still
        # between B1 (225kV) and B2 (225kV) -- to actually test the cross-voltage
        # guard we need T1's OWN two sides to differ, so build a dedicated net.
        net = pn.create_empty()
        net.create_substations(id="S1")
        net.create_voltage_levels(substation_id="S1", id="VL1", nominal_v=225.0,
                                  topology_kind="BUS_BREAKER")
        net.create_voltage_levels(substation_id="S1", id="VL2", nominal_v=90.0,
                                  topology_kind="BUS_BREAKER")
        net.create_buses(voltage_level_id="VL1", id="B1")
        net.create_buses(voltage_level_id="VL2", id="B2")
        net.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1",
                              target_p=10.0, target_v=225.0,
                              min_p=0.0, max_p=100.0, voltage_regulator_on=True)
        net.create_loads(id="LD1", voltage_level_id="VL2", bus_id="B2", p0=5.0, q0=1.0)
        net.create_2_windings_transformers(id="T1", voltage_level1_id="VL1", bus1_id="B1",
                                           voltage_level2_id="VL2", bus2_id="B2",
                                           r=0.0, x=0.0, g=0.0, b=0.0,
                                           rated_u1=225.0, rated_u2=90.0)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(net, sort_index=False,
                                        fuse_zero_impedance_branches=True)
        trafo = next(t for t in model.get_trafos() if t.name == "T1")
        self.assertTrue(trafo.connected1, "a genuine (even if lossless) "
                        "step-up/down transformer must never be fused")
        self.assertTrue(trafo.connected2)


class TestBusFusionFlagOffUnchanged(unittest.TestCase):
    """With the flag omitted (default), behaviour must be identical to before
    this feature existed."""

    def test_default_is_false_and_unchanged(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        net = pn.create_ieee14()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model_default = init_from_pypowsybl(net, sort_index=False)
            model_explicit_off = init_from_pypowsybl(net, sort_index=False,
                                                     fuse_zero_impedance_branches=False)
        self.assertEqual(model_default.get_bus_vn_kv().shape[0],
                         model_explicit_off.get_bus_vn_kv().shape[0])
        Vdc1, ok1 = _dc_pf_ok(model_default)
        Vdc2, ok2 = _dc_pf_ok(model_explicit_off)
        self.assertEqual(ok1, ok2)
        if ok1:
            np.testing.assert_array_almost_equal(Vdc1, Vdc2)


if __name__ == "__main__":
    unittest.main()
