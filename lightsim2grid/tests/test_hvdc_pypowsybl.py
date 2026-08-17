# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""HVDC parity against PowSyBl OpenLoadFlow through the pypowsybl converter.

Covers the four regimes: angle-droop linear, fixed setpoint (VSC), fixed setpoint
(LCC inverter side) and forced saturation. Skipped when pypowsybl is unavailable.
"""

import unittest
import numpy as np

try:
    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _build_net(max_p=300.0, droop_enabled=True, lcc=False, target_p=50.0):
    n = pp.network.create_empty()
    n.create_substations(id=["S1", "S2"])
    n.create_voltage_levels(id=["VL1", "VL2"], substation_id=["S1", "S2"],
                            topology_kind=["BUS_BREAKER"] * 2, nominal_v=[400.0, 400.0])
    n.create_buses(id=["B1", "B2"], voltage_level_id=["VL1", "VL2"])
    n.create_lines(id="L", voltage_level1_id="VL1", voltage_level2_id="VL2",
                   bus1_id="B1", bus2_id="B2", r=1.0, x=20.0, g1=0.0, b1=0.0, g2=0.0, b2=0.0)
    n.create_vsc_converter_stations(id="VSC1", voltage_level_id="VL1", bus_id="B1", loss_factor=1.1,
                                    voltage_regulator_on=True, target_v=400.0, target_q=0.0)
    if lcc:
        n.create_lcc_converter_stations(id="LCC2", voltage_level_id="VL2", bus_id="B2",
                                        loss_factor=1.1, power_factor=0.9)
        cs2 = "LCC2"
    else:
        n.create_vsc_converter_stations(id="VSC2", voltage_level_id="VL2", bus_id="B2", loss_factor=1.1,
                                        voltage_regulator_on=True, target_v=400.0, target_q=0.0)
        cs2 = "VSC2"
    n.create_hvdc_lines(id="HVDC", converter_station1_id="VSC1", converter_station2_id=cs2,
                        r=1.0, nominal_v=400.0, converters_mode="SIDE_1_RECTIFIER_SIDE_2_INVERTER",
                        target_p=target_p, max_p=max_p)
    n.create_generators(id="G", voltage_level_id="VL1", bus_id="B1",
                        target_p=100.0, target_q=0.0, target_v=400.0,
                        voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=100.0, q0=0.0)
    if droop_enabled:
        n.create_extensions("hvdcAngleDroopActivePowerControl", id="HVDC",
                            p0=20.0, droop=180.0, enabled=True)
    return n


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestHvdcPypowsybl(unittest.TestCase):
    def setUp(self):
        self.params = pp.loadflow.Parameters(distributed_slack=False, use_reactive_limits=False)
        # OLF carries its own floating-point noise in the loss formulas: ~1e-6 is the
        # achievable parity, not 1e-8.
        self.tol_v = 5e-6
        self.tol_a = 5e-6

    def _run_ls(self, n):
        n.per_unit = False
        model = init_from_pypowsybl(n, gen_slack_id="G", sort_index=True)
        nb_bus = len(model.get_bus_status())
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, 1e-10)
        self.assertGreater(V.shape[0], 0, "lightsim2grid diverged")
        return model, V

    def _compare(self, n, name):
        ref = pp.loadflow.run_ac(n, parameters=self.params)
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        n.per_unit = False
        buses_si = n.get_buses()
        vlevels = n.get_voltage_levels()
        model, V = self._run_ls(n)
        for i, (bus_id, row) in enumerate(buses_si.iterrows()):
            nominal_v = vlevels.loc[row["voltage_level_id"], "nominal_v"]
            v_olf = row["v_mag"] / nominal_v
            a_olf = np.deg2rad(row["v_angle"])
            self.assertLess(abs(abs(V[i]) - v_olf), self.tol_v, f"{name}: V at {bus_id}")
            self.assertLess(abs(np.angle(V[i]) - a_olf), self.tol_a, f"{name}: angle at {bus_id}")
        return model

    def test_droop_linear(self):
        model = self._compare(_build_net(max_p=300.0), "droop-linear")
        h = model.get_dclines()[0]
        self.assertTrue(h.droop_enabled)
        self.assertEqual(model.get_status_droop_hvdc(0), 0)

    def test_fixed_vsc(self):
        model = self._compare(_build_net(droop_enabled=False), "fixed-vsc")
        self.assertFalse(model.get_dclines()[0].droop_enabled)

    def test_fixed_lcc(self):
        model = self._compare(_build_net(droop_enabled=False, lcc=True), "fixed-lcc")
        self.assertEqual(model.get_dclines()[0].station2.converter_type, 1)

    def test_saturated_droop(self):
        # OLF saturates the droop at max_p; we reproduce the single saturation switch the
        # python outer loop would perform and re-check parity.
        n = _build_net(max_p=10.0)
        ref = pp.loadflow.run_ac(n, parameters=self.params)
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        n.per_unit = False
        buses_si = n.get_buses()
        vlevels = n.get_voltage_levels()
        model, _ = self._run_ls(n)
        h = model.get_dclines()[0]
        raw = h.droop_p0_mw + h.droop_k_mw_per_rad * np.deg2rad(h.res_theta1_deg - h.res_theta2_deg)
        self.assertGreater(raw, h.pmax_1to2_mw)  # the linear solve overshoots the cap
        model.set_status_droop_hvdc(0, 1)
        nb_bus = len(model.get_bus_status())
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        for i, (bus_id, row) in enumerate(buses_si.iterrows()):
            nominal_v = vlevels.loc[row["voltage_level_id"], "nominal_v"]
            self.assertLess(abs(abs(V[i]) - row["v_mag"] / nominal_v), 1e-6, f"saturated V {bus_id}")
            # In the saturated regime the reference angle comes from OLF's droop-saturation +
            # loss formulas, whose last digits depend on the pypowsybl/OLF build: locally the
            # agreement is ~1e-13, but other environments (CI) drift to ~1e-6. Allow that band
            # for the angle so the test tracks parity rather than a single build's rounding.
            self.assertLess(abs(np.angle(V[i]) - np.deg2rad(row["v_angle"])), 5e-6, f"saturated ang {bus_id}")


if __name__ == "__main__":
    unittest.main()
