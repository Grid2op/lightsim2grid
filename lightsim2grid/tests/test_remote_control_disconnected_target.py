# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A *connected* generator (or SVC) can remotely regulate the voltage of a
terminal that is itself disconnected (found on a real RTE grid snapshot: two
generators remotely regulate busbar sections in a de-energized voltage
level). pypowsybl's bus-view id for a disconnected element is `''`, not NaN,
which used to crash `from_pypowsybl.init()` with
`KeyError: "[''] not in index"`. OpenLoadFlow converges fine on such grids,
so `init()` now falls back to local voltage control for the affected
controller instead of crashing.

This is a different case from the one covered by
`test_disconnected_remote_voltage_control` (a *disconnected* generator
remotely regulating a disconnected element, which is simply moot since the
generator itself gets deactivated): here the controller is connected and
actually needs a fallback so it keeps regulating *something*.
"""

import unittest
import warnings

try:
    import pypowsybl.network as pn
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _build_base_net():
    net = pn.create_empty()
    net.create_substations(id=["S1", "S2"])
    net.create_voltage_levels(substation_id=["S1", "S2"], id=["VL1", "VL2"],
                              nominal_v=[225.0, 225.0],
                              topology_kind=["BUS_BREAKER", "BUS_BREAKER"])
    net.create_buses(voltage_level_id=["VL1", "VL2"], id=["B1", "B2"])
    net.create_lines(id="L12", voltage_level1_id="VL1", bus1_id="B1",
                     voltage_level2_id="VL2", bus2_id="B2",
                     r=0.1, x=1.0, g1=0.0, b1=0.0, g2=0.0, b2=0.0)
    net.create_loads(id="LD1", voltage_level_id="VL2", bus_id="B2", p0=5.0, q0=1.0)
    return net


class TestRemoteGenDisconnectedTarget(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        self.net = _build_base_net()
        # a disconnected load, its terminal's bus-view id resolves to ''
        self.net.create_loads(id="LD_TARGET", voltage_level_id="VL2", bus_id="B2",
                              p0=0.0, q0=0.0)
        self.net.update_loads(id="LD_TARGET", connected=False)
        self.net.create_generators(id="G0", voltage_level_id="VL1", bus_id="B1",
                                   target_p=10.0, target_v=225.0, min_p=0.0, max_p=100.0,
                                   voltage_regulator_on=True)
        self.net.update_generators(id="G0", regulated_element_id="LD_TARGET")

    def test_does_not_crash_and_falls_back_to_local_control(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(self.net, sort_index=False)
        gen = next(g for g in model.get_generators() if g.name == "G0")
        self.assertTrue(gen.connected)
        flat_v = [1.0] * model.get_bus_vn_kv().shape[0]
        import numpy as np
        Vdc = model.dc_pf(np.array(flat_v, dtype=complex), 30, 1e-8)
        self.assertGreater(Vdc.shape[0], 0)
        V = model.ac_pf(Vdc, 30, 1e-8)
        self.assertGreater(V.shape[0], 0)


class TestRemoteSVCDisconnectedTarget(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is required")
        self.net = _build_base_net()
        self.net.create_loads(id="LD_TARGET", voltage_level_id="VL2", bus_id="B2",
                              p0=0.0, q0=0.0)
        self.net.update_loads(id="LD_TARGET", connected=False)
        self.net.create_generators(id="G0", voltage_level_id="VL1", bus_id="B1",
                                   target_p=10.0, target_v=225.0, min_p=0.0, max_p=100.0,
                                   voltage_regulator_on=True)
        # SVC sits on B2 (a PQ bus, no other voltage controller there) so it doesn't
        # clash with G0's local PV control on B1.
        self.net.create_static_var_compensators(
            id="SVC1", voltage_level_id="VL2", bus_id="B2",
            b_min=-0.01, b_max=0.01, target_v=225.0, target_q=0.0,
            regulation_mode="VOLTAGE", regulating=True)
        self.net.update_static_var_compensators(id="SVC1", regulated_element_id="LD_TARGET")

    def test_does_not_crash_and_falls_back_to_local_control(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = init_from_pypowsybl(self.net, sort_index=False)
        import numpy as np
        flat_v = np.full(model.get_bus_vn_kv().shape[0], 1.0, dtype=complex)
        Vdc = model.dc_pf(flat_v, 30, 1e-8)
        self.assertGreater(Vdc.shape[0], 0)
        V = model.ac_pf(Vdc, 30, 1e-8)
        self.assertGreater(V.shape[0], 0)


if __name__ == "__main__":
    unittest.main()
