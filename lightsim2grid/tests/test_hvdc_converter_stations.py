# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Converter station personalities (VSC / LCC) on a fixed-setpoint HVDC line."""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc, recv_mw  # noqa: E402


class TestHvdcConverterStations(unittest.TestCase):
    def setUp(self):
        self.b1, self.b2 = 3, 9

    def _solve(self, model, nb_bus):
        Vinit = np.ones(nb_bus, dtype=np.complex128)
        V = model.ac_pf(Vinit, 30, 1e-10)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged")
        return V

    def test_fixed_setpoint_flow_and_losses(self):
        # side 1 rectifier (converters_mode == 0): bus1 sends p_setpoint, bus2 receives
        # the post-loss power. res_p uses the generator convention (injection at bus).
        lf1, lf2, psp = 0.02, 0.03, 50.0
        net, model = make_case14_hvdc(self.b1, self.b2,
                                      loss_factor1=lf1, loss_factor2=lf2, p_setpoint=psp)
        self._solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p1_mw, -psp, places=8)
        self.assertAlmostEqual(h.res_p2_mw, recv_mw(psp, lf1, lf2), places=8)

    def test_reversed_converters_mode(self):
        # side 2 rectifier (converters_mode == 1): direction is mirrored.
        psp = 50.0
        net, model = make_case14_hvdc(self.b1, self.b2, converters_mode=1, p_setpoint=psp)
        self._solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p2_mw, -psp, places=8)
        self.assertAlmostEqual(h.res_p1_mw, recv_mw(psp, 0.0, 0.0), places=8)

    def test_vsc_regulating_is_pv(self):
        # a regulating VSC behaves like a generator: it holds its voltage setpoint.
        vm2 = 1.03
        net, model = make_case14_hvdc(self.b1, self.b2, vreg2=True, vm2=vm2, p_setpoint=50.0)
        V = self._solve(model, net.bus.shape[0])
        self.assertAlmostEqual(abs(V[self.b2]), vm2, places=8)

    def test_vsc_regulator_off_constant_q(self):
        # a non-regulating VSC injects a constant reactive power, independent of voltage.
        q2 = 15.0
        net, model = make_case14_hvdc(self.b1, self.b2, vreg2=False, q2=q2, p_setpoint=50.0)
        self._solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_q2_mvar, q2, places=8)
        self.assertAlmostEqual(h.station2.res_q_mvar, q2, places=8)

    def test_lcc_constant_pq(self):
        # an LCC station is a constant-PQ load: Q = -|P| * tan(acos(power_factor)).
        pf2, psp = 0.9, 50.0
        net, model = make_case14_hvdc(self.b1, self.b2, type2=1, power_factor2=pf2, p_setpoint=psp)
        self._solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        st2 = h.station2
        self.assertEqual(st2.converter_type, 1)  # LCC
        expected_q = -abs(h.res_p2_mw) * np.tan(np.arccos(pf2))
        self.assertAlmostEqual(h.res_q2_mvar, expected_q, places=8)
        self.assertAlmostEqual(st2.res_q_mvar, expected_q, places=8)

    def test_lcc_q_voltage_independent(self):
        # the LCC reactive consumption must not depend on the bus voltage: forcing the
        # sending VSC to a different voltage leaves Q (which only depends on |P|) untouched.
        pf2, psp = 0.92, 60.0
        q_ref = None
        for vm1 in (1.0, 1.06):
            net, model = make_case14_hvdc(self.b1, self.b2, type2=1, power_factor2=pf2,
                                          p_setpoint=psp, vreg1=True, vm1=vm1)
            self._solve(model, net.bus.shape[0])
            h = model.get_dclines()[0]
            q = h.res_q2_mvar
            if q_ref is None:
                q_ref = q
            else:
                self.assertAlmostEqual(q, q_ref, places=8)
        self.assertAlmostEqual(q_ref, -psp * np.tan(np.arccos(pf2)), places=8)

    def test_loss_factor_roundtrip(self):
        # both station loss factors compound multiplicatively on a lossless DC line.
        for lf1, lf2 in [(0.0, 0.0), (0.011, 0.0), (0.0, 0.05), (0.02, 0.04)]:
            net, model = make_case14_hvdc(self.b1, self.b2,
                                          loss_factor1=lf1, loss_factor2=lf2, p_setpoint=40.0)
            self._solve(model, net.bus.shape[0])
            h = model.get_dclines()[0]
            self.assertAlmostEqual(h.res_p2_mw, recv_mw(40.0, lf1, lf2), places=8,
                                   msg=f"loss factors {lf1}/{lf2}")


if __name__ == "__main__":
    unittest.main()
