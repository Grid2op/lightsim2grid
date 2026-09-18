# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Angle-droop HVDC in the DC powerflow (B-matrix + RHS terms)."""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc, recv_mw, droop_raw_mw  # noqa: E402


class TestHvdcDc(unittest.TestCase):
    def setUp(self):
        self.b1, self.b2 = 3, 9
        self.p0 = 10.0
        self.droop = 5.0
        self.lf = 0.011
        self.pmax = 300.0

    def _model(self, **kw):
        defaults = dict(loss_factor1=self.lf, loss_factor2=self.lf,
                        droop_enabled=True, p0=self.p0, droop_mw_per_deg=self.droop,
                        pmax12=self.pmax, pmax21=self.pmax, p_setpoint=20.0)
        defaults.update(kw)
        return make_case14_hvdc(self.b1, self.b2, **defaults)

    def _dc_solve(self, model, nb_bus):
        V = model.dc_pf(np.ones(nb_bus, dtype=np.complex128), 1, 1.0)
        self.assertGreater(V.shape[0], 0, "dc_pf diverged with droop hvdc")
        return V

    def test_dc_linear_droop_flow(self):
        # DC is lossless: p1 = -raw, p2 = +raw (no receive-side loss in DC, mirroring
        # _pyloadflow dc.py).
        net, model = self._model()
        self._dc_solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        raw = droop_raw_mw(self.p0, self.droop, h.res_theta1_deg, h.res_theta2_deg)
        self.assertAlmostEqual(h.res_p1_mw, -raw, places=8)
        self.assertAlmostEqual(h.res_p2_mw, raw, places=8)

    def test_dc_rhs_p0_term(self):
        # with a zero slope (k == 0) the droop reduces to a pure RHS injection of +/-p0:
        # the DC line carries exactly p0 regardless of the angles.
        net, model = self._model(droop_mw_per_deg=0.0, p0=15.0)
        self._dc_solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p1_mw, -15.0, places=8)
        self.assertAlmostEqual(h.res_p2_mw, 15.0, places=8)

    def test_dc_b_term_couples_angles(self):
        # with p0 == 0 the flow is driven purely by the k * (theta1 - theta2) B-matrix
        # term: the solved flow must equal that closed form (and be non-zero).
        net, model = self._model(p0=0.0, droop_mw_per_deg=5.0)
        self._dc_solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        raw = droop_raw_mw(0.0, 5.0, h.res_theta1_deg, h.res_theta2_deg)
        self.assertGreater(abs(raw), 1.0)  # the two ends are genuinely coupled
        self.assertAlmostEqual(h.res_p1_mw, -raw, places=8)
        self.assertAlmostEqual(h.res_p2_mw, raw, places=8)

    def test_dc_saturated_is_fixed_injection(self):
        # saturated lines arrive through Sbus as a fixed injection: bus1 sends pmax, bus2
        # receives the post-loss power (the DC linear term is lossless, but a saturated
        # line behaves like a fixed-setpoint line, which does carry the converter losses).
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        self._dc_solve(model, nb_bus)
        model.set_status_droop_hvdc(0, 1)
        self._dc_solve(model, nb_bus)
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p1_mw, -self.pmax, places=8)
        self.assertAlmostEqual(h.res_p2_mw, recv_mw(self.pmax, self.lf, self.lf), places=8)

    def test_dc_flip_back_consistency(self):
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        V0 = self._dc_solve(model, nb_bus).copy()
        model.set_status_droop_hvdc(0, 1)
        self._dc_solve(model, nb_bus)
        model.set_status_droop_hvdc(0, 0)
        V2 = self._dc_solve(model, nb_bus)
        self.assertLess(np.abs(V2 - V0).max(), 1e-10)


if __name__ == "__main__":
    unittest.main()
