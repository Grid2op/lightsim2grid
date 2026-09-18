# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Pickle / copy round-trip of a grid carrying converter stations + an HVDC line."""

import os
import sys
import pickle
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc  # noqa: E402


class TestHvdcPickle(unittest.TestCase):
    def _grids(self):
        # a fixed LCC-terminated line and a VSC-VSC angle-droop line
        yield "lcc_fixed", make_case14_hvdc(3, 9, type2=1, power_factor2=0.9,
                                            loss_factor1=0.02, loss_factor2=0.03, p_setpoint=50.0)
        yield "vsc_droop", make_case14_hvdc(3, 9, loss_factor1=0.011, loss_factor2=0.011,
                                            droop_enabled=True, p0=10.0, droop_mw_per_deg=5.0,
                                            p_setpoint=20.0)

    def _assert_same_solve(self, ref_model, new_model, net):
        nb = net.bus.shape[0]
        Vref = ref_model.ac_pf(np.ones(nb, dtype=np.complex128), 30, 1e-10)
        Vnew = new_model.ac_pf(np.ones(nb, dtype=np.complex128), 30, 1e-10)
        self.assertGreater(Vnew.shape[0], 0)
        self.assertLess(np.abs(Vref - Vnew).max(), 1e-12)

    def _assert_line_preserved(self, ref_model, new_model):
        h0 = ref_model.get_dclines()[0]
        h1 = new_model.get_dclines()[0]
        self.assertEqual(h1.droop_enabled, h0.droop_enabled)
        self.assertAlmostEqual(h1.droop_p0_mw, h0.droop_p0_mw)
        self.assertAlmostEqual(h1.droop_k_mw_per_rad, h0.droop_k_mw_per_rad)
        self.assertAlmostEqual(h1.p_setpoint_mw, h0.p_setpoint_mw)
        self.assertEqual(h1.converters_mode, h0.converters_mode)
        self.assertEqual(h1.station1.converter_type, h0.station1.converter_type)
        self.assertEqual(h1.station2.converter_type, h0.station2.converter_type)
        self.assertAlmostEqual(h1.station2.power_factor, h0.station2.power_factor)

    def test_pickle_roundtrip(self):
        for name, (net, model) in self._grids():
            with self.subTest(grid=name):
                model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)
                restored = pickle.loads(pickle.dumps(model))
                self._assert_line_preserved(model, restored)
                self._assert_same_solve(model, restored, net)

    def test_copy_roundtrip(self):
        for name, (net, model) in self._grids():
            with self.subTest(grid=name):
                model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)
                cpy = model.copy()
                self._assert_line_preserved(model, cpy)
                self._assert_same_solve(model, cpy, net)


if __name__ == "__main__":
    unittest.main()
