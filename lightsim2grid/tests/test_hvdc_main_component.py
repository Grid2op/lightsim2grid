# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""HVDC line across the boundary of ``consider_only_main_component``.

An HVDC line bridges two AC grids that are NOT synchronous: ``get_graph`` adds no
edge for it, so the connectivity BFS of ``consider_only_main_component`` really
splits the grid into *synchronous* components. When the two converters land in
different components only one of them is in the main (solved) one.

OpenLoadFlow keeps that in-main converter as a fixed boundary injection (it still
imports / exports its scheduled HVDC power). lightsim2grid used to deactivate the
WHOLE line as soon as one side was out of the main component, silently dropping
that injection (hundreds of MW on real RTE grids -- see HVDC_OLF_FINDINGS.md).
These tests pin the fixed behaviour: keep the line connected, keep the in-main
converter injecting, open only the out-of-main converter.
"""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc  # noqa: E402

# bus 7 of case14 is a leaf (degree 1, no load): isolating it leaves the rest of the
# grid connected and well-posed, so it is the perfect "other synchronous component".
_LEAF_BUS = 7


class TestHvdcMainComponent(unittest.TestCase):
    def _isolate_leaf(self, model):
        """Open every branch touching ``_LEAF_BUS`` so the bus (and the converter
        wired to it) becomes its own connected component."""
        for line in model.get_lines():
            if _LEAF_BUS in (line.bus1_id, line.bus2_id):
                model.deactivate_powerline(line.id)
        for trafo in model.get_trafos():
            if _LEAF_BUS in (trafo.bus1_id, trafo.bus2_id):
                model.deactivate_trafo(trafo.id)

    def test_one_converter_out_of_main_keeps_injection(self):
        # side 1 (bus 3, rectifier) is in the main component, side 2 (bus 7) is
        # isolated. The line must stay connected and side 1 must keep injecting.
        psp = 30.0
        net, model = make_case14_hvdc(3, _LEAF_BUS, converters_mode=0, p_setpoint=psp)
        self._isolate_leaf(model)
        model.consider_only_main_component()

        h = model.get_dclines()[0]
        self.assertTrue(h.connected_global, "the HVDC line must stay connected_global")
        self.assertTrue(h.connected1, "the in-main converter (side 1) must stay connected")
        self.assertFalse(h.connected2, "the out-of-main converter (side 2) must be opened")

        V = model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-10)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged after isolating the HVDC far end")
        # the in-main rectifier still draws its full setpoint (generator convention),
        # exactly as if the line were fully connected: the injection is preserved.
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p1_mw, -psp, places=6)

    def test_both_converters_in_main_unchanged(self):
        # control case: both converters in the main component -> nothing is opened.
        net, model = make_case14_hvdc(3, 9, converters_mode=0, p_setpoint=30.0)
        model.consider_only_main_component()
        h = model.get_dclines()[0]
        self.assertTrue(h.connected_global)
        self.assertTrue(h.connected1)
        self.assertTrue(h.connected2)


if __name__ == "__main__":
    unittest.main()
