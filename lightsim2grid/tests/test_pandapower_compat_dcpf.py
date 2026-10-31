# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings
import numpy as np

try:
    import pandapower as pp
    import pandapower.networks as pn
    from pandapower.pd2ppc import _pd2ppc
    from pandapower.pypower.makeBdc import makeBdc
    from pandapower.pypower.makeSbus import makeSbus
    from pandapower.pypower.dcpf import dcpf as pp_dcpf
    from pandapower.pypower.idx_bus import BUS_TYPE, REF, PV, PQ, VA, GS
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False

from lightsim2grid.pandapower_compat import dcpf


@unittest.skipUnless(PP_AVAILABLE, "pandapower is not installed")
class TestPandapowerCompatDcpf(unittest.TestCase):
    """lightsim2grid.pandapower_compat.dcpf gives the angles pandapower's own dcpf does."""

    def _inputs(self, net):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(net)
            _, ppci = _pd2ppc(net)
        bus, branch, base_mva = ppci["bus"], ppci["branch"], ppci["baseMVA"]
        B, _, Pbusinj, _ = makeBdc(bus, branch)[:4]
        Pbus = makeSbus(base_mva, bus, ppci["gen"]).real - Pbusinj - bus[:, GS] / base_mva
        ref = np.flatnonzero(bus[:, BUS_TYPE] == REF)
        pv = np.flatnonzero(bus[:, BUS_TYPE] == PV)
        pq = np.flatnonzero(bus[:, BUS_TYPE] == PQ)
        Va0 = np.deg2rad(bus[:, VA])
        return B, Pbus, Va0, ref, pv, pq

    def _check(self, net):
        B, Pbus, Va0, ref, pv, pq = self._inputs(net)
        Va_pp = pp_dcpf(B, Pbus, Va0, ref, pv, pq)
        Va_ls = dcpf(B, Pbus, Va0, ref, pv, pq)
        self.assertEqual(Va_ls.shape, Va_pp.shape)
        np.testing.assert_allclose(Va_ls, Va_pp, atol=1e-8, rtol=0.)

    def test_case14(self):
        self._check(pn.case14())

    def test_case118(self):
        self._check(pn.case118())


if __name__ == "__main__":
    unittest.main()
