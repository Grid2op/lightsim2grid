# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Every connected N-1 row of a grid with a remote voltage regulator converges, whatever
``nb_thread`` (the bordered voltage rows put zeros on the diagonal, so a KLU
refactorization could hit a zero pivot or reuse a poor one from another row)."""

import unittest
import numpy as np

try:
    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestContingencyAnalysisRemoteRegulator(unittest.TestCase):
    def setUp(self):
        net = pp.network.create_ieee14()
        loads = net.get_loads(attributes=["voltage_level_id"])
        load = loads.index[loads["voltage_level_id"] == "VL12"][0]
        net.update_generators(id="B6-G", regulated_element_id=load, target_v=14.2,
                              min_q=-9999., max_q=9999.)
        self.grid = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=True)
        self.V0 = self.grid.ac_pf(np.ones(self.grid.total_bus(), dtype=complex), 30, 1e-11)
        assert self.V0.shape[0] > 0

    def _single_solves(self):
        n_lines = len(self.grid.get_lines())
        res = {}
        for b in range(n_lines + len(self.grid.get_trafos())):
            g = self.grid.copy()
            if b < n_lines:
                g.deactivate_powerline(b)
            else:
                g.deactivate_trafo(b - n_lines)
            res[b] = g.ac_pf(self.V0.copy(), 30, 1e-11).shape[0] > 0
        return res

    def test_connected_rows_converge(self):
        single = self._single_solves()
        for nb_thread in (1, 2, 4, 8):
            with self.subTest(nb_thread=nb_thread):
                ca = ContingencyAnalysisCPP(self.grid)
                ca.nb_thread = nb_thread
                ca.add_all_n1()
                ca.compute(self.V0.copy(), 30, 1e-11)
                connected = np.asarray(ca.is_grid_connected_after_contingency()).ravel().astype(bool)
                conv = np.asarray(ca.converged_mask()).ravel()
                for b, ok in single.items():
                    if connected[b] and ok:
                        self.assertTrue(conv[b], f"branch {b} did not converge")

    def test_single_contingency_alone(self):
        ca = ContingencyAnalysisCPP(self.grid)
        ca.add_n1(5)
        ca.compute(self.V0.copy(), 30, 1e-11)
        self.assertTrue(np.asarray(ca.converged_mask()).all())


if __name__ == "__main__":
    unittest.main()
