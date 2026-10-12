# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A DC contingency analysis on several threads.

In DC a contingency is applied by editing the solver's own reduced matrix, not a copy
of Ybus. A worker thread's solver is built fresh, without that matrix, and editing it
crashed the process (segmentation fault) as soon as `nb_thread > 1`. These tests run
every outage of a small grid on several threads and compare it with one thread, and
with a plain `dc_pf` on a copy of the grid with the branch switched off.
"""

import unittest
import warnings

import numpy as np

from lightsim2grid.algorithm import AlgorithmType
from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP

try:
    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False


@unittest.skipUnless(PP_AVAILABLE, "needs pandapower")
class TestDCContingencyAnalysisThreads(unittest.TestCase):
    def setUp(self):
        warnings.filterwarnings("ignore")
        self.grid = init_from_pandapower(pn.case14())
        self.n_line = len(self.grid.get_lines())
        self.n_branch = self.n_line + len(self.grid.get_trafos())
        self.v0 = np.ones(self.grid.total_bus(), dtype=complex)
        available = ContingencyAnalysisCPP(self.grid).available_default_algorithms()
        self.algos = [a for a in (AlgorithmType.DC_SparseLU, AlgorithmType.DC_KLU) if a in available]

    def _analysis(self, algo, nb_thread, handle_disconnected_grid=False):
        sa = ContingencyAnalysisCPP(self.grid)
        sa.change_algorithm(algo)
        sa.handle_disconnected_grid = handle_disconnected_grid
        sa.add_multiple_n1(list(range(self.n_branch)))
        sa.nb_thread = nb_thread
        return sa

    @staticmethod
    def _results(sa):
        sa.compute_power_flows()
        rows = np.argsort([int(d[0]) for d in sa.my_defaults()])
        return (np.asarray(sa.converged_mask(), dtype=bool)[rows],
                np.asarray(sa.get_voltages())[rows],
                np.asarray(sa.get_power_flows())[rows])

    def _reference(self, algo, branch_id):
        grid = self.grid.copy()
        grid.change_algorithm(algo)
        if branch_id < self.n_line:
            grid.deactivate_powerline(branch_id)
        else:
            grid.deactivate_trafo(branch_id - self.n_line)
        V = grid.dc_pf(self.v0, 10, 1e-8)
        flows = np.r_[grid.get_line_res1()[0], grid.get_trafo_res1()[0]]
        return V, flows

    def test_threads_match_one_thread_and_dc_pf(self):
        for algo in self.algos:
            for hdg in (False, True):
                sa1 = self._analysis(algo, 1, hdg)
                sa1.compute(self.v0, 10, 1e-8)
                conv1, V1, P1 = self._results(sa1)
                self.assertGreater(conv1.sum(), 0)
                for nb_thread in (2, 4):
                    with self.subTest(algo=algo, handle_disconnected_grid=hdg, nb_thread=nb_thread):
                        sa = self._analysis(algo, nb_thread, hdg)
                        sa.compute(self.v0, 10, 1e-8)
                        conv, V, P = self._results(sa)
                        np.testing.assert_array_equal(conv, conv1)
                        np.testing.assert_allclose(V[conv], V1[conv], rtol=1e-10, atol=1e-10)
                        np.testing.assert_allclose(P[conv], P1[conv], rtol=1e-10, atol=1e-8)
                        # every worker's first outage is the one that went wrong: check
                        # them all against an independent solve, not just against 1 thread
                        for b in np.where(conv)[0]:
                            V_ref, P_ref = self._reference(algo, int(b))
                            if V_ref.shape[0] == 0:
                                continue   # islanding outage, only solved in handle_disconnected_grid mode
                            P_ref[b] = 0.
                            np.testing.assert_allclose(np.angle(V[b]), np.angle(V_ref), atol=1e-10)
                            np.testing.assert_allclose(P[b], P_ref, atol=1e-8)

    def test_kept_workers_give_the_same_answer(self):
        for algo in self.algos:
            for reuse in (True, False):
                with self.subTest(algo=algo, reuse_base_case=reuse):
                    sa = self._analysis(algo, 3)
                    sa.reuse_base_case = reuse
                    sa.compute(self.v0, 10, 1e-8)
                    conv, V, P = self._results(sa)
                    sa.compute(self.v0, 10, 1e-8)
                    self.assertEqual(sa.thread_algos_were_reused(), reuse)
                    conv2, V2, P2 = self._results(sa)
                    np.testing.assert_array_equal(conv2, conv)
                    np.testing.assert_allclose(V2[conv], V[conv], rtol=1e-10, atol=1e-10)
                    np.testing.assert_allclose(P2[conv], P[conv], rtol=1e-10, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
