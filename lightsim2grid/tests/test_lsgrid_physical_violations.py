# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``LSGrid.get_physical_violations``: the batch algorithms' ``compute_physical_violations``
checks, on the grid's own last powerflow. It must report exactly what a batch reports on
its base case, and what a hand-made comparison of the results against the limits says.
"""

import unittest
import warnings
import numpy as np
import pandapower.networks as pn

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from lightsim2grid.gridmodel import init_from_pandapower
    from lightsim2grid.lightsim2grid_cpp import ContingencyAnalysisCPP
    from lightsim2grid.algorithm import AlgorithmType

_MAX_IT = 30
_TOL = 1e-10


def _key(v):
    return (str(v.element_type), int(v.element_id), str(v.violation_type))


class TestLSGridPhysicalViolations(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.net = pn.case14()
            self.grid = init_from_pandapower(self.net)
        gens = self.grid.get_generators()
        self.n_gen = len(gens)
        for g in gens:
            self.grid.add_gen_slackbus(g.id, 1.)
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)

    def _batch_n_case(self, dc=False):
        """the same checks, through a batch (one contingency, whatever it is) base case"""
        ca = ContingencyAnalysisCPP(self.grid)
        if dc:
            ca.change_algorithm(AlgorithmType.DC_SparseLU)
        ca.add_n1(0)
        ca.compute_physical_violations = True
        ca.compute(1. * self.V0, _MAX_IT, _TOL)
        return ca.get_physical_violations_n()

    def test_raises_before_any_powerflow(self):
        with self.assertRaises(RuntimeError):
            self.grid.get_physical_violations()

    def test_no_limits_no_violation(self):
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        # case14's generators regulate within their Q limits: nothing to report
        viol = [v for v in self.grid.get_physical_violations() if "_P" in str(v.violation_type)]
        self.assertEqual(viol, [])

    def test_gen_p_limits_match_hand_check_and_batch(self):
        # tight p limits around the targets: the slack share pushes some past them
        targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        max_p = targets + 0.5
        min_p = targets - 0.5
        self.grid.set_gen_p_limits(min_p, max_p)
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        viol = self.grid.get_physical_violations()
        p_viol = [v for v in viol if str(v.element_type).endswith("GENERATOR")]
        self.assertGreater(len(p_viol), 0, "this test wants a generator past its limit")
        # by hand: res_p vs the limits
        res_p = np.array([g.res_p_mw for g in self.grid.get_generators()])
        expected = {}
        for g in range(self.n_gen):
            if res_p[g] > max_p[g] + 1e-4:
                expected[g] = ("HIGH_P", max_p[g])
            elif res_p[g] < min_p[g] - 1e-4:
                expected[g] = ("LOW_P", min_p[g])
        self.assertEqual({v.element_id for v in p_viol}, set(expected))
        for v in p_viol:
            kind, limit = expected[v.element_id]
            self.assertIn(kind, str(v.violation_type))
            self.assertAlmostEqual(v.limit, limit, places=9)
            self.assertAlmostEqual(v.value, res_p[v.element_id], places=6)
        # ... and the batch says the same on its base case
        batch = self._batch_n_case()
        self.assertEqual(sorted(_key(v) for v in viol), sorted(_key(v) for v in batch))
        for a, b in zip(sorted(viol, key=_key), sorted(batch, key=_key)):
            self.assertAlmostEqual(a.value, b.value, places=9)

    def test_dc(self):
        targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        self.grid.set_gen_p_limits(targets - 0.5, targets + 0.5)
        with self.assertRaises(RuntimeError):
            self.grid.get_physical_violations(ac=False)  # no DC powerflow yet
        V = self.grid.dc_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        viol = self.grid.get_physical_violations(ac=False)
        res_p = np.array([g.res_p_mw for g in self.grid.get_generators()])
        expected = {g for g in range(self.n_gen) if abs(res_p[g] - targets[g]) > 0.5 + 1e-4}
        self.assertEqual({v.element_id for v in viol}, expected)
        batch = self._batch_n_case(dc=True)
        self.assertEqual(sorted(_key(v) for v in viol), sorted(_key(v) for v in batch))

    def test_tolerance(self):
        targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        self.grid.set_gen_p_limits(targets - 0.5, targets + 0.5)
        self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(len(self.grid.get_physical_violations()), 0)
        # a huge tolerance forgives everything
        self.assertEqual(self.grid.get_physical_violations(tol_mva=1e3), [])

    def test_saturated_units_after_redistribution_not_reported(self):
        # a unit the slack pre-pass clamped at max_p left the slack: it sits exactly at
        # its bound and is not reported (it takes no share of the losses any more)
        targets = np.array([g.target_p_mw for g in self.grid.get_generators()])
        max_p = np.full(self.n_gen, np.inf)
        max_p[1] = targets[1] + 1.
        self.grid.set_gen_p_limits(np.full(self.n_gen, -np.inf), max_p)
        report = self.grid.redistribute_active_power(30.)
        self.assertEqual(report.nb_saturated, 1)
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual([v for v in self.grid.get_physical_violations() if "_P" in str(v.violation_type)], [])
        self.assertAlmostEqual(self.grid.get_generators()[1].res_p_mw, max_p[1], places=6)


if __name__ == "__main__":
    unittest.main()
