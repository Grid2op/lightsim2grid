# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""``LSGrid.get_violations``: the batch algorithms' ``compute_limit_violations`` (operational
limits: bus voltages, branch currents on both sides) on the grid's own last powerflow. It
must report exactly what a batch reports on its base case, and what a hand-made comparison
of the results against the limits says.
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
    return (str(v.element_type), int(v.element_id), int(v.side), str(v.violation_type))


class TestLSGridViolations(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.net = pn.case14()
            self.grid = init_from_pandapower(self.net)
        self.n_bus = self.grid.get_bus_vn_kv().shape[0]
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())
        self.V0 = np.ones(self.n_bus, dtype=complex)
        self.vn = np.asarray(self.grid.get_bus_vn_kv())

    def _set_limits(self, v_margin=0.03, i_ratio=0.8):
        """voltage limits at +/- v_margin around nominal, current limits at i_ratio of the
        base-case currents (so that some sides are above their limit right away)"""
        self.grid.set_bus_voltage_limits(self.vn * (1. - v_margin), self.vn * (1. + v_margin))
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        assert V.shape[0] > 0
        i1 = np.array([l.res_a1_ka for l in self.grid.get_lines()]) * 1000.
        i2 = np.array([l.res_a2_ka for l in self.grid.get_lines()]) * 1000.
        self.line_lim1, self.line_lim2 = i_ratio * i1 / 1000., i_ratio * i2 / 1000.  # kA
        self.grid.set_line_current_limit_side1(self.line_lim1)
        self.grid.set_line_current_limit_side2(self.line_lim2)
        t1 = np.array([t.res_a1_ka for t in self.grid.get_trafos()]) * 1000.
        t2 = np.array([t.res_a2_ka for t in self.grid.get_trafos()]) * 1000.
        self.trafo_lim1, self.trafo_lim2 = 1.2 * t1 / 1000., 1.2 * t2 / 1000.  # never violated
        self.grid.set_trafo_current_limit_side1(self.trafo_lim1)
        self.grid.set_trafo_current_limit_side2(self.trafo_lim2)

    def _batch_n_case(self, threshold=1., dc=False):
        ca = ContingencyAnalysisCPP(self.grid, True)  # compute_limit_violations first
        if dc:
            ca.change_algorithm(AlgorithmType.DC_SparseLU)
        ca.add_n1(0)
        ca.violation_threshold = threshold
        ca.compute(1. * self.V0, _MAX_IT, _TOL)
        return ca.get_violations_n()

    def _hand_check(self, threshold=1.):
        """what the results say, compared to the limits by hand"""
        expected = set()
        vm = np.abs(self.grid.get_V()) * self.vn  # get_V is in solver order == grid order here (all buses solved)
        vmin, vmax = self.vn * 0.97, self.vn * 1.03
        for b in range(self.n_bus):
            anchor = min(max(self.vn[b], vmin[b]), vmax[b])
            low = vmin[b] + (1. - threshold) * (anchor - vmin[b])
            high = vmax[b] - (1. - threshold) * (vmax[b] - anchor)
            if vm[b] <= low:
                expected.add(("ViolationElementType.BUS", b, 0, "LimitViolationType.LOW_VOLTAGE"))
            elif vm[b] >= high:
                expected.add(("ViolationElementType.BUS", b, 0, "LimitViolationType.HIGH_VOLTAGE"))
        for l in self.grid.get_lines():
            if l.res_a1_ka >= threshold * self.line_lim1[l.id]:
                expected.add(("ViolationElementType.LINE", l.id, 1, "LimitViolationType.CURRENT"))
            if l.res_a2_ka >= threshold * self.line_lim2[l.id]:
                expected.add(("ViolationElementType.LINE", l.id, 2, "LimitViolationType.CURRENT"))
        return expected

    def test_raises_before_any_powerflow(self):
        with self.assertRaises(RuntimeError):
            self.grid.get_violations()

    def test_no_limits_no_violation(self):
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(self.grid.get_violations(), [])

    def test_bad_threshold(self):
        self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        for thr in (0., -1., 1.5):
            with self.assertRaises(RuntimeError):
                self.grid.get_violations(thr)

    def test_matches_hand_check_and_batch(self):
        self._set_limits()
        V = self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        viol = self.grid.get_violations()
        keys = {_key(v) for v in viol}
        self.assertGreater(sum(1 for k in keys if "CURRENT" in k[3]), 0, "this test wants current violations")
        self.assertGreater(sum(1 for k in keys if "VOLTAGE" in k[3]), 0, "this test wants voltage violations")
        self.assertEqual(keys, self._hand_check())
        # values and limits: what the grid computed
        for v in viol:
            if "CURRENT" in str(v.violation_type):
                el = self.grid.get_lines()[v.element_id]
                # currents are reported in kA, like the limits they are compared to
                self.assertAlmostEqual(v.value, (el.res_a1_ka if v.side == 1 else el.res_a2_ka), places=9)
                self.assertAlmostEqual(v.limit, (self.line_lim1 if v.side == 1 else self.line_lim2)[v.element_id], places=12)
            else:
                self.assertAlmostEqual(v.value, abs(self.grid.get_V()[v.element_id]) * self.vn[v.element_id], places=9)
        # ... and the batch says the same on its base case
        batch = self._batch_n_case()
        self.assertEqual(sorted(_key(v) for v in viol), sorted(_key(v) for v in batch))
        for a, b in zip(sorted(viol, key=_key), sorted(batch, key=_key)):
            self.assertAlmostEqual(a.value, b.value, places=9)
            self.assertAlmostEqual(a.limit, b.limit, places=9)

    def test_threshold(self):
        self._set_limits(i_ratio=1.05)  # every current just BELOW its limit
        self.grid.ac_pf(1. * self.V0, _MAX_IT, _TOL)
        current_at_1 = [v for v in self.grid.get_violations(1.) if "CURRENT" in str(v.violation_type)]
        self.assertEqual(current_at_1, [])
        current_at_09 = [v for v in self.grid.get_violations(0.9) if "CURRENT" in str(v.violation_type)]
        self.assertGreater(len(current_at_09), 0)
        self.assertEqual({_key(v) for v in self.grid.get_violations(0.9)}, self._hand_check(0.9))
        self.assertEqual(sorted(_key(v) for v in self.grid.get_violations(0.9)),
                         sorted(_key(v) for v in self._batch_n_case(0.9)))

    def test_dc(self):
        self._set_limits()
        with self.assertRaises(RuntimeError):
            self.grid.get_violations(ac=False)  # no DC powerflow yet
        V = self.grid.dc_pf(1. * self.V0, _MAX_IT, _TOL)
        self.assertGreater(V.shape[0], 0)
        viol = self.grid.get_violations(ac=False)
        # DC: the currents come from P only, the magnitudes are the set-points / the seed
        batch = self._batch_n_case(dc=True)
        self.assertEqual(sorted(_key(v) for v in viol), sorted(_key(v) for v in batch))
        for a, b in zip(sorted(viol, key=_key), sorted(batch, key=_key)):
            self.assertAlmostEqual(a.value, b.value, places=9)

    def test_divergence_sentinel(self):
        self._set_limits()
        # an absurd load makes the AC powerflow diverge
        self.grid.change_p_load(0, 1e5)
        V = self.grid.ac_pf(1. * self.V0, 5, _TOL)
        self.assertEqual(V.shape[0], 0, "this test needs a diverging powerflow")
        viol = self.grid.get_violations()
        self.assertEqual(len(viol), 1)
        self.assertIn("GRID", str(viol[0].element_type))
        self.assertIn("DIVERGENCE", str(viol[0].violation_type))


if __name__ == "__main__":
    unittest.main()
