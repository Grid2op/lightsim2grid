# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Angle-droop ("AC emulation") HVDC in the Newton-Raphson solver."""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _aux_make_hvdc import make_case14_hvdc, recv_mw, droop_raw_mw, MW_PER_DEG_TO_MW_PER_RAD  # noqa: E402


class TestHvdcDroop(unittest.TestCase):
    def setUp(self):
        # buses 3 and 9 are NOT adjacent in case14: a good probe for the structural
        # cross-bus Jacobian entries the droop term must inject.
        self.b1, self.b2 = 3, 9
        self.p0 = 10.0
        self.droop = 5.0      # MW / deg
        self.lf = 0.011
        self.pmax = 300.0

    def _model(self, **kw):
        b1 = kw.pop("b1", self.b1)
        b2 = kw.pop("b2", self.b2)
        defaults = dict(loss_factor1=self.lf, loss_factor2=self.lf,
                        droop_enabled=True, p0=self.p0, droop_mw_per_deg=self.droop,
                        pmax12=self.pmax, pmax21=self.pmax, p_setpoint=20.0)
        defaults.update(kw)
        return make_case14_hvdc(b1, b2, **defaults)

    def _solve(self, model, nb_bus, tol=1e-10):
        V = model.ac_pf(np.ones(nb_bus, dtype=np.complex128), 30, tol)
        self.assertGreater(V.shape[0], 0, "ac_pf diverged with droop hvdc")
        return V

    def _expected_flows(self, raw):
        if raw >= 0:
            return -raw, recv_mw(raw, self.lf, self.lf)
        return recv_mw(-raw, self.lf, self.lf), raw

    def test_linear_flow_matches_closed_form(self):
        net, model = self._model()
        self._solve(model, net.bus.shape[0])
        h = model.get_dclines()[0]
        raw = droop_raw_mw(self.p0, self.droop, h.res_theta1_deg, h.res_theta2_deg)
        p1_exp, p2_exp = self._expected_flows(raw)
        self.assertAlmostEqual(h.res_p1_mw, p1_exp, places=8)
        self.assertAlmostEqual(h.res_p2_mw, p2_exp, places=8)
        self.assertEqual(model.get_status_droop_hvdc(0), 0)

    def test_kirchhoff_balance(self):
        # the converged solution must balance once the droop flows are accounted for.
        net, model = self._model()
        V = self._solve(model, net.bus.shape[0])
        mismatch = model.check_solution(V, False)
        self.assertLess(np.abs(np.real(mismatch)).max(), 1e-6)
        self.assertLess(np.abs(np.imag(mismatch)).max(), 1e-6)

    def test_structural_cross_entries(self):
        # buses 3 and 9 are not adjacent: without the droop FeatureSink entries the
        # (p_row(b1), theta_col(b2)) slots would be absent and NR would crawl/diverge.
        net, model = self._model()
        self._solve(model, net.bus.shape[0])
        solver = model.get_solver()
        J = solver.get_J()
        theta_cols = solver.get_theta_to_J_col()
        c1, c2 = theta_cols[self.b1], theta_cols[self.b2]
        self.assertGreaterEqual(c1, 0)
        self.assertGreaterEqual(c2, 0)
        Jd = J.todense()
        # the off-diagonal coupling between the two HVDC ends must be present
        self.assertNotEqual(Jd[c1, c2], 0.0)
        self.assertNotEqual(Jd[c2, c1], 0.0)

    def test_jacobian_quadratic_convergence(self):
        # an exact analytic Jacobian => Newton converges in a handful of iterations.
        # a wrong droop Jacobian would still reach the same fixed point but only after
        # many more iterations (or diverge): assert a tight iteration budget.
        net, model = self._model()
        self._solve(model, net.bus.shape[0])
        nb_iter = model.get_solver().get_nb_iter()
        self.assertLessEqual(nb_iter, 6, f"too many NR iterations ({nb_iter})")

    def test_saturated_regime_pins_flow(self):
        # force the +1 saturation: the flow is pinned at pmax_1to2, the receive side
        # gets the post-loss power, and the four droop Jacobian entries become 0.
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        self._solve(model, nb_bus)
        model.set_status_droop_hvdc(0, 1)
        self.assertEqual(model.get_status_droop_hvdc(0), 1)
        self._solve(model, nb_bus)
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p1_mw, -self.pmax, places=8)
        self.assertAlmostEqual(h.res_p2_mw, recv_mw(self.pmax, self.lf, self.lf), places=8)

    def test_saturated_jacobian_entries_zero(self):
        # buses 3 and 9 are not adjacent, so their only Jacobian coupling comes from the
        # droop term. In the linear regime the cross entries are the droop slopes; in the
        # saturated regime the injection is constant => those slopes (and hence the cross
        # entries) are exactly zero, even though the sparsity slots are kept.
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        self._solve(model, nb_bus)
        solver = model.get_solver()
        theta_cols = solver.get_theta_to_J_col()
        c1, c2 = theta_cols[self.b1], theta_cols[self.b2]
        J_lin = np.asarray(solver.get_J().todense())
        self.assertGreater(abs(J_lin[c1, c2]) + abs(J_lin[c2, c1]), 0.0)

        model.set_status_droop_hvdc(0, 1)
        self._solve(model, nb_bus)
        J_sat = np.asarray(model.get_solver().get_J().todense())
        self.assertEqual(J_sat[c1, c2], 0.0)
        self.assertEqual(J_sat[c2, c1], 0.0)

    def test_pattern_stable_across_status_flips(self):
        # the KLU symbolic factorization is reused across status flips: declaring the
        # droop entries regardless of regime keeps J's sparsity pattern (nnz) constant.
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        self._solve(model, nb_bus)
        nnz_lin = model.get_solver().get_J().nnz
        model.set_status_droop_hvdc(0, 1)
        self._solve(model, nb_bus)
        nnz_sat = model.get_solver().get_J().nnz
        self.assertEqual(nnz_lin, nnz_sat)

    def test_flip_back_consistency(self):
        # linear -> saturated -> linear returns to the original solution bit-for-bit.
        net, model = self._model()
        nb_bus = net.bus.shape[0]
        V0 = self._solve(model, nb_bus).copy()
        model.set_status_droop_hvdc(0, 1)
        self._solve(model, nb_bus)
        model.set_status_droop_hvdc(0, 0)
        V2 = self._solve(model, nb_bus)
        self.assertLess(np.abs(V2 - V0).max(), 1e-8)

    def test_saturated_reverse_direction(self):
        # -1 saturation pins the reverse flow at pmax_2to1. A full reverse 300 MW is not
        # feasible on these two case14 buses, so cap the reverse direction at 50 MW and
        # warm-start from the linear solution.
        pmax21 = 50.0
        net, model = self._model(pmax21=pmax21)
        nb_bus = net.bus.shape[0]
        V = self._solve(model, nb_bus)
        model.set_status_droop_hvdc(0, -1)
        V2 = model.ac_pf(V.copy(), 30, 1e-10)
        self.assertGreater(V2.shape[0], 0, "reverse-saturated ac_pf diverged")
        h = model.get_dclines()[0]
        self.assertAlmostEqual(h.res_p2_mw, -pmax21, places=8)
        self.assertAlmostEqual(h.res_p1_mw, recv_mw(pmax21, self.lf, self.lf), places=8)

    def test_slack_end_droop_still_balances(self):
        # one HVDC end sitting on the slack bus: its angle is the fixed reference, the
        # droop is still evaluated from it (no theta unknown there), and the network
        # must balance. case14's slack is bus 0.
        net, model = self._model(b1=0, b2=9)
        nb_bus = net.bus.shape[0]
        V = self._solve(model, nb_bus)
        h = model.get_dclines()[0]
        raw = droop_raw_mw(self.p0, self.droop, h.res_theta1_deg, h.res_theta2_deg)
        p1_exp, p2_exp = self._expected_flows(raw)
        self.assertAlmostEqual(h.res_p1_mw, p1_exp, places=6)
        self.assertAlmostEqual(h.res_p2_mw, p2_exp, places=6)
        mismatch = model.check_solution(V, False)
        self.assertLess(np.abs(np.real(mismatch)).max(), 1e-6)

    def test_fdpf_guard_raises(self):
        # the droop term only lives in the NR system: selecting a non-NR algorithm with
        # an active droop line must be rejected rather than silently ignored.
        from lightsim2grid.solver import AlgorithmType
        net, model = self._model()
        model.change_algorithm(AlgorithmType.GaussSeidel)
        with self.assertRaises(Exception):
            model.ac_pf(np.ones(net.bus.shape[0], dtype=np.complex128), 30, 1e-8)


if __name__ == "__main__":
    unittest.main()
