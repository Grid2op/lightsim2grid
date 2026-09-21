# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""ContingencyAnalysisCPP has the same C++ class behind it as the three sweeps.

Its python interface was spelled out by hand rather than sharing the sweeps' binder,
so every feature added to that binder silently skipped it: base-case reuse and the
whole reverse-mode block existed in C++ and could not be reached from python at all.

The first test here is the guard against that happening again -- it compares the two
interfaces and insists the only differences are the ones the C++ really does gate
off. The rest check that what is now reachable actually works, rather than merely
being callable.
"""

import unittest
import warnings

import numpy as np

from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP, InjectionSweepCPP,
                                             ScenarioSweepCPP)

try:
    import pandapower.networks as pn
    from lightsim2grid.network import init_from_pandapower
    PP_AVAILABLE = True
except ImportError:
    PP_AVAILABLE = False

# The six the C++ genuinely gates off for this instantiation (SbusPolicy::NOOP): a
# contingency analysis varies the topology, not the injections, so it has no per-row
# setter and no aggregate status. Everything else the sweeps expose, it must expose.
VARY_ONLY = {"modify_gen_p", "modify_sgen_p", "modify_load_p", "modify_load_q",
             "modify_gen_v", "get_status"}


@unittest.skipUnless(PP_AVAILABLE, "needs pandapower")
class TestContingencyAnalysisBindings(unittest.TestCase):
    N_CONT = 6

    def setUp(self):
        warnings.filterwarnings("ignore")

    @staticmethod
    def _grid(load_scale=1.0, gen_v_delta=0.0):
        net = pn.case14()
        grid = init_from_pandapower(net)
        if load_scale != 1.0:
            for i, p in enumerate(np.asarray(grid.get_load_target_p())):
                grid.change_p_load(i, float(p * load_scale))
        if gen_v_delta != 0.0:
            for g in grid.get_generators():
                grid.change_v_gen(g.id, float(g.target_vm_pu + gen_v_delta))
        return grid

    def _analysis(self, grid, keep_jacobian=False):
        sa = ContingencyAnalysisCPP(grid)
        for c in range(self.N_CONT):
            sa.add_n1(c)
        if keep_jacobian:
            sa.keep_jacobian = True
        return sa

    @staticmethod
    def _v0(grid):
        return np.ones(grid.total_bus(), dtype=complex)

    # ------------------------------------------------------------------ the guard
    def test_it_exposes_everything_the_sweeps_do(self):
        def api(cls):
            return {x for x in dir(cls) if not x.startswith("_")}
        shared = api(InjectionSweepCPP) & api(ScenarioSweepCPP)
        missing = shared - api(ContingencyAnalysisCPP)
        self.assertEqual(
            missing, VARY_ONLY,
            "ContingencyAnalysisCPP's interface drifted from the sweeps'. Anything "
            "shared belongs in bind_batch_shared(), which both call; only the "
            "SbusPolicy::Vary members may be missing here.")

    # --------------------------------------------------------- base-case reuse
    def test_base_case_reuse(self):
        grid = self._grid()
        sa = self._analysis(grid)
        self.assertTrue(sa.reuse_base_case)          # on by default

        sa.compute(self._v0(grid), 30, 1e-10)
        self.assertFalse(sa.base_case_was_reused())  # nothing to reuse yet
        first = np.array(sa.get_voltages())
        analyze_after_first = sa.get_linear_solver_stats().nb_analyze
        self.assertGreaterEqual(analyze_after_first, 1)

        sa.compute(self._v0(grid), 30, 1e-10)
        self.assertTrue(sa.base_case_was_reused())
        # the symbolic factorization of the first call is still the one in use
        self.assertEqual(sa.get_linear_solver_stats().nb_analyze, analyze_after_first)
        # not bit-identical: the kept base case maps this call's Vinit onto the kept
        # labelling instead of rebuilding it, which lands a ulp away
        np.testing.assert_allclose(np.array(sa.get_voltages()), first, rtol=1e-12, atol=1e-12)

    def test_registering_a_contingency_drops_the_base_case(self):
        grid = self._grid()
        sa = self._analysis(grid)
        sa.compute(self._v0(grid), 30, 1e-10)
        sa.compute(self._v0(grid), 30, 1e-10)
        self.assertTrue(sa.base_case_was_reused())

        sa.add_n1(self.N_CONT)   # what the graph walk settled no longer describes this batch
        sa.compute(self._v0(grid), 30, 1e-10)
        self.assertFalse(sa.base_case_was_reused())

        fresh = self._analysis(self._grid())
        fresh.add_n1(self.N_CONT)
        fresh.compute(self._v0(grid), 30, 1e-10)
        np.testing.assert_allclose(np.array(sa.get_voltages()),
                                   np.array(fresh.get_voltages()), rtol=1e-10, atol=1e-10)

    def test_reuse_can_be_turned_off(self):
        grid = self._grid()
        sa = self._analysis(grid)
        sa.reuse_base_case = False
        sa.compute(self._v0(grid), 30, 1e-10)
        first = np.array(sa.get_voltages())
        n0 = sa.get_linear_solver_stats().nb_analyze
        sa.compute(self._v0(grid), 30, 1e-10)
        self.assertFalse(sa.base_case_was_reused())
        self.assertGreater(sa.get_linear_solver_stats().nb_analyze, n0)
        np.testing.assert_allclose(np.array(sa.get_voltages()), first, rtol=1e-10, atol=1e-10)

    # ------------------------------------------------------------- the adjoint
    @staticmethod
    def _loss(V):
        n_bus = V.shape[1]
        wr = 1. + 0.5 * np.arange(n_bus)
        wi = -0.75 + 0.25 * np.arange(n_bus)
        return float((wr * V.real + wi * V.imag).sum())

    def _loss_of_run(self, **grid_kwargs):
        grid = self._grid(**grid_kwargs)
        sa = self._analysis(grid)
        sa.compute(self._v0(grid), 30, 1e-12)
        return self._loss(np.array(sa.get_voltages()))

    def test_the_adjoint_is_the_gradient_a_finite_difference_measures(self):
        """lambda here is the gradient with respect to the per-unit injection of every
        bus, in every contingency at once. The reference is a central finite difference
        of the same loss with the whole analysis re-run on a perturbed grid."""
        grid = self._grid()
        sa = self._analysis(grid, keep_jacobian=True)
        sa.compute(self._v0(grid), 30, 1e-12)
        self.assertGreater(sa.adjoint_memory_bytes(), 0)

        V = np.array(sa.get_voltages())
        n_row, n_bus = V.shape
        dim_J = sa.dim_J()
        theta_col = np.asarray(sa.get_theta_col_of_bus())
        vm_col = np.asarray(sa.get_vm_col_of_bus())
        p_row = np.asarray(sa.get_p_row_of_bus())

        wr = 1. + 0.5 * np.arange(n_bus)
        wi = -0.75 + 0.25 * np.arange(n_bus)
        gV = wr + 1j * wi

        xbar = np.zeros((n_row, dim_J))
        for b in range(n_bus):
            v = V[:, b]
            if theta_col[b] >= 0:
                xbar[:, theta_col[b]] = (np.conj(v) * gV[b]).imag
            if vm_col[b] >= 0:
                with np.errstate(invalid="ignore", divide="ignore"):
                    xbar[:, vm_col[b]] = np.nan_to_num((np.conj(v / np.abs(v)) * gV[b]).real)

        lam = np.asarray(sa.solve_JT(np.ascontiguousarray(xbar)))
        self.assertEqual(lam.shape, (n_row, dim_J))
        self.assertTrue(all(sa.adjoint_row_ok()))

        # scale every load by (1 +- eps): the loss then moves by the sum over load
        # buses of the gradient times that bus' own injection
        sn_mva = grid.get_sn_mva()
        load_p = np.asarray(grid.get_load_target_p())
        load_bus = np.asarray(grid.get_loads().get_bus_id())
        eps = 1e-6
        fd = (self._loss_of_run(load_scale=1. + eps)
              - self._loss_of_run(load_scale=1. - eps)) / (2. * eps)

        # d loss / d scale = sum over rows and loads of -lambda[p_row[bus]] * p / sn_mva
        analytic = 0.
        for el, bus in enumerate(load_bus):
            if bus < 0 or p_row[bus] < 0:
                continue
            analytic += float((-lam[:, p_row[bus]] * load_p[el] / sn_mva).sum())
        self.assertAlmostEqual(analytic, fd, places=4,
                               msg=f"adjoint {analytic} vs finite difference {fd}")

    def test_the_gen_v_gradient_is_reachable_and_right(self):
        grid = self._grid()
        sa = self._analysis(grid, keep_jacobian=True)
        sa.compute(self._v0(grid), 30, 1e-12)

        V = np.array(sa.get_voltages())
        n_row, n_bus = V.shape
        dim_J = sa.dim_J()
        theta_col = np.asarray(sa.get_theta_col_of_bus())
        vm_col = np.asarray(sa.get_vm_col_of_bus())
        wr = 1. + 0.5 * np.arange(n_bus)
        wi = -0.75 + 0.25 * np.arange(n_bus)
        gV = wr + 1j * wi

        xbar = np.zeros((n_row, dim_J))
        for b in range(n_bus):
            v = V[:, b]
            if theta_col[b] >= 0:
                xbar[:, theta_col[b]] = (np.conj(v) * gV[b]).imag
            if vm_col[b] >= 0:
                xbar[:, vm_col[b]] = (np.conj(v / np.abs(v)) * gV[b]).real

        lam = np.asarray(sa.solve_JT(np.ascontiguousarray(xbar)))
        grad = np.asarray(sa.gen_v_indirect_grad(np.ascontiguousarray(lam)))
        target = np.asarray(sa.get_gen_v_target_bus())
        live = np.flatnonzero(target >= 0)
        self.assertGreater(live.size, 0)
        for g in live:                                  # + the direct half
            b = target[g]
            grad[:, g] += (np.conj(V[:, b] / np.abs(V[:, b])) * gV[b]).real

        # every generator's set-point moved together, so the finite difference is the
        # sum of the live columns
        eps = 1e-6
        fd = (self._loss_of_run(gen_v_delta=eps)
              - self._loss_of_run(gen_v_delta=-eps)) / (2. * eps)
        self.assertAlmostEqual(float(grad.sum()), fd, places=3,
                               msg=f"gen_v gradient {grad.sum()} vs finite difference {fd}")
        # a generator whose set-point never reaches the solve keeps a zero column
        for g in np.flatnonzero(target < 0):
            np.testing.assert_array_equal(grad[:, g], np.zeros(n_row))


if __name__ == "__main__":
    unittest.main()
