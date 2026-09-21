# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Continuation powerflow (ContinuationSweepCPP + the ContinuationPowerFlow wrapper).

The oracles used here are deliberately independent of the continuation itself:

- every traced point must equal a plain ``LSGrid.ac_pf`` run at the same injections
  (``test_every_point_matches_a_direct_solve`` -- this is the correctness test, everything
  else is a property of it);
- the nose lambda must equal a brute-force bisection on the loading factor, i.e. the
  largest factor for which an ordinary powerflow still converges;
- the whole curve must cost exactly ONE symbolic factorization, which is the reason a
  continuation belongs in the batch layer at all.

Two traps have their own tests because getting them wrong produces a plausible-looking
curve rather than an error: the sign of the load direction (loads are stamped into Sbus
with a MINUS, so increasing them must LOWER the voltages), and a zero direction (which
would otherwise walk lambda along a curve on which nothing moves and report success).
"""

import unittest
import warnings

import numpy as np

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    import pandapower.networks as pn

from lightsim2grid import AlgorithmType
from lightsim2grid.network import init_from_pandapower
from lightsim2grid.continuationPowerflow import ContinuationPowerFlow, run_cpf


def _case14():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(pn.case14())


class TestContinuationPowerFlow(unittest.TestCase):
    def setUp(self):
        self.grid = _case14()
        self.nb_bus = self.grid.total_bus()
        self.base_load_p = np.array(self.grid.get_load_target_p(), dtype=float)
        self.base_load_q = np.array([el.target_q_mvar for el in self.grid.get_loads()], dtype=float)
        self.base_gen_p = np.array(self.grid.get_gen_target_p(), dtype=float)
        self.gen_is_slack = np.array([el.is_slack for el in self.grid.get_generators()], dtype=bool)
        self.nb_load = self.base_load_p.shape[0]
        self.nb_gen = self.base_gen_p.shape[0]
        self.v_init = np.full(self.nb_bus, 1.04, dtype=complex)

    # ---- helpers -----------------------------------------------------------
    def _solve_at(self, lam, k, alpha=None, beta=None, scale_q=True, max_iter=100, tol=1e-10):
        """A plain powerflow at the injections the continuation should have at `lam`."""
        if alpha is None:
            alpha = np.ones(self.nb_load)
        if beta is None:
            beta = np.ones(self.nb_gen)
        grid = _case14()
        for l_id, val in enumerate(self.base_load_p * (1.0 + alpha * lam * (k - 1.0))):
            grid.change_p_load(l_id, val)
        if scale_q:
            for l_id, val in enumerate(self.base_load_q * (1.0 + alpha * lam * (k - 1.0))):
                grid.change_q_load(l_id, val)
        for g_id, val in enumerate(self.base_gen_p * (1.0 + beta * lam * (k - 1.0))):
            grid.change_p_gen(g_id, val)
        return grid.ac_pf(self.v_init, max_iter, tol)

    # ---- correctness -------------------------------------------------------
    def test_every_point_matches_a_direct_solve(self):
        """Each traced (lam, V) is a genuine solution of the powerflow at that lambda."""
        k = 4.0
        res = ContinuationPowerFlow(_case14()).run(loading_factor=k)
        self.assertTrue(res.success, res.msg)
        self.assertGreater(res.lam.size, 20)
        # Every 20th point. The tolerance is keyed on tangent_lam, which IS the
        # conditioning of the point: near the nose the Jacobian is nearly singular and
        # the reference solve is just as ill-conditioned as the traced one, so demanding
        # 1e-8 there would be testing the arithmetic, not the algorithm.
        for i in list(range(1, res.lam.size - 1, 20)):
            V = self._solve_at(res.lam[i], k)
            self.assertGreater(V.shape[0], 0, f"the reference solve diverged at lam={res.lam[i]}")
            atol = 1e-8 if res.tangent_lam[i] > 1e-3 else 1e-4
            np.testing.assert_allclose(V, res.V[i], atol=atol,
                                       err_msg=f"traced point {i} (lam={res.lam[i]}, "
                                               f"tangent_lam={res.tangent_lam[i]:.2e}) is not a solution")

    def test_nose_matches_a_bisection_on_the_loading_factor(self):
        """
        The nose is the largest lambda for which a powerflow still converges. Found here
        without any continuation at all, by bisection, and compared with what the CPF
        reports.
        """
        k = 4.0
        res = ContinuationPowerFlow(_case14()).run(loading_factor=k)
        self.assertTrue(res.success, res.msg)

        def solvable(lam):
            return self._solve_at(lam, k, max_iter=100, tol=1e-8).shape[0] > 0

        lo, hi = 0.0, 3.0
        self.assertTrue(solvable(lo))
        self.assertFalse(solvable(hi))
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if solvable(mid):
                lo = mid
            else:
                hi = mid
        self.assertAlmostEqual(res.lam_max, lo, delta=1e-4 * max(1.0, lo))

    def test_one_analyze_for_the_whole_curve(self):
        """
        The entire premise of putting a continuation in the batch layer: one symbolic
        factorization, however many points the curve has.
        """
        cpf = ContinuationPowerFlow(_case14())
        res = cpf.run(loading_factor=4.0)
        self.assertGreater(res.lam.size, 50)
        stats = cpf.cpp.get_linear_solver_stats()
        self.assertEqual(stats.nb_analyze, 1)
        self.assertGreater(stats.nb_refactorize, res.lam.size)

    def test_loading_lowers_the_voltages(self):
        """
        The load direction's sign. Loads are stamped into Sbus with a MINUS, so a
        direction that increases them must DECREASE the voltages -- get this backwards
        and the continuation happily traces the grid unloading instead.
        """
        res = ContinuationPowerFlow(_case14()).run(loading_factor=4.0)
        alive = np.flatnonzero(res.Vm[0] > 0.0)
        worst = alive[np.argmax(res.Vm[0, alive] - res.Vm[-1, alive])]
        self.assertLess(res.Vm[-1, worst], res.Vm[0, worst])
        # and it gets there monotonically, not by wandering
        self.assertTrue(np.all(np.diff(res.Vm[:, worst]) <= 1e-9))

    def test_tangent_lam_falls_to_zero_at_the_nose(self):
        res = ContinuationPowerFlow(_case14()).run(loading_factor=4.0)
        self.assertTrue(res.success, res.msg)
        # strictly positive throughout (this parameterisation cannot make it change
        # sign), and collapsing towards 0 at the end
        self.assertTrue(np.all(res.tangent_lam[:-1] > 0.0))
        self.assertLess(res.tangent_lam[-2], res.tangent_lam[0])
        self.assertLess(res.tangent_lam[-2], 1e-4)

    # ---- steering ----------------------------------------------------------
    def test_all_ones_steering_is_the_default(self):
        a = ContinuationPowerFlow(_case14()).run(loading_factor=3.0)
        b = ContinuationPowerFlow(_case14()).run(loading_factor=3.0,
                                                 load_steering=np.ones(self.nb_load))
        np.testing.assert_array_equal(a.lam, b.lam)
        np.testing.assert_array_equal(a.V, b.V)

    def test_steering_one_load_moves_only_that_bus(self):
        load_id = 3
        alpha = np.zeros(self.nb_load)
        alpha[load_id] = 1.0
        cpf = ContinuationPowerFlow(_case14())
        cpf.run(loading_factor=3.0, load_steering=alpha, gen_steering=0.0)
        direction = np.array(cpf.cpp.get_direction_solver())
        moved = np.flatnonzero(np.abs(direction) > 1e-12)
        self.assertEqual(moved.size, 1)
        # case14 has no disconnected bus, so solver and grid ids coincide here
        self.assertEqual(int(moved[0]), self.grid.get_loads()[load_id].bus_id)

    def test_steered_curve_matches_a_direct_solve(self):
        """The steering reaches the solver: the traced points are solutions of the
        STEERED problem, not of the uniform one."""
        k = 3.0
        alpha = np.zeros(self.nb_load)
        alpha[[2, 5]] = [1.0, 0.4]
        beta = np.zeros(self.nb_gen)
        res = ContinuationPowerFlow(_case14()).run(loading_factor=k, load_steering=alpha,
                                                   gen_steering=0.0)
        for i in [1, res.lam.size // 2]:
            V = self._solve_at(res.lam[i], k, alpha=alpha, beta=beta)
            self.assertGreater(V.shape[0], 0)
            np.testing.assert_allclose(V, res.V[i], atol=1e-8)

    def test_zero_coefficient_load_does_not_move(self):
        """A load steered with 0 keeps its base value all along the curve."""
        alpha = np.ones(self.nb_load)
        alpha[4] = 0.0
        cpf = ContinuationPowerFlow(_case14())
        target = cpf._build_target(3.0, alpha, None, True, None)
        self.assertAlmostEqual(target["load_p"][4], self.base_load_p[4])
        self.assertAlmostEqual(target["load_q"][4], self.base_load_q[4])

    def test_gen_steering_default_scales_every_generator_slack_included(self):
        """
        Slack machines are NOT a special case, unlike in MATPOWER's target-case rule --
        see the class docstring. The two tests below pin the two facts that justify it.
        """
        cpf = ContinuationPowerFlow(_case14())
        target = cpf._build_target(3.0, None, None, True, None)
        np.testing.assert_allclose(target["gen_p"], 3.0 * self.base_gen_p)

    def test_single_slack_machine_setpoint_is_inert(self):
        """
        Why excluding the slack would be a no-op rather than a correction, with ONE slack:
        that bus has no active-power equation, so the machine's target_p never reaches the
        mismatch. Scaling it leaves the solution bit-identical -- and the machine's actual
        output unchanged, because that output is what balances the grid, not what was asked
        of it.
        """
        slack_id = int(np.flatnonzero(self.gen_is_slack)[0])
        self.assertGreater(self.base_gen_p[slack_id], 0.0)  # or the test proves nothing

        grid_a = _case14()
        V_a = grid_a.ac_pf(self.v_init, 30, 1e-11)
        grid_b = _case14()
        grid_b.change_p_gen(slack_id, 2.0 * self.base_gen_p[slack_id])
        V_b = grid_b.ac_pf(self.v_init, 30, 1e-11)

        np.testing.assert_array_equal(V_a, V_b)
        self.assertAlmostEqual(grid_a.get_generators()[slack_id].res_p_mw,
                               grid_b.get_generators()[slack_id].res_p_mw, places=9)

    def test_distributed_slack_participant_setpoint_is_not_inert(self):
        """
        ... and why the exclusion cannot simply be kept anyway: with a DISTRIBUTED slack a
        participant's target_p does enter its bus' equation, so excluding every slack
        machine would silently freeze real generation.
        """
        slack_id = int(np.flatnonzero(self.gen_is_slack)[0])

        def solve(scale):
            grid = _case14()
            grid.add_gen_slackbus(1, 0.5)  # a second participant -> distributed slack
            if scale:
                grid.change_p_gen(slack_id, 2.0 * self.base_gen_p[slack_id])
            return grid.ac_pf(self.v_init, 30, 1e-11)

        V_a, V_b = solve(False), solve(True)
        self.assertGreater(np.max(np.abs(V_a - V_b)), 1e-3)

    def test_direction_the_slack_absorbs_stops_instead_of_being_traced(self):
        """
        Steering ONLY the single slack machine is a legal thing to ask for -- it moves a
        real input, just not one the powerflow reads -- so it is not refused. What it must
        not do is look like a result: before this was handled, lambda marched to 20.8 over
        a thousand points on which |V| changed by 3e-15, and the run reported success.
        """
        beta = self.gen_is_slack.astype(float)
        res = ContinuationPowerFlow(_case14()).run(loading_factor=2.0,
                                                   load_steering=np.zeros(self.nb_load),
                                                   gen_steering=beta)
        self.assertFalse(res.success)
        self.assertEqual(res.lam.size, 1)      # the base case, and nothing past it
        self.assertEqual(res.lam_max, 0.0)
        self.assertIn("absorbed by the slack", res.msg)

    def test_gen_steering_zero_leaves_generation_fixed(self):
        cpf = ContinuationPowerFlow(_case14())
        target = cpf._build_target(3.0, None, 0.0, True, None)
        np.testing.assert_allclose(target["gen_p"], self.base_gen_p)

    def test_scale_q_false_keeps_reactive_load(self):
        cpf = ContinuationPowerFlow(_case14())
        target = cpf._build_target(3.0, None, None, False, None)
        self.assertNotIn("load_q", target)

    def test_slack_supplied_margin_is_smaller_than_the_shared_one(self):
        """
        Sanity on the modelling choice: holding generation fixed makes the slack supply
        the whole increase, which collapses the grid earlier than sharing it out.
        """
        shared = ContinuationPowerFlow(_case14()).run(loading_factor=4.0)
        slack_only = ContinuationPowerFlow(_case14()).run(loading_factor=4.0, gen_steering=0.0)
        self.assertLess(slack_only.lam_max, shared.lam_max)

    # ---- explicit direction ------------------------------------------------
    def test_explicit_direction(self):
        delta = np.zeros(self.nb_load)
        delta[1] = 50.0  # +50 MW on load 1, nothing else moves
        cpf = ContinuationPowerFlow(_case14())
        cpf.run(direction={"load_p": delta}, stop_at_lam=1.0)
        direction = np.array(cpf.cpp.get_direction_solver())
        moved = np.flatnonzero(np.abs(direction) > 1e-12)
        self.assertEqual(moved.size, 1)
        self.assertEqual(int(moved[0]), self.grid.get_loads()[1].bus_id)
        # a load is a NEGATIVE injection: +50 MW of load is -50 MW at the bus
        sn_mva = self.grid.get_sn_mva()
        self.assertAlmostEqual(direction[moved[0]].real, -50.0 / sn_mva)

    def test_direction_and_steering_are_mutually_exclusive(self):
        with self.assertRaises(ValueError):
            ContinuationPowerFlow(_case14()).run(direction={"load_p": np.zeros(self.nb_load)},
                                                 load_steering=np.ones(self.nb_load))

    # ---- stopping ----------------------------------------------------------
    def test_stop_at_lam_lands_exactly_on_it(self):
        res = ContinuationPowerFlow(_case14()).run(loading_factor=2.0, stop_at_lam=0.5)
        self.assertTrue(res.success, res.msg)
        self.assertAlmostEqual(res.lam[-1], 0.5, places=12)
        # and that point is a real solution
        V = self._solve_at(0.5, 2.0)
        np.testing.assert_allclose(V, res.V[-1], atol=1e-8)

    def test_max_steps_stops_the_curve(self):
        res = ContinuationPowerFlow(_case14()).run(loading_factor=4.0, max_steps=5)
        self.assertFalse(res.success)
        self.assertEqual(res.lam.size, 6)  # the base case + 5 steps
        self.assertIn("maximum number of steps", res.msg)

    def test_adapt_step_reaches_the_same_nose_in_fewer_points(self):
        fixed = ContinuationPowerFlow(_case14()).run(loading_factor=4.0, adapt_step=False)
        adapt = ContinuationPowerFlow(_case14()).run(loading_factor=4.0, adapt_step=True)
        self.assertTrue(adapt.success, adapt.msg)
        self.assertAlmostEqual(adapt.lam_max, fixed.lam_max, delta=1e-3)
        self.assertLess(adapt.lam.size, fixed.lam.size)

    def test_exact_tangent_reaches_the_same_nose(self):
        loose = ContinuationPowerFlow(_case14()).run(loading_factor=4.0)
        exact = ContinuationPowerFlow(_case14()).run(loading_factor=4.0, exact_tangent=True)
        self.assertTrue(exact.success, exact.msg)
        self.assertAlmostEqual(exact.lam_max, loose.lam_max, delta=1e-3)

    # ---- refusals ----------------------------------------------------------
    def test_zero_direction_is_refused(self):
        """Not a degenerate curve, a meaningless run -- see the module docstring."""
        with self.assertRaises(RuntimeError) as ctx:
            ContinuationPowerFlow(_case14()).run(loading_factor=3.0,
                                                 load_steering=np.zeros(self.nb_load),
                                                 gen_steering=0.0)
        self.assertIn("direction is zero", str(ctx.exception))

    def test_non_newton_algorithm_is_refused(self):
        cpf = ContinuationPowerFlow(_case14(), algorithm=AlgorithmType.GaussSeidel)
        with self.assertRaises(RuntimeError) as ctx:
            cpf.run(loading_factor=2.0)
        self.assertIn("Newton-Raphson", str(ctx.exception))

    def test_dc_algorithm_is_refused(self):
        cpf = ContinuationPowerFlow(_case14(), algorithm=AlgorithmType.DC_SparseLU)
        with self.assertRaises(RuntimeError) as ctx:
            cpf.run(loading_factor=2.0)
        self.assertIn("AC", str(ctx.exception))

    def test_loading_factor_must_exceed_one(self):
        for bad in (1.0, 0.5, -1.0, np.nan):
            with self.assertRaises(ValueError):
                ContinuationPowerFlow(_case14()).run(loading_factor=bad)

    def test_steering_out_of_range_is_refused(self):
        for bad in (np.full(self.nb_load, 1.5), np.full(self.nb_load, -0.1)):
            with self.assertRaises(ValueError):
                ContinuationPowerFlow(_case14()).run(loading_factor=2.0, load_steering=bad)

    def test_steering_wrong_size_is_refused(self):
        with self.assertRaises(ValueError):
            ContinuationPowerFlow(_case14()).run(loading_factor=2.0,
                                                 load_steering=np.ones(self.nb_load + 1))

    def test_multithreading_is_refused(self):
        """The points are chained, so there is nothing to split."""
        cpf = ContinuationPowerFlow(_case14())
        with self.assertRaises(RuntimeError):
            cpf.cpp.set_nb_thread(4)

    # ---- API ---------------------------------------------------------------
    def test_run_cpf_helper(self):
        res = run_cpf(_case14(), loading_factor=2.0, stop_at_lam=0.25)
        self.assertTrue(res.success, res.msg)
        self.assertAlmostEqual(res.lam[-1], 0.25, places=12)
        self.assertEqual(res.V.shape, (res.lam.size, self.nb_bus))
        self.assertEqual(res.Vm.shape, res.V.shape)
        self.assertEqual(res.Va_deg.shape, res.V.shape)

    def test_results_survive_a_second_run(self):
        """The result object holds copies, not views into the C++ buffers."""
        cpf = ContinuationPowerFlow(_case14())
        first = cpf.run(loading_factor=2.0, stop_at_lam=0.3)
        lam_before = first.lam.copy()
        V_before = first.V.copy()
        cpf.run(loading_factor=4.0)
        np.testing.assert_array_equal(first.lam, lam_before)
        np.testing.assert_array_equal(first.V, V_before)


class TestContinuationPowerFlowMultiSlack(unittest.TestCase):
    """
    The tangent's right-hand side is built from the NRLedger, so it must be right for a
    distributed slack too -- where the augmented Jacobian carries an extra column (the
    absorbed slack power) that does not depend on lambda.
    """

    def setUp(self):
        self.grid = _case14()
        self.v_init = np.full(self.grid.total_bus(), 1.04, dtype=complex)

    def _distributed(self):
        grid = _case14()
        # every generator participates in the slack
        nb_gen = len(grid.get_generators())
        grid.update_slack_weights(np.ones(nb_gen, dtype=bool))
        return grid

    def test_distributed_slack_curve_matches_a_direct_solve(self):
        k = 3.0
        cpf = ContinuationPowerFlow(self._distributed())
        res = cpf.run(loading_factor=k)
        self.assertTrue(res.success, res.msg)
        self.assertEqual(cpf.cpp.get_linear_solver_stats().nb_analyze, 1)

        base_load_p = np.array(self.grid.get_load_target_p(), dtype=float)
        base_load_q = np.array([el.target_q_mvar for el in self.grid.get_loads()], dtype=float)
        base_gen_p = np.array(self.grid.get_gen_target_p(), dtype=float)
        # gen_steering's default is all ones -- slack participants included, which is
        # exactly what matters here: with a distributed slack their setpoints are real
        # inputs, so freezing them would trace a different curve (see
        # test_distributed_slack_participant_setpoint_is_not_inert).
        beta = np.ones(base_gen_p.shape[0])

        for i in [1, res.lam.size // 2]:
            lam = res.lam[i]
            grid = self._distributed()
            for l_id, val in enumerate(base_load_p * (1.0 + lam * (k - 1.0))):
                grid.change_p_load(l_id, val)
            for l_id, val in enumerate(base_load_q * (1.0 + lam * (k - 1.0))):
                grid.change_q_load(l_id, val)
            for g_id, val in enumerate(base_gen_p * (1.0 + beta * lam * (k - 1.0))):
                grid.change_p_gen(g_id, val)
            V = grid.ac_pf(self.v_init, 100, 1e-10)
            self.assertGreater(V.shape[0], 0, f"reference solve diverged at lam={lam}")
            np.testing.assert_allclose(V, res.V[i], atol=1e-8)


if __name__ == "__main__":
    unittest.main()
