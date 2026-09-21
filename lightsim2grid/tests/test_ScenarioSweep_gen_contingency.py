# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
ScenarioSweep's third contingency axis: generator disconnection
(``set_contingency_gens``, see batch_algorithm/BaseBatchSweep.hpp).

The oracle throughout is a plain, one-off powerflow: for each row, a fresh copy of the
grid with that row's generators actually deactivated, solved by ``ac_pf``. A sweep row
must land on the same voltages -- including when the row's contingency takes the LAST
locally voltage-regulating generator off a bus, which turns that bus from PV to PQ and
is the whole point of the feature.

The second thing pinned here is the reason the feature exists at all: doing that must
NOT cost a symbolic re-factorization per row. ``test_single_symbolic_analysis`` asserts
the linear solver analyzes once for the whole sweep and refactorizes afterwards -- if a
row ever raises ``has_pv_changed()``, that count is what catches it.
"""

import copy
import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.scenarioSweep import ScenarioSweep, ScenarioSweepCPP


class TestScenarioSweepGenContingency(unittest.TestCase):
    def setUp(self):
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(),
                                    param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.Vinit = self.env.backend.V
        self.max_it = self.env.backend.max_it
        self.tol = self.env.backend.tol
        self.n_gen = self.env.n_gen
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())
        # case14: gens 2 and 3 share bus 5, so masking one of them leaves the bus PV and
        # masking both turns it PQ -- the two cases that must not be confused.
        self.bus_of_gen = [g.bus_id for g in self.grid.get_generators()]

        # real per-row injections: without them every row starts already converged
        # (Vinit IS the base solution and nothing varies), so the solver never builds a
        # Jacobian and a comparison proves nothing about this feature.
        data = self.env.chronics_handler.real_data.data
        self.nb_steps = 8
        self.gen_p = 1.0 * data.prod_p[:self.nb_steps]
        self.load_p = 1.0 * data.load_p[:self.nb_steps]
        self.load_q = 1.0 * data.load_q[:self.nb_steps]

    def tearDown(self):
        self.env.close()

    # ------------------------------------------------------------------ helpers
    def _reference_V(self, gens_off, lines_off=()):
        """One-off powerflow on a copy of the grid with those elements really removed."""
        grid = copy.deepcopy(self.grid)
        for gen_id in gens_off:
            grid.deactivate_gen(int(gen_id))
        for line_id in lines_off:
            grid.deactivate_powerline(int(line_id))
        V = grid.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
        return V

    def _sweep(self, gen_mask, line_mask=None):
        sweep = ScenarioSweepCPP(self.grid)
        sweep.set_contingency_gens(gen_mask)
        if line_mask is not None:
            sweep.set_contingency_lines(line_mask)
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        return sweep

    def _solved_buses(self):
        """Grid bus ids that are actually part of the solved system.

        Everything else -- an unused second bus of a substation, say -- is left at
        exactly complex 0 by the sweep but carries Vinit in a one-off ``ac_pf``, so
        comparing those columns would compare two different conventions, not results.
        """
        return np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _assert_row_matches(self, sweep, row, gens_off, lines_off=()):
        ref = self._reference_V(gens_off, lines_off)
        self.assertGreater(ref.shape[0], 0,
                           f"the reference powerflow itself diverged for row {row}")
        got = sweep.get_voltages()[row]
        buses = self._solved_buses()
        self.assertGreater(buses.size, 0, "no solved bus to compare")
        np.testing.assert_allclose(got[buses], ref[buses], rtol=1e-8, atol=1e-8,
                                   err_msg=f"row {row}, generators {list(gens_off)} off")

    # ------------------------------------------------------------------ tests
    def test_gen_n1_matches_one_off_powerflow(self):
        """One row per generator, each disconnecting it, against a real one-off solve."""
        gen_mask = np.zeros((self.n_gen, self.n_gen), dtype=bool)
        np.fill_diagonal(gen_mask, True)
        # gen 5 carries the whole slack weight on case14: taking it out is a separate
        # test (test_slack_gen_off), not this one.
        rows = [g for g in range(self.n_gen) if not self.grid.get_generators()[g].is_slack]
        gen_mask = gen_mask[rows]
        sweep = self._sweep(gen_mask)
        for row, gen_id in enumerate(rows):
            with self.subTest(gen=gen_id):
                self.assertTrue(sweep.converged_mask()[row], f"row {row} did not converge")
                self._assert_row_matches(sweep, row, [gen_id])

    def test_bus_stays_pv_while_one_gen_remains(self):
        """Gens 2 and 3 share bus 5: one off keeps it PV, both off turns it PQ."""
        shared = [g for g in range(self.n_gen)
                  if self.bus_of_gen.count(self.bus_of_gen[g]) > 1]
        self.assertGreaterEqual(len(shared), 2, "this test needs a bus with two generators")
        g_a, g_b = shared[0], shared[1]
        gen_mask = np.zeros((3, self.n_gen), dtype=bool)
        gen_mask[0, g_a] = True             # bus keeps a controller -> stays PV
        gen_mask[1, g_b] = True             # ditto
        gen_mask[2, [g_a, g_b]] = True      # bus loses both -> becomes PQ
        sweep = self._sweep(gen_mask)
        for row, gens_off in enumerate([[g_a], [g_b], [g_a, g_b]]):
            with self.subTest(row=row):
                self.assertTrue(sweep.converged_mask()[row])
                self._assert_row_matches(sweep, row, gens_off)

        # and the distinction is real: the bus' magnitude is pinned in rows 0/1 and free
        # in row 2, so row 2 must differ from the others at that bus
        bus = self.bus_of_gen[g_a]
        v = sweep.get_voltages()
        self.assertAlmostEqual(abs(v[0][bus]), abs(v[1][bus]), places=8,
                               msg="both rows keep a controller: |V| should hold at the setpoint")
        self.assertNotAlmostEqual(abs(v[0][bus]), abs(v[2][bus]), places=6,
                                  msg="losing every controller must free the bus' magnitude")

    def _with_injections(self, sweep):
        sweep.modify_gen_p(self.gen_p)
        sweep.modify_load_p(self.load_p)
        sweep.modify_load_q(self.load_q)
        return sweep

    def test_empty_mask_is_bit_identical(self):
        """An all-False mask must reproduce the plain, contingency-free sweep exactly.

        Run over varying injections on purpose, so every row really solves: an
        all-defaults sweep starts at its own answer and converges in zero iterations,
        which would compare two Jacobians neither of which was ever built.
        """
        with_mask = ScenarioSweepCPP(self.grid)
        self._with_injections(with_mask)
        with_mask.set_contingency_gens(np.zeros((self.nb_steps, self.n_gen), dtype=bool))
        with_mask.compute(1.0 * self.Vinit, self.max_it, self.tol)

        plain = ScenarioSweepCPP(self.grid)
        self._with_injections(plain)
        plain.set_contingency_lines(np.zeros((self.nb_steps, self.n_line), dtype=bool))
        plain.compute(1.0 * self.Vinit, self.max_it, self.tol)

        self.assertGreater(with_mask.get_linear_solver_stats().nb_factorize, 0,
                           "the rows must actually solve for this comparison to mean anything")
        np.testing.assert_array_equal(with_mask.get_voltages(), plain.get_voltages())

    def test_matches_one_off_powerflow_with_injections(self):
        """The real case: this row's own injections AND this row's own generator out."""
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack]
        gen_mask = np.zeros((self.nb_steps, self.n_gen), dtype=bool)
        for row in range(self.nb_steps):
            gen_mask[row, non_slack[row % len(non_slack)]] = True

        sweep = ScenarioSweepCPP(self.grid)
        self._with_injections(sweep)
        sweep.set_contingency_gens(gen_mask)
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)

        buses = self._solved_buses()
        for row in range(self.nb_steps):
            gen_id = non_slack[row % len(non_slack)]
            with self.subTest(row=row, gen=gen_id):
                self.assertTrue(sweep.converged_mask()[row])
                grid = copy.deepcopy(self.grid)
                # update_* take an (all-True mask, values) pair, values as float32
                grid.update_gens_p(np.ones(self.n_gen, dtype=bool),
                                   self.gen_p[row].astype(np.float32))
                grid.update_loads_p(np.ones(self.env.n_load, dtype=bool),
                                    self.load_p[row].astype(np.float32))
                grid.update_loads_q(np.ones(self.env.n_load, dtype=bool),
                                    self.load_q[row].astype(np.float32))
                grid.deactivate_gen(int(gen_id))
                ref = grid.ac_pf(1.0 * self.Vinit, self.max_it, self.tol)
                self.assertGreater(ref.shape[0], 0, "the reference itself diverged")
                np.testing.assert_allclose(sweep.get_voltages()[row][buses], ref[buses],
                                           rtol=1e-8, atol=1e-8)

    def test_combined_with_line_contingency(self):
        """A row disconnecting both a generator and a powerline."""
        gen_id = next(g for g in range(self.n_gen)
                      if not self.grid.get_generators()[g].is_slack)
        line_id = 0
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)
        line_mask = np.zeros((2, self.n_line), dtype=bool)
        gen_mask[1, gen_id] = True
        line_mask[1, line_id] = True
        sweep = self._sweep(gen_mask, line_mask)
        self._assert_row_matches(sweep, 0, [])
        self._assert_row_matches(sweep, 1, [gen_id], [line_id])

    def test_slack_gen_off(self):
        """Disconnecting a slack-participating generator re-weights the slack."""
        slack_gens = [g for g in range(self.n_gen) if self.grid.get_generators()[g].is_slack]
        self.assertTrue(slack_gens, "this test needs a generator carrying the slack")
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, slack_gens[0]] = True
        sweep = self._sweep(gen_mask)
        # the reference keeps the same reference bus (deactivate_gen does not move it),
        # which is exactly the contract: the slack bus set is a property of the batch
        self.assertTrue(sweep.converged_mask()[0])
        self._assert_row_matches(sweep, 0, [slack_gens[0]])

    def test_single_symbolic_analysis(self):
        """The point of the feature: N rows, ONE symbolic analysis."""
        nb_rows = 12
        rng = np.random.default_rng(0)
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack]
        gen_mask = np.zeros((nb_rows, self.n_gen), dtype=bool)
        for row in range(nb_rows):
            gen_mask[row, rng.choice(non_slack)] = True
        sweep = self._sweep(gen_mask)
        stats = sweep.get_linear_solver_stats()
        self.assertEqual(stats.nb_analyze, 1,
                         f"expected a single symbolic analysis for the whole sweep, "
                         f"got {stats.nb_analyze} -- a row is changing the sparsity pattern")
        self.assertGreater(stats.nb_refactorize, nb_rows,
                           "the rows should be running on refactorizations")

    def test_one_analysis_per_algorithm_when_threaded(self):
        """Multi-threaded: one analyze per WORKER, not one per row and not one overall.

        Each worker owns its own algorithm, so the count that stays flat is per algorithm.
        The sum alone cannot tell "every algorithm analyzed once" from "one algorithm
        analyzed four times", which is why the per-algorithm breakdown is checked instead.
        """
        nb_rows, nb_thread = 24, 4
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack]
        gen_mask = np.zeros((nb_rows, self.n_gen), dtype=bool)
        for row in range(nb_rows):
            gen_mask[row, non_slack[row % len(non_slack)]] = True

        sweep = ScenarioSweepCPP(self.grid)
        sweep.nb_thread = nb_thread
        sweep.set_contingency_gens(gen_mask)
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)

        per_algo = sweep.get_linear_solver_stats_per_algo()
        self.assertEqual(len(per_algo), 1 + nb_thread,
                         "expected the member algorithm plus one per worker")
        for i, st in enumerate(per_algo):
            with self.subTest(algo=i):
                self.assertLessEqual(st.nb_analyze, 1,
                                     f"algorithm {i} analyzed {st.nb_analyze} times; "
                                     f"each one should analyze at most once")
        self.assertEqual(sweep.get_linear_solver_stats().nb_analyze,
                         sum(st.nb_analyze for st in per_algo),
                         "the aggregate must be the sum of the per-algorithm counts")
        self.assertGreater(sum(st.nb_refactorize for st in per_algo), nb_rows,
                           "the rows should be running on refactorizations")

    def test_remote_controller_is_refused(self):
        """Remote voltage control is out of scope, and must say so rather than mislead."""
        grid = copy.deepcopy(self.grid)
        gen_id = next(g for g in range(self.n_gen)
                      if not grid.get_generators()[g].is_slack)
        own_bus = grid.get_generators()[gen_id].bus_id
        other_bus = next(b for b in self.bus_of_gen if b != own_bus)
        grid.set_gen_regulated_bus(gen_id, int(other_bus))
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, gen_id] = True
        sweep = ScenarioSweepCPP(grid)
        sweep.set_contingency_gens(gen_mask)
        with self.assertRaises(RuntimeError) as ctx:
            sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)
        self.assertIn("remote", str(ctx.exception).lower())

    def test_multithreaded_matches_single_threaded(self):
        """Each worker owns its algorithm, so the switchable set must reach all of them.

        The per-row helpers are const and touch no shared mutable state (which is why
        ``get_slack_weights_solver_without`` is const and does not refresh the container's
        weight cache) -- this is what pins that.
        """
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack]
        nb_rows = 16
        gen_mask = np.zeros((nb_rows, self.n_gen), dtype=bool)
        for row in range(nb_rows):
            gen_mask[row, non_slack[row % len(non_slack)]] = True

        single = self._sweep(gen_mask)
        multi = ScenarioSweepCPP(self.grid)
        multi.nb_thread = 4
        multi.set_contingency_gens(gen_mask)
        multi.compute(1.0 * self.Vinit, self.max_it, self.tol)

        np.testing.assert_array_equal(single.converged_mask(), multi.converged_mask())
        np.testing.assert_allclose(single.get_voltages(), multi.get_voltages(),
                                   rtol=1e-10, atol=1e-10)

    def test_with_handle_disconnected_grid(self):
        """The masked code path takes the generator contingencies too.

        ``handle_disconnected_grid`` runs its own per-row loop (_run_range_masked), so a
        generator contingency has to be applied there as well -- including on rows whose
        line contingency also strands part of the grid, where the bus mask and the PV
        pinning have to compose.
        """
        non_slack = [g for g in range(self.n_gen)
                     if not self.grid.get_generators()[g].is_slack]
        gen_id = non_slack[0]
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)
        line_mask = np.zeros((2, self.n_line), dtype=bool)
        gen_mask[0, gen_id] = True
        gen_mask[1, gen_id] = True
        line_mask[1, 0] = True

        sweep = ScenarioSweepCPP(self.grid)
        sweep.compute_limit_violations = True   # must be set first, see the C++ note
        sweep.handle_disconnected_grid = True
        sweep.set_contingency_gens(gen_mask)
        sweep.set_contingency_lines(line_mask)
        sweep.compute(1.0 * self.Vinit, self.max_it, self.tol)

        self.assertTrue(sweep.converged_mask()[0])
        self._assert_row_matches(sweep, 0, [gen_id])
        if sweep.converged_mask()[1]:
            self._assert_row_matches(sweep, 1, [gen_id], [0])

    def test_python_wrapper_validates_shape(self):
        sweep = ScenarioSweep(self.env)
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_gens(np.zeros((2, self.n_gen + 1), dtype=bool))
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_gens(np.zeros(self.n_gen, dtype=bool))
        sweep.set_contingency_gens(np.zeros((2, self.n_gen), dtype=bool))
        # the row count is now locked, like every other setter
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_lines(np.zeros((3, self.n_line), dtype=bool))


if __name__ == "__main__":
    unittest.main()


class TestGenContingencyWithRemoteControlAcrossTrafo(unittest.TestCase):
    """Interaction with remote voltage control, which lives in a different part of J.

    A generator regulating a bus on the far side of a transformer is a ``VoltageControl``
    controller, not a PV bus: it owns a Q column and a setpoint row of its own, and
    disconnecting a transformer on the path can strand it. That case is handled and tested
    on the C++ side (``src/tests/test_batch_voltage_control.cpp``, "a contingency stranding
    a lone controller's own bus falls back to plain PQ"). What is new here is the
    combination: the generator-contingency axis has to compose with it, including on a row
    that disconnects a transformer AND a generator at once.

    Uses the same fixture as ``test_voltage_control_batch.py`` -- pandapower case14 with
    generator 3 (on bus 7) regulating bus 9 -- so the two files describe the same grid.
    """

    GEN_REMOTE = 3     # on bus 7 ...
    REG_BUS = 9        # ... regulating bus 9, remotely
    TRAFO_ON_PATH = 4  # buses 6-8: on the way, and the row still solves without it
    TRAFO_BEHIND = 3   # buses 6-7: the controller's OWN bus is behind this one
    MAX_IT, TOL = 30, 1e-11

    @classmethod
    def setUpClass(cls):
        try:
            import pandapower.networks as pn
        except ImportError:
            raise unittest.SkipTest("pandapower is needed for case14")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            cls._net = pn.case14()

    def _grid(self, trafo_off=None, gens_off=()):
        # imported here, not stored on the class: a plain function assigned to a class
        # attribute becomes a bound method, and `self` would be passed as its first arg
        from lightsim2grid.network import init_from_pandapower
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(self._net)
        grid.set_gen_regulated_bus(self.GEN_REMOTE, self.REG_BUS)
        if trafo_off is not None:
            grid.deactivate_trafo(int(trafo_off))
        for g in gens_off:
            grid.deactivate_gen(int(g))
        grid.tell_solver_need_reset()
        return grid

    def setUp(self):
        self.grid = self._grid()
        self.Vinit = np.ones(self.grid.total_bus(), dtype=complex)
        self.n_gen = len(self.grid.get_generators())
        self.n_trafo = len(self.grid.get_trafos())
        self.buses = np.asarray(self.grid.id_ac_solver_to_me(), dtype=int)

    def _reference(self, trafo_off, gens_off):
        return self._grid(trafo_off, gens_off).ac_pf(1.0 * self.Vinit, self.MAX_IT, self.TOL)

    def _sweep(self, trafo_mask, gen_mask):
        sweep = ScenarioSweepCPP(self._grid())
        sweep.set_contingency_trafos(trafo_mask)
        sweep.set_contingency_gens(gen_mask)
        sweep.compute(1.0 * self.Vinit, self.MAX_IT, self.TOL)
        return sweep

    def test_trafo_on_the_control_path_still_matches(self):
        """Disconnect a transformer between the controller and the bus it regulates."""
        trafo_mask = np.zeros((2, self.n_trafo), dtype=bool)
        trafo_mask[1, self.TRAFO_ON_PATH] = True
        gen_mask = np.zeros((2, self.n_gen), dtype=bool)   # axis on, but empty
        sweep = self._sweep(trafo_mask, gen_mask)

        for row, trafo in enumerate([None, self.TRAFO_ON_PATH]):
            with self.subTest(row=row, trafo=trafo):
                ref = self._reference(trafo, [])
                self.assertGreater(ref.shape[0], 0, "the reference itself diverged")
                self.assertTrue(sweep.converged_mask()[row])
                np.testing.assert_allclose(sweep.get_voltages()[row][self.buses],
                                           ref[self.buses], rtol=1e-8, atol=1e-8)

    def test_trafo_and_generator_contingency_on_the_same_row(self):
        """Both at once: a transformer out AND a bus turned PV -> PQ by a lost generator."""
        candidates = [g for g in range(self.n_gen)
                      if g != self.GEN_REMOTE
                      and self._reference(self.TRAFO_ON_PATH, [g]).shape[0] > 0]
        self.assertTrue(candidates, "no generator gives a converging reference")
        gen_id = candidates[0]

        trafo_mask = np.zeros((1, self.n_trafo), dtype=bool)
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        trafo_mask[0, self.TRAFO_ON_PATH] = True
        gen_mask[0, gen_id] = True
        sweep = self._sweep(trafo_mask, gen_mask)

        self.assertTrue(sweep.converged_mask()[0])
        ref = self._reference(self.TRAFO_ON_PATH, [gen_id])
        np.testing.assert_allclose(sweep.get_voltages()[0][self.buses], ref[self.buses],
                                   rtol=1e-8, atol=1e-8)

    def test_stranding_the_controller_agrees_with_the_one_off_solve(self):
        """The controller's own bus goes behind a disconnected transformer.

        Whatever the batch does here it must do what a one-off solve of the same case
        does -- the point is that adding the generator-contingency axis does not make the
        two disagree.
        """
        trafo_mask = np.zeros((1, self.n_trafo), dtype=bool)
        trafo_mask[0, self.TRAFO_BEHIND] = True
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        sweep = self._sweep(trafo_mask, gen_mask)

        ref = self._reference(self.TRAFO_BEHIND, [])
        ref_converged = ref.shape[0] > 0
        self.assertEqual(bool(sweep.converged_mask()[0]), ref_converged,
                         "batch and one-off solve disagree on whether this case solves")
        if ref_converged:
            np.testing.assert_allclose(sweep.get_voltages()[0][self.buses], ref[self.buses],
                                       rtol=1e-8, atol=1e-8)

    def test_masking_the_remote_controller_itself_is_still_refused(self):
        """The scope guard must not be weakened by the transformer cases above."""
        gen_mask = np.zeros((1, self.n_gen), dtype=bool)
        gen_mask[0, self.GEN_REMOTE] = True
        sweep = ScenarioSweepCPP(self._grid())
        sweep.set_contingency_gens(gen_mask)
        with self.assertRaises(RuntimeError) as ctx:
            sweep.compute(1.0 * self.Vinit, self.MAX_IT, self.TOL)
        self.assertIn("remote", str(ctx.exception).lower())
