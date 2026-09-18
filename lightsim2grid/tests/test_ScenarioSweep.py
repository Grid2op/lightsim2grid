# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
ScenarioSweepCPP is the 4th instantiation of the same C++ template as TimeSeriesCPP /
InjectionSweepCPP / ContingencyAnalysisCPP (see batch_algorithm/BaseBatchSweep.hpp): it
varies BOTH the injection and a contingency per row, independently, row-aligned.

The tests below use the two simpler algorithms as correctness oracles: a row that only
sets injections (contingency mask all-False) must match the equivalent InjectionSweep
row bit for bit; a row that only sets a contingency (injections left at the grid's own
state) must match the equivalent ContingencyAnalysis row. They also pin the API
separation from ContingencyAnalysis (set_contingency_lines/trafos vs add_n1/add_nk are
deliberately two different APIs) and the row-count locking / fail-fast validation.
"""

import unittest
import warnings

import numpy as np
import grid2op
from grid2op.Parameters import Parameters

from lightsim2grid import LightSimBackend
from lightsim2grid.scenarioSweep import ScenarioSweep, ScenarioSweepCPP
from lightsim2grid.injectionSweep import InjectionSweepCPP
from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP
from test_ContingencyAnalysis_limit_violations import _set_tight_limits


class TestScenarioSweepCPP(unittest.TestCase):
    def setUp(self):
        env_name = "l2rpn_case14_sandbox"
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(env_name, backend=LightSimBackend(), param=param, test=True)
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.Vinit = self.env.backend.V
        self.max_it = self.env.backend.max_it
        self.tol = self.env.backend.tol

        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())

        data = self.env.chronics_handler.real_data.data
        self.nb_steps = 8
        self.gen_p = 1.0 * data.prod_p[:self.nb_steps]
        self.load_p = 1.0 * data.load_p[:self.nb_steps]
        self.load_q = 1.0 * data.load_q[:self.nb_steps]
        self.sgen_p = np.zeros((self.nb_steps, 0))

        self.no_line_mask = np.zeros((self.nb_steps, self.n_line), dtype=bool)
        self.no_trafo_mask = np.zeros((self.nb_steps, self.n_trafo), dtype=bool)

    def tearDown(self):
        self.env.close()

    def _compute_injection_sweep(self):
        computer = InjectionSweepCPP(self.grid)
        status = computer.compute_Vs(self.gen_p, self.sgen_p, self.load_p, self.load_q,
                                     self.Vinit, self.max_it, self.tol)
        assert status == 1, f"InjectionSweep diverged after {computer.nb_solved()} step(s)"
        return 1.0 * computer.get_voltages()

    def _compute_contingency_analysis(self, line_id):
        computer = ContingencyAnalysisCPP(self.grid)
        computer.add_n1(line_id)
        computer.compute(self.Vinit, self.max_it, self.tol)
        return 1.0 * computer.get_voltages()

    def test_injection_only_row_matches_injection_sweep(self):
        """a row with no contingency (mask all-False) must match InjectionSweep exactly"""
        Vs_inj = self._compute_injection_sweep()

        computer = ScenarioSweepCPP(self.grid)
        computer.modify_gen_p(self.gen_p)
        computer.modify_sgen_p(self.sgen_p)
        computer.modify_load_p(self.load_p)
        computer.modify_load_q(self.load_q)
        computer.set_contingency_lines(self.no_line_mask)
        computer.set_contingency_trafos(self.no_trafo_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1
        Vs = 1.0 * computer.get_voltages()

        assert np.array_equal(Vs, Vs_inj), (
            f"max deviation {np.max(np.abs(Vs - Vs_inj))} vs InjectionSweep")

    def test_contingency_only_row_matches_contingency_analysis(self):
        """a row with injections left at grid defaults must match ContingencyAnalysis"""
        line_id = 3
        Vs_ca = self._compute_contingency_analysis(line_id)

        computer = ScenarioSweepCPP(self.grid)
        line_mask = np.zeros((1, self.n_line), dtype=bool)
        line_mask[0, line_id] = True
        computer.set_contingency_lines(line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1
        Vs = 1.0 * computer.get_voltages()

        assert np.max(np.abs(Vs - Vs_ca)) <= 1e-8, (
            f"max deviation {np.max(np.abs(Vs - Vs_ca))} vs ContingencyAnalysis")

    def test_unset_axis_defaults_to_grid_state(self):
        """modify_gen_p/load_q never called: they must default to the grid's own current
        target, broadcast across every row -- same result as passing that same constant
        value explicitly. (Using self.gen_p/self.load_q -- the time-varying chronics --
        here would NOT be a fair comparison: the default is a constant broadcast of the
        grid's *current* state, not a replay of the chronics.)"""
        gen_p_const = np.tile(self.grid.get_gen_target_p(), (self.nb_steps, 1))
        load_q_const = np.tile([l.target_q_mvar for l in self.grid.get_loads()], (self.nb_steps, 1))

        computer_explicit = ScenarioSweepCPP(self.grid)
        computer_explicit.modify_gen_p(gen_p_const)
        computer_explicit.modify_load_p(self.load_p)
        computer_explicit.modify_load_q(load_q_const)
        computer_explicit.compute(self.Vinit, self.max_it, self.tol)
        Vs_explicit = 1.0 * computer_explicit.get_voltages()

        computer_default = ScenarioSweepCPP(self.grid)
        computer_default.modify_load_p(self.load_p)  # only lock the row count via load_p
        computer_default.compute(self.Vinit, self.max_it, self.tol)
        Vs_default = 1.0 * computer_default.get_voltages()

        assert np.array_equal(Vs_explicit, Vs_default)

    def test_row_count_lock_fails_fast(self):
        computer = ScenarioSweepCPP(self.grid)
        computer.modify_load_p(self.load_p)  # locks nb_steps
        with self.assertRaises(RuntimeError) as ctx:
            computer.modify_gen_p(self.gen_p[:3])
        assert "ScenarioSweep" in str(ctx.exception)

    def test_compute_raises_if_nothing_set(self):
        computer = ScenarioSweepCPP(self.grid)
        with self.assertRaises(RuntimeError):
            computer.compute(self.Vinit, self.max_it, self.tol)

    def test_api_separation_from_contingency_analysis(self):
        """set_contingency_lines/trafos and add_n1/add_nk are deliberately two different
        APIs -- this guard fails if that ever regresses"""
        sweep = ScenarioSweepCPP(self.grid)
        ca = ContingencyAnalysisCPP(self.grid)

        assert hasattr(sweep, "set_contingency_lines")
        assert hasattr(sweep, "set_contingency_trafos")
        assert not hasattr(sweep, "add_n1")
        assert not hasattr(sweep, "add_nk")
        assert not hasattr(sweep, "my_defaults")

        assert hasattr(ca, "add_n1")
        assert hasattr(ca, "add_nk")
        assert not hasattr(ca, "set_contingency_lines")
        assert not hasattr(ca, "set_contingency_trafos")

        # legacy bundled compute_Vs: ScenarioSweep never had it
        assert not hasattr(sweep, "compute_Vs")
        assert not hasattr(sweep, "get_sbuses")

    def test_no_converged_or_converged_n_on_scenario_sweep(self):
        """handle_disconnected_grid / compute_limit_violations / get_violations /
        get_violations_n are now shared with ContingencyAnalysis (same names); but
        converged() / converged_n() deliberately stay ContingencyAnalysis-only -- a
        non-converged row's get_violations() entry already carries a GRID-type
        NOT_SIMULATED/DIVERGENCE sentinel, so a separate convergence flag on
        ScenarioSweep would be redundant. This is its own test, separate from the
        general API-separation sweep above, so a future regression here is easy to
        diagnose on its own."""
        sweep = ScenarioSweepCPP(self.grid)
        ca = ContingencyAnalysisCPP(self.grid)

        assert hasattr(sweep, "handle_disconnected_grid")
        assert hasattr(sweep, "compute_limit_violations")
        assert hasattr(sweep, "violation_threshold")
        assert hasattr(sweep, "get_violations")
        assert hasattr(sweep, "get_violations_n")
        assert not hasattr(sweep, "converged")
        assert not hasattr(sweep, "converged_n")

        assert hasattr(ca, "converged")
        assert hasattr(ca, "converged_n")

    def test_thread_independence(self):
        computer1 = ScenarioSweepCPP(self.grid)
        computer1.modify_load_p(self.load_p)
        computer1.set_contingency_lines(self.no_line_mask)
        computer1.compute(self.Vinit, self.max_it, self.tol)
        Vs1 = 1.0 * computer1.get_voltages()

        for nb_thread in (2, 4):
            computer = ScenarioSweepCPP(self.grid)
            computer.nb_thread = nb_thread
            computer.modify_load_p(self.load_p)
            computer.set_contingency_lines(self.no_line_mask)
            computer.compute(self.Vinit, self.max_it, self.tol)
            Vs = 1.0 * computer.get_voltages()
            assert np.array_equal(Vs, Vs1), (
                f"results differ with nb_thread={nb_thread} (max deviation "
                f"{np.max(np.abs(Vs - Vs1))})")

    def test_dc(self):
        from lightsim2grid.algorithm import AlgorithmType
        dc_algo = AlgorithmType.DC_SparseLU
        computer = ScenarioSweepCPP(self.grid)
        computer.change_algorithm(dc_algo)
        computer.modify_load_p(self.load_p)
        computer.set_contingency_lines(self.no_line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1

    def test_get_violations_requires_compute_limit_violations(self):
        """mirrors ContingencyAnalysisCPP's own test_default_is_off_and_raises"""
        computer = ScenarioSweepCPP(self.grid)
        computer.modify_load_p(self.load_p)
        computer.set_contingency_lines(self.no_line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        with self.assertRaises(RuntimeError):
            computer.get_violations()
        with self.assertRaises(RuntimeError):
            computer.get_violations_n()

    def test_violations_match_contingency_analysis(self):
        """a row with only a contingency (injections at grid defaults) must match the
        equivalent ContingencyAnalysis row, violations included"""
        line_id = 3
        ca = ContingencyAnalysisCPP(self.grid, True)
        ca.add_n1(line_id)
        ca.compute(self.Vinit, self.max_it, self.tol)
        ca_violations = ca.get_violations()[0]

        computer = ScenarioSweepCPP(self.grid)
        computer.compute_limit_violations = True
        line_mask = np.zeros((1, self.n_line), dtype=bool)
        line_mask[0, line_id] = True
        computer.set_contingency_lines(line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        ss_violations = computer.get_violations()[0]

        assert len(ss_violations) == len(ca_violations)
        for v_ss, v_ca in zip(ss_violations, ca_violations):
            assert v_ss.element_type == v_ca.element_type
            assert v_ss.violation_type == v_ca.violation_type
            assert v_ss.element_id == v_ca.element_id

    def test_violations_match_contingency_analysis_with_handle_disconnected_grid(self):
        """same as above, with handle_disconnected_grid=True on both sides -- also
        exercises the mask-mode _status fix (a compute() that only succeeds via
        masking must still report get_status() == 1, not 0)"""
        line_id = 3
        ca = ContingencyAnalysisCPP(self.grid, True)
        ca.handle_disconnected_grid = True
        ca.add_n1(line_id)
        ca.compute(self.Vinit, self.max_it, self.tol)
        ca_violations = ca.get_violations()[0]
        ca_voltages = 1.0 * ca.get_voltages()

        computer = ScenarioSweepCPP(self.grid)
        computer.compute_limit_violations = True
        computer.handle_disconnected_grid = True
        line_mask = np.zeros((1, self.n_line), dtype=bool)
        line_mask[0, line_id] = True
        computer.set_contingency_lines(line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1
        ss_violations = computer.get_violations()[0]
        ss_voltages = 1.0 * computer.get_voltages()

        assert np.max(np.abs(ss_voltages - ca_voltages)) <= 1e-8
        assert len(ss_violations) == len(ca_violations)

    def test_n_divergence_stamps_violations_not_empty_lists(self):
        """a diverging pre-batch "n" powerflow (an impossible tol forces this, with
        just enough max_iter to actually attempt and fail rather than never running)
        must stamp a GRID/DIVERGENCE sentinel into get_violations_n() and every row
        of get_violations(), not leave them as empty lists indistinguishable from
        "converged, no violations" -- there is no converged()/converged_n() escape
        hatch on ScenarioSweep to fall back on."""
        from lightsim2grid.lightsim2grid_cpp import ViolationElementType, LimitViolationType
        computer = ScenarioSweepCPP(self.grid)
        computer.compute_limit_violations = True
        computer.modify_load_p(self.load_p)
        computer.set_contingency_lines(self.no_line_mask)
        computer.compute(self.Vinit, max_iter=1, tol=1e-300)
        assert computer.get_status() == 0

        violations_n = computer.get_violations_n()
        assert len(violations_n) == 1
        assert violations_n[0].element_type == ViolationElementType.GRID
        assert violations_n[0].violation_type == LimitViolationType.DIVERGENCE

        for row in computer.get_violations():
            assert len(row) == 1
            assert row[0].element_type == ViolationElementType.GRID
            assert row[0].violation_type == LimitViolationType.DIVERGENCE

    def _pv_gen(self):
        gens = self.grid.get_generators()
        pv_gen_id = next(i for i, g in enumerate(gens) if not g.is_slack)
        return gens, pv_gen_id, gens[pv_gen_id].bus_id

    def test_modify_gen_v(self):
        """per-step generator target voltage magnitude, plain (non-masked) path."""
        gens, pv_gen_id, bus_id = self._pv_gen()
        nb_steps = 3
        targets = np.array([1.00, 1.03, 0.97])
        gen_v = np.tile(np.array([g.target_vm_pu for g in gens]), (nb_steps, 1))
        gen_v[:, pv_gen_id] = targets

        computer = ScenarioSweepCPP(self.grid)
        computer.modify_gen_p(self.gen_p[:nb_steps])
        computer.modify_load_p(self.load_p[:nb_steps])
        computer.modify_load_q(self.load_q[:nb_steps])
        computer.modify_gen_v(gen_v)
        computer.set_contingency_lines(self.no_line_mask[:nb_steps])
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1
        vm = np.abs(np.array(computer.get_voltages())[:, bus_id])
        np.testing.assert_allclose(vm, targets, atol=1e-6)

    def test_modify_gen_v_through_masked_path(self):
        """same as above, but with handle_disconnected_grid=True and an active (but
        non-islanding, on this meshed grid) contingency -- exercises
        _run_range_masked, the only per-step solve path that does not go through
        _run_one_step."""
        gens, pv_gen_id, bus_id = self._pv_gen()
        nb_steps = 3
        targets = np.array([1.00, 1.03, 0.97])
        gen_v = np.tile(np.array([g.target_vm_pu for g in gens]), (nb_steps, 1))
        gen_v[:, pv_gen_id] = targets

        computer = ScenarioSweepCPP(self.grid)
        computer.handle_disconnected_grid = True
        computer.modify_gen_p(self.gen_p[:nb_steps])
        computer.modify_load_p(self.load_p[:nb_steps])
        computer.modify_load_q(self.load_q[:nb_steps])
        computer.modify_gen_v(gen_v)
        line_mask = np.zeros((nb_steps, self.n_line), dtype=bool)
        line_mask[:, 0] = True  # a contingency that does not island the grid here
        computer.set_contingency_lines(line_mask)
        computer.compute(self.Vinit, self.max_it, self.tol)
        assert computer.get_status() == 1
        vm = np.abs(np.array(computer.get_voltages())[:, bus_id])
        np.testing.assert_allclose(vm, targets, atol=1e-6)


class TestScenarioSweepGrid2op(unittest.TestCase):
    """the grid2op-level wrapper (lightsim2grid.scenarioSweep.ScenarioSweep)"""

    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", backend=LightSimBackend(), test=True)
        self.env.set_id(0)
        self.obs = self.env.reset()

    def tearDown(self):
        self.env.close()

    def test_behaviour(self):
        sweep = ScenarioSweep(self.env)
        assert isinstance(sweep.computer, ScenarioSweepCPP)

        n_sim = 4
        load_p = np.tile(self.obs.load_p, (n_sim, 1))
        n_line = len(self.env.backend._grid.get_lines())
        line_mask = np.zeros((n_sim, n_line), dtype=bool)
        line_mask[2:, 0] = True

        sweep.modify_load_p(load_p)
        sweep.set_contingency_lines(line_mask)
        sweep.compute()
        Vs = sweep.get_voltages()
        Ps = sweep.compute_P()
        As = sweep.compute_A()
        assert Vs.shape[0] == n_sim
        assert Ps.shape[0] == n_sim
        assert As.shape[0] == n_sim
        # the disconnected line's flow is cleaned (reported as 0) for the masked rows
        assert np.all(As[2:, 0] == 0.0)
        assert np.all(As[:2, 0] != 0.0)

    def test_get_flows_before_compute_raises(self):
        sweep = ScenarioSweep(self.env)
        with self.assertRaises(RuntimeError):
            sweep.compute_A()
        with self.assertRaises(RuntimeError):
            sweep.compute_P()
        with self.assertRaises(RuntimeError):
            sweep.get_voltages()

    def test_nb_thread_validated(self):
        sweep = ScenarioSweep(self.env)
        with self.assertRaises(ValueError):
            sweep.nb_thread = 1.5

    def test_column_count_validated(self):
        sweep = ScenarioSweep(self.env)
        with self.assertRaises(RuntimeError):
            sweep.modify_load_p(np.zeros((3, self.env.n_load + 1)))
        with self.assertRaises(RuntimeError):
            sweep.set_contingency_lines(np.zeros((3, 2), dtype=bool))

    def test_run(self):
        """ScenarioSweep.run() -- reuses ContingencyResult/PreContingencyResult/
        SecurityAnalysisResult from contingencyAnalysis.py, same as
        ContingencyAnalysis.run() returns."""
        from lightsim2grid.scenarioSweep import SecurityAnalysisResult, ContingencyResult, PreContingencyResult

        _set_tight_limits(self.env.backend._grid)

        n_sim = 3
        load_p = np.tile(self.obs.load_p, (n_sim, 1))
        n_line = len(self.env.backend._grid.get_lines())
        line_mask = np.zeros((n_sim, n_line), dtype=bool)
        line_mask[1:, 2] = True

        sweep = ScenarioSweep(self.env)
        sweep.compute_limit_violations = True
        sweep.modify_load_p(load_p)
        sweep.set_contingency_lines(line_mask)
        res = sweep.run()

        assert isinstance(res, SecurityAnalysisResult)
        assert isinstance(res.pre_contingency_result, PreContingencyResult)
        assert len(res.post_contingency_results) == n_sim
        # tight limits guarantee at least one violation somewhere
        assert len(res.pre_contingency_result.limit_violations) > 0
        for row_id, row in enumerate(res.post_contingency_results):
            assert isinstance(row, ContingencyResult)
            assert row.contingency_name is None
            assert len(row.limit_violations) > 0
        assert res.post_contingency_results[0].element_ids == []
        assert res.post_contingency_results[1].element_ids == [2]
        assert res.post_contingency_results[2].element_ids == [2]

    def test_run_requires_the_flag(self):
        sweep = ScenarioSweep(self.env)
        sweep.modify_load_p(np.tile(self.obs.load_p, (2, 1)))
        with self.assertRaises(RuntimeError):
            sweep.run()


if __name__ == "__main__":
    unittest.main()
