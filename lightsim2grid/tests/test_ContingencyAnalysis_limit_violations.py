# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the (opt-in) limit-violation reporting of `ContingencyAnalysis` /
`ContingencyAnalysisCPP`: `compute_limit_violations`, `converged` / `converged_n`,
`get_violations` / `get_violations_n`, and the python-side `run` / `run_ac` / `run_dc`.
"""

import unittest
import warnings
import tempfile
import os
import numpy as np
import pandapower as pp
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (ContingencyAnalysisCPP,
                                              ViolationElementType,
                                              LimitViolationType)
from lightsim2grid.algorithm import AlgorithmType


def _set_tight_limits(grid, cur_ka=1e-6):
    """makes every bus voltage and every line/trafo current violate: guarantees at least one
    HIGH_VOLTAGE and one CURRENT violation everywhere the corresponding element is used.

    ``vmax == vn`` is violated by every bus at or above its nominal voltage, which on the
    l2rpn_case14_sandbox fixture is all of them (no bus sits below nominal there)."""
    vn = grid.get_bus_vn_kv()
    nb_bus = vn.shape[0]
    grid.set_bus_voltage_limits(np.full(nb_bus, np.nan), 1. * vn)
    n_line = len(grid.get_lines())
    n_trafo = len(grid.get_trafos())
    grid.set_line_current_limit_side1(np.full(n_line, cur_ka))
    grid.set_line_current_limit_side2(np.full(n_line, np.nan))
    grid.set_trafo_current_limit_side1(np.full(n_trafo, cur_ka))
    grid.set_trafo_current_limit_side2(np.full(n_trafo, np.nan))


class TestFlagGating(unittest.TestCase):
    """the whole feature must be opt-in: default False, and toggling it must not affect
    the pre-existing get_flows / get_voltages API at all."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_default_is_off_and_raises(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid)
        assert SA.compute_limit_violations is False
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        for fun_name in ("converged", "get_violations", "converged_n", "get_violations_n"):
            with self.assertRaises(RuntimeError):
                getattr(SA, fun_name)()

    def test_on_at_construction_works(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        assert SA.compute_limit_violations is True
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        SA.converged()
        SA.get_violations()
        SA.converged_n()
        SA.get_violations_n()

    def test_flows_unaffected_by_flag(self):
        lid_cont = [0, 1, 2, 3]
        SA_off = ContingencyAnalysisCPP(self.env.backend._grid, False)
        SA_off.add_multiple_n1(lid_cont)
        SA_off.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        flows_off = SA_off.compute_flows()

        SA_on = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA_on.add_multiple_n1(lid_cont)
        SA_on.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        flows_on = SA_on.compute_flows()

        assert np.array_equal(flows_off, flows_on, equal_nan=True)

    def test_setter_toggle_clears_and_raises_when_off(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        assert len(SA.get_violations()) == len(SA.my_defaults())

        SA.compute_limit_violations = False
        assert len(SA.my_defaults()) == 0, "toggling the flag should clear() the object"
        with self.assertRaises(RuntimeError):
            SA.get_violations()


class TestConvergedFlag(unittest.TestCase):
    """the `converged` flag must exactly reflect which contingencies actually produced a
    result (and get_violations() must stay empty, not fabricated, for the others)."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_converged_matches_nan_heuristic(self):
        # 17 converges, 18 diverges (splits the grid -> a pre-check skips it before the solver
        # is ever invoked, hence NOT_SIMULATED), 19 converges -- see
        # test_SecurityAnlysis_cpp.test_compute_nonconnected_graph for the reference behaviour
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        lid_cont = [17, 18, 19]
        SA.add_multiple_n1(lid_cont)
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        res_v = SA.get_voltages()
        res_flows = SA.compute_flows()
        converged = SA.converged()
        violations = SA.get_violations()
        for cont_id in range(len(lid_cont)):
            nan_heuristic_converged = np.all(np.isfinite(res_flows[cont_id])) and \
                np.max(np.abs(res_v[cont_id])) > 0.
            assert bool(converged[cont_id]) == nan_heuristic_converged, \
                f"converged() disagrees with the nan/0 heuristic for contingency {cont_id}"
            if not converged[cont_id]:
                viol = violations[cont_id]
                assert len(viol) == 1, \
                    "a non-converged contingency must report exactly one GRID violation, " \
                    f"got {viol}"
                assert viol[0].element_type == ViolationElementType.GRID
                assert viol[0].violation_type == LimitViolationType.NOT_SIMULATED, \
                    "contingency 18 splits the grid: a pre-check should skip it before the " \
                    "solver is invoked"


class TestBaseCaseDivergenceRaises(unittest.TestCase):
    """the pre-contingency ("n") case is special: unlike a single contingency, it cannot
    silently be reported as "no result available" (every contingency is solved relative to
    it), so a non-converging base case must raise instead."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_n_divergence_raises(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        nb_bus = self.env.backend._grid.get_bus_vn_kv().shape[0]
        # deliberately bad starting point + a single NR iteration: the base case cannot converge
        bad_vinit = np.full(nb_bus, 0.5, dtype=complex)
        with self.assertRaises(RuntimeError):
            SA.compute(bad_vinit, 1, self.env.backend.tol)


class TestThreadIndependence(unittest.TestCase):
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        _set_tight_limits(self.env.backend._grid)

    def tearDown(self):
        self.env.close()

    def test_converged_and_violation_counts_match(self):
        def _run(nb_thread):
            SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
            SA.nb_thread = nb_thread
            SA.add_all_n1()
            SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
            return SA.converged(), SA.get_violations()

        conv1, viol1 = _run(1)
        conv4, viol4 = _run(4)
        assert conv1 == conv4
        assert [len(v) for v in viol1] == [len(v) for v in viol4]


class TestBasicViolations(unittest.TestCase):
    """end-to-end: with deliberately tight limits, every kind of LimitViolationType /
    ViolationElementType must show up, with the correct fields."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        _set_tight_limits(self.env.backend._grid)

    def tearDown(self):
        self.env.close()

    def test_base_case_and_per_contingency_violations(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)

        assert SA.converged_n() is True
        viol_n = SA.get_violations_n()
        assert len(viol_n) > 0
        types_seen = {(v.element_type, v.violation_type) for v in viol_n}
        assert (ViolationElementType.BUS, LimitViolationType.HIGH_VOLTAGE) in types_seen
        assert (ViolationElementType.LINE, LimitViolationType.CURRENT) in types_seen
        assert (ViolationElementType.TRAFO, LimitViolationType.CURRENT) in types_seen
        sub_names = self.env.backend._grid.get_substation_names()
        for v in viol_n:
            assert np.isfinite(v.value)
            assert np.isfinite(v.limit)
            if v.violation_type == LimitViolationType.LOW_VOLTAGE:
                assert v.value <= v.limit
            else:
                assert v.value >= v.limit
            if v.element_type == ViolationElementType.BUS:
                # `name` is the *substation* the violating bus belongs to (no per-bus name
                # exists in LSGrid), not the bus id itself
                assert v.name == sub_names[v.element_id]

        violations = SA.get_violations()
        assert sum(len(v) for v in violations) > 0

    def test_disconnected_branch_of_its_own_contingency_is_skipped(self):
        # a contingency must never report a CURRENT violation on the very branch it disconnects
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_n1(0)
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        viol0 = SA.get_violations()[0]
        assert not any(v.element_type == ViolationElementType.LINE and v.element_id == 0
                       for v in viol0)


class TestMaskedBusExclusion(unittest.TestCase):
    """with `handle_disconnected_grid=True`, a masked (forced-to-0V) bus must never be
    reported as a LOW_VOLTAGE violation."""
    def setUp(self):
        self.max_it = 30
        self.tol = 1e-8
        net = pp.create_empty_network(sn_mva=1.)
        for _ in range(4):
            pp.create_bus(net, vn_kv=20.)
        pp.create_ext_grid(net, 0, vm_pu=1.0)
        pp.create_line_from_parameters(net, 0, 1, length_km=1., r_ohm_per_km=0.1,
                                       x_ohm_per_km=0.3, c_nf_per_km=0., max_i_ka=1.)  # line 0
        pp.create_line_from_parameters(net, 1, 2, length_km=1., r_ohm_per_km=0.1,
                                       x_ohm_per_km=0.3, c_nf_per_km=0., max_i_ka=1.)  # line 1
        pp.create_line_from_parameters(net, 2, 3, length_km=1., r_ohm_per_km=0.1,
                                       x_ohm_per_km=0.3, c_nf_per_km=0., max_i_ka=1.)  # line 2 (splits)
        pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
        pp.create_load(net, 2, p_mw=1.0, q_mvar=0.2)
        pp.create_load(net, 3, p_mw=1.0, q_mvar=0.2)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(net)
        self.V0 = np.ones(self.grid.get_bus_vn_kv().shape[0], dtype=complex)
        # any positive vmin_kv would make the masked (0V) bus 3 violate if it were checked
        nb_bus = self.grid.get_bus_vn_kv().shape[0]
        self.grid.set_bus_voltage_limits(np.full(nb_bus, 1.), np.full(nb_bus, np.nan))

    def test_masked_bus_not_reported(self):
        SA = ContingencyAnalysisCPP(self.grid, True)
        SA.add_n1(2)  # disconnect line 2 -> isolates bus 3
        SA.handle_disconnected_grid = True
        SA.compute(1. * self.V0, self.max_it, self.tol)
        v = SA.get_voltages()
        assert v[0, 3] == 0., "sanity check: bus 3 should be masked (0V)"
        assert SA.converged()[0] is True
        viol0 = SA.get_violations()[0]
        assert not any(v.element_type == ViolationElementType.BUS and v.element_id == 3
                       for v in viol0), "the masked bus must be excluded from voltage checks"


class TestDCPhaseShiftSide2(unittest.TestCase):
    """side-2 current computation (only exercised by limit-violation checking, the existing
    batch API only ever computes side 1): build a grid with a phase-shifting transformer (so
    dc_x_tau_shift != 0) and check both sides trigger sensible, finite, positive violations
    for a very tight limit, in both AC and DC."""
    def _build(self):
        case = pn.case14()
        hv_bus, lv_bus = 0, 2
        pp.create_transformer_from_parameters(
            case, hv_bus=hv_bus, lv_bus=lv_bus, sn_mva=9900.0,
            vn_hv_kv=case.bus.iloc[hv_bus]["vn_kv"], vn_lv_kv=case.bus.iloc[lv_bus]["vn_kv"],
            i0_percent=0.0, vk_percent=2070.288000, vkr_percent=0.0,
            shift_degree=-10.0, pfe_kw=0., tap_side="hv")
        case.name = "case14_shift"
        with tempfile.TemporaryDirectory() as path:
            case_name = os.path.join(path, "this_case.json")
            pp.to_json(case, case_name)
            from lightsim2grid import LightSimBackend
            backend = LightSimBackend()
            type(backend)._clear_grid_dependant_class_attributes()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                type(backend).env_name = case.name
                backend.load_grid(case_name)
                backend.assert_grid_correct()
        backend._init_pp_backend._grid.ext_grid["in_service"] = False
        backend._init_pp_backend._grid.ext_grid.loc[0, "in_service"] = True
        backend._init_pp_backend._grid.gen["slack"] = False
        return backend

    def _check_both_sides_trigger(self, is_dc):
        backend = self._build()
        conv, exc_ = backend.runpf(is_dc=is_dc)
        assert conv, f"reference powerflow did not converge: {exc_}"
        grid = backend._grid
        trafos = grid.get_trafos()
        idx = len(trafos) - 1  # the phase-shifting trafo we just added

        n_trafo = len(trafos)
        lim1 = np.full(n_trafo, np.nan)
        lim2 = np.full(n_trafo, np.nan)
        lim1[idx] = 1e-9  # any nonzero current on this trafo will violate
        lim2[idx] = 1e-9
        grid.set_trafo_current_limit_side1(lim1)
        grid.set_trafo_current_limit_side2(lim2)

        SA = ContingencyAnalysisCPP(grid, True)
        if is_dc:
            SA.change_algorithm(AlgorithmType.DC_KLU)
        SA.add_all_n1()
        SA.compute(backend.V.copy(), 10, 1e-8)

        assert SA.converged_n() is True
        viol_n = {v.side: v for v in SA.get_violations_n() if v.element_id == idx}
        assert 1 in viol_n, "side 1 of the phase-shifting trafo should violate"
        assert 2 in viol_n, "side 2 of the phase-shifting trafo should violate"
        for side, v in viol_n.items():
            assert v.element_type == ViolationElementType.TRAFO
            assert v.violation_type == LimitViolationType.CURRENT
            assert np.isfinite(v.value) and v.value > 0.

    def test_ac(self):
        self._check_both_sides_trigger(is_dc=False)

    def test_dc(self):
        self._check_both_sides_trigger(is_dc=True)


class TestPythonWrapperRun(unittest.TestCase):
    """end-to-end test of the python-side `ContingencyAnalysis.run` API (pypowsybl-like
    shape: res.pre_contingency_result.limit_violations / res.post_contingency_results)."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        _set_tight_limits(self.env.backend._grid)

    def tearDown(self):
        self.env.close()

    def test_run_requires_the_flag(self):
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env)
        sa.add_single_contingency(0)
        with self.assertRaises(RuntimeError):
            sa.run()

    def test_run_result_shape_order_and_names(self):
        from lightsim2grid.contingencyAnalysis import (ContingencyAnalysis, SecurityAnalysisResult,
                                                         PreContingencyResult, ContingencyResult)
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_single_contingency(4, name="fourth")
        sa.add_single_contingency(0, name="first")
        sa.add_single_contingency(7)  # no name

        res = sa.run()
        assert isinstance(res, SecurityAnalysisResult)
        assert isinstance(res.pre_contingency_result, PreContingencyResult)
        assert res.pre_contingency_result.converged is True
        assert len(res.pre_contingency_result.limit_violations) > 0

        assert isinstance(res.post_contingency_results, list)
        assert len(res.post_contingency_results) == 3
        # order must match insertion order (4, 0, 7), not the c++-internal order
        expected_ids = [[4], [0], [7]]
        expected_names = ["fourth", "first", None]
        for cont, exp_ids, exp_name in zip(res.post_contingency_results, expected_ids, expected_names):
            assert isinstance(cont, ContingencyResult)
            assert cont.element_ids == exp_ids
            assert cont.element_names == [str(self.env.name_line[el_id]) for el_id in exp_ids]
            assert cont.contingency_name == exp_name
            assert cont.converged is True
            assert len(cont.limit_violations) > 0

    def test_run_ac_and_run_dc(self):
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_single_contingency(0)
        res_ac = sa.run_ac()
        assert len(res_ac.post_contingency_results) == 1

        # switching family clears the registered contingencies (documented behaviour,
        # mirrors change_algorithm) -- switch first, then (re-)register, then run
        sa._change_algorithm_family(want_dc=True)
        sa.add_single_contingency(0)
        res_dc = sa.run_dc()
        assert len(res_dc.post_contingency_results) == 1


class TestViolationThreshold(unittest.TestCase):
    """`violation_threshold`: default value, validation range, and the fact that the
    default (1.) reproduces the legacy, threshold-less behaviour."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_default_is_one(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        assert SA.violation_threshold == 1.0

    def test_default_reproduces_legacy_behaviour(self):
        # the legacy checks are strict (`>` / `<`); with threshold == 1. they become
        # `>=` / `<=`, which only differ if a value lands *exactly* on its limit. With
        # the deliberately extreme limits of `_set_tight_limits` nothing does, so the
        # reported violations must be exactly the ones the pre-threshold code reported.
        _set_tight_limits(self.env.backend._grid)
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        viol_n = SA.get_violations_n()
        assert len(viol_n) > 0
        types_seen = {(v.element_type, v.violation_type) for v in viol_n}
        assert (ViolationElementType.BUS, LimitViolationType.HIGH_VOLTAGE) in types_seen
        assert (ViolationElementType.LINE, LimitViolationType.CURRENT) in types_seen
        assert (ViolationElementType.TRAFO, LimitViolationType.CURRENT) in types_seen

    def test_rejects_out_of_range(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        for bad in (0.0, -0.5, 1.0001, 2.0, float("nan")):
            with self.assertRaises(RuntimeError):
                SA.violation_threshold = bad
            assert SA.violation_threshold == 1.0, "a rejected value must not be applied"

    def test_accepts_boundary_one(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.violation_threshold = 0.5
        SA.violation_threshold = 1.0  # the upper bound is inclusive
        assert SA.violation_threshold == 1.0

    def test_python_wrapper_validates_and_raises_value_error(self):
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        for bad in (0.0, -1.0, 1.5, "not_a_number", float("nan")):
            with self.assertRaises(ValueError):
                sa.violation_threshold = bad
        sa.violation_threshold = 0.9
        assert sa.violation_threshold == 0.9


class TestViolationThresholdBoundary(unittest.TestCase):
    """the actual semantics of `violation_threshold`, one "boundary straddling" case per
    violation type: the limit is placed 1% away from the value actually reached, so
    threshold == 1. reports nothing and a threshold below the crossover reports a violation.

    The reference value is not recomputed here (which would only test that two
    reimplementations of the same formula agree): it is harvested from a first run with a
    deliberately unreachable limit, so it is *by construction* the exact value the checker
    itself computes.
    """
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        self.grid = self.env.backend._grid
        self.nb_bus = self.grid.get_bus_vn_kv().shape[0]
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())
        self._no_limits()

    def tearDown(self):
        self.env.close()

    def _no_limits(self):
        """clears every limit, so each test only sets the single one it is about"""
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, np.nan), np.full(self.nb_bus, np.nan))
        self.grid.set_line_current_limit_side1(np.full(self.n_line, np.nan))
        self.grid.set_line_current_limit_side2(np.full(self.n_line, np.nan))
        self.grid.set_trafo_current_limit_side1(np.full(self.n_trafo, np.nan))
        self.grid.set_trafo_current_limit_side2(np.full(self.n_trafo, np.nan))

    def _violations_n(self, threshold):
        SA = ContingencyAnalysisCPP(self.grid, True)
        SA.violation_threshold = threshold
        SA.add_n1(0)
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        return SA.get_violations_n()

    def test_current_threshold_straddle(self):
        line_id = 0
        # harvest the exact current the checker computes on this line (unreachable limit)
        lim1 = np.full(self.n_line, np.nan)
        lim1[line_id] = 1e-9
        self.grid.set_line_current_limit_side1(lim1)
        harvested = [v for v in self._violations_n(1.0)
                     if v.element_type == ViolationElementType.LINE and v.element_id == line_id]
        assert len(harvested) == 1
        amps = harvested[0].value

        # place the limit 1% above it: no violation at threshold == 1.
        lim1[line_id] = amps * 1.01
        self.grid.set_line_current_limit_side1(lim1)

        def violates(th):
            return any(v.element_type == ViolationElementType.LINE and v.element_id == line_id
                       and v.side == 1 for v in self._violations_n(th))

        assert not violates(1.0), "amps < limit: nothing to report at the nominal threshold"
        # threshold * limit == 0.98 * 1.01 * amps == 0.9898 * amps <= amps -> violates
        assert violates(0.98)

    def _harvest_bus_voltage(self):
        """returns (bus_id, v_kv, vnom_kv) for the bus furthest ABOVE its nominal voltage, as
        computed by the checker itself. `vmax == vn` is the tightest upper bound a grid may
        declare (the band must bracket vn), and every bus at or above nominal violates it, so
        it doubles as a probe reporting each such bus's exact voltage."""
        vn = self.grid.get_bus_vn_kv()
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, np.nan), 1. * vn)
        harvested = {v.element_id: v.value for v in self._violations_n(1.0)
                     if v.element_type == ViolationElementType.BUS}
        assert len(harvested) > 0
        bus_id = max(harvested, key=lambda b: harvested[b] / vn[b])
        assert harvested[bus_id] > vn[bus_id], "expected a bus above its nominal voltage"
        return bus_id, harvested[bus_id], vn[bus_id]

    def test_high_voltage_lerp_crossover_at_half(self):
        """the effective bound is the linear interpolation
        `threshold * v_max + (1 - threshold) * v_nom`. Placing v_max at exactly TWICE the
        bus's excess over nominal therefore puts the crossover at exactly threshold == 0.5,
        which is what this pins down (a rule that merely got "stricter somehow" would not
        land there)."""
        bus_id, v0, vnom = self._harvest_bus_voltage()
        excess = v0 - vnom
        vmax = np.full(self.nb_bus, np.nan)
        vmax[bus_id] = vnom + 2. * excess   # -> high_eff(th) == vnom + th * 2 * excess
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, np.nan), vmax)

        def violates(th):
            return any(v.element_type == ViolationElementType.BUS and v.element_id == bus_id
                       and v.violation_type == LimitViolationType.HIGH_VOLTAGE
                       for v in self._violations_n(th))

        assert not violates(1.0), "v is below v_max: nothing to report at the nominal threshold"
        assert not violates(0.55), "high_eff has not come down to v yet"
        assert violates(0.5), "high_eff == vnom + excess == v exactly, and the test is `>=`"
        assert violates(0.45)

    def test_bus_above_nominal_is_never_low_voltage(self):
        """the anti-saturation property: because the low bound is interpolated towards
        v_nom, it can never go past it, so a bus sitting ABOVE its nominal voltage is never
        reported as LOW_VOLTAGE -- however small the threshold gets. The previous
        bound-scaling rule (`v <= v_min / threshold`) flagged every such bus instead."""
        bus_id, v0, vnom = self._harvest_bus_voltage()
        assert v0 > vnom
        vmin = np.full(self.nb_bus, np.nan)
        vmin[bus_id] = vnom * 0.95   # a perfectly ordinary lower bound, below nominal
        self.grid.set_bus_voltage_limits(vmin, np.full(self.nb_bus, np.nan))

        for th in (1.0, 0.9, 0.5, 0.1, 1e-9):
            assert not any(v.element_type == ViolationElementType.BUS and v.element_id == bus_id
                           and v.violation_type == LimitViolationType.LOW_VOLTAGE
                           for v in self._violations_n(th)), \
                f"a bus above nominal must never be LOW_VOLTAGE (threshold={th})"

    def test_bounds_never_cross_so_low_and_high_are_exclusive(self):
        """with a band bracketing v_nom (the only kind a grid may declare) the two effective
        bounds each converge towards v_nom from their own side, so they can never cross: no
        bus is ever reported as both too low and too high, however small the threshold."""
        vn = self.grid.get_bus_vn_kv()
        self.grid.set_bus_voltage_limits(0.999 * vn, 1.001 * vn)  # deliberately very narrow

        for th in (1.0, 0.99, 0.9, 0.5, 0.2, 1e-9):
            per_bus = {}
            for v in self._violations_n(th):
                if v.element_type != ViolationElementType.BUS:
                    continue
                per_bus.setdefault(v.element_id, []).append(v.violation_type)
            for bus_id, kinds in per_bus.items():
                assert len(kinds) == 1, (
                    f"bus {bus_id} reported {len(kinds)} voltage violations at threshold "
                    f"{th} ({kinds}): the effective bounds must never cross")

    def test_thread_independence_under_non_default_threshold(self):
        # `_violation_threshold_` is read directly (not passed down) by every worker
        # thread: check the results really are thread-count independent
        _set_tight_limits(self.grid)

        def _run(nb_thread):
            SA = ContingencyAnalysisCPP(self.grid, True)
            SA.violation_threshold = 0.5
            SA.nb_thread = nb_thread
            SA.add_all_n1()
            SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
            return SA.converged(), SA.get_violations()

        conv1, viol1 = _run(1)
        conv4, viol4 = _run(4)
        assert conv1 == conv4
        assert [len(v) for v in viol1] == [len(v) for v in viol4]


class TestViolationThresholdLowVoltage(unittest.TestCase):
    """the LOW_VOLTAGE side of the interpolation needs buses actually sagging BELOW their
    nominal voltage (the l2rpn_case14_sandbox fixture has none), so this uses a small radial
    feeder whose far end droops under load."""
    def setUp(self):
        self.max_it, self.tol = 30, 1e-8
        net = pp.create_empty_network(sn_mva=1.)
        for _ in range(4):
            pp.create_bus(net, vn_kv=20.)
        pp.create_ext_grid(net, 0, vm_pu=1.0)
        for i in range(3):
            pp.create_line_from_parameters(net, i, i + 1, length_km=1., r_ohm_per_km=0.1,
                                           x_ohm_per_km=0.3, c_nf_per_km=0., max_i_ka=1.)
        for b in (1, 2, 3):
            pp.create_load(net, b, p_mw=1.0, q_mvar=0.2)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(net)
        self.nb_bus = self.grid.get_bus_vn_kv().shape[0]
        self.V0 = np.ones(self.nb_bus, dtype=complex)

    def _violations_n(self, threshold):
        SA = ContingencyAnalysisCPP(self.grid, True)
        SA.violation_threshold = threshold
        SA.add_n1(0)
        SA.compute(1. * self.V0, self.max_it, self.tol)
        return SA.get_violations_n()

    def test_low_voltage_lerp_crossover_at_half(self):
        """mirror of `test_high_voltage_lerp_crossover_at_half`: the effective bound is
        `threshold * v_min + (1 - threshold) * v_nom`, so placing v_min at exactly TWICE the
        bus's sag below nominal puts the crossover at exactly threshold == 0.5."""
        # harvest the sagging bus: `vmin == vn` is the tightest lower bound a grid may
        # declare, and every bus at or below nominal violates it, so it reports their values
        vn = self.grid.get_bus_vn_kv()
        self.grid.set_bus_voltage_limits(1. * vn, np.full(self.nb_bus, np.nan))
        harvested = {v.element_id: v.value for v in self._violations_n(1.0)
                     if v.element_type == ViolationElementType.BUS}
        bus_id = min(harvested, key=lambda b: harvested[b] / vn[b])
        v0, vnom = harvested[bus_id], vn[bus_id]
        assert v0 < vnom, "the far end of the feeder should sag below nominal"

        sag = vnom - v0
        vmin = np.full(self.nb_bus, np.nan)
        vmin[bus_id] = vnom - 2. * sag   # -> low_eff(th) == vnom - th * 2 * sag
        self.grid.set_bus_voltage_limits(vmin, np.full(self.nb_bus, np.nan))

        def violates(th):
            return any(v.element_type == ViolationElementType.BUS and v.element_id == bus_id
                       and v.violation_type == LimitViolationType.LOW_VOLTAGE
                       for v in self._violations_n(th))

        assert not violates(1.0), "v is above v_min: nothing to report at the nominal threshold"
        assert not violates(0.55), "low_eff has not come up to v yet"
        assert violates(0.5), "low_eff == vnom - sag == v exactly, and the test is `<=`"
        assert violates(0.45)


    def test_low_and_high_checks_are_independent(self):
        """each voltage check is anchored on v_nom, NOT on the opposite bound, so the
        LOW_VOLTAGE verdict is a function of (v_min, v_nom, threshold) alone: setting v_max,
        changing it, or leaving it at NaN can never alter it -- and symmetrically. This grid
        is the one with buses on both sides of nominal, so both directions are exercised."""
        vn = self.grid.get_bus_vn_kv()
        nan = np.full(self.nb_bus, np.nan)

        def sets(vmin_pu, vmax_pu, th):
            self.grid.set_bus_voltage_limits(vmin_pu * vn if vmin_pu else nan,
                                             vmax_pu * vn if vmax_pu else nan)
            viol = self._violations_n(th)
            low = {v.element_id for v in viol if v.element_type == ViolationElementType.BUS
                   and v.violation_type == LimitViolationType.LOW_VOLTAGE}
            high = {v.element_id for v in viol if v.element_type == ViolationElementType.BUS
                    and v.violation_type == LimitViolationType.HIGH_VOLTAGE}
            return low, high

        for th in (1.0, 0.9, 0.5):
            # v_min fixed at nominal (it catches the sagging buses), v_max varied
            ref_low, _ = sets(1.0, None, th)
            assert len(ref_low) > 0, "the fixture should produce some LOW_VOLTAGE violations"
            for vmax_pu in (1.05, 1.20, 5.00):
                low, _ = sets(1.0, vmax_pu, th)
                assert low == ref_low, \
                    f"LOW_VOLTAGE changed when v_max was set to {vmax_pu} pu (threshold={th})"

            # ... and symmetrically: v_max fixed at nominal, v_min varied
            _, ref_high = sets(None, 1.0, th)
            assert len(ref_high) > 0, "the slack bus sits exactly at nominal, so it violates"
            for vmin_pu in (0.95, 0.90, 0.10):
                _, high = sets(vmin_pu, 1.0, th)
                assert high == ref_high, \
                    f"HIGH_VOLTAGE changed when v_min was set to {vmin_pu} pu (threshold={th})"

    def test_fully_degenerate_band_reports_low_voltage(self):
        """the one -- measure-zero -- exception to "a bus is never both too low and too
        high": the band collapsed onto nominal (v_min == v_nom == v_max) makes both effective
        bounds equal to v_nom, so a bus sitting EXACTLY at nominal satisfies `v <= low_eff`
        and `v >= high_eff` at once. It is still reported once, as LOW_VOLTAGE (the checks
        are an `if low ... else if high` chain). Pinned here so the precedence is a decision
        rather than an accident."""
        vn = self.grid.get_bus_vn_kv()
        self.grid.set_bus_voltage_limits(1. * vn, 1. * vn)
        kinds = {}
        for v in self._violations_n(1.0):
            if v.element_type == ViolationElementType.BUS:
                kinds.setdefault(v.element_id, []).append(v.violation_type)
        # bus 0 is the slack, held exactly at 1.0 pu == v_nom
        assert kinds[0] == [LimitViolationType.LOW_VOLTAGE]
        for bus_id, k in kinds.items():
            assert len(k) == 1, f"bus {bus_id} was reported twice ({k})"

class TestBandNotBracketingNominal(unittest.TestCase):
    """A band is NOT required to bracket the bus nominal voltage: real IIDM data violates it
    (a 380 kV level with an operating range of [390, 450] kV is ordinary on the European
    400 kV network), so such a grid must still load and still be analysable. The threshold
    copes by clamping its interpolation anchor into [vmin, vmax]."""
    def setUp(self):
        self.max_it, self.tol = 30, 1e-8
        net = pp.create_empty_network(sn_mva=1.)
        for _ in range(2):
            pp.create_bus(net, vn_kv=380.)
        pp.create_ext_grid(net, 0, vm_pu=1.08)   # operated well above nominal, as such a level is
        pp.create_line_from_parameters(net, 0, 1, length_km=1., r_ohm_per_km=0.01,
                                       x_ohm_per_km=0.03, c_nf_per_km=0., max_i_ka=10.)
        pp.create_load(net, 1, p_mw=1.0, q_mvar=0.2)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.grid = init_from_pandapower(net)
        self.nb_bus = self.grid.get_bus_vn_kv().shape[0]
        self.vn = self.grid.get_bus_vn_kv()
        self.V0 = np.ones(self.nb_bus, dtype=complex)

    def _violations_n(self, threshold):
        SA = ContingencyAnalysisCPP(self.grid, True)
        SA.violation_threshold = threshold
        SA.add_n1(0)
        SA.compute(1. * self.V0, self.max_it, self.tol)
        return SA.get_violations_n()

    def _bus_kinds(self, threshold):
        out = {}
        for v in self._violations_n(threshold):
            if v.element_type == ViolationElementType.BUS:
                out.setdefault(v.element_id, []).append(v.violation_type)
        return out

    def test_such_a_band_is_accepted_and_round_trips(self):
        # a realistic shape: nominal 380, operating range [390, 450]
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, 390.), np.full(self.nb_bus, 450.))
        self.grid.check_grid()
        assert np.allclose(self.grid.get_bus_vmin_kv(), 390.)
        cls = type(self.grid)
        restored = cls.__new__(cls)
        restored.__setstate__(self.grid.__getstate__())
        restored.check_grid()
        assert np.allclose(restored.get_bus_vmin_kv(), 390.)

    def test_bounds_still_never_cross(self):
        """with the anchor clamped to vmin (= 390 here) the low bound has no room to move and
        the high bound interpolates down towards it, so they converge instead of crossing."""
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, 390.), np.full(self.nb_bus, 450.))
        for th in (1.0, 0.9, 0.5, 0.1, 1e-9):
            for bus_id, kinds in self._bus_kinds(th).items():
                assert len(kinds) == 1, \
                    f"bus {bus_id} reported {kinds} at threshold {th}: bounds must not cross"

    def test_low_bound_is_pinned_when_the_anchor_clamps_to_it(self):
        """anchor == vmin, so `low_eff = vmin + (1-th) * (anchor - vmin) == vmin` for every
        threshold: the low side simply has no margin left to give (it is already at its own
        bound), rather than drifting the wrong way."""
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, 390.), np.full(self.nb_bus, 450.))
        # bus 1 sags just below the slack; both sit above 390, so neither is LOW at th == 1
        assert not any(k == [LimitViolationType.LOW_VOLTAGE] for k in self._bus_kinds(1.0).values())
        # ... and no threshold can make them LOW either, since low_eff stays at 390
        for th in (0.9, 0.5, 1e-9):
            for bus_id, kinds in self._bus_kinds(th).items():
                assert LimitViolationType.LOW_VOLTAGE not in kinds, \
                    f"bus {bus_id} became LOW_VOLTAGE at threshold {th}; low_eff must stay at vmin"

    def test_high_bound_still_tightens_monotonically(self):
        """the high side keeps a real interval to shrink ([anchor, vmax] == [390, 450]), so
        lowering the threshold still reports more, never fewer, violations."""
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, 390.), np.full(self.nb_bus, 450.))
        seen = None
        for th in (1.0, 0.9, 0.5, 0.2, 1e-9):
            high = {b for b, k in self._bus_kinds(th).items()
                    if LimitViolationType.HIGH_VOLTAGE in k}
            assert seen is None or seen <= high, \
                f"a bus stopped being HIGH_VOLTAGE when the threshold fell to {th}"
            seen = high
        assert len(seen) > 0, "at a vanishing threshold the whole band should be violated"

    def test_inconsistent_band_still_raises(self):
        """vmin > vmax remains a genuine input error (unlike a band that merely fails to
        bracket nominal), and is still reported as one rather than as an arbitrary violation."""
        self.grid.set_bus_voltage_limits(np.full(self.nb_bus, 450.), np.full(self.nb_bus, 390.))
        with self.assertRaises(RuntimeError):
            self._violations_n(1.0)


class TestViolationThresholdInvalidation(unittest.TestCase):
    """*lowering* the threshold makes every check stricter, so results computed under the
    previous (looser) one would silently under-report: they must be dropped. Raising it
    back up must not drop anything (a stricter result is a superset)."""
    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})
        _set_tight_limits(self.env.backend._grid)

    def tearDown(self):
        self.env.close()

    def test_lowering_clears_results_but_keeps_contingencies(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        nb_cont = len(SA.my_defaults())
        assert nb_cont > 0
        assert len(SA.get_violations()) == nb_cont
        assert len(SA.get_violations_n()) > 0

        SA.violation_threshold = 0.99

        assert len(SA.my_defaults()) == nb_cont, "the contingencies themselves must be kept"
        assert len(SA.get_violations()) == 0, "results computed at the looser threshold must be dropped"
        assert len(SA.converged()) == 0
        assert len(SA.get_violations_n()) == 0

        # ... and recomputing (no need to re-add anything) brings them back
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        assert len(SA.get_violations()) == nb_cont

    def test_raising_does_not_clear(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.violation_threshold = 0.5
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        nb_cont = len(SA.my_defaults())
        nb_viol_n = len(SA.get_violations_n())
        assert nb_cont > 0 and nb_viol_n > 0

        SA.violation_threshold = 0.8  # loosening: a stricter result is still a valid superset

        assert len(SA.my_defaults()) == nb_cont
        assert len(SA.get_violations()) == nb_cont
        assert len(SA.get_violations_n()) == nb_viol_n

    def test_setting_the_same_value_does_not_clear(self):
        SA = ContingencyAnalysisCPP(self.env.backend._grid, True)
        SA.add_all_n1()
        SA.compute(self.env.backend.V, self.env.backend.max_it, self.env.backend.tol)
        nb_cont = len(SA.my_defaults())
        SA.violation_threshold = 1.0  # same as the current value: not a tightening
        assert len(SA.get_violations()) == nb_cont

    def test_python_wrapper_reruns_after_tightening(self):
        # without the python-side cache reset mirroring the c++ `clear_results_only()`,
        # `__computed` would stay True while the c++ results are gone, and this second
        # `run()` would read empty result arrays instead of recomputing
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_all_n1_contingencies()
        res1 = sa.run()
        nb_cont = len(res1.post_contingency_results)
        assert nb_cont > 0

        sa.violation_threshold = 0.9
        res2 = sa.run()
        assert len(res2.post_contingency_results) == nb_cont
        assert len(res2.pre_contingency_result.limit_violations) > 0
        # stricter threshold => at least as many violations as before, on every contingency
        for c1, c2 in zip(res1.post_contingency_results, res2.post_contingency_results):
            assert c1.element_ids == c2.element_ids
            if c1.converged and c2.converged:
                assert len(c2.limit_violations) >= len(c1.limit_violations)


if __name__ == "__main__":
    unittest.main()
