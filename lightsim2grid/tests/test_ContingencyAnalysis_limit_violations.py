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


def _set_tight_limits(grid, vmin_kv=1e6, cur_ka=1e-6):
    """makes every bus voltage and every line/trafo current violate: guarantees at least
    one LOW_VOLTAGE and one CURRENT violation everywhere the corresponding element is used."""
    nb_bus = grid.get_bus_vn_kv().shape[0]
    grid.set_bus_voltage_limits(np.full(nb_bus, vmin_kv), np.full(nb_bus, np.nan))
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
        assert (ViolationElementType.BUS, LimitViolationType.LOW_VOLTAGE) in types_seen
        assert (ViolationElementType.LINE, LimitViolationType.CURRENT) in types_seen
        assert (ViolationElementType.TRAFO, LimitViolationType.CURRENT) in types_seen
        sub_names = self.env.backend._grid.get_substation_names()
        for v in viol_n:
            assert np.isfinite(v.value)
            assert np.isfinite(v.limit)
            if v.violation_type == LimitViolationType.LOW_VOLTAGE:
                assert v.value < v.limit
            else:
                assert v.value > v.limit
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


if __name__ == "__main__":
    unittest.main()
