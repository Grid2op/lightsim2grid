# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the (opt-in) reactive-capability reporting of the batch algorithms:
`compute_bus_q_violations` / `bus_q_violation_tol_mvar` / `get_bus_q_violations` /
`get_bus_q_violations_n`, on the C++ classes and on the python wrappers
(`TimeSerie`, `InjectionSweep`, `ContingencyAnalysis`, `ScenarioSweep`).

The reported VALUES are pinned against a single-shot `ac_pf` here as they are in
`src/tests/test_batch_bus_q_violations.cpp`; what this file adds is the python layer --
the properties, what they invalidate, and the `bus_q_violations` field of the
`run()` result.
"""

import unittest
import warnings
import numpy as np
import pandapower.networks as pn

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import (TimeSeriesCPP, ContingencyAnalysisCPP,
                                             ScenarioSweepCPP,
                                             ViolationCategory, ViolationElementType,
                                             LimitViolationType, violation_category)
from lightsim2grid.algorithm import AlgorithmType


def _case14_tight_q(max_q_mvar=5.):
    """case14, with every generator's reactive range narrowed to +/- `max_q_mvar` -- so the
    buses they hold ask for far more reactive power than they own. Reactive limits have no
    influence at all on the solution (nothing enforces them), so this changes only what is
    reported, never what is solved."""
    net = pn.case14()
    net.gen["min_q_mvar"] = -max_q_mvar
    net.gen["max_q_mvar"] = max_q_mvar
    if "min_q_mvar" in net.ext_grid:
        net.ext_grid["min_q_mvar"] = -max_q_mvar
        net.ext_grid["max_q_mvar"] = max_q_mvar
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pandapower(net)


def _reference_bus_q(grid):
    """what a single ac_pf publishes, summed per bus over the voltage-regulating generators:
    exactly what `get_bus_q_violations` reports a value for"""
    v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
    grid.ac_pf(v_init, 20, 1e-11)
    per_bus = {}
    for gen in grid.get_generators():
        if not gen.connected or not gen.voltage_regulator_on:
            continue
        per_bus[gen.bus_id] = per_bus.get(gen.bus_id, 0.) + gen.res_q_mvar
    return per_bus


def _one_row_time_series(grid):
    """a TimeSeriesCPP with a single row: the grid's own generator set-points, so the row
    solves the very state `ac_pf` solves"""
    ts = TimeSeriesCPP(grid)
    ts.compute_bus_q_violations = True
    ts.bus_q_violation_tol_mvar = 0.
    gen_p = np.array([[gen.target_p_mw for gen in grid.get_generators()]])
    ts.modify_gen_p(gen_p)
    return ts


class TestBusQViolationsCpp(unittest.TestCase):
    """the C++ classes, through their bindings"""

    def test_default_is_off_and_raises(self):
        grid = _case14_tight_q()
        for cls in (TimeSeriesCPP, ContingencyAnalysisCPP, ScenarioSweepCPP):
            algo = cls(grid)
            assert algo.compute_bus_q_violations is False, cls.__name__
            assert algo.bus_q_violation_tol_mvar == 1e-4, cls.__name__
            with self.assertRaises(RuntimeError):
                algo.get_bus_q_violations()
            with self.assertRaises(RuntimeError):
                algo.get_bus_q_violations_n()

    def test_reports_what_ac_pf_publishes_per_bus(self):
        grid = _case14_tight_q()
        expected = _reference_bus_q(grid)
        assert len(expected) > 1, "the fixture is only meaningful with several regulated buses"

        ts = _one_row_time_series(grid)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]

        viols = ts.get_bus_q_violations()
        assert len(viols) == 1, "one row in, one row out"
        assert len(viols[0]) > 0, "+/- 5 MVAr per machine cannot hold case14's voltages"
        for v in viols[0]:
            assert v.element_type == ViolationElementType.BUS
            assert v.violation_type in (LimitViolationType.LOW_Q, LimitViolationType.HIGH_Q)
            # NOT an operational limit: the whole point of the category
            assert v.category == ViolationCategory.PHYSICAL
            assert v.element_id in expected, f"bus {v.element_id} holds no regulating generator"
            self.assertAlmostEqual(v.value, expected[v.element_id], places=5)
            # the limit is the SUM over that bus' machines, so a multiple of 5 MVAr here
            self.assertAlmostEqual(abs(v.limit) % 5., 0., places=9)
        # the base ("n") case solves the same state, so it reports the same thing
        assert len(ts.get_bus_q_violations_n()) == len(viols[0])

    def test_a_wide_tolerance_hides_everything(self):
        grid = _case14_tight_q()
        ts = _one_row_time_series(grid)
        ts.bus_q_violation_tol_mvar = 1e6
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert len(ts.get_bus_q_violations()[0]) == 0
        assert len(ts.get_bus_q_violations_n()) == 0

    def test_wide_limits_report_nothing(self):
        # the same grid, with the reactive ranges pandapower's case14 actually carries
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(pn.case14())
        ts = _one_row_time_series(grid)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]
        assert len(ts.get_bus_q_violations()[0]) == 0

    def test_dc_is_refused(self):
        grid = _case14_tight_q()
        ts = _one_row_time_series(grid)
        ts.change_algorithm(AlgorithmType.DC_SparseLU)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        with self.assertRaises(RuntimeError):
            ts.compute(v_init, 20, 1e-11)

    def test_independent_of_compute_limit_violations(self):
        # the two flags are separate opt-ins, and a reactive violation never lands in
        # get_violations() (nor a current one in get_bus_q_violations())
        grid = _case14_tight_q()
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_bus_q_violations = True
        ca.bus_q_violation_tol_mvar = 0.
        ca.add_n1(0)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ca.compute(v_init, 20, 1e-11)
        assert ca.compute_limit_violations is False
        with self.assertRaises(RuntimeError):
            ca.get_violations()
        assert len(ca.get_bus_q_violations()) == 1
        assert len(ca.get_bus_q_violations()[0]) > 0
        for v in ca.get_bus_q_violations()[0]:
            assert violation_category(v.violation_type) == ViolationCategory.PHYSICAL

    def test_toggling_the_flag_keeps_the_registered_contingencies(self):
        # unlike compute_limit_violations, whose setter clear()s the whole object
        grid = _case14_tight_q()
        ca = ContingencyAnalysisCPP(grid)
        ca.add_n1(0)
        ca.add_n1(1)
        ca.compute_bus_q_violations = True
        assert len(ca.my_defaults()) == 2
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ca.compute(v_init, 20, 1e-11)
        assert len(ca.get_bus_q_violations()) == 2


class TestSvcCapabilityFromPython(unittest.TestCase):
    """the one family whose capability is not already in MVAr: an SVC's `b_min` / `b_max` are
    a susceptance range, and the reported `limit` is that range at the solved voltage. A user
    must be able to reproduce it from the published data, which is what this pins."""

    def test_limit_is_b_times_v_squared(self):
        from lightsim2grid.lightsim2grid_cpp import LSGrid, SvcContainer
        sn_mva, vn_kv, b_max = 100., 138., 0.2
        grid = LSGrid()
        grid.set_sn_mva(sn_mva)
        grid.set_init_vm_pu(1.0)
        grid.init_bus(2, 1, np.full(2, vn_kv), 0, 0)
        grid.init_powerlines(np.array([0.01]), np.array([0.1]), np.zeros(1, dtype=complex),
                             np.array([0]), np.array([1]))
        grid.init_loads(np.array([40.]), np.array([30.]), np.array([1]))
        grid.init_generators(np.array([0.]), np.array([1.02]), np.array([-1000.]),
                             np.array([1000.]), np.array([0]))
        grid.add_gen_slackbus(0, 1.)
        # a voltage-mode SVC holding the load bus, with a susceptance range too small for it
        grid.init_svcs([int(SvcContainer.RegulationMode.VOLTAGE)], np.array([1.05]), np.array([0.]),
                       np.array([0.]), np.array([-b_max]), np.array([b_max]),
                       np.array([1]), np.array([1]))
        grid.tell_solver_need_reset()

        ts = TimeSeriesCPP(grid)
        ts.compute_bus_q_violations = True
        ts.bus_q_violation_tol_mvar = 0.
        ts.modify_gen_p(np.array([[gen.target_p_mw for gen in grid.get_generators()]]))
        v_init = np.full(grid.total_bus(), 1.0 + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]

        viols = [v for v in ts.get_bus_q_violations()[0] if v.element_id == 1]
        assert len(viols) == 1, "the SVC cannot hold 1.05 pu with 0.2 pu of susceptance"
        viol = viols[0]
        assert viol.violation_type == LimitViolationType.HIGH_Q
        assert viol.category == ViolationCategory.PHYSICAL

        # reproduce the limit from the published data alone: b_max . |V|^2 . sn_mva. The
        # batch works on a private COPY of the grid, so this object has published nothing
        # yet -- solve the same state on it to get the SVC's own results.
        grid.ac_pf(v_init, 20, 1e-11)
        svc = grid.get_svcs()[0]
        v_pu = svc.res_v_kv / vn_kv
        self.assertAlmostEqual(viol.limit, svc.b_max * v_pu ** 2 * sn_mva, places=6)
        # ... and the value against what the SVC itself published
        self.assertAlmostEqual(viol.value, svc.res_q_mvar, places=5)


class TestBusQViolationsWrapper(unittest.TestCase):
    """the python wrappers: the properties they expose, what they invalidate, and the
    `bus_q_violations` field of the `run()` result"""

    def setUp(self):
        import grid2op
        from lightsim2grid import LightSimBackend
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.env.reset(seed=0, options={"time serie id": 0})

    def tearDown(self):
        self.env.close()

    def test_time_serie_properties(self):
        from lightsim2grid import TimeSerie
        ts = TimeSerie(self.env)
        assert ts.compute_bus_q_violations is False
        assert ts.bus_q_violation_tol_mvar == 1e-4
        with self.assertRaises(RuntimeError):
            ts.get_bus_q_violations()
        with self.assertRaises(ValueError):
            ts.compute_bus_q_violations = "yes"
        with self.assertRaises(ValueError):
            ts.bus_q_violation_tol_mvar = "tight"
        with self.assertRaises(RuntimeError):
            ts.bus_q_violation_tol_mvar = -1.  # rejected C++-side

        ts.compute_bus_q_violations = True
        ts.compute_V(scenario_id=0)
        viols = ts.get_bus_q_violations()
        assert len(viols) == ts.computer.get_voltages().shape[0]
        for row in viols:
            for v in row:
                assert v.element_type == ViolationElementType.BUS
                assert v.category == ViolationCategory.PHYSICAL
        ts.get_bus_q_violations_n()
        ts.close()

    def test_injection_sweep_inherits_it(self):
        from lightsim2grid import InjectionSweep
        sweep = InjectionSweep(self.env)
        assert sweep.compute_bus_q_violations is False
        sweep.compute_bus_q_violations = True
        sweep.compute_V(scenario_id=0)
        assert len(sweep.get_bus_q_violations()) == sweep.computer.get_voltages().shape[0]
        sweep.close()

    def test_contingency_analysis_run_carries_the_field(self):
        from lightsim2grid import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_single_contingency(0)
        sa.add_single_contingency(1)

        # off: the field is there and empty, so a caller can read it unconditionally
        res = sa.run()
        assert res.pre_contingency_result.bus_q_violations == []
        for cont in res.post_contingency_results:
            assert cont.bus_q_violations == []

        # on: the contingencies survive the flag change (unlike compute_limit_violations)
        sa.compute_bus_q_violations = True
        assert len(sa.computer.my_defaults()) == 2
        res = sa.run()
        assert len(res.post_contingency_results) == 2
        for cont in res.post_contingency_results:
            for v in cont.bus_q_violations:
                assert v.category == ViolationCategory.PHYSICAL
            # the two lists never mix
            for v in cont.limit_violations:
                assert v.violation_type not in (LimitViolationType.LOW_Q,
                                                LimitViolationType.HIGH_Q)
        sa.close()

    def test_scenario_sweep_run_carries_the_field(self):
        from lightsim2grid import ScenarioSweep
        sweep = ScenarioSweep(self.env)
        sweep.compute_limit_violations = True
        sweep.compute_bus_q_violations = True
        n_line = len(self.env.backend._grid.get_lines())
        mask = np.zeros((2, n_line), dtype=bool)
        mask[1, 0] = True
        sweep.set_contingency_lines(mask)
        res = sweep.run()
        assert len(res.post_contingency_results) == 2
        for cont in res.post_contingency_results:
            for v in cont.bus_q_violations:
                assert v.element_type == ViolationElementType.BUS
                assert v.category == ViolationCategory.PHYSICAL
        assert sweep.get_bus_q_violations() is not None
        sweep.get_bus_q_violations_n()
        sweep.close()


if __name__ == "__main__":
    unittest.main()
