# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Tests for the (opt-in) physical-limit reporting of the batch algorithms:
`compute_physical_violations` / `physical_violation_tol_mva` / `get_physical_violations` /
`get_physical_violations_n`, on the C++ classes and on the python wrappers
(`TimeSerie`, `InjectionSweep`, `ContingencyAnalysis`, `ScenarioSweep`).

The reported VALUES are pinned against a single-shot `ac_pf` here as they are in
`src/tests/test_batch_physical_violations.cpp`; what this file adds is the python layer --
the properties, what they invalidate, and the `physical_violations` field of the
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
    exactly what `get_physical_violations` reports a value for"""
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
    ts.compute_physical_violations = True
    ts.physical_violation_tol_mva = 0.
    gen_p = np.array([[gen.target_p_mw for gen in grid.get_generators()]])
    ts.modify_gen_p(gen_p)
    return ts


class TestPhysicalViolationsCpp(unittest.TestCase):
    """the C++ classes, through their bindings"""

    def test_default_is_off_and_raises(self):
        grid = _case14_tight_q()
        for cls in (TimeSeriesCPP, ContingencyAnalysisCPP, ScenarioSweepCPP):
            algo = cls(grid)
            assert algo.compute_physical_violations is False, cls.__name__
            assert algo.physical_violation_tol_mva == 1e-4, cls.__name__
            with self.assertRaises(RuntimeError):
                algo.get_physical_violations()
            with self.assertRaises(RuntimeError):
                algo.get_physical_violations_n()

    def test_reports_what_ac_pf_publishes_per_bus(self):
        grid = _case14_tight_q()
        expected = _reference_bus_q(grid)
        assert len(expected) > 1, "the fixture is only meaningful with several regulated buses"

        ts = _one_row_time_series(grid)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]

        viols = ts.get_physical_violations()
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
        assert len(ts.get_physical_violations_n()) == len(viols[0])

    def test_a_wide_tolerance_hides_everything(self):
        grid = _case14_tight_q()
        ts = _one_row_time_series(grid)
        ts.physical_violation_tol_mva = 1e6
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert len(ts.get_physical_violations()[0]) == 0
        assert len(ts.get_physical_violations_n()) == 0

    def test_wide_limits_report_nothing(self):
        # the same grid, with the reactive ranges pandapower's case14 actually carries
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(pn.case14())
        ts = _one_row_time_series(grid)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]
        assert len(ts.get_physical_violations()[0]) == 0

    def test_dc_reports_no_reactive_violation_and_does_not_raise(self):
        # a DC powerflow has no reactive power at all, so the reactive half is not
        # applicable rather than missing -- and the hvdc active-power half (nothing here to
        # trigger it on case14) still runs, so compute() must not raise
        grid = _case14_tight_q()
        ts = _one_row_time_series(grid)
        ts.change_algorithm(AlgorithmType.DC_SparseLU)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]
        assert len(ts.get_physical_violations()[0]) == 0

    def test_independent_of_compute_limit_violations(self):
        # the two flags are separate opt-ins, and a reactive violation never lands in
        # get_violations() (nor a current one in get_physical_violations())
        grid = _case14_tight_q()
        ca = ContingencyAnalysisCPP(grid)
        ca.compute_physical_violations = True
        ca.physical_violation_tol_mva = 0.
        ca.add_n1(0)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ca.compute(v_init, 20, 1e-11)
        assert ca.compute_limit_violations is False
        with self.assertRaises(RuntimeError):
            ca.get_violations()
        assert len(ca.get_physical_violations()) == 1
        assert len(ca.get_physical_violations()[0]) > 0
        for v in ca.get_physical_violations()[0]:
            assert violation_category(v.violation_type) == ViolationCategory.PHYSICAL

    def test_toggling_the_flag_keeps_the_registered_contingencies(self):
        # unlike compute_limit_violations, whose setter clear()s the whole object
        grid = _case14_tight_q()
        ca = ContingencyAnalysisCPP(grid)
        ca.add_n1(0)
        ca.add_n1(1)
        ca.compute_physical_violations = True
        assert len(ca.my_defaults()) == 2
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ca.compute(v_init, 20, 1e-11)
        assert len(ca.get_physical_violations()) == 2


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
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.modify_gen_p(np.array([[gen.target_p_mw for gen in grid.get_generators()]]))
        v_init = np.full(grid.total_bus(), 1.0 + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]

        viols = [v for v in ts.get_physical_violations()[0] if v.element_id == 1]
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


class TestStorageCapabilityFromPython(unittest.TestCase):
    """a storage unit regulating its own bus (``init_storages_full``) holds it like a local
    generator: its reactive range counts, and the reported value is the reactive power it
    produced -- GENERATOR convention, whereas ``StorageInfo.res_q_mvar`` is in load
    convention."""
    BUS = 9  # a load bus of case14, no generator on it

    def _grid(self, max_q_mvar):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(pn.case14())
        grid.init_storages_full(np.array([0.]), np.array([0.]), [True], np.array([1.035]),
                                np.array([-max_q_mvar]), np.array([max_q_mvar]),
                                np.array([self.BUS], dtype=np.int32))
        grid.tell_solver_need_reset()
        return grid

    def _bus_viols(self, grid):
        ts = _one_row_time_series(grid)
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        ts.compute(v_init, 20, 1e-11)
        assert ts.converged_mask()[0]
        return [v for v in ts.get_physical_violations()[0] if v.element_id == self.BUS]

    def test_limits_of_a_regulating_storage_unit(self):
        # what the unit has to produce, from a single solve with a range wide enough
        grid = self._grid(1e3)
        assert len(self._bus_viols(grid)) == 0
        v_init = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
        grid.ac_pf(v_init, 20, 1e-11)
        q_prod = -grid.get_storages()[0].res_q_mvar
        assert abs(q_prod) > 1., "the fixture needs the unit to actually produce or absorb"

        half = 0.5 * abs(q_prod)
        viols = self._bus_viols(self._grid(half))
        assert len(viols) == 1, "half of what it produces cannot hold the bus"
        viol = viols[0]
        assert viol.element_type == ViolationElementType.BUS
        assert viol.category == ViolationCategory.PHYSICAL
        if q_prod > 0.:
            assert viol.violation_type == LimitViolationType.HIGH_Q
            self.assertAlmostEqual(viol.limit, half, places=9)
        else:
            assert viol.violation_type == LimitViolationType.LOW_Q
            self.assertAlmostEqual(viol.limit, -half, places=9)
        self.assertAlmostEqual(viol.value, q_prod, places=5)


class TestHvdcPFromPython(unittest.TestCase):
    """the active-power half: an angle-droop hvdc line beyond what its converters can
    transmit. `status_droop` is an input of the solve, so nothing saturates the droop -- see
    `LSGrid.set_status_droop_hvdc`, whose documentation says the saturation belongs between
    two solves, which is the outer loop this detects."""

    @staticmethod
    def _droop_grid(pmax_1to2, pmax_2to1=1000., p0_mw=30., k_mw_per_deg=400.):
        from lightsim2grid.lightsim2grid_cpp import LSGrid
        grid = LSGrid()
        grid.set_sn_mva(100.)
        grid.set_init_vm_pu(1.0)
        grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
        grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                             np.array([0, 1, 2]), np.array([1, 2, 3]))
        grid.init_loads(np.array([80.]), np.array([60.]), np.array([3]))
        grid.init_generators(np.array([0.]), np.array([1.02]), np.array([-1e3]),
                             np.array([1e3]), np.array([0]))
        grid.add_gen_slackbus(0, 1.)
        one = np.ones(1)
        # bus 1 -> bus 3, no converter or dc-line losses, so the flow leaving side 1 is
        # exactly p0 + k.(theta1 - theta3). NB: the slope argument is MW per DEGREE.
        grid.init_hvdc_lines(np.array([1]), np.array([3]), [0], [0], np.zeros(1), np.zeros(1),
                             [False], [False], one, one, np.zeros(1), np.zeros(1),
                             np.full(1, -1e3), np.full(1, 1e3), np.full(1, -1e3), np.full(1, 1e3),
                             one, one, [0], np.zeros(1), np.zeros(1), np.zeros(1),
                             [True], np.array([p0_mw]), np.array([k_mw_per_deg]),
                             np.array([pmax_1to2]), np.array([pmax_2to1]))
        grid.tell_solver_need_reset()
        return grid

    def _one_row(self, grid):
        ts = TimeSeriesCPP(grid)
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.modify_gen_p(np.array([[gen.target_p_mw for gen in grid.get_generators()]]))
        v_init = np.full(grid.total_bus(), 1.0 + 0j)
        ts.compute(v_init, 30, 1e-11)
        assert ts.converged_mask()[0]
        return ts

    def test_reports_the_flow_ac_pf_publishes(self):
        ref = self._droop_grid(pmax_1to2=1000.)
        ref.ac_pf(np.full(ref.total_bus(), 1.0 + 0j), 30, 1e-11)
        hv = ref.get_dclines()[0]
        # generator convention: side 1 DRAWS this much from the AC grid
        p_flow = -hv.res_p1_mw
        assert p_flow > 1., "the fixture must actually flow 1 -> 2"
        self.assertAlmostEqual(p_flow, hv.res_p2_mw, places=6)  # lossless here

        pmax = p_flow - 5.
        grid = self._droop_grid(pmax_1to2=pmax)
        grid.set_dcline_names(["dc_link"])
        viols = self._one_row(grid).get_physical_violations()[0]
        assert len(viols) == 1
        assert viols[0].element_type == ViolationElementType.HVDC
        assert viols[0].element_id == 0
        assert viols[0].side == 1  # the flow leaves side 1
        assert viols[0].violation_type == LimitViolationType.HIGH_P
        assert viols[0].category == ViolationCategory.PHYSICAL
        self.assertAlmostEqual(viols[0].value, p_flow, places=6)
        self.assertAlmostEqual(viols[0].limit, pmax, places=9)
        assert viols[0].name == "dc_link"

    def test_within_the_limit_reports_nothing(self):
        grid = self._droop_grid(pmax_1to2=1000.)
        assert len(self._one_row(grid).get_physical_violations()[0]) == 0

    def test_an_already_saturated_droop_is_not_reported(self):
        # status_droop != 0: the solver pins the flow at the very limit this would compare
        # against, so there is nothing left to detect
        grid = self._droop_grid(pmax_1to2=5.)
        grid.set_status_droop_hvdc(0, 1)
        assert len(self._one_row(grid).get_physical_violations()[0]) == 0

    def test_it_works_in_dc_too(self):
        # the hvdc half needs only the bus angles, which a DC powerflow solves
        ref = self._droop_grid(pmax_1to2=5.)
        ref.change_algorithm(AlgorithmType.DC_SparseLU)
        ref.dc_pf(np.full(ref.total_bus(), 1.0 + 0j), 30, 1e-11)
        p_dc = -ref.get_dclines()[0].res_p1_mw
        assert p_dc > 5.

        grid = self._droop_grid(pmax_1to2=5.)
        ts = TimeSeriesCPP(grid)
        ts.change_algorithm(AlgorithmType.DC_SparseLU)
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.modify_gen_p(np.array([[0.]]))
        ts.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        viols = ts.get_physical_violations()[0]
        assert len(viols) == 1
        assert viols[0].violation_type == LimitViolationType.HIGH_P
        self.assertAlmostEqual(viols[0].value, p_dc, places=6)


class TestGenPFromPython(unittest.TestCase):
    """the distributed slack's own half: a generator whose converged active power -- its
    target plus its share of the imbalance -- left `min_p_mw` / `max_p_mw`. The slack is
    solved inside the Jacobian by participation factors that know nothing about limits, so
    nothing stops it; that is what OpenLoadFlow's `DistributedSlack` outer loop re-shares.

    Per machine, where the reactive check is per bus: the active split is not a convention,
    it is the participation factors the caller chose."""

    @staticmethod
    def _slack_grid(w0=1., w1=1.):
        """the 4-bus radial feeder 0-1-2-3 with the 80 MW / 60 MVAr load on bus 3, and the
        slack SHARED between gen 0 (bus 0, target 0 MW) and gen 1 (bus 1, target 10 MW)."""
        from lightsim2grid.lightsim2grid_cpp import LSGrid
        grid = LSGrid()
        grid.set_sn_mva(100.)
        grid.set_init_vm_pu(1.0)
        grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
        grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                             np.array([0, 1, 2]), np.array([1, 2, 3]))
        grid.init_loads(np.array([80.]), np.array([60.]), np.array([3]))
        grid.init_generators(np.array([0., 10.]), np.array([1.02, 1.05]),
                             np.full(2, -1e3), np.full(2, 1e3), np.array([0, 1]))
        grid.add_gen_slackbus(0, w0)
        grid.add_gen_slackbus(1, w1)
        grid.tell_solver_need_reset()
        return grid

    @staticmethod
    def _reference_gen_p(grid):
        """each generator's converged active power as a single ac_pf publishes it: target
        plus its share of the distributed slack, the very number the check re-derives"""
        grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        return np.array([gen.res_p_mw for gen in grid.get_generators()])

    def _one_row(self, grid):
        ts = TimeSeriesCPP(grid)
        ts.compute_physical_violations = True
        ts.physical_violation_tol_mva = 0.
        ts.modify_gen_p(np.array([[gen.target_p_mw for gen in grid.get_generators()]]))
        ts.compute(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-11)
        assert ts.converged_mask()[0]
        return ts

    def test_reports_the_power_ac_pf_publishes(self):
        p_ref = self._reference_gen_p(self._slack_grid())
        assert p_ref[1] > 10., "gen 1 must actually be pushed above its target"

        pmax = p_ref[1] - 5.
        grid = self._slack_grid()
        grid.set_gen_p_limits(np.array([np.nan, np.nan]), np.array([np.nan, pmax]))
        grid.set_gen_names(["slack_unit", "shared_unit"])
        viols = self._one_row(grid).get_physical_violations()[0]
        assert len(viols) == 1
        assert viols[0].element_type == ViolationElementType.GENERATOR
        assert viols[0].element_id == 1
        assert viols[0].side == 0
        assert viols[0].violation_type == LimitViolationType.HIGH_P
        assert viols[0].category == ViolationCategory.PHYSICAL
        self.assertAlmostEqual(viols[0].value, p_ref[1], places=6)
        self.assertAlmostEqual(viols[0].limit, pmax, places=9)
        assert viols[0].name == "shared_unit"

    def test_below_min_p_is_low_p(self):
        p_ref = self._reference_gen_p(self._slack_grid())
        pmin = p_ref[0] + 5.
        grid = self._slack_grid()
        grid.set_gen_p_limits(np.array([pmin, np.nan]), np.array([np.nan, np.nan]))
        viols = self._one_row(grid).get_physical_violations()[0]
        assert len(viols) == 1
        assert viols[0].element_id == 0
        assert viols[0].violation_type == LimitViolationType.LOW_P
        assert viols[0].category == ViolationCategory.PHYSICAL
        self.assertAlmostEqual(viols[0].value, p_ref[0], places=6)

    def test_within_the_limits_reports_nothing(self):
        p_ref = self._reference_gen_p(self._slack_grid())
        grid = self._slack_grid()
        grid.set_gen_p_limits(np.array([-1e3, -1e3]), p_ref + 5.)
        assert len(self._one_row(grid).get_physical_violations()[0]) == 0

    def test_a_grid_without_limits_reports_nothing(self):
        grid = self._slack_grid()
        for gen in grid.get_generators():
            assert np.isnan(gen.min_p_mw)
            assert np.isnan(gen.max_p_mw)
        assert len(self._one_row(grid).get_physical_violations()[0]) == 0

    def test_each_machine_at_its_own_share(self):
        # the participation factors are what the split follows -- a 1:3 share moves both
        # reported values, each still matching what ac_pf gives that machine
        p_even = self._reference_gen_p(self._slack_grid(1., 1.))
        p_ref = self._reference_gen_p(self._slack_grid(1., 3.))
        assert abs(p_ref[1] - p_even[1]) > 1.

        grid = self._slack_grid(1., 3.)
        grid.set_gen_p_limits(np.full(2, np.nan), p_ref - 1.)
        viols = self._one_row(grid).get_physical_violations()[0]
        assert len(viols) == 2
        by_id = {v.element_id: v for v in viols}
        self.assertAlmostEqual(by_id[0].value, p_ref[0], places=6)
        self.assertAlmostEqual(by_id[1].value, p_ref[1], places=6)

    def test_a_non_participant_is_never_reported(self):
        # it keeps its target exactly, so a violation there is an input error rather than
        # something the solve produced
        grid = self._slack_grid()
        grid.remove_gen_slackbus(1)
        grid.set_gen_p_limits(np.full(2, np.nan), np.array([np.nan, 1.]))
        assert len(self._one_row(grid).get_physical_violations()[0]) == 0

    def test_the_limits_are_optional_and_droppable(self):
        grid = self._slack_grid()
        grid.set_gen_p_limits(np.array([0., 0.]), np.array([1., 2.]))
        assert grid.get_generators()[1].max_p_mw == 2.
        grid.set_gen_p_limits(np.array([]), np.array([]))
        assert np.isnan(grid.get_generators()[1].max_p_mw)

    def test_pandapower_columns_are_threaded_through(self):
        # `min_p_mw` / `max_p_mw` on `net.gen` reach `GenInfo`; the generators
        # `_aux_add_slack` appends for an ext_grid have no pandapower row, hence no limit
        net = pn.case14()
        net.gen["min_p_mw"] = 1. + np.arange(net.gen.shape[0])
        net.gen["max_p_mw"] = 100. + np.arange(net.gen.shape[0])
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(net)
        gens = grid.get_generators()
        nb_pp_gen = net.gen.shape[0]
        for gen_id in range(nb_pp_gen):
            self.assertAlmostEqual(gens[gen_id].min_p_mw, 1. + gen_id, places=9)
            self.assertAlmostEqual(gens[gen_id].max_p_mw, 100. + gen_id, places=9)
        for gen_id in range(nb_pp_gen, len(gens)):
            assert np.isnan(gens[gen_id].min_p_mw)
            assert np.isnan(gens[gen_id].max_p_mw)

    def test_no_pandapower_column_means_no_limit(self):
        net = pn.case14()
        net.gen["min_p_mw"] = np.nan
        net.gen["max_p_mw"] = np.nan
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(net)
        for gen in grid.get_generators():
            assert np.isnan(gen.min_p_mw)
            assert np.isnan(gen.max_p_mw)


class TestPhysicalViolationsWrapper(unittest.TestCase):
    """the python wrappers: the properties they expose, what they invalidate, and the
    `physical_violations` field of the `run()` result"""

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
        assert ts.compute_physical_violations is False
        assert ts.physical_violation_tol_mva == 1e-4
        with self.assertRaises(RuntimeError):
            ts.get_physical_violations()
        with self.assertRaises(ValueError):
            ts.compute_physical_violations = "yes"
        with self.assertRaises(ValueError):
            ts.physical_violation_tol_mva = "tight"
        with self.assertRaises(RuntimeError):
            ts.physical_violation_tol_mva = -1.  # rejected C++-side

        ts.compute_physical_violations = True
        ts.compute_V(scenario_id=0)
        viols = ts.get_physical_violations()
        assert len(viols) == ts.computer.get_voltages().shape[0]
        for row in viols:
            for v in row:
                assert v.element_type == ViolationElementType.BUS
                assert v.category == ViolationCategory.PHYSICAL
        ts.get_physical_violations_n()
        ts.close()

    def test_injection_sweep_inherits_it(self):
        from lightsim2grid import InjectionSweep
        sweep = InjectionSweep(self.env)
        assert sweep.compute_physical_violations is False
        sweep.compute_physical_violations = True
        sweep.compute_V(scenario_id=0)
        assert len(sweep.get_physical_violations()) == sweep.computer.get_voltages().shape[0]
        sweep.close()

    def test_contingency_analysis_run_carries_the_field(self):
        from lightsim2grid import ContingencyAnalysis
        sa = ContingencyAnalysis(self.env, compute_limit_violations=True)
        sa.add_single_contingency(0)
        sa.add_single_contingency(1)

        # off: the field is there and empty, so a caller can read it unconditionally
        res = sa.run()
        assert res.pre_contingency_result.physical_violations == []
        for cont in res.post_contingency_results:
            assert cont.physical_violations == []

        # on: the contingencies survive the flag change (unlike compute_limit_violations)
        sa.compute_physical_violations = True
        assert len(sa.computer.my_defaults()) == 2
        res = sa.run()
        assert len(res.post_contingency_results) == 2
        for cont in res.post_contingency_results:
            for v in cont.physical_violations:
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
        sweep.compute_physical_violations = True
        n_line = len(self.env.backend._grid.get_lines())
        mask = np.zeros((2, n_line), dtype=bool)
        mask[1, 0] = True
        sweep.set_contingency_lines(mask)
        res = sweep.run()
        assert len(res.post_contingency_results) == 2
        for cont in res.post_contingency_results:
            for v in cont.physical_violations:
                assert v.element_type == ViolationElementType.BUS
                assert v.category == ViolationCategory.PHYSICAL
        assert sweep.get_physical_violations() is not None
        sweep.get_physical_violations_n()
        sweep.close()


if __name__ == "__main__":
    unittest.main()
