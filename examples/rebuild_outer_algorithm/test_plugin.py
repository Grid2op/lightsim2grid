#!/usr/bin/env python3
"""Tests of the NROuterRebuild_* plugin: OpenLoadFlow's outer loops on a Newton whose Jacobian
is rebuilt, and analyzed again, whenever the active equations change (see
RebuildOuterInner.hpp).

Build the plugin first (see RebuildOuter_plugin.cpp), then run:
    python examples/rebuild_outer_algorithm/test_plugin.py

NROuter_* solves the same equations on one fixed pattern, so both must reach the same
voltages, after the same outer iterations; only what the linear solver does differs.
"""
import platform
import pathlib
import sys
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import AlgorithmType, ContingencyAnalysisCPP
from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import (OuterLoopStatus, DistributedSlack, ReactiveLimits, TransformerVoltageControl,
                                     ShuntVoltageControl, PhaseControl, HvdcAcEmulationLimits)

try:
    import pandas as pd
    import pypowsybl as pp
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

MAX_IT = 30
TOL = 1e-10
ATOL_V = 1e-8


def find_plugin():
    build = pathlib.Path(__file__).parent / "build"
    if platform.system() == "Windows":
        candidates = [build / "Release" / "rebuild_outer_algorithm.dll", build / "rebuild_outer_algorithm.dll"]
    else:
        candidates = [build / "librebuild_outer_algorithm.so"]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(f"Plugin not found (tried {[str(c) for c in candidates]}). "
                            "Build it first (see RebuildOuter_plugin.cpp).")


lightsim2grid.load_algorithm_plugin(find_plugin())
_names = lightsim2grid.network.LSGrid().available_algorithm_names()
if "NROuterRebuild_KLU" in _names:
    REBUILD, FIXED = "NROuterRebuild_KLU", "NROuter_KLU"
else:
    REBUILD, FIXED = "NROuterRebuild_SparseLU", "NROuter_SparseLU"


def _pandapower(net, algo):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(net)
    grid.change_algorithm(algo)
    return grid


def _case14(algo, distributed=True):
    """case14, its generators short of reactive power (ReactiveLimits switches buses PQ) and,
    if `distributed`, sharing the slack in proportion to their maximum output"""
    net = pn.case14()
    net.gen["min_q_mvar"] = -15.
    net.gen["max_q_mvar"] = 15.
    grid = _pandapower(net, algo)
    nb_gen = len(grid.get_generators())
    grid.set_gen_p_limits(np.zeros(nb_gen), np.full(nb_gen, 300.))
    if distributed:
        grid.set_gen_can_participate_slack([True] * nb_gen, np.full(nb_gen, 300.))
    return grid


def _loops(grid, *loops):
    grid.clear_outer_loops()
    for loop in loops:
        grid.add_outer_loop(loop)
    return grid


def _solve(grid, v0=1.04):
    return grid.ac_pf(np.full(grid.total_bus(), v0 + 0j), MAX_IT, TOL)


def _stats(grid):
    return grid.get_algo().get_outer_loop_stats(), grid.get_algo().get_linear_solver_stats()


class TestRegistration(unittest.TestCase):
    def test_registered_as_an_outer_loop_algorithm(self):
        for name in ("NROuterRebuild_SparseLU", "NROuterRebuild_KLU"):
            if name not in _names:
                continue
            grid = _case14(name)
            self.assertEqual(grid.get_algo_type(), AlgorithmType.Custom)
            self.assertTrue(grid.get_algo().supports_outer_loops())


class _SameAsFixedPattern:
    """the rebuilt Newton reaches what the fixed-pattern one does, after the same outer
    iterations; `_make(algo)` builds the grid, its loops included"""

    def _make(self, algo):
        raise NotImplementedError

    def _check(self, v0=1.04):
        fixed, rebuilt = self._make(FIXED), self._make(REBUILD)
        V_fixed, V_rebuilt = _solve(fixed, v0), _solve(rebuilt, v0)
        self.assertGreater(V_fixed.shape[0], 0, "the fixed-pattern solve failed")
        self.assertGreater(V_rebuilt.shape[0], 0, f"the rebuilt solve failed: {rebuilt.get_algo().get_error()}")
        np.testing.assert_allclose(V_rebuilt, V_fixed, rtol=0., atol=ATOL_V)
        np.testing.assert_allclose(np.stack(rebuilt.get_trafo_res1()[:2]), np.stack(fixed.get_trafo_res1()[:2]),
                                   rtol=0., atol=1e-5)
        (out_f, lin_f), (out_r, lin_r) = _stats(fixed), _stats(rebuilt)
        self.assertEqual(out_r.status, OuterLoopStatus.STABLE)
        self.assertEqual(out_r.nb_outer_iterations, out_f.nb_outer_iterations)
        self.assertEqual(lin_f.nb_analyze, 1)
        return out_r, lin_r


class TestCase118(_SameAsFixedPattern, unittest.TestCase):
    def _make(self, algo):
        return _loops(_pandapower(pn.case118(), algo), DistributedSlack(), ReactiveLimits())

    def test_same_solution(self):
        out, lin = self._check()
        self.assertGreater(out.nb_outer_iterations, 0)
        # ReactiveLimits switched buses: a new pattern, analyzed again
        self.assertGreater(lin.nb_analyze, 1)


class TestCase14(_SameAsFixedPattern, unittest.TestCase):
    def _make(self, algo):
        return _loops(_case14(algo), DistributedSlack(), ReactiveLimits())

    def test_same_solution(self):
        self._check()


class TestCostTiers(unittest.TestCase):
    """what the linear solver does is what OpenLoadFlow's JacobianMatrix does"""

    def test_values_only_is_a_refactorization(self):
        # DistributedSlack only moves targets: one analysis, then refactorizations
        grid = _loops(_case14(REBUILD), DistributedSlack())
        self.assertGreater(_solve(grid).shape[0], 0)
        out, lin = _stats(grid)
        self.assertGreater(out.nb_outer_iterations, 0)
        self.assertEqual(lin.nb_analyze, 1)
        self.assertEqual(lin.nb_factorize - lin.nb_fallback_factorize, 1)
        self.assertGreater(lin.nb_refactorize, 0)

    def test_each_switch_is_one_analysis(self):
        # ReactiveLimits only: every re-solve follows a PV / PQ switch
        grid = _loops(_case14(REBUILD, distributed=False), ReactiveLimits())
        self.assertGreater(_solve(grid).shape[0], 0)
        out, lin = _stats(grid)
        self.assertGreater(out.nb_outer_iterations, 0)
        self.assertEqual(lin.nb_analyze, 1 + out.nb_outer_iterations)

    def test_a_new_solve_starts_from_the_grid_structure(self):
        # the second solve starts from the grid again, not from where the first one left the
        # buses: one analysis for that, then one per switch
        grid = _loops(_case14(REBUILD, distributed=False), ReactiveLimits())
        _solve(grid)
        first = _stats(grid)[1].nb_analyze
        _solve(grid)
        out, lin = _stats(grid)
        self.assertEqual(lin.nb_analyze - first, out.nb_outer_iterations + 1)


    def test_hvdc_regime_is_a_factorization(self):
        # the hvdc line saturates: its terms switch, on the pattern analyzed already (the
        # helper is the source tree's, next to this example)
        sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "lightsim2grid" / "tests"))
        try:
            from _aux_make_hvdc import make_case14_hvdc
        except ImportError:
            self.skipTest("not run from a source tree")

        def make(algo):
            _, grid = make_case14_hvdc(3, 9, loss_factor1=0.011, loss_factor2=0.011, droop_enabled=True, p0=10.,
                                       droop_mw_per_deg=5., pmax12=15., pmax21=300.)
            grid.change_algorithm(algo)
            return _loops(grid, HvdcAcEmulationLimits())
        fixed, rebuilt = make(FIXED), make(REBUILD)
        np.testing.assert_allclose(_solve(rebuilt), _solve(fixed), rtol=0., atol=ATOL_V)
        out, lin = _stats(rebuilt)
        self.assertEqual(out.nb_outer_iterations, 1)
        self.assertEqual(lin.nb_analyze, 1)
        self.assertEqual(lin.nb_factorize - lin.nb_fallback_factorize, 2)


class TestContingencyAnalysis(unittest.TestCase):
    """a batch: every row reaches what the fixed-pattern batch does; opening a branch is a
    numeric factorization, a switch a new analysis"""

    N_ROWS = 8

    def _batch(self, algo, *loops):
        grid = _loops(_case14(algo), *(loops or (DistributedSlack(), ReactiveLimits())))
        batch = ContingencyAnalysisCPP(grid)
        for line_id in range(self.N_ROWS):
            batch.add_n1(line_id)
        batch.compute(np.full(grid.total_bus(), 1.04 + 0j), MAX_IT, TOL)
        return batch

    def test_rows_equal_the_fixed_pattern_ones(self):
        fixed, rebuilt = self._batch(FIXED), self._batch(REBUILD)
        V_fixed, V_rebuilt = fixed.get_voltages(), rebuilt.get_voltages()
        self.assertEqual(V_rebuilt.shape, V_fixed.shape)
        np.testing.assert_allclose(V_rebuilt, V_fixed, rtol=0., atol=ATOL_V)
        res_f, res_r = fixed.get_outer_results(), rebuilt.get_outer_results()
        for i, (f, r) in enumerate(zip(res_f, res_r)):
            with self.subTest(row=i):
                self.assertEqual(r.converged, f.converged)
                self.assertEqual(r.nb_outer_iterations, f.nb_outer_iterations)
        lin_f, lin_r = fixed.get_linear_solver_stats(), rebuilt.get_linear_solver_stats()
        self.assertEqual(lin_f.nb_analyze, 1)
        self.assertGreater(lin_r.nb_analyze, 1)

    def test_opening_a_branch_is_a_factorization(self):
        # DistributedSlack switches nothing: one analysis for the batch (its base solve), then
        # each row factorizes again on it, since its branch is out
        batch = self._batch(REBUILD, DistributedSlack())
        self.assertTrue(all(r.converged for r in batch.get_outer_results()))
        lin = batch.get_linear_solver_stats()
        self.assertEqual(lin.nb_analyze, 1)
        self.assertEqual(lin.nb_factorize - lin.nb_fallback_factorize, 1 + self.N_ROWS)


def _ieee14_tap():
    """pypowsybl's IEEE 14 with a ratio tap changer on T4-9-1 regulating its low voltage side"""
    net = pp.network.create_ieee14()
    rtc = pd.DataFrame.from_records(
        index="id", columns=["id", "target_deadband", "target_v", "oltc", "low_tap", "tap", "regulating", "regulated_side"],
        data=[("T4-9-1", 0.1, 12.6, True, 0, 10, False, "TWO")])
    steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho"],
                                      data=[("T4-9-1", 0., 0., 0., 0., 0.9 + 0.01 * k) for k in range(21)])
    net.create_ratio_tap_changers(rtc, steps)
    net.update_ratio_tap_changers(id="T4-9-1", regulating=True)
    return net


def _ieee14_shunts(bus="B9", target_kv=14.):
    """pypowsybl's IEEE 14 with two linear shunts regulating `bus`"""
    net = pp.network.create_ieee14()
    vl = "VL" + bus[1:]
    sh = pd.DataFrame.from_records(
        index="id", columns=["id", "model_type", "section_count", "target_v", "target_deadband", "voltage_level_id",
                             "bus_id", "connectable_bus_id"],
        data=[("SH-A", "LINEAR", 3, target_kv, 0., vl, bus, bus),
              ("SH-B", "LINEAR", 1, target_kv, 0., vl, bus, bus)])
    model = pd.DataFrame.from_records(index="id", columns=["id", "g_per_section", "b_per_section", "max_section_count"],
                                      data=[("SH-A", 0., 0.02, 10), ("SH-B", 0., 0.05, 4)])
    net.create_shunt_compensators(sh, model)
    net.update_shunt_compensators(id=["SH-A", "SH-B"], voltage_regulation_on=[True, True])
    return net


def _ieee14_pst(target_mw=30.):
    """pypowsybl's IEEE 14 with a phase tap changer on T5-6-1 limiting its active power flow"""
    pst = "T5-6-1"
    net = pp.network.create_ieee14()
    ptc = pd.DataFrame.from_records(
        index="id", columns=["id", "target_deadband", "regulation_mode", "low_tap", "tap", "regulating", "regulated_side"],
        data=[(pst, 0., "ACTIVE_POWER_CONTROL", -10, 0, False, "ONE")])
    steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho", "alpha"],
                                      data=[(pst, 0., 0., 0., 0.5 * a, 1., float(a)) for a in range(-10, 11)])
    net.create_phase_tap_changers(ptc, steps)
    net.update_phase_tap_changers(id=pst, regulation_value=target_mw)
    net.update_phase_tap_changers(id=pst, regulating=True)
    return net


def _ieee14_two_stars(target_kv=12.6):
    """pypowsybl's IEEE 14 with two three-winding transformers from B4, each modelled as a star:
    a fictitious bus STARk and three two-winding legs, HVk from B4, LVLk to a load bus of its
    own and MVk to B9. Both MV legs have a ratio tap changer regulating B9."""
    net = pp.network.create_ieee14()
    sub = net.get_voltage_levels().loc["VL4", "substation_id"]
    vls, buses, loads, legs = [], [], [], []
    for k in (1, 2):
        vls += [(f"STAR{k}", sub, "BUS_BREAKER", 132.), (f"LV{k}", sub, "BUS_BREAKER", 20.)]
        buses += [(f"STAR{k}_0", f"STAR{k}"), (f"LV{k}_0", f"LV{k}")]
        loads += [(f"LD{k}", f"LV{k}", f"LV{k}_0", 10., 2.5)]
        legs += [(f"HV{k}", "VL4", "B4", f"STAR{k}", f"STAR{k}_0", 132., 132., 0.9, 17.),
                 (f"LVL{k}", f"STAR{k}", f"STAR{k}_0", f"LV{k}", f"LV{k}_0", 132., 20., 0.02, 0.4),
                 (f"MV{k}", f"STAR{k}", f"STAR{k}_0", "VL9", "B9", 132., 12., 0.008, 0.15)]
    net.create_voltage_levels(pd.DataFrame.from_records(
        vls, index="id", columns=["id", "substation_id", "topology_kind", "nominal_v"]))
    net.create_buses(pd.DataFrame.from_records(buses, index="id", columns=["id", "voltage_level_id"]))
    net.create_loads(pd.DataFrame.from_records(loads, index="id", columns=["id", "voltage_level_id", "bus_id", "p0", "q0"]))
    net.create_2_windings_transformers(pd.DataFrame.from_records(
        legs, index="id", columns=["id", "voltage_level1_id", "bus1_id", "voltage_level2_id", "bus2_id",
                                   "rated_u1", "rated_u2", "r", "x"]).assign(g=0., b=0.))
    rtc = pd.DataFrame.from_records(
        index="id", columns=["id", "target_deadband", "target_v", "oltc", "low_tap", "tap", "regulating", "regulated_side"],
        data=[("MV1", 0., target_kv, True, 0, 8, False, "TWO"), ("MV2", 0., target_kv, True, 0, 7, False, "TWO")])
    steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho"],
                                      data=[(f"MV{k}", 0., 0., 0., 0., 0.88 + 0.01 * s) for k in (1, 2) for s in range(25)])
    net.create_ratio_tap_changers(rtc, steps)
    net.update_ratio_tap_changers(id=["MV1", "MV2"], regulating=[True, True])
    return net


def _pypowsybl(net, algo, reactive_limits=False):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                   olf_rules=OlfLoadingParameters(reactive_limits=reactive_limits))
    grid.change_algorithm(algo)
    return grid


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestTransformerVoltageControl(_SameAsFixedPattern, unittest.TestCase):
    def _make(self, algo):
        return _loops(_pypowsybl(_ieee14_tap(), algo, reactive_limits=True), ReactiveLimits(),
                      TransformerVoltageControl(use_initial_tap_position=True, max_controlled_nominal_voltage=120.))

    def test_same_solution(self):
        self._check(1.0)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestShuntVoltageControl(_SameAsFixedPattern, unittest.TestCase):
    def _make(self, algo):
        return _loops(_pypowsybl(_ieee14_shunts(), algo), ShuntVoltageControl())

    def test_same_solution(self):
        self._check(1.0)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestPhaseControl(_SameAsFixedPattern, unittest.TestCase):
    def _make(self, algo):
        return _loops(_pypowsybl(_ieee14_pst(), algo), PhaseControl())

    def test_same_solution(self):
        self._check(1.0)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestRowAfterATrippedRatioControl(unittest.TestCase):
    """a batch row reaches what the fixed-pattern batch does whatever the rows before it: here a
    row islanding a load bus, then one tripping a regulating MV leg, then one where both MV legs
    regulate B9 again"""

    ROWS = ("LVL1", "MV1", "HV2")

    def _batch(self, algo):
        grid = _loops(_pypowsybl(_ieee14_two_stars(), algo),
                      TransformerVoltageControl(use_initial_tap_position=True, max_controlled_nominal_voltage=120.))
        branches = [el.name for el in grid.get_lines()] + [el.name for el in grid.get_trafos()]
        rows = [branches.index(name) for name in self.ROWS]
        batch = ContingencyAnalysisCPP(grid)
        batch.handle_disconnected_grid = True  # the LVL1 row islands LV1
        batch.add_multiple_n1(rows)
        self.assertEqual([int(c[0]) for c in batch.my_defaults()], rows)
        batch.compute(np.full(grid.total_bus(), 1. + 0j), MAX_IT, TOL)
        return batch

    def test_rows_equal_the_fixed_pattern_ones(self):
        fixed, rebuilt = self._batch(FIXED), self._batch(REBUILD)
        V_fixed, V_rebuilt = fixed.get_voltages(), rebuilt.get_voltages()
        res_f, res_r = fixed.get_outer_results(), rebuilt.get_outer_results()
        for i, name in enumerate(self.ROWS):
            with self.subTest(row=name):
                self.assertTrue(res_f[i].converged)
                self.assertTrue(res_r[i].converged)
                self.assertEqual(sorted(res_r[i].ratio_tap), sorted(res_f[i].ratio_tap))
                np.testing.assert_allclose(V_rebuilt[i], V_fixed[i], rtol=0., atol=ATOL_V)


if __name__ == "__main__":
    unittest.main()
