# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The non-structural setters of LSGrid (issue #299): every value `init_*` takes and
that does not change the set of elements can be changed on a kept grid.

Each test checks the two things a caller that updates a grid in place relies on:

* the setter means what the `init_*` argument means: the kept grid, changed and solved
  again from its previous solution, lands where a grid built from scratch with that
  value lands;
* the setter raises the flags its change needs: the kept grid lands where a copy of it
  (whose solver starts from nothing) lands from the same starting point. A forgotten
  flag shows up here as an answer read off a stale cache.

The grid is the exotic-elements case (IEEE 14 + an SVC, a storage unit, three hvdc
lines -- an LCC one, an angle-droop VSC one, a plain VSC one -- and a phase shifter).
"""

import unittest
from unittest import mock

import numpy as np

from lightsim2grid.network import LSGrid
from lightsim2grid.tests import _exotic_elements_case_ls as case

MAX_IT = 30
TOL = 1e-11
ATOL_SAME_GRID = 1e-8   # two solves of the same grid from two starting points
ATOL_SAME_START = 1e-9  # the same grid, the same starting point, a cold vs a kept solver

# argument positions of the python-facing init functions, as the fixture calls them
GEN_VREG, GEN_MIN_Q, GEN_MAX_Q = 3, 4, 5
SVC_MODE, SVC_VM, SVC_Q, SVC_SLOPE, SVC_BMIN, SVC_BMAX = 0, 1, 2, 3, 4, 5
(HVDC_LF1, HVDC_LF2, HVDC_VREG1, HVDC_Q1, HVDC_MINQ1, HVDC_MAXQ1, HVDC_PF1, HVDC_PF2,
 HVDC_MODE, HVDC_P, HVDC_R, HVDC_VNOM, HVDC_DROOP, HVDC_P0, HVDC_K, HVDC_PMAX12, HVDC_PMAX21) = (
    4, 5, 6, 10, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)

LCC_LINE, DROOP_LINE, VSC_LINE = 0, 1, 2


def _set(args, pos, el_id, value):
    """args[pos][el_id] = value, on a copy of that argument (list or array)"""
    arg = args[pos]
    arg = list(arg) if isinstance(arg, list) else np.array(arg, copy=True)
    arg[el_id] = value
    args[pos] = arg
    return args


def _build(**patches):
    """The fixture grid, with some `init_*` calls patched: `patches[name](self, orig, args)`
    makes the call in place of `orig(self, *args)`."""
    if not patches:
        return case.build_exotic_elements_case_grid()
    methods = {}
    for name, fun in patches.items():
        orig = getattr(LSGrid, name)

        def meth(self, *args, _orig=orig, _fun=fun):
            return _fun(self, _orig, list(args))
        methods[name] = meth
    patched = type("_PatchedLSGrid", (LSGrid,), methods)
    with mock.patch.object(case, "LSGrid", patched):
        return case.build_exotic_elements_case_grid()


def _args_patch(fun):
    """a patch that only rewrites the arguments"""
    return lambda self, orig, args: orig(self, *fun(args))


def _solve(grid, V0=None):
    if V0 is None:
        V0 = np.full(grid.total_bus(), grid.get_init_vm_pu() + 0j)
    V = grid.ac_pf(np.array(V0, copy=True), MAX_IT, TOL)
    assert V.shape[0] > 0, "the powerflow did not converge"
    return V


class TestSettersInPlace(unittest.TestCase):

    def _check(self, change, reference):
        """`change(grid)` applied to a solved grid, compared to `reference` (a grid built
        with the value at init, or None to skip that half) and to a copy of the changed
        grid solved from the same point."""
        kept = _build()
        V_before = _solve(kept)
        change(kept)
        cold = kept.copy()
        V_kept = _solve(kept, V_before)
        V_cold = _solve(cold, V_before)
        np.testing.assert_allclose(V_kept, V_cold, rtol=0, atol=ATOL_SAME_START,
                                   err_msg="the kept grid read a stale cache (a flag is missing)")
        # a stale Jacobian (one missing an entry) still converges to the right point, in
        # more iterations: the iteration counts tell it apart
        self.assertEqual(kept.get_algo().get_nb_iter(), cold.get_algo().get_nb_iter(),
                         "the kept grid solved with a stale Jacobian (a flag is missing)")
        if reference is not None:
            V_ref = _solve(reference)
            np.testing.assert_allclose(V_kept, V_ref, rtol=0, atol=ATOL_SAME_GRID,
                                       err_msg="the setter does not mean what the init argument means")
        return kept, V_before, V_kept

    # ---- generators ----------------------------------------------------------------
    def test_gen_voltage_regulator_off(self):
        ref = _build(init_generators_full=_args_patch(lambda a: _set(a, GEN_VREG, 1, False)))
        kept, V_before, V_kept = self._check(lambda g: g.set_gen_voltage_regulator_on(1, False), ref)
        self.assertFalse(kept.get_generators()[1].voltage_regulator_on)
        # it now injects its reactive setpoint: its bus no longer sits at its target
        self.assertGreater(abs(abs(V_kept[6]) - abs(V_before[6])), 1e-6)
        # and back on
        kept.set_gen_voltage_regulator_on(1, True)
        np.testing.assert_allclose(_solve(kept, V_kept), V_before, rtol=0, atol=ATOL_SAME_GRID)

    def test_gen_q_limits(self):
        ref = _build(init_generators_full=_args_patch(
            lambda a: _set(_set(a, GEN_MIN_Q, 1, -20.), GEN_MAX_Q, 1, 30.)))
        kept, _, _ = self._check(lambda g: g.set_gen_q_limits(1, -20., 30.), ref)
        gen = kept.get_generators()[1]
        self.assertEqual((gen.min_q_mvar, gen.max_q_mvar), (-20., 30.))
        with self.assertRaises(RuntimeError):
            kept.set_gen_q_limits(1, 10., -10.)
        with self.assertRaises(IndexError):
            kept.set_gen_q_limits(99, -1., 1.)

    # ---- storage units -------------------------------------------------------------
    def test_storage_voltage_regulation(self):
        def full(self, orig, args):
            p, q, bus = args
            return self.init_storages_full(p, q, [True], np.array([1.02]), np.array([-30.]),
                                           np.array([30.]), bus)
        ref = _build(init_storages=full)

        def change(grid):
            grid.change_v_storage(0, 1.02)
            grid.set_storage_q_limits(0, -30., 30.)
            grid.set_storage_voltage_regulator_on(0, True)
        kept, _, V_kept = self._check(change, ref)
        sto = kept.get_storages()[0]
        self.assertTrue(sto.voltage_regulator_on)
        self.assertEqual((sto.min_q_mvar, sto.max_q_mvar), (-30., 30.))
        self.assertAlmostEqual(abs(V_kept[sto.bus_id]), 1.02, places=8)

    # ---- SVC -----------------------------------------------------------------------
    def test_svc_target_v(self):
        ref = _build(init_svcs=_args_patch(lambda a: _set(a, SVC_VM, 0, 1.03)))
        kept, _, _ = self._check(lambda g: g.change_v_svc(0, 1.03), ref)
        self.assertEqual(kept.get_svcs()[0].target_vm_pu, 1.03)

    def test_svc_reactive_power_mode(self):
        ref = _build(init_svcs=_args_patch(lambda a: _set(_set(a, SVC_MODE, 0, 2), SVC_Q, 0, 7.5)))

        def change(grid):
            grid.change_q_svc(0, 7.5)
            grid.set_svc_regulation_mode(0, 2)
        kept, _, _ = self._check(change, ref)
        svc = kept.get_svcs()[0]
        self.assertEqual(svc.regulation_mode, 2)
        self.assertAlmostEqual(svc.res_q_mvar, 7.5)

    def test_svc_off_and_back(self):
        ref = _build(init_svcs=_args_patch(lambda a: _set(a, SVC_MODE, 0, 0)))
        kept, V_before, V_kept = self._check(lambda g: g.set_svc_regulation_mode(0, 0), ref)
        kept.set_svc_regulation_mode(0, 1)
        np.testing.assert_allclose(_solve(kept, V_kept), V_before, rtol=0, atol=ATOL_SAME_GRID)
        with self.assertRaises(RuntimeError):
            kept.set_svc_regulation_mode(0, 3)

    def test_svc_slope_and_b_limits(self):
        ref = _build(init_svcs=_args_patch(
            lambda a: _set(_set(_set(a, SVC_SLOPE, 0, 0.02), SVC_BMIN, 0, -0.05), SVC_BMAX, 0, 0.04)))

        def change(grid):
            grid.set_svc_slope(0, 0.02)
            grid.set_svc_b_limits(0, -0.05, 0.04)
        kept, _, _ = self._check(change, ref)
        svc = kept.get_svcs()[0]
        self.assertEqual((svc.slope_pu, svc.b_min, svc.b_max), (0.02, -0.05, 0.04))
        with self.assertRaises(RuntimeError):
            kept.set_svc_b_limits(0, 1., -1.)

    # ---- hvdc: converter stations --------------------------------------------------
    def test_station_reactive_mode(self):
        def patch(a):
            a = _set(a, HVDC_VREG1, VSC_LINE, False)
            return _set(a, HVDC_Q1, VSC_LINE, 12.)
        ref = _build(init_hvdc_lines=_args_patch(patch))

        def change(grid):
            grid.change_hvdc_station_q(VSC_LINE, 1, 12.)
            grid.set_hvdc_station_voltage_regulator_on(VSC_LINE, 1, False)
        kept, _, _ = self._check(change, ref)
        station = kept.get_dclines()[VSC_LINE].station1
        self.assertFalse(station.voltage_regulator_on)
        self.assertEqual(station.target_q_mvar, 12.)

    def test_station_q_limits(self):
        ref = _build(init_hvdc_lines=_args_patch(
            lambda a: _set(_set(a, HVDC_MINQ1, VSC_LINE, -40.), HVDC_MAXQ1, VSC_LINE, 45.)))
        kept, _, _ = self._check(lambda g: g.set_hvdc_station_q_limits(VSC_LINE, 1, -40., 45.), ref)
        station = kept.get_dclines()[VSC_LINE].station1
        self.assertEqual((station.min_q_mvar, station.max_q_mvar), (-40., 45.))

    def test_station_errors(self):
        grid = _build()
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_station_voltage_regulator_on(LCC_LINE, 1, True)   # an LCC
        with self.assertRaises(RuntimeError):
            grid.change_hvdc_station_q(LCC_LINE, 2, 1.)
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_station_q_limits(VSC_LINE, 3, -1., 1.)             # no side 3
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_loss_factors(VSC_LINE, 0.01, 1.5)
        # the first side was not left half-written
        self.assertEqual(grid.get_dclines()[VSC_LINE].station1.loss_factor,
                         _build().get_dclines()[VSC_LINE].station1.loss_factor)

    def test_loss_factors(self):
        ref = _build(init_hvdc_lines=_args_patch(
            lambda a: _set(_set(a, HVDC_LF1, VSC_LINE, 0.02), HVDC_LF2, VSC_LINE, 0.03)))
        kept, _, _ = self._check(lambda g: g.set_hvdc_loss_factors(VSC_LINE, 0.02, 0.03), ref)
        self._assert_same_station_targets(kept, ref, VSC_LINE)

    def test_lcc_power_factors(self):
        ref = _build(init_hvdc_lines=_args_patch(
            lambda a: _set(_set(a, HVDC_PF1, LCC_LINE, 0.8), HVDC_PF2, LCC_LINE, 0.85)))
        kept, _, _ = self._check(lambda g: g.set_hvdc_power_factors(LCC_LINE, 0.8, 0.85), ref)
        self._assert_same_station_targets(kept, ref, LCC_LINE)

    # ---- hvdc: the line ------------------------------------------------------------
    def _assert_same_station_targets(self, grid, ref, hvdc_id):
        """bit for bit: the setter derives the station targets as init does"""
        mine, theirs = grid.get_dclines()[hvdc_id], ref.get_dclines()[hvdc_id]
        for side in ("station1", "station2"):
            self.assertEqual(getattr(mine, side).target_p_mw, getattr(theirs, side).target_p_mw)
            self.assertEqual(getattr(mine, side).target_q_mvar, getattr(theirs, side).target_q_mvar)

    def test_active_power_setpoint_iidm(self):
        for mode, p in ((1, 5.0), (0, 0.0), (1, 0.0), (0, 3.25)):
            with self.subTest(mode=mode, p=p):
                ref = _build(init_hvdc_lines=_args_patch(
                    lambda a: _set(_set(a, HVDC_MODE, VSC_LINE, mode), HVDC_P, VSC_LINE, p)))
                kept, _, _ = self._check(lambda g: g.set_hvdc_active_power_setpoint(VSC_LINE, mode, p), ref)
                line = kept.get_dclines()[VSC_LINE]
                # kept exactly, 0 MW and the converters mode included
                self.assertEqual((line.converters_mode, line.p_setpoint_mw), (mode, p))
                self._assert_same_station_targets(kept, ref, VSC_LINE)
        grid = _build()
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_active_power_setpoint(VSC_LINE, 2, 1.)
        with self.assertRaises(RuntimeError):
            grid.set_hvdc_active_power_setpoint(VSC_LINE, 0, -1.)

    def test_dc_line_parameters(self):
        ref = _build(init_hvdc_lines=_args_patch(
            lambda a: _set(_set(a, HVDC_R, VSC_LINE, 3.), HVDC_VNOM, VSC_LINE, 10.)))
        kept, _, _ = self._check(lambda g: g.set_hvdc_dc_line_parameters(VSC_LINE, 3., 10.), ref)
        self.assertEqual((kept.get_dclines()[VSC_LINE].r_ohm, kept.get_dclines()[VSC_LINE].nominal_v_kv), (3., 10.))
        self._assert_same_station_targets(kept, ref, VSC_LINE)

    def test_droop_parameters(self):
        ref = _build(init_hvdc_lines=_args_patch(
            lambda a: _set(_set(a, HVDC_P0, DROOP_LINE, 4.), HVDC_K, DROOP_LINE, 2.5)))
        kept, V_before, V_kept = self._check(lambda g: g.set_hvdc_droop_parameters(DROOP_LINE, 4., 2.5), ref)
        line = kept.get_dclines()[DROOP_LINE]
        self.assertEqual(line.droop_p0_mw, 4.)
        self.assertGreater(np.max(np.abs(V_kept - V_before)), 1e-6)   # the change is seen

    def test_droop_pmax(self):
        def change(grid):
            grid.set_hvdc_droop_pmax(DROOP_LINE, 0.5, 0.5)
            grid.set_status_droop_hvdc(DROOP_LINE, 1)     # saturated 1 -> 2: injects its pmax
        kept, V_before, V_kept = self._check(change, None)
        self.assertEqual(kept.get_dclines()[DROOP_LINE].pmax_1to2_mw, 0.5)
        self.assertGreater(np.max(np.abs(V_kept - V_before)), 1e-6)

    def test_droop_disable_and_enable(self):
        ref = _build(init_hvdc_lines=_args_patch(lambda a: _set(a, HVDC_DROOP, DROOP_LINE, False)))
        kept, V_before, V_kept = self._check(lambda g: g.set_hvdc_droop_enabled(DROOP_LINE, False), ref)
        self.assertFalse(kept.get_dclines()[DROOP_LINE].droop_enabled)
        kept.set_hvdc_droop_enabled(DROOP_LINE, True)
        np.testing.assert_allclose(_solve(kept, V_kept), V_before, rtol=0, atol=ATOL_SAME_GRID)

    def test_droop_enable_on_another_line(self):
        def patch(a):
            a = _set(a, HVDC_DROOP, VSC_LINE, True)
            a = _set(a, HVDC_P0, VSC_LINE, 1.5)
            return _set(a, HVDC_K, VSC_LINE, 2.)
        ref = _build(init_hvdc_lines=_args_patch(patch))

        def change(grid):
            grid.set_hvdc_droop_parameters(VSC_LINE, 1.5, 2.)
            grid.set_hvdc_droop_enabled(VSC_LINE, True)
        self._check(change, ref)
        with self.assertRaises(RuntimeError):
            _build().set_hvdc_droop_enabled(LCC_LINE, True)       # an LCC line

    # ---- branches ------------------------------------------------------------------
    def test_one_powerline_matches_whole_vectors(self):
        line_id, r, x, h1, h2 = 3, 0.25, 0.21, 0.01j, 0.02j

        def whole(grid):
            lines = grid.get_lines()
            rs = np.array([l.r_pu for l in lines]); xs = np.array([l.x_pu for l in lines])
            h1s = np.array([l.h1_pu for l in lines]); h2s = np.array([l.h2_pu for l in lines])
            rs[line_id], xs[line_id], h1s[line_id], h2s[line_id] = r, x, h1, h2
            grid.update_powerlines_parameters(rs, xs, h1s, h2s)
        ref = _build()
        whole(ref)
        self._check(lambda g: g.update_one_powerline_parameters(line_id, r, x, h1, h2), ref)

    def test_one_trafo_matches_whole_vectors(self):
        # trafo 0 is the phase shifter whose impedance follows its shift angle
        for trafo_id in (0, 2):
            with self.subTest(trafo_id=trafo_id):
                r, x, h1, h2 = 0.01, 0.3, 0.001j, 0.002j

                def whole(grid):
                    trafos = grid.get_trafos()
                    rs = np.array([t.r_pu for t in trafos]); xs = np.array([t.x_pu for t in trafos])
                    h1s = np.array([t.h1_pu for t in trafos]); h2s = np.array([t.h2_pu for t in trafos])
                    rs[trafo_id], xs[trafo_id], h1s[trafo_id], h2s[trafo_id] = r, x, h1, h2
                    grid.update_trafos_parameters_full(rs, xs, h1s, h2s)
                ref = _build()
                whole(ref)
                self._check(lambda g: g.update_one_trafo_parameters_full(trafo_id, r, x, h1, h2), ref)


if __name__ == "__main__":
    unittest.main()
