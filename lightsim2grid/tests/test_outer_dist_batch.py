# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The NROuterDist_* algorithms (the outer loops around a distributed-slack Newton) in the batch
classes: each row equals the single solve of its own grid, with no more symbolic analyses and
factorizations than the plain Newton-Raphson batch. The equivalence harness and the rows are
test_outer_batch's; here one unit has little room above its target, so that the DistributedSlack
loop caps it in the rows and takes it out of the Newton's share."""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import test_outer_batch as tob  # noqa: E402

try:
    from lightsim2grid.algorithm import NROuterDist_KLU  # noqa: F401
    DIST = "NROuterDist_KLU"
except ImportError:
    DIST = "NROuterDist_SparseLU"

# case14's generator 1 (bus 2, 0 MW): its share of the losses is more than this
CAPPED_GEN, CAPPED_MAX_P = 1, 1.


class _Dist:
    ALGO = DIST

    def _grid(self, distributed=True):
        grid = tob._grid(distributed)
        nb_gen = len(grid.get_generators())
        max_p = np.full(nb_gen, 300.)
        max_p[CAPPED_GEN] = CAPPED_MAX_P
        grid.set_gen_p_limits(np.zeros(nb_gen), max_p)
        grid.change_algorithm(self.ALGO)
        return grid


class TestDistContingencyAnalysis(_Dist, tob.TestOuterContingencyAnalysis):
    def test_the_unit_is_capped(self):
        grid = self._grid()
        self.assertGreater(grid.ac_pf(tob._v0(grid), tob.MAX_IT, tob.TOL).size, 0)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertTrue(any(d.action == "PARK" and d.element_id == CAPPED_GEN for d in stats.decisions))
        p, _, _ = grid.get_gen_res()
        self.assertAlmostEqual(p[CAPPED_GEN], CAPPED_MAX_P, places=6)


class TestDistScenarioSweep(_Dist, tob.TestOuterScenarioSweep):
    pass


class TestDistInjectionSweep(_Dist, tob.TestOuterInjectionSweep):
    pass


class TestDistTimeSeries(_Dist, tob.TestOuterTimeSeries):
    pass


class TestDistContingencyAnalysisFromN(_Dist, tob.TestOuterContingencyAnalysisFromN):
    pass


class TestDistInjectionSweepFromN(_Dist, tob.TestOuterInjectionSweepFromN):
    pass


class TestDistResumeKeepsTheParkedUnits(_Dist, unittest.TestCase):
    """a solve resumed from another -- a batch row from its "n" solve, a time series' next step --
    starts with the units that solve took out of the Newton's share still out: the capped unit is
    not pushed past its limit again. DistributedSlack is the grid's only loop, so that any outer
    iteration is its own."""

    def _grid(self, distributed=True):
        from lightsim2grid.algorithm import DistributedSlack
        grid = super()._grid(distributed)
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack())
        return grid

    def _n_grid(self):
        grid = self._grid()
        v_n = grid.ac_pf(tob._v0(grid), tob.MAX_IT, tob.TOL)
        self.assertGreater(v_n.size, 0)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertTrue(any(d.action == "PARK" and d.element_id == CAPPED_GEN for d in stats.decisions))
        return grid, v_n

    @staticmethod
    def _more_load(grid, factor=1.05):
        for load_id, p in enumerate(np.asarray(grid.get_load_target_p())):
            grid.change_p_load(load_id, factor * float(p))

    def test_a_resumed_solve_that_changes_nothing(self):
        grid, v_n = self._n_grid()
        tob._resume_on(grid)
        V = grid.ac_pf(v_n, tob.MAX_IT, tob.TOL)
        self.assertGreater(V.size, 0)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.nb_outer_iterations, 0)
        self.assertFalse([d for d in stats.decisions if d.action in ("PARK", "UNPARK")])
        p, _, _ = grid.get_gen_res()
        self.assertAlmostEqual(p[CAPPED_GEN], CAPPED_MAX_P, places=6)
        np.testing.assert_allclose(V, v_n, rtol=0., atol=1e-10)

    def test_a_resumed_solve_with_more_load(self):
        """a time series' next step takes more load: the capped unit stays out, the others share
        the rest, the step needs no outer iteration"""
        grid, v_n = self._n_grid()
        tob._resume_on(grid)
        self._more_load(grid)
        V = grid.ac_pf(v_n, tob.MAX_IT, tob.TOL)
        self.assertGreater(V.size, 0)
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.nb_outer_iterations, 0)
        p, _, _ = grid.get_gen_res()
        self.assertAlmostEqual(p[CAPPED_GEN], CAPPED_MAX_P, places=6)
        # the same as solving that step afresh, where the loop caps the unit again
        fresh = self._grid()
        self._more_load(fresh)
        V_fresh = fresh.ac_pf(tob._v0(fresh), tob.MAX_IT, tob.TOL)
        self.assertGreater(fresh.get_algo().get_outer_loop_stats().nb_outer_iterations, 0)
        np.testing.assert_allclose(V, V_fresh, rtol=0., atol=1e-9)
        np.testing.assert_allclose(p, fresh.get_gen_res()[0], rtol=0., atol=1e-6)

    def test_rows_from_n_with_more_load(self):
        """the rows of a batch seeded from its "n" solve: the same, from the "n" solve's loops"""
        grid, v_n = self._n_grid()
        tob._resume_on(grid)
        self._more_load(grid)
        V_ref = grid.ac_pf(v_n, tob.MAX_IT, tob.TOL)
        self.assertGreater(V_ref.size, 0)
        p_ref = np.array(grid.get_gen_res()[0])
        base = self._grid()
        sweep = tob.InjectionSweepCPP(base)
        sweep.modify_load_p(np.outer(np.ones(tob.N_ROWS), 1.05 * np.asarray(base.get_load_target_p())))
        sweep.init_from_n_powerflow = True
        sweep.change_algorithm(self.ALGO)
        sweep.compute_gen_results = True
        sweep.compute(tob._v0(base), tob.MAX_IT, tob.TOL)
        gens = sweep.get_gen_results()
        for i, (res, V) in enumerate(zip(sweep.get_outer_results(), sweep.get_voltages())):
            with self.subTest(row=i):
                self.assertTrue(res.converged)
                self.assertEqual(res.nb_outer_iterations, 0)
                np.testing.assert_allclose(V, V_ref, rtol=0., atol=1e-10)
                np.testing.assert_allclose(gens[i, :, 0], p_ref, rtol=0., atol=1e-6)

class TestDistLostPowerParksTheUnitsItSaturates(unittest.TestCase):
    """a contingency row that loses a unit (redistribute_slack): the units the share of the lost
    power puts at a bound, before the row's first solve, are out of the Newton's share from the
    start. On a feeder that loses the unit next to its load the losses rise, the Newton would push
    them past that bound, and the loop would re-solve to cap them."""

    # buses 0-1-2-3 in a line, the load at 3: units A (70 MW) and C (20 MW, at most 25) at bus 1,
    # B (80 MW, tripped) at bus 2, the external grid at bus 0
    UNIT_C, UNIT_B, C_MAX_P = 1, 2, 25.

    def _grid(self):
        import warnings
        import pandapower as pp
        from lightsim2grid.network import init_from_pandapower
        from lightsim2grid.algorithm import DistributedSlack
        net = pp.create_empty_network()
        buses = [pp.create_bus(net, vn_kv=138.) for _ in range(4)]
        for i in range(3):
            pp.create_line(net, buses[i], buses[i + 1], length_km=40., std_type="149-AL1/24-ST1A 110.0")
        pp.create_ext_grid(net, buses[0], vm_pu=1.02)
        pp.create_gen(net, buses[1], p_mw=70., vm_pu=1.01)
        pp.create_gen(net, buses[1], p_mw=20., vm_pu=1.01)
        pp.create_gen(net, buses[2], p_mw=80., vm_pu=1.01)
        pp.create_load(net, buses[3], p_mw=170., q_mvar=20.)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pandapower(net)
        nb_gen = len(grid.get_generators())
        max_p = np.full(nb_gen, 300.)
        max_p[self.UNIT_C] = self.C_MAX_P
        grid.set_gen_p_limits(np.zeros(nb_gen), max_p)
        grid.set_gen_can_participate_slack([True] * nb_gen, np.full(nb_gen, 300.))
        grid.clear_outer_loops()
        grid.add_outer_loop(DistributedSlack())
        grid.change_algorithm(DIST)
        return grid

    def test_no_outer_iteration(self):
        base = self._grid()
        v0 = np.full(base.total_bus(), 1.02 + 0j)
        self.assertGreater(base.ac_pf(v0, tob.MAX_IT, tob.TOL).size, 0)
        p_n = np.array(base.get_gen_res()[0])
        self.assertLess(p_n[self.UNIT_C], self.C_MAX_P)  # not at its limit in the "n" solve
        ca = tob.ContingencyAnalysisCPP(base)
        ca.handle_disconnected_grid = True
        ca.init_from_n_powerflow = True
        ca.redistribute_slack = True
        ca.add_contingency_gen([self.UNIT_B])
        ca.change_algorithm(DIST)
        ca.compute_gen_results = True
        ca.compute(v0, tob.MAX_IT, tob.TOL)
        res = ca.get_outer_results()[0]
        self.assertTrue(res.converged)
        self.assertEqual(res.nb_outer_iterations, 0)
        p = ca.get_gen_results()[0, :, 0]
        self.assertAlmostEqual(p[self.UNIT_B], 0., places=9)
        self.assertAlmostEqual(p[self.UNIT_C], self.C_MAX_P, places=6)
        # the two others shared the rest alike (same weight), more than the lost power alone
        others = [g for g in range(len(p)) if g not in (self.UNIT_B, self.UNIT_C)]
        shares = p[others] - p_n[others]
        np.testing.assert_allclose(shares, shares[0], rtol=0., atol=1e-6)
        self.assertGreater(shares.sum() + (self.C_MAX_P - p_n[self.UNIT_C]), p_n[self.UNIT_B])

if __name__ == "__main__":
    unittest.main()
