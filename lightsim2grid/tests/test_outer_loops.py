# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The NROuter_* algorithms (OpenLoadFlow's outer loops around a single-slack Newton), seen
from python. The driver itself is tested in C++ (src/tests/test_outer_loop_driver.cpp) with
scripted loops; each loop has its own test file."""

import copy
import unittest
import warnings

import numpy as np
import pandapower.networks as pn

from lightsim2grid.network import init_from_pandapower
from lightsim2grid.algorithm import ErrorType, OuterLoopStatus

try:
    from lightsim2grid.algorithm import NROuter_KLU  # noqa: F401
    OUTER, SING = "NROuter_KLU", "NRSing_KLU"
except ImportError:
    OUTER, SING = "NROuter_SparseLU", "NRSing_SparseLU"


def _grid(algo):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(pn.case118())
    grid.change_algorithm(algo)
    return grid


def _solve(grid):
    return grid.ac_pf(np.full(grid.total_bus(), 1.04 + 0j), 30, 1e-8)


class TestNROuter(unittest.TestCase):
    def test_without_loop_it_is_the_single_slack_newton(self):
        V_ref = _solve(_grid(SING))
        grid = _grid(OUTER)
        grid.clear_outer_loops()
        V = _solve(grid)
        self.assertEqual(V.shape, V_ref.shape)
        self.assertTrue(np.array_equal(V, V_ref))
        stats = grid.get_algo().get_outer_loop_stats()
        self.assertEqual(stats.status, OuterLoopStatus.STABLE)
        self.assertEqual(stats.nb_outer_iterations, 0)
        self.assertEqual(len(stats.nr_iterations), 1)
        self.assertTrue(grid.get_algo().supports_outer_loops())

    def test_one_analysis_over_several_solves(self):
        grid = _grid(OUTER)
        for _ in range(3):
            self.assertGreater(_solve(grid).shape[0], 0)
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

    def test_loop_list(self):
        grid = _grid(OUTER)
        default = grid.get_outer_loops()
        grid.clear_outer_loops()
        self.assertEqual(grid.get_outer_loops(), [])
        grid.reset_outer_loops()
        self.assertEqual([l.name() for l in grid.get_outer_loops()], [l.name() for l in default])
        with self.assertRaises(RuntimeError):
            grid.add_outer_loop(None)

    def test_driver_parameters_round_trip(self):
        grid = _grid(OUTER)
        cfg = grid.get_ac_algo_config()
        # the Newton's 4 + 6 parameters, then the driver's 2 + 3
        self.assertEqual(len(cfg.int_params), 6)
        self.assertEqual(len(cfg.real_params), 9)
        self.assertEqual(cfg.int_params[4], 30)
        self.assertAlmostEqual(cfg.real_params[6], 0.8)
        self.assertAlmostEqual(cfg.real_params[7], 1.2)
        self.assertAlmostEqual(cfg.real_params[8], 180.)
        ip = list(cfg.int_params)
        ip[4] = 12
        cfg.int_params = ip
        grid.set_ac_algo_config(cfg)
        self.assertEqual(grid.get_ac_algo_config().int_params[4], 12)
        # a copy keeps the algorithm and its configuration
        grid2 = grid.copy()
        self.assertTrue(grid2.get_algo().supports_outer_loops())
        self.assertEqual(grid2.get_ac_algo_config().int_params[4], 12)
        self.assertEqual(copy.deepcopy(grid).get_ac_algo_config().int_params[4], 12)

    def test_unrealistic_state(self):
        grid = _grid(OUTER)
        cfg = grid.get_ac_algo_config()
        rp = list(cfg.real_params)
        rp[6] = 1.2  # nothing is realistic any more ...
        rp[7] = 1.3
        rp[8] = 0.   # ... at any nominal voltage
        cfg.real_params = rp
        grid.set_ac_algo_config(cfg)
        self.assertEqual(_solve(grid).shape[0], 0)
        self.assertEqual(grid.get_algo().get_error(), ErrorType.UnrealisticState)
        self.assertTrue(grid.get_algo().get_outer_loop_stats().unrealistic_state)

    def test_standalone_algorithm(self):
        from lightsim2grid import algorithm
        algo = getattr(algorithm, OUTER)()
        algo.max_outer_iterations = 7
        self.assertEqual(algo.max_outer_iterations, 7)
        with self.assertRaises(RuntimeError):
            algo.min_realistic_voltage = 2.  # above the max
        self.assertEqual(algo.get_outer_loop_stats().nb_outer_iterations, 0)


if __name__ == "__main__":
    unittest.main()
