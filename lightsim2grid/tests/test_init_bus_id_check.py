# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A grid built by hand with `LSGrid.init_*` is checked before it is solved (issue #259).

A bus id outside [0, total_bus()) used to be stored as given and only dereferenced by
the powerflow, which then crashed, returned an empty V, or silently ignored the
element. The init_* methods still store what they are given (the grid is incomplete
between two of them), but the first ac_pf / dc_pf / check_solution after one of them,
and a batch algorithm built on such a grid, run check_grid() first: the bad id is
refused with an IndexError (C++ std::out_of_range) instead.
"""

import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid, ContingencyAnalysisCPP, TimeSeriesCPP

N_BUS = 4
# beyond the grid, exactly total_bus(), a negative id that is not -1
BAD_IDS = (99, N_BUS, -2)
CASES = ("load", "gen", "shunt", "line_1", "line_2", "trafo_1", "trafo_2")


def _grid(bad=None, case=None):
    """The 4-bus grid of the issue; the bus of one element of `case` is `bad`."""
    def pick(name, default):
        return bad if case == name else default

    g = LSGrid()
    g.set_sn_mva(100.)
    g.init_bus(N_BUS, 1, np.array([138., 138., 138., 20.]), 3, 1)
    g.init_powerlines(np.array([.01, .02, .015]), np.array([.1, .15, .12]),
                      np.array([.02j, .03j, .025j]),
                      np.array([0, 1, pick("line_1", 0)]), np.array([1, 2, pick("line_2", 2)]))
    g.init_trafo(np.array([.005]), np.array([.08]), np.zeros(1, complex), np.ones(1), np.zeros(1),
                 [True], np.array([pick("trafo_1", 2)]), np.array([pick("trafo_2", 3)]), False)
    g.init_loads(np.array([60., 30.]), np.array([15., 10.]), np.array([2, pick("load", 3)]))
    g.init_generators(np.array([0., 40.]), np.array([1.02, 1.01]), np.array([-1e3, -1e3]),
                      np.array([1e3, 1e3]), np.array([0, pick("gen", 1)]))
    g.init_shunt(np.array([0.]), np.array([1.]), np.array([pick("shunt", 2)]))
    g.add_gen_slackbus(0, 1.)
    return g


def _ac(g):
    return g.ac_pf(np.ones(N_BUS, complex), 20, 1e-8)


class TestCheckGridAfterInit(unittest.TestCase):
    def test_valid_grid_converges(self):
        g = _grid()
        self.assertEqual(_ac(g).shape[0], N_BUS)
        self.assertEqual(_ac(g).shape[0], N_BUS)  # and again, the check is not in the way

    def test_issue_cases_refused_by_ac_pf(self):
        for case in CASES:
            for bad in BAD_IDS:
                with self.subTest(case=case, bad=bad):
                    g = _grid(bad, case)  # the init_* still accept it
                    with self.assertRaises(IndexError) as ctx:
                        _ac(g)
                    self.assertIn(f"bus id {bad} which is out of range [0, {N_BUS})", str(ctx.exception))
                    # still refused at the next call: the failed check does not clear the flag
                    with self.assertRaises(IndexError):
                        _ac(g)

    def test_dc_pf_and_check_solution(self):
        for case in CASES:
            with self.subTest(case=case):
                with self.assertRaises(IndexError):
                    _grid(99, case).dc_pf(np.ones(N_BUS, complex), 1, 1e-8)
                with self.assertRaises(IndexError):
                    _grid(99, case).check_solution(np.ones(N_BUS, complex), False)

    def test_batch_algorithms(self):
        for cls in (ContingencyAnalysisCPP, TimeSeriesCPP):
            with self.subTest(cls=cls.__name__):
                cls(_grid())
                with self.assertRaises(IndexError):
                    cls(_grid(99, "load"))

    def test_init_after_a_solve_is_checked_again(self):
        g = _grid()
        self.assertEqual(_ac(g).shape[0], N_BUS)
        g.init_loads(np.array([1., 2., 3.]), np.zeros(3), np.array([0, 1, 7]))
        with self.assertRaises(IndexError) as ctx:
            _ac(g)
        self.assertIn("id 2 is on bus id 7", str(ctx.exception))
        # fixed by another init_*: solves again
        g.init_loads(np.array([1., 2., 3.]), np.zeros(3), np.array([0, 1, 2]))
        self.assertEqual(_ac(g).shape[0], N_BUS)

    def test_substation_ids_left_from_a_smaller_container(self):
        # set_load_to_subid sized for 2 loads, then init_loads with more: init does not
        # touch the substation ids, and the automatic check_grid must not read past
        # their end (update_topo is what rejects the mismatch)
        g = _grid()
        g.set_load_to_subid(np.array([2, 3]))
        g.init_loads(np.ones(5), np.zeros(5), np.array([2, 3, 2, 3, 2]))
        self.assertEqual(_ac(g).shape[0], N_BUS)

    def test_copy_keeps_the_pending_check(self):
        g = _grid(99, "load")
        with self.assertRaises(IndexError):
            _ac(g.copy())

    def test_deactivated_bus_id_then_deactivate(self):
        # what the pypowsybl converter does for a disconnected element: -1 at init,
        # deactivated right after, before any powerflow
        g = _grid()
        g.init_loads(np.array([60., 30.]), np.array([15., 10.]), np.array([2, -1]))
        g.deactivate_load(1)
        self.assertEqual(_ac(g).shape[0], N_BUS)

    def test_deactivated_bus_id_left_connected(self):
        # -1 at init and never deactivated: connected on no bus, refused
        g = _grid()
        g.init_loads(np.array([60., 30.]), np.array([15., 10.]), np.array([2, -1]))
        with self.assertRaises(RuntimeError):
            _ac(g)


if __name__ == "__main__":
    unittest.main()
