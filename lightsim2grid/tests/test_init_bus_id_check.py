# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The `LSGrid.init_*` methods reject a bus id outside [0, total_bus()) (issue #259).

Such an id used to be stored as given and only dereferenced by the powerflow, which
then crashed, returned an empty V, or silently ignored the element. It must now raise
at init (`std::out_of_range`, an `IndexError` in Python), naming the element and the
id, and leave the container as it was. -1 (`_deactivated_bus_id`) stays accepted:
the pypowsybl converter passes it for a disconnected element and deactivates it
right after.
"""

import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid

N_BUS = 4
# beyond the grid, exactly total_bus(), a negative id that is not -1
BAD_IDS = (99, N_BUS, -2)


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
    return g


class TestInitBusIdCheck(unittest.TestCase):
    def test_valid_grid_converges(self):
        g = _grid()
        g.add_gen_slackbus(0, 1.)
        V = g.ac_pf(np.ones(N_BUS, complex), 20, 1e-8)
        self.assertEqual(V.shape[0], N_BUS)

    def test_issue_cases_raise_at_init(self):
        for case in ("load", "gen", "shunt", "line_1", "line_2", "trafo_1", "trafo_2"):
            for bad in BAD_IDS:
                with self.subTest(case=case, bad=bad):
                    with self.assertRaises(IndexError) as ctx:
                        _grid(bad, case)
                    msg = str(ctx.exception)
                    self.assertIn(f" is {bad},", msg)
                    self.assertIn(f"< {N_BUS}", msg)

    def test_message_names_the_element(self):
        g = _grid()
        with self.assertRaises(IndexError) as ctx:
            g.init_loads(np.array([1., 2., 3.]), np.zeros(3), np.array([0, 1, 7]))
        msg = str(ctx.exception)
        self.assertIn("LSGrid::init_loads", msg)
        self.assertIn("`loads_bus_id` of element 2 is 7", msg)

    def test_failed_init_leaves_container_unchanged(self):
        g = _grid()
        with self.assertRaises(IndexError):
            g.init_loads(np.array([1., 2., 3.]), np.zeros(3), np.array([0, 1, 7]))
        self.assertEqual(len(g.get_loads()), 2)
        g.add_gen_slackbus(0, 1.)
        V = g.ac_pf(np.ones(N_BUS, complex), 20, 1e-8)
        self.assertEqual(V.shape[0], N_BUS)

    def test_deactivated_bus_id_accepted(self):
        g = _grid()
        g.init_loads(np.array([60., 30.]), np.array([15., 10.]), np.array([2, -1]))
        g.deactivate_load(1)
        g.add_gen_slackbus(0, 1.)
        V = g.ac_pf(np.ones(N_BUS, complex), 20, 1e-8)
        self.assertEqual(V.shape[0], N_BUS)

    def test_other_inits(self):
        g = _grid()
        one = np.ones(1)
        zero = np.zeros(1)
        for bad in BAD_IDS:
            with self.subTest(bad=bad):
                b = np.array([bad])
                with self.assertRaises(IndexError):
                    g.init_sgens(one, zero, zero, one, -one, one, b)
                with self.assertRaises(IndexError):
                    g.init_storages(one, zero, b)
                with self.assertRaises(IndexError):  # svc bus
                    g.init_svcs([0], one, zero, zero, -one, one, np.array([2]), b)
                with self.assertRaises(IndexError):  # regulated bus
                    g.init_svcs([0], one, zero, zero, -one, one, b, np.array([2]))
                with self.assertRaises(IndexError):
                    g.init_dclines(np.array([0]), b, one, zero, zero, one, one,
                                   -one, one, -one, one)

    def test_init_before_init_bus(self):
        g = LSGrid()
        with self.assertRaises(IndexError) as ctx:
            g.init_loads(np.ones(1), np.zeros(1), np.array([0]))
        self.assertIn("init_bus", str(ctx.exception))
        # empty containers are fine without any bus
        g.init_loads(np.zeros(0), np.zeros(0), np.zeros(0, dtype=int))


if __name__ == "__main__":
    unittest.main()
