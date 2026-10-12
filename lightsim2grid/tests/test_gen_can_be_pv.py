# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The per-generator ``can_be_pv`` flag (``LSGrid.set_gen_can_be_pv`` / ``GenInfo.can_be_pv``):
a hint that a PQ machine is one an outer loop pinned at a reactive limit. False by default,
never read by a powerflow, kept by ``copy``, pickle and the binary format."""

import os
import pickle
import tempfile
import unittest
import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid


def _feeder():
    """a 4-bus radial feeder 0-1-2-3, gen 0 (PV, slack) and gen 1 (PQ at its min_q)"""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
    grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                         np.array([0, 1, 2]), np.array([1, 2, 3]))
    grid.init_loads(np.array([80.]), np.array([20.]), np.array([3]))
    grid.init_generators_full(np.array([0., 10.]), np.array([1.02, 1.05]), np.array([0., -20.]),
                              [True, False], np.array([-1e3, -20.]), np.array([1e3, 20.]),
                              np.array([0, 1]))
    grid.add_gen_slackbus(0, 1.)
    grid.tell_solver_need_reset()
    return grid


class TestGenCanBePv(unittest.TestCase):
    def test_default_false(self):
        grid = _feeder()
        self.assertEqual([g.can_be_pv for g in grid.get_generators()], [False, False])

    def test_setter_and_info(self):
        grid = _feeder()
        grid.set_gen_can_be_pv(np.array([False, True]))
        self.assertEqual([g.can_be_pv for g in grid.get_generators()], [False, True])
        # a list works too, and the flag is not tied to the regulation flag
        grid.set_gen_can_be_pv([True, True])
        self.assertEqual([g.can_be_pv for g in grid.get_generators()], [True, True])
        self.assertEqual([g.voltage_regulator_on for g in grid.get_generators()], [True, False])

    def test_wrong_size_is_refused(self):
        grid = _feeder()
        with self.assertRaises(RuntimeError):
            grid.set_gen_can_be_pv(np.array([True, False, True]))
        self.assertEqual([g.can_be_pv for g in grid.get_generators()], [False, False])

    def test_powerflow_unchanged(self):
        grid = _feeder()
        v0 = np.full(grid.total_bus(), 1.0 + 0j)
        V_a = grid.ac_pf(1. * v0, 30, 1e-11)
        self.assertGreater(V_a.shape[0], 0)
        grid.set_gen_can_be_pv(np.array([False, True]))
        grid.tell_solver_need_reset()
        V_b = grid.ac_pf(1. * v0, 30, 1e-11)
        np.testing.assert_array_equal(V_a, V_b)

    def test_copy_pickle_binary(self):
        grid = _feeder()
        grid.set_gen_can_be_pv(np.array([False, True]))
        self.assertEqual([g.can_be_pv for g in grid.copy().get_generators()], [False, True])
        self.assertEqual([g.can_be_pv for g in pickle.loads(pickle.dumps(grid)).get_generators()],
                         [False, True])
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            loaded = LSGrid.load_binary(path)
        self.assertEqual([g.can_be_pv for g in loaded.get_generators()], [False, True])

    def test_reinit_resets(self):
        grid = _feeder()
        grid.set_gen_can_be_pv(np.array([True, True]))
        grid.init_generators_full(np.array([0., 10.]), np.array([1.02, 1.05]), np.array([0., -20.]),
                                  [True, False], np.array([-1e3, -20.]), np.array([1e3, 20.]),
                                  np.array([0, 1]))
        self.assertEqual([g.can_be_pv for g in grid.get_generators()], [False, False])


if __name__ == "__main__":
    unittest.main()
