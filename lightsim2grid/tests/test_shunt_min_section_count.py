# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`LSGrid.set_shunt_sections` forwards `min_section_count`: the fewest sections the
ShuntVoltageControl outer loop may leave on (1 for OpenLoadFlow with
allowNonLinearShuntZeroSection=false)."""

import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid
from lightsim2grid.algorithm import ShuntVoltageControl

N_SECTIONS = 4


def _grid(min_section_count=None, target_vm_pu=0.80):
    """buses 0-1-2 in a line, the slack generator at 0, a load and a shunt of `N_SECTIONS`
    sections of 5 MVar at 2. The target is below the voltage the load leaves, so the loop
    wants no section on."""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.)
    grid.init_bus(3, 1, np.full(3, 138.), 0, 0)
    grid.init_powerlines(np.full(2, .01), np.full(2, .1), np.zeros(2, complex), np.array([0, 1]), np.array([1, 2]))
    grid.init_loads(np.array([60.]), np.array([40.]), np.array([2]))
    grid.init_shunt(np.zeros(1), np.zeros(1), np.array([2]))
    grid.init_generators(np.zeros(1), np.ones(1), np.full(1, -1000.), np.full(1, 1000.), np.array([0]))
    grid.add_gen_slackbus(0, 1.)
    p = [0.] * N_SECTIONS
    q = [-5. * k for k in range(1, N_SECTIONS + 1)]
    if min_section_count is None:
        grid.set_shunt_sections(0, 2, p, q)
    else:
        grid.set_shunt_sections(0, 2, p, q, min_section_count=min_section_count)
    grid.set_shunt_section_regulation(0, True, target_vm_pu, 0., 2)
    grid.change_algorithm("NROuter_SparseLU")
    grid.clear_outer_loops()
    grid.add_outer_loop(ShuntVoltageControl())
    return grid


def _solve(grid):
    V = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 30, 1e-8)
    assert V.shape[0] > 0
    return grid.get_shunts()[0]


class TestShuntMinSectionCount(unittest.TestCase):
    def test_default_is_zero_and_the_loop_may_switch_everything_off(self):
        info = _solve(_grid())
        self.assertEqual(info.min_section_count, 0)
        self.assertEqual(info.res_section_count, 0)

    def test_min_section_count_is_kept_by_the_loop(self):
        info = _solve(_grid(min_section_count=1))
        self.assertEqual(info.min_section_count, 1)
        self.assertEqual(info.res_section_count, 1)

    def test_min_section_count_is_clamped_to_the_table(self):
        self.assertEqual(_grid(min_section_count=N_SECTIONS + 5).get_shunts()[0].min_section_count, N_SECTIONS)
        self.assertEqual(_grid(min_section_count=-3).get_shunts()[0].min_section_count, 0)


if __name__ == "__main__":
    unittest.main()
