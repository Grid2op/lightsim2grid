# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The pandapower converter keeps the in_service status of the generators when it appends
one generator per ext_grid (no generator stands on the slack bus): the offline ones stay
offline, and the powerflow matches pandapower's."""

import unittest
import warnings

import numpy as np
import pandapower as pp

from lightsim2grid.network import init_from_pandapower


class TestExtGridKeepsGenStatus(unittest.TestCase):
    def _net(self, statuses):
        net = pp.create_empty_network()
        pp.create_buses(net, 3, 110)
        pp.create_ext_grid(net, 0, vm_pu=1.02)
        pp.create_line_from_parameters(net, 0, 1, 10, .1, .2, 0, 1)
        pp.create_line_from_parameters(net, 1, 2, 10, .1, .2, 0, 1)
        pp.create_load(net, 1, 10, 5)
        pp.create_load(net, 2, 1, .5)
        for j, status in enumerate(statuses):
            pp.create_gen(net, 2 - j, p_mw=5, vm_pu=1.10 - .09 * j,
                          in_service=status, min_q_mvar=-100, max_q_mvar=100)
        pp.runpp(net, numba=False, lightsim2grid=False, tolerance_mva=1e-8)
        return net

    def _check(self, statuses):
        net = self._net(statuses)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = init_from_pandapower(net)
        gens = model.get_generators()
        self.assertEqual([gens[j].connected for j in range(len(statuses))], statuses)
        self.assertTrue(gens[len(statuses)].connected)  # the generator added for the ext_grid
        expected = net._ppc["internal"]["V"]
        got = model.ac_pf(expected.copy(), 30, 1e-8)
        self.assertEqual(len(got), len(expected))
        self.assertLess(float(np.max(np.abs(got - expected))), 1e-8)

    def test_offline_generator_stays_offline(self):
        self._check([False])

    def test_online_generator_stays_online(self):
        self._check([True])

    def test_mixed_generator_statuses(self):
        self._check([False, True])


if __name__ == "__main__":
    unittest.main()
