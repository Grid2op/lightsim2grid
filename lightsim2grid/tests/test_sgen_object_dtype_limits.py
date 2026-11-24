# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`init_from_pandapower` accepts static generators whose limit columns have object dtype.

SimBench grids (and any net read from a file that stores missing limits as ``None``) carry
``min_p_mw``, ``max_p_mw``, ``min_q_mvar`` and ``max_q_mvar`` of the ``sgen`` table as
``object`` columns. ``np.isfinite`` refuses an object array, so the converter raised a
``TypeError`` before it built the grid. The generator converter already casts the same columns
to float.
"""

import unittest
import warnings

import numpy as np
import pandapower as pp

from lightsim2grid.network import init_from_pandapower


def _net(dtype):
    net = pp.create_empty_network()
    b0 = pp.create_bus(net, vn_kv=20.)
    b1 = pp.create_bus(net, vn_kv=20.)
    pp.create_ext_grid(net, b0)
    pp.create_line_from_parameters(net, b0, b1, 1., 0.1, 0.1, 10., 0.5)
    pp.create_load(net, b1, p_mw=1., q_mvar=0.2)
    pp.create_sgen(net, b1, p_mw=0.5, q_mvar=0.)

    for column in ("min_p_mw", "max_p_mw", "min_q_mvar", "max_q_mvar"):
        net.sgen[column] = np.array([None], dtype=dtype)

    return net


class TestSgenObjectDtypeLimits(unittest.TestCase):
    def test_object_dtype_none_limits(self):
        net = _net(object)
        self.assertEqual(net.sgen["min_p_mw"].dtype, object)
        pp.runpp(net)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid = init_from_pandapower(net)

        v = grid.ac_pf(np.ones(grid.total_bus(), dtype=complex), 20, 1e-8)
        self.assertEqual(v.shape[0], 2)
        np.testing.assert_allclose(np.abs(v), net.res_bus["vm_pu"].to_numpy(), atol=1e-6)


if __name__ == "__main__":
    unittest.main()
