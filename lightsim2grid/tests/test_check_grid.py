# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Tests for ``LSGrid.check_grid()`` (the whole-grid consistency validator) from
Python.

``check_grid()`` range-checks every index the grid carries (bus / substation /
topology-vector) and the finiteness of the physical inputs; it is called
automatically when a grid is loaded (pickle / binary format, via ``set_state``)
and by every grid loader (pandapower / pypowsybl / matpower / powermodels). The
deep ``set_state`` poison cases are covered by the C++ suite (test_check_grid.cpp,
also run under ASan/UBSan/valgrind); here we check the Python-visible behaviour:

  * a real grid loaded from pandapower passes (the loader calls check_grid, so a
    false positive would make the loader itself raise);
  * check_grid() is reachable from Python and rejects an inconsistent
    topology-vector with a clean exception (``IndexError`` for an out-of-range
    index -- C++ ``std::out_of_range`` --, ``RuntimeError`` for a duplicate --
    C++ ``std::runtime_error``), never a crash.
"""

import unittest
import warnings
import numpy as np

import pandapower.networks as pn
from lightsim2grid.network import init_from_pandapower


# an id far above any element count, but still inside the C `int` range
BIG = 10 ** 6


class TestCheckGrid(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            net = pn.case14()
            self.grid = init_from_pandapower(net)  # the loader calls check_grid()
        self.n_load = len(self.grid.get_loads())
        self.assertGreater(self.n_load, 1)  # case14 has several loads

    def test_loader_accepts_valid_grid(self):
        # init_from_pandapower already ran check_grid() in setUp without raising;
        # calling it again explicitly must also succeed and return None.
        self.assertIsNone(self.grid.check_grid())

    def test_duplicated_pos_topo_vect_is_rejected(self):
        # put every load at the same position in the topology vector: the collected
        # pos_topo_vect entries can then no longer be a permutation of [0, dim_topo).
        self.grid.set_load_pos_topo_vect(np.zeros(self.n_load, dtype=np.int32))
        with self.assertRaises(RuntimeError):
            self.grid.check_grid()

    def test_out_of_range_pos_topo_vect_is_rejected(self):
        # distinct positions, but one is far out of range.
        pos = np.arange(self.n_load, dtype=np.int32)
        pos[-1] = BIG
        self.grid.set_load_pos_topo_vect(pos)
        with self.assertRaises(IndexError):
            self.grid.check_grid()

    def test_negative_pos_topo_vect_is_rejected(self):
        pos = np.arange(self.n_load, dtype=np.int32)
        pos[0] = -1
        self.grid.set_load_pos_topo_vect(pos)
        with self.assertRaises(IndexError):
            self.grid.check_grid()

    def test_pickle_roundtrip_of_valid_grid(self):
        # pickling goes through set_state on load, which runs check_grid(); a valid
        # grid must round-trip cleanly.
        import pickle
        restored = pickle.loads(pickle.dumps(self.grid))
        self.assertIsNone(restored.check_grid())


if __name__ == "__main__":
    unittest.main()
