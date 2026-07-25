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
        self.n_gen = len(self.grid.get_generators())
        self.n_sub = self.grid.get_n_sub()
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

    # ------------------------------------------------------------------
    # substation ids (`*_to_subid`), same treatment as pos_topo_vect above
    # ------------------------------------------------------------------
    def test_out_of_range_subid_is_rejected(self):
        # one element sits on a substation id far above the substation count.
        for name, n in [("set_load_to_subid", self.n_load),
                        ("set_gen_to_subid", self.n_gen)]:
            with self.subTest(name):
                sub = np.zeros(n, dtype=np.int32)
                sub[-1] = BIG
                getattr(self.grid, name)(sub)
                with self.assertRaises(IndexError):
                    self.grid.check_grid()
                # restore a valid mapping so the next subTest starts clean
                getattr(self.grid, name)(np.zeros(n, dtype=np.int32))

    def test_negative_subid_is_rejected(self):
        for name, n in [("set_load_to_subid", self.n_load),
                        ("set_gen_to_subid", self.n_gen)]:
            with self.subTest(name):
                sub = np.zeros(n, dtype=np.int32)
                sub[0] = -1
                getattr(self.grid, name)(sub)
                with self.assertRaises(IndexError):
                    self.grid.check_grid()
                getattr(self.grid, name)(np.zeros(n, dtype=np.int32))

    def test_subid_just_out_of_range_is_rejected(self):
        # off-by-one: `n_sub` itself is the first invalid substation id.
        sub = np.zeros(self.n_load, dtype=np.int32)
        sub[-1] = self.n_sub
        self.grid.set_load_to_subid(sub)
        with self.assertRaises(IndexError):
            self.grid.check_grid()

    def test_valid_subid_is_accepted(self):
        # every element on the last valid substation: in range, must not raise.
        self.grid.set_load_to_subid(np.full(self.n_load, self.n_sub - 1, dtype=np.int32))
        self.assertIsNone(self.grid.check_grid())

    # ------------------------------------------------------------------
    # solver policy parameters (the `lightsim2grid.algorithm` / newtonpf entry
    # point, and the AlgoConfig carried by the serialized grid state)
    # ------------------------------------------------------------------
    def test_refactor_every_n_zero_is_rejected(self):
        # `refactor_every_n` is used as `iter % refactor_every_n`: 0 used to be an
        # integer division by zero, ie a SIGFPE that killed the interpreter.
        from lightsim2grid.algorithm import NR_SparseLU
        solver = NR_SparseLU()
        for bad in (0, -1):
            with self.subTest(bad):
                with self.assertRaises(RuntimeError):
                    solver.set_refactor_every_n(bad)
        solver.set_refactor_every_n(3)  # a sane value still works
        self.assertEqual(solver.get_refactor_every_n(), 3)

    def test_line_search_params_are_validated(self):
        from lightsim2grid.algorithm import NR_SparseLU
        solver = NR_SparseLU()
        for name, bad_values in [("set_ls_rho", (0.0, 1.0, 1.5, float("nan"))),
                                 ("set_ls_c", (0.0, 1.0, float("inf"))),
                                 ("set_ls_max_iter", (-1,)),
                                 ("set_max_dVa", (0.0, -1.0, float("nan"))),
                                 ("set_max_dVm", (0.0, float("inf")))]:
            for bad in bad_values:
                with self.subTest(f"{name}({bad})"):
                    with self.assertRaises(RuntimeError):
                        getattr(solver, name)(bad)

    def test_poisoned_algo_config_is_rejected(self):
        # AlgoConfig is part of the serialized grid state, so this value can come
        # straight from a pickle or a binary file.
        cfg = self.grid.get_ac_algo_config()
        int_params = list(cfg.int_params)
        int_params[3] = 0  # refactor_every_n
        cfg.int_params = int_params
        with self.assertRaises(RuntimeError):
            self.grid.set_ac_algo_config(cfg)

    def test_every_registered_solver_name_is_valid(self):
        # A solver name is persisted in every saved grid and re-selects the solver
        # on load, so it is restricted to [A-Za-z_][A-Za-z0-9_.]{0,63} (rejecting
        # non-ASCII homoglyphs of a built-in name, control characters, and
        # unbounded lengths). Guards against a built-in ever gaining a bad name.
        import re
        pattern = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]{0,63}$")
        names = self.grid.available_algorithm_names()
        self.assertGreater(len(names), 0)
        for name in names:
            with self.subTest(name):
                self.assertRegex(name, pattern)

    def test_load_binary_without_algorithm(self):
        # the escape hatch for a grid saved with a solver that is not available
        # here (typically a plugin that has not been loaded): the grid data loads,
        # only the solver selection is skipped.
        import tempfile, os
        from lightsim2grid.lightsim2grid_cpp import LSGrid
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "grid.lsb")
            self.grid.save_binary(path)
            loaded = LSGrid.load_binary_without_algorithm(path)
            self.assertIsNone(loaded.check_grid())
            self.assertEqual(len(loaded.get_loads()), self.n_load)
            # and it is usable: a powerflow runs with the default solver
            V = loaded.ac_pf(np.ones(loaded.total_bus(), dtype=complex), 20, 1e-8)
            self.assertGreater(V.shape[0], 0)

    def test_pickle_roundtrip_of_valid_grid(self):
        # pickling goes through set_state on load, which runs check_grid(); a valid
        # grid must round-trip cleanly.
        import pickle
        restored = pickle.loads(pickle.dumps(self.grid))
        self.assertIsNone(restored.check_grid())


if __name__ == "__main__":
    unittest.main()
