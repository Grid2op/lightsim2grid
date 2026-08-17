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
        # a negative position is rejected by the setter itself (the one bound it
        # can check locally, unlike the upper bound -- see
        # test_out_of_range_pos_topo_vect_is_rejected below, which still needs
        # check_grid() since that bound needs cross-container context), so
        # check_grid() never even gets a chance to see it here.
        pos = np.arange(self.n_load, dtype=np.int32)
        pos[0] = -1
        with self.assertRaises(IndexError):
            self.grid.set_load_pos_topo_vect(pos)

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
        # same reasoning as test_negative_pos_topo_vect_is_rejected above: the
        # setter rejects a negative id immediately (n_sub, the upper bound, is
        # not available to it), so check_grid() is never reached here.
        for name, n in [("set_load_to_subid", self.n_load),
                        ("set_gen_to_subid", self.n_gen)]:
            with self.subTest(name):
                sub = np.zeros(n, dtype=np.int32)
                sub[0] = -1
                with self.assertRaises(IndexError):
                    getattr(self.grid, name)(sub)
                # the rejected call never mutated the container, but re-apply a
                # valid mapping anyway so the next subTest starts from a known state
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


class TestSubstationAndMappingConsistency(unittest.TestCase):
    """
    The substation container and the bus-id mapping vectors.

    ``SubstationContainer`` is the root of the grid's index space: ``nb_bus()``
    is the bound ``check_grid()`` validates every element's bus id against, and
    ``bus_status_`` is the vector those same ids are used to index. They used to
    be restored from a pickle / binary file with no cross-check at all, so a file
    whose two lengths disagreed passed ``check_grid()`` and then corrupted the
    heap on the next ``init_bus_status()``. Same idea for ``_ls_to_orig``, whose
    *values* are used as indices into the reverse mapping it allocates.

    These are built by hand (no pandapower) so the checks are exercised directly.
    """

    def _grid(self, n_sub=3, n_busbar=2):
        from lightsim2grid.lightsim2grid_cpp import LSGrid
        grid = LSGrid()
        grid.init_bus(n_sub, n_busbar,
                      np.array([138.0] * (n_sub * n_busbar)), 0, 0)
        grid.init_loads(np.array([1.0]), np.array([0.5]),
                        np.array([0], dtype=np.int32))
        grid.init_generators(np.array([1.0]), np.array([1.0]),
                             np.array([-10.0]), np.array([10.0]),
                             np.array([1], dtype=np.int32))
        return grid

    def test_hand_built_grid_is_consistent(self):
        # the new substation / mapping checks must not reject a legitimate grid
        # (sub_vn_kv and the substation names are optional and normally empty).
        self.assertIsNone(self._grid().check_grid())

    def test_substation_info_without_names(self):
        # `sub_names_` is optional -- the pandapower / matpower / powermodels
        # loaders never set it -- and SubstationInfo used to read sub_names_[id]
        # unconditionally, ie read a std::string out of bounds. Iterating the
        # substations of such a grid used to segfault.
        grid = self._grid()
        self.assertEqual(list(grid.get_substation_names()), [])
        for sub in grid.get_substations():
            self.assertEqual(sub.name, "")
            self.assertEqual(sub.vn_kv, 138.0)

    def test_negative_ls_to_orig_is_rejected(self):
        # the reverse mapping is sized from max(abs(values)) and then indexed with
        # the values themselves: -5 sizes it from 5 and writes at index -5.
        grid = self._grid()
        bad = np.full(grid.total_bus(), -1, dtype=np.int32)
        bad[0] = -5
        with self.assertRaises(IndexError):
            grid._ls_to_orig = bad

    def test_absurd_ls_to_orig_is_rejected(self):
        grid = self._grid()
        bad = np.full(grid.total_bus(), -1, dtype=np.int32)
        bad[0] = np.iinfo(np.int32).max
        with self.assertRaises(IndexError):
            grid._ls_to_orig = bad

    def test_valid_ls_to_orig_is_accepted(self):
        grid = self._grid()
        ok = np.arange(grid.total_bus(), dtype=np.int32)
        grid._ls_to_orig = ok
        np.testing.assert_array_equal(np.asarray(grid._ls_to_orig), ok)
        self.assertIsNone(grid.check_grid())

    def test_orig_to_ls_is_the_inverse_of_ls_to_orig(self):
        # orig bus 0 has no lightsim counterpart; orig 1..6 map to ls 0..5.
        grid = self._grid()
        orig_to_ls = np.array([-1, 0, 1, 2, 3, 4, 5], dtype=np.int32)
        grid._orig_to_ls = orig_to_ls
        # the inverse: lightsim bus i corresponds to original bus i + 1
        np.testing.assert_array_equal(np.asarray(grid._ls_to_orig),
                                      np.arange(1, 7, dtype=np.int32))

    def test_out_of_range_orig_to_ls_is_rejected(self):
        grid = self._grid()
        bad = np.arange(grid.total_bus(), dtype=np.int32)
        bad[0] = BIG  # used to index a total_bus()-long _ls_to_orig
        with self.assertRaises(IndexError):
            grid._orig_to_ls = bad

    def test_update_topo_rejects_short_arrays(self):
        # each container indexes has_changed / new_values by position in the
        # topology vector with an unchecked Eigen operator(): a short array was
        # simply read past its end.
        grid = self._grid()
        grid.set_load_pos_topo_vect(np.array([0], dtype=np.int32))
        grid.set_gen_pos_topo_vect(np.array([1], dtype=np.int32))
        with self.assertRaises(RuntimeError):
            grid.update_topo(np.zeros(0, dtype=bool), np.zeros(0, dtype=np.int32))
        # correctly-sized (dim_topo == 2 here: one load + one gen), no-op call
        grid.update_topo(np.zeros(2, dtype=bool), np.ones(2, dtype=np.int32))

    def test_update_slack_weights_rejects_short_array(self):
        grid = self._grid()  # exactly one generator
        with self.assertRaises(RuntimeError):
            grid.update_slack_weights(np.zeros(0, dtype=bool))
        grid.update_slack_weights(np.ones(1, dtype=bool))

    def test_poisoned_substation_state_is_rejected(self):
        # go through the real pickle path: tamper with the substation block of a
        # grid's state, then load it. Every case must raise, never corrupt memory.
        import copyreg
        import io
        import pickle
        from lightsim2grid.lightsim2grid_cpp import LSGrid

        SUBSTATION_ID = 7
        # SubstationContainer::StateRes = (n_sub, nmax_busbar, sub_vn_kv,
        #                                  bus_status, bus_vn_kv, sub_names,
        #                                  bus_vmin_kv, bus_vmax_kv)
        def dump_with(grid, mutate):
            outer = grid.__getstate__()
            inner = list(outer[3])
            sub = list(inner[SUBSTATION_ID])
            mutate(sub)
            inner[SUBSTATION_ID] = tuple(sub)
            state = (outer[0], outer[1], outer[2], tuple(inner))

            class _P(pickle.Pickler):
                def reducer_override(self, obj):
                    if isinstance(obj, LSGrid):
                        return (copyreg.__newobj__, (LSGrid,), state)
                    return NotImplemented

            buf = io.BytesIO()
            _P(buf, protocol=4).dump(grid)
            return buf.getvalue()

        def set_bus_status_short(sub):
            sub[3] = [True]

        def set_nsub_zero(sub):
            sub[0] = 0

        def set_nsub_overflow(sub):
            sub[0] = 100000
            sub[1] = 100000

        def set_bus_vn_kv_short(sub):
            sub[4] = [138.0]

        for name, mutate in [("bus_status too short", set_bus_status_short),
                             ("n_sub == 0", set_nsub_zero),
                             ("n_sub * nmax overflows", set_nsub_overflow),
                             ("bus_vn_kv too short", set_bus_vn_kv_short)]:
            with self.subTest(name):
                blob = dump_with(self._grid(), mutate)
                with self.assertRaises(RuntimeError):
                    pickle.loads(blob)

    def test_pickle_roundtrip_of_hand_built_grid(self):
        import pickle
        grid = self._grid()
        restored = pickle.loads(pickle.dumps(grid))
        self.assertIsNone(restored.check_grid())
        self.assertEqual(restored.total_bus(), grid.total_bus())


if __name__ == "__main__":
    unittest.main()
