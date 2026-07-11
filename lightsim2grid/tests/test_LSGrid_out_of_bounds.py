# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Regression tests for the memory-safety fixes on the Python-exposed `LSGrid` /
`TimeSeriesCPP` methods.

Before the fix, passing an out-of-range or negative element id (or a mis-shaped
array) to these methods reached unchecked C++ indexing (``operator[]`` / Eigen
``operator()`` / ``.col()``, whose bounds asserts are compiled out in release
builds) and could read or write out of bounds. They must now raise a clean Python
exception instead:

  * ``IndexError``   for out-of-range / negative element ids (``_check_in_range``
    and the ``LSGrid``-level bus check throw ``std::out_of_range``);
  * ``RuntimeError`` for array shape / length mismatches (``check_size``,
    ``update_continuous_values`` and ``TimeSeries::compute_Vs`` throw
    ``std::runtime_error``).
"""

import unittest
import warnings
import numpy as np

import pandapower.networks as pn
from lightsim2grid.network import init_from_pandapower

try:
    from lightsim2grid.lightsim2grid_cpp import TimeSeriesCPP
    TimeSeriesCPP_AVAILABLE = True
except ImportError:
    TimeSeriesCPP_AVAILABLE = False


# an id far above any element count, but still inside the C `int` range
BIG = 10 ** 6


class TestOutOfBoundsBindings(unittest.TestCase):
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            net = pn.case14()
            self.grid = init_from_pandapower(net)

        # element counts as seen by the *gridmodel* (e.g. the pandapower ext_grid
        # is modelled as an extra generator, so these differ from len(net.gen)).
        self.n_gen = len(self.grid.get_generators())
        self.n_sgen = len(self.grid.get_static_generators())
        self.n_load = len(self.grid.get_loads())
        self.n_shunt = len(self.grid.get_shunts())
        self.n_storage = len(self.grid.get_storages())
        self.n_line = len(self.grid.get_lines())
        self.n_trafo = len(self.grid.get_trafos())
        self.n_bus = self.grid.total_bus()

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _assert_raises_each(self, exc, calls):
        """`calls` is an iterable of (description, callable); each must raise `exc`."""
        for descr, fun in calls:
            with self.subTest(descr):
                with self.assertRaises(exc):
                    fun()

    # ------------------------------------------------------------------
    # out-of-bounds READS: deactivate / reactivate / change_bus family
    # ------------------------------------------------------------------
    def test_deactivate_reactivate_out_of_bounds(self):
        # one-side elements (loads, gens, sgens, storages, shunts, svc) and
        # two-side branches (powerlines, trafos). Empty containers (e.g. sgen /
        # storage / svc in case14) still raise, so no special-casing is needed.
        families = ["load", "gen", "sgen", "storage", "shunt", "svc",
                    "powerline", "trafo"]
        calls = []
        for fam in families:
            for bad in (-1, BIG):
                calls.append((f"deactivate_{fam}({bad})",
                              lambda fam=fam, bad=bad: getattr(self.grid, f"deactivate_{fam}")(bad)))
                calls.append((f"reactivate_{fam}({bad})",
                              lambda fam=fam, bad=bad: getattr(self.grid, f"reactivate_{fam}")(bad)))
        self._assert_raises_each(IndexError, calls)

    def test_change_bus_out_of_bounds(self):
        valid_bus = 1  # any in-range gridmodel bus, to isolate the element-id check
        calls = []
        # one-side "change_bus_<fam>(el_id, new_bus)"
        for fam in ["load", "gen", "sgen", "storage", "shunt", "svc"]:
            for bad in (-1, BIG):
                calls.append((f"change_bus_{fam}({bad})",
                              lambda fam=fam, bad=bad: getattr(self.grid, f"change_bus_{fam}")(bad, valid_bus)))
        # two-side "change_bus{1,2}_<fam>(el_id, new_bus)"
        for fam in ["powerline", "trafo"]:
            for side in (1, 2):
                for bad in (-1, BIG):
                    calls.append((f"change_bus{side}_{fam}({bad})",
                                  lambda fam=fam, side=side, bad=bad:
                                  getattr(self.grid, f"change_bus{side}_{fam}")(bad, valid_bus)))
        self._assert_raises_each(IndexError, calls)

    # ------------------------------------------------------------------
    # out-of-bounds WRITES
    # ------------------------------------------------------------------
    def test_update_slack_weights_by_id_out_of_bounds(self):
        calls = [
            ("update_slack_weights_by_id([BIG])",
             lambda: self.grid.update_slack_weights_by_id(np.array([BIG], dtype=np.int32))),
            ("update_slack_weights_by_id([-1])",
             lambda: self.grid.update_slack_weights_by_id(np.array([-1], dtype=np.int32))),
        ]
        self._assert_raises_each(IndexError, calls)

    def test_set_gen_regulated_bus_out_of_bounds(self):
        calls = [
            # bad gen_id (valid bus)
            ("set_gen_regulated_bus(-1, 0)", lambda: self.grid.set_gen_regulated_bus(-1, 0)),
            ("set_gen_regulated_bus(BIG, 0)", lambda: self.grid.set_gen_regulated_bus(BIG, 0)),
            # valid gen_id, bad bus_id
            ("set_gen_regulated_bus(0, -1)", lambda: self.grid.set_gen_regulated_bus(0, -1)),
            ("set_gen_regulated_bus(0, BIG)", lambda: self.grid.set_gen_regulated_bus(0, BIG)),
        ]
        self._assert_raises_each(IndexError, calls)

    def test_change_ratio_shift_trafo_out_of_bounds(self):
        calls = []
        for name, arg in [("change_ratio_trafo", 1.0),
                          ("change_shift_trafo", 0.1),
                          ("change_shift_trafo_deg", 1.0)]:
            for bad in (-1, BIG):
                calls.append((f"{name}({bad})",
                              lambda name=name, arg=arg, bad=bad: getattr(self.grid, name)(bad, arg)))
        self._assert_raises_each(IndexError, calls)

    # ------------------------------------------------------------------
    # out-of-bounds READ (getter)
    # ------------------------------------------------------------------
    def test_get_status_droop_hvdc_out_of_bounds(self):
        # case14 has no hvdc line; the empty-container bounds check still fires.
        calls = [
            ("get_status_droop_hvdc(-1)", lambda: self.grid.get_status_droop_hvdc(-1)),
            ("get_status_droop_hvdc(BIG)", lambda: self.grid.get_status_droop_hvdc(BIG)),
        ]
        self._assert_raises_each(IndexError, calls)

    # ------------------------------------------------------------------
    # array length / shape mismatches
    # ------------------------------------------------------------------
    def test_update_continuous_values_length_mismatch(self):
        # for each grid2op fast-update method, `new_values` must have the same
        # length as `has_changed`; a shorter `new_values` used to over-read.
        specs = [
            ("update_gens_p", self.n_gen),
            ("update_gens_v", self.n_gen),
            ("update_sgens_p", self.n_sgen),
            ("update_loads_p", self.n_load),
            ("update_loads_q", self.n_load),
            ("update_storages_p", self.n_storage),
        ]
        calls = []
        for name, n in specs:
            # has_changed longer than new_values by one -> mismatch (works for n == 0 too)
            has_changed = np.zeros(n + 1, dtype=bool)
            has_changed[-1] = True
            new_values = np.zeros(n, dtype=np.float32)
            calls.append((f"{name}(len mismatch)",
                          lambda name=name, hc=has_changed, nv=new_values:
                          getattr(self.grid, name)(hc, nv)))
        self._assert_raises_each(RuntimeError, calls)

    def test_set_pos_topo_vect_and_subid_size_mismatch(self):
        specs = [
            ("set_gen_pos_topo_vect", self.n_gen),
            ("set_load_pos_topo_vect", self.n_load),
            ("set_storage_pos_topo_vect", self.n_storage),
            ("set_gen_to_subid", self.n_gen),
            ("set_load_to_subid", self.n_load),
            ("set_storage_to_subid", self.n_storage),
            ("set_shunt_to_subid", self.n_shunt),
        ]
        calls = []
        for name, n in specs:
            wrong = np.zeros(n + 1, dtype=np.int32)  # one element too many
            calls.append((f"{name}(size mismatch)",
                          lambda name=name, arr=wrong: getattr(self.grid, name)(arr)))
        self._assert_raises_each(RuntimeError, calls)

    @unittest.skipUnless(TimeSeriesCPP_AVAILABLE, "TimeSeriesCPP not available")
    def test_compute_Vs_shape_mismatch(self):
        ts = TimeSeriesCPP(self.grid)
        nb_steps = 3
        Vinit = np.ones(self.n_bus, dtype=complex)

        def mat(rows, cols):
            return np.zeros((rows, cols), dtype=np.float64)

        # correctly-shaped matrices (baseline used to build the mis-shaped variants)
        gen_p = mat(nb_steps, self.n_gen)
        sgen_p = mat(nb_steps, self.n_sgen)
        load_p = mat(nb_steps, self.n_load)
        load_q = mat(nb_steps, self.n_load)

        calls = [
            # wrong number of columns (fewer gens than the grid has)
            ("gen_p wrong #cols",
             lambda: ts.compute_Vs(mat(nb_steps, self.n_gen - 1), sgen_p, load_p, load_q, Vinit, 10, 1e-8)),
            ("load_p wrong #cols",
             lambda: ts.compute_Vs(gen_p, sgen_p, mat(nb_steps, self.n_load - 1), load_q, Vinit, 10, 1e-8)),
            # mismatched number of rows (time steps) between matrices
            ("load_q wrong #rows",
             lambda: ts.compute_Vs(gen_p, sgen_p, load_p, mat(nb_steps - 1, self.n_load), Vinit, 10, 1e-8)),
        ]
        self._assert_raises_each(RuntimeError, calls)

    # ------------------------------------------------------------------
    # regression: the same methods still work for valid inputs
    # ------------------------------------------------------------------
    def test_valid_calls_still_work(self):
        # these must NOT raise
        self.grid.deactivate_load(0)
        self.grid.reactivate_load(0)
        self.grid.change_ratio_trafo(0, 1.0)
        self.grid.change_shift_trafo(0, 0.0)
        self.grid.set_gen_regulated_bus(0, 0)

        has_changed = np.zeros(self.n_gen, dtype=bool)
        new_values = np.zeros(self.n_gen, dtype=np.float32)
        self.grid.update_gens_p(has_changed, new_values)

        if TimeSeriesCPP_AVAILABLE:
            ts = TimeSeriesCPP(self.grid)
            nb_steps = 3
            Vinit = np.ones(self.n_bus, dtype=complex)
            res = ts.compute_Vs(np.zeros((nb_steps, self.n_gen), dtype=np.float64),
                                 np.zeros((nb_steps, self.n_sgen), dtype=np.float64),
                                 np.zeros((nb_steps, self.n_load), dtype=np.float64),
                                 np.zeros((nb_steps, self.n_load), dtype=np.float64),
                                 Vinit, 10, 1e-8)
            self.assertIsNotNone(res)


if __name__ == "__main__":
    unittest.main()
