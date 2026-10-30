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

import gc
import os
import sys
import unittest
import warnings
import numpy as np

import pandapower.networks as pn
from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import LSGrid, RegulationMode
from lightsim2grid.algorithm import PhaseControl

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
        # gen 0 sits on bus 1; bus 9 is a PQ load bus, so this is a voltage-control
        # configuration v1 actually supports.
        #
        # Bus 0 would ALSO be an in-range id -- which is all this method itself
        # bounds-checks -- but it is case14's slack (the pandapower ext_grid), so it
        # owns no Vm unknown and LSGrid::fill_voltage_control_solver_data rejects the
        # grid on its `is_pq[reg_solver]` check. Note the criterion is whether the
        # regulated bus' magnitude is already DETERMINED (slack, or a PV bus pinned by
        # a local regulator), not whether a generator happens to sit there: several
        # generators may perfectly well share one regulated bus, each with its own Q
        # unknown plus a sharing row. Either way it is not the argument check under
        # test here, and it makes the compute_Vs call below raise -- which used to slip
        # through only because the batch algorithms silently dropped voltage control.
        self.grid.set_gen_regulated_bus(0, 9)

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


class TestConfigGettersOwnTheirData(unittest.TestCase):
    """The reactive ranges and raw target Q a converter sets on a grid, read back: arrays that
    stay valid once the grid is gone, or once a setter gave it another size. Run under ASan
    (sanitizers.yml), a view into the grid's buffer is a heap-use-after-free here."""

    @staticmethod
    def _grid():
        """case14 with an hvdc line: without one, the station ranges are empty arrays, which
        own their (empty) buffer whatever the getter returns, and check nothing"""
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _aux_make_hvdc import make_case14_hvdc
        _, grid = make_case14_hvdc()
        return grid

    def test_arrays_outlive_the_grid(self):
        grid = self._grid()
        n_gen = len(grid.get_generators())
        n_hvdc = len(grid.get_dclines())
        self.assertGreater(n_hvdc, 0)
        grid.set_gen_reactive_range_max(np.arange(n_gen, dtype=float))
        grid.set_gen_raw_target_q(np.full(n_gen, 3.))
        grid.set_station_reactive_range_max(np.full(n_hvdc, 7.), np.full(n_hvdc, 9.))
        arrays = [grid.get_gen_reactive_range_max(), grid.get_gen_raw_target_q(),
                  grid.get_station_reactive_range_max(1), grid.get_station_reactive_range_max(2)]
        for arr in arrays:
            # its own buffer, or a reference to whatever holds it
            self.assertTrue(arr.flags.owndata or arr.base is not None)
        del grid
        gc.collect()
        np.testing.assert_array_equal(arrays[0], np.arange(n_gen, dtype=float))
        np.testing.assert_array_equal(arrays[1], np.full(n_gen, 3.))
        np.testing.assert_array_equal(arrays[2], np.full(n_hvdc, 7.))
        np.testing.assert_array_equal(arrays[3], np.full(n_hvdc, 9.))

    def test_arrays_outlive_a_resize(self):
        grid = self._grid()
        n_gen = len(grid.get_generators())
        n_hvdc = len(grid.get_dclines())
        self.assertGreater(n_hvdc, 0)
        grid.set_gen_reactive_range_max(np.arange(n_gen, dtype=float))
        grid.set_gen_raw_target_q(np.full(n_gen, 3.))
        grid.set_station_reactive_range_max(np.full(n_hvdc, 7.), np.full(n_hvdc, 9.))
        ranges = grid.get_gen_reactive_range_max()
        raw_q = grid.get_gen_raw_target_q()
        station1 = grid.get_station_reactive_range_max(1)
        station2 = grid.get_station_reactive_range_max(2)
        # "none": the grid's buffers are freed
        grid.set_gen_reactive_range_max(np.zeros(0))
        grid.set_gen_raw_target_q(np.zeros(0))
        grid.set_station_reactive_range_max(np.zeros(0), np.zeros(0))
        np.testing.assert_array_equal(ranges, np.arange(n_gen, dtype=float))
        np.testing.assert_array_equal(raw_q, np.full(n_gen, 3.))
        np.testing.assert_array_equal(station1, np.full(n_hvdc, 7.))
        np.testing.assert_array_equal(station2, np.full(n_hvdc, 9.))


class TestOuterTapsOfReplacedTransformers(unittest.TestCase):
    """The checks after an NROuter_* solve read the taps its loops left, one per transformer of
    that solve. ``init_trafo`` then replaces the transformers: the checks must refuse those taps,
    not read them in the tap tables it freed (a heap-use-after-free under ASan, then an id check
    raising after the read: the message shows the taps are refused before any read)."""

    @staticmethod
    def _grid():
        """a triangle 0-1-2: lines 0-1 and 1-2, a phase shifter 0-2 (-10 .. 10 degrees, 21 taps)
        regulating its active power, loads at 1 and 2, the slack generator at 0 -- the grid of
        src/tests/test_phase_control_loop.cpp"""
        grid = LSGrid()
        grid.set_sn_mva(100.)
        grid.set_init_vm_pu(1.)
        grid.init_bus(3, 1, np.full(3, 138.), 0, 0)
        grid.init_powerlines(np.full(2, .01), np.full(2, .08), np.zeros(2, complex), np.array([0, 1]), np.array([1, 2]))
        grid.init_trafo(np.array([.005]), np.array([.06]), np.zeros(1, complex), np.ones(1), np.zeros(1),
                        [False], np.array([0]), np.array([2]), False)
        grid.init_loads(np.array([60., 90.]), np.array([20., 30.]), np.array([1, 2]))
        grid.init_generators(np.zeros(1), np.full(1, 1.02), np.full(1, -1000.), np.full(1, 1000.), np.array([0]))
        grid.add_gen_slackbus(0, 1.)
        alpha = [-10. + k for k in range(21)]
        zero = [0.] * 21
        grid.set_trafo_phase_tap_changer(0, 0, 10, [1.] * 21, alpha, zero, zero, zero, zero)
        return grid

    def test_checks_after_init_trafo_raise(self):
        ref = self._grid()
        ref.change_algorithm("NRSing_SparseLU")
        self.assertEqual(ref.ac_pf(np.ones(3, dtype=complex), 30, 1e-12).size, 3)
        p0 = ref.get_trafo_res1()[0][0]

        grid = self._grid()
        grid.set_trafo_phase_tap_regulation(0, RegulationMode.ACTIVE_POWER, True, p0 + 5., 0., 1)
        grid.change_algorithm("NROuter_SparseLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(PhaseControl())
        self.assertEqual(grid.ac_pf(np.ones(3, dtype=complex), 30, 1e-12).size, 3)
        grid.get_violations()
        grid.get_physical_violations()

        grid.init_trafo(np.zeros(0), np.zeros(0), np.zeros(0, complex), np.zeros(0), np.zeros(0),
                        [], np.zeros(0, dtype=int), np.zeros(0, dtype=int), False)
        with self.assertRaisesRegex(IndexError, "were the transformers replaced"):
            grid.get_violations()
        with self.assertRaisesRegex(IndexError, "were the transformers replaced"):
            grid.get_physical_violations()


if __name__ == "__main__":
    unittest.main()
