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

The classes at the end check the arrays read back from a grid or an algorithm instead: none
may read memory its owner has freed.
"""

import gc
import os
import re
import sys
import unittest
import warnings
import weakref
import numpy as np
import scipy.sparse as sp

import pandapower as pp
import pandapower.networks as pn
from lightsim2grid.network import init_from_pandapower
from lightsim2grid.lightsim2grid_cpp import LSGrid, RegulationMode, NR_SparseLU
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


def _solve(grid, ac=True):
    V0 = np.ones(grid.total_bus(), dtype=complex)
    V = grid.ac_pf(V0, 10, 1e-8) if ac else grid.dc_pf(V0, 10, 1e-8)
    if V.size == 0:
        raise AssertionError(f"the {'AC' if ac else 'DC'} powerflow did not converge")


def _solved_grid():
    """case14 with one of each element it lacks (a static generator, a storage unit, an hvdc
    line), solved in AC, in DC, then in AC again: every result getter has something to show,
    and the AC and the DC solver both have their caches filled"""
    net = pn.case14()
    pp.create_sgen(net, 4, p_mw=5., q_mvar=1.)
    pp.create_storage(net, 5, p_mw=3., max_e_mwh=10.)
    pp.create_dcline(net, 3, 9, p_mw=20., loss_percent=1., loss_mw=.5, vm_from_pu=1.01, vm_to_pu=1.01)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        grid = init_from_pandapower(net)
    _solve(grid, ac=True)
    _solve(grid, ac=False)
    _solve(grid, ac=True)
    return grid


def _arrays_of(res):
    """the numpy arrays in what a getter returned: the array itself, those of a tuple, or those
    of a scipy sparse matrix"""
    if isinstance(res, np.ndarray):
        return [res]
    if isinstance(res, tuple):
        return [arr for arr in res if isinstance(arr, np.ndarray)]
    if sp.issparse(res):
        return [res.data, res.indices, res.indptr]
    return []


def _memory_holder(arr):
    """What keeps the memory of `arr` alive: the array itself when it owns its data, else the
    object its `.base` chain ends on (a capsule, the grid...). None for a view whose chain ends
    on nothing: memory that nothing keeps alive."""
    while isinstance(arr, np.ndarray):
        if arr.flags.owndata:
            return arr
        if arr.base is None:
            return None
        arr = arr.base
    return arr


def _is_snapshot(arr):
    """memory of its own: held by the array, or by the capsule pybind11 hands a result over in"""
    holder = _memory_holder(arr)
    return isinstance(holder, np.ndarray) or type(holder).__name__ == "PyCapsule"


def _array_getters(obj):
    """{name: what it returned} for the get_* / id_* methods of `obj` that take no argument and
    return arrays"""
    res = {}
    for name in dir(obj):
        if not name.startswith(("get_", "id_")):
            continue
        try:
            got = getattr(obj, name)()
        except Exception:  # takes arguments, or does not apply to this grid or algorithm
            continue
        if _arrays_of(got):
            res[name] = got
    return res


_CONTAINERS = ("get_generators", "get_loads", "get_static_generators", "get_storages", "get_shunts",
               "get_lines", "get_trafos", "get_dclines")


def _grid_array_getters(grid):
    """_array_getters of the grid and of its element containers (as ``get_lines().get_bus_id_side_1``)"""
    res = _array_getters(grid)
    for cont in _CONTAINERS:
        for name, got in _array_getters(getattr(grid, cont)()).items():
            res[f"{cont}().{name}"] = got
    return res


def _data_ptr(arr):
    return arr.__array_interface__["data"][0]


def _reuse_freed_memory():
    """allocate blocks of every small size and fill them: freed memory read afterwards no longer
    shows what it held"""
    junk = [np.full(n, -1.2345e300) for n in range(1, 2000)]
    del junk


def _standalone_solver(grid):
    """an algorithm solving the grid's AC problem on its own, as pandapower_compat uses one"""
    solver = NR_SparseLU()
    conv = solver.compute_pf(grid.get_Ybus_solver().copy(),
                             np.abs(grid.get_V_solver()).astype(complex),
                             1.0 * grid.get_Sbus_solver(),
                             np.array(grid.get_slack_ids_solver(), dtype=np.int32),
                             1.0 * grid.get_slack_weights_solver(),
                             np.array(grid.get_pv_solver(), dtype=np.int32),
                             np.array(grid.get_pq_solver(), dtype=np.int32),
                             10, 1e-8)
    if not conv:
        raise AssertionError("the algorithm alone did not converge")
    return solver


class TestResultViewsKeepTheirGrid(unittest.TestCase):
    """LSGrid's result getters (get_gen_res, get_line_res1_full, get_load_theta,
    get_gen_target_p...) return views of the grid's own buffers, which show its next powerflow.
    A view must keep its grid alive, and the grid must keep the buffers it handed out: the views
    did not keep the grid, and a diverged powerflow replaced every result buffer. Either way an
    array read freed memory afterwards (a heap-use-after-free under ASan, sanitizers.yml)."""

    maxDiff = None  # the lists of getters at fault, in full

    def test_every_view_keeps_its_grid_alive(self):
        """through `.base`, where numpy (and the reader) can see what holds the memory"""
        grid = _solved_grid()
        got = {nm: res for nm, res in _grid_array_getters(grid).items() if not nm.endswith("_solver")}
        # the sweep found the result getters, each with something to show
        for nm in ("get_gen_res", "get_line_res1_full", "get_trafo_res2", "get_sgens_res",
                   "get_storages_res_full", "get_dcline_res1_full", "get_load_theta",
                   "get_gen_target_p", "get_all_shunt_buses", "get_lines().get_yac_eff_11",
                   "get_dclines().get_bus_id_side_1"):
            self.assertIn(nm, got)
            self.assertTrue(all(arr.size > 0 for arr in _arrays_of(got[nm])), nm)
        dangling = sorted(nm for nm, res in got.items()
                          if any(_memory_holder(arr) is None for arr in _arrays_of(res)))
        self.assertEqual(dangling, [], "these views do not keep their grid alive")

    def test_a_view_outlives_the_function_that_made_its_grid(self):
        """the pattern of the helper that made test_outer_redistribute_slack.py fail on one CI
        runner only: the grid is local, the array is returned"""
        def gen_p():
            grid = _solved_grid()
            p = grid.get_gen_res()[0]
            return p, p.copy(), weakref.ref(grid)

        p, expected, grid_alive = gen_p()
        gc.collect()
        self.assertIsNotNone(grid_alive(), "the array does not keep its grid alive")
        _reuse_freed_memory()
        np.testing.assert_array_equal(p, expected)
        # and only the array keeps it: the grid goes with it
        del p
        gc.collect()
        self.assertIsNone(grid_alive())

    def test_a_control_keeps_its_grid_alive(self):
        grid = _solved_grid()
        controls = [grid.get_ac_algo_controler(), grid.get_dc_algo_controler()]
        grid_alive = weakref.ref(grid)
        del grid
        gc.collect()
        self.assertIsNotNone(grid_alive(), "a control does not keep its grid alive")
        for ctrl in controls:
            self.assertIn(ctrl.need_reset_solver(), (True, False))
        del controls, ctrl
        gc.collect()
        self.assertIsNone(grid_alive())

    def test_views_stay_the_grid_buffers(self):
        """After any ordinary call -- powerflows, a change of injection, of topology, one bus
        less, a divergence, a reset, another algorithm -- what a view shows is still a buffer
        the grid holds, of the same size. Only pointers are compared: nothing reads a buffer
        that may be gone. The `_internal` views of the solver side are left out: they promise
        nothing of the kind (TestInternalTwinsAreViews)."""
        grid = _solved_grid()
        held = {}
        for nm, res in _grid_array_getters(grid).items():
            if nm.endswith("_internal"):
                continue
            for idx, arr in enumerate(_arrays_of(res)):
                if not _is_snapshot(arr):
                    held[nm, idx] = (_data_ptr(arr), arr.shape)
        self.assertIn(("get_gen_res", 0), held)
        self.assertIn(("get_trafos().get_yac_eff_22", 0), held)

        V0 = np.ones(grid.total_bus(), dtype=complex)
        load_p = grid.get_loads()[0].target_p_mw
        trafo_7 = [tr.id for tr in grid.get_trafos() if {tr.bus1_id, tr.bus2_id} == {6, 7}]
        gen_7 = [gen.id for gen in grid.get_generators() if gen.bus_id == 7]
        gen_2 = [gen.id for gen in grid.get_generators() if gen.bus_id == 2]
        self.assertEqual((len(trafo_7), len(gen_7), len(gen_2)), (1, 1, 1))
        nb_bus = grid.nb_connected_bus()

        def one_bus_less():
            # bus 7 hangs on one transformer and holds one generator: both off, it is gone
            grid.deactivate_trafo(trafo_7[0])
            grid.deactivate_gen(gen_7[0])

        def divergence():
            grid.change_p_load(0, 1e5)
            self.assertEqual(grid.ac_pf(V0, 10, 1e-8).size, 0)
            grid.change_p_load(0, load_p)

        steps = [("powerflows", lambda: None),
                 ("a change of load", lambda: grid.change_p_load(0, 1.1 * load_p)),
                 ("a line disconnected", lambda: grid.deactivate_powerline(5)),
                 ("a generator disconnected", lambda: grid.deactivate_gen(gen_2[0])),
                 ("one bus less", one_bus_less),
                 ("a divergence", divergence),
                 ("a reset", grid.tell_solver_need_reset),
                 ("another algorithm", lambda: grid.change_algorithm("NR_SparseLU"))]
        moved = {}
        for label, step in steps:
            step()
            _solve(grid, ac=True)
            _solve(grid, ac=False)
            fresh = _grid_array_getters(grid)
            for (nm, idx), (ptr, shape) in held.items():
                arr = _arrays_of(fresh[nm])[idx] if nm in fresh else None
                if arr is None or _data_ptr(arr) != ptr or arr.shape != shape:
                    moved.setdefault(label, []).append(nm)
                    if arr is not None:
                        held[nm, idx] = (_data_ptr(arr), arr.shape)  # report each step's own
        self.assertEqual(grid.nb_connected_bus(), nb_bus - 1)
        self.assertEqual({label: sorted(set(nms)) for label, nms in moved.items()}, {},
                         "buffers the grid replaced while a view of them was out")

    def test_a_divergence_resets_the_views_in_place(self):
        """a diverged powerflow leaves no result: NaN, rather than the previous solve's values"""
        grid = _solved_grid()
        V0 = np.ones(grid.total_bus(), dtype=complex)
        results = {nm: _arrays_of(res) for nm, res in _array_getters(grid).items()
                   if re.search(r"_res\d?(_full)?$|_theta\d?$", nm)}
        for nm in ("get_gen_res", "get_line_res1_full", "get_dcline_res2_full", "get_trafo_theta1"):
            self.assertIn(nm, results)
        load_p = grid.get_loads()[0].target_p_mw
        grid.change_p_load(0, 1e5)
        self.assertEqual(grid.ac_pf(V0, 10, 1e-8).size, 0)
        not_reset = sorted(nm for nm, arrs in results.items()
                           if not all(np.isnan(arr).all() for arr in arrs))
        self.assertEqual(not_reset, [])
        # the next powerflow writes into the same arrays
        grid.change_p_load(0, load_p)
        _solve(grid)
        for nm, arrs in results.items():
            for arr, fresh in zip(arrs, _arrays_of(getattr(grid, nm)())):
                np.testing.assert_array_equal(arr, fresh)
                self.assertTrue(np.isfinite(arr).all(), nm)


class TestSolverGettersAreSnapshots(unittest.TestCase):
    """What the solver side holds -- its bus labelling, pv / pq split, slack, Sbus, Ybus,
    voltages -- is rebuilt by later powerflows (a change of topology, a reset, another
    algorithm), which free the buffers it was in: a view of them would read freed memory,
    whatever keeps the grid alive. The `*_solver` getters of LSGrid and the voltages of an
    algorithm are copies, taken when called; their `_internal` twins are the views, at the
    caller's risk (TestInternalTwinsAreViews)."""

    maxDiff = None  # the lists of getters at fault, in full

    def test_solver_getters_own_their_data(self):
        grid = _solved_grid()
        got = {nm: res for nm, res in _array_getters(grid).items() if nm.endswith("_solver")}
        for nm in ("get_V_solver", "get_Va_solver", "get_Vm_solver", "get_Sbus_solver",
                   "get_dcSbus_solver", "get_pv_solver", "get_pq_solver", "get_slack_ids_solver",
                   "get_slack_weights_solver", "get_Ybus_solver", "get_dcYbus_solver"):
            self.assertIn(nm, got)
        views = sorted(nm for nm, res in got.items()
                       if not all(_is_snapshot(arr) for arr in _arrays_of(res)))
        self.assertEqual(views, [], "these solver getters return views")

    def test_solver_snapshots_keep_their_values(self):
        grid = _solved_grid()
        got = {nm: _arrays_of(res) for nm, res in _array_getters(grid).items() if nm.endswith("_solver")}
        expected = {nm: [arr.copy() for arr in arrs] for nm, arrs in got.items()}
        grid.deactivate_powerline(5)
        grid.deactivate_gen([gen.id for gen in grid.get_generators() if gen.bus_id == 2][0])
        grid.change_p_load(0, 1.1 * grid.get_loads()[0].target_p_mw)
        grid.tell_solver_need_reset()
        _solve(grid, ac=True)
        _solve(grid, ac=False)
        # the solver side did change
        self.assertFalse(np.array_equal(grid.get_Sbus_solver(), expected["get_Sbus_solver"][0]))
        changed = sorted(nm for nm, arrs in got.items()
                         if not all(np.array_equal(arr, exp, equal_nan=True)
                                    for arr, exp in zip(arrs, expected[nm])))
        self.assertEqual(changed, [], "these solver getters follow the grid")

    def test_algorithm_getters_own_their_data(self):
        grid = _solved_grid()
        for label, algo in (("get_algo()", grid.get_algo()), ("get_dc_algo()", grid.get_dc_algo()),
                            ("NR_SparseLU", _standalone_solver(grid))):
            got = _array_getters(algo)
            self.assertIn("get_V", got, label)
            views = sorted(nm for nm, res in got.items()
                           if not nm.endswith("_internal")
                           and not all(_is_snapshot(arr) for arr in _arrays_of(res)))
            self.assertEqual(views, [], f"{label}: these getters return views")

    def test_voltages_outlive_their_algorithm(self):
        """the pandapower_compat pattern: an algorithm made for one call, its voltages returned"""
        grid = _solved_grid()

        def solve_alone():
            solver = _standalone_solver(grid)
            got = solver.get_Va(), solver.get_Vm(), solver.get_V()
            return got, [arr.copy() for arr in got]

        got, expected = solve_alone()
        gc.collect()
        self.assertTrue(all(_is_snapshot(arr) for arr in got), "the voltages are a view of the algorithm")
        _reuse_freed_memory()
        for arr, exp in zip(got, expected):
            np.testing.assert_array_equal(arr, exp)


# the solver getters that copy, each with a twin `<name>_internal` that does not
_GRID_INTERNAL_TWINS = ("get_V_solver", "get_Va_solver", "get_Vm_solver", "get_Sbus_solver",
                        "get_dcSbus_solver", "get_Ybus_solver", "get_dcYbus_solver",
                        "get_pv_solver", "get_pq_solver", "get_ac_pv_solver", "get_dc_pv_solver",
                        "get_ac_pq_solver", "get_dc_pq_solver", "get_slack_ids_solver",
                        "get_slack_ids_dc_solver", "get_slack_weights_solver",
                        "get_ac_slack_weights_solver", "get_dc_slack_weights_solver")
_ALGO_INTERNAL_TWINS = ("get_V", "get_Va", "get_Vm")
# the Jacobian's, on a Newton-Raphson only
_GRID_JACOBIAN_TWIN = "get_J_solver"
_ALGO_JACOBIAN_TWIN = "get_J"


class TestInternalTwinsAreViews(unittest.TestCase):
    """Each solver getter that copies has a twin, `<name>_internal`, that does not: a read only
    view of the buffer the solver holds, for code that measured the copy. It keeps its owner
    alive, and is valid only until a later call rebuilds that buffer, which nothing checks: at
    the caller's risk (docs/security.rst). Hence the sweeps above leave the twins out."""

    maxDiff = None

    def _check_twins(self, owner, names):
        nb_checked = 0
        for nm in names:
            with self.subTest(getter=nm):
                self.assertTrue(hasattr(owner, nm + "_internal"), f"{nm} has no _internal twin")
                safe = getattr(owner, nm)()
                view = getattr(owner, nm + "_internal")()
                again = getattr(owner, nm + "_internal")()
                self.assertIs(type(view), type(safe))
                self.assertEqual(len(_arrays_of(view)), len(_arrays_of(safe)))
                for arr_safe, arr, arr_again in zip(_arrays_of(safe), _arrays_of(view), _arrays_of(again)):
                    # what the copy shows...
                    np.testing.assert_array_equal(arr, arr_safe)
                    if arr.size == 0:
                        continue  # numpy gives an empty array a buffer of its own
                    nb_checked += 1
                    # ... read in the owner's own buffer, call after call, and read only
                    self.assertIs(_memory_holder(arr), owner)
                    self.assertTrue(np.shares_memory(arr, arr_again))
                    self.assertFalse(arr.flags.writeable)
        self.assertGreater(nb_checked, len(names) // 2)

    def test_grid_solver_getters(self):
        grid = _solved_grid()
        self._check_twins(grid, _GRID_INTERNAL_TWINS)

    def test_algorithm_voltages(self):
        grid = _solved_grid()
        for label, algo in (("get_algo()", grid.get_algo()), ("get_dc_algo()", grid.get_dc_algo()),
                            ("NR_SparseLU", _standalone_solver(grid))):
            with self.subTest(algo=label):
                self._check_twins(algo, _ALGO_INTERNAL_TWINS)

    def test_jacobian(self):
        """the Jacobian a Newton-Raphson holds, through the grid, its algorithm, an algorithm
        alone, and an NROuter_* (which forwards to the Newton-Raphson inside it)"""
        grid = _solved_grid()
        self._check_twins(grid, [_GRID_JACOBIAN_TWIN])
        with self.subTest(algo="get_algo()"):
            self._check_twins(grid.get_algo(), [_ALGO_JACOBIAN_TWIN])
        with self.subTest(algo="NR_SparseLU"):
            self._check_twins(_standalone_solver(grid), [_ALGO_JACOBIAN_TWIN])
        with self.subTest(algo="NROuter_SparseLU"):
            outer = _solved_grid()
            outer.change_algorithm("NROuter_SparseLU")
            outer.clear_outer_loops()
            _solve(outer)
            self._check_twins(outer, [_GRID_JACOBIAN_TWIN])
            self._check_twins(outer.get_algo(), [_ALGO_JACOBIAN_TWIN])

    def test_every_internal_method_has_its_twin(self):
        """no `_internal` name but those: each is the view of a getter that copies"""
        grid = _solved_grid()
        for label, owner, twins in (("LSGrid", grid, _GRID_INTERNAL_TWINS + (_GRID_JACOBIAN_TWIN,)),
                                    ("get_algo()", grid.get_algo(), _ALGO_INTERNAL_TWINS + (_ALGO_JACOBIAN_TWIN,)),
                                    ("NR_SparseLU", NR_SparseLU(), _ALGO_INTERNAL_TWINS + (_ALGO_JACOBIAN_TWIN,))):
            with self.subTest(owner=label):
                internal = sorted(nm for nm in dir(owner) if nm.endswith("_internal") and not nm.startswith("_"))
                self.assertEqual(internal, sorted(nm + "_internal" for nm in twins))

    def test_a_view_keeps_its_grid_alive(self):
        def view_of_a_local_grid():
            grid = _solved_grid()
            return grid.get_V_solver_internal(), grid.get_V_solver(), weakref.ref(grid)

        view, expected, grid_alive = view_of_a_local_grid()
        gc.collect()
        self.assertIsNotNone(grid_alive(), "the view does not keep its grid alive")
        np.testing.assert_array_equal(view, expected)
        del view
        gc.collect()
        self.assertIsNone(grid_alive())


if __name__ == "__main__":
    unittest.main()
