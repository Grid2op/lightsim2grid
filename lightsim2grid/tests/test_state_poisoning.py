# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Python-visible behaviour of the "poisoned state" rejections.

A pickle (or a binary file: both go through ``LSGrid::set_state``) is only
length-checked while it is being read. Whatever it declares afterwards is copied
straight into the C++ containers, and release wheels are built ``-O3 -DNDEBUG``,
so an inconsistency that nothing rejects becomes an out-of-bounds read or write
on the next powerflow rather than an error. Every case below segfaulted the
interpreter before its fix.

The deep ``set_state`` cases are covered exhaustively (and under
valgrind/ASan) by the C++ suite, ``src/tests/test_state_poisoning.cpp``; what is
checked here is that the same poisoned states raise a clean Python exception
through the real ``pickle.loads`` path, and that the checks reachable without any
serialization at all (an out-of-range SVC regulated bus set by a grid loader)
raise too.
"""

import pickle
import unittest

import numpy as np

from lightsim2grid.lightsim2grid_cpp import LSGrid, PandaPowerConverter, ContingencyAnalysisCPP
from lightsim2grid.tests._exotic_elements_case_ls import build_exotic_elements_case_grid


# state indices, mirroring LSGrid::StateRes (see LSGrid.hpp)
LINE_ID = 8
SHUNT_ID = 9
TRAFO_ID = 10
HVDC_ID = 15
SVC_ID = 16

BIG = 10 ** 7


def _rebuild(state):
    """Reconstruct an LSGrid from a (possibly tampered) pickle state, exactly the
    way ``pickle.loads`` does -- ``__setstate__`` on an already-built object is a
    no-op for a pybind11 ``py::pickle`` type."""
    obj = LSGrid.__new__(LSGrid)
    obj.__setstate__(state)
    return obj


def _set_at(state, path, value):
    """Return a copy of the nested state tuple with `path` replaced by `value`."""
    if not path:
        return value
    head, rest = path[0], path[1:]
    as_list = list(state)
    as_list[head] = _set_at(as_list[head], rest, value)
    return tuple(as_list)


def _get_at(state, path):
    for i in path:
        state = state[i]
    return state


class TestPoisonedState(unittest.TestCase):
    """The grid carries a phase-shifting transformer, an SVC and 3 HVDC lines, so
    every container the cases below poison is non-empty."""

    def setUp(self):
        self.grid = build_exotic_elements_case_grid()
        self.state = self.grid.__getstate__()
        # (version_major, version_medium, version_minor, LSGrid::StateRes)
        self.inner = (3,)

    def _poison(self, path, value):
        return _set_at(self.state, self.inner + tuple(path), value)

    # ------------------------------------------------------------------
    # sanity: the untouched state round-trips
    # ------------------------------------------------------------------
    def test_valid_state_round_trips(self):
        grid = pickle.loads(pickle.dumps(self.grid))
        self.assertIsNone(grid.check_grid())
        v0 = np.ones(grid.total_bus(), dtype=complex)
        self.assertNotEqual(len(grid.ac_pf(v0, 20, 1e-7)), 0)

    # ------------------------------------------------------------------
    # TrafoContainer: the (alpha -> r/x correction) tables. They are indexed by
    # transformer id inside set_state itself (_update_model_coeffs), i.e. before
    # check_grid() gets a chance to run.
    # ------------------------------------------------------------------
    def test_rx_corr_tables_missing_is_rejected(self):
        # shift-dependent r/x enabled (it is, in this grid) but no table at all
        self.assertTrue(_get_at(self.state, (3, TRAFO_ID, 5)))
        bad = self._poison((TRAFO_ID, 8), ())
        bad = _set_at(bad, (3, TRAFO_ID, 9), ())
        with self.assertRaises(RuntimeError):
            _rebuild(bad)

    def test_rx_corr_tables_too_short_is_rejected(self):
        alpha = _get_at(self.state, (3, TRAFO_ID, 8))
        self.assertGreater(len(alpha), 1)
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((TRAFO_ID, 8), tuple(alpha)[:-1]))

    def test_rx_corr_alpha_and_pct_of_different_lengths_is_rejected(self):
        alpha = list(_get_at(self.state, (3, TRAFO_ID, 8)))
        pct = list(_get_at(self.state, (3, TRAFO_ID, 9)))
        alpha[0] = (-0.5, 0.0, 0.5)
        pct[0] = (1.0, 2.0)  # one correction short
        bad = self._poison((TRAFO_ID, 8), tuple(alpha))
        bad = _set_at(bad, (3, TRAFO_ID, 9), tuple(pct))
        with self.assertRaises(RuntimeError):
            _rebuild(bad)

    # ------------------------------------------------------------------
    # TwoSidesContainer: status_global_. `nb()` is side_1_.nb(), so nothing tied
    # status_global_'s length to the number of elements indexing it.
    # ------------------------------------------------------------------
    def test_empty_status_global_on_lines_is_rejected(self):
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((LINE_ID, 0, 0, 3), ()))

    def test_short_status_global_on_lines_is_rejected(self):
        status = _get_at(self.state, (3, LINE_ID, 0, 0, 3))
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((LINE_ID, 0, 0, 3), tuple(status)[:-1]))

    def test_wrong_status_global_on_trafos_is_rejected(self):
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((TRAFO_ID, 0, 0, 3), ()))

    def test_wrong_status_global_on_hvdc_is_rejected(self):
        status = _get_at(self.state, (3, HVDC_ID, 0, 3))
        with self.assertRaises(RuntimeError):
            _rebuild(self._poison((HVDC_ID, 0, 3), tuple(status) + (True,)))

    # ------------------------------------------------------------------
    # SvcContainer: the regulated bus id, used as an index by the powerflow.
    # ------------------------------------------------------------------
    def test_out_of_range_svc_regulated_bus_is_rejected(self):
        reg = _get_at(self.state, (3, SVC_ID, 6))
        self.assertGreater(len(reg), 0)
        with self.assertRaises(IndexError):
            _rebuild(self._poison((SVC_ID, 6), (BIG,) * len(reg)))

    def test_negative_svc_regulated_bus_is_rejected(self):
        reg = _get_at(self.state, (3, SVC_ID, 6))
        with self.assertRaises(IndexError):
            _rebuild(self._poison((SVC_ID, 6), (-5,) * len(reg)))

    def test_minus_one_svc_regulated_bus_is_accepted(self):
        reg = _get_at(self.state, (3, SVC_ID, 6))
        grid = _rebuild(self._poison((SVC_ID, 6), (-1,) * len(reg)))
        self.assertIsNone(grid.check_grid())

    def test_out_of_range_svc_regulated_bus_from_init_svcs_is_rejected(self):
        # the same hole with no serialization involved: this is what
        # init_from_pypowsybl does with the bus ids it reads from the source grid
        grid = build_exotic_elements_case_grid()
        grid.init_svcs([1],
                       np.array([1.0]), np.array([0.0]), np.array([0.0]),
                       np.array([-0.03]), np.array([0.03]),
                       np.array([BIG], dtype=np.int32),   # regulated bus, out of range
                       np.array([5], dtype=np.int32))
        with self.assertRaises(IndexError):
            grid.check_grid()

    # ------------------------------------------------------------------
    # pos_topo_vect on an element that is not in the grid2op topology vector.
    # check_grid() proved the collected positions were a permutation of [0, K)
    # with K counting shunts / sgens / svcs / hvdc too, while update_topo() sizes
    # its caller arrays without them -- so a validated load position could still
    # be past the end of the array update_topo() indexes with it.
    # ------------------------------------------------------------------
    def test_pos_topo_vect_on_a_shunt_is_rejected(self):
        nb_shunt = len(_get_at(self.state, (3, SHUNT_ID, 0, 0, 1)))
        self.assertGreater(nb_shunt, 0)
        bad = self._poison((SHUNT_ID, 0, 0, 5), True)
        bad = _set_at(bad, (3, SHUNT_ID, 0, 0, 6), tuple(range(nb_shunt)))
        with self.assertRaises(RuntimeError):
            _rebuild(bad)

    def test_update_topo_stays_in_bounds_with_a_full_topo_vect(self):
        # the legitimate arrangement the check above protects: positions only on
        # the containers update_topo() drives, forming a permutation of
        # [0, dim_topo).
        grid = build_exotic_elements_case_grid()
        n_load = len(grid.get_loads())
        n_gen = len(grid.get_generators())
        n_sto = len(grid.get_storages())
        n_line = len(grid.get_lines())
        n_trafo = len(grid.get_trafos())
        dim_topo = n_load + n_gen + n_sto + 2 * n_line + 2 * n_trafo
        pos = np.arange(dim_topo, dtype=np.int32)
        off = 0
        for setter, nb in ((grid.set_load_pos_topo_vect, n_load),
                           (grid.set_gen_pos_topo_vect, n_gen),
                           (grid.set_storage_pos_topo_vect, n_sto),
                           (grid.set_line_pos1_topo_vect, n_line),
                           (grid.set_line_pos2_topo_vect, n_line),
                           (grid.set_trafo_pos1_topo_vect, n_trafo),
                           (grid.set_trafo_pos2_topo_vect, n_trafo)):
            setter(pos[off:off + nb])
            off += nb
        self.assertIsNone(grid.check_grid())
        grid.update_topo(np.zeros(dim_topo, dtype=bool), np.ones(dim_topo, dtype=np.int32))


class TestDegenerateGridScalars(unittest.TestCase):
    """``sn_mva`` is the base power of the whole per-unit system and ``init_vm_pu``
    the flat-start voltage magnitude. Neither was validated anywhere, and a
    degenerate value does not make the powerflow fail -- it makes it quietly
    wrong. With ``sn_mva = nan`` the DC powerflow used to report CONVERGENCE and
    return NaN branch flows: the built-in solvers' own finiteness guards see a
    perfectly finite ``Va`` (Bbus does not involve sn_mva), the NaN only appears
    later when the results are scaled back to MW, and the size/finiteness check on
    the solver output is deliberately reserved for external solvers."""

    BAD = (0.0, -100.0, float("nan"), float("inf"))

    def test_set_sn_mva_rejects_degenerate_values(self):
        grid = build_exotic_elements_case_grid()
        for v in self.BAD:
            with self.subTest(sn_mva=v), self.assertRaises(RuntimeError):
                grid.set_sn_mva(v)

    def test_set_init_vm_pu_rejects_degenerate_values(self):
        grid = build_exotic_elements_case_grid()
        for v in self.BAD:
            with self.subTest(init_vm_pu=v), self.assertRaises(RuntimeError):
                grid.set_init_vm_pu(v)

    def test_check_grid_rejects_degenerate_sn_mva_from_a_state(self):
        # set_state assigns sn_mva straight from the serialized state, bypassing
        # the setter: check_grid() is the backstop for pickle / the binary format
        state = build_exotic_elements_case_grid().__getstate__()
        for v in self.BAD:
            with self.subTest(sn_mva=v), self.assertRaises(RuntimeError):
                _rebuild(_set_at(state, (3, 5), v))

    def test_check_grid_rejects_degenerate_init_vm_pu_from_a_state(self):
        state = build_exotic_elements_case_grid().__getstate__()
        for v in self.BAD:
            with self.subTest(init_vm_pu=v), self.assertRaises(RuntimeError):
                _rebuild(_set_at(state, (3, 4), v))

    def test_good_scalars_still_accepted(self):
        grid = build_exotic_elements_case_grid()
        grid.set_sn_mva(100.0)
        grid.set_init_vm_pu(1.04)
        self.assertIsNone(grid.check_grid())
        v0 = np.ones(grid.total_bus(), dtype=complex)
        self.assertNotEqual(len(grid.ac_pf(v0, 20, 1e-7)), 0)

    def test_powermodels_loader_rejects_a_degenerate_baseMVA(self):
        # init_from_matpower routes through this loader, which does
        # `set_sn_mva(float(network["baseMVA"]))` straight from the source file
        import copy
        import json
        import os
        try:
            from lightsim2grid.network import init_from_powermodels
        except ImportError:
            self.skipTest("powermodels loader not available")
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "pf_delta_case14.json")
        if not os.path.exists(path):
            self.skipTest("pf_delta_case14.json fixture not available")
        with open(path) as f:
            network = json.load(f)["network"]
        for v in (0.0, -100.0, float("nan")):
            bad = copy.deepcopy(network)
            bad["baseMVA"] = v
            with self.subTest(baseMVA=v), self.assertRaises(RuntimeError):
                init_from_powermodels(bad)


class TestPandaPowerConverterInputSizes(unittest.TestCase):
    """``PandaPowerConverter`` walks its inputs together element by element with
    unchecked accessors, and combines them in coefficient-wise Eigen expressions
    sized from one of them. Eigen's own size assertions are compiled out of the
    release wheels, so mismatched lengths used to write past the end of the
    shorter arrays -- observed as heap corruption (``free(): invalid next size``),
    not an exception."""

    def setUp(self):
        self.conv = PandaPowerConverter()
        self.conv.set_f_hz(50.0)
        self.conv.set_sn_mva(100.0)
        self.long = np.ones(200)
        self.one = np.ones(1)

    def _trafo_args(self, first):
        return (first, self.one, self.one, [True],
                self.one, self.one, self.one, self.one, self.one, self.one, self.one)

    def test_get_trafo_param_pp2_rejects_mismatched_sizes(self):
        with self.assertRaises(RuntimeError):
            self.conv.get_trafo_param_pp2(*self._trafo_args(self.long))

    def test_get_trafo_param_pp3_rejects_mismatched_sizes(self):
        with self.assertRaises(RuntimeError):
            self.conv.get_trafo_param_pp3(*self._trafo_args(self.long), True)

    def test_check_init_rejects_a_nan_sn_mva_or_f_hz(self):
        # the guard was `sn_mva_ <= 0.`, and every comparison with NaN is false, so a
        # NaN went straight through -- and it is a divisor in every method here
        conv = PandaPowerConverter()
        conv.set_sn_mva(float("nan"))
        conv.set_f_hz(50.0)
        with self.assertRaises(RuntimeError):
            conv.get_line_param(self.one, self.one, self.one, self.one, self.one, self.one)
        conv2 = PandaPowerConverter()
        conv2.set_sn_mva(100.0)
        conv2.set_f_hz(float("nan"))
        with self.assertRaises(RuntimeError):
            conv2.get_line_param(self.one, self.one, self.one, self.one, self.one, self.one)

    def test_get_trafo_param_rejects_a_short_bool_vector(self):
        args = (self.one, self.one, self.one, [],  # is_tap_hv_side empty
                self.one, self.one, self.one, self.one, self.one, self.one, self.one)
        with self.assertRaises(RuntimeError):
            self.conv.get_trafo_param_pp2(*args)

    def test_get_line_param_rejects_mismatched_sizes(self):
        with self.assertRaises(RuntimeError):
            self.conv.get_line_param(self.long, self.one, self.one, self.one,
                                     self.one, self.one)

    def test_get_line_param_legacy_rejects_mismatched_sizes(self):
        with self.assertRaises(RuntimeError):
            self.conv.get_line_param_legacy(self.one, self.one, self.one, self.long,
                                            self.one, self.one)

    def test_malformed_pandapower_net_raises_instead_of_corrupting_the_heap(self):
        # the same hole, reached through init_from_pandapower: the loader looks the
        # transformer / line end voltages up with `pp_net.bus.loc[pp_net.trafo["hv_bus"]]`,
        # so a net whose bus index carries a DUPLICATE label hands the converter more
        # nominal voltages than there are branches.
        try:
            import pandas as pd
            import pandapower.networks as pn
            from lightsim2grid.network import init_from_pandapower
        except ImportError:
            self.skipTest("pandapower is not installed")
        net = pn.case14()
        net.bus = pd.concat([net.bus, net.bus.loc[[3]].copy()])  # bus label 3 twice
        with self.assertRaises(RuntimeError):
            init_from_pandapower(net)

    def test_matching_sizes_still_work(self):
        n = 3
        v = np.ones(n)
        r, x, h = self.conv.get_trafo_param_pp2(v, v, v, [True] * n, v, v, v, v, v, v, v)
        self.assertEqual(r.shape[0], n)
        lr, lx, h_or, h_ex = self.conv.get_line_param(v, v, v, v, v, v)
        self.assertEqual(lr.shape[0], n)


class TestBusVoltageLimitSymmetry(unittest.TestCase):
    """bus_vmin_kv_ / bus_vmax_kv_ are each optional (empty or complete), but they
    are consumed *together* -- ContingencyAnalysis::check_bus_voltage_violations
    loops to vmin.size() and then reads vmax(grid_id). A state with one set and the
    other empty passes the per-field checks yet reads past the end of the empty one
    (segfault before the fix). SUBSTATION_ID sub-state layout: n_sub(0), nmax(1),
    sub_vn_kv(2), bus_status(3), bus_vn_kv(4), sub_names(5), bus_vmin(6), bus_vmax(7).
    """

    SUBSTATION_ID = 7
    VMIN_FIELD = 6
    VMAX_FIELD = 7

    def _grid_with_limits(self):
        grid = build_exotic_elements_case_grid()
        nb_bus = len(grid.get_bus_vn_kv())
        grid.set_bus_voltage_limits(np.full(nb_bus, 90.0), np.full(nb_bus, 160.0))
        return grid

    def test_both_set_round_trips(self):
        grid = pickle.loads(pickle.dumps(self._grid_with_limits()))
        self.assertIsNone(grid.check_grid())
        self.assertEqual(len(grid.get_bus_vmin_kv()), len(grid.get_bus_vmax_kv()))

    def test_vmax_empty_while_vmin_set_is_rejected(self):
        state = self._grid_with_limits().__getstate__()
        bad = _set_at(state, (3, self.SUBSTATION_ID, self.VMAX_FIELD), ())
        with self.assertRaises(RuntimeError):
            _rebuild(bad)

    def test_vmin_empty_while_vmax_set_is_rejected(self):
        state = self._grid_with_limits().__getstate__()
        bad = _set_at(state, (3, self.SUBSTATION_ID, self.VMIN_FIELD), ())
        with self.assertRaises(RuntimeError):
            _rebuild(bad)

    def test_contingency_voltage_violations_stay_in_bounds(self):
        # the honest arrangement the check protects: with both vectors present the
        # inline voltage-violation check runs without over-reading.
        grid = self._grid_with_limits()
        sa = ContingencyAnalysisCPP(grid, True)  # compute_limit_violations=True
        sa.add_n1(0)
        sa.compute(np.ones(grid.total_bus(), dtype=complex), 20, 1e-7)
        self.assertTrue(sa.converged_n())


class TestPosTopoVectSetterBounds(unittest.TestCase):
    """set_*_pos_topo_vect() only length-checks its argument; the *values* are
    range-checked in check_grid(), which the setters do not call. update_topo()
    then indexes its caller arrays with the stored positions (unchecked Eigen
    operator()), so a position written straight through a setter -- bypassing
    check_grid -- read out of bounds before the fix (segfault)."""

    def _grid_with_topo(self):
        grid = build_exotic_elements_case_grid()
        n_load = len(grid.get_loads())
        n_gen = len(grid.get_generators())
        n_sto = len(grid.get_storages())
        n_line = len(grid.get_lines())
        n_trafo = len(grid.get_trafos())
        dim_topo = n_load + n_gen + n_sto + 2 * n_line + 2 * n_trafo
        pos = np.arange(dim_topo, dtype=np.int32)
        off = 0
        for setter, nb in ((grid.set_load_pos_topo_vect, n_load),
                           (grid.set_gen_pos_topo_vect, n_gen),
                           (grid.set_storage_pos_topo_vect, n_sto),
                           (grid.set_line_pos1_topo_vect, n_line),
                           (grid.set_line_pos2_topo_vect, n_line),
                           (grid.set_trafo_pos1_topo_vect, n_trafo),
                           (grid.set_trafo_pos2_topo_vect, n_trafo)):
            setter(pos[off:off + nb])
            off += nb
        return grid, dim_topo, n_load

    def test_out_of_range_position_from_setter_is_rejected_by_update_topo(self):
        grid, dim_topo, n_load = self._grid_with_topo()
        self.assertIsNone(grid.check_grid())  # honest layout validates
        # poison a load position AFTER validation; the setter accepts it silently
        bad = np.arange(n_load, dtype=np.int32)
        bad[0] = 1 << 28
        grid.set_load_pos_topo_vect(bad)
        has_changed = np.ones(dim_topo, dtype=bool)
        new_values = np.ones(dim_topo, dtype=np.int32)
        with self.assertRaises((IndexError, RuntimeError)):
            grid.update_topo(has_changed, new_values)


class TestContingencyResultsMatchDefaults(unittest.TestCase):
    """compute() sizes the result matrices with one row per registered contingency.
    Adding / removing contingencies afterwards without recomputing left the flow
    computation (clean_flows) walking the new, longer contingency set while indexing
    the old, shorter matrices -- an out-of-bounds write / heap corruption before the
    fix. compute_flows() / compute_power_flows() now refuse the stale state."""

    def _sa(self):
        grid = build_exotic_elements_case_grid()
        sa = ContingencyAnalysisCPP(grid)
        sa.add_n1(0)
        sa.compute(np.ones(grid.total_bus(), dtype=complex), 20, 1e-7)
        return sa

    def test_compute_flows_after_adding_contingency_is_rejected(self):
        sa = self._sa()
        for i in range(1, 5):
            sa.add_n1(i)
        with self.assertRaises(RuntimeError):
            sa.compute_flows()

    def test_compute_power_flows_after_adding_contingency_is_rejected(self):
        sa = self._sa()
        for i in range(1, 5):
            sa.add_n1(i)
        with self.assertRaises(RuntimeError):
            sa.compute_power_flows()

    def test_flows_are_fine_without_mutation(self):
        sa = self._sa()
        flows = sa.compute_flows()
        self.assertEqual(flows.shape[0], 1)

    def test_recomputing_clears_the_staleness(self):
        sa = self._sa()
        sa.add_n1(1)
        grid = build_exotic_elements_case_grid()
        sa.compute(np.ones(grid.total_bus(), dtype=complex), 20, 1e-7)
        flows = sa.compute_flows()  # no longer stale
        self.assertEqual(flows.shape[0], 2)


class TestSubObjectGettersKeepGridAlive(unittest.TestCase):
    """LSGrid's sub-object getters (get_lines / get_generators / get_bus_vn_kv /
    get_V_solver / ...) return a reference / numpy view into grid-owned memory. They
    were bound with return_value_policy::reference, which does NOT keep the parent
    grid alive, so the wrapper dangled once the grid was collected (segfault / heap
    disclosure). reference_internal ties their lifetime to the grid."""

    def test_container_getter_outlives_grid(self):
        lines = build_exotic_elements_case_grid().get_lines()  # grid is a temporary
        import gc
        gc.collect()
        _ = [bytearray(b"\x41" * 8192) for _ in range(2000)]
        gc.collect()
        # reading through it must still work (the grid is kept alive)
        buses = np.asarray(lines.get_bus_id_side_1())
        self.assertTrue(np.all(buses >= 0))

    def test_eigen_view_getter_outlives_grid(self):
        vn = build_exotic_elements_case_grid().get_bus_vn_kv()  # numpy view
        import gc
        gc.collect()
        _ = [bytearray(b"\x41" * 8192) for _ in range(2000)]
        gc.collect()
        self.assertIsNotNone(vn.base)          # a parent object keeps it alive
        self.assertTrue(np.all(np.asarray(vn) > 0.))


if __name__ == "__main__":
    unittest.main()
