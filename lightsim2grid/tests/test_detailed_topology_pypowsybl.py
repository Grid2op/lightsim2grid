# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`init_from_pypowsybl(detailed_topology=...)`: the switches of a pypowsybl
grid, read into the model, must put every element exactly where pypowsybl's
bus view puts it -- when the grid is loaded and after any switch is operated
on both sides."""

import os
import pickle
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd
import pypowsybl.network as pn

from lightsim2grid.network import init_from_pypowsybl, LSGrid
from lightsim2grid.network.compare_lsgrid import compare_network_input
from lightsim2grid.network.from_pypowsybl._aux_add_detailed_topology import _aux_scan_detailed_topology
from lightsim2grid.elements import SwitchKind


def _load(net, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pypowsybl(net, **kwargs)


def _flat(grid):
    return np.full(grid.total_bus(), grid.get_init_vm_pu(), dtype=complex)


def _element_buses(grid, net):
    """{(kind, pypowsybl id): (connected, bus)} as the model sees it and as
    pypowsybl's bus view sees it. The model's elements are in sorted-id order
    (sort_index=True), like the frames here."""
    ls, pp = {}, {}
    for kind, els, df in [("gen", grid.get_generators(), net.get_generators()),
                          ("load", grid.get_loads(), net.get_loads()),
                          ("shunt", grid.get_shunts(), net.get_shunt_compensators()),
                          ("svc", grid.get_svcs(), net.get_static_var_compensators()),
                          ("storage", grid.get_storages(), net.get_batteries())]:
        df = df.sort_index()
        assert len(els) == df.shape[0], kind
        for el, (pid, row) in zip(els, df.iterrows()):
            ls[(kind, pid)] = (bool(el.connected), int(el.bus_id))
            pp[(kind, pid)] = (bool(row["connected"]), str(row["bus_id"]))
    for kind, els, df in [("line", grid.get_lines(), net.get_lines()),
                          ("trafo", grid.get_trafos(), net.get_2_windings_transformers())]:
        df = df.sort_index()
        assert len(els) == df.shape[0], kind
        for el, (pid, row) in zip(els, df.iterrows()):
            ls[(kind + "1", pid)] = (bool(el.connected1), int(el.bus1_id))
            ls[(kind + "2", pid)] = (bool(el.connected2), int(el.bus2_id))
            pp[(kind + "1", pid)] = (bool(row["connected1"]), str(row["bus1_id"]))
            pp[(kind + "2", pid)] = (bool(row["connected2"]), str(row["bus2_id"]))
    df_dc = net.get_hvdc_lines().sort_index()
    stations = pd.concat([net.get_vsc_converter_stations(), net.get_lcc_converter_stations()])
    for el, (pid, row) in zip(grid.get_dclines(), df_dc.iterrows()):
        for side, st_id in ((1, row["converter_station1_id"]), (2, row["converter_station2_id"])):
            st = stations.loc[st_id]
            ls[(f"hvdc{side}", pid)] = (bool(getattr(el, f"connected{side}")), int(getattr(el, f"bus{side}_id")))
            pp[(f"hvdc{side}", pid)] = (bool(st["connected"]), str(st["bus_id"]))
    return ls, pp


def _check_same_partition(test, grid, net, what="", synched=False):
    """Every element is connected on both sides or on neither, and the buses
    of the model are the buses of pypowsybl, up to a renaming.

    pypowsybl keeps a branch half-open (one end connected, the other not);
    lightsim2grid does too with ``keep_half_open_lines=True``, but by default
    (``synched``) it mirrors the two ends and a branch with one end off is off
    on both -- the rule the topology vector already applies."""
    ls, pp = _element_buses(grid, net)
    ls2pp, pp2ls = {}, {}
    for key in ls:
        kind, pid = key
        expected = pp[key][0]
        if synched and kind[:-1] in ("line", "trafo"):
            expected = pp[(kind[:-1] + "1", pid)][0] and pp[(kind[:-1] + "2", pid)][0]
        test.assertEqual(ls[key][0], expected, f"{what}: {key} connected in the model: {ls[key][0]}, expected {expected} (pypowsybl: {pp[key][0]})")
        if not ls[key][0]:
            continue
        b_ls, b_pp = ls[key][1], pp[key][1]
        test.assertEqual(ls2pp.setdefault(b_ls, b_pp), b_pp, f"{what}: model bus {b_ls} holds two pypowsybl buses ({key})")
        test.assertEqual(pp2ls.setdefault(b_pp, b_ls), b_ls, f"{what}: pypowsybl bus {b_pp} is split in the model ({key})")
    return ls2pp


def _operable(grid):
    return [sw for sw in grid.get_switches() if sw.kind != SwitchKind.INTERNAL_CONNECTION]


class TestFourSubstations(unittest.TestCase):
    """The node-breaker example grid (5 voltage levels, 59 switches, 6 busbar
    sections, HVDC lines, an SVC, a shunt). Two synchronous islands, so the
    main-component clean-up is off: every element is compared."""

    def setUp(self):
        self.net = pn.create_four_substations_node_breaker_network()
        # half-open branches kept, as pypowsybl keeps them: the per-side comparison is exact
        self.grid = _load(self.net, detailed_topology=True, only_main_component=False, keep_half_open_lines=True)

    def test_declared(self):
        grid, net = self.grid, self.net
        assert grid.has_detailed_topology()
        assert len(grid.get_switches()) == net.get_switches().shape[0]
        assert len(grid.get_busbar_sections()) == net.get_busbar_sections().shape[0]
        assert len(grid.get_substations()) == net.get_voltage_levels().shape[0]
        # switch and busbar-section names are the pypowsybl ids
        assert sorted(sw.name for sw in grid.get_switches()) == sorted(net.get_switches().index)
        assert sorted(b.name for b in grid.get_busbar_sections()) == sorted(net.get_busbar_sections().index)
        # each switch belongs to its voltage level, in voltage-level order
        vls = list(net.get_voltage_levels().sort_index().index)
        sw = net.get_switches()
        for s in grid.get_switches():
            assert vls[s.sub_id] == sw.loc[s.name, "voltage_level_id"]
            assert s.open == bool(sw.loc[s.name, "open"])
        # every terminal of a node-breaker grid is described
        assert all(el.node_id >= 0 for el in grid.get_generators())
        assert all(el.node1_id >= 0 and el.node2_id >= 0 for el in grid.get_lines())
        assert all(el.node1_id >= 0 and el.node2_id >= 0 for el in grid.get_dclines())
        assert all(el.node_id >= 0 for el in grid.get_svcs())

    def test_initial_state_matches_bus_view(self):
        _check_same_partition(self, self.grid, self.net, "initial state")

    def test_busbar_sections_match(self):
        ls2pp = _check_same_partition(self, self.grid, self.net)
        pp_bbs = self.net.get_busbar_sections(all_attributes=True)
        for b in self.grid.get_busbar_sections():
            assert b.connected == bool(pp_bbs.loc[b.name, "connected"]), b.name
            if b.connected and b.bus_id in ls2pp:
                assert ls2pp[b.bus_id] == pp_bbs.loc[b.name, "bus_id"], b.name

    def test_capacity(self):
        scan = _aux_scan_detailed_topology(self.net, True)
        buses_per_vl = self.net.get_buses().groupby("voltage_level_id").size().max()
        expected = max(int(scan.bound_per_vl.max()), int(buses_per_vl))
        assert self.grid.get_substations()[0].nb_max_busbars == expected
        # the bound is what the rule allows: sections + branch ends, capped by the feeders
        assert expected >= buses_per_vl
        # a capacity below it is refused
        with self.assertRaises(RuntimeError):
            _load(self.net, detailed_topology=True, n_busbar_per_sub=1)

    def test_random_switches_follow_pypowsybl(self):
        grid, net = self.grid, self.net
        rng = np.random.default_rng(0)
        switches = _operable(grid)
        state = {sw.id: sw.open for sw in switches}
        for it in range(50):
            sw = switches[int(rng.integers(len(switches)))]
            new_open = not state[sw.id]
            assert grid.set_switch_open(sw.id, new_open)
            net.update_switches(id=sw.name, open=new_open)
            state[sw.id] = new_open
            _check_same_partition(self, grid, net, f"iteration {it}, switch {sw.name} open={new_open}")

    def test_mirrored_ends_by_default(self):
        """Without keep_half_open_lines, a branch whose one end loses its bus
        goes off on both ends (as with the topology vector); pypowsybl keeps
        the other end. The comparison knows the rule."""
        net = pn.create_four_substations_node_breaker_network()
        grid = _load(net, detailed_topology=True, only_main_component=False)
        rng = np.random.default_rng(3)
        switches = _operable(grid)
        state = {sw.id: sw.open for sw in switches}
        for it in range(30):
            sw = switches[int(rng.integers(len(switches)))]
            new_open = not state[sw.id]
            grid.set_switch_open(sw.id, new_open)
            net.update_switches(id=sw.name, open=new_open)
            state[sw.id] = new_open
            _check_same_partition(self, grid, net, f"iteration {it}, switch {sw.name} open={new_open}", synched=True)
        # and one concrete case: the far end's disconnector of a line
        net = pn.create_four_substations_node_breaker_network()
        grid = _load(net, detailed_topology=True, only_main_component=False)
        sw = [s for s in grid.get_switches() if s.name == "S3VL1_BBS_LINES3S4_DISCONNECTOR"][0]
        grid.set_switch_open(sw.id, True)
        net.update_switches(id=sw.name, open=True)
        line = [el for el in grid.get_lines() if el.name == "LINE_S3S4"][0]
        assert not line.connected_global and not line.connected1 and not line.connected2
        pp = net.get_lines().loc["LINE_S3S4"]
        assert bool(pp["connected1"]) != bool(pp["connected2"])  # half-open in pypowsybl

    def test_bulk_update(self):
        grid, net = self.grid, self.net
        switches = _operable(grid)
        chosen = switches[::7]
        has_changed = np.zeros(len(grid.get_switches()), dtype=bool)
        new_open = np.zeros(len(grid.get_switches()), dtype=bool)
        for sw in chosen:
            has_changed[sw.id] = True
            new_open[sw.id] = not sw.open
            net.update_switches(id=sw.name, open=not sw.open)
        grid.update_switches(has_changed, new_open)
        _check_same_partition(self, grid, net, "bulk update")

    def test_pickle_binary_then_switch(self):
        grid, net = self.grid, self.net
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "grid.lsb")
            grid.save_binary(path)
            others = [pickle.loads(pickle.dumps(grid)), LSGrid.load_binary(path), grid.copy()]
        coupler = [sw for sw in grid.get_switches() if sw.name == "S1VL2_COUPLER"][0]
        assert not coupler.open
        net.update_switches(id="S1VL2_COUPLER", open=True)
        for g in [grid] + others:
            assert g.set_switch_open(coupler.id, True)
            _check_same_partition(self, g, net, "coupler opened after a round trip")
            assert g.get_node_bus().tolist() == grid.get_node_bus().tolist()


class TestMetrixSixBuses(unittest.TestCase):
    """The other node-breaker example (6 voltage levels, 141 switches, 12
    busbar sections), a single synchronous component that solves."""

    def setUp(self):
        self.net = pn.create_metrix_tutorial_six_buses_network()
        self.grid = _load(self.net, detailed_topology=True, only_main_component=False, keep_half_open_lines=True)

    def test_initial_state_matches_bus_view(self):
        assert self.grid.has_detailed_topology()
        assert len(self.grid.get_switches()) == self.net.get_switches().shape[0]
        _check_same_partition(self, self.grid, self.net, "initial state")
        assert len(self.grid.ac_pf(_flat(self.grid), 30, 1e-8)) > 0

    def test_random_switches_follow_pypowsybl(self):
        grid, net = self.grid, self.net
        rng = np.random.default_rng(1)
        switches = _operable(grid)
        state = {sw.id: sw.open for sw in switches}
        for it in range(30):
            sw = switches[int(rng.integers(len(switches)))]
            new_open = not state[sw.id]
            assert grid.set_switch_open(sw.id, new_open)
            net.update_switches(id=sw.name, open=new_open)
            state[sw.id] = new_open
            _check_same_partition(self, grid, net, f"iteration {it}, switch {sw.name} open={new_open}")

    def test_results_equal_a_fresh_import(self):
        """After switches moved, the projected grid solves to the same
        element results as the modified network imported from scratch."""
        grid, net = self.grid, self.net
        rng = np.random.default_rng(2)
        switches = [sw for sw in _operable(grid) if not sw.open]
        compared = 0
        for it in range(6):
            sw = switches[int(rng.integers(len(switches)))]
            grid.set_switch_open(sw.id, True)
            net.update_switches(id=sw.name, open=True)
            _check_same_partition(self, grid, net, f"iteration {it}")
            V = grid.ac_pf(_flat(grid), 30, 1e-8)
            fresh = _load(net, detailed_topology=False, only_main_component=False, keep_half_open_lines=True)
            V_fresh = fresh.ac_pf(_flat(fresh), 30, 1e-8)
            assert (len(V) > 0) == (len(V_fresh) > 0), f"iteration {it}: one converged, not the other"
            if len(V) == 0:
                continue
            compared += 1
            for get in ("get_generators", "get_loads", "get_lines", "get_trafos"):
                for a, b in zip(getattr(grid, get)(), getattr(fresh, get)()):
                    if hasattr(a, "connected"):
                        assert a.connected == b.connected, (get, a.name)
                    else:
                        assert (a.connected1, a.connected2) == (b.connected1, b.connected2), (get, a.name)
                    for attr in ("res_p_mw", "res_q_mvar", "res_p1_mw", "res_q1_mvar", "res_p2_mw"):
                        if hasattr(a, attr):
                            assert abs(getattr(a, attr) - getattr(b, attr)) < 1e-5, (get, a.name, attr)
            # the busbar sections read their bus' voltage
            for b in grid.get_busbar_sections():
                if b.connected:
                    assert b.has_res
                    assert abs(b.res_v_kv - grid.get_Vm()[b.bus_id] * grid.get_bus_vn_kv()[b.bus_id]) < 1e-9
        assert compared > 0


class TestBusBreakerGrids(unittest.TestCase):
    def test_auto_leaves_a_bus_breaker_grid_alone(self):
        net = pn.create_ieee14()
        plain = _load(net)
        auto = _load(net, detailed_topology="auto")
        assert not auto.has_detailed_topology()
        assert compare_network_input(plain, auto) == {}

    def test_synthetic_node_model(self):
        """A bus-breaker grid asked for its detailed topology gets one node per
        bus, one breaker per terminal, and behaves like the plain import."""
        net = pn.create_ieee14()
        plain = _load(net)
        grid = _load(net, detailed_topology=True)
        assert grid.has_detailed_topology()
        nb_terminals = (net.get_generators().shape[0] + net.get_loads().shape[0]
                        + net.get_shunt_compensators().shape[0]
                        + 2 * net.get_lines().shape[0] + 2 * net.get_2_windings_transformers().shape[0])
        assert len(grid.get_switches()) == nb_terminals
        assert all(sw.kind == SwitchKind.BREAKER for sw in grid.get_switches())
        assert len(grid.get_busbar_sections()) == net.get_bus_breaker_view_buses().shape[0]
        _check_same_partition(self, grid, net, "ieee14")
        # one bus per voltage level here: the numbering is the plain one
        V = grid.ac_pf(_flat(grid), 30, 1e-8)
        V_plain = plain.ac_pf(_flat(plain), 30, 1e-8)
        assert len(V) > 0 and len(V_plain) > 0
        assert np.abs(V[:14] - V_plain[:14]).max() < 1e-10
        # and a load's breaker opened is the load disconnected
        sw = [s for s in grid.get_switches() if s.name == "B2-L"][0]
        grid.set_switch_open(sw.id, True)
        net.update_switches(id="B2-L", open=True) if "B2-L" in net.get_switches().index else net.disconnect("B2-L")
        _check_same_partition(self, grid, net, "load B2-L disconnected")

    def test_refusals(self):
        net = pn.create_four_substations_node_breaker_network()
        with self.assertRaises(RuntimeError):
            _load(net, detailed_topology=True, buses_for_sub=True)
        with self.assertRaises(RuntimeError):
            _load(net, detailed_topology=True, fuse_zero_impedance_branches=True)
        with self.assertRaises(RuntimeError):
            _load(net, detailed_topology=True, convert_dangling_lines=True)
        with self.assertRaises(ValueError):
            _load(net, detailed_topology="yes")
        # "auto" with a refused option silently stays off
        assert not _load(net, detailed_topology="auto", buses_for_sub=True).has_detailed_topology()


class TestBackendPassThrough(unittest.TestCase):
    def test_loader_kwargs(self):
        try:
            from lightsim2grid import LightSimBackend
            import grid2op  # noqa: F401
        except ImportError:
            self.skipTest("grid2op is not installed")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = os.path.join(dir_path, "case_14_iidm")
        backend = LightSimBackend(loader_method="pypowsybl",
                                  loader_kwargs={"detailed_topology": True, "use_buses_for_sub": False})
        type(backend)._clear_grid_dependant_class_attributes()
        backend.set_env_name("case_14_iidm_detailed_topology")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            backend.load_grid(path, "grid.xiidm")
        assert backend._grid.has_detailed_topology()
        assert len(backend._grid.get_switches()) > 0


if __name__ == "__main__":
    unittest.main()
