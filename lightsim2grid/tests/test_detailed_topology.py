# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The python surface of the detailed topology (switches inside each substation),
on a grid built by hand: no pypowsybl, no grid2op. The C++ suite pins the
semantics (src/tests/test_switch_projection.cpp and friends); this checks that
they are reachable, typed and picklable from python."""

import pickle
import unittest

import numpy as np

from lightsim2grid.network import LSGrid
from lightsim2grid.elements import (SwitchKind, SwitchContainer, SwitchInfo,
                                    BusbarSectionContainer, BusbarSectionInfo,
                                    SubstationTopology)


def make_switch_grid():
    """The fixture of src/tests/test_switch_projection.cpp: three substations of
    two busbars, a slack generator in substation 0, lines 0-1, 1-2 and a second
    0-1, a load in substation 2. Substation 0 has two sections A (node 0) and B
    (node 1) with an open coupler (switch 0); the generator (node 2) and line 0's
    bay (node 3, terminal node 4 behind an internal connection) can reach either
    section (switches 1-4); line 2's end (node 5) is on A for good (switch 6).
    Substations 1 and 2 have one section each and a closed breaker per feeder
    (switches 7-9 and 10-11)."""
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(3, 2, np.full(6, 138.), 0, 0)
    grid.init_powerlines(np.array([0.01, 0.01, 0.02]),
                         np.array([0.1, 0.1, 0.2]),
                         np.zeros(3, dtype=complex),
                         np.array([0, 1, 0], dtype=np.int32),
                         np.array([1, 2, 1], dtype=np.int32))
    grid.init_loads(np.array([50.]), np.array([10.]), np.array([2], dtype=np.int32))
    grid.init_generators(np.array([0.]), np.array([1.02]), np.array([-1000.]), np.array([1000.]),
                         np.array([0], dtype=np.int32))
    grid.add_gen_slackbus(0, 1.)
    grid.set_load_to_subid(np.array([2], dtype=np.int32))
    grid.set_gen_to_subid(np.array([0], dtype=np.int32))
    grid.set_line_to_sub1_id(np.array([0, 1, 0], dtype=np.int32))
    grid.set_line_to_sub2_id(np.array([1, 2, 1], dtype=np.int32))

    B, D, I = int(SwitchKind.BREAKER), int(SwitchKind.DISCONNECTOR), int(SwitchKind.INTERNAL_CONNECTION)
    grid.init_detailed_topology(
        np.array([6, 4, 3], dtype=np.int32),
        np.array([0, 0, 1, 2], dtype=np.int32), np.array([0, 1, 0, 0], dtype=np.int32),
        np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2], dtype=np.int32),
        np.array([0, 2, 2, 3, 3, 4, 5, 1, 2, 3, 1, 2], dtype=np.int32),
        np.array([1, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 0], dtype=np.int32),
        [B, D, D, D, D, I, D, B, B, B, B, B],
        [True, False, True, False, True, False, False, False, False, False, False, False],
        [True, False, False, False, False, False, False, True, True, True, True, True])
    grid.set_load_to_node_id(np.array([1], dtype=np.int32))
    grid.set_gen_to_node_id(np.array([2], dtype=np.int32))
    grid.set_line_to_node1_id(np.array([4, 2, 5], dtype=np.int32))
    grid.set_line_to_node2_id(np.array([1, 2, 3], dtype=np.int32))
    grid.set_switch_names([f"sw{i}" for i in range(12)])
    grid.set_busbar_section_names(["A", "B", "C", "D"])
    return grid


def make_plain_grid():
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(2, 1, np.full(2, 138.), 0, 0)
    grid.init_powerlines(np.array([0.01]), np.array([0.1]), np.zeros(1, dtype=complex),
                         np.array([0], dtype=np.int32), np.array([1], dtype=np.int32))
    grid.init_loads(np.array([50.]), np.array([10.]), np.array([1], dtype=np.int32))
    grid.init_generators(np.array([0.]), np.array([1.02]), np.array([-1000.]), np.array([1000.]),
                         np.array([0], dtype=np.int32))
    grid.add_gen_slackbus(0, 1.)
    return grid


def flat_start(grid):
    return np.full(grid.total_bus(), grid.get_init_vm_pu(), dtype=complex)


class TestDetailedTopologyAPI(unittest.TestCase):
    def setUp(self):
        self.grid = make_switch_grid()

    def test_declared(self):
        assert self.grid.has_detailed_topology()
        subs = self.grid.get_substations()
        assert subs[0].nb_nodes == 6
        assert subs[0].nb_switches == 7
        assert subs[0].nb_busbar_sections == 2
        assert subs[1].first_node == 6
        assert subs[1].first_switch == 7
        assert subs[2].first_busbar_section == 3

    def test_switches(self):
        switches = self.grid.get_switches()
        assert isinstance(switches, SwitchContainer)
        assert len(switches) == 12
        all_sw = list(switches)
        assert len(all_sw) == 12
        assert [sw.id for sw in all_sw] == list(range(12))
        sw7 = switches[7]
        assert isinstance(sw7, SwitchInfo)
        assert sw7.name == "sw7"
        assert sw7.sub_id == 1
        assert sw7.voltage_level_id == 1
        assert sw7.local_id == 0
        assert (sw7.node1, sw7.node2) == (1, 0)
        assert sw7.kind == SwitchKind.BREAKER
        assert not sw7.open
        assert sw7.retained
        assert switches[5].kind == SwitchKind.INTERNAL_CONNECTION
        assert switches[0].open
        with self.assertRaises(ValueError):  # out of bound: std::range_error
            switches[12]

    def test_busbar_sections(self):
        sections = self.grid.get_busbar_sections()
        assert isinstance(sections, BusbarSectionContainer)
        assert len(sections) == 4
        assert [bbs.name for bbs in sections] == ["A", "B", "C", "D"]
        b = sections[1]
        assert isinstance(b, BusbarSectionInfo)
        assert b.sub_id == 0
        assert b.local_id == 1
        assert b.node == 1
        # B is isolated: not a bus
        assert not b.connected
        assert b.bus_id == -1
        assert not b.has_res
        a = sections[0]
        assert a.connected
        assert a.bus_id == 0

    def test_node_ids_on_the_elements(self):
        assert self.grid.get_loads()[0].node_id == 1
        assert self.grid.get_generators()[0].node_id == 2
        assert self.grid.get_lines()[0].node1_id == 4
        assert self.grid.get_lines()[0].node2_id == 1
        assert self.grid.get_lines()[2].node1_id == 5
        node_bus = self.grid.get_node_bus()
        assert node_bus.shape == (13,)
        assert node_bus.tolist() == [0, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2]

    def test_substation_topology(self):
        topo = self.grid.get_substation_topology(0)
        assert isinstance(topo, SubstationTopology)
        assert topo.nb_nodes() == 6
        assert topo.nb_switches() == 7
        assert topo.nb_busbar_sections() == 2
        assert topo.nb_buses() == 1
        assert topo.labels_ready()
        assert topo.node_bus().tolist() == [1, -1, 1, 1, 1, 1]
        assert topo.sw_kind(5) == SwitchKind.INTERNAL_CONNECTION
        assert topo.sw_name(0) == "sw0"
        assert topo.bbs_name(1) == "B"
        assert topo.is_open(0)
        assert topo.bbs_bus(1) == -1
        with self.assertRaises(RuntimeError):
            make_plain_grid().get_substation_topology(0)

    def test_set_switch_open(self):
        grid = self.grid
        # gen and line 0's bay to section B: B becomes bus 3, A keeps line 2
        assert grid.set_switch_open(1, True)
        assert grid.set_switch_open(2, False)
        assert grid.set_switch_open(3, True)
        assert grid.set_switch_open(4, False)
        assert grid.get_generators()[0].bus_id == 3
        assert grid.get_lines()[0].bus1_id == 3
        assert grid.get_lines()[2].bus1_id == 0
        assert grid.get_busbar_sections()[1].bus_id == 3
        assert grid.get_substation_topology(0).nb_buses() == 2
        # already there: nothing moves
        assert not grid.set_switch_open(1, True)
        V = grid.ac_pf(flat_start(grid), 20, 1e-8)
        assert len(V) > 0
        sections = grid.get_busbar_sections()
        assert sections[1].has_res
        assert abs(sections[1].res_v_kv - grid.get_Vm()[3] * 138.) < 1e-9
        # closing the coupler merges everything back on bus 0 (section A first)
        assert grid.set_switch_open(0, False)
        assert grid.get_generators()[0].bus_id == 0
        assert grid.get_substation_topology(0).nb_buses() == 1

    def test_update_switches(self):
        grid = self.grid
        has_changed = np.zeros(12, dtype=bool)
        new_open = np.zeros(12, dtype=bool)
        has_changed[[1, 2, 3, 4]] = True
        new_open[[1, 3]] = True
        grid.update_switches(has_changed, new_open)
        assert grid.get_generators()[0].bus_id == 3
        assert grid.get_lines()[0].bus1_id == 3
        # the load's breaker: disconnected, its section still a bus (line 1's end)
        has_changed[:] = False
        has_changed[10] = True
        new_open[10] = True
        grid.update_switches(has_changed, new_open)
        assert not grid.get_loads()[0].connected
        assert grid.get_busbar_sections()[3].connected
        with self.assertRaises(RuntimeError):
            grid.update_switches(np.zeros(3, dtype=bool), np.zeros(3, dtype=bool))

    def test_errors(self):
        with self.assertRaises(RuntimeError):  # an internal connection
            self.grid.set_switch_open(5, True)
        with self.assertRaises(IndexError):
            self.grid.set_switch_open(12, True)
        with self.assertRaises(RuntimeError):  # unsorted substations
            make_plain_grid().init_detailed_topology(
                np.array([2, 2], dtype=np.int32),
                np.array([0, 1], dtype=np.int32), np.array([0, 0], dtype=np.int32),
                np.array([1, 0], dtype=np.int32), np.array([0, 0], dtype=np.int32),
                np.array([1, 1], dtype=np.int32),
                [int(SwitchKind.BREAKER)] * 2, [False, False], [True, True])

    def test_pickle_and_copy(self):
        grid = self.grid
        grid.set_switch_open(1, True)
        grid.set_switch_open(2, False)
        for other in (pickle.loads(pickle.dumps(grid)), grid.copy()):
            assert other.has_detailed_topology()
            assert len(other.get_switches()) == 12
            assert other.get_switches()[1].open
            assert other.get_switches()[7].name == "sw7"
            assert other.get_busbar_sections()[3].name == "D"
            assert other.get_generators()[0].node_id == 2
            assert other.get_node_bus().tolist() == grid.get_node_bus().tolist()
            # and it keeps behaving the same
            assert other.set_switch_open(0, False)
            assert grid.copy().set_switch_open(0, False)

    def test_project_switches_overwrites_the_escape_hatch(self):
        grid = self.grid
        grid.deactivate_load(0)
        assert not grid.get_loads()[0].connected
        grid.project_switches()
        assert grid.get_loads()[0].connected
        assert grid.get_loads()[0].bus_id == 2


class TestNoDetailedTopology(unittest.TestCase):
    def test_plain_grid(self):
        grid = make_plain_grid()
        assert not grid.has_detailed_topology()
        assert len(grid.get_switches()) == 0
        assert len(grid.get_busbar_sections()) == 0
        assert list(grid.get_switches()) == []
        assert grid.get_node_bus().shape == (0,)
        assert grid.get_loads()[0].node_id == -1
        assert grid.get_lines()[0].node1_id == -1
        sub0 = grid.get_substations()[0]
        assert sub0.nb_nodes == 0
        assert sub0.nb_switches == 0
        assert sub0.first_node == -1
        with self.assertRaises(RuntimeError):
            grid.project_switches()
        with self.assertRaises(RuntimeError):
            grid.set_switch_open(0, True)
        restored = pickle.loads(pickle.dumps(grid))
        assert not restored.has_detailed_topology()
        assert len(restored.ac_pf(flat_start(restored), 20, 1e-8)) > 0

    def test_switch_kind_values(self):
        assert int(SwitchKind.BREAKER) == 0
        assert int(SwitchKind.DISCONNECTOR) == 1
        assert int(SwitchKind.LOAD_BREAK_SWITCH) == 2
        assert int(SwitchKind.INTERNAL_CONNECTION) == 3


if __name__ == "__main__":
    unittest.main()
