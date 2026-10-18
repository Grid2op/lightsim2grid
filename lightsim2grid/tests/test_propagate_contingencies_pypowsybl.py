# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""`propagate_contingencies`: a contingency also takes down what stands on the same piece of
a node-breaker substation, up to the breakers.

The four-substations network has a breaker on every feeder, so nothing propagates on it as
delivered. The tests hang extra equipment on the feeder of LINE_S2S3 (node 6 of S2VL1, behind
the breaker between nodes 5 and 6; the busbar section is node 0) to make it do so.
"""

import unittest

try:
    import pypowsybl.network as pn
    from lightsim2grid.network import (ContingencyPropagator, propagate_contingencies,
                                       group_equivalent_contingencies)
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False


def _net_with_stuff_on_the_feeder(kind="DISCONNECTOR", is_open=False):
    """LOAD_FEEDER hangs on the line's terminal through one switch; LOAD_BEHIND_BREAKER
    hangs on the busbar side of the line's breaker (node 5)."""
    net = pn.create_four_substations_node_breaker_network()
    net.create_loads(id=["LOAD_FEEDER", "LOAD_BEHIND_BREAKER"], voltage_level_id=["S2VL1", "S2VL1"],
                     node=[20, 21], p0=[1., 1.], q0=[0., 0.])
    net.create_switches(id=["SW_FEEDER", "SW_BEHIND"], voltage_level_id=["S2VL1", "S2VL1"],
                        node1=[6, 5], node2=[20, 21], kind=[kind, "DISCONNECTOR"],
                        open=[is_open, False])
    return net


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestPropagateContingencies(unittest.TestCase):
    def test_nothing_to_propagate_with_a_breaker_on_every_feeder(self):
        net = pn.create_four_substations_node_breaker_network()
        res = propagate_contingencies(net, [["LINE_S2S3"], ["GH1"], ["LINE_S3S4", "LD5"]])
        self.assertEqual(res, [["LINE_S2S3"], ["GH1"], ["LINE_S3S4", "LD5"]])

    def test_equipment_on_the_feeder_goes_down_with_the_line(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(), [["LINE_S2S3"]])
        # the load behind the breaker is on the busbar side: it stays
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER"]])

    def test_the_other_side_of_the_line_is_walked_too(self):
        net = _net_with_stuff_on_the_feeder()
        net.create_loads(id="LOAD_S3", voltage_level_id="S3VL1", node=30, p0=1., q0=0.)
        net.create_switches(id="SW_S3", voltage_level_id="S3VL1", node1=2, node2=30,
                            kind="DISCONNECTOR", open=False)
        res = propagate_contingencies(net, [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER", "LOAD_S3"]])

    def test_a_load_break_switch_is_not_a_breaker(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(kind="LOAD_BREAK_SWITCH"), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER"]])

    def test_a_breaker_stops_the_walk(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(kind="BREAKER"), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3"]])

    def test_an_open_switch_stops_the_walk(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(is_open=True), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3"]])

    def test_tripping_the_load_does_not_take_the_line(self):
        # the walk from the load goes through the disconnector to node 6, where the line stands:
        # the line is equipment on the way, so it is tripped (powsybl does the same)
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(), [["LOAD_FEEDER"]])
        self.assertEqual(res, [["LOAD_FEEDER", "LINE_S2S3"]])

    def test_order_duplicates_and_independence_between_contingencies(self):
        net = _net_with_stuff_on_the_feeder()
        res = propagate_contingencies(net, [["LINE_S2S3", "LOAD_FEEDER", "LINE_S2S3"], ["GH1"], []])
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER"], ["GH1"], []])

    def test_hvdc_line_and_stations(self):
        net = pn.create_four_substations_node_breaker_network()
        self.assertEqual(propagate_contingencies(net, [["HVDC1"]]), [["HVDC1", "VSC1", "VSC2"]])
        self.assertEqual(propagate_contingencies(net, [["VSC1"]]), [["VSC1", "HVDC1", "VSC2"]])

    def test_busbar_section_is_walked_through(self):
        # tripping a busbar section opens its breakers and trips what is on the way, here nothing
        # since every feeder has a breaker: the section does not pull the other sections in
        net = pn.create_four_substations_node_breaker_network()
        self.assertEqual(propagate_contingencies(net, [["S1VL2_BBS1"]]), [["S1VL2_BBS1"]])

    def test_bus_breaker_network_does_not_propagate(self):
        net = pn.create_eurostag_tutorial_example1_network()
        names = list(net.get_lines().index[:1]) + list(net.get_generators().index[:1])
        self.assertEqual(propagate_contingencies(net, [names]), [names])

    def test_internal_connection(self):
        # a load joined to the line's terminal by an internal connection (no switch on it)
        net = pn.create_four_substations_node_breaker_network()
        net.create_loads(id="LOAD_IC", voltage_level_id="S2VL1", node=22, p0=1., q0=0.)
        net.create_internal_connections(voltage_level_id="S2VL1", node1=6, node2=22)
        expected = [["LINE_S2S3", "LOAD_IC"]]
        self.assertEqual(propagate_contingencies(net, [["LINE_S2S3"]]), expected)
        self.assertEqual(propagate_contingencies(net, [["LINE_S2S3"]], internal_connections=True), expected)
        self.assertEqual(propagate_contingencies(net, [["LINE_S2S3"]], internal_connections=False),
                         [["LINE_S2S3"]])

    def test_group_equivalent_contingencies(self):
        net = _net_with_stuff_on_the_feeder()
        propagated = propagate_contingencies(net, [["LINE_S2S3"], ["GH1"], ["LOAD_FEEDER"], ["GH1"]])
        # tripping the line or the load on its feeder ends in the same grid
        unique, inverse = group_equivalent_contingencies(propagated)
        self.assertEqual(unique, [["LINE_S2S3", "LOAD_FEEDER"], ["GH1"]])
        self.assertEqual(inverse, [0, 1, 0, 1])
        self.assertEqual([unique[i] for i in inverse], [propagated[0], propagated[1], propagated[0], propagated[1]])

    def test_group_ignores_order_and_keeps_the_first(self):
        unique, inverse = group_equivalent_contingencies([["B", "A"], ["A", "B"], ["A"], []])
        self.assertEqual(unique, [["B", "A"], ["A"], []])
        self.assertEqual(inverse, [0, 0, 1, 2])
        self.assertEqual(group_equivalent_contingencies([]), ([], []))

    def test_unknown_name_and_bad_input(self):
        net = pn.create_four_substations_node_breaker_network()
        with self.assertRaises(ValueError):
            propagate_contingencies(net, [["NOT_AN_ELEMENT"]])
        with self.assertRaises(ValueError):  # a switch is not an element to trip
            propagate_contingencies(net, [["S2VL1_LINES2S3_BREAKER"]])
        with self.assertRaises(TypeError):
            propagate_contingencies(net, ["LINE_S2S3"])



@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestContingencyPropagator(unittest.TestCase):
    def test_add_in_several_goes_propagate_and_group(self):
        propagator = ContingencyPropagator(_net_with_stuff_on_the_feeder())
        self.assertEqual(propagator.propagate_contingencies(), [])
        self.assertEqual(propagator.group_contingencies(), ([], []))
        propagator.add_contingencies([["LINE_S2S3"], ["GH1"]])
        self.assertEqual(propagator.propagate_contingencies(), [["LINE_S2S3", "LOAD_FEEDER"], ["GH1"]])
        propagator.add_contingencies([["LOAD_FEEDER", "LOAD_FEEDER"], ["GH1"]])
        self.assertEqual(propagator.contingencies, [["LINE_S2S3"], ["GH1"], ["LOAD_FEEDER"], ["GH1"]])
        self.assertEqual(propagator.propagate_contingencies(),
                         [["LINE_S2S3", "LOAD_FEEDER"], ["GH1"], ["LOAD_FEEDER", "LINE_S2S3"], ["GH1"]])
        unique, inverse = propagator.group_contingencies()
        self.assertEqual(unique, [["LINE_S2S3", "LOAD_FEEDER"], ["GH1"]])
        self.assertEqual(inverse, [0, 1, 0, 1])

    def test_same_result_as_the_function(self):
        net = _net_with_stuff_on_the_feeder()
        conts = [["LINE_S2S3"], ["HVDC1"], ["VSC1"], ["LOAD_FEEDER"], []]
        propagator = ContingencyPropagator(net)
        propagator.add_contingencies(conts)
        self.assertEqual(propagator.propagate_contingencies(), propagate_contingencies(net, conts))

    def test_returned_lists_are_copies(self):
        propagator = ContingencyPropagator(_net_with_stuff_on_the_feeder())
        propagator.add_contingencies([["LINE_S2S3"]])
        propagator.propagate_contingencies()[0].append("X")
        propagator.contingencies[0].append("X")
        self.assertEqual(propagator.propagate_contingencies(), [["LINE_S2S3", "LOAD_FEEDER"]])
        self.assertEqual(propagator.contingencies, [["LINE_S2S3"]])

    def test_bad_input_adds_nothing_and_clear(self):
        propagator = ContingencyPropagator(pn.create_four_substations_node_breaker_network())
        propagator.add_contingencies([["GH1"]])
        with self.assertRaises(ValueError):
            propagator.add_contingencies([["LINE_S2S3"], ["NOT_AN_ELEMENT"]])
        with self.assertRaises(TypeError):
            propagator.add_contingencies(["LINE_S2S3"])
        self.assertEqual(propagator.contingencies, [["GH1"]])
        propagator.clear()
        self.assertEqual(propagator.contingencies, [])
        self.assertEqual(propagator.propagate_contingencies(), [])


if __name__ == "__main__":
    unittest.main()
