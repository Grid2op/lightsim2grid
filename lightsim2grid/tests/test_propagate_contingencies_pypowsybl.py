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


def _net_with_a_transformer_bay(line_has_breaker):
    """A substation S5: S5VL1 holds a busbar section, the bay of line L5 (to S3VL1, with or
    without a breaker at S5VL1) and the bay of transformer T5 (behind a breaker); S5VL2, the other
    side of T5, is connected to S2VL1 by line L6, and L7 joins S3VL1 and S2VL1."""
    net = pn.create_four_substations_node_breaker_network()
    net.create_substations(id="S5", country="FR")
    net.create_voltage_levels(id=["S5VL1", "S5VL2"], substation_id=["S5", "S5"],
                              topology_kind=["NODE_BREAKER"] * 2, nominal_v=[400., 400.],
                              high_voltage_limit=[440., 440.], low_voltage_limit=[360., 360.])
    net.create_busbar_sections(id=["S5VL1_BBS", "S5VL2_BBS"], voltage_level_id=["S5VL1", "S5VL2"], node=[0, 0])
    switches = [("S5VL1", "S5_L5_D", 0, 4 if line_has_breaker else 1, "DISCONNECTOR")]
    if line_has_breaker:
        switches.append(("S5VL1", "S5_L5_B", 4, 1, "BREAKER"))
    switches += [("S5VL1", "S5_T5_D", 0, 2, "DISCONNECTOR"), ("S5VL1", "S5_T5_B", 2, 3, "BREAKER"),
                 ("S5VL2", "S5b_T5_D", 0, 2, "DISCONNECTOR"), ("S5VL2", "S5b_T5_B", 2, 1, "BREAKER"),
                 ("S5VL2", "S5b_L6_D", 0, 4, "DISCONNECTOR"), ("S5VL2", "S5b_L6_B", 4, 3, "BREAKER"),
                 ("S3VL1", "S3_L5_D", 0, 60, "DISCONNECTOR"), ("S3VL1", "S3_L5_B", 60, 61, "BREAKER"),
                 ("S2VL1", "S2_L6_D", 0, 60, "DISCONNECTOR"), ("S2VL1", "S2_L6_B", 60, 61, "BREAKER"),
                 ("S3VL1", "S3_L7_D", 0, 62, "DISCONNECTOR"), ("S3VL1", "S3_L7_B", 62, 63, "BREAKER"),
                 ("S2VL1", "S2_L7_D", 0, 62, "DISCONNECTOR"), ("S2VL1", "S2_L7_B", 62, 63, "BREAKER")]
    net.create_switches(id=[s[1] for s in switches], voltage_level_id=[s[0] for s in switches],
                        node1=[s[2] for s in switches], node2=[s[3] for s in switches],
                        kind=[s[4] for s in switches], open=[False] * len(switches))
    net.create_lines(id=["L5", "L6", "L7"], voltage_level1_id=["S5VL1", "S5VL2", "S3VL1"], node1=[1, 3, 63],
                     voltage_level2_id=["S3VL1", "S2VL1", "S2VL1"], node2=[61, 61, 63],
                     r=[1.] * 3, x=[10.] * 3, g1=[0.] * 3, b1=[0.] * 3, g2=[0.] * 3, b2=[0.] * 3)
    net.create_2_windings_transformers(id="T5", voltage_level1_id="S5VL1", node1=3, voltage_level2_id="S5VL2",
                                       node2=1, rated_u1=400., rated_u2=400., r=0.1, x=10., g=0., b=0.)
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

    def test_equipment_tripped_on_the_way_is_walked_out_of(self):
        # LINE_X hangs on the feeder of LINE_S2S3 (a breakerless tee), LOAD_X on the far end of
        # LINE_X: tripping either end of the chain takes the whole chain down
        net = pn.create_four_substations_node_breaker_network()
        net.create_switches(id=["SW_X2", "SW_X3"], voltage_level_id=["S2VL1", "S3VL1"],
                            node1=[6, 31], node2=[23, 32], kind=["DISCONNECTOR", "DISCONNECTOR"],
                            open=[False, False])
        net.create_lines(id="LINE_X", voltage_level1_id="S2VL1", node1=23, voltage_level2_id="S3VL1",
                         node2=31, r=0.1, x=1., g1=0., b1=0., g2=0., b2=0.)
        net.create_loads(id="LOAD_X", voltage_level_id="S3VL1", node=32, p0=1., q0=0.)
        res = propagate_contingencies(net, [["LINE_S2S3"], ["LOAD_X"], ["LINE_X"]])
        self.assertEqual(res, [["LINE_S2S3", "LINE_X", "LOAD_X"],
                               ["LOAD_X", "LINE_X", "LINE_S2S3"],
                               ["LINE_X", "LINE_S2S3", "LOAD_X"]])
        unique, inverse = group_equivalent_contingencies(res)
        self.assertEqual(inverse, [0, 0, 0])

    def test_a_load_break_switch_is_not_a_breaker(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(kind="LOAD_BREAK_SWITCH"), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER"]])

    def test_a_breaker_stops_the_walk(self):
        # the walk stops at the load's breaker, which opens: the load is then left alone behind
        # it, so it is disconnected, but not walked out of
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(kind="BREAKER"), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3", "LOAD_FEEDER"]])
        # what tells "stopped" from "walked through": the load still fed from the busbar
        # section (node 0) behind the breaker stays. Walked through, the walk would trip it
        net = _net_with_stuff_on_the_feeder(kind="BREAKER")
        net.create_switches(id="SW_BBS_FEEDER", voltage_level_id="S2VL1", node1=0, node2=20,
                            kind="DISCONNECTOR", open=False)
        self.assertEqual(propagate_contingencies(net, [["LINE_S2S3"]]), [["LINE_S2S3"]])

    def test_an_open_switch_stops_the_walk(self):
        res = propagate_contingencies(_net_with_stuff_on_the_feeder(is_open=True), [["LINE_S2S3"]])
        self.assertEqual(res, [["LINE_S2S3"]])

    def test_tripping_the_load_takes_the_line(self):
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
        # tripping a busbar section opens the breakers around it: each feeder of that section is
        # left alone behind its breaker, so it is disconnected; the coupler leads to the other
        # section, still energized, so that section and its feeders stay
        net = pn.create_four_substations_node_breaker_network()
        self.assertEqual(propagate_contingencies(net, [["S1VL2_BBS1"]]),
                         [["S1VL2_BBS1", "SHUNT", "GH3", "GH2", "GH1", "TWT"]])

    def test_equipment_left_alone_behind_an_opened_breaker(self):
        # S5VL1: a busbar section, the bay of L5 and the bay of T5 (behind its breaker); T5's
        # other side is connected through L6, so nothing is islanded (OpenLoadFlow's security
        # analysis reports the same disconnected elements on this network)
        for line_has_breaker, expected in (
                (True, [["L5"], ["T5"], ["S5VL1_BBS", "T5", "L5"]]),
                # L5 has no breaker at S5VL1: its walk takes the busbar section, whose breakers
                # open, which leaves T5 alone in its bay
                (False, [["L5", "T5"], ["T5"], ["S5VL1_BBS", "L5", "T5"]])):
            with self.subTest(line_has_breaker=line_has_breaker):
                net = _net_with_a_transformer_bay(line_has_breaker)
                self.assertEqual(propagate_contingencies(net, [["L5"], ["T5"], ["S5VL1_BBS"]]), expected)

    def test_opened_sides(self):
        # what is cut off behind an opened breaker is disconnected at that end only (as
        # OpenLoadFlow opens the side of each terminal it disconnects); what the walk meets, and
        # the contingency itself, at all ends
        for line_has_breaker, expected in (
                (True, [{}, {}, {"T5": 1, "L5": 1}]),
                (False, [{"T5": 1}, {}, {"T5": 1}])):
            with self.subTest(line_has_breaker=line_has_breaker):
                net = _net_with_a_transformer_bay(line_has_breaker)
                propagator = ContingencyPropagator(net)
                propagator.add_contingencies([["L5"], ["T5"], ["S5VL1_BBS"]])
                self.assertEqual(propagator.opened_sides(), expected)
        # the line itself, met by the walk: opened at both ends
        net = _net_with_stuff_on_the_feeder()
        propagator = ContingencyPropagator(net)
        propagator.add_contingencies([["LOAD_FEEDER"]])
        self.assertEqual(propagator.propagate_contingencies(), [["LOAD_FEEDER", "LINE_S2S3"]])
        self.assertEqual(propagator.opened_sides(), [{}])

    def test_a_branch_left_at_a_dead_end(self):
        # L5 has no breaker at S5VL1: its walk takes the busbar section and opens T5's breaker.
        # With T5 further behind it (a disconnector, a fictitious breaker or an internal
        # connection in between), powsybl only opens the breaker: T5 stays in service, fed from
        # S5VL2, with a bus of its own at S5VL1 (OpenLoadFlow's security analysis reports it so)
        for link in ("DISCONNECTOR", "BREAKER", "INTERNAL_CONNECTION"):
            with self.subTest(link=link):
                net = _net_with_a_transformer_bay(line_has_breaker=False)
                net.remove_elements(["T5"])
                if link == "INTERNAL_CONNECTION":
                    net.create_internal_connections(voltage_level_id="S5VL1", node1=3, node2=30)
                else:
                    net.create_switches(id="S5_T5_X", voltage_level_id="S5VL1", node1=3, node2=30, kind=link,
                                        open=False, fictitious=(link == "BREAKER"))
                net.create_2_windings_transformers(id="T5", voltage_level1_id="S5VL1", node1=30,
                                                   voltage_level2_id="S5VL2", node2=1, rated_u1=400.,
                                                   rated_u2=400., r=0.1, x=10., g=0., b=0.)
                propagator = ContingencyPropagator(net)
                propagator.add_contingencies([["L5"], ["T5"]])
                self.assertEqual(propagator.propagate_contingencies(), [["L5"], ["T5"]])
                self.assertEqual(propagator.opened_sides(), [{}, {}])
                self.assertEqual(propagator.dead_ends(), [{"T5": 1}, {}])
                # simulated opened at that side, which tells it apart from L5 alone
                self.assertEqual(group_equivalent_contingencies([["L5"], ["L5"]], propagator.dead_ends()[:1] + [{}]),
                                 ([["L5"], ["L5"]], [0, 1]))
                self.assertEqual(propagator.group_contingencies(), ([["L5"], ["T5"]], [0, 1]))
        # an open switch on T5's node does not count: T5 is still right behind its breaker
        net = _net_with_a_transformer_bay(line_has_breaker=False)
        net.create_loads(id="LD_X", voltage_level_id="S5VL1", node=31, p0=1., q0=0.)
        net.create_switches(id="S5_X", voltage_level_id="S5VL1", node1=3, node2=31, kind="DISCONNECTOR", open=True)
        propagator = ContingencyPropagator(net)
        propagator.add_contingencies([["L5"]])
        self.assertEqual(propagator.propagate_contingencies(), [["L5", "T5"]])
        self.assertEqual(propagator.opened_sides(), [{"T5": 1}])
        self.assertEqual(propagator.dead_ends(), [{}])

    def test_group_tells_the_opened_sides_apart(self):
        unique, inverse = group_equivalent_contingencies([["A", "B"], ["B", "A"], ["A", "B"]],
                                                         [{"B": 1}, {"B": 1}, {}])
        self.assertEqual(unique, [["A", "B"], ["A", "B"]])
        self.assertEqual(inverse, [0, 0, 1])

    def test_a_fictitious_breaker_is_walked_through(self):
        # powsybl does not open a fictitious breaker: L5's is then as if L5 had none
        net = _net_with_a_transformer_bay(line_has_breaker=True)
        net.update_switches(id="S5_L5_B", fictitious=True)
        self.assertEqual(propagate_contingencies(net, [["L5"]]), [["L5", "T5"]])

    def test_what_is_left_behind_an_opened_breaker(self):
        # L5 has no breaker at S5VL1, so the walk takes the busbar section and opens T5's breaker
        # and the breakers of the other feeders
        net = _net_with_a_transformer_bay(line_has_breaker=False)
        # a load next to T5, behind its breaker: still fed through T5, both stay
        net.create_loads(id="LD_T5", voltage_level_id="S5VL1", node=10, p0=1., q0=0.)
        net.create_switches(id="S5_LDT5_D", voltage_level_id="S5VL1", node1=3, node2=10,
                            kind="DISCONNECTOR", open=False)
        self.assertEqual(propagate_contingencies(net, [["L5"]]), [["L5"]])

        # a feeder holding only injections: nothing feeds them any more, they are disconnected
        net = _net_with_a_transformer_bay(line_has_breaker=False)
        net.create_switches(id=["S5_INJ_B", "S5_INJ_D1", "S5_INJ_D2"], voltage_level_id=["S5VL1"] * 3,
                            node1=[0, 20, 20], node2=[20, 21, 22], kind=["BREAKER", "DISCONNECTOR", "DISCONNECTOR"],
                            open=[False] * 3)
        net.create_loads(id="LD_INJ", voltage_level_id="S5VL1", node=21, p0=1., q0=0.)
        net.create_generators(id="G_INJ", voltage_level_id="S5VL1", node=22, target_p=1., min_p=0., max_p=10.,
                              target_v=400., voltage_regulator_on=False, target_q=0.)
        self.assertEqual(sorted(propagate_contingencies(net, [["L5"]])[0]), ["G_INJ", "L5", "LD_INJ", "T5"])

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


def _documented_n1_example():
    """The code of the full N-1 example of docs/network.rst (the block that follows "A full
    N-1 of the lines and transformers"), or None when the docs are not next to the package."""
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "network.rst")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    start = next(i for i, line in enumerate(lines) if "A full N-1 of the lines and transformers" in line)
    start = next(i for i in range(start, len(lines)) if lines[i].strip() == ".. code-block:: python") + 1
    code = []
    for line in lines[start:]:
        if line.strip() and not line.startswith("    "):
            break
        code.append(line[4:])
    return "\n".join(code)


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestDocumentedN1Example(unittest.TestCase):
    def test_the_example_runs_and_opens_one_end_where_the_propagation_does(self):
        import numpy as np
        from lightsim2grid.network import init_from_pypowsybl
        from lightsim2grid.contingencyAnalysis import ContingencyAnalysisCPP
        code = _documented_n1_example()
        if code is None:
            self.skipTest("docs/network.rst is not available")
        # L5 has no breaker at S5VL1: its contingency cuts T5 off at side 1 only, and T5 is
        # given a charging admittance so that this differs from tripping it
        net = _net_with_a_transformer_bay(line_has_breaker=False)
        net.update_2_windings_transformers(id="T5", b=1e-4)
        propagator = ContingencyPropagator(net)
        propagator.add_contingencies([["L5"]])
        self.assertEqual(propagator.opened_sides(), [{"T5": 1}])
        grid = init_from_pypowsybl(net)
        v_init = np.ones(grid.total_bus(), dtype=complex)
        namespace = {"net": net, "v_init": v_init, "max_iter": 10, "tol": 1e-8}
        exec(compile(code, "docs/network.rst", "exec"), namespace)
        V = namespace["V"]
        conts = [[name] for name in net.get_lines().index] + \
                [[name] for name in net.get_2_windings_transformers().index]
        self.assertEqual(V.shape[0], len(conts))

        branch_id = {name: i for i, name in enumerate(list(grid.get_line_names()) + list(grid.get_trafo_names()))}

        def ref_V(**kwargs):
            ca = ContingencyAnalysisCPP(grid)
            ca.handle_disconnected_grid = True
            ca.add_contingency(**kwargs)
            ca.compute(v_init, 10, 1e-8)
            return ca.get_voltages()[0]

        row_l5 = conts.index(["L5"])
        opened = ref_V(branch_ids=[branch_id["L5"]], branch_ids_side1=[branch_id["T5"]])
        tripped = ref_V(branch_ids=[branch_id["L5"], branch_id["T5"]])
        self.assertGreater(np.abs(opened - tripped).max(), 1e-6)
        np.testing.assert_allclose(V[row_l5], opened, atol=1e-10, rtol=0.)
        # every other row is a plain branch trip
        for k, cont in enumerate(conts):
            if k == row_l5:
                continue
            np.testing.assert_allclose(V[k], ref_V(branch_ids=[branch_id[cont[0]]]), atol=1e-10, rtol=0.,
                                       err_msg=str(cont))


if __name__ == "__main__":
    unittest.main()
