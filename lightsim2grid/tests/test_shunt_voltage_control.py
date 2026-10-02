# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's ShuntVoltageControl outer loop (WITH_GENERATOR_VOLTAGE_CONTROL) against
OpenLoadFlow itself: two switched shunts added to a bus of pypowsybl's IEEE 14-bus grid,
regulating its voltage."""

import unittest
import warnings

import numpy as np

try:
    import pandas as pd
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
    from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
    from lightsim2grid.network.from_pypowsybl._result_network import LightsimResultNetwork
    from lightsim2grid.algorithm import ShuntVoltageControl, TransformerVoltageControl
    from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationCategory, ViolationElementType
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

SLACK_GEN = "B1-G"
SHUNTS = ["SH-A", "SH-B"]


def _net(bus, target_kv, sections=(3, 1), with_trafo_control=False):
    """IEEE 14 with two linear shunts on `bus` (10 sections of 0.02 S, 4 of 0.05 S), regulating
    it at `target_kv`; with `with_trafo_control`, T4-9-1 also regulates bus 9 at 12.6 kV."""
    net = pp.network.create_ieee14()
    vl = "VL" + bus[1:]
    sh = pd.DataFrame.from_records(
        index="id", columns=["id", "model_type", "section_count", "target_v", "target_deadband", "voltage_level_id",
                             "bus_id", "connectable_bus_id"],
        data=[("SH-A", "LINEAR", sections[0], target_kv, 0., vl, bus, bus),
              ("SH-B", "LINEAR", sections[1], target_kv, 0., vl, bus, bus)])
    model = pd.DataFrame.from_records(index="id", columns=["id", "g_per_section", "b_per_section", "max_section_count"],
                                      data=[("SH-A", 0., 0.02, 10), ("SH-B", 0., 0.05, 4)])
    net.create_shunt_compensators(sh, model)
    net.update_shunt_compensators(id=SHUNTS, voltage_regulation_on=[True, True])
    if with_trafo_control:
        rtc = pd.DataFrame.from_records(
            index="id", columns=["id", "target_deadband", "target_v", "oltc", "low_tap", "tap", "regulating", "regulated_side"],
            data=[("T4-9-1", 0.1, 12.6, True, 0, 10, False, "TWO")])
        steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho"],
                                          data=[("T4-9-1", 0., 0., 0., 0., 0.9 + 0.01 * k) for k in range(21)])
        net.create_ratio_tap_changers(rtc, steps)
        net.update_ratio_tap_changers(id="T4-9-1", regulating=True)
    return net


def _olf(net, with_trafo_control):
    from _olf_reference import reference_parameters
    loops = ("TransformerVoltageControl," if with_trafo_control else "") + "ShuntVoltageControl"
    params = reference_parameters(
        provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[SLACK_GEN, "bus_id"],
                  "newtonRaphsonConvEpsPerEq": "1e-10", "outerLoopNames": loops},
        distributed_slack=False, use_reactive_limits=False, read_slack_bus=False,
        shunt_compensator_voltage_control_on=True, transformer_voltage_control_on=with_trafo_control)
    res = lf.run_ac(net, params)[0]
    assert res.status.name == "CONVERGED", res.status_text
    return [int(c) for c in net.get_shunt_compensators(all_attributes=True).loc[SHUNTS, "solved_section_count"]]


def _grid(net):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pypowsybl(net, gen_slack_id=SLACK_GEN, sort_index=False, buses_for_sub=False,
                                   olf_rules=OlfLoadingParameters(reactive_limits=False))


class TestShuntVoltageControl(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _compare(self, bus, target_kv, sections=(3, 1), with_trafo_control=False):
        net = _net(bus, target_kv, sections, with_trafo_control)
        olf_counts = _olf(net, with_trafo_control)

        grid = _grid(_net(bus, target_kv, sections, with_trafo_control))
        grid.change_algorithm("NROuter_KLU")
        grid.clear_outer_loops()
        if with_trafo_control:
            grid.add_outer_loop(TransformerVoltageControl())
        grid.add_outer_loop(ShuntVoltageControl())
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().status.name, "STABLE")
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

        res = {s.name: s for s in grid.get_shunts()}
        self.assertEqual([res[s].res_section_count for s in SHUNTS], olf_counts)
        self.assertEqual([res[s].section_count for s in SHUNTS], list(sections))  # inputs untouched
        olf_v = iidm_bus_voltages(net)
        ls = LightsimResultNetwork(grid, net).get_buses()
        vm = ls["v_mag"] / ls["voltage_level_id"].map(net.get_voltage_levels()["nominal_v"])
        common = olf_v.index.intersection(vm.dropna().index)
        self.assertLess((olf_v.loc[common, "vm_pu"] - vm.loc[common]).abs().max(), 1e-9)
        return olf_counts

    def test_sections(self):
        # the susceptance shared out, the largest shunt first
        self.assertEqual(self._compare("B9", 12.8), [2, 0])
        self.assertEqual(self._compare("B9", 13.2), [6, 3])
        self.assertEqual(self._compare("B9", 12.4), [0, 0])
        self.assertEqual(self._compare("B10", 12.8, sections=(5, 0)), [2, 1])

    def test_hidden_by_a_generator(self):
        # bus 6 is a generator's: the shunts keep their sections, not even reshuffled
        self.assertEqual(self._compare("B6", 12.8, sections=(0, 2)), [0, 2])

    def test_hidden_by_a_transformer(self):
        # a transformer regulating the same bus takes precedence
        self.assertEqual(self._compare("B9", 13.2, sections=(0, 2), with_trafo_control=True), [0, 2])

    def test_detection(self):
        grid = _grid(_net("B9", 13.2))
        grid.change_algorithm("NRSing_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(ShuntVoltageControl())
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        viols = [v for v in grid.get_physical_violations(True, 1e-3, 1e-4)
                 if v.violation_type == LimitViolationType.SHUNT_VOLTAGE_CONTROL]
        self.assertEqual(len(viols), 1)
        self.assertEqual(viols[0].element_type, ViolationElementType.BUS)
        self.assertEqual(viols[0].category, ViolationCategory.CONTROL)
        self.assertAlmostEqual(viols[0].limit, 13.2, places=9)

    def test_not_in_the_default_list(self):
        grid = _grid(_net("B9", 13.2))
        self.assertNotIn("ShuntVoltageControl", [l.name() for l in grid.get_outer_loops()])


if __name__ == "__main__":
    unittest.main()
