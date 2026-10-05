# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's TransformerVoltageControl outer loop (AFTER_GENERATOR_VOLTAGE_CONTROL)
against OpenLoadFlow itself: a ratio tap changer added to a transformer of pypowsybl's IEEE
14-bus grid, regulating the voltage of its low voltage side. The grid's generators of the low
voltage buses are frozen while the transformer acts, which the solved tap depends on."""

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
    from lightsim2grid.algorithm import TransformerVoltageControl, ReactiveLimits
    from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationCategory, ViolationElementType
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

SLACK_GEN = "B1-G"


def _net(trafo, target_kv, deadband_kv, side="TWO"):
    """IEEE 14 with a ratio tap changer on `trafo`: rho 0.9 .. 1.1 in 21 positions, at 1.0,
    regulating its side `side`."""
    net = pp.network.create_ieee14()
    rtc = pd.DataFrame.from_records(
        index="id", columns=["id", "target_deadband", "target_v", "oltc", "low_tap", "tap", "regulating", "regulated_side"],
        data=[(trafo, deadband_kv, target_kv, True, 0, 10, False, side)])
    steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho"],
                                      data=[(trafo, 0., 0., 0., 0., 0.9 + 0.01 * k) for k in range(21)])
    net.create_ratio_tap_changers(rtc, steps)
    net.update_ratio_tap_changers(id=trafo, regulating=True)
    return net


def _olf(net, reactive_limits):
    from _olf_reference import reference_parameters
    loops = ("ReactiveLimits," if reactive_limits else "") + "TransformerVoltageControl"
    params = reference_parameters(
        provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[SLACK_GEN, "bus_id"],
                  "newtonRaphsonConvEpsPerEq": "1e-10", "outerLoopNames": loops},
        distributed_slack=False, use_reactive_limits=reactive_limits, read_slack_bus=False,
        transformer_voltage_control_on=True)
    res = lf.run_ac(net, params)[0]
    assert res.status.name == "CONVERGED", res.status_text
    return res


def _grid(net, reactive_limits=False):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pypowsybl(net, gen_slack_id=SLACK_GEN, sort_index=False, buses_for_sub=False,
                                   olf_rules=OlfLoadingParameters(reactive_limits=reactive_limits))


class TestTransformerVoltageControl(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    def _compare(self, trafo, target_kv, deadband_kv, reactive_limits=False, side="TWO"):
        net = _net(trafo, target_kv, deadband_kv, side)
        _olf(net, reactive_limits)
        olf_tap = int(net.get_ratio_tap_changers(all_attributes=True).loc[trafo, "solved_tap_position"])

        grid = _grid(_net(trafo, target_kv, deadband_kv, side), reactive_limits)
        grid.change_algorithm("NROuter_KLU")
        grid.clear_outer_loops()
        if reactive_limits:
            grid.add_outer_loop(ReactiveLimits())
        grid.add_outer_loop(TransformerVoltageControl())
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().status.name, "STABLE")
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

        tr = next(el for el in grid.get_trafos() if el.name == trafo)
        self.assertEqual(tr.res_ratio_tap_position, olf_tap)
        self.assertEqual(tr.ratio_tap_position, 10)  # the input is not modified
        olf_v = iidm_bus_voltages(net)
        ls = LightsimResultNetwork(grid, net).get_buses()
        vm = ls["v_mag"] / ls["voltage_level_id"].map(net.get_voltage_levels()["nominal_v"])
        common = olf_v.index.intersection(vm.dropna().index)
        self.assertLess((olf_v.loc[common, "vm_pu"] - vm.loc[common]).abs().max(), 1e-10)
        return olf_tap

    def test_interior_tap(self):
        # the generators of the low voltage buses frozen meanwhile decide where it lands
        self.assertEqual(self._compare("T4-9-1", 12.6, 0.1), 11)
        self.assertEqual(self._compare("T4-9-1", 12.4, 0.05), 6)
        self.assertEqual(self._compare("T4-9-1", 12.5, 0.02), 8)

    def test_extreme_tap(self):
        # beyond the range: rounded to the last tap, then once more the whole procedure
        self.assertEqual(self._compare("T4-9-1", 12.0, 0.1), 0)

    def test_with_reactive_limits(self):
        self.assertEqual(self._compare("T4-9-1", 12.3, 0.05, reactive_limits=True), 3)

    def test_hidden_by_a_generator(self):
        # the bus a generator holds: the transformer does not regulate it
        self.assertEqual(self._compare("T5-6-1", 12.6, 0.1), 10)

    def test_detection(self):
        grid = _grid(_net("T4-9-1", 12.6, 0.1))
        grid.change_algorithm("NRSing_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(TransformerVoltageControl())
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        viols = [v for v in grid.get_physical_violations(True, 1e-3, 1e-4)
                 if v.violation_type == LimitViolationType.TRANSFORMER_VOLTAGE_DEADBAND]
        self.assertEqual(len(viols), 1)
        v = viols[0]
        self.assertEqual(v.element_type, ViolationElementType.BUS)
        self.assertEqual(v.category, ViolationCategory.CONTROL)
        self.assertAlmostEqual(v.limit, 12.6, places=9)
        self.assertGreater(abs(v.value - v.limit), 0.05)
        # hidden by a generator: nothing reported
        grid = _grid(_net("T5-6-1", 12.6, 0.1))
        grid.change_algorithm("NRSing_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(TransformerVoltageControl())
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10).shape[0], 0)
        self.assertFalse([v for v in grid.get_physical_violations(True, 1e-3, 1e-4)
                          if v.violation_type == LimitViolationType.TRANSFORMER_VOLTAGE_DEADBAND])

    def test_parameters_and_default_list(self):
        loop = TransformerVoltageControl(use_initial_tap_position=False, max_controlled_nominal_voltage=-1.,
                                         min_target_deadband_kv=0.2)
        self.assertFalse(loop.use_initial_tap_position)
        self.assertEqual(loop.max_controlled_nominal_voltage, -1.)
        self.assertEqual(loop.min_target_deadband_kv, 0.2)
        with self.assertRaises(AttributeError):  # set at construction only
            loop.use_initial_tap_position = True
        grid = _grid(_net("T4-9-1", 12.6, 0.1))
        self.assertNotIn("TransformerVoltageControl", [l.name() for l in grid.get_outer_loops()])


if __name__ == "__main__":
    unittest.main()
