# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's PhaseControl outer loop (CONTINUOUS_WITH_DISCRETISATION) against
OpenLoadFlow itself: a phase shifter added to a meshed transformer of pypowsybl's IEEE 14-bus
grid, regulating an active power or limiting a current."""

import unittest
import warnings

import numpy as np

try:
    import pandas as pd
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
    from lightsim2grid.network.from_pypowsybl._result_network import LightsimResultNetwork
    from lightsim2grid.algorithm import PhaseControl
    from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationCategory, ViolationElementType
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

PST = "T5-6-1"  # in a loop of the grid: its flow follows its shift
SLACK_GEN = "B1-G"


def _net(mode, target, regulating=True):
    """IEEE 14 with a phase tap changer on PST: -10 .. 10 degrees, its reactance changing with
    the step, at 0 degree."""
    net = pp.network.create_ieee14()
    ptc = pd.DataFrame.from_records(
        index="id", columns=["id", "target_deadband", "regulation_mode", "low_tap", "tap", "regulating", "regulated_side"],
        data=[(PST, 0., mode, -10, 0, False, "ONE")])
    steps = pd.DataFrame.from_records(index="id", columns=["id", "b", "g", "r", "x", "rho", "alpha"],
                                      data=[(PST, 0., 0., 0., 0.5 * a, 1., float(a)) for a in range(-10, 11)])
    net.create_phase_tap_changers(ptc, steps)
    net.update_phase_tap_changers(id=PST, regulation_value=target)
    if regulating:
        net.update_phase_tap_changers(id=PST, regulating=True)
    return net


def _olf(net, phase_control):
    from _olf_reference import reference_parameters
    params = reference_parameters(
        provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[SLACK_GEN, "bus_id"],
                  "newtonRaphsonConvEpsPerEq": "1e-12", "outerLoopNames": "PhaseControl" if phase_control else ""},
        distributed_slack=False, use_reactive_limits=False, read_slack_bus=False,
        phase_shifter_regulation_on=phase_control)
    res = lf.run_ac(net, params)[0]
    assert res.status.name == "CONVERGED", res.status_text
    return net.get_2_windings_transformers().loc[PST]


def _grid(net):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pypowsybl(net, gen_slack_id=SLACK_GEN, sort_index=False, buses_for_sub=False)


class TestPhaseControl(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    @staticmethod
    def _initial_flow():
        return _olf(_net("ACTIVE_POWER_CONTROL", 0., regulating=False), phase_control=False)

    def _compare(self, mode, target):
        net = _net(mode, target)
        olf = _olf(net, phase_control=True)
        olf_tap = net.get_phase_tap_changers(all_attributes=True).loc[PST, "solved_tap_position"]

        grid = _grid(_net(mode, target))
        grid.change_algorithm("NROuter_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(PhaseControl())
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-12)
        self.assertGreater(V.shape[0], 0)
        self.assertEqual(grid.get_algo().get_outer_loop_stats().status.name, "STABLE")
        self.assertEqual(grid.get_algo().get_linear_solver_stats().nb_analyze, 1)

        pst = next(el for el in grid.get_trafos() if el.name == PST)
        self.assertEqual(pst.res_phase_tap_position, olf_tap)
        self.assertEqual(pst.phase_tap_position, 0)  # the input is not modified
        self.assertAlmostEqual(pst.res_p1_mw, olf["p1"], places=8)
        olf_v = iidm_bus_voltages(net)
        ls = LightsimResultNetwork(grid, net).get_buses()
        vm = ls["v_mag"] / ls["voltage_level_id"].map(net.get_voltage_levels()["nominal_v"])
        common = olf_v.index.intersection(vm.dropna().index)
        self.assertLess((olf_v.loc[common, "vm_pu"] - vm.loc[common]).abs().max(), 1e-10)
        return olf_tap

    def test_active_power(self):
        p1 = self._initial_flow()["p1"]
        # rounded to an inside tap, up and down
        self.assertEqual(self._compare("ACTIVE_POWER_CONTROL", p1 + 7.), 3)
        self.assertEqual(self._compare("ACTIVE_POWER_CONTROL", p1 - 12.3), -5)
        # beyond the last tap
        self.assertEqual(self._compare("ACTIVE_POWER_CONTROL", p1 + 30.), 10)

    def test_current_limiter(self):
        i1 = self._initial_flow()["i1"]
        self.assertEqual(self._compare("CURRENT_LIMITER", 0.8 * i1), -5)
        self.assertEqual(self._compare("CURRENT_LIMITER", 0.95 * i1), -2)
        # below its limit: nothing moves
        self.assertEqual(self._compare("CURRENT_LIMITER", 1.1 * i1), 0)

    def test_detection(self):
        # a plain solve reports what the loop would act on, as a CONTROL
        flow = self._initial_flow()
        for mode, target, kind, value in (("ACTIVE_POWER_CONTROL", flow["p1"] + 7., LimitViolationType.PHASE_CONTROL_P,
                                           flow["p1"]),
                                          ("CURRENT_LIMITER", 0.8 * flow["i1"], LimitViolationType.PHASE_LIMITER_CURRENT,
                                           flow["i1"])):
            grid = _grid(_net(mode, target))
            grid.change_algorithm("NRSing_KLU")
            grid.clear_outer_loops()
            grid.add_outer_loop(PhaseControl())
            self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-12).shape[0], 0)
            viols = [v for v in grid.get_physical_violations(True, 1e-3, 0.) if v.violation_type == kind]
            self.assertEqual(len(viols), 1)
            v = viols[0]
            self.assertEqual(v.element_type, ViolationElementType.TRAFO)
            self.assertEqual(v.category, ViolationCategory.CONTROL)
            self.assertEqual(v.side, 1)
            self.assertAlmostEqual(v.value, value, places=6)
            self.assertAlmostEqual(v.limit, target, places=9)
        # within its limit: nothing
        grid = _grid(_net("CURRENT_LIMITER", 1.1 * flow["i1"]))
        grid.change_algorithm("NRSing_KLU")
        grid.clear_outer_loops()
        grid.add_outer_loop(PhaseControl())
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-12).shape[0], 0)
        self.assertFalse([v for v in grid.get_physical_violations(True, 1e-3, 0.)
                          if v.violation_type == LimitViolationType.PHASE_LIMITER_CURRENT])

    def test_not_in_the_default_list(self):
        grid = _grid(_net("ACTIVE_POWER_CONTROL", 0.))
        self.assertNotIn("PhaseControl", [loop.name() for loop in grid.get_outer_loops()])


if __name__ == "__main__":
    unittest.main()
