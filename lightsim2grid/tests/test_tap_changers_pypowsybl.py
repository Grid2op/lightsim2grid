# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The tap changers and shunt sections read from pypowsybl, against OpenLoadFlow: the pi model
at the taps, a tap or a section count moved on both sides, and what they regulate."""

import unittest
import warnings

import numpy as np

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
    from lightsim2grid.network.from_pypowsybl._result_network import LightsimResultNetwork
    from lightsim2grid.lightsim2grid_cpp import RegulationMode
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

SLACK_GEN = "GH1"


def _net():
    return pp.network.create_four_substations_node_breaker_network()


def _grid(net):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        return init_from_pypowsybl(net, gen_slack_id=SLACK_GEN, sort_index=False, buses_for_sub=False)


class TestTapChangersPypowsybl(unittest.TestCase):
    def setUp(self):
        if not HAS_PYPOWSYBL:
            self.skipTest("pypowsybl is not installed")

    @staticmethod
    def _trafo(grid, name):
        return next(el for el in grid.get_trafos() if el.name == name)

    def _max_dvm(self, net, grid):
        """OpenLoadFlow on `net` (no outer loop, slack on SLACK_GEN's bus) against `grid`, on
        the synchronous area of SLACK_GEN: the network has two, joined by hvdc lines, and the
        transformer is not in the largest one (OpenLoadFlow's main component)."""
        from _olf_reference import reference_parameters
        params = reference_parameters(
            provider={"slackBusSelectionMode": "NAME", "slackBusesIds": net.get_generators().loc[SLACK_GEN, "bus_id"],
                      "newtonRaphsonConvEpsPerEq": "1e-12", "outerLoopNames": ""},
            distributed_slack=False, use_reactive_limits=False, read_slack_bus=False,
            component_mode=lf.ComponentMode.ALL_CONNECTED)
        self.assertEqual(lf.run_ac(net, params)[0].status.name, "CONVERGED")
        grid.change_algorithm("NRSing_KLU")
        self.assertGreater(grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-12).shape[0], 0)
        olf = iidm_bus_voltages(net)
        ls = LightsimResultNetwork(grid, net).get_buses()
        vm = (ls["v_mag"] / ls["voltage_level_id"].map(net.get_voltage_levels()["nominal_v"])).dropna()
        vm = vm[vm > 0.]  # the other synchronous area, not solved here
        common = olf.dropna().index.intersection(vm.index)
        self.assertGreaterEqual(len(common), 2)
        return float((olf.loc[common, "vm_pu"] - vm.loc[common]).abs().max())

    def test_read(self):
        net = _net()
        grid = _grid(net)
        twt = self._trafo(grid, "TWT")
        rtc = net.get_ratio_tap_changers().loc["TWT"]
        ptc = net.get_phase_tap_changers(all_attributes=True).loc["TWT"]
        self.assertTrue(twt.has_ratio_tap_changer)
        self.assertEqual((twt.ratio_tap_position, twt.ratio_low_tap, twt.ratio_high_tap),
                         (rtc["tap"], rtc["low_tap"], rtc["high_tap"]))
        self.assertTrue(twt.has_phase_tap_changer)
        self.assertEqual(twt.phase_tap_position, ptc["tap"])
        # at the taps, the same pi model as before (pypowsybl's rho and alpha)
        net_pu = _net()
        net_pu.per_unit = True
        pu = net_pu.get_2_windings_transformers(all_attributes=True).loc["TWT"]
        self.assertAlmostEqual(twt.ratio, pu["rho"], places=12)
        self.assertAlmostEqual(twt.shift_rad, pu["alpha"], places=12)
        step = net.get_phase_tap_changer_steps().loc[("TWT", int(ptc["tap"]))]
        self.assertAlmostEqual(twt.x_pu, pu["x"] * (1. + step["x"] / 100.), places=12)
        # what they regulate
        self.assertEqual(twt.ratio_regulation_mode, RegulationMode.VOLTAGE)
        self.assertTrue(twt.ratio_regulating)
        self.assertAlmostEqual(twt.ratio_target, rtc["target_v"] / 225., places=12)
        self.assertGreaterEqual(twt.ratio_regulated, 0)
        self.assertEqual(twt.phase_regulation_mode, RegulationMode.CURRENT_LIMITER)
        self.assertFalse(twt.phase_regulating)
        self.assertEqual(twt.phase_regulated, 1)
        # the shunt's sections
        shunt = next(el for el in grid.get_shunts() if el.name == "SHUNT")
        self.assertTrue(shunt.has_sections)
        self.assertEqual((shunt.section_count, shunt.max_section_count), (1, 1))

    def test_same_as_olf(self):
        net = _net()
        self.assertLess(self._max_dvm(net, _grid(net)), 1e-9)

    def test_moved_taps_same_as_olf(self):
        # both changers moved, on OpenLoadFlow's network and on the grid built before
        net = _net()
        grid = _grid(net)
        tid = self._trafo(grid, "TWT").id
        net.update_ratio_tap_changers(id="TWT", tap=2, regulating=False)
        net.update_phase_tap_changers(id="TWT", tap=10)
        grid.change_trafo_ratio_tap(tid, 2)
        grid.change_trafo_phase_tap(tid, 10)
        self.assertLess(self._max_dvm(net, grid), 1e-9)
        self.assertEqual(grid.closest_trafo_phase_tap(tid, self._trafo(grid, "TWT").shift_rad), 10)

    def test_ratio_step_corrections(self):
        # a ratio tap step correcting r, x, g and b: OpenLoadFlow applies it, so does the grid
        net = _net()
        net.update_ratio_tap_changer_steps(id="TWT", position=1, r=10., x=-20., g=5., b=30.)
        self.assertLess(self._max_dvm(net, _grid(net)), 1e-9)

    def test_moved_sections_same_as_olf(self):
        net = _net()
        grid = _grid(net)
        sid = next(el.id for el in grid.get_shunts() if el.name == "SHUNT")
        net.update_shunt_compensators(id="SHUNT", section_count=0)
        grid.change_shunt_section_count(sid, 0)
        self.assertEqual(grid.get_shunts()[sid].target_q_mvar, 0.)
        self.assertLess(self._max_dvm(net, grid), 1e-9)

    def test_copy(self):
        grid = _grid(_net())
        tid = self._trafo(grid, "TWT").id
        grid.change_trafo_phase_tap(tid, 3)
        copy = grid.copy()
        self.assertEqual(self._trafo(copy, "TWT").phase_tap_position, 3)
        self.assertAlmostEqual(self._trafo(copy, "TWT").x_pu, self._trafo(grid, "TWT").x_pu, places=15)


if __name__ == "__main__":
    unittest.main()
