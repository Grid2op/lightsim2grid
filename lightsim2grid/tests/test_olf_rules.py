# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's network-loading rules (``from_pypowsybl/_olf_rules.py``), one rule at a
time on pypowsybl's IEEE 14-bus network, then end to end: a grid built with
``init_from_pypowsybl(olf_rules=True)`` solves to OpenLoadFlow's answer where the raw grid
does not."""

import unittest
import warnings

import numpy as np

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl, OlfLoadingParameters
    from lightsim2grid.network.from_pypowsybl import _olf_rules
    from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages, lightsim_bus_to_iidm
    PP_OK = True
except ImportError:
    PP_OK = False


def _net():
    net = pp.network.create_ieee14()
    # every unit regulates locally, with a positive minimum P and dispatched above it (the
    # three synchronous condensers of the case, at 0 MW, would otherwise be "not started")
    gen = net.get_generators()
    target_p = gen["target_p"].where(gen["target_p"] > 0., 10.)
    net.update_generators(id=list(gen.index), min_p=[1.] * len(gen), max_p=[500.] * len(gen),
                          target_p=list(target_p))
    return net


def _discards(net, params=None):
    params = OlfLoadingParameters() if params is None else params
    return _olf_rules.generator_voltage_control(net, params=params)


@unittest.skipIf(not PP_OK, "pypowsybl is not installed")
class TestOlfRules(unittest.TestCase):
    def test_nothing_discarded_on_the_plain_grid(self):
        self.assertFalse(_discards(_net())["discarded"].any())

    def test_not_started(self):
        net = _net()
        # 0.005 MW is "zero" for OpenLoadFlow's voltage control (its epsilon is per unit)
        net.update_generators(id=["B3-G", "B6-G"], target_p=[0., 0.005])
        net.update_generators(id="B8-G", target_p=0., min_p=0.)  # allowed to sit at 0 MW
        res = _discards(net)
        self.assertTrue(res.loc["B3-G", "not_started"])
        self.assertTrue(res.loc["B6-G", "not_started"])
        self.assertFalse(res.loc["B8-G", "discarded"])
        # a condenser and (FORCED mode) a fictitious unit are never "not started" (the rules
        # are functions of the generator frame: pypowsybl does not let `condenser` be edited)
        gen = _olf_rules.generators(net)
        gen.loc["B3-G", "condenser"] = True
        gen.loc["B6-G", "fictitious"] = True
        res = _olf_rules.generator_voltage_control(net, gen=gen)
        self.assertFalse(res.loc["B3-G", "discarded"])
        self.assertFalse(res.loc["B6-G", "discarded"])
        res = _olf_rules.generator_voltage_control(
            net, gen=gen, params=OlfLoadingParameters(fictitious_voltage_control_forced=False))
        self.assertTrue(res.loc["B6-G", "discarded"])
        # and the rule can be switched off
        self.assertFalse(_discards(net, OlfLoadingParameters(zero_mw_target_not_started=False)).loc["B3-G", "discarded"])

    def test_reactive_range(self):
        net = _net()
        net.update_generators(id="B8-G", min_q=-0.4, max_q=0.4)  # a 0.8 MVar range
        self.assertTrue(_discards(net).loc["B8-G", "reactive_range_too_small"])
        # without reactive limits, the range is not looked at
        self.assertFalse(_discards(net, OlfLoadingParameters(reactive_limits=False)).loc["B8-G", "discarded"])

    def test_implausible_target_v(self):
        net = _net()
        nominal = net.get_voltage_levels().loc[net.get_generators().loc["B2-G", "voltage_level_id"], "nominal_v"]
        net.update_generators(id="B2-G", target_v=1.25 * nominal)
        self.assertTrue(_discards(net).loc["B2-G", "target_v_implausible"])
        # no check below the nominal voltage threshold
        self.assertFalse(_discards(net, OlfLoadingParameters(
            min_nominal_voltage_target_voltage_check=nominal + 1.)).loc["B2-G", "discarded"])

    def test_inconsistent_controls(self):
        net = _net()
        gen = net.get_generators()
        vl = gen.loc["B2-G", "voltage_level_id"]
        bus = net.get_bus_breaker_view_buses().query("voltage_level_id == @vl").index[0]
        target_v = gen.loc["B2-G", "target_v"]
        nominal = net.get_voltage_levels().loc[vl, "nominal_v"]
        net.create_generators(id="B2-G2", voltage_level_id=vl, bus_id=bus, target_p=10., min_p=1., max_p=50.,
                              target_v=target_v + 0.005 * nominal, voltage_regulator_on=True, target_q=0.)
        net.update_generators(id="B2-G2", min_q=-50., max_q=50.)
        # targets 0.005 pu apart: consistent
        self.assertFalse(_discards(net)["discarded"].any())
        # 0.02 pu apart: both lose their voltage control
        net.update_generators(id="B2-G2", target_v=target_v + 0.02 * nominal)
        res = _discards(net)
        self.assertTrue(res.loc["B2-G", "inconsistent_controls"])
        self.assertTrue(res.loc["B2-G2", "inconsistent_controls"])
        self.assertFalse(_discards(net, OlfLoadingParameters(
            disable_inconsistent_voltage_controls=False))["discarded"].any())
        # a unit discarded for another reason no longer counts on the bus
        net.update_generators(id="B2-G2", target_p=0.)
        res = _discards(net)
        self.assertTrue(res.loc["B2-G2", "not_started"])
        self.assertFalse(res.loc["B2-G", "discarded"])

    def test_target_q_clamped(self):
        net = _net()
        net.update_generators(id="B8-G", voltage_regulator_on=False, target_q=80., min_q=-6., max_q=24.)
        tq = _olf_rules.generator_target_q(net)
        self.assertAlmostEqual(tq.loc["B8-G"], 24.)
        tq = _olf_rules.generator_target_q(net, params=OlfLoadingParameters(force_target_q_in_reactive_limits=False))
        self.assertAlmostEqual(tq.loc["B8-G"], 80.)

    def test_init_applies_the_rules(self):
        net = _net()
        net.update_generators(id="B3-G", target_p=0.)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            raw = init_from_pypowsybl(net, gen_slack_id="B1-G")
            ruled = init_from_pypowsybl(net, gen_slack_id="B1-G", olf_rules=True)
        gen_ids = list(net.get_generators().index)
        i = gen_ids.index("B3-G")
        self.assertTrue(raw.get_generators()[i].voltage_regulator_on)
        self.assertFalse(ruled.get_generators()[i].voltage_regulator_on)
        with self.assertRaises(TypeError):
            init_from_pypowsybl(net, gen_slack_id="B1-G", olf_rules="yes")


@unittest.skipIf(not PP_OK, "pypowsybl is not installed")
class TestOlfRulesAgainstOLF(unittest.TestCase):
    """The rules, end to end: OpenLoadFlow run without any outer loop, its loading rules
    explicitly on (they differ between pypowsybl builds), and lightsim2grid built with the
    same rules land on the same voltages."""

    @staticmethod
    def _modified_net():
        net = _net()
        # a not-started unit, and two units of one bus with inconsistent targets
        net.update_generators(id="B3-G", target_p=0.)
        gen = net.get_generators()
        vl = gen.loc["B6-G", "voltage_level_id"]
        bus = net.get_bus_breaker_view_buses().query("voltage_level_id == @vl").index[0]
        nominal = net.get_voltage_levels().loc[vl, "nominal_v"]
        net.create_generators(id="B6-G2", voltage_level_id=vl, bus_id=bus, target_p=5., min_p=1., max_p=50.,
                              target_v=gen.loc["B6-G", "target_v"] + 0.03 * nominal,
                              voltage_regulator_on=True, target_q=2.)
        return net

    @staticmethod
    def _olf_params(slack_bus):
        params = lf.Parameters(distributed_slack=False, use_reactive_limits=False,
                               read_slack_bus=False, twt_split_shunt_admittance=True)
        wanted = {"slackBusSelectionMode": "NAME", "slackBusesIds": slack_bus,
                  "outerLoopNames": "", "newtonRaphsonConvEpsPerEq": "1e-10",
                  "generatorsWithZeroMwTargetAreNotStarted": "true",
                  "disableInconsistentVoltageControls": "true"}
        known = set(lf.get_provider_parameters_names())
        missing = {k for k in wanted if k not in known}
        if missing:
            raise unittest.SkipTest(f"this pypowsybl's OpenLoadFlow lacks {sorted(missing)}")
        params.provider_parameters = wanted
        return params

    def _solve_ls(self, olf_rules):
        net = self._modified_net()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                       olf_rules=olf_rules)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
        self.assertGreater(V.shape[0], 0)
        to_iidm = lightsim_bus_to_iidm(grid, net)
        return {to_iidm[i]: abs(V[i]) for i in range(V.shape[0]) if i in to_iidm}

    def test_same_voltages_as_olf(self):
        net = self._modified_net()
        slack_bus = net.get_generators().loc["B1-G", "bus_id"]
        res = lf.run_ac(net, self._olf_params(slack_bus))
        self.assertEqual(res[0].status, lf.ComponentStatus.CONVERGED)
        olf_vm = iidm_bus_voltages(net)["vm_pu"]

        ruled = self._solve_ls(OlfLoadingParameters(reactive_limits=False))
        dvm = max(abs(ruled[b] - olf_vm[b]) for b in ruled if b in olf_vm.index)
        self.assertLess(dvm, 1e-8)

        # read as written, the grid has two set-points on one bus, which lightsim2grid refuses
        with self.assertRaises(RuntimeError):
            self._solve_ls(False)


if __name__ == "__main__":
    unittest.main()
