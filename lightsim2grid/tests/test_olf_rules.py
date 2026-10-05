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

import os
import sys
import unittest
import warnings

import numpy as np
import pandas as pd

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

    def test_curve_extrapolated_limits_crossing(self):
        # below the curve's first point, the extrapolated limits cross: both are their mean,
        # as powsybl-core (hence OpenLoadFlow) gives them
        net = _net()
        net.create_curve_reactive_limits(id=["B8-G", "B8-G"], p=[5., 50.], min_q=[-1., 4.], max_q=[-1., 5.])
        net.update_generators(id="B8-G", target_p=0.)
        lo, hi = _olf_rules.generator_limits_at_target_p(net, _olf_rules.generators(net).loc[["B8-G"]])
        # at 0 MW: min -1 - 5/45 * 5, max -1 - 6/45 * 5, crossed
        mean = ((-1. - 5. / 45. * 5.) + (-1. - 6. / 45. * 5.)) / 2.
        self.assertAlmostEqual(lo["B8-G"], mean)
        self.assertAlmostEqual(hi["B8-G"], mean)

    def test_vsc_station_target_p(self):
        # what OpenLoadFlow makes each station inject: its published P (load convention)
        net = pp.network.create_four_substations_node_breaker_network()
        target = _olf_rules.vsc_station_target_p(net)
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _olf_reference import reference_parameters
        lf.run_ac(net, reference_parameters(provider={"newtonRaphsonConvEpsPerEq": "1e-10"},
                                            distributed_slack=False, use_reactive_limits=False))
        published = net.get_vsc_converter_stations()["p"].dropna()  # not every station has one
        self.assertIn("VSC2", published.index)  # the inverter: the losses of both ends and the line
        for station in published.index:
            self.assertAlmostEqual(target[station], -published[station], places=6)

    def test_curve_limits_at_p(self):
        net = _net()
        net.create_curve_reactive_limits(id=["B8-G"] * 3, p=[0., 50., 100.], min_q=[-10., -20., -10.],
                                         max_q=[10., 30., 10.])
        lo, hi = _olf_rules.curve_limits_at_p(net, ["B8-G"], [25.], [np.nan], [np.nan])
        self.assertAlmostEqual(lo[0], -15.)  # interpolated inside the curve
        self.assertAlmostEqual(hi[0], 20.)
        lo, hi = _olf_rules.curve_limits_at_p(net, ["B8-G"], [120.], [np.nan], [np.nan])
        self.assertAlmostEqual(lo[0], -6.)   # extrapolated past its end
        self.assertAlmostEqual(hi[0], 2.)
        lo, hi = _olf_rules.curve_limits_at_p(net, ["B8-G"], [150.], [np.nan], [np.nan])
        self.assertAlmostEqual(lo[0], -5.)   # crossed (0 and -10): both their mean
        self.assertAlmostEqual(hi[0], -5.)

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

    def test_participation(self):
        net = _net()
        gen = net.get_generators()
        w = _olf_rules.generator_participation_weight(net, gen)
        # without the extension: every unit, with the default droop
        self.assertTrue(np.allclose(w.to_numpy(), gen["max_p"].to_numpy() / 4.))
        p = OlfLoadingParameters()
        # one condition of checkActivePowerControl at a time
        cases = {
            "zero target": dict(target_p=[0.], min_p=[-10.], max_p=[100.]),
            "implausible max_p": dict(target_p=[50.], min_p=[0.], max_p=[20000.]),
            "above its range": dict(target_p=[120.], min_p=[0.], max_p=[100.]),
            "below its range": dict(target_p=[5.], min_p=[10.], max_p=[100.]),
            "degenerate range": dict(target_p=[50.], min_p=[50.], max_p=[50.]),
        }
        for label, c in cases.items():
            w = _olf_rules.participation_weight(c["target_p"], c["min_p"], c["max_p"], [True], [np.nan],
                                                [np.nan], [np.nan], p)
            self.assertEqual(w[0], 0., label)
        # a negative target takes part (a pump), with the same key
        w = _olf_rules.participation_weight([-30.], [-100.], [100.], [True], [np.nan], [np.nan], [np.nan], p)
        self.assertAlmostEqual(w[0], 25.)
        # the extension: its flag, its droop (0 excludes), its target range
        w = _olf_rules.participation_weight([50.] * 4, [0.] * 4, [100.] * 4, [False, True, True, True],
                                            [np.nan, 0., 2., np.nan], [np.nan, np.nan, np.nan, 60.],
                                            [np.nan] * 4, p)
        self.assertTrue(np.allclose(w, [0., 0., 50., 0.]))
        # without active limits, the range is not looked at
        w = _olf_rules.participation_weight([120.], [0.], [100.], [True], [np.nan], [np.nan], [np.nan],
                                            OlfLoadingParameters(use_active_limits=False))
        self.assertAlmostEqual(w[0], 25.)

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
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _olf_reference import reference_parameters
        wanted = {"slackBusSelectionMode": "NAME", "slackBusesIds": slack_bus,
                  "outerLoopNames": "", "newtonRaphsonConvEpsPerEq": "1e-10",
                  "generatorsWithZeroMwTargetAreNotStarted": "true",
                  "disableInconsistentVoltageControls": "true"}
        known = set(lf.get_provider_parameters_names())
        missing = {k for k in wanted if k not in known}
        if missing:
            raise unittest.SkipTest(f"this pypowsybl's OpenLoadFlow lacks {sorted(missing)}")
        return reference_parameters(provider=wanted, distributed_slack=False, use_reactive_limits=False,
                                    read_slack_bus=False)

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

    def test_branches_on_same_bus_as_olf(self):
        net = _net_with_self_loops()
        slack_bus = net.get_generators().loc["B1-G", "bus_id"]
        res = lf.run_ac(net, self._olf_params(slack_bus))
        self.assertEqual(res[0].status, lf.ComponentStatus.CONVERGED)
        olf_vm = iidm_bus_voltages(net)["vm_pu"]

        def solve(olf_rules):
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                grid = init_from_pypowsybl(_net_with_self_loops(), gen_slack_id="B1-G", sort_index=False,
                                           buses_for_sub=False, olf_rules=olf_rules)
            V = grid.ac_pf(np.full(grid.total_bus(), 1.0 + 0j), 30, 1e-10)
            self.assertGreater(V.shape[0], 0)
            to_iidm = lightsim_bus_to_iidm(grid, net)
            return max(abs(abs(V[i]) - olf_vm[to_iidm[i]]) for i in range(V.shape[0])
                       if i in to_iidm and to_iidm[i] in olf_vm.index)

        self.assertLess(solve(OlfLoadingParameters(reactive_limits=False)), 1e-8)
        # kept, the phase shifter carries a flow around itself and the voltages move
        self.assertGreater(solve(False), 1e-4)


def _net_with_self_loops():
    """``_net`` plus a phase-shifting transformer and a line, each with both ends on the
    same bus (and a branch between two buses, to check it is left alone)"""
    net = _net()
    vl = "VL4"
    bus = net.get_bus_breaker_view_buses().query("voltage_level_id == @vl").index[0]
    net.create_2_windings_transformers(id="T-self", voltage_level1_id=vl, bus1_id=bus, voltage_level2_id=vl,
                                       bus2_id=bus, rated_u1=135., rated_u2=135., r=0.5, x=10., g=1e-6, b=-1e-5)
    net.create_phase_tap_changers(
        pd.DataFrame.from_records(index="id", data=[{"id": "T-self", "target_deadband": 0.,
                                                     "regulation_mode": "CURRENT_LIMITER", "low_tap": 0, "tap": 0}]),
        pd.DataFrame.from_records(index="id", data=[{"id": "T-self", "b": 0., "g": 0., "r": 0., "x": 0.,
                                                     "rho": 1., "alpha": 3.}]))
    net.create_lines(id="L-self", voltage_level1_id=vl, bus1_id=bus, voltage_level2_id=vl, bus2_id=bus,
                     r=1., x=10., g1=0., b1=1e-3, g2=0., b2=1e-3)
    return net


@unittest.skipIf(not PP_OK, "pypowsybl is not installed")
class TestOlfBranchOnSameBus(unittest.TestCase):
    def test_rule(self):
        net = _net_with_self_loops()
        lines = net.get_lines()
        trafos = net.get_2_windings_transformers()
        self.assertEqual(list(lines.index[_olf_rules.branch_on_same_bus(lines)]), ["L-self"])
        self.assertEqual(list(trafos.index[_olf_rules.branch_on_same_bus(trafos)]), ["T-self"])
        # one end open: no longer on the same bus
        net.update_lines(id="L-self", connected2=False)
        self.assertFalse(_olf_rules.branch_on_same_bus(net.get_lines()).any())

    def test_init_disconnects_them(self):
        net = _net_with_self_loops()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            raw = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=False)
            ruled = init_from_pypowsybl(net, gen_slack_id="B1-G", sort_index=False, olf_rules=True)
        line_ids = list(net.get_lines().index)
        trafo_ids = list(net.get_2_windings_transformers().index)
        for grid, expected in ((raw, True), (ruled, False)):
            lines = grid.get_lines()
            trafos = grid.get_trafos()
            self.assertEqual(lines[line_ids.index("L-self")].connected_global, expected)
            self.assertEqual(trafos[trafo_ids.index("T-self")].connected_global, expected)
            # every other branch untouched
            self.assertTrue(all(el.connected_global for i, el in enumerate(lines) if line_ids[i] != "L-self"))
            self.assertTrue(all(el.connected_global for i, el in enumerate(trafos) if trafo_ids[i] != "T-self"))


def _net_with_monitors():
    """IEEE 14 with stand-by SVCs: alone on a load bus (SVC9), next to a regulating
    generator (SVC8), and two on one bus (SVC10a, SVC10b)."""
    net = _net()
    nominal_v = net.get_voltage_levels()["nominal_v"]
    for svc_id, bus in (("SVC9", "B9"), ("SVC8", "B8"), ("SVC10a", "B10"), ("SVC10b", "B10")):
        vn = float(nominal_v.loc["VL" + bus[1:]])
        net.create_static_var_compensators(id=svc_id, voltage_level_id="VL" + bus[1:], bus_id=bus,
                                           connectable_bus_id=bus, b_min=-0.01, b_max=0.01,
                                           regulation_mode="VOLTAGE", target_v=vn, target_q=0.,
                                           regulating=True)
        net.create_extensions("standbyAutomaton", id=svc_id, b0=0., standby=True,
                              low_voltage_threshold=0.9 * vn, high_voltage_threshold=1.1 * vn,
                              low_voltage_setpoint=0.95 * vn, high_voltage_setpoint=1.05 * vn)
    return net


@unittest.skipIf(not PP_OK, "pypowsybl is not installed")
class TestOlfVoltageControllers(unittest.TestCase):
    """The voltage-control rules over every kind of unit (``voltage_controllers``)."""

    def test_monitors(self):
        vc = _olf_rules.voltage_controllers(_net_with_monitors())
        self.assertTrue(vc.loc["SVC9", "monitor"])
        self.assertFalse(vc.loc["SVC9", "discarded"])
        # next to a regulating unit, the monitor is switched off; the unit keeps regulating
        self.assertTrue(vc.loc["SVC8", "monitor_with_regulator"])
        self.assertFalse(vc.loc["SVC8", "monitor"])
        self.assertFalse(vc.loc["B8-G", "discarded"])
        # two on a bus: both regulate
        for svc_id in ("SVC10a", "SVC10b"):
            self.assertFalse(vc.loc[svc_id, "monitor"])
            self.assertFalse(vc.loc[svc_id, "discarded"])
        # without svcVoltageMonitoring, a standby SVC is an ordinary regulator
        vc = _olf_rules.voltage_controllers(_net_with_monitors(),
                                            OlfLoadingParameters(svc_voltage_monitoring=False))
        self.assertFalse(vc["monitor"].any())
        self.assertFalse(vc.loc["SVC9", "discarded"])
        # SVC8 then regulates its bus with the generator there, at another target
        self.assertTrue(vc.loc["SVC8", "inconsistent_controls"])

    def test_monitors_in_the_grid(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(_net_with_monitors(), olf_rules=True)
        svcs = {svc.name: svc for svc in grid.get_svcs()}
        # idle (off) and flagged standby with its thresholds, in pu of its own nominal voltage
        self.assertEqual(svcs["SVC9"].regulation_mode, 0)
        self.assertTrue(svcs["SVC9"].standby)
        self.assertAlmostEqual(svcs["SVC9"].standby_low_vm_pu, 0.9)
        self.assertAlmostEqual(svcs["SVC9"].standby_high_vm_pu, 1.1)
        self.assertEqual(svcs["SVC8"].regulation_mode, 0)
        self.assertFalse(svcs["SVC8"].standby)
        self.assertEqual(svcs["SVC10a"].regulation_mode, 1)

    def test_svc_reactive_range(self):
        net = _net_with_monitors()
        # (b_max - b_min) * V^2 below 1 MVar at the voltage the snapshot stores
        net.update_static_var_compensators(id="SVC10a", b_min=-1e-6, b_max=1e-6)
        vc = _olf_rules.voltage_controllers(net, OlfLoadingParameters(svc_voltage_monitoring=False))
        self.assertTrue(vc.loc["SVC10a", "reactive_range_too_small"])
        self.assertFalse(vc.loc["SVC10b", "discarded"])

    def test_kinds_judged_together(self):
        net = pp.network.create_four_substations_node_breaker_network()
        vc = _olf_rules.voltage_controllers(net)
        self.assertEqual(set(vc["kind"]), {"generator", "svc", "vsc"})
        self.assertFalse(vc["discarded"].any())
        # a VSC station regulating its bus with a target 0.025 pu away from the generators'
        net.update_vsc_converter_stations(id="VSC1", target_v=410.)
        vc = _olf_rules.voltage_controllers(net)
        for unit in ("GH1", "GH2", "GH3", "VSC1"):
            self.assertTrue(vc.loc[unit, "inconsistent_controls"], unit)
        self.assertTrue(_olf_rules.generator_voltage_control(net).loc["GH1", "discarded"])
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            grid = init_from_pypowsybl(net, olf_rules=True)
        gens = {gen.name: gen for gen in grid.get_generators()}
        self.assertFalse(gens["GH1"].voltage_regulator_on)
        self.assertTrue(gens["GTH2"].voltage_regulator_on)


if __name__ == "__main__":
    unittest.main()
