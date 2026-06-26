# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Phase 4 reference matrix: remote generator voltage control + Static Var
Compensators (SVC) through the bordered ``VoltageControl`` NR extension, validated
against PowSyBl OpenLoadFlow (OLF) to 1e-6 through the pypowsybl converter.

OLF's slack is pinned (by voltage level NAME) to the lightsim2grid slack so the
angle reference matches; reactive limits stay OFF (the feature has no Q-limit outer
loop) and the remote-control / slope provider parameters are enabled. Skipped when
pypowsybl is unavailable.
"""

import unittest
import numpy as np

try:
    import pypowsybl as pp
    from lightsim2grid.network import init_from_pypowsybl
    from lightsim2grid.network.from_pypowsybl import bake_outer_loops
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

SN_MVA = 100.0
# SvcContainer.RegulationMode codes mirrored from C++
OFF, VOLTAGE, REACTIVE_POWER = 0, 1, 2


def _vl(ids, nominal_v=400.0):
    return dict(id=ids, substation_id=[f"S{i}" for i in range(len(ids))],
                topology_kind=["BUS_BREAKER"] * len(ids), nominal_v=[nominal_v] * len(ids))


def _star(extra_q_g1=None, extra_q_g2=None, g1_reg="LD", g2_reg=None,
          g1_tv=405.0, g2_tv=405.0):
    """A 4-bus star: slack gen G0 at B0, controllers G1@B1 / (optional) G2@B3, all
    radial to the load bus B2. Returns the network."""
    n = pp.network.create_empty()
    subs = ["S0", "S1", "S2"] + (["S3"] if g2_reg is not None else [])
    vls = ["VL0", "VL1", "VL2"] + (["VL3"] if g2_reg is not None else [])
    n.create_substations(id=subs)
    n.create_voltage_levels(id=vls, substation_id=subs,
                            topology_kind=["BUS_BREAKER"] * len(vls), nominal_v=[400.0] * len(vls))
    buses = ["B0", "B1", "B2"] + (["B3"] if g2_reg is not None else [])
    n.create_buses(id=buses, voltage_level_id=vls)
    lines_v1 = ["VL0", "VL1"] + (["VL3"] if g2_reg is not None else [])
    lines_b1 = ["B0", "B1"] + (["B3"] if g2_reg is not None else [])
    lid = ["L02", "L12"] + (["L32"] if g2_reg is not None else [])
    nl = len(lid)
    n.create_lines(id=lid, voltage_level1_id=lines_v1, voltage_level2_id=["VL2"] * nl,
                   bus1_id=lines_b1, bus2_id=["B2"] * nl,
                   r=[1.0] * nl, x=[20.0] * nl, g1=[0.0] * nl, b1=[0.0] * nl, g2=[0.0] * nl, b2=[0.0] * nl)
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=120.0, q0=60.0)
    n.create_generators(id="G0", voltage_level_id="VL0", bus_id="B0", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=0.0, target_v=g1_tv, voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    if g2_reg is not None:
        n.create_generators(id="G2", voltage_level_id="VL3", bus_id="B3", target_p=60.0,
                            target_q=0.0, target_v=g2_tv, voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    if extra_q_g1 is not None:
        ids = ["G1"] + (["G2"] if g2_reg is not None else [])
        mins = [extra_q_g1[0]] + ([extra_q_g2[0]] if g2_reg is not None else [])
        maxs = [extra_q_g1[1]] + ([extra_q_g2[1]] if g2_reg is not None else [])
        n.create_minmax_reactive_limits(id=ids, min_q=mins, max_q=maxs)
    reg = {"G1": g1_reg}
    if g2_reg is not None:
        reg["G2"] = g2_reg
    n.update_generators(id=list(reg.keys()), regulated_element_id=list(reg.values()))
    return n


def _svc_net(mode=VOLTAGE, target_v=398.0, target_q=0.0, slope=0.0, remote=False, with_svc=True):
    """3 radial buses: slack gen G@B0 -- B1 (SVC) -- B2 (load). The SVC regulates its
    own bus B1 (local) or the remote PQ bus B2; B0 (slack) is never the SVC target."""
    n = pp.network.create_empty()
    n.create_substations(id=["S0", "S1", "S2"])
    n.create_voltage_levels(id=["VL0", "VL1", "VL2"], substation_id=["S0", "S1", "S2"],
                            topology_kind=["BUS_BREAKER"] * 3, nominal_v=[400.0] * 3)
    n.create_buses(id=["B0", "B1", "B2"], voltage_level_id=["VL0", "VL1", "VL2"])
    n.create_lines(id=["L01", "L12"], voltage_level1_id=["VL0", "VL1"], voltage_level2_id=["VL1", "VL2"],
                   bus1_id=["B0", "B1"], bus2_id=["B1", "B2"],
                   r=[1.0, 1.0], x=[20.0, 20.0], g1=[0.0, 0.0], b1=[0.0, 0.0], g2=[0.0, 0.0], b2=[0.0, 0.0])
    n.create_generators(id="G", voltage_level_id="VL0", bus_id="B0", target_p=120.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_loads(id=["LD1", "LD2"], voltage_level_id=["VL1", "VL2"],
                   bus_id=["B1", "B2"], p0=[40.0, 60.0], q0=[20.0, 30.0])
    if with_svc:
        mode_str = {VOLTAGE: "VOLTAGE", REACTIVE_POWER: "REACTIVE_POWER", OFF: "VOLTAGE"}[mode]
        n.create_static_var_compensators(id="SVC1", voltage_level_id="VL1", bus_id="B1",
                                         regulation_mode=mode_str, regulating=(mode != OFF),
                                         target_v=target_v, target_q=target_q, b_min=-1.0, b_max=1.0)
        if remote:
            n.update_static_var_compensators(id="SVC1", regulated_element_id="LD2")  # regulate B2
        if slope != 0.0:
            n.create_extensions("voltagePerReactivePowerControl", id="SVC1", slope=slope)
    return n


def _dist_slack_pq_net():
    """3 radial buses for a distributed slack with a PQ participant: G0@B0 is a
    voltage-regulating gen, G1@B1 is a NON-regulating gen (voltage_regulator_on=
    False, e.g. baked to PQ at its Q-limit). Both join the distributed slack, so
    G1's bus is a slack bus that is NOT pinned by a local PV gen: it must keep a
    free Vm unknown + Q equation, otherwise its magnitude is frozen at the init
    value. The load exceeds generation so the slack actually distributes."""
    n = pp.network.create_empty()
    n.create_substations(id=["S0", "S1", "S2"])
    n.create_voltage_levels(id=["VL0", "VL1", "VL2"], substation_id=["S0", "S1", "S2"],
                            topology_kind=["BUS_BREAKER"] * 3, nominal_v=[400.0] * 3)
    n.create_buses(id=["B0", "B1", "B2"], voltage_level_id=["VL0", "VL1", "VL2"])
    n.create_lines(id=["L02", "L12"], voltage_level1_id=["VL0", "VL1"], voltage_level2_id=["VL2", "VL2"],
                   bus1_id=["B0", "B1"], bus2_id=["B2", "B2"],
                   r=[1.0, 1.0], x=[20.0, 20.0], g1=[0.0, 0.0], b1=[0.0, 0.0], g2=[0.0, 0.0], b2=[0.0, 0.0])
    n.create_loads(id="LD", voltage_level_id="VL2", bus_id="B2", p0=180.0, q0=60.0)
    n.create_generators(id="G0", voltage_level_id="VL0", bus_id="B0", target_p=60.0,
                        target_q=0.0, target_v=400.0, voltage_regulator_on=True, max_p=1000.0, min_p=0.0)
    n.create_generators(id="G1", voltage_level_id="VL1", bus_id="B1", target_p=60.0,
                        target_q=10.0, target_v=400.0, voltage_regulator_on=False, max_p=1000.0, min_p=0.0)
    return n


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
class TestVoltageControlPypowsybl(unittest.TestCase):
    def setUp(self):
        self.tol = 1e-6

    def _params(self, slack_vl, **extra):
        prov = {"slackBusSelectionMode": "NAME", "slackBusesIds": slack_vl,
                "voltageRemoteControl": "true",
                "generatorReactivePowerRemoteControl": "true",
                "voltagePerReactivePowerControl": "true"}
        prov.update(extra)
        return pp.loadflow.Parameters(distributed_slack=False, use_reactive_limits=False,
                                      provider_parameters=prov)

    def _run_ls(self, n, gen_slack_id, dc=False):
        n.per_unit = False
        model = init_from_pypowsybl(n, gen_slack_id=gen_slack_id, sort_index=True)
        nb = len(model.get_bus_status())
        V0 = np.ones(nb, dtype=np.complex128)
        V = (model.dc_pf if dc else model.ac_pf)(V0, 30, 1e-11)
        self.assertGreater(V.shape[0], 0, "lightsim2grid diverged")
        return model, V

    def _assert_bus_parity(self, n, V, name, check_mag=True, check_angle=True):
        n.per_unit = False
        buses = n.get_buses()
        vls = n.get_voltage_levels()
        for i, (bid, row) in enumerate(buses.iterrows()):
            vn = vls.loc[row["voltage_level_id"], "nominal_v"]
            if check_mag:
                self.assertLess(abs(abs(V[i]) - row["v_mag"] / vn), self.tol, f"{name}: |V| at {bid}")
            if check_angle:
                self.assertLess(abs(np.angle(V[i]) - np.deg2rad(row["v_angle"])), self.tol,
                                f"{name}: angle at {bid}")

    def _compare(self, n, slack_vl, gen_slack_id, name, **extra):
        ref = pp.loadflow.run_ac(n, parameters=self._params(slack_vl, **extra))
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED, name)
        model, V = self._run_ls(n, gen_slack_id)
        self._assert_bus_parity(n, V, name)
        return model

    # ----- remote generators ------------------------------------------------
    def test_remote_gen_single(self):
        n = _star(g2_reg=None, g1_reg="LD", g1_tv=405.0)
        self._compare(n, "VL0", "G0", "remote-gen-single")

    def test_remote_gen_two_share(self):
        # two gens, different reactive ranges, regulating the SAME remote bus B2.
        # OLF shares Q proportionally to (qmax-qmin) under Q_EQUAL_PROPORTION.
        n = _star(extra_q_g1=(-50.0, 50.0), extra_q_g2=(-100.0, 200.0),
                  g1_reg="LD", g2_reg="LD", g1_tv=405.0, g2_tv=405.0)
        model = self._compare(n, "VL0", "G0", "remote-gen-share",
                              reactivePowerDispatchMode="Q_EQUAL_PROPORTION")
        # reactive sharing ratio Q1/Q2 must equal w1/w2 = 100/300
        q = {g.name: g.res_q_mvar for g in model.get_generators()}
        self.assertAlmostEqual(q["G1"] / q["G2"], 100.0 / 300.0, places=4)

    def test_remote_gen_on_slack(self):
        # A remote controller whose OWN bus is the slack. The active-power slack
        # role and the reactive/voltage role are decoupled: the slack bus is given
        # a free Vm unknown + a Q equation (MultiSlack) so the remote control
        # attaches, while its angle stays the reference. (This used to be rejected,
        # cf the former test_remote_gen_on_pv_bus_raises.) Matches OLF when OLF's
        # slack is pinned to the same voltage level.
        n = _star(g2_reg=None, g1_reg="LD", g1_tv=405.0)
        self._compare(n, "VL1", "G1", "remote-gen-on-slack")

    def test_remote_gen_on_slack_baked(self):
        # The same controller-on-slack case that test_remote_gen_on_pv_bus_raises
        # rejects: bake_outer_loops rewrites G1's remote control to local control at
        # its solved terminal voltage, so it loads and reproduces the (remote-control)
        # OLF solution. OLF's slack is pinned to G1's voltage level so the angle
        # references coincide; baking does not re-solve, so the network still carries
        # the with-loops solution that _assert_bus_parity compares against.
        n = _star(g2_reg=None, g1_reg="LD", g1_tv=405.0)
        n.per_unit = False
        ref = pp.loadflow.run_ac(n, parameters=self._params("VL1"))
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        bake_outer_loops(n)
        # G1 is now a LOCAL controller of its own terminal
        self.assertEqual(n.get_generators().loc["G1", "regulated_element_id"], "G1")
        _, V = self._run_ls(n, "G1")  # used to raise; now supported via baking
        self._assert_bus_parity(n, V, "remote-gen-on-slack-baked")

    def test_dist_slack_pq_participant(self):
        # Regression: a distributed-slack participant whose generator is PQ
        # (voltage_regulator_on=False) gets the active-power slack role on its P
        # row, but its magnitude must still be solved through a Q equation -- it is
        # NOT pinned by a local PV gen. Before MultiSlack handled this, the bus was
        # Vm-fixed at the init value (1.0), which on a real grid froze the magnitude
        # of every non-regulating distributed-slack gen. Compare magnitudes to OLF's
        # distributed slack (proportional to Pmax), with G0's voltage level as the
        # angle reference.
        n = _dist_slack_pq_net()
        n.per_unit = False
        params = pp.loadflow.Parameters(
            distributed_slack=True, use_reactive_limits=False,
            balance_type=pp.loadflow.BalanceType.PROPORTIONAL_TO_GENERATION_P_MAX,
            provider_parameters={"slackBusSelectionMode": "NAME", "slackBusesIds": "VL0"})
        ref = pp.loadflow.run_ac(n, parameters=params)
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        self.assertGreater(abs(ref[0].distributed_active_power), 1.0)  # slack really distributes

        model, V = self._run_ls(n, {"G0": 1.0, "G1": 1.0})
        # magnitude parity with OLF (the frozen bug showed up on |V|; the angle
        # datum of a distributed slack is reference-dependent, so check |V| only)
        self._assert_bus_parity(n, V, "dist-slack-pq", check_angle=False)
        # G1's bus (VL1) must be solved away from the frozen init value, and G1
        # keeps its fixed reactive setpoint (PQ behaviour)
        vl1 = list(n.get_buses().index).index(
            n.get_buses().index[n.get_buses()["voltage_level_id"] == "VL1"][0])
        self.assertGreater(abs(abs(V[vl1]) - 1.0), 1e-3, "VL1 magnitude was frozen at init")
        q = {g.name: g.res_q_mvar for g in model.get_generators()}
        self.assertAlmostEqual(q["G1"], 10.0, places=4)

    # ----- SVC --------------------------------------------------------------
    def test_svc_local_voltage(self):
        self._compare(_svc_net(mode=VOLTAGE, target_v=398.0), "VL0", "G", "svc-local-V")

    def test_svc_remote_voltage(self):
        self._compare(_svc_net(mode=VOLTAGE, target_v=402.0, remote=True),
                      "VL0", "G", "svc-remote-V")

    def test_svc_local_slope(self):
        self._compare(_svc_net(mode=VOLTAGE, target_v=398.0, slope=0.05),
                      "VL0", "G", "svc-local-slope")

    def test_svc_remote_slope(self):
        # Remote SVC voltage control WITH a slope is underspecified in OLF: OLF holds
        # the regulated bus at exactly vset (slope ignored for remote control), so we
        # cannot assert OLF parity here. Instead we check our own (self-consistent)
        # slope law Vm(reg) = vset - s_pu * Q_pu holds at the regulated bus B2.
        n = _svc_net(mode=VOLTAGE, target_v=402.0, slope=0.05, remote=True)
        model, V = self._run_ls(n, "G")
        reg_bus_view = 2  # B2 == VL2_0 is bus index 2
        s_pu = 0.05 * SN_MVA / 400.0
        q_pu = model.get_svcs()[0].res_q_mvar / SN_MVA
        self.assertLess(abs(abs(V[reg_bus_view]) - (402.0 / 400.0 - s_pu * q_pu)), self.tol,
                        "remote SVC slope law Vm(reg)=vset-s*Q broken")

    def test_svc_reactive_power(self):
        n = _svc_net(mode=REACTIVE_POWER, target_q=30.0)
        model = self._compare(n, "VL0", "G", "svc-Q")
        # OLF terminal q is load-sign; lightsim2grid res_q is injection -> opposite sign
        n.per_unit = False
        q_olf = n.get_static_var_compensators().loc["SVC1", "q"]
        self.assertAlmostEqual(model.get_svcs()[0].res_q_mvar, -q_olf, places=4)
        self.assertAlmostEqual(model.get_svcs()[0].res_q_mvar, -30.0, places=4)

    def test_svc_off_is_noop(self):
        n_off = _svc_net(mode=OFF)
        n_plain = _svc_net(with_svc=False)
        _, V_off = self._run_ls(n_off, "G")
        _, V_plain = self._run_ls(n_plain, "G")
        self.assertLess(np.abs(V_off - V_plain).max(), 1e-10)

    def test_mixed_gen_svc_same_bus_raises(self):
        # v1: an SVC may not share its regulated bus with another controller
        n = _star(g2_reg=None, g1_reg="LD", g1_tv=405.0)  # G1 regulates B2
        n.create_static_var_compensators(id="SVC1", voltage_level_id="VL2", bus_id="B2",
                                         regulation_mode="VOLTAGE", regulating=True,
                                         target_v=405.0, b_min=-1.0, b_max=1.0)
        with self.assertRaises(Exception):
            init_from_pypowsybl(n, gen_slack_id="G0", sort_index=True).ac_pf(
                np.ones(3, dtype=np.complex128), 30, 1e-10)

    # ----- backward compatibility & DC -------------------------------------
    def test_backward_compat_no_feature(self):
        # a plain grid (no remote control, no SVC) matches OLF and is structurally
        # unchanged: declaring the (empty) SVC container leaves J sparsity intact
        n = _svc_net(with_svc=False)
        self._compare(n, "VL0", "G", "backward-compat")

    def test_dc_remote_gen_angles(self):
        # DC: voltage control contributes nothing; angles still match OLF DC
        # (OLF leaves v_mag NaN in DC, so only angles are compared).
        n = _star(g2_reg=None, g1_reg="LD", g1_tv=405.0)
        ref = pp.loadflow.run_dc(n, parameters=self._params("VL0"))
        self.assertEqual(ref[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        _, V = self._run_ls(n, "G0", dc=True)
        self._assert_bus_parity(n, V, "dc-remote-gen", check_mag=False, check_angle=True)


if __name__ == "__main__":
    unittest.main()
