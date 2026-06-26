# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Validation of ``bake_outer_loops`` + ``compare_baked`` against PowSyBl Open Load
Flow.

Two groups of tests:

* ``test_olf_*`` -- pure pypowsybl: solve with outer loops, bake, solve
  loop-free, and check the baked loop-free solve reproduces the with-loops
  result. They exercise the bake on generators, VSC, SVC, shunts, and ratio +
  phase tap changers.
* ``test_ls_*`` -- check lightsim2grid reproduces the OLF loop-free result,
  including through line / transformer outages.
"""

import unittest

import numpy as np

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import (
        init_from_pypowsybl,
        bake_outer_loops,
        compare_baked,
        get_pypowsybl_loopfree_parameters,
    )
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

from global_var_tests import (
    CURRENT_PYPOW_VERSION,
    VERSION_PHASESHIFT_OK_PYPOW,
)


# Solver-tolerance thresholds.
TOL_VM_PU = 1e-3
TOL_VA_DEG = 5e-2
TOL_VM_KV = 1e-2  # for the pure-OLF kV-space check


# Default OLF provider parameters of pypowsybl 1.15.0 in the project's reference
# venv (``venv_ls``). They are pinned EXPLICITLY here, never read from the
# installed build's defaults: some pypowsybl builds ship the same version string
# but different defaults (e.g. ``extrapolateReactiveLimits``,
# ``stateVectorScalingMode``, ``maxRealisticVoltage``), which would otherwise make
# the with-loops solve -- and hence the "baked grid is inert" check -- behave
# differently from one environment to the next. Pinning every field makes these
# tests reproducible across builds.
_REF_PROVIDER_PARAMS = {
    'maxVoltageMismatch': '1.0E-4', 'generatorVoltageControlMinNominalVoltage': '-1.0',
    'startWithFrozenACEmulation': 'false', 'networkCacheEnabled': 'false',
    'reactiveRangeCheckMode': 'MAX', 'maxVoltageChangeStateVectorScalingMaxDphi': '0.17453292519943295',
    'maxRatioMismatch': '1.0E-5', 'areaInterchangeControl': 'false', 'loadPowerFactorConstant': 'false',
    'acSolverType': 'NEWTON_RAPHSON', 'actionableTransformersIds': '',
    'maxVoltageChangeStateVectorScalingMaxDv': '0.1', 'incrementalShuntControlOuterLoopMaxSectionShift': '3',
    'maxSusceptanceMismatch': '1.0E-4', 'maxNewtonKrylovIterations': '100',
    'reactivePowerDispatchMode': 'Q_EQUAL_PROPORTION', 'phaseShifterControlMode': 'CONTINUOUS_WITH_DISCRETISATION',
    'extrapolateReactiveLimits': 'false', 'asymmetrical': 'false', 'slackBusPMaxMismatch': '1.0',
    'maxActivePowerMismatch': '0.01', 'disableVoltageControlOfGeneratorsOutsideActivePowerLimits': 'false',
    'maxSlackBusCount': '1', 'mostMeshedSlackBusSelectorMaxNominalVoltagePercentile': '95.0',
    'maxNewtonRaphsonIterations': '15', 'minPlausibleTargetVoltage': '0.8', 'secondaryVoltageControl': 'false',
    'useActiveLimits': 'true', 'stateVectorScalingMode': 'NONE', 'useLoadModel': 'false',
    'voltageRemoteControlRobustMode': 'true', 'maxOuterLoopIterations': '20',
    'generatorReactivePowerRemoteControl': 'false', 'newtonRaphsonConvEpsPerEq': '1.0E-4',
    'lowImpedanceBranchMode': 'REPLACE_BY_ZERO_IMPEDANCE_LINE', 'actionableSwitchesIds': '',
    'maxRealisticVoltage': '2.0', 'fictitiousGeneratorVoltageControlCheckMode': 'FORCED',
    'maxReactivePowerMismatch': '0.01', 'transformerReactivePowerControl': 'false', 'minRealisticVoltage': '0.5',
    'fixVoltageTargets': 'false', 'acDcNetwork': 'false', 'simulateAutomationSystems': 'false',
    'forceTargetQInReactiveLimits': 'false', 'plausibleActivePowerLimit': '10000.0',
    'voltagePerReactivePowerControl': 'false', 'dcApproximationType': 'IGNORE_R', 'linePerUnitMode': 'IMPEDANCE',
    'referenceBusSelectionMode': 'FIRST_SLACK', 'newtonRaphsonStoppingCriteriaType': 'UNIFORM_CRITERIA',
    'disableInconsistentVoltageControls': 'false', 'generatorsWithZeroMwTargetAreNotStarted': 'true',
    'svcVoltageMonitoring': 'true', 'alwaysUpdateNetwork': 'false', 'slackBusSelectionMode': 'MOST_MESHED',
    'voltageInitModeOverride': 'NONE', 'slackBusCountryFilter': '', 'minNominalVoltageTargetVoltageCheck': '20.0',
    'areaInterchangePMaxMismatch': '2.0', 'maxPlausibleTargetVoltage': '1.2', 'lowImpedanceThreshold': '1.0E-8',
    'transformerVoltageControlUseInitialTapPosition': 'false', 'newtonKrylovLineSearch': 'false',
    'maxAngleMismatch': '1.0E-5', 'transformerVoltageControlMode': 'INCREMENTAL_VOLTAGE_CONTROL',
    'voltageRemoteControl': 'true', 'slackBusesIds': '', 'reportedFeatures': '',
    'areaInterchangeControlAreaType': 'ControlArea', 'shuntVoltageControlMode': 'WITH_GENERATOR_VOLTAGE_CONTROL',
    'lineSearchStateVectorScalingStepFold': '1.3333333333333333', 'reactiveLimitsMaxPqPvSwitch': '3',
    'incrementalTransformerRatioTapControlOuterLoopMaxTapShift': '3', 'lineSearchStateVectorScalingMaxIteration': '10',
    'slackDistributionFailureBehavior': 'FAIL', 'minNominalVoltageRealisticVoltageCheck': '0.0',
    'writeReferenceTerminals': 'true', 'voltageTargetPriorities': 'VOLTAGE_SOURCE_CONVERTER,GENERATOR,TRANSFORMER,SHUNT',
}


def _with_loops_params():
    # Every field is pinned to the ``venv_ls`` (pypowsybl 1.15.0) default so the
    # outer-loop solve is identical across pypowsybl builds (see
    # ``_REF_PROVIDER_PARAMS``). The single intentional deviation from that default
    # is ``twt_split_shunt_admittance=True`` (default is False): it matches the
    # transformer model the loop-free / lightsim2grid solve uses (see
    # get_pypowsybl_loopfree_parameters), which isolates the outer-loop effect that
    # baking neutralizes.
    return lf.Parameters(
        voltage_init_mode=lf.VoltageInitMode.UNIFORM_VALUES,
        transformer_voltage_control_on=False,
        use_reactive_limits=True,
        phase_shifter_regulation_on=False,
        twt_split_shunt_admittance=True,  # intentional deviation; see above
        shunt_compensator_voltage_control_on=False,
        read_slack_bus=True,
        write_slack_bus=True,
        distributed_slack=True,
        balance_type=lf.BalanceType.PROPORTIONAL_TO_GENERATION_P_MAX,
        dc_use_transformer_ratio=True,
        countries_to_balance=[],
        component_mode=lf.ComponentMode.MAIN_CONNECTED,
        dc_power_factor=1.0,
        hvdc_ac_emulation=True,
        dc=False,
        provider_parameters=dict(_REF_PROVIDER_PARAMS),
    )


def ieee14_with_qbind():
    """IEEE-14 with one generator's Q range tightened so the reactive-limit
    outer loop is forced to switch it PV->PQ. Exercises the freeze logic."""
    n = pp.network.create_ieee14()
    n.update_generators(id="B2-G", max_q=10.0, min_q=-10.0)
    return n


def ieee14_forced_pv_pq():
    """IEEE-14 where one generator is forced PV->PQ at its *upper* Q limit and
    another at its *lower* limit, by tightening each limit past the value that
    generator reaches with limits wide open. Exercises both freeze directions
    with the switch happening naturally during the loadflow."""
    n = pp.network.create_ieee14()
    n.update_generators(
        id=list(n.get_generators().index),
        min_q=[-9999] * 5, max_q=[9999] * 5,
    )
    # B3-G unconstrained Q_gen ~= 25.08 -> cap below to bind at max.
    n.update_generators(id="B3-G", max_q=20.08)
    # B6-G unconstrained Q_gen ~= 12.73 -> raise min above to bind at min.
    n.update_generators(id="B6-G", min_q=18.0)
    return n


def four_substations():
    """Node-breaker grid with VSC + LCC HVDC, an SVC regulating voltage, a
    shunt, and ratio + phase tap changers. Exercises every bake path."""
    return pp.network.create_four_substations_node_breaker_network()


def _olf_roundtrip_max_dev(network_factory):
    """Return (max |dV| kV, max |dAngle| deg) between OLF-with-loops and the
    baked OLF-loop-free solve. Pure pypowsybl; no lightsim2grid. This is the
    'without outer loop' test: the baked grid solved with every loop disabled
    must reproduce the original with-loops operating point."""
    with_loops = _with_loops_params()
    loop_free = get_pypowsybl_loopfree_parameters()

    n_ref = network_factory()
    lf.run_ac(n_ref, with_loops)
    ref = n_ref.get_buses()[["v_mag", "v_angle"]].copy()

    n_baked = network_factory()
    lf.run_ac(n_baked, with_loops)
    bake_outer_loops(n_baked)
    res = lf.run_ac(n_baked, loop_free)
    assert res[0].status == pp.loadflow.ComponentStatus.CONVERGED
    baked = n_baked.get_buses()[["v_mag", "v_angle"]]

    cmp = ref.join(baked, lsuffix="_r", rsuffix="_b")
    return (
        (cmp["v_mag_r"] - cmp["v_mag_b"]).abs().max(),
        (cmp["v_angle_r"] - cmp["v_angle_b"]).abs().max(),
    )


def _control_snapshot(n):
    """State of everything an AC outer loop can mutate. Used to detect, in a
    fully OLF-version-independent way, whether any outer loop took an action:
    if it did, at least one of these changes."""
    g = n.get_generators(attributes=["voltage_regulator_on", "q", "target_p"])
    snap = {
        "reg": g["voltage_regulator_on"].copy(),
        "q": g["q"].copy(),
        "target_p": g["target_p"].copy(),
    }
    for name, getter, col in [
        ("rtc", "get_ratio_tap_changers", "tap"),
        ("ptc", "get_phase_tap_changers", "tap"),
        ("shunt", "get_shunt_compensators", "section_count"),
        ("svc", "get_static_var_compensators", "regulation_mode"),
    ]:
        df = getattr(n, getter)(attributes=[col])
        snap[name] = df[col].copy() if len(df) else None
    return snap


def _state_max_change(a, b):
    """Largest discrete/continuous change between two control snapshots."""
    flips = int((a["reg"] != b["reg"]).sum())
    dq = float((a["q"] - b["q"]).abs().max())
    dtp = float((a["target_p"] - b["target_p"]).abs().max())
    discrete = flips
    for k in ("rtc", "ptc", "shunt", "svc"):
        if a[k] is not None and b[k] is not None:
            discrete += int((a[k] != b[k]).sum())
    return discrete, dq, dtp


@unittest.skipUnless(HAS_PYPOWSYBL, "pypowsybl is not installed")
@unittest.skipUnless(
    HAS_PYPOWSYBL and CURRENT_PYPOW_VERSION >= VERSION_PHASESHIFT_OK_PYPOW,
    "pypowsybl too old (no solved_tap_position / phase-shifter support)",
)
class TestOlfBake(unittest.TestCase):
    # -----------------------------------------------------------------
    # Pure-OLF tests (no lightsim2grid needed)
    # -----------------------------------------------------------------
    def _assert_baked_is_inert(self, network_factory, tol_kv=TOL_VM_KV, tol_q=1e-2):
        """OLF-version-independent proof that no outer loop acts on the baked grid.

        Two independent angles, neither parsing report text:

        (A) Solution agreement: on the baked grid, solving WITH all outer loops
            and solving loop-free give the same bus voltages. If any loop had
            acted in the with-loops run, that action would be absent from the
            loop-free run and the two would differ.

        (B) State invariance: re-running WITH loops on the baked grid flips no
            regulation flag, moves no tap / section / SVC mode, and leaves Q
            unchanged.
        """
        with_loops = _with_loops_params()
        loop_free = get_pypowsybl_loopfree_parameters()

        # (A) with-loops vs loop-free on the baked grid
        n_with = network_factory()
        lf.run_ac(n_with, with_loops)
        bake_outer_loops(n_with)
        res_with = lf.run_ac(n_with, with_loops)
        v_with = n_with.get_buses()[["v_mag", "v_angle"]].copy()

        n_free = network_factory()
        lf.run_ac(n_free, with_loops)
        bake_outer_loops(n_free)
        lf.run_ac(n_free, loop_free)
        v_free = n_free.get_buses()[["v_mag", "v_angle"]]

        cmp = v_with.join(v_free, lsuffix="_w", rsuffix="_f")
        self.assertLess(
            (cmp["v_mag_w"] - cmp["v_mag_f"]).abs().max(), tol_kv,
            "with-loops and loop-free disagree on baked grid -> a loop acted",
        )

        # (B) state invariance across the with-loops re-run
        n_state = network_factory()
        lf.run_ac(n_state, with_loops)
        bake_outer_loops(n_state)
        before = _control_snapshot(n_state)
        lf.run_ac(n_state, with_loops)
        after = _control_snapshot(n_state)
        discrete, dq, _dtp = _state_max_change(before, after)
        self.assertEqual(discrete, 0, f"{discrete} controllers changed state on baked grid")
        self.assertLess(dq, tol_q, f"Q moved by {dq} on baked re-run -> a loop acted")

        return res_with[0].status

    def test_olf_ieee14_reactive_limit_roundtrip(self):
        """WITHOUT outer loops: baked grid solved loop-free reproduces the
        original with-loops result."""
        dvm, dva = _olf_roundtrip_max_dev(ieee14_with_qbind)
        self.assertLess(dvm, TOL_VM_KV)
        self.assertLess(dva, 1e-2)

    def test_olf_four_substations_roundtrip(self):
        """WITHOUT outer loops, on the VSC+LCC HVDC / SVC / shunt / tap grid."""
        dvm, dva = _olf_roundtrip_max_dev(four_substations)
        self.assertLess(dvm, TOL_VM_KV)
        self.assertLess(dva, 1e-2)

    def test_olf_pv_pq_switch_both_directions(self):
        """The headline case: a generator that hits max_q and one that hits
        min_q are frozen to fixed-Q (PQ) at the binding value; generators inside
        their limits stay PV. The baked loop-free solve reproduces the
        with-limits one."""
        n_ref = ieee14_forced_pv_pq()
        lf.run_ac(n_ref, _with_loops_params())
        q_ref = n_ref.get_generators(attributes=["q"])["q"]

        n = ieee14_forced_pv_pq()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        g = n.get_generators(attributes=["voltage_regulator_on", "target_q"])

        # Bound generators frozen to PQ at the limit.
        self.assertFalse(g.loc["B3-G", "voltage_regulator_on"])
        self.assertLess(abs(g.loc["B3-G", "target_q"] - 20.08), 1e-1)
        self.assertFalse(g.loc["B6-G", "voltage_regulator_on"])
        self.assertLess(abs(g.loc["B6-G", "target_q"] - 18.0), 1e-1)
        # Unbound generators still PV.
        self.assertTrue(g.loc[["B1-G", "B2-G", "B8-G"], "voltage_regulator_on"].all())

        # Loop-free re-solve reproduces the with-limits reactive powers.
        res = lf.run_ac(n, get_pypowsybl_loopfree_parameters())
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        q_redo = n.get_generators(attributes=["q"])["q"]
        self.assertLess((q_redo - q_ref).abs().max(), 1e-2)

    def test_olf_ieee14_baked_inert(self):
        """WITH outer loops: they trigger nothing on the baked grid (robust
        check, no reporter)."""
        status = self._assert_baked_is_inert(ieee14_forced_pv_pq)
        self.assertEqual(status, pp.loadflow.ComponentStatus.CONVERGED)

    def test_olf_four_substations_baked_inert(self):
        """WITH outer loops on the HVDC/SVC/shunt/tap grid: nothing triggers."""
        status = self._assert_baked_is_inert(four_substations)
        self.assertEqual(status, pp.loadflow.ComponentStatus.CONVERGED)

    def test_olf_unbaked_loops_do_act_control(self):
        """Control: on the UN-baked grid the outer loops genuinely act, so the
        inertness checks above are not vacuously passing."""
        n_with = ieee14_forced_pv_pq()
        lf.run_ac(n_with, _with_loops_params())
        v_with = n_with.get_buses()[["v_mag"]].copy()
        n_free = ieee14_forced_pv_pq()
        lf.run_ac(n_free, get_pypowsybl_loopfree_parameters())
        v_free = n_free.get_buses()[["v_mag"]]
        cmp = v_with.join(v_free, lsuffix="_w", rsuffix="_f")
        self.assertGreater((cmp["v_mag_w"] - cmp["v_mag_f"]).abs().max(), 1e-2)

    # -----------------------------------------------------------------
    # Supplementary reporter-based check (OLF-version dependent: it parses
    # report text, whose wording can change between PowSyBl releases). The
    # robust checks above are the primary guarantees; this one is skipped
    # rather than failed if the marker strings ever stop matching.
    # -----------------------------------------------------------------
    _OUTER_LOOP_ACTION_MARKERS = (
        "Outer loop iteration",
        "PV -> PQ",
        "PQ -> PV",
        "switched",
    )

    def _outer_loop_acted(self, report_text):
        return any(m in report_text for m in self._OUTER_LOOP_ACTION_MARKERS)

    def test_olf_baked_reporter_shows_no_action(self):
        """Supplementary: a node reporter shows no outer-loop action on the
        baked grid. First asserts the markers DO fire on the un-baked grid; if
        they don't (e.g. PowSyBl reworded the report), the whole check is
        skipped, since the robust tests already cover the invariant."""
        n_orig = ieee14_forced_pv_pq()
        rep_orig = pp.report.ReportNode()
        lf.run_ac(n_orig, _with_loops_params(), report_node=rep_orig)
        if not self._outer_loop_acted(str(rep_orig)):
            self.skipTest("report markers not present in this PowSyBl version")

        n = ieee14_forced_pv_pq()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        rep = pp.report.ReportNode()
        lf.run_ac(n, _with_loops_params(), report_node=rep)
        self.assertFalse(
            self._outer_loop_acted(str(rep)),
            "reporter shows an outer loop acted on the baked grid:\n" + str(rep),
        )

    # -----------------------------------------------------------------
    # lightsim2grid agreement tests
    # -----------------------------------------------------------------
    def test_ls_no_outage_agrees(self):
        r = compare_baked(ieee14_with_qbind, slack_gen_id="B1-G")
        self.assertLess(r.max_dvm_pu, TOL_VM_PU)
        self.assertLess(r.max_dva_deg_offset_removed, TOL_VA_DEG)

    def test_ls_single_line_outage_agrees(self):
        r = compare_baked(
            ieee14_with_qbind, slack_gen_id="B1-G", line_outages=["L1-2-1"]
        )
        self.assertLess(r.max_dvm_pu, TOL_VM_PU)
        self.assertLess(r.max_dva_deg_offset_removed, TOL_VA_DEG)

    def test_ls_transformer_outage_agrees(self):
        r = compare_baked(
            ieee14_with_qbind, slack_gen_id="B1-G", trafo_outages=["T4-7-1"]
        )
        self.assertLess(r.max_dvm_pu, TOL_VM_PU)
        self.assertLess(r.max_dva_deg_offset_removed, TOL_VA_DEG)

    def test_ls_two_line_outage_agrees(self):
        """Disconnecting L1-2-1 and L7-9-1 strands the VL7/VL8 corner (VL8
        carries generator B8-G and connects out only through line L7-8 into the
        injection-free junction VL7). With the bus mapping taken from
        ``grid._ls_to_orig`` the two engines agree to solver tolerance here too,
        including on the junction bus."""
        r = compare_baked(
            ieee14_with_qbind,
            slack_gen_id="B1-G",
            line_outages=["L1-2-1", "L7-9-1"],
        )
        self.assertLess(r.max_dvm_pu, TOL_VM_PU)
        self.assertLess(r.max_dva_deg_offset_removed, TOL_VA_DEG)

    def test_ls_unbaked_disagrees(self):
        """Control: WITHOUT baking, OLF (full outer loops) and lightsim2grid
        must disagree after a line outage -- that disagreement is the whole
        reason baking exists."""
        LINE = "L1-2-1"
        n_olf = ieee14_with_qbind()
        n_olf.update_lines(id=LINE, connected1=False, connected2=False)
        lf.run_ac(n_olf, _with_loops_params())
        b = n_olf.get_buses().join(
            n_olf.get_voltage_levels()[["nominal_v"]], on="voltage_level_id"
        )
        olf_vm = np.sort((b["v_mag"] / b["nominal_v"]).to_numpy())

        n_ls = ieee14_with_qbind()
        grid = init_from_pypowsybl(
            n_ls, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False
        )
        grid.deactivate_powerline(list(n_ls.get_lines().index).index(LINE))
        V = grid.ac_pf(np.full(grid.total_bus(), 1.06 + 0j), 20, 1e-10)
        ls_vm = np.sort(np.abs(V))

        spread = np.max(np.abs(olf_vm - ls_vm))
        self.assertGreater(spread, 1e-2, f"expected disagreement, got {spread:.2e}")


if __name__ == "__main__":
    unittest.main()
