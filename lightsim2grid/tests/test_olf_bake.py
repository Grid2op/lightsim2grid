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
import pandas as pd

try:
    import pypowsybl as pp
    import pypowsybl.loadflow as lf
    from lightsim2grid.network import (
        init_from_pypowsybl,
        bake_outer_loops,
        compare_baked,
        remove_outer_loops,
    )
    HAS_PYPOWSYBL = True
except ImportError:
    HAS_PYPOWSYBL = False

from lightsim2grid.network.from_pypowsybl._olf_const import _Q_SATURATED_HELD_TOL_MVAR
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
    # ``_REF_PROVIDER_PARAMS``). The corresponding loop-free reference used
    # throughout this file is ``remove_outer_loops(_with_loops_params())``, not a
    # separately-defaulted factory: deriving it from this same object guarantees
    # every field it does not touch (twt_split_shunt_admittance, component_mode,
    # use_reactive_limits, voltage_init_mode, ...) is identical between the two
    # runs, isolating the outer-loop effect that baking neutralizes.
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


def ieee14_curve_reactive_range_too_small():
    """IEEE-14 with B2-G's reactive limits replaced by a degenerate CURVE (a
    fixed +/-0.3 MVAr range, well under OLF's 1 MVar plausibility floor)
    instead of the default MIN_MAX box. Mirrors real generators (small
    run-of-river hydro units, etc.) that OLF silently treats as PQ regardless
    of ``voltage_regulator_on``.

    The range is kept comfortably nonzero (0.6 MVAr) so the *ordinary*
    Q-at-limit saturation freeze (whose tolerance blows up and would
    coincidentally catch an exactly-zero range) does not also fire here --
    this exercises ``_bake_generator_voltage_control_discards`` alone. It is
    also the only fixture in this file with CURVE (rather than MIN_MAX)
    reactive limits, exercising ``_generator_max_reactive_range``'s
    curve-points code path.
    """
    n = pp.network.create_ieee14()
    n.create_curve_reactive_limits(
        id=["B2-G", "B2-G"], p=[0.0, 100.0], min_q=[-0.3, -0.3], max_q=[0.3, 0.3]
    )
    return n


def ieee14_curve_reactive_range_too_small_zero_target_q():
    """Same too-small (+/-0.3 MVAr) CURVE range as
    ``ieee14_curve_reactive_range_too_small``, but with B2-G's raw target_q
    pinned to 0 -- i.e. *inside* the tiny box -- instead of IEEE-14's default
    42.4 MVAr. Used only to isolate ``_bake_generator_voltage_control_discards``
    from the ordinary Q-at-limit saturation freeze: with a target_q outside the
    box (the other fixture), that pre-existing freeze also fires on its own
    (realized Q lands far outside the box either way) and would mask the flag
    having any effect; with target_q inside the box, only the new too-small-
    range check discards the generator, so disabling it is observable."""
    n = ieee14_curve_reactive_range_too_small()
    n.update_generators(id="B2-G", target_q=0.0)
    return n


def ieee14_implausible_target_v():
    """IEEE-14 with B2-G's target_v set far outside OLF's plausible target-
    voltage window (0.8-1.2 pu of nominal): 50 kV on a 135 kV bus is ~0.37 pu."""
    n = pp.network.create_ieee14()
    n.update_generators(id="B2-G", target_v=50.0)
    return n


def star_shared_remote_control_one_switched():
    """Two generators on DIFFERENT buses, G1@B1 and G2@B3, remotely regulating the
    same load bus B2 at 405 kV. G1's upper reactive limit (5 MVAr) is well under its
    share of the group's reactive output, so OLF's reactive-limit loop switches it
    to PQ at 5 MVAr and G2 alone keeps holding B2 -- whose voltage therefore still
    reads as on target, for G1 as much as for G2. How much each one injects sets
    the voltage of its own bus (B1, B3), so a bake that keeps G1 regulating hands
    the loop-free solve a different split and moves B3 by ~2 kV."""
    from test_voltage_control_pypowsybl import _star
    return _star(extra_q_g1=(-100.0, 5.0), extra_q_g2=(-100.0, 200.0),
                 g1_reg="LD", g2_reg="LD", g1_tv=405.0, g2_tv=405.0)


def _olf_roundtrip_max_dev(network_factory):
    """Return (max |dV| kV, max |dAngle| deg) between OLF-with-loops and the
    baked OLF-loop-free solve. Pure pypowsybl; no lightsim2grid. This is the
    'without outer loop' test: the baked grid solved with every loop disabled
    must reproduce the original with-loops operating point."""
    with_loops = _with_loops_params()
    loop_free = remove_outer_loops(with_loops)

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
        loop_free = remove_outer_loops(with_loops)

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
        res = lf.run_ac(n, remove_outer_loops(_with_loops_params()))
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        q_redo = n.get_generators(attributes=["q"])["q"]
        self.assertLess((q_redo - q_ref).abs().max(), 1e-2)

    def test_bake_returns_the_pinned_generators_and_init_flags_them(self):
        """`bake_outer_loops` returns the generators it froze AT A REACTIVE LIMIT (and
        nothing on a second, idempotent bake); `init_from_pypowsybl(can_be_pv=...)` flags
        exactly those, and the result network shows the flag."""
        from lightsim2grid.network.from_pypowsybl import LightsimResultNetwork
        n = ieee14_forced_pv_pq()
        lf.run_ac(n, _with_loops_params())
        pinned = bake_outer_loops(n)
        self.assertEqual(set(pinned), {"B3-G", "B6-G"})
        self.assertEqual(len(bake_outer_loops(n)), 0)  # already baked: nothing left at a limit

        grid = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                   can_be_pv=pinned)
        flags = {g.name: g.can_be_pv for g in grid.get_generators()}
        self.assertEqual(flags, {"B1-G": False, "B2-G": False, "B3-G": True, "B6-G": True,
                                 "B8-G": False})
        grid.ac_pf(np.full(grid.total_bus(), 1.06 + 0j), 20, 1e-10)
        res = LightsimResultNetwork(grid, n).get_generators()
        self.assertEqual(res.loc[["B3-G", "B6-G"], "can_be_pv"].tolist(), [True, True])
        self.assertFalse(res.loc[["B1-G", "B2-G", "B8-G"], "can_be_pv"].any())
        # the other ways of saying it
        grid_b = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                     can_be_pv=pd.Series(True, index=["B6-G"]))
        self.assertEqual([g.can_be_pv for g in grid_b.get_generators()],
                         [g.name == "B6-G" for g in grid_b.get_generators()])
        # nothing flagged by default, and an unknown id is refused
        grid_c = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False)
        self.assertFalse(any(g.can_be_pv for g in grid_c.get_generators()))
        with self.assertRaises(ValueError):
            init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                can_be_pv=["B3-G", "NOT-A-GEN"])

    def test_pinned_generator_released_is_reported(self):
        """The baked grid pins B3-G as PQ at its max_q. Drop the load on its bus: the
        voltage rises above the target it would hold, which OLF's loop (on the raw grid)
        answers by keeping the unit PV inside its range -- and which lightsim2grid, unable
        to release it, reports as HIGH_VOLTAGE_AT_MAX_Q on that generator. B6-G, pinned at
        its min_q with a voltage ABOVE its target, is released by neither."""
        from lightsim2grid.lightsim2grid_cpp import LimitViolationType, ViolationElementType
        n = ieee14_forced_pv_pq()
        lf.run_ac(n, _with_loops_params())
        pinned = bake_outer_loops(n)

        # untouched: the baked grid is the reference solve, nothing to release
        grid = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                   can_be_pv=pinned)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.06 + 0j), 20, 1e-10)
        self.assertGreater(V.shape[0], 0)
        release = [v for v in grid.get_physical_violations()
                   if v.violation_type in (LimitViolationType.LOW_VOLTAGE_AT_MIN_Q,
                                           LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)]
        self.assertEqual(release, [])

        # the load on bus 3 gone
        n.update_loads(id="B3-L", p0=0., q0=0.)
        grid = init_from_pypowsybl(n, gen_slack_id="B1-G", sort_index=False, buses_for_sub=False,
                                   can_be_pv=pinned)
        V = grid.ac_pf(np.full(grid.total_bus(), 1.06 + 0j), 20, 1e-10)
        self.assertGreater(V.shape[0], 0)
        release = [v for v in grid.get_physical_violations()
                   if v.violation_type in (LimitViolationType.LOW_VOLTAGE_AT_MIN_Q,
                                           LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)]
        self.assertEqual([v.name for v in release], ["B3-G"])
        viol = release[0]
        self.assertEqual(viol.element_type, ViolationElementType.GENERATOR)
        self.assertEqual(viol.violation_type, LimitViolationType.HIGH_VOLTAGE_AT_MAX_Q)
        target_v = n.get_generators(attributes=["target_v"]).loc["B3-G", "target_v"]
        self.assertAlmostEqual(viol.limit, target_v, places=6)
        self.assertGreater(viol.value, viol.limit + 1.)  # kV, well above the target

        # the reference, outer loops on, on the raw perturbed grid: B3-G regulates again,
        # strictly inside its range, and its bus sits on target; B6-G stays at its limit
        n_ref = ieee14_forced_pv_pq()
        n_ref.update_loads(id="B3-L", p0=0., q0=0.)
        lf.run_ac(n_ref, _with_loops_params())
        g = n_ref.get_generators(attributes=["bus_id", "q", "target_v"])
        v = n_ref.get_buses(attributes=["v_mag"])["v_mag"]
        q_gen_b3 = -g.loc["B3-G", "q"]  # generator convention
        self.assertLess(q_gen_b3, 20.08 - 0.5)
        self.assertAlmostEqual(v[g.loc["B3-G", "bus_id"]], g.loc["B3-G", "target_v"], places=3)
        self.assertAlmostEqual(-g.loc["B6-G", "q"], 18.0, places=3)

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
        lf.run_ac(n_free, remove_outer_loops(_with_loops_params()))
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
    # Voltage-control discards beyond Q-limit saturation: too-small
    # reactive range and implausible target_v
    # (_bake_generator_voltage_control_discards).
    # -----------------------------------------------------------------
    def test_olf_shared_remote_control_switched_member_frozen(self):
        """A member of a shared remote voltage control switched to PQ at its limit
        is frozen even though the group's regulated bus is still held (by the other
        member), and the loop-free solve then reproduces the with-loops voltages at
        the controller buses too. Regression: the "target held" test is read at the
        regulated bus, so it used to keep every member of such a group regulating."""
        n = star_shared_remote_control_one_switched()
        lf.run_ac(n, _with_loops_params())
        gen = n.get_generators(attributes=["q", "voltage_regulator_on"])
        self.assertAlmostEqual(-gen.at["G1", "q"], 5.0, places=6, msg="fixture: G1 is not switched at its limit")
        bake_outer_loops(n)
        gen = n.get_generators(attributes=["target_q", "voltage_regulator_on"])
        self.assertFalse(gen.at["G1", "voltage_regulator_on"])
        self.assertAlmostEqual(gen.at["G1", "target_q"], 5.0, places=6)
        self.assertTrue(gen.at["G2", "voltage_regulator_on"])

        dvm, dva = _olf_roundtrip_max_dev(star_shared_remote_control_one_switched)
        self.assertLess(dvm, TOL_VM_KV)
        self.assertLess(dva, 1e-2)

    def test_olf_saturated_held_unit_frozen_on_request(self):
        """A unit whose target OLF held while its Q sits exactly at a limit stays PV
        by default; bake_saturated_voltage_control=True freezes it at that limit. The
        loop-free solve reproduces the with-loops voltages either way."""
        loop_free = remove_outer_loops(_with_loops_params())
        for flag in (False, True):
            n = pp.network.create_ieee14()
            n.update_generators(id=list(n.get_generators().index), min_q=[-9999] * 5, max_q=[9999] * 5)
            lf.run_ac(n, _with_loops_params())
            ref = n.get_buses()[["v_mag", "v_angle"]].copy()
            # the upper limit put where B3-G already sits: target held, Q at max_q
            q_b3 = -n.get_generators(attributes=["q"]).at["B3-G", "q"]
            n.update_generators(id="B3-G", max_q=q_b3)
            bake_outer_loops(n, bake_saturated_voltage_control=flag)
            gen = n.get_generators(attributes=["target_q", "voltage_regulator_on"])
            self.assertEqual(gen.at["B3-G", "voltage_regulator_on"], not flag)
            if flag:
                self.assertAlmostEqual(gen.at["B3-G", "target_q"], q_b3, places=9)
            # no other unit is at a limit: nothing else frozen
            others = gen.drop(index="B3-G")
            self.assertTrue(others["voltage_regulator_on"].all())

            res = lf.run_ac(n, loop_free)
            self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
            cmp = ref.join(n.get_buses()[["v_mag", "v_angle"]], lsuffix="_r", rsuffix="_b")
            self.assertLess((cmp["v_mag_r"] - cmp["v_mag_b"]).abs().max(), TOL_VM_KV)
            self.assertLess((cmp["v_angle_r"] - cmp["v_angle_b"]).abs().max(), 1e-2)

    def test_olf_held_unit_with_headroom_not_frozen(self):
        """bake_saturated_voltage_control only freezes a held unit at its limit to
        _Q_SATURATED_HELD_TOL_MVAR: one still inside the relative tolerance of its limit but
        with more headroom than that keeps regulating, as OLF does."""
        for headroom_mvar, frozen in [(0.5 * _Q_SATURATED_HELD_TOL_MVAR, True),
                                      (5. * _Q_SATURATED_HELD_TOL_MVAR, False)]:
            n = pp.network.create_ieee14()
            n.update_generators(id=list(n.get_generators().index), min_q=[-9999] * 5, max_q=[9999] * 5)
            lf.run_ac(n, _with_loops_params())
            q_b3 = -n.get_generators(attributes=["q"]).at["B3-G", "q"]
            # a Q range wide enough for the relative tolerance to cover the headroom
            n.update_generators(id="B3-G", min_q=q_b3 - 1000., max_q=q_b3 + headroom_mvar)
            bake_outer_loops(n, bake_saturated_voltage_control=True)
            reg = n.get_generators(attributes=["voltage_regulator_on"]).at["B3-G", "voltage_regulator_on"]
            self.assertEqual(reg, not frozen, f"headroom {headroom_mvar} MVAr")

    def test_hit_qlimit_picks_the_nearer_limit(self):
        from lightsim2grid.network.from_pypowsybl._olf_bake import _hit_qlimit
        import pandas as pd
        df = pd.DataFrame({"min_q": [-10.0, -10.0, np.nan, -4.0], "max_q": [5.0, 5.0, 3.0, np.nan]},
                          index=["at_max", "at_min", "no_min", "no_max"])
        q_gen = pd.Series([4.9995, -9.998, 2.999, -3.99], index=df.index)
        hit = _hit_qlimit(df, q_gen)
        self.assertEqual(hit.to_dict(), {"at_max": 5.0, "at_min": -10.0, "no_min": 3.0, "no_max": -4.0})

    def test_olf_switched_unit_baked_at_its_limit_not_reported_q(self):
        """A unit switched at its Q limit injects that limit, but OLF does not always
        report it (on real grid snapshots it re-splits a bus' reactive target among the
        units of the bus when writing results). The bake must write the limit: here the
        reported q of one of two switched units is overwritten slightly inside its limit
        after the solve, standing for such a report."""
        n = pp.network.create_ieee14()
        g0 = n.get_generators(attributes=["bus_breaker_bus_id", "voltage_level_id", "target_v"]).loc["B2-G"]
        n.create_generators(id="B2-G2", voltage_level_id=g0["voltage_level_id"], bus_id=g0["bus_breaker_bus_id"],
                            target_p=10.0, target_q=0.0, target_v=g0["target_v"], voltage_regulator_on=True,
                            max_p=100.0, min_p=0.0)
        n.create_minmax_reactive_limits(id=["B2-G", "B2-G2"], min_q=[-10.0, -10.0], max_q=[5.0, 3.0])
        lf.run_ac(n, _with_loops_params())
        q = n.get_generators(attributes=["q"])["q"]
        self.assertAlmostEqual(-q["B2-G"], 5.0, places=6, msg="fixture: B2-G not switched at its limit")
        self.assertAlmostEqual(-q["B2-G2"], 3.0, places=6, msg="fixture: B2-G2 not switched at its limit")
        # B2-G2 reported 5e-4 MVAr inside the limit it injects (a misreport: baked at the
        # limit); B2-G reported 0.05 MVAr inside its limit, within the relative at-limit
        # tolerance (0.075) but too far for a misreport: a unit that settled there
        # injects what it reports, baked as reported
        n.update_generators(id=["B2-G", "B2-G2"], q=[-4.95, -2.9995])
        bake_outer_loops(n)
        gen = n.get_generators(attributes=["target_q", "voltage_regulator_on"])
        self.assertFalse(gen.at["B2-G", "voltage_regulator_on"])
        self.assertFalse(gen.at["B2-G2", "voltage_regulator_on"])
        self.assertEqual(gen.at["B2-G2", "target_q"], 3.0)
        self.assertAlmostEqual(gen.at["B2-G", "target_q"], 4.95, places=9)

    @staticmethod
    def _ieee14_below_curve():
        """IEEE-14 with B3-G (dispatched at 0 MW) on a reactive capability curve that
        only covers 50-100 MW (max_q 20 -> 21), every other limit wide open. OLF's
        ``extrapolateReactiveLimits`` extends the curve's first segment down to 0 MW:
        max_q 19, where pypowsybl's max_q_at_p clamps to 20."""
        n = pp.network.create_ieee14()
        n.update_generators(id=list(n.get_generators().index), min_q=[-9999] * 5, max_q=[9999] * 5)
        n.create_curve_reactive_limits(id=["B3-G", "B3-G"], p=[50.0, 100.0],
                                       min_q=[-9999.0, -9999.0], max_q=[20.0, 21.0])
        return n

    @staticmethod
    def _extrapolating_params():
        params = _with_loops_params()
        prov = dict(params.provider_parameters)
        prov["extrapolateReactiveLimits"] = "true"
        params.provider_parameters = prov
        return params

    def test_extrapolate_curve_limits(self):
        from lightsim2grid.network.from_pypowsybl._olf_bake import _extrapolate_curve_limits
        n = self._ieee14_below_curve()
        n.create_curve_reactive_limits(id=["B2-G", "B2-G"], p=[0.0, 100.0], min_q=[-30.0, -40.0], max_q=[30.0, 50.0])
        lf.run_ac(n, self._extrapolating_params())
        gen = n.get_generators(attributes=["p", "min_q_at_p", "max_q_at_p"])
        self.assertAlmostEqual(gen.at["B3-G", "max_q_at_p"], 20.0, places=9, msg="pypowsybl no longer clamps")
        out = _extrapolate_curve_limits(n, gen)
        self.assertAlmostEqual(out.at["B3-G", "max_q_at_p"], 19.0, places=9)
        self.assertAlmostEqual(out.at["B3-G", "min_q_at_p"], -9999.0, places=9)
        # inside its curve: untouched; without a curve: untouched
        self.assertEqual(out.at["B2-G", "max_q_at_p"], gen.at["B2-G", "max_q_at_p"])
        self.assertEqual(out.at["B6-G", "max_q_at_p"], gen.at["B6-G", "max_q_at_p"])

    def test_olf_switched_below_its_curve_baked_at_extrapolated_limit(self):
        """OLF switches B3-G at its extrapolated limit (19 MVAr); its reported q, nudged
        5e-4 MVAr away, must be baked at 19 -- not kept as reported because the clamped
        limit (20) is a whole MVAr away -- and the loop-free solve must reproduce the
        voltages. With ``extrapolate_reactive_limits=False`` the reported value is kept."""
        with_loops = self._extrapolating_params()
        n_ref = self._ieee14_below_curve()
        lf.run_ac(n_ref, with_loops)
        self.assertAlmostEqual(-n_ref.get_generators().at["B3-G", "q"], 19.0, places=6,
                               msg="fixture: OLF no longer switches B3-G at the extrapolated limit")
        ref = n_ref.get_buses()[["v_mag", "v_angle"]].copy()

        for extrapolate, expected in ((True, 19.0), (False, 18.9995)):
            n = self._ieee14_below_curve()
            lf.run_ac(n, with_loops)
            n.update_generators(id="B3-G", q=-18.9995)
            bake_outer_loops(n, extrapolate_reactive_limits=extrapolate)
            gen = n.get_generators(attributes=["target_q", "voltage_regulator_on"])
            self.assertFalse(gen.at["B3-G", "voltage_regulator_on"], f"extrapolate={extrapolate}")
            self.assertAlmostEqual(gen.at["B3-G", "target_q"], expected, places=9, msg=f"extrapolate={extrapolate}")
            if extrapolate:
                res = lf.run_ac(n, remove_outer_loops(with_loops))
                self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
                cmp = ref.join(n.get_buses()[["v_mag", "v_angle"]], lsuffix="_r", rsuffix="_b")
                self.assertLess((cmp["v_mag_r"] - cmp["v_mag_b"]).abs().max(), TOL_VM_KV)
                self.assertLess((cmp["v_angle_r"] - cmp["v_angle_b"]).abs().max(), 1e-2)

    def test_olf_pq_target_q_forced_in_limits_baked(self):
        """A non-regulating generator whose target_q lies outside its reactive limits
        injects the limit under OLF's ``forceTargetQInReactiveLimits`` (the default of
        recent OLF, pinned off in _REF_PROVIDER_PARAMS, hence turned on here). The
        loop-free solve, reactive limits off, would inject the raw target: the bake must
        write the clamped value into target_q."""
        def factory():
            n = pp.network.create_ieee14()
            n.update_generators(id="B8-G", voltage_regulator_on=False, target_q=0.0)
            n.create_minmax_reactive_limits(id="B8-G", min_q=5.0, max_q=20.0)
            return n

        with_loops = _with_loops_params()
        prov = dict(with_loops.provider_parameters)
        prov["forceTargetQInReactiveLimits"] = "true"
        with_loops.provider_parameters = prov
        loop_free = remove_outer_loops(with_loops)

        n_ref = factory()
        lf.run_ac(n_ref, with_loops)
        self.assertAlmostEqual(-n_ref.get_generators().at["B8-G", "q"], 5.0, places=6,
                               msg="fixture: OLF does not clamp target_q into the limits")
        ref = n_ref.get_buses()[["v_mag", "v_angle"]].copy()

        n = factory()
        lf.run_ac(n, with_loops)
        bake_outer_loops(n)
        self.assertAlmostEqual(n.get_generators().at["B8-G", "target_q"], 5.0, places=6)
        res = lf.run_ac(n, loop_free)
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        cmp = ref.join(n.get_buses()[["v_mag", "v_angle"]], lsuffix="_r", rsuffix="_b")
        self.assertLess((cmp["v_mag_r"] - cmp["v_mag_b"]).abs().max(), TOL_VM_KV)
        self.assertLess((cmp["v_angle_r"] - cmp["v_angle_b"]).abs().max(), 1e-2)

    def test_olf_reactive_range_too_small_sharing_a_held_bus_frozen(self):
        """A unit with too small a reactive range on the bus of a controller OLF keeps:
        the bus is held (by that other unit), but OLF discarded the small one all the
        same. Bake must freeze it at its realized q -- left regulating, it would add its
        range to the bus' capability -- and leave the other one regulating; the baked
        loop-free re-solve reproduces the reference."""
        def make():
            n = pp.network.create_ieee14()
            n.create_generators(id="B2-SMALL", voltage_level_id="VL2", bus_id="B2",
                                target_p=5.0, target_q=0.1, target_v=141.075,
                                voltage_regulator_on=True, min_p=0.0, max_p=20.0)
            n.update_generators(id="B2-SMALL", min_q=-0.3, max_q=0.3)
            return n
        n_ref = make()
        lf.run_ac(n_ref, _with_loops_params())
        q_ref = n_ref.get_generators(attributes=["q"])["q"]
        # OLF did not regulate it: its raw target_q, not a share of the bus' reactive power
        self.assertAlmostEqual(-q_ref["B2-SMALL"], 0.1, places=6)

        n = make()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        g = n.get_generators(attributes=["voltage_regulator_on", "target_q"])
        self.assertFalse(g.loc["B2-SMALL", "voltage_regulator_on"])
        self.assertAlmostEqual(g.loc["B2-SMALL", "target_q"], 0.1, places=6)
        self.assertTrue(g.loc["B2-G", "voltage_regulator_on"])

        res = lf.run_ac(n, remove_outer_loops(_with_loops_params()))
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        q_redo = n.get_generators(attributes=["q"])["q"]
        self.assertLess((q_redo - q_ref).abs().max(), 1e-2)

    def test_olf_reactive_range_too_small_frozen(self):
        """A CURVE-kind generator with a sub-1-MVar reactive range is not
        actually voltage-controlled by OLF: its realized Q sits far outside
        the tiny +/-0.3 MVAr box the curve declares, proving OLF fell back to
        the generator's raw (unconfined) target_q rather than confining it
        through voltage control. Bake must freeze it to fixed-Q at that
        realized q, and the baked loop-free re-solve must reproduce it."""
        n_ref = ieee14_curve_reactive_range_too_small()
        lf.run_ac(n_ref, _with_loops_params())
        # not actually voltage-controlled: q falls way outside the +/-0.3 box
        self.assertGreater(
            abs(n_ref.get_generators(attributes=["q"]).loc["B2-G", "q"]), 1.0
        )

        n = ieee14_curve_reactive_range_too_small()
        lf.run_ac(n, _with_loops_params())
        q_ref = n.get_generators(attributes=["q"]).loc["B2-G", "q"]
        bake_outer_loops(n)
        g = n.get_generators(attributes=["voltage_regulator_on", "target_q"])
        self.assertFalse(g.loc["B2-G", "voltage_regulator_on"])
        self.assertLess(abs(g.loc["B2-G", "target_q"] - (-q_ref)), 1e-2)
        # other generators (plain MIN_MAX, ample range) stay untouched
        self.assertTrue(g.loc[["B1-G", "B3-G", "B6-G", "B8-G"], "voltage_regulator_on"].all())

        res = lf.run_ac(n, remove_outer_loops(_with_loops_params()))
        self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
        q_redo = n.get_generators(attributes=["q"]).loc["B2-G", "q"]
        self.assertLess(abs(q_redo - q_ref), 1e-2)

    def test_olf_reactive_range_too_small_flag_off(self):
        """``bake_generator_voltage_control_discards=False`` leaves the
        too-small-range generator regulating voltage. Uses the target_q=0
        variant so the pre-existing Q-at-limit saturation freeze -- which
        would otherwise also catch this generator on its own, masking the
        flag -- does not fire (see the fixture's docstring)."""
        n = ieee14_curve_reactive_range_too_small_zero_target_q()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n, bake_generator_voltage_control_discards=False)
        self.assertTrue(
            n.get_generators(attributes=["voltage_regulator_on"]).loc["B2-G", "voltage_regulator_on"]
        )

    def test_olf_reactive_range_too_small_zero_target_q_frozen(self):
        """Same as ``test_olf_reactive_range_too_small_frozen`` but with the
        flag on (default): confirms the new check alone -- independent of
        the ordinary saturation freeze -- discards this generator too."""
        n = ieee14_curve_reactive_range_too_small_zero_target_q()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        self.assertFalse(
            n.get_generators(attributes=["voltage_regulator_on"]).loc["B2-G", "voltage_regulator_on"]
        )

    def test_olf_implausible_target_v_frozen(self):
        """A generator with target_v far outside OLF's plausible window is
        frozen to fixed-Q, at the realized q, the same way."""
        n_ref = ieee14_implausible_target_v()
        lf.run_ac(n_ref, _with_loops_params())
        q_ref = n_ref.get_generators(attributes=["q"]).loc["B2-G", "q"]
        b_ref = n_ref.get_buses(attributes=["v_mag"])
        bus_id = n_ref.get_generators(attributes=["bus_id"]).loc["B2-G", "bus_id"]
        self.assertGreater(abs(b_ref.loc[bus_id, "v_mag"] - 50.0), 1.0)

        n = ieee14_implausible_target_v()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        g = n.get_generators(attributes=["voltage_regulator_on", "target_q"])
        self.assertFalse(g.loc["B2-G", "voltage_regulator_on"])
        self.assertLess(abs(g.loc["B2-G", "target_q"] - (-q_ref)), 1e-2)
        self.assertTrue(g.loc[["B1-G", "B3-G", "B6-G", "B8-G"], "voltage_regulator_on"].all())

    def test_olf_implausible_target_v_flag_off(self):
        n = ieee14_implausible_target_v()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n, bake_generator_voltage_control_discards=False)
        self.assertTrue(
            n.get_generators(attributes=["voltage_regulator_on"]).loc["B2-G", "voltage_regulator_on"]
        )

    # -----------------------------------------------------------------
    # Active-power (slack-distribution) participation zeroing
    # (_bake_active_power_control_participation).
    # -----------------------------------------------------------------
    def test_olf_active_power_control_participation_excluded(self):
        """B3-G/B6-G/B8-G sit at target_p=0 MW by default on IEEE-14, with
        min_p < 0 (so they are NOT also frozen by the "not started" voltage
        rule -- this isolates the active-power-only exclusion); B1-G gets an
        implausible max_p. All four get participate=False written into the
        network's own activePowerControl extension; the untouched B2-G does
        not get an extension entry at all."""
        n = pp.network.create_ieee14()
        n.update_generators(id="B1-G", max_p=20000.0)
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n)
        apc = n.get_extensions("activePowerControl")
        for gid in ["B1-G", "B3-G", "B6-G", "B8-G"]:
            self.assertIn(gid, apc.index, f"{gid} should have an activePowerControl entry")
            self.assertFalse(bool(apc.loc[gid, "participate"]), f"{gid} should not participate")
        self.assertNotIn("B2-G", apc.index)

    def test_olf_active_power_round_off_snapped_on_bound(self):
        """OLF writes a unit's p back with a round-off, so a unit dispatched at its max_p
        can come back a hair above it (seen on real grid snapshots): written as is into
        target_p, OLF's checkActivePowerControl then takes the unit out of the slack of
        every solve of the baked network. The bake puts it back on the bound; a realized
        p further off than the round-off tolerance is written as is."""
        for p_over, expected in [(1e-11, 40.), (1., 41.)]:
            n = pp.network.create_ieee14()
            n.update_generators(id="B2-G", max_p=40.)  # dispatched at its max_p
            lf.run_ac(n, _with_loops_params())
            n.update_generators(id="B2-G", p=-(40. + p_over))
            bake_outer_loops(n)
            self.assertEqual(n.get_generators().loc["B2-G", "target_p"], expected)

    def test_snap_realized_into_target_range(self):
        from lightsim2grid.network.from_pypowsybl._olf_bake import _snap_realized_into_target_range
        min_tp = np.array([10., 10., 10., 10., 10., 10.])
        max_tp = np.array([40., 40., 40., 40., 40., 40.])
        target_p = np.array([40., 10., 45., 25., 40., 40.])  # 3rd: outside before the bake
        realized = np.array([40. + 1e-13, 10. - 1e-13, 45. + 1e-13, 25., 41., np.nan])
        res = _snap_realized_into_target_range(realized, target_p, min_tp, max_tp)
        np.testing.assert_array_equal(res, [40., 10., 45. + 1e-13, 25., 41., np.nan])

    def test_olf_active_power_control_participation_flag_off(self):
        """``bake_active_power_control_participation=False`` creates no
        activePowerControl extension at all."""
        n = pp.network.create_ieee14()
        lf.run_ac(n, _with_loops_params())
        bake_outer_loops(n, bake_active_power_control_participation=False)
        apc = n.get_extensions("activePowerControl")
        self.assertEqual(len(apc), 0)

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

    def test_copy_parameters_keeps_component_mode(self):
        """The deprecated ``connected_component_mode`` alias reads None for
        MAIN_SYNCHRONOUS and writes None back as ALL_CONNECTED; copying it after
        ``component_mode`` (hash-seed dependent) used to widen the solve to every
        island, and a warm-started solve then read an unsolved island's NaN V."""
        from lightsim2grid.network.from_pypowsybl._olf_params import _copy_parameters
        modes = [lf.Parameters().component_mode]  # the version's own default
        modes += [lf.ComponentMode.MAIN_CONNECTED, lf.ComponentMode.ALL_CONNECTED]
        if hasattr(lf.ComponentMode, "MAIN_SYNCHRONOUS"):
            modes.append(lf.ComponentMode.MAIN_SYNCHRONOUS)
        for mode in modes:
            params = lf.Parameters(component_mode=mode)
            self.assertEqual(_copy_parameters(params).component_mode, mode)
            self.assertEqual(remove_outer_loops(params).component_mode, mode)


    def test_svc_near_its_limit_follows_the_generator_rule(self):
        """An SVC within the saturation tolerance of its limit, whose regulated bus the
        reference solve held at its target, was regulating: left in VOLTAGE mode (frozen
        only with ``bake_saturated_voltage_control``). One that could not hold its target
        is frozen, as a generator."""
        def _baked_mode(b_max_mvar, **kwargs):
            n = four_substations()
            # the SVC holds its bus at 400 kV producing a bit more than 12.5 MVAr
            n.update_static_var_compensators(id="SVC", b_max=b_max_mvar / 400. ** 2)
            lf.run_ac(n, _with_loops_params())
            bake_outer_loops(n, **kwargs)
            return n.get_static_var_compensators().loc["SVC", "regulation_mode"]

        # a hair inside its limit, far within the relative tolerance, target held
        self.assertEqual(_baked_mode(12.9), "VOLTAGE")
        # ... but with more headroom than _Q_SATURATED_HELD_TOL_MVAR: still regulating
        self.assertEqual(_baked_mode(12.9, bake_saturated_voltage_control=True), "VOLTAGE")
        # saturated for real: the target is not held
        self.assertEqual(_baked_mode(10.), "REACTIVE_POWER")

    def test_saturated_held_svc_frozen_on_request(self):
        """With ``bake_saturated_voltage_control``, a held SVC is frozen only when its Q sits
        at its limit to ``_Q_SATURATED_HELD_TOL_MVAR``."""
        n = four_substations()
        lf.run_ac(n, _with_loops_params())
        q_svc = -n.get_static_var_compensators().loc["SVC", "q"]  # generator convention
        self.assertGreater(q_svc, 1.)
        for headroom_mvar, mode in [(0.5 * _Q_SATURATED_HELD_TOL_MVAR, "REACTIVE_POWER"),
                                    (5. * _Q_SATURATED_HELD_TOL_MVAR, "VOLTAGE")]:
            n = four_substations()
            lf.run_ac(n, _with_loops_params())
            v_kv = n.get_buses().loc[n.get_static_var_compensators().loc["SVC", "bus_id"], "v_mag"]
            n.update_static_var_compensators(id="SVC", b_max=(q_svc + headroom_mvar) / v_kv ** 2)
            bake_outer_loops(n, bake_saturated_voltage_control=True)
            self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], mode,
                             f"headroom {headroom_mvar} MVAr")

    def test_unit_clamped_inside_a_regulating_bus_not_frozen(self):
        """GB2 shares bus B with GB, both holding the remote load bus: OLF splits the bus'
        reactive power among them and clamps GB2 at its small limit, the bus (GB2 too)
        still regulating. Neither freezing rule freezes it."""
        def build():
            n = pp.network.create_empty()
            n.create_substations(id=["S"])
            n.create_voltage_levels(id=["VR", "VA", "VB", "VS"], substation_id=["S"] * 4,
                                    topology_kind=["BUS_BREAKER"] * 4, nominal_v=[225., 20., 20., 225.])
            n.create_buses(id=["R", "A", "B", "SB"], voltage_level_id=["VR", "VA", "VB", "VS"])
            n.create_lines(id="LSR", voltage_level1_id="VS", voltage_level2_id="VR", bus1_id="SB", bus2_id="R",
                           r=1., x=10., g1=0., b1=0., g2=0., b2=0.)
            for t, vl, b in (("TA", "VA", "A"), ("TB", "VB", "B")):
                n.create_2_windings_transformers(id=t, voltage_level1_id="VR", bus1_id="R", voltage_level2_id=vl,
                                                 bus2_id=b, rated_u1=225., rated_u2=20., r=0.5, x=3., g=0., b=0.)
            n.create_generators(id="GS", voltage_level_id="VS", bus_id="SB", target_p=0., target_q=0.,
                                target_v=230., voltage_regulator_on=True, min_p=0., max_p=1000.)
            for gid, vl, b in (("GA", "VA", "A"), ("GB", "VB", "B"), ("GB2", "VB", "B")):
                n.create_generators(id=gid, voltage_level_id=vl, bus_id=b, target_p=20., target_q=0.,
                                    target_v=232., voltage_regulator_on=True, min_p=0., max_p=50.)
            n.create_minmax_reactive_limits(id=["GA", "GB", "GB2"], min_q=[-100., -100., -30.],
                                            max_q=[100., 100., 2.])
            n.create_loads(id="LR", voltage_level_id="VR", bus_id="R", p0=60., q0=-40.)
            n.update_generators(id=["GA", "GB", "GB2"], regulated_element_id=["LR"] * 3)
            # equal keys: inside bus B, GB2 is asked as much as GB, beyond its range
            n.create_extensions("coordinatedReactiveControl", generator_id=["GA", "GB", "GB2"],
                                q_percent=[50.] * 3)
            return n
        params = lf.Parameters(distributed_slack=False)
        n = build()
        self.assertEqual(lf.run_ac(n, params)[0].status, lf.ComponentStatus.CONVERGED)
        q = -n.get_generators().loc["GB2", "q"]
        self.assertAlmostEqual(q, 2., places=6)  # clamped at its limit by the dispatch
        for saturated in (False, True):
            n = build()
            lf.run_ac(n, params)
            bake_outer_loops(n, bake_saturated_voltage_control=saturated)
            self.assertTrue(n.get_generators().loc["GB2", "voltage_regulator_on"], f"saturated={saturated}")

    def test_svc_clamped_inside_a_regulating_bus_not_frozen(self):
        """The SVC shares bus B with GB, both holding the remote load bus: OLF splits the
        bus' reactive power among them by reactive range and clamps the SVC at its small
        limit, the bus (the SVC too) still regulating. As for a generator, neither freezing
        rule freezes it."""
        def build():
            n = pp.network.create_empty()
            n.create_substations(id=["S"])
            n.create_voltage_levels(id=["VR", "VA", "VB", "VS"], substation_id=["S"] * 4,
                                    topology_kind=["BUS_BREAKER"] * 4, nominal_v=[225., 20., 20., 225.])
            n.create_buses(id=["R", "A", "B", "SB"], voltage_level_id=["VR", "VA", "VB", "VS"])
            n.create_lines(id="LSR", voltage_level1_id="VS", voltage_level2_id="VR", bus1_id="SB", bus2_id="R",
                           r=1., x=10., g1=0., b1=0., g2=0., b2=0.)
            for t, vl, b in (("TA", "VA", "A"), ("TB", "VB", "B")):
                n.create_2_windings_transformers(id=t, voltage_level1_id="VR", bus1_id="R", voltage_level2_id=vl,
                                                 bus2_id=b, rated_u1=225., rated_u2=20., r=0.5, x=3., g=0., b=0.)
            n.create_generators(id="GS", voltage_level_id="VS", bus_id="SB", target_p=0., target_q=0.,
                                target_v=230., voltage_regulator_on=True, min_p=0., max_p=1000.)
            for gid, vl, b in (("GA", "VA", "A"), ("GB", "VB", "B")):
                n.create_generators(id=gid, voltage_level_id=vl, bus_id=b, target_p=20., target_q=0.,
                                    target_v=232., voltage_regulator_on=True, min_p=0., max_p=50.)
            n.create_minmax_reactive_limits(id=["GA", "GB"], min_q=[-100., -100.], max_q=[100., 100.])
            # a small range, mostly on the absorbing side: by range, it is asked beyond its max
            n.create_static_var_compensators(id="SVC", voltage_level_id="VB", bus_id="B",
                                             b_min=-30. / 20. ** 2, b_max=2. / 20. ** 2,
                                             regulation_mode="VOLTAGE", target_v=232., target_q=0.,
                                             regulating=True)
            n.create_loads(id="LR", voltage_level_id="VR", bus_id="R", p0=60., q0=40.)
            n.update_generators(id=["GA", "GB"], regulated_element_id=["LR"] * 2)
            n.update_static_var_compensators(id="SVC", regulated_element_id="LR")
            return n
        params = lf.Parameters(distributed_slack=False)
        n = build()
        self.assertEqual(lf.run_ac(n, params)[0].status, lf.ComponentStatus.CONVERGED)
        svc = n.get_static_var_compensators().loc["SVC"]
        v_kv = n.get_buses().loc[svc["bus_id"], "v_mag"]
        self.assertAlmostEqual(-svc["q"], svc["b_max"] * v_kv ** 2, places=6)  # clamped at its limit
        q_gb = -n.get_generators().loc["GB", "q"]
        self.assertTrue(-100. + 1. < q_gb < 100. - 1.)  # GB, on the same bus, has room
        for saturated in (False, True):
            n = build()
            lf.run_ac(n, params)
            bake_outer_loops(n, bake_saturated_voltage_control=saturated)
            self.assertEqual(n.get_static_var_compensators().loc["SVC", "regulation_mode"], "VOLTAGE",
                             f"saturated={saturated}")

    def test_olf_hvdc_saturated_at_its_limit_baked(self):
        """An angle-droop hvdc line OLF saturated at its operator range is baked as a fixed
        setpoint at that limit, in the direction it flowed: the loop-free re-solve then
        reproduces the reference. A line in its linear regime keeps its droop."""
        from test_hvdc_pypowsybl import _build_net, add_operator_range
        loop_free = remove_outer_loops(_with_loops_params())
        # p0 sets the direction of the flow; (12, 8) MW is well below what the droop asks
        for p0, opr, mode, limit in [(20.0, (12.0, 8.0), "SIDE_1_RECTIFIER_SIDE_2_INVERTER", 12.0),
                                     (-200.0, (12.0, 8.0), "SIDE_1_INVERTER_SIDE_2_RECTIFIER", 8.0),
                                     (20.0, (500.0, 500.0), None, None)]:
            def make():
                n = add_operator_range(_build_net(max_p=300.0), *opr)
                n.update_extensions("hvdcAngleDroopActivePowerControl", id="HVDC", p0=p0)
                return n
            n_ref = make()
            lf.run_ac(n_ref, _with_loops_params())
            p_ref = n_ref.get_vsc_converter_stations()["p"]

            n = make()
            lf.run_ac(n, _with_loops_params())
            bake_outer_loops(n)
            droop_on = bool(n.get_extensions("hvdcAngleDroopActivePowerControl").loc["HVDC", "enabled"])
            hvdc = n.get_hvdc_lines().loc["HVDC"]
            if limit is None:
                self.assertTrue(droop_on, "a line in its linear regime must keep its droop")
                continue
            self.assertFalse(droop_on)
            self.assertEqual(hvdc["converters_mode"], mode)
            self.assertAlmostEqual(hvdc["target_p"], limit, places=6)
            self.assertAlmostEqual(p_ref.max(), limit, places=6)  # what the rectifier took
            res = lf.run_ac(n, loop_free)
            self.assertEqual(res[0].status, pp.loadflow.ComponentStatus.CONVERGED)
            p_baked = n.get_vsc_converter_stations()["p"]
            self.assertLess((p_baked - p_ref).abs().max(), 1e-4)


if __name__ == "__main__":
    unittest.main()
