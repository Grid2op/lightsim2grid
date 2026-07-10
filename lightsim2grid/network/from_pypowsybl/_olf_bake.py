# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Bake PowSyBl Open Load Flow (OLF) outer-loop results into an IIDM network so
that a subsequent *outer-loop-free* power flow -- in OLF itself or in
lightsim2grid -- reproduces the converged state exactly.

Why this exists
---------------
lightsim2grid solves a single power-flow problem: it does not run the discrete
"outer loops" that OLF runs (distributed slack, reactive-limit PV<->PQ
switching, transformer/shunt voltage control, phase control, ...). Those loops
change the *inputs* of the power flow between inner solves: a generator that
hits its Q limit becomes a fixed-Q (PQ) injection, a tap changer settles on a
discrete position, the slack mismatch is spread over participating units, etc.

If you take the raw network, run OLF *with* outer loops, and then hand the same
raw network to lightsim2grid, the two engines disagree -- not because the
solvers differ, but because they are solving different problems.

``bake_outer_loops`` rewrites the network's input setpoints to the values the
outer loops settled on, and disables the corresponding regulation, so the
problem becomes a plain power flow. After baking, OLF (loop-free) and
lightsim2grid agree to solver tolerance, and they keep agreeing through
topology changes (line/transformer outages) applied identically to both.

The loop-free OLF parameters that reproduce the baked operating point are
available as :func:`lightsim2grid.network.get_pypowsybl_loopfree_parameters`.

What gets baked
---------------
* Discrete tap / section positions (ratio tap, phase tap, shunt section):
  the solved position (``solved_tap_position`` / ``solved_section_count``) is
  copied into the input position and regulation is switched off. The copy is
  guarded: it only happens where the changer was regulating AND a solved value
  exists, because ``solved_tap_position`` is NaN for changers that did not take
  part in an outer loop -- copying NaN would destroy the input tap.
* Generators dispatched at (approximately) 0 MW with a strictly positive
  ``min_p``: OLF never sets these up as voltage-controlling PV buses in the
  first place (``generatorsWithZeroMwTargetAreNotStarted``, "not started"),
  regardless of their ``voltage_regulator_on``/``target_v`` attributes. Frozen
  to fixed-Q (their realized reactive output) before the reactive-limit switch
  below runs. A generator with ``min_p <= 0`` (legitimately allowed to sit at
  0 MW) is left alone.
* PV -> PQ reactive-limit switches (generators, VSC converter stations): only
  the units whose realized reactive power actually sits at a Q limit are frozen
  to fixed-Q at that value, with voltage regulation switched off. Units still
  inside their limits keep voltage control -- their equation type is unchanged.
  Units with NaN (unlimited) reactive limits are never frozen, since NaN
  comparisons are False.
* Generators OLF's own voltage-control consistency checks would discard for a
  reason other than "not started": too small a reactive range (default
  ``reactiveRangeCheckMode``, widest ``max_q - min_q`` below 1 MVar) or an
  implausible ``target_v`` (outside 0.8-1.2 pu of nominal, on buses above 20 kV).
  Frozen to fixed-Q the same way as the "not started" case above, since OLF
  itself never actually voltage-controlled them either.
* Active-power redistribution from distributed slack / area interchange: the
  realized P is written back into the target P of generators and batteries
  (and optionally loads, if the slack was distributed on load).
* Generators OLF's own ``checkActivePowerControl`` would exclude from
  slack-distribution participation (dispatched at ~0 MW, an implausible
  ``max_p``, ``target_p`` outside ``[min_p, max_p]``, or a degenerate P range):
  their ``activePowerControl`` extension's ``participate`` flag is set to
  ``False`` (created if absent), so neither a subsequent pypowsybl OLF re-solve
  nor lightsim2grid's own distributed slack -- which reads that same flag --
  puts slack mismatch back onto them.
* Static var compensators whose realized reactive output sits at (or beyond)
  their voltage-dependent susceptance envelope (``Q(V) = b * V^2``, ``b`` in
  ``b_min``..``b_max``, recomputed in MVAr at the SVC's own solved terminal
  voltage) are frozen to fixed-Q (``REACTIVE_POWER`` mode), mirroring the
  generator/VSC reactive-limit switch above. An SVC still comfortably inside
  its envelope is left regulating (its target reproduces the OLF result
  exactly, since it isn't saturated).
* Static var compensators carrying a "standby automaton" (``standby=True``):
  OLF's own outer loop starts these as a fixed ``b0`` shunt (no voltage
  control) and only switches them to active voltage control once their bus
  voltage crosses a low/high threshold. The static ``regulating`` /
  ``target_v`` attributes do *not* reflect this -- they read as if the SVC
  were always actively regulating. Resolved before the reactive-limit
  switch above: an SVC that stayed within its deadband is frozen to fixed-Q
  (0 when, as commonly configured, ``b0 == 0``); one that crossed a
  threshold is switched to ``VOLTAGE`` mode targeting the corresponding
  ``low_voltage_setpoint`` / ``high_voltage_setpoint`` (then still eligible
  for the ordinary saturation freeze).

Verified against pypowsybl 1.15.0. The OLF-internal round trip (solve with
outer loops -> bake -> solve loop-free) reproduces bus voltages to ~2e-4 kV /
~1e-5 deg, both on IEEE-14 (with a forced reactive-limit switch) and on the
four-substations node-breaker network, which carries VSC and LCC HVDC, an SVC
regulating voltage, a shunt, and ratio + phase tap changers.
"""

import numpy as np
import pandas as pd


# Tolerance (MVAr) for deciding a reactive injection sits "at" its Q limit: the
# larger of a small absolute floor (catches exact/near-exact hits, dominated by
# float rounding) and a fraction of the unit's own Q range (catches OLF's discrete
# reactive-limit outer loop settling a hair inside the limit -- e.g. because P kept
# shifting slightly in later outer-loop iterations after the PV->PQ switch already
# happened). A fixed absolute tolerance alone is either too tight for large-range
# units (missing genuine saturation) or too loose for small-range ones (freezing
# units that still have real headroom).
_Q_LIMIT_TOL_ABS = 1e-3
_Q_LIMIT_TOL_REL = 0.005  # 0.5% of (qmax - qmin)


def _q_limit_tol(qmin, qmax):
    rng = np.asarray(qmax) - np.asarray(qmin)
    # An unset reactive-limit box defaults to +/-Double.MAX in pypowsybl (an
    # "unbounded" sentinel, not a real 3.6e308 MVAr range): a relative tolerance on
    # that would swallow every generator. Treat absurdly large / non-finite ranges
    # as "no limit" -- no relative bonus, just the absolute floor -- rather than
    # let them dominate the max().
    rng = np.where(np.isfinite(rng) & (rng < 1e6), rng, 0.0)
    return np.maximum(_Q_LIMIT_TOL_ABS, _Q_LIMIT_TOL_REL * rng)
            
            
def _keep_only_main_comp(df_el, df_bus):
    """
    keep only element (modeled in df_el) that are on the main component => bus_els["synchronous_component"] == 0
    
    This does not deactivate anything.
    """
    mask_conn = df_el["connected"]
    bus_els = df_bus.loc[df_el.loc[mask_conn, "bus_id"]]
    mask_main = (bus_els["synchronous_component"] == 0).to_numpy()
    df_el = df_el.loc[df_el.loc[mask_conn][mask_main].index]
    return df_el


def _reactive_limits(df: pd.DataFrame):
    """Return (qmin, qmax) in *generator* convention.

    Prefer the P-dependent capability-curve limits (``min_q_at_p`` /
    ``max_q_at_p``) when present and finite; fall back to the fixed
    ``min_q`` / ``max_q`` box otherwise.
    """
    if "min_q_at_p" in df.columns:
        qmin = df["min_q_at_p"].where(df["min_q_at_p"].notna(), df["min_q"])
        qmax = df["max_q_at_p"].where(df["max_q_at_p"].notna(), df["max_q"])
    else:
        qmin, qmax = df["min_q"], df["max_q"]
    return qmin, qmax


def _bound_at_qlimit(df: pd.DataFrame, regulating: pd.Series):
    """Boolean mask of rows that regulate voltage *and* sit at a Q limit,
    plus the realized reactive power in generator convention.

    ``df`` must carry result column ``q`` (load convention) and the reactive
    limit columns. ``regulating`` is the per-row voltage-regulation flag.
    """
    q_gen = -df["q"]  # result column is load convention; flip to generator
    qmin, qmax = _reactive_limits(df)
    tol = _q_limit_tol(qmin, qmax)
    at_max = regulating & (q_gen >= qmax - tol)
    at_min = regulating & (q_gen <= qmin + tol)
    return (at_max | at_min), q_gen


def bake_outer_loops(
    network,
    bake_taps: bool = True,
    bake_reactive_limits: bool = True,
    bake_generator_voltage_control_discards: bool = True,
    bake_active_power: bool = True,
    bake_active_power_control_participation: bool = True,
    bake_remote_voltage_control: bool = False,
    balance_on_loads: bool = False,
    load_power_factor_constant: bool = False,
    keep_only_main_comp: bool=True
):
    """Rewrite ``network`` input setpoints to the converged outer-loop state.

    Call this on a network that has just been solved by OLF *with* outer loops.
    Afterwards the network represents a plain power-flow problem: a loop-free
    OLF run (see :func:`get_pypowsybl_loopfree_parameters`) or a lightsim2grid
    run (via :func:`init_from_pypowsybl`) will reproduce the same operating
    point.

    Parameters
    ----------
    network
        A pypowsybl network, freshly solved with the outer loops enabled.
    bake_taps
        Copy solved ratio/phase tap positions and shunt sections into the
        input positions and disable their regulation.
    bake_reactive_limits
        Freeze generators / VSC stations that hit a Q limit to fixed-Q (PQ).
    bake_generator_voltage_control_discards
        Freeze generators OLF's own voltage-control consistency checks would
        discard for a reason other than "not started": too small a reactive
        range, or an implausible ``target_v`` (see
        :func:`_bake_generator_voltage_control_discards`). Also gated by
        ``bake_reactive_limits`` -- has no effect if that is off.
    bake_active_power
        Write realized active power back into generator/battery target P
        (and load p0/q0 if ``balance_on_loads``).
    bake_active_power_control_participation
        Zero out (``activePowerControl`` extension ``participate=False``)
        slack-distribution participation for generators OLF's own
        ``checkActivePowerControl`` would exclude (see
        :func:`_bake_active_power_control_participation`). Also gated by
        ``bake_active_power`` -- has no effect if that is off.
    bake_remote_voltage_control
        Rewrite remote voltage control to local control at the solved terminal
        voltage (see :func:`_bake_remote_voltage_control`). Needed so that
        remote-regulating generators can sit on a (distributed) slack bus, which
        lightsim2grid v1 does not otherwise support.
    balance_on_loads
        Set if the slack was distributed on loads (BalanceType
        PROPORTIONAL_TO_LOAD / CONFORM_LOAD).
    load_power_factor_constant
        Mirror OLF's ``loadPowerFactorConstant``: also rewrite load q0 so the
        power factor is preserved.
    keep_only_main_comp
        Only elements of the main connected component are updated (True by default)

    Notes
    -----
    Operates in place and is idempotent on an already-baked network. The
    voltage-regulation flag is *not* used as a switch signal: OLF does not flip
    it in IIDM, so PV->PQ is detected from the realized Q sitting at a limit.
    """
    if bake_taps:
        _bake_taps_and_sections(network, keep_only_main_comp)
    if bake_reactive_limits:
        _bake_reactive_limit_switches(
            network, keep_only_main_comp, bake_generator_voltage_control_discards
        )
    if bake_active_power:
        _bake_active_power(
            network,
            balance_on_loads=balance_on_loads,
            load_power_factor_constant=load_power_factor_constant,
            keep_only_main_comp=keep_only_main_comp,
            bake_active_power_control_participation=bake_active_power_control_participation
        )
    if bake_remote_voltage_control:
        _bake_remote_voltage_control(network, keep_only_main_comp)


def _bake_taps_and_sections(network, keep_only_main_comp=True):
    # Ratio tap changers (transformer voltage control outer loop).
    df_bus = network.get_buses(attributes=["synchronous_component"])
    rtc = network.get_ratio_tap_changers(
        attributes=["tap", "solved_tap_position", "regulating"]
    )
    if len(rtc):
        # Only adopt the solved position where the changer was actually
        # regulating AND a solved position exists. ``solved_tap_position`` is
        # NaN when the changer did not participate in an outer loop; copying it
        # blindly would destroy the input tap.
        keep = rtc["regulating"] & rtc["solved_tap_position"].notna()
        if keep.any():
            upd = pd.DataFrame(index=rtc.index[keep])
            upd["tap"] = rtc["solved_tap_position"][keep].astype(int)
            # upd["regulating"] = False
            network.update_ratio_tap_changers(upd)

    # Phase tap changers (phase-shifter regulation outer loop).
    ptc = network.get_phase_tap_changers(
        attributes=["tap", "solved_tap_position", "regulating"]
    )
    if len(ptc):
        keep = ptc["regulating"] & ptc["solved_tap_position"].notna()
        if keep.any():
            upd = pd.DataFrame(index=ptc.index[keep])
            upd["tap"] = ptc["solved_tap_position"][keep].astype(int)
            # upd["regulating"] = False
            network.update_phase_tap_changers(upd)

    # Shunt compensators (shunt voltage control outer loop).
    sh = network.get_shunt_compensators(
        attributes=[
            "section_count", 
            "solved_section_count", 
            "voltage_regulation_on",
            "connected", 
            "bus_id"
        ]
    )
    if keep_only_main_comp:
        sh = _keep_only_main_comp(sh, df_bus)
    if len(sh):
        keep = sh["voltage_regulation_on"] & sh["solved_section_count"].notna()
        if keep.any():
            upd = pd.DataFrame(index=sh.index[keep])
            upd["section_count"] = sh["solved_section_count"][keep].astype(int)
            upd["voltage_regulation_on"] = False
            network.update_shunt_compensators(upd)


# Tolerance (MW) for deciding a generator's target_p is "zero": mirrors OLF's own
# POWER_EPSILON_SI = 1e-4 MW (AbstractLfGenerator.checkIfGeneratorStartedForVoltageControl).
_ZERO_P_TOL = 1e-4

# OLF's own PlausibleValues.MIN_REACTIVE_RANGE (1 MVar): with the default
# reactiveRangeCheckMode ("MAX"), a generator whose widest reactive range across its
# whole active-power range (or its fixed min_q/max_q box, for a MIN_MAX-kind generator)
# falls below this is discarded from voltage control
# (AbstractLfGenerator.checkIfReactiveRangesAreLargeEnoughForVoltageControl).
_MIN_REACTIVE_RANGE_MVAR = 1.0

# OLF's own defaults (LfNetworkParameters): a generator's targetV, expressed in per
# unit of its (regulated bus) nominal voltage, is discarded from voltage control if
# outside this range -- but only on buses above the nominal-voltage floor below (this
# is minNominalVoltageTargetVoltageCheck, distinct from the *realistic-voltage-check*
# floor that PARAMS_STANDARD sets to 180 kV; this one is left at its own OLF default).
_MIN_PLAUSIBLE_TARGET_V_PU = 0.8
_MAX_PLAUSIBLE_TARGET_V_PU = 1.2
_MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV = 20.0

# OLF's own plausibleActivePowerLimit default (MW): a generator whose maxP exceeds
# this is discarded from active-power (slack) participation regardless of any other
# parameter (AbstractLfGenerator.checkActivePowerControl).
_MAX_PLAUSIBLE_ACTIVE_POWER_MW = 10000.0


def _bake_generator_not_started(network, keep_only_main_comp=True):
    """Freeze a voltage-regulating generator dispatched at (approximately) 0 MW,
    with a strictly positive minimum active power, to fixed-Q (PQ) -- mirroring
    PowSyBl OpenLoadFlow's own generator setup rule (default-on parameter
    ``generatorsWithZeroMwTargetAreNotStarted``,
    ``AbstractLfGenerator.checkIfGeneratorStartedForVoltageControl``): such a
    generator is declared unable to legitimately run at 0 MW (``min_p > 0``), so
    OLF treats it as "not started" and never sets it up as a voltage-controlling
    PV bus in the first place -- regardless of its ``voltage_regulator_on`` /
    ``target_v`` attributes, which read as if it were actively regulating. A
    generator with ``min_p <= 0`` (legitimately allowed to sit at 0 MW, e.g. a
    curtailed renewable) is NOT affected -- OLF keeps it on voltage control.

    ``lightsim2grid``'s own C++ analogue (``GeneratorContainer::is_pseudo_off``,
    gated by ``turnedoff_no_pv`` / ``LightSimBackend(turned_off_pv=False)``) does
    not check ``min_p`` yet (see the CHANGELOG.rst TODO) and is off by default
    anyway (Grid2Op semantics: an agent redispatching a generator to ~0 MW
    mid-episode should not silently also drop its voltage support) -- so this
    OLF-specific case is instead resolved here, on the pypowsybl-loading path
    only, the same way ``_bake_svc_standby`` resolves the SVC standby automaton.

    Runs before the reactive-limit-switch check below: a generator resolved to PQ
    here is already frozen and skipped there (it only looks at
    ``voltage_regulator_on == True``).
    """
    df_bus = network.get_buses(attributes=["synchronous_component"])
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "target_p", "min_p", "q", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    not_started = (
        gen["voltage_regulator_on"]
        & (gen["target_p"].abs() < _ZERO_P_TOL)
        & (gen["min_p"] > _ZERO_P_TOL)
    )
    if not not_started.any():
        return
    upd = pd.DataFrame(index=gen.index[not_started])
    # result "q" is load/receptor convention; Generator.target_q is generator convention
    upd["target_q"] = -gen["q"].to_numpy()[not_started]
    upd["voltage_regulator_on"] = False
    network.update_generators(upd)


def _generator_max_reactive_range(network, gen_index):
    """OLF's own default ``reactiveRangeCheckMode`` is ``MAX``: the widest
    ``max_q - min_q`` across the whole active-power range of the reactive
    capability curve for a CURVE-kind generator, or simply ``max_q - min_q``
    for a MIN_MAX (fixed box) generator. Returns a ``pandas.Series`` of ranges
    (MVAr), indexed like ``gen_index``, ``NaN`` for any id not found.
    """
    rng = pd.Series(np.nan, index=gen_index)
    if not len(gen_index):
        return rng
    box = network.get_generators(attributes=["reactive_limits_kind", "min_q", "max_q"])
    box = box.loc[box.index.intersection(gen_index)]
    is_box = box["reactive_limits_kind"] == "MIN_MAX"
    rng.loc[box.index[is_box]] = (box["max_q"] - box["min_q"])[is_box]
    curve_ids = box.index[~is_box]
    if len(curve_ids):
        pts = network.get_reactive_capability_curve_points()
        # pts is indexed on a (id, num) MultiIndex -- intersecting it directly against
        # a flat Index of generator ids matches nothing, since a tuple never equals a
        # bare id; filter on the first level instead.
        pts = pts.loc[pts.index.get_level_values(0).isin(curve_ids)]
        if len(pts):
            pt_range = pts["max_q"] - pts["min_q"]
            rng.update(pt_range.groupby(level=0).max())
    return rng


def _bake_generator_voltage_control_discards(network, keep_only_main_comp=True):
    """Freeze a voltage-regulating generator to fixed-Q when PowSyBl OLF's own
    consistency checks (``AbstractLfGenerator.checkVoltageControlConsistency``)
    would have discarded it from voltage control for a reason other than "not
    started" (handled separately by :func:`_bake_generator_not_started`, which
    must run before this so its frozen generators are excluded here):

    * too small a reactive range (``checkIfReactiveRangesAreLargeEnoughForVoltageControl``,
      see :func:`_generator_max_reactive_range` and ``_MIN_REACTIVE_RANGE_MVAR``);
    * an implausible target voltage (``checkTargetV``, see
      ``_MIN_PLAUSIBLE_TARGET_V_PU`` / ``_MAX_PLAUSIBLE_TARGET_V_PU``).

    Like ``_bake_generator_not_started``, these generators never actually reached
    voltage control in the reference (outer-loop) run: OLF falls back to a fixed-Q
    injection at the raw IIDM ``target_q`` (``LfGeneratorImpl.getTargetQ()``), *not*
    at any computed value -- so the realized ``q`` from the reference solve already
    equals that fixed output, and freezing to it here is exact. Without this,
    lightsim2grid -- which has no equivalent of these OLF-specific plausibility
    screens -- would treat such a generator as a normal PV bus after baking.
    """
    df_bus = network.get_buses(attributes=["synchronous_component"])
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "target_v", "q", "voltage_level_id", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    reg = gen[gen["voltage_regulator_on"]]
    if not len(reg):
        return

    nominal_v = network.get_voltage_levels(attributes=["nominal_v"])["nominal_v"]
    # approximation: uses the generator's own terminal nominal voltage, not the
    # (possibly remote) regulated bus -- same base-case approximation already
    # documented for _bake_remote_voltage_control.
    gen_nominal_v = nominal_v.reindex(reg["voltage_level_id"]).to_numpy()

    max_range = _generator_max_reactive_range(network, reg.index).to_numpy()
    too_small_range = max_range < _MIN_REACTIVE_RANGE_MVAR

    target_v_pu = reg["target_v"].to_numpy() / gen_nominal_v
    implausible_v = (
        (gen_nominal_v > _MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV)
        & ((target_v_pu < _MIN_PLAUSIBLE_TARGET_V_PU) | (target_v_pu > _MAX_PLAUSIBLE_TARGET_V_PU))
    )

    mask = too_small_range | implausible_v
    if not mask.any():
        return
    upd = pd.DataFrame(index=reg.index[mask])
    upd["target_q"] = -reg["q"].to_numpy()[mask]
    upd["voltage_regulator_on"] = False
    network.update_generators(upd)


def _bake_reactive_limit_switches(
    network,
    keep_only_main_comp=True,
    bake_generator_voltage_control_discards=True):
    _bake_generator_not_started(network, keep_only_main_comp)
    if bake_generator_voltage_control_discards:
        _bake_generator_voltage_control_discards(network, keep_only_main_comp)
    df_bus = network.get_buses(attributes=["synchronous_component"])
    gen = network.get_generators(
        attributes=[
            "voltage_regulator_on", "q",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id"
        ]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    mask, q_gen = _bound_at_qlimit(gen, gen["voltage_regulator_on"])
    if mask.any():
        upd = pd.DataFrame(index=gen.index[mask])
        upd["target_q"] = q_gen[mask]
        upd["voltage_regulator_on"] = False
        network.update_generators(upd)
        
    vsc = network.get_vsc_converter_stations(
        attributes=[
            "voltage_regulator_on", "q",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id"
        ]
    )
    if keep_only_main_comp:
        vsc = _keep_only_main_comp(vsc, df_bus)
    if len(vsc):
        mask, q_gen = _bound_at_qlimit(vsc, vsc["voltage_regulator_on"])
        if mask.any():
            upd = pd.DataFrame(index=vsc.index[mask])
            upd["target_q"] = q_gen[mask]
            upd["voltage_regulator_on"] = False
            network.update_vsc_converter_stations(upd)

    _bake_svc_standby(network, keep_only_main_comp)
    _bake_svc_saturation(network, keep_only_main_comp)


def _bake_svc_standby(network, keep_only_main_comp=True):
    """Resolve an SVC's "standby automaton" (PowSyBl OLF's ``MonitoringVoltageOuterLoop``)
    to the state the outer loop actually settled on.

    A standby SVC (``standbyAutomaton`` extension, ``standby=True``) starts the outer loop
    as a fixed-susceptance shunt (``b0``, no voltage control) rather than a normal
    voltage-regulating device -- regardless of what its static ``regulating`` /
    ``regulation_mode`` / ``target_v`` attributes say. It only switches to active PV control
    -- targeting ``low_voltage_setpoint`` / ``high_voltage_setpoint`` -- once its own bus
    voltage crosses ``low_voltage_threshold`` / ``high_voltage_threshold``. Inside the
    deadband it stays a plain ``b0`` shunt: the realized ``q`` there already equals
    ``b0 * v_mag_kv ** 2`` (0 when, as commonly configured, ``b0 == 0``), so freezing it to
    that realized value -- like a saturated SVC -- reproduces the deadband state exactly.

    Runs before ``_bake_svc_saturation``: an SVC resolved to VOLTAGE mode here (crossed a
    threshold) is still eligible for the ordinary saturation freeze; an SVC resolved to
    REACTIVE_POWER here is already frozen and the saturation check leaves it alone (it only
    looks at ``regulation_mode == "VOLTAGE"``).
    """
    try:
        automaton = network.get_extensions("standbyAutomaton")
    except Exception:  # noqa: BLE001 - extension unsupported / absent on old pypowsybl
        return
    automaton = automaton[automaton["standby"]]
    if not len(automaton):
        return
    df_bus = network.get_buses(attributes=["v_mag", "synchronous_component"])
    svc = network.get_static_var_compensators(attributes=["connected", "bus_id"])
    svc = svc.loc[svc.index.intersection(automaton.index)]
    if keep_only_main_comp:
        svc = _keep_only_main_comp(svc, df_bus)
    automaton = automaton.loc[svc.index]
    if not len(automaton):
        return

    v_kv = df_bus.loc[svc["bus_id"].values, "v_mag"].to_numpy()
    above_high = v_kv > automaton["high_voltage_threshold"].to_numpy()
    below_low = v_kv < automaton["low_voltage_threshold"].to_numpy()

    active = above_high | below_low
    if active.any():
        upd = pd.DataFrame(index=automaton.index[active])
        upd["target_v"] = np.where(
            above_high[active],
            automaton["high_voltage_setpoint"].to_numpy()[active],
            automaton["low_voltage_setpoint"].to_numpy()[active],
        )
        upd["regulation_mode"] = "VOLTAGE"
        network.update_static_var_compensators(upd)

    idle_ids = automaton.index[~active]
    if len(idle_ids):
        # realized "q" (receptor convention, matches target_q) already equals b0 * v^2
        svc_q = network.get_static_var_compensators(attributes=["q"]).loc[idle_ids, "q"]
        upd = pd.DataFrame(index=idle_ids)
        upd["target_q"] = svc_q
        upd["regulation_mode"] = "REACTIVE_POWER"
        network.update_static_var_compensators(upd)


def _bake_svc_saturation(network, keep_only_main_comp=True):
    """Freeze a VOLTAGE-mode SVC whose realized reactive output sits at (or beyond)
    its voltage-dependent susceptance envelope to fixed-Q (REACTIVE_POWER mode),
    mirroring ``_bake_reactive_limit_switches`` for generators/VSC stations above.

    Unlike a generator's fixed Q box, an SVC's reactive range is
    ``Q(V) = b * V^2`` (``b`` in ``b_min``..``b_max``, in Siemens): recompute
    ``qmin``/``qmax`` in MVAr at the SVC's own solved terminal voltage before
    comparing against the realized ``q``.
    """
    df_bus = network.get_buses(attributes=["v_mag", "synchronous_component"])
    svc = network.get_static_var_compensators(
        attributes=["regulating", "regulation_mode", "q", "b_min", "b_max", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        svc = _keep_only_main_comp(svc, df_bus)
    is_voltage = svc["regulating"] & (svc["regulation_mode"] == "VOLTAGE")
    svc = svc[is_voltage]
    if not len(svc):
        return
    v_kv = df_bus.loc[svc["bus_id"].values, "v_mag"].to_numpy()
    # SVC "q" (like generators') is the terminal/receptor-convention result: flip to
    # generator/injection convention to compare against the susceptance envelope.
    q_gen = -svc["q"].to_numpy()
    qmax = svc["b_max"].to_numpy() * v_kv ** 2  # Q[MVAr] = B[S] * V[kV]^2
    qmin = svc["b_min"].to_numpy() * v_kv ** 2
    tol = _q_limit_tol(qmin, qmax)
    mask = (q_gen >= qmax - tol) | (q_gen <= qmin + tol)
    if mask.any():
        upd = pd.DataFrame(index=svc.index[mask])
        # unlike Generator.target_q (already generator convention), StaticVarCompensator
        # .target_q is receptor convention like its "q" result column -- write the raw
        # (unflipped) realized value.
        upd["target_q"] = svc["q"].to_numpy()[mask]
        upd["regulation_mode"] = "REACTIVE_POWER"
        network.update_static_var_compensators(upd)


def _bake_remote_voltage_control(network, keep_only_main_comp=True):
    """Rewrite *remote* voltage control into *local* control at the solved terminal.

    A generator regulating a bus other than its own terminal (``regulated_element_id``
    resolving to a different bus -- e.g. holding the 400 kV grid-connection point
    across its step-up transformer) is rewritten to regulate its OWN terminal at that
    terminal's solved magnitude, and the remote regulation is cleared so the converter
    treats it as local control.

    At the converged operating point this is exact: the generator's own terminal
    already sits at ``v_mag`` with the realized reactive output, so fixing it locally
    reproduces the same fixed point (same V everywhere, same Q; the formerly-regulated
    remote bus still lands on its solved value). A subsequent loop-free solve -- OLF
    (whose loop-free parameters disable remote control anyway) or lightsim2grid --
    reproduces the baked state.

    Why it exists: lightsim2grid v1 cannot host a remote voltage controller on a slack
    bus, and the default distributed slack puts many remote-regulating generators on
    slack buses. Making every controller local sidesteps that.

    .. warning::
        This is a *base-case* faithful approximation. Under a topology change the
        generator then holds its own terminal magnitude instead of the remote bus, so
        the post-contingency reactive behaviour differs from true remote control.

    Operates in place; idempotent (an already-local generator is left untouched).
    """
    df_bus = network.get_buses(attributes=["v_mag", "synchronous_component"])
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "regulated_element_id", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if not len(gen):
        return
    # "remote" matches the converter's own test (see _from_pypowsybl.init): a non-empty
    # regulated element that is not the generator's own id.
    reg = gen["regulated_element_id"].fillna("")
    remote = gen["voltage_regulator_on"] & gen["connected"] & (reg != "") & (reg != gen.index)
    gen = gen[remote]
    if not len(gen):
        return
    upd = pd.DataFrame(index=gen.index)
    upd["target_v"] = df_bus.loc[gen["bus_id"].values, "v_mag"].values  # own terminal, kV
    upd["regulated_element_id"] = gen.index  # regulate own terminal -> local control
    network.update_generators(upd)


def _bake_active_power_control_participation(network, gen):
    """Zero out active-power (slack-distribution) participation for generators that
    PowSyBl OLF's own ``AbstractLfGenerator.checkActivePowerControl`` would exclude:

    * dispatched at (approximately) zero MW (``POWER_EPSILON_SI``, regardless of
      ``min_p`` -- unlike :func:`_bake_generator_not_started`, this check has no
      ``min_p > 0`` guard);
    * an implausible ``max_p`` (``_MAX_PLAUSIBLE_ACTIVE_POWER_MW``, always checked);
    * ``target_p`` outside ``[min_target_p, max_target_p]``, or a degenerate
      (``max_target_p == min_target_p``) range -- OLF's own defaults have
      ``useActiveLimits`` on, so these are always checked too. ``min_target_p`` /
      ``max_target_p`` default to ``min_p`` / ``max_p`` unless the generator's own
      ``activePowerControl`` extension overrides them, mirroring
      ``ActivePowerControlHelper``.

    ``gen`` must be the pre-bake generator frame (``p``, ``target_p``, ``min_p``,
    ``max_p``), read *before* :func:`_bake_active_power` overwrites ``target_p``
    with the realized dispatch -- the discard decision is about what OLF actually
    used to build the reference solve, not the post-bake value.

    Writes ``participate=False`` into the network's own ``activePowerControl``
    extension (created where absent) rather than any lightsim2grid-side state:
    lightsim2grid's own ``_default_distributed_slack`` (mirroring the same OLF
    rule) already reads exactly this extension's ``participate`` flag when
    computing slack weights, so this single write keeps both lightsim2grid and any
    subsequent pypowsybl OLF re-solve from putting slack mismatch back onto a
    generator OLF itself excluded from participating.
    """
    apc = network.get_extensions("activePowerControl")
    min_target_p = gen["min_p"].to_numpy(copy=True)
    max_target_p = gen["max_p"].to_numpy(copy=True)
    if len(apc):
        common = gen.index.intersection(apc.index)
        if len(common):
            pos = gen.index.get_indexer(common)
            if "min_target_p" in apc.columns:
                v = apc.loc[common, "min_target_p"].to_numpy()
                ok = ~np.isnan(v)
                min_target_p[pos[ok]] = v[ok]
            if "max_target_p" in apc.columns:
                v = apc.loc[common, "max_target_p"].to_numpy()
                ok = ~np.isnan(v)
                max_target_p[pos[ok]] = v[ok]

    target_p = gen["target_p"].to_numpy()
    max_p = gen["max_p"].to_numpy()

    zero_target = np.abs(target_p) < _ZERO_P_TOL
    maxp_not_plausible = max_p > _MAX_PLAUSIBLE_ACTIVE_POWER_MW
    outside_limits = (target_p > max_target_p) | (target_p < min_target_p)
    degenerate_range = (max_target_p - min_target_p) < _ZERO_P_TOL

    excluded = zero_target | maxp_not_plausible | outside_limits | degenerate_range
    if not excluded.any():
        return

    excluded_ids = gen.index[excluded]
    already_apc = excluded_ids.intersection(apc.index) if len(apc) else excluded_ids[:0]
    new_apc = excluded_ids.difference(already_apc)

    if len(already_apc):
        network.update_extensions(
            "activePowerControl", pd.DataFrame({"participate": False}, index=already_apc)
        )
    if len(new_apc):
        network.create_extensions(
            "activePowerControl", pd.DataFrame({"participate": False}, index=new_apc)
        )


def _bake_active_power(
    network,
    balance_on_loads,
    load_power_factor_constant,
    keep_only_main_comp=True,
    bake_active_power_control_participation=True
):
    df_bus = network.get_buses(attributes=["synchronous_component"])

    gen = network.get_generators(attributes=["p", "target_p", "min_p", "max_p", "connected", "bus_id"])
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if len(gen):
        if bake_active_power_control_participation:
            _bake_active_power_control_participation(network, gen)
        # result p is load convention; target_p is generator convention
        network.update_generators(
            pd.DataFrame({"target_p": -gen["p"]}, index=gen.index)
        )

    bat = network.get_batteries(attributes=["p", "connected", "bus_id"])
    if keep_only_main_comp:
        bat = _keep_only_main_comp(bat, df_bus)
    if len(bat):
        network.update_batteries(
            pd.DataFrame({"target_p": -bat["p"]}, index=bat.index)
        )

    if balance_on_loads:
        load = network.get_loads(attributes=["p", "q", "connected", "bus_id"])
        if keep_only_main_comp:
            load = _keep_only_main_comp(load, df_bus)
        if len(load):
            upd = pd.DataFrame(index=load.index)
            upd["p0"] = load["p"]  # load convention both sides, no flip
            if load_power_factor_constant:
                upd["q0"] = load["q"]
            network.update_loads(upd)
