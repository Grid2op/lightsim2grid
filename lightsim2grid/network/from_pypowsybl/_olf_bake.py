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
change the input of the inner solves (*eg* the input of the inner Newton
Raphson) between iterations: a generator that hits its Q limit becomes a
fixed-Q (PQ) injection, a tap changer settles on a discrete position, the
slack mismatch is spread over participating units, etc.

If you take the raw network, run OLF *with* outer loops, and then hand the same
raw network to lightsim2grid, the two engines disagree -- not because the
solvers differ, but because they are solving different problems.

``bake_outer_loops`` rewrites the network's input setpoints to the values the
outer loops settled on, and disables the corresponding regulation, so the
problem becomes a plain power flow. After baking, OLF (loop-free) and
lightsim2grid *should* agree to solver tolerance, and keep agreeing through
topology changes (line/transformer outages) applied identically to both --
this does not always hold in practice (see e.g. multi-root / basin-boundary
cases where the two engines land on different, both-valid, roots).

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
  to fixed-Q at that limit, with voltage regulation switched off. The limit, not
  the reported ``q``, when the two are within 1e-3 MVAr: OLF can report a different
  split of a bus' reactive power among its units than the one it injected; a unit
  further from its limit injects what it reports (see :func:`_baked_q_at_limit`). Units still
  inside their limits keep voltage control -- their equation type is unchanged.
  Units with NaN (unlimited) reactive limits are never frozen, since NaN
  comparisons are False.
* For generators, every one of the voltage-control rules above and below is
  subordinate to what the reference solve shows OLF actually did: a generator
  whose regulated bus sits exactly at its ``target_v`` (to 1e-8 pu, see
  ``_TARGET_V_HELD_TOL_PU``; OLF pins a PV magnitude to machine precision) was
  voltage-controlled and is never frozen, whatever the rules say; a regulating
  generator whose target was NOT held was not controlling (switched, not
  started, discarded, merged, below OLF's ``generatorVoltageControlMinNominalVoltage``,
  ...) and is frozen to its realized Q in a final catch-all pass
  (:func:`_bake_generator_voltage_control_not_held`). On real grid snapshots
  this catch-all is what decides the bulk of the disagreements: it un-freezes
  the generators the rules alone froze although OLF held their target (mostly a
  target-voltage check that used the generator's own terminal nominal voltage
  instead of the regulated bus's) and freezes the ones they left regulating
  although OLF did not. One opt-in exception:
  ``bake_saturated_voltage_control=True`` also freezes a held unit whose Q sits
  at a limit, at that limit.
* Generators OLF's own voltage-control consistency checks would discard for a
  reason other than "not started": too small a reactive range (default
  ``reactiveRangeCheckMode``, widest ``max_q - min_q`` below 1 MVar) or an
  implausible ``target_v`` (outside 0.8-1.2 pu of nominal, on buses above 20 kV).
  Frozen to fixed-Q the same way as the "not started" case above, since OLF
  itself never actually voltage-controlled them either.
* A remotely-regulating generator whose own bus already hosts another
  connected, locally-regulating generator: the shared bus's voltage is already
  pinned by the local unit, so the remote unit's reactive injection cannot
  independently move any voltage (including the one it targets). Frozen to
  fixed-Q the same way, mirroring OLF's own handling of this configuration
  (see :func:`_bake_remote_control_bus_conflicts`).
* The reactive target of a non-regulating generator, where the reference solve injected
  something else: OLF's default ``forceTargetQInReactiveLimits`` clamps it into the
  reactive capability at the generator's active power, a clamp the loop-free parameters
  (reactive limits off) and lightsim2grid do not apply. The realized value is written
  into ``target_q`` (see :func:`_bake_generator_target_q_forced_in_limits`).
* Active-power redistribution from distributed slack / area interchange: the
  realized P is written back into the target P of generators and batteries
  (and optionally loads, if the slack was distributed on load).
* Generators OLF's own ``checkActivePowerControl`` would exclude from
  slack-distribution participation (dispatched at ~0 MW, an implausible
  ``max_p``, ``target_p`` outside ``[min_p, max_p]``, or a degenerate P range):
  their ``activePowerControl`` extension's ``participate`` flag is set to
  ``False`` (created if absent), so neither a subsequent pypowsybl OLF re-solve
  nor lightsim2grid's own distributed slack -- which reads that same flag --
  puts slack mismatch back onto them. So are the generators OLF *capped* while
  distributing the reference mismatch (realized dispatch at ``max_p`` for a
  positive mismatch, at ``min_p`` for a negative one): they took no share, and
  a static distribution that gives them one is wrong by exactly that share.
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

# every OLF-mirrored constant and reading-back tolerance lives in one module, see there
from ._olf_const import (
    _MAX_PLAUSIBLE_ACTIVE_POWER_MW,
    _MAX_PLAUSIBLE_TARGET_V_PU,
    _MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV,
    _MIN_PLAUSIBLE_TARGET_V_PU,
    _MIN_REACTIVE_RANGE_MVAR,
    _Q_LIMIT_TOL_ABS,
    _Q_LIMIT_TOL_REL,
    _TARGET_Q_TOL_MVAR,
    _TARGET_V_HELD_TOL_PU,
    _ZERO_P_TOL,
)


def _q_limit_tol(qmin, qmax):
    # An unset reactive-limit box defaults to +/-Double.MAX in pypowsybl (an
    # "unbounded" sentinel, not a real 3.6e308 MVAr range): subtracting the two
    # overflows to +inf, which is exactly the "no limit" case filtered out below --
    # expected and already handled, so silence the resulting RuntimeWarning.
    with np.errstate(over="ignore"):
        rng = np.asarray(qmax) - np.asarray(qmin)
    # A relative tolerance on that overflowed/absurd range would swallow every
    # generator. Treat absurdly large / non-finite ranges as "no limit" -- no
    # relative bonus, just the absolute floor -- rather than let them dominate max().
    rng = np.where(np.isfinite(rng) & (rng < 1e6), rng, 0.0)
    return np.maximum(_Q_LIMIT_TOL_ABS, _Q_LIMIT_TOL_REL * rng)
            
            
def _get_buses(network):
    """The bus frame every step of the bake reads (solved ``v_mag``, ``voltage_level_id``,
    ``synchronous_component``, plus the ``nominal_v`` of its voltage level).

    Fetched once per bake and handed down: baking rewrites input setpoints only, never
    the solved state nor the topology, so the frame stays valid from the first step to
    the last -- and each fetch is a full round trip to the Java network.
    """
    buses = network.get_buses(attributes=["v_mag", "voltage_level_id", "synchronous_component"])
    nominal_v = network.get_voltage_levels(attributes=["nominal_v"])["nominal_v"]
    buses["nominal_v"] = nominal_v.reindex(buses["voltage_level_id"].to_numpy()).to_numpy(float)
    return buses


def _keep_only_main_comp(df_el, df_bus):
    """
    keep only element (modeled in df_el) that are on the main component => bus_els["synchronous_component"] == 0
    
    This does not deactivate anything.
    """
    comp = df_el["bus_id"].map(df_bus["synchronous_component"])
    return df_el[df_el["connected"].to_numpy(bool) & (comp == 0).to_numpy()]


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


def _hit_qlimit(df: pd.DataFrame, q_gen: pd.Series):
    """The reactive limit (generator convention, same limits as
    :func:`_bound_at_qlimit`) each row's realized output ``q_gen`` sits at: the
    nearer of the two, a missing limit counting as infinitely far.

    A unit OLF switched at its limit injects exactly that limit, but the ``q`` it
    reports is not always it: OLF spreads the reactive target of a bus over the
    generators of that bus again when writing the results, so two units of one bus
    can each be reported a hair away from what they injected, the bus no longer
    balancing. Baking the reported value would then move the injection."""
    qmin, qmax = _reactive_limits(df)
    to_min = (q_gen - qmin).abs().fillna(np.inf)
    to_max = (qmax - q_gen).abs().fillna(np.inf)
    return pd.Series(np.where(to_min < to_max, qmin, qmax), index=df.index)


def _baked_q_at_limit(df: pd.DataFrame, q_gen: pd.Series):
    """The fixed reactive output to bake for a unit frozen at its Q limit: the limit
    (see :func:`_hit_qlimit`) when the reported ``q`` is within ``_Q_LIMIT_TOL_ABS`` of
    it, the reported ``q`` otherwise.

    Only a report that close is the misreport of a limit injection: on real grid
    snapshots the bus balance is off by exactly the reported discrepancy, which is
    what identifies it. The relative tolerance of :func:`_bound_at_qlimit` also
    catches units a visible distance from their limit -- settled inside it, or
    injecting a ``target_q`` beyond it -- and those do inject what they report
    (their bus balances), so the reported value is what must be baked for them."""
    limit = _hit_qlimit(df, q_gen)
    return limit.where((limit - q_gen).abs() <= _Q_LIMIT_TOL_ABS, q_gen)


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
    keep_only_main_comp: bool=True,
    extrapolate_reactive_limits: bool = True,
    bake_saturated_voltage_control: bool = False,
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
        slack-distribution participation for generators and batteries OLF's own
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
    extrapolate_reactive_limits
        Mirror OLF's ``extrapolateReactiveLimits`` provider parameter (on by default):
        a unit whose active power lies outside the P range of its reactive capability
        curve has limits linearly extrapolated from the curve's end segment, where
        pypowsybl's ``min_q_at_p`` / ``max_q_at_p`` clamp to the end point. Set it to
        what the reference solve used, or a unit switched at an extrapolated limit
        reads as a visible distance from its limit (see
        :func:`_extrapolate_curve_limits`).
    bake_saturated_voltage_control
        Also freeze, at its limit, a generator whose reactive output sits at a Q
        limit although the reference solve held its target voltage (a PV unit
        exactly saturated). The base case is unchanged, but any change of the grid
        asking it for more reactive power would make OLF switch it to PQ anyway;
        kept regulating, it reports a reactive-limit violation for almost every
        contingency. ``False`` (default) keeps it regulating, as OLF did.

    Notes
    -----
    Operates in place and is idempotent on an already-baked network. The
    voltage-regulation flag is *not* used as a switch signal: OLF does not flip
    it in IIDM, so PV->PQ is detected from the realized Q sitting at a limit.
    """
    df_bus = _get_buses(network)
    if bake_taps:
        _bake_taps_and_sections(network, keep_only_main_comp, df_bus)
    if bake_reactive_limits:
        _bake_reactive_limit_switches(
            network, keep_only_main_comp, bake_generator_voltage_control_discards,
            extrapolate_reactive_limits, bake_saturated_voltage_control, df_bus
        )
    if bake_active_power:
        _bake_active_power(
            network,
            balance_on_loads=balance_on_loads,
            load_power_factor_constant=load_power_factor_constant,
            keep_only_main_comp=keep_only_main_comp,
            bake_active_power_control_participation=bake_active_power_control_participation,
            df_bus=df_bus
        )
    if bake_remote_voltage_control:
        _bake_remote_voltage_control(network, keep_only_main_comp, df_bus)


def _bake_taps_and_sections(network, keep_only_main_comp=True, df_bus=None):
    # Ratio tap changers (transformer voltage control outer loop).
    df_bus = _get_buses(network) if df_bus is None else df_bus
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


def _resolve_regulated_bus(network, own_bus, regulated_element_id):
    """Bus-view id of the bus each element regulates, as a ``pandas.Series`` indexed
    like ``own_bus``: the element's own bus where it regulates locally, or where the
    remote element cannot be resolved.

    Element-type agnostic (a generator's ``regulated_element_id`` column and a
    battery's ``voltageRegulation`` extension say the same thing), so both go through
    it -- see :func:`_generator_regulated_bus` and :func:`_bake_battery_voltage_control`.
    """
    rel = regulated_element_id.reindex(own_bus.index).fillna("").astype(str)
    reg_bus = own_bus.copy()
    remote = (rel != "") & (rel != rel.index.to_series())
    if not remote.any():
        return reg_bus
    from ._aux_common import _aux_regulated_bus_view_ids
    try:
        resolved = _aux_regulated_bus_view_ids(network, rel[remote].to_numpy())
    except RuntimeError:
        return reg_bus
    resolved = pd.Series(resolved, index=rel.index[remote])
    reg_bus.loc[remote] = resolved.where(resolved != "", own_bus[remote])
    return reg_bus


def _generator_regulated_bus(network, gen):
    """Bus-view id of the bus each generator in ``gen`` regulates (its own bus for
    a local controller), as a ``pandas.Series`` indexed like ``gen``.

    Uses pypowsybl's own ``regulated_bus_id`` when the installed version provides
    it (it resolves a ``regulated_element_id`` that is a busbar section, a
    transformer terminal, another generator, ...), else falls back to resolving
    ``regulated_element_id`` through the bus-connected elements
    (:func:`_resolve_regulated_bus`), then to the generator's own bus.
    """
    own = gen["bus_id"]
    try:
        rb = network.get_generators(attributes=["regulated_bus_id"])["regulated_bus_id"]
        rb = rb.reindex(gen.index)
        return rb.where(rb.notna() & (rb != ""), own)
    except Exception:
        pass
    rel = gen["regulated_element_id"] if "regulated_element_id" in gen.columns else pd.Series("", index=gen.index)
    return _resolve_regulated_bus(network, own, rel)


def _target_v_held(network, reg_bus, target_v, df_bus=None):
    """Boolean ``pandas.Series`` (indexed like ``reg_bus``): the reference solve held
    ``target_v`` at the bus of ``reg_bus`` -- the empirical signature of a unit OLF
    actually voltage-controlled (see ``_TARGET_V_HELD_TOL_PU``). ``False`` where the
    regulated bus has no solved voltage.

    The decision is about the regulated bus, not about the element holding it, so
    generators and batteries share it (their regulation flag and target live in
    different places, what to do with the answer is all that differs).
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
    buses = df_bus.reindex(reg_bus.to_numpy())
    v = buses["v_mag"].to_numpy(float)
    nom = buses["nominal_v"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        dv = np.abs(v - np.asarray(target_v, dtype=float)) / nom
    return pd.Series(np.isfinite(dv) & (dv < _TARGET_V_HELD_TOL_PU), index=reg_bus.index)


def _generator_target_v_held(network, gen, df_bus=None, reg_bus=None):
    """:func:`_target_v_held` for the generators of ``gen``, which must carry
    ``target_v``, ``bus_id`` and (for the fallback resolution) ``regulated_element_id``.
    ``reg_bus``, when given, is :func:`_generator_regulated_bus` already resolved."""
    if reg_bus is None:
        reg_bus = _generator_regulated_bus(network, gen)
    return _target_v_held(network, reg_bus.reindex(gen.index), gen["target_v"], df_bus)


def _generator_regulated_nominal_v(network, gen, df_bus=None, reg_bus=None):
    """Nominal voltage (kV) of the bus each generator in ``gen`` regulates, as a
    numpy array aligned with ``gen`` (NaN where unresolved)."""
    if reg_bus is None:
        reg_bus = _generator_regulated_bus(network, gen)
    df_bus = _get_buses(network) if df_bus is None else df_bus
    return df_bus["nominal_v"].reindex(reg_bus.reindex(gen.index).to_numpy()).to_numpy(float)


def _bake_generator_not_started(network, keep_only_main_comp=True, held=None, df_bus=None):
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

    ``held`` (optional, see :func:`_generator_target_v_held`) is the set of
    generators whose target voltage the reference solve actually held: those are
    never frozen here, whatever this rule says -- the result is the authority on
    what OLF did.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
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
    if held is not None:
        not_started &= ~held.reindex(gen.index).fillna(False).astype(bool)
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


def _bake_generator_voltage_control_discards(network, keep_only_main_comp=True, held=None, df_bus=None,
                                              reg_bus=None):
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

    The target-voltage plausibility is judged in per unit of the **regulated**
    bus's nominal voltage (OLF's ``checkTargetV`` does the same): a 24 kV unit
    regulating the 400 kV side of its step-up transformer legitimately carries a
    ~400 kV ``target_v``, which read as ~16 pu of its own terminal and got it
    wrongly frozen by an earlier version of this check.

    ``held`` (optional, see :func:`_generator_target_v_held`): generators whose
    target the reference solve actually held are never frozen here.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "target_v", "q", "voltage_level_id", "connected", "bus_id",
                    "regulated_element_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    reg = gen[gen["voltage_regulator_on"]]
    if not len(reg):
        return

    reg_nominal_v = _generator_regulated_nominal_v(network, reg, df_bus, reg_bus)

    max_range = _generator_max_reactive_range(network, reg.index).to_numpy()
    too_small_range = max_range < _MIN_REACTIVE_RANGE_MVAR

    with np.errstate(invalid="ignore", divide="ignore"):
        target_v_pu = reg["target_v"].to_numpy() / reg_nominal_v
    implausible_v = (
        (reg_nominal_v > _MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV)
        & ((target_v_pu < _MIN_PLAUSIBLE_TARGET_V_PU) | (target_v_pu > _MAX_PLAUSIBLE_TARGET_V_PU))
    )

    mask = too_small_range | implausible_v
    if held is not None:
        mask &= ~held.reindex(reg.index).fillna(False).to_numpy(bool)
    if not mask.any():
        return
    upd = pd.DataFrame(index=reg.index[mask])
    upd["target_q"] = -reg["q"].to_numpy()[mask]
    upd["voltage_regulator_on"] = False
    network.update_generators(upd)


def _bake_remote_control_bus_conflicts(network, keep_only_main_comp=True, held=None, df_bus=None):
    """Freeze a *remotely*-regulating generator to fixed-Q when its OWN bus
    already hosts another connected, *locally* voltage-regulating generator.

    Physical reason: if a co-located generator locally pins this bus's voltage
    (holding it at its own target, within its own reactive limits), that bus's
    voltage is not a free network unknown from the remote controller's point of
    view -- AC power flow couples buses purely through voltage phasors, so how
    the *combined* reactive injection at that bus splits between the two
    co-located generators has *zero* effect on any other bus's voltage,
    including the one the second generator is trying to remotely control. Its
    own voltage target is therefore unreachable regardless of its dispatched Q.

    Verified on a real grid carrying this exact configuration: solving the
    baked (loop-free) network with OLF, with BOTH generators still nominally
    regulating, reproduces almost exactly the voltage obtained by freezing the
    remote one to its realized Q (much closer than OLF gets to the -- provably
    unreachable -- raw remote target itself), showing OLF itself discards the
    remote controller rather than actually achieving its target through free
    reactive power.

    lightsim2grid's own C++ voltage-control extension has no equivalent
    discard: a lone remote controller sharing a bus with a local one produces a
    genuinely singular Jacobian (the remote controller's own reactive-injection
    unknown ends up with no equation to pair with) which either crashes
    (``ErrorType.SolverFactor``, if the shared bus is classified as a
    distributed-slack participant) or raises ``LSGrid::fill_voltage_control_
    solver_data``'s own "not supported in v1" error (if it lands on an ordinary
    PV bus) -- so this must be resolved before ``init()`` runs.

    Runs after :func:`_bake_generator_not_started` and
    :func:`_bake_generator_voltage_control_discards`: a co-located generator
    already frozen there (e.g. "not started", or discarded for too small a
    reactive range) no longer counts as "locally regulating", so a remote
    generator sharing its bus is correctly left alone once the true conflict is
    resolved.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "regulated_element_id", "q", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    reg = gen["regulated_element_id"].fillna("")
    is_remote = gen["voltage_regulator_on"] & (reg != "") & (reg != gen.index)
    is_local = gen["voltage_regulator_on"] & ~is_remote
    if not is_remote.any() or not is_local.any():
        return
    local_buses = set(gen.loc[is_local, "bus_id"])
    conflict = is_remote & gen["bus_id"].isin(local_buses)
    if held is not None:
        conflict &= ~held.reindex(gen.index).fillna(False).astype(bool)
    if not conflict.any():
        return
    upd = pd.DataFrame(index=gen.index[conflict])
    upd["target_q"] = -gen["q"].to_numpy()[conflict]
    upd["voltage_regulator_on"] = False
    network.update_generators(upd)


def _bake_generator_voltage_control_not_held(network, keep_only_main_comp=True, held=None, df_bus=None):
    """Catch-all: freeze to fixed-Q (its realized reactive output) every
    voltage-regulating generator whose target voltage the reference solve did
    NOT hold (see :func:`_generator_target_v_held`), whatever the reason --
    a reactive-limit switch the tolerance-based rule missed, a generator on a
    bus below ``generatorVoltageControlMinNominalVoltage``, a voltage-control
    discard for a rule not reproduced here, a controller OLF merged into another
    one's ... The realized ``q`` is exactly what OLF injected, so the freeze is
    exact; and a generator OLF *did* control keeps regulating. Only generators
    with a solved (finite) reactive output are touched.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "target_v", "q", "connected", "bus_id", "regulated_element_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if held is None:
        held = _generator_target_v_held(network, gen, df_bus)
    mask = (gen["voltage_regulator_on"] & gen["q"].notna()
            & ~held.reindex(gen.index).fillna(False).astype(bool))
    if not mask.any():
        return
    upd = pd.DataFrame(index=gen.index[mask])
    upd["target_q"] = -gen["q"].to_numpy()[mask]
    upd["voltage_regulator_on"] = False
    network.update_generators(upd)


def _bake_generator_target_q_forced_in_limits(network, keep_only_main_comp=True, df_bus=None):
    """Write the realized reactive output into ``target_q`` for every connected,
    non-regulating generator whose reference solve did not inject its ``target_q``.

    OLF's default ``forceTargetQInReactiveLimits`` clamps a PQ generator's target into
    its reactive capability at its active power (the case seen on real grid snapshots:
    a unit at 0 MW whose curve starts above 0 MVAr, with ``target_q = 0``). The loop-free
    parameters turn the reactive limits off, which turns the clamp off with them, and
    lightsim2grid never clamps: without this both inject the raw target, and the
    voltages of the low-voltage side it feeds are visibly off. The realized value is what
    OLF injected, so the rewrite is exact, and it also covers any other reason a PQ
    injection moved.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen = network.get_generators(attributes=["voltage_regulator_on", "target_q", "q", "connected", "bus_id"])
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    q_gen = -gen["q"]  # result column is load convention; target_q is generator convention
    moved = (~gen["voltage_regulator_on"].astype(bool) & q_gen.notna()
             & ((q_gen - gen["target_q"]).abs() > _TARGET_Q_TOL_MVAR))
    if not moved.any():
        return
    network.update_generators(pd.DataFrame({"target_q": q_gen[moved].to_numpy()}, index=gen.index[moved]))


def _bake_battery_voltage_control(network, keep_only_main_comp=True, df_bus=None):
    """:func:`_bake_generator_voltage_control_not_held`, for the batteries carrying an
    IIDM ``voltageRegulation`` extension (OLF runs them as PV; the converter models them
    as voltage-regulating storage units, see
    ``_aux_add_storage._aux_battery_voltage_regulation``): a battery whose target voltage
    the reference solve did not hold (switched at a reactive limit, discarded, ...) gets
    the extension switched off and its realized reactive output as fixed ``target_q``.

    The decision itself -- :func:`_target_v_held` on the regulated bus -- is the shared
    one. What does not merge is the plumbing around it: a battery's regulation flag,
    target voltage and regulated element live in an extension frame rather than in
    columns of ``get_batteries()``, and switching it off is an ``update_extensions`` call
    rather than an ``update_batteries`` one.
    """
    try:
        vr = network.get_extensions("voltageRegulation")
    except Exception:
        return
    if vr is None or not len(vr) or "voltage_regulator_on" not in vr.columns:
        return
    bat = network.get_batteries(attributes=["q", "connected", "bus_id"])
    if keep_only_main_comp:
        df_bus = _get_buses(network) if df_bus is None else df_bus
        bat = _keep_only_main_comp(bat, df_bus)
    vr = vr.reindex(bat.index)
    on = vr["voltage_regulator_on"].fillna(False).astype(bool)
    if not on.any():
        return
    ids = bat.index[on.to_numpy()]
    reg = vr.loc[ids, "regulated_element_id"] if "regulated_element_id" in vr.columns \
        else pd.Series("", index=ids)
    reg_bus = _resolve_regulated_bus(network, bat.loc[ids, "bus_id"], reg)
    held = _target_v_held(network, reg_bus, vr.loc[ids, "target_v"], df_bus).to_numpy()
    freeze = ~held & bat.loc[ids, "q"].notna().to_numpy()
    if not freeze.any():
        return
    frozen = ids[freeze]
    network.update_extensions("voltageRegulation", pd.DataFrame({"voltage_regulator_on": False}, index=frozen))
    network.update_batteries(pd.DataFrame({"target_q": -bat.loc[frozen, "q"].to_numpy()}, index=frozen))


def _switched_group_members(network, gen, q_gen, reg_bus=None):
    """Boolean ``pandas.Series`` (indexed like ``gen``): a voltage-regulating
    generator whose realized reactive output sits *exactly* at a limit (absolute
    tolerance only, a PQ unit injects its limit to the digit) while at least one
    other generator regulating the same bus is still inside its range -- the
    signature of a unit OLF switched out of a shared control group, whose regulated
    bus nevertheless reads as held."""
    regulating = gen["voltage_regulator_on"].astype(bool)
    qmin, qmax = _reactive_limits(gen)
    exact = regulating & ((q_gen >= qmax - _Q_LIMIT_TOL_ABS) | (q_gen <= qmin + _Q_LIMIT_TOL_ABS))
    if not exact.any():
        return exact
    if reg_bus is None:
        reg_bus = _generator_regulated_bus(network, gen)
    reg_bus = reg_bus.reindex(gen.index)
    n_free = (regulating & ~exact).astype(int).groupby(reg_bus).transform("sum")
    return exact & (n_free > 0)


def _extrapolate_curve_limits(network, df: pd.DataFrame):
    """Return ``df`` with ``min_q_at_p`` / ``max_q_at_p`` replaced, for every row whose
    realized active power (``-p``) lies outside the P range of its reactive capability
    curve, by the linear extrapolation of the curve's end segment on that side.

    That is what OLF evaluates with its default ``extrapolateReactiveLimits``, where
    pypowsybl's columns clamp to the end point. On real grid snapshots, units running
    below the P range of their curve are switched at the *extrapolated* limit, which sits
    a visible distance from the clamped one -- far enough that the freeze below would
    otherwise bake the wrong value. Rows inside their curve, without a curve, or with
    fewer than two points are untouched.
    """
    if not len(df) or "p" not in df.columns or "min_q_at_p" not in df.columns:
        return df
    try:
        pts = network.get_reactive_capability_curve_points()
    except Exception:  # noqa: BLE001 - no curve support in this pypowsybl
        return df
    pts = pts.loc[pts.index.get_level_values(0).isin(df.index)]
    if not len(pts):
        return df
    # one vectorised pass over every curve (a per-curve loop dominates the whole bake on
    # a large grid): sort the points by (element, p), then read each curve's first two
    # and last two points by position
    ids = pts.index.get_level_values(0).to_numpy()
    p_pt = pts["p"].to_numpy(float)
    codes = pd.factorize(ids)[0]
    order = np.lexsort((p_pt, codes))  # stable: points of equal p keep their curve order
    ids, codes, p_pt = ids[order], codes[order], p_pt[order]
    qmin_pt = pts["min_q"].to_numpy(float)[order]
    qmax_pt = pts["max_q"].to_numpy(float)[order]
    start = np.flatnonzero(np.r_[True, codes[1:] != codes[:-1]])
    count = np.diff(np.r_[start, len(ids)])
    start, count = start[count >= 2], count[count >= 2]
    if not len(start):
        return df
    el_ids = ids[start]
    # result column is load convention; curves are in generator convention
    power = -df["p"].reindex(el_ids).to_numpy(float)
    below = power < p_pt[start]
    above = power > p_pt[start + count - 1]
    # the end segment on the side the unit lies (NaN power is neither below nor above)
    i1 = np.where(below, start, start + count - 2)
    i2 = i1 + 1
    p1, p2 = p_pt[i1], p_pt[i2]
    keep = (below | above) & (p2 != p1)
    if not keep.any():
        return df
    i1, i2, p1, p2, power = i1[keep], i2[keep], p1[keep], p2[keep], power[keep]
    df = df.copy()
    rows = df.index.get_indexer(el_ids[keep])
    for col, q in (("min_q_at_p", qmin_pt), ("max_q_at_p", qmax_pt)):
        vals = df[col].to_numpy(float, copy=True)
        vals[rows] = q[i1] + (q[i2] - q[i1]) * (power - p1) / (p2 - p1)
        df[col] = vals
    return df


def _bake_reactive_limit_switches(
    network,
    keep_only_main_comp=True,
    bake_generator_voltage_control_discards=True,
    extrapolate_reactive_limits=True,
    bake_saturated_voltage_control=False,
    df_bus=None):
    # What OLF actually did with each generator's voltage control, read off the
    # reference solve itself: every rule below defers to it (a generator whose
    # target was held is never frozen, one whose target was not is always frozen
    # in the end), the rules only document *why* OLF dropped a control.
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen0 = network.get_generators(
        attributes=["voltage_regulator_on", "target_v", "connected", "bus_id", "regulated_element_id"])
    if keep_only_main_comp:
        gen0 = _keep_only_main_comp(gen0, df_bus)
    # resolved once: which bus a generator regulates is not something the bake changes
    reg_bus = _generator_regulated_bus(network, gen0)
    held = _generator_target_v_held(network, gen0, df_bus, reg_bus) & gen0["voltage_regulator_on"]

    # first, while "not regulating" still means "PQ in the reference solve": the
    # freezes below turn regulation off on units whose reported q is not what they
    # inject (see _hit_qlimit), which this must not write back
    _bake_generator_target_q_forced_in_limits(network, keep_only_main_comp, df_bus)
    _bake_generator_not_started(network, keep_only_main_comp, held, df_bus)
    if bake_generator_voltage_control_discards:
        _bake_generator_voltage_control_discards(network, keep_only_main_comp, held, df_bus, reg_bus)
    _bake_remote_control_bus_conflicts(network, keep_only_main_comp, held, df_bus)
    gen = network.get_generators(
        attributes=[
            "voltage_regulator_on", "q", "p",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id", "regulated_element_id"
        ]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if extrapolate_reactive_limits:
        gen = _extrapolate_curve_limits(network, gen)
    mask, q_gen = _bound_at_qlimit(gen, gen["voltage_regulator_on"])
    # a unit sitting at its Q limit while OLF still held its target is a PV bus
    # exactly saturated, not a switched one: the relative tolerance of
    # _q_limit_tol exists to catch a switched unit that settled a hair *inside*
    # its limit, and that unit's voltage is free, hence off target.
    # "held" is a property of the regulated bus, so it cannot tell apart the
    # members of a group sharing one: a member exactly at its limit while another
    # member is still inside its range was switched (OLF drops it from the group
    # and the others keep the target), and the split it leaves behind decides the
    # voltages behind each unit's step-up transformer.
    # bake_saturated_voltage_control freezes the exactly saturated ones too.
    if not bake_saturated_voltage_control:
        mask &= ~held.reindex(gen.index).fillna(False).astype(bool) | _switched_group_members(network, gen, q_gen, reg_bus)
    if mask.any():
        upd = pd.DataFrame(index=gen.index[mask])
        # the limit it was switched at rather than a misreported q (see _baked_q_at_limit)
        upd["target_q"] = _baked_q_at_limit(gen, q_gen)[mask]
        upd["voltage_regulator_on"] = False
        network.update_generators(upd)
    if bake_generator_voltage_control_discards:
        # everything OLF dropped for a reason the rules above do not spell out
        _bake_generator_voltage_control_not_held(network, keep_only_main_comp, held, df_bus)
    _bake_battery_voltage_control(network, keep_only_main_comp, df_bus)

    vsc = network.get_vsc_converter_stations(
        attributes=[
            "voltage_regulator_on", "q", "p",
            "min_q", "max_q", "min_q_at_p", "max_q_at_p",
            "connected", "bus_id"
        ]
    )
    if keep_only_main_comp:
        vsc = _keep_only_main_comp(vsc, df_bus)
    if extrapolate_reactive_limits:
        vsc = _extrapolate_curve_limits(network, vsc)
    if len(vsc):
        mask, q_gen = _bound_at_qlimit(vsc, vsc["voltage_regulator_on"])
        if mask.any():
            upd = pd.DataFrame(index=vsc.index[mask])
            upd["target_q"] = _baked_q_at_limit(vsc, q_gen)[mask]
            upd["voltage_regulator_on"] = False
            network.update_vsc_converter_stations(upd)

    _bake_svc_standby(network, keep_only_main_comp, df_bus)
    _bake_svc_saturation(network, keep_only_main_comp, df_bus)


def _bake_svc_standby(network, keep_only_main_comp=True, df_bus=None):
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
    df_bus = _get_buses(network) if df_bus is None else df_bus
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


def _bake_svc_saturation(network, keep_only_main_comp=True, df_bus=None):
    """Freeze a VOLTAGE-mode SVC whose realized reactive output sits at (or beyond)
    its voltage-dependent susceptance envelope to fixed-Q (REACTIVE_POWER mode),
    mirroring ``_bake_reactive_limit_switches`` for generators/VSC stations above.

    Unlike a generator's fixed Q box, an SVC's reactive range is
    ``Q(V) = b * V^2`` (``b`` in ``b_min``..``b_max``, in Siemens): recompute
    ``qmin``/``qmax`` in MVAr at the SVC's own solved terminal voltage before
    comparing against the realized ``q``.
    """
    df_bus = _get_buses(network) if df_bus is None else df_bus
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


def _bake_remote_voltage_control(network, keep_only_main_comp=True, df_bus=None):
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
    df_bus = _get_buses(network) if df_bus is None else df_bus
    gen = network.get_generators(
        attributes=["voltage_regulator_on", "regulated_element_id", "connected", "bus_id"]
    )
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    if not len(gen):
        return
    # "remote" matches the converter's own test (see _aux_add_generators.py): a non-empty
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


def _bake_active_power_control_participation(network, gen, bat=None):
    """Zero out active-power (slack-distribution) participation for generators (and
    batteries, see below) that PowSyBl OLF's own
    ``AbstractLfGenerator.checkActivePowerControl`` would exclude:

    * dispatched at (approximately) zero MW (``POWER_EPSILON_SI``, regardless of
      ``min_p`` -- unlike :func:`_bake_generator_not_started`, this check has no
      ``min_p > 0`` guard);
    * an implausible ``max_p`` (``_MAX_PLAUSIBLE_ACTIVE_POWER_MW``, always checked);
    * ``target_p`` outside ``[min_target_p, max_target_p]``, or a degenerate
      (``max_target_p == min_target_p``) range -- OLF's own defaults have
      ``useActiveLimits`` on, so these are always checked too. ``min_target_p`` /
      ``max_target_p`` default to ``min_p`` / ``max_p`` unless the generator's own
      ``activePowerControl`` extension overrides them, mirroring
      ``ActivePowerControlHelper``;
    * **capped** in the reference distribution: OLF spreads the mismatch by the
      sharing key, caps every generator that would cross ``max_target_p`` (for a
      positive mismatch; ``min_target_p`` for a negative one) at that limit and
      redistributes the rest over the others, so a generator whose realized
      dispatch sits at the limit on the side the mismatch pushes it took no share
      at all. That capping is a dynamic effect no static weight reproduces; the
      sign of the reference mismatch -- the total realized minus dispatched
      generation -- is known here, so the generators OLF capped are excluded, and
      the static distribution lightsim2grid then derives from the flag is the one
      OLF effectively used -- on real grid snapshots the units capped this way
      carry a large enough share of the raw key that ignoring them visibly skews
      the distribution. This only holds for a mismatch of the same sign as the
      reference one -- the common case, a loss increase or a load pick-up -- and
      is documented as such in ``_default_distributed_slack``.

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

    Batteries (``bat``, the same pre-bake frame) follow the same rules -- OLF runs
    ``checkActivePowerControl`` on them too -- read off their own extension
    (:func:`~._aux_battery_apc.battery_active_power_control`), with one more detail of
    the capping: a unit is never pushed across 0 MW, so a charging battery is capped at
    0 by a positive mismatch (and a discharging one by a negative mismatch). The sign of
    the reference mismatch is taken over the generators and the batteries together.

    Returns the batteries OLF capped that this pypowsybl cannot mark in their extension
    (<= 1.16.1 rejects a battery id there): the caller pins their active range on their
    realized dispatch (``min_p = max_p``), a degenerate range both OLF and lightsim2grid
    exclude from the slack.
    """
    # not at the top: _aux_battery_apc reads its OLF constants from this module
    from ._aux_battery_apc import (
        _pypowsybl_exposes_battery_apc,
        battery_active_power_control,
        olf_participation_weight,
        olf_target_p_range,
    )
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

    # capped in the reference distribution (see the docstring): realized generation
    # (result "p" is load convention) at the limit on the side the mismatch pushed.
    participate = np.ones(len(gen), bool)
    if len(apc) and "participate" in apc.columns:
        participate = apc["participate"].reindex(gen.index).fillna(True).to_numpy(bool)
    realized = -gen["p"].to_numpy()
    solved = np.isfinite(realized)
    cand = solved & participate & ~excluded

    # the batteries: OLF's rule, read off their own extension, and a range cut at 0 on
    # the side opposite their dispatch (OLF keeps the sign of a unit when it caps it)
    has_bat = bat is not None and len(bat) > 0
    if has_bat:
        b_participate, b_droop, b_min_tp, b_max_tp = battery_active_power_control(network, bat)
        b_target_p = bat["target_p"].to_numpy(float)
        b_min_p = bat["min_p"].to_numpy(float)
        b_max_p = bat["max_p"].to_numpy(float)
        b_weight = olf_participation_weight(b_target_p, b_min_p, b_max_p,
                                            b_participate, b_droop, b_min_tp, b_max_tp)
        b_min_tp, b_max_tp = olf_target_p_range(b_min_p, b_max_p, b_min_tp, b_max_tp)
        b_max_tp = np.where(b_target_p < 0., np.minimum(b_max_tp, 0.), b_max_tp)
        b_min_tp = np.where(b_target_p < 0., b_min_tp, np.maximum(b_min_tp, 0.))
        b_realized = -bat["p"].to_numpy(float)
        b_cand = np.isfinite(b_realized) & (b_weight > 0.)
        b_capped = np.zeros(len(bat), dtype=bool)

    # the sign of the reference mismatch, over every participant
    mismatch = float(np.sum((realized - target_p)[cand]))
    if has_bat:
        mismatch += float(np.sum((b_realized - b_target_p)[b_cand]))
    if mismatch > _ZERO_P_TOL:
        excluded |= cand & (realized >= max_target_p - _ZERO_P_TOL)
        if has_bat:
            b_capped = b_cand & (b_realized >= b_max_tp - _ZERO_P_TOL)
    elif mismatch < -_ZERO_P_TOL:
        excluded |= cand & (realized <= min_target_p + _ZERO_P_TOL)
        if has_bat:
            b_capped = b_cand & (b_realized <= b_min_tp + _ZERO_P_TOL)

    excluded_ids = gen.index[excluded]
    pinned_batteries = bat.index[b_capped] if has_bat else gen.index[:0]
    if has_bat and _pypowsybl_exposes_battery_apc():
        # this pypowsybl writes the extension on a battery as on a generator
        excluded_ids = excluded_ids.append(pinned_batteries)
        pinned_batteries = pinned_batteries[:0]

    if len(excluded_ids):
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
    return pinned_batteries


def _bake_active_power(
    network,
    balance_on_loads,
    load_power_factor_constant,
    keep_only_main_comp=True,
    bake_active_power_control_participation=True,
    df_bus=None
):
    df_bus = _get_buses(network) if df_bus is None else df_bus

    gen = network.get_generators(attributes=["p", "target_p", "min_p", "max_p", "connected", "bus_id"])
    if keep_only_main_comp:
        gen = _keep_only_main_comp(gen, df_bus)
    # the batteries' frame is read before any target_p is overwritten too: they take part
    # in the slack (see _bake_active_power_control_participation)
    bat = network.get_batteries(attributes=["p", "target_p", "min_p", "max_p", "connected", "bus_id"])
    if keep_only_main_comp:
        bat = _keep_only_main_comp(bat, df_bus)
    pinned_batteries = bat.index[:0]
    if bake_active_power_control_participation and (len(gen) or len(bat)):
        pinned_batteries = _bake_active_power_control_participation(network, gen, bat)

    if len(gen):
        # result p is load convention; target_p is generator convention
        network.update_generators(
            pd.DataFrame({"target_p": -gen["p"]}, index=gen.index)
        )

    if len(bat):
        network.update_batteries(
            pd.DataFrame({"target_p": -bat["p"]}, index=bat.index)
        )
        if len(pinned_batteries):
            # a battery OLF capped, on a pypowsybl that cannot mark its extension: a
            # degenerate active range at its realized dispatch takes it out of the slack,
            # in OLF (checkActivePowerControl) and in lightsim2grid alike
            realized = -bat.loc[pinned_batteries, "p"]
            network.update_batteries(
                pd.DataFrame({"min_p": realized, "max_p": realized}, index=pinned_batteries)
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
