# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""OpenLoadFlow's network-loading rules, as pure functions of the pypowsybl network.

When OpenLoadFlow builds its own network from an IIDM one (``LfNetworkLoaderImpl``,
``AbstractLfGenerator``, ``LfGeneratorImpl``), it decides which units actually regulate a
voltage and what a non-regulating unit injects -- before any powerflow, and whatever the
outer loops do afterwards. This module is the one place lightsim2grid reproduces those
decisions:

* ``init_from_pypowsybl(..., olf_rules=True)`` applies them while building the grid;
* :func:`bake_outer_loops` reads them to tell a unit OpenLoadFlow never let regulate from one
  an outer loop switched.

Nothing here modifies the network, and nothing reads a powerflow result: every rule is a
function of the input data only. The reference is OpenLoadFlow 2.3.0 with the default
parameters of the pypowsybl used for the comparisons (``OlfLoadingParameters``), see
``docs/dev_notes/outer_loops_fixed_sparsity.md``.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ._olf_const import (
    _MAX_PLAUSIBLE_ACTIVE_POWER_MW,
    _MAX_PLAUSIBLE_TARGET_V_PU,
    _MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV,
    _MIN_PLAUSIBLE_TARGET_V_PU,
    _MIN_REACTIVE_RANGE_MVAR,
    _NOT_STARTED_P_TOL_MW,
    _OLF_DEFAULT_DROOP,
    _TARGET_V_EPSILON_PU,
    _ZERO_P_TOL,
)


@dataclass(frozen=True)
class OlfLoadingParameters:
    """The OpenLoadFlow parameters the loading rules depend on, with their OpenLoadFlow
    names in the comments and the defaults of the reference pypowsybl."""

    #: ``useReactiveLimits``: off, no unit is discarded for its reactive range and no
    #: target Q is clamped
    reactive_limits: bool = True
    #: ``generatorsWithZeroMwTargetAreNotStarted``
    zero_mw_target_not_started: bool = True
    #: ``fictitiousGeneratorVoltageControlCheckMode`` = FORCED: a fictitious unit keeps its
    #: voltage control even when it looks not started
    fictitious_voltage_control_forced: bool = True
    #: ``disableInconsistentVoltageControls``
    disable_inconsistent_voltage_controls: bool = True
    #: ``forceTargetQInReactiveLimits`` (only with ``reactive_limits``)
    force_target_q_in_reactive_limits: bool = True
    #: ``extrapolateReactiveLimits``: a capability curve is extrapolated past its ends
    extrapolate_reactive_limits: bool = True
    #: ``minPlausibleTargetVoltage`` / ``maxPlausibleTargetVoltage``, pu
    min_plausible_target_v: float = _MIN_PLAUSIBLE_TARGET_V_PU
    max_plausible_target_v: float = _MAX_PLAUSIBLE_TARGET_V_PU
    #: ``minNominalVoltageTargetVoltageCheck``, kV: no plausibility check below it
    min_nominal_voltage_target_voltage_check: float = _MIN_NOMINAL_V_FOR_TARGET_V_CHECK_KV
    #: ``useActiveLimits``: a unit outside its active range does not take part in the slack
    use_active_limits: bool = True
    #: ``plausibleActivePowerLimit``, MW: a unit with a larger max P never takes part in it
    plausible_active_power_limit: float = _MAX_PLAUSIBLE_ACTIVE_POWER_MW
    #: ``disableVoltageControlOfGeneratorsOutsideActivePowerLimits`` (with ``use_active_limits``):
    #: a unit whose target P is outside its active range does not regulate
    disable_voltage_control_outside_active_limits: bool = False
    #: ``svcVoltageMonitoring``: an SVC whose ``standbyAutomaton`` says ``standby`` is a voltage
    #: monitor (idle until the VoltageMonitoring loop switches it on), and the automaton's
    #: ``b0`` is a fixed susceptance of every SVC carrying it
    svc_voltage_monitoring: bool = True
    #: ``voltagePerReactivePowerControl``: off, an SVC's slope is ignored
    voltage_per_reactive_power_control: bool = False


GEN_ATTRIBUTES = ["target_p", "min_p", "target_q", "target_v", "voltage_regulator_on",
                  "regulated_element_id", "regulated_bus_id", "bus_id", "connected",
                  "condenser", "fictitious", "reactive_limits_kind", "min_q", "max_q",
                  "min_q_at_target_p", "max_q_at_target_p"]


def generators(network):
    """The generator frame every rule of this module reads."""
    return network.get_generators(attributes=GEN_ATTRIBUTES)


# ---------------------------------------------------------------------------------------
# reactive capability
# ---------------------------------------------------------------------------------------

def curve_limits_at(network, ids, power_mw, min_q, max_q):
    """Reactive limits (generator convention, MVar) of the units ``ids`` at the active
    powers ``power_mw`` (generator convention), the capability curve being extrapolated
    past its ends by its end segment, as OpenLoadFlow does with ``extrapolateReactiveLimits``
    (limits extrapolated across each other are both their mean).

    ``min_q`` / ``max_q`` are the limits pypowsybl reports at those powers -- clamped to the
    end points of the curve -- and are returned as is for a unit inside the P range of its
    curve, without a curve, or with a curve of fewer than two points. Returns two numpy
    arrays aligned on ``ids``.
    """
    ids = pd.Index(ids)
    qmin = np.asarray(min_q, dtype=float).copy()
    qmax = np.asarray(max_q, dtype=float).copy()
    power = np.asarray(power_mw, dtype=float)
    if not len(ids):
        return qmin, qmax
    try:
        pts = network.get_reactive_capability_curve_points()
    except Exception:  # noqa: BLE001 - no curve support in this pypowsybl
        return qmin, qmax
    pts = pts.loc[pts.index.get_level_values(0).isin(ids)]
    if not len(pts):
        return qmin, qmax
    # one vectorised pass over every curve (a per-curve loop dominates the whole bake on
    # a large grid): sort the points by (element, p), then read each curve's first two
    # and last two points by position
    pt_ids = pts.index.get_level_values(0).to_numpy()
    p_pt = pts["p"].to_numpy(float)
    codes = pd.factorize(pt_ids)[0]
    order = np.lexsort((p_pt, codes))  # stable: points of equal p keep their curve order
    pt_ids, codes, p_pt = pt_ids[order], codes[order], p_pt[order]
    qmin_pt = pts["min_q"].to_numpy(float)[order]
    qmax_pt = pts["max_q"].to_numpy(float)[order]
    start = np.flatnonzero(np.r_[True, codes[1:] != codes[:-1]])
    count = np.diff(np.r_[start, len(pt_ids)])
    start, count = start[count >= 2], count[count >= 2]
    if not len(start):
        return qmin, qmax
    el_ids = pt_ids[start]
    rows = ids.get_indexer(el_ids)
    p_el = power[rows]
    below = p_el < p_pt[start]
    above = p_el > p_pt[start + count - 1]
    # the end segment on the side the unit lies (NaN power is neither below nor above)
    i1 = np.where(below, start, start + count - 2)
    i2 = i1 + 1
    p1, p2 = p_pt[i1], p_pt[i2]
    keep = (below | above) & (p2 != p1)
    if not keep.any():
        return qmin, qmax
    i1, i2, p1, p2, p_el, rows = i1[keep], i2[keep], p1[keep], p2[keep], p_el[keep], rows[keep]
    lo = qmin_pt[i1] + (qmin_pt[i2] - qmin_pt[i1]) * (p_el - p1) / (p2 - p1)
    hi = qmax_pt[i1] + (qmax_pt[i2] - qmax_pt[i1]) * (p_el - p1) / (p2 - p1)
    # extrapolated limits that cross are both their mean (powsybl-core's
    # ReactiveCapabilityCurveUtil, which OpenLoadFlow reads them through)
    crossed = lo > hi
    mean = (lo + hi) / 2.
    lo = np.where(crossed, mean, lo)
    hi = np.where(crossed, mean, hi)
    qmin[rows] = lo
    qmax[rows] = hi
    return qmin, qmax


def curve_limits_at_p(network, ids, power_mw, min_q, max_q):
    """Reactive limits (MVar) of the units ``ids`` at the active powers ``power_mw``, read on
    their capability curve -- interpolated inside it, extrapolated past its ends by its end
    segment, crossed extrapolated limits being their mean (powsybl-core's
    ReactiveCapabilityCurveUtil, as OpenLoadFlow reads them). ``min_q`` / ``max_q`` are kept
    for a unit without a curve of at least two points. Two numpy arrays aligned on ``ids``."""
    ids = pd.Index(ids)
    qmin = np.asarray(min_q, dtype=float).copy()
    qmax = np.asarray(max_q, dtype=float).copy()
    power = np.asarray(power_mw, dtype=float)
    if not len(ids):
        return qmin, qmax
    try:
        pts = network.get_reactive_capability_curve_points()
    except Exception:  # noqa: BLE001 - no curve support in this pypowsybl
        return qmin, qmax
    pts = pts.loc[pts.index.get_level_values(0).isin(ids)]
    for el_id, curve in pts.groupby(level=0):
        curve = curve.sort_values("p")
        if len(curve) < 2:
            continue
        k = ids.get_loc(el_id)
        p = power[k]
        if not np.isfinite(p):
            continue
        p_pt = curve["p"].to_numpy(float)
        # the segment p lies on, or the end one it extrapolates
        i = int(np.clip(np.searchsorted(p_pt, p) - 1, 0, len(p_pt) - 2))
        p1, p2 = p_pt[i], p_pt[i + 1]
        if p2 == p1:
            continue
        t = (p - p1) / (p2 - p1)
        lo_pt = curve["min_q"].to_numpy(float)
        hi_pt = curve["max_q"].to_numpy(float)
        lo = lo_pt[i] + (lo_pt[i + 1] - lo_pt[i]) * t
        hi = hi_pt[i] + (hi_pt[i + 1] - hi_pt[i]) * t
        if lo > hi:
            lo = hi = (lo + hi) / 2.
        qmin[k], qmax[k] = lo, hi
    return qmin, qmax


def vsc_station_target_p(network):
    """The active power (MW, generator convention) each VSC station injects, as OpenLoadFlow
    loads it (powsybl-core's HvdcUtils.getConverterStationTargetP): the rectifier draws the
    hvdc line's set-point, the inverter gives it back after the losses of both converters and
    of the line (R . pDc^2 / V^2). 0 when the other station is disconnected. A
    ``pandas.Series`` indexed by station id."""
    vsc = network.get_vsc_converter_stations(attributes=["loss_factor", "connected"])
    res = pd.Series(0., index=vsc.index)
    hvdc = network.get_hvdc_lines(attributes=["converter_station1_id", "converter_station2_id",
                                              "converters_mode", "target_p", "r", "nominal_v"])
    for _, line in hvdc.iterrows():
        s1, s2 = line["converter_station1_id"], line["converter_station2_id"]
        if s1 not in vsc.index or s2 not in vsc.index:
            continue
        if not (bool(vsc.loc[s1, "connected"]) and bool(vsc.loc[s2, "connected"])):
            continue
        rect, inv = (s1, s2) if str(line["converters_mode"]) == "SIDE_1_RECTIFIER_SIDE_2_INVERTER" else (s2, s1)
        setpoint = abs(float(line["target_p"]))
        p_dc1 = setpoint * (1. - float(vsc.loc[rect, "loss_factor"]) / 100.)
        v = float(line["nominal_v"])
        p_dc2 = p_dc1 - (float(line["r"]) * p_dc1 * p_dc1 / (v * v) if v > 0. else 0.)
        res[rect] = -setpoint
        res[inv] = p_dc2 * (1. - float(vsc.loc[inv, "loss_factor"]) / 100.)
    return res


def generator_limits_at_target_p(network, gen, params=OlfLoadingParameters()):
    """Reactive limits (MVar, generator convention) of each generator of ``gen`` at its
    target P, as OpenLoadFlow's ``getMinQ`` / ``getMaxQ`` evaluate them: the fixed box of a
    MIN_MAX unit, its capability curve at target P otherwise (extrapolated with
    ``extrapolate_reactive_limits``). Two ``pandas.Series`` indexed like ``gen``."""
    qmin = gen["min_q_at_target_p"].where(gen["min_q_at_target_p"].notna(), gen["min_q"])
    qmax = gen["max_q_at_target_p"].where(gen["max_q_at_target_p"].notna(), gen["max_q"])
    if params.extrapolate_reactive_limits:
        lo, hi = curve_limits_at(network, gen.index, gen["target_p"].to_numpy(float),
                                 qmin.to_numpy(float), qmax.to_numpy(float))
        qmin = pd.Series(lo, index=gen.index)
        qmax = pd.Series(hi, index=gen.index)
    return qmin, qmax


def generator_max_reactive_range(network, gen_index, box=None):
    """OpenLoadFlow's ``reactiveRangeCheckMode`` = MAX: the widest ``max_q - min_q`` across
    the whole active-power range of the capability curve of a CURVE-kind generator, simply
    ``max_q - min_q`` for a MIN_MAX (fixed box) one. A ``pandas.Series`` of ranges (MVar)
    indexed like ``gen_index``, NaN for an id not found. ``box`` (``reactive_limits_kind``,
    ``min_q``, ``max_q``) is the generators' by default; a battery's or a VSC station's frame
    gives theirs, their capability curves being read the same way."""
    rng = pd.Series(np.nan, index=gen_index)
    if not len(gen_index):
        return rng
    if box is None:
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


# ---------------------------------------------------------------------------------------
# voltage control of the generators
# ---------------------------------------------------------------------------------------

def generator_not_started(gen, params=OlfLoadingParameters()):
    """``AbstractLfGenerator.checkIfGeneratorStartedForVoltageControl``: a unit dispatched
    at zero whose minimum active power is not, unless it is a condenser or (in FORCED mode)
    a fictitious unit. OpenLoadFlow compares per-unit values with its epsilon there (see
    ``_NOT_STARTED_P_TOL_MW``)."""
    if not params.zero_mw_target_not_started:
        return pd.Series(False, index=gen.index)
    forced = gen["condenser"].fillna(False).astype(bool)
    if params.fictitious_voltage_control_forced:
        forced |= gen["fictitious"].fillna(False).astype(bool)
    return (gen["target_p"].abs() < _NOT_STARTED_P_TOL_MW) & (gen["min_p"] > _NOT_STARTED_P_TOL_MW) & ~forced


def generator_reactive_range_too_small(network, gen, params=OlfLoadingParameters()):
    """``checkIfReactiveRangesAreLargeEnoughForVoltageControl`` (MAX mode): a reactive range
    below ``PlausibleValues.MIN_REACTIVE_RANGE``. Only with reactive limits."""
    if not params.reactive_limits:
        return pd.Series(False, index=gen.index)
    return generator_max_reactive_range(network, gen.index) < _MIN_REACTIVE_RANGE_MVAR


def generator_regulated_bus(gen):
    """Bus-view id of the bus each generator regulates (pypowsybl's ``regulated_bus_id``,
    which resolves any regulating terminal), "" where it cannot be resolved: the regulated
    element is disconnected."""
    reg = gen["regulated_bus_id"].fillna("").astype(str)
    return reg


def generator_target_v_implausible(gen, reg_nominal_v, params=OlfLoadingParameters()):
    """``AbstractLfGenerator.checkTargetV``: the target voltage, in per unit of the
    REGULATED bus' nominal voltage, outside the plausible band -- only on a regulated bus
    above ``min_nominal_voltage_target_voltage_check``."""
    nom = np.asarray(reg_nominal_v, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        target_pu = gen["target_v"].to_numpy(float) / nom
    bad = ((nom > params.min_nominal_voltage_target_voltage_check)
           & ((target_pu < params.min_plausible_target_v) | (target_pu > params.max_plausible_target_v)))
    return pd.Series(bad, index=gen.index)


def generator_inconsistent_controls(gen, kept, reg_bus, target_v_pu, params=OlfLoadingParameters()):
    """``LfNetworkLoaderImpl.checkAndCreateVoltageControl`` with
    ``disableInconsistentVoltageControls``: on a bus whose voltage-controlling units (those
    of ``kept``) regulate different buses, or the same bus with targets more than
    ``TARGET_V_EPSILON`` apart from the first unit's, every one of them loses its voltage
    control. The first unit of a bus is the first in ``gen``'s order."""
    res = pd.Series(False, index=gen.index)
    if not params.disable_inconsistent_voltage_controls or not kept.any():
        return res
    frame = pd.DataFrame({"bus": gen["bus_id"], "reg": reg_bus, "tv": target_v_pu})[kept.to_numpy(bool)]
    first = frame.groupby("bus", sort=False).transform("first")
    several_regulated = frame.groupby("bus", sort=False)["reg"].transform("nunique") > 1
    # the target is compared only between units regulating the first unit's bus
    same_reg = frame["reg"] == first["reg"]
    far = same_reg & ((frame["tv"] - first["tv"]).abs() > _TARGET_V_EPSILON_PU)
    bad_bus = (several_regulated | far).groupby(frame["bus"], sort=False).transform("any")
    res.loc[frame.index[bad_bus.to_numpy(bool)]] = True
    return res


def generator_voltage_control(network, gen=None, params=OlfLoadingParameters(), bus_nominal_v=None):
    """Which voltage-regulating generators OpenLoadFlow lets regulate, and why not the others.

    Returns a ``pandas.DataFrame`` indexed like ``gen`` (by default every generator) with one
    boolean column per reason a regulating unit is discarded, in the order OpenLoadFlow
    checks them, and ``discarded`` (any of them). A unit not regulating in the first place
    reads False everywhere. ``bus_nominal_v`` maps a bus-view id to its nominal voltage (kV),
    read from the network when not given. The generators' rows of
    :func:`voltage_controllers`, so a generator sharing its bus with another kind of
    regulating unit is judged together with it, as OpenLoadFlow does.
    """
    gen = generators(network) if gen is None else gen
    res = voltage_controllers(network, params, gen=gen, bus_nominal_v=bus_nominal_v)
    res = res[res["kind"] == "generator"].drop(columns=["kind", "monitor"])
    return res.reindex(gen.index).fillna(False).astype(bool)


# ---------------------------------------------------------------------------------------
# voltage control of every generator-like unit
# ---------------------------------------------------------------------------------------

#: the generator-like units OpenLoadFlow sets up the same way
#: (``AbstractLfGenerator.setVoltageControl``). OpenLoadFlow orders the units of a bus as
#: IIDM visits its terminals, which the dataframes do not give: they are taken by kind, in
#: this order. That order only decides which unit's target the others of a bus are compared
#: with (``disableInconsistentVoltageControls``).
UNIT_KINDS = ("generator", "battery", "svc", "vsc")

#: the reasons a regulating unit is discarded, in the order OpenLoadFlow checks them
DISCARD_REASONS = ("reactive_range_too_small", "not_started", "outside_active_limits",
                   "regulated_bus_unresolved", "target_v_implausible", "monitor_with_regulator",
                   "inconsistent_controls")


def _extension(network, name):
    try:
        ext = network.get_extensions(name)
    except Exception:  # noqa: BLE001 - extension unsupported / absent on this pypowsybl
        return None
    return ext if ext is not None and len(ext) else None


def _box(frame):
    """The reactive box frame :func:`generator_max_reactive_range` reads (MIN_MAX when the
    kind is not reported)."""
    kind = frame["reactive_limits_kind"] if "reactive_limits_kind" in frame else "MIN_MAX"
    return pd.DataFrame({"reactive_limits_kind": kind, "min_q": frame["min_q"],
                         "max_q": frame["max_q"]}, index=frame.index)


def _units(network, gen):
    """One frame over every generator-like unit of ``network`` (``gen``: the generators'),
    with what the voltage-control rules read: ``kind``, ``bus_id``, ``regulating``,
    ``regulated_bus_id`` ("" if unresolved), ``target_v`` (kV), ``target_p`` / ``min_p`` /
    ``max_p`` (MW), ``forced`` (exempt from the not-started rule), ``range_q`` (MVar, the
    reactive range of ``reactiveRangeCheckMode`` = MAX, NaN where unknown) and ``standby``
    (an SVC whose automaton is in standby)."""
    parts = []
    # generators
    g = gen
    max_p = network.get_generators(attributes=["max_p"])["max_p"].reindex(g.index)
    parts.append(pd.DataFrame({
        "kind": "generator", "bus_id": g["bus_id"],
        "regulating": g["voltage_regulator_on"].fillna(False).astype(bool) & g["connected"].astype(bool),
        "regulated_bus_id": generator_regulated_bus(g), "target_v": g["target_v"],
        "target_p": g["target_p"], "min_p": g["min_p"], "max_p": max_p,
        "forced": g["condenser"].fillna(False).astype(bool), "fictitious": g["fictitious"].fillna(False).astype(bool),
        "range_q": generator_max_reactive_range(network, g.index), "standby": False}, index=g.index))
    # batteries: regulating through the voltageRegulation extension, their own bus
    b = network.get_batteries(all_attributes=True)
    if len(b):
        ext = _extension(network, "voltageRegulation")
        on = pd.Series(False, index=b.index)
        tv = pd.Series(np.nan, index=b.index)
        if ext is not None and "voltage_regulator_on" in ext.columns:
            ext = ext.reindex(b.index)
            on = ext["voltage_regulator_on"].fillna(False).astype(bool)
            tv = ext["target_v"].astype(float)
        parts.append(pd.DataFrame({
            "kind": "battery", "bus_id": b["bus_id"], "regulating": on & b["connected"].astype(bool),
            "regulated_bus_id": b["bus_id"].fillna("").astype(str), "target_v": tv,
            "target_p": b["target_p"], "min_p": b["min_p"], "max_p": b["max_p"],
            "forced": False, "fictitious": False,
            "range_q": generator_max_reactive_range(network, b.index, _box(b)), "standby": False}, index=b.index))
    # SVCs: the range of B over [b_min, b_max] at the voltage the snapshot stores (NaN: none)
    sv = network.get_static_var_compensators(all_attributes=True)
    if len(sv):
        regulating = (sv["regulation_mode"].astype(str) == "VOLTAGE") & sv["connected"].astype(bool)
        if "regulating" in sv.columns:
            regulating &= sv["regulating"].fillna(False).astype(bool)
        v_kv = sv["bus_id"].map(network.get_buses(attributes=["v_mag"])["v_mag"]).astype(float)
        standby = pd.Series(False, index=sv.index)
        sa = _extension(network, "standbyAutomaton")
        if sa is not None and "standby" in sa.columns:
            standby = sa["standby"].reindex(sv.index).fillna(False).astype(bool)
        parts.append(pd.DataFrame({
            "kind": "svc", "bus_id": sv["bus_id"], "regulating": regulating,
            "regulated_bus_id": sv["regulated_bus_id"].fillna("").astype(str), "target_v": sv["target_v"],
            "target_p": 0., "min_p": -np.inf, "max_p": np.inf, "forced": False, "fictitious": False,
            "range_q": (sv["b_max"] - sv["b_min"]) * v_kv ** 2, "standby": standby}, index=sv.index))
    # VSC stations: local control only; their active range is the hvdc line's
    vs = network.get_vsc_converter_stations(all_attributes=True)
    if len(vs):
        hvdc = network.get_hvdc_lines(attributes=["converter_station1_id", "converter_station2_id",
                                                  "max_p", "target_p"])
        line_max_p = pd.concat([hvdc.set_index("converter_station1_id")["max_p"],
                                hvdc.set_index("converter_station2_id")["max_p"]])
        line_target_p = pd.concat([hvdc.set_index("converter_station1_id")["target_p"],
                                   hvdc.set_index("converter_station2_id")["target_p"]])
        mp = line_max_p.reindex(vs.index).astype(float).fillna(np.inf)
        parts.append(pd.DataFrame({
            "kind": "vsc", "bus_id": vs["bus_id"],
            "regulating": vs["voltage_regulator_on"].fillna(False).astype(bool) & vs["connected"].astype(bool),
            "regulated_bus_id": vs["bus_id"].fillna("").astype(str), "target_v": vs["target_v"],
            # the size of the transfer: what the active range is checked against
            "target_p": line_target_p.reindex(vs.index).astype(float).abs(), "min_p": -mp, "max_p": mp,
            "forced": False, "fictitious": False,
            "range_q": generator_max_reactive_range(network, vs.index, _box(vs)), "standby": False}, index=vs.index))
    units = pd.concat(parts)
    units["regulated_bus_id"] = units["regulated_bus_id"].fillna("").astype(str)
    return units


def voltage_controllers(network, params=OlfLoadingParameters(), gen=None, bus_nominal_v=None):
    """Which voltage-regulating units -- generators, batteries, SVCs, VSC stations --
    OpenLoadFlow lets regulate, which SVCs it makes voltage monitors, and why.

    ``AbstractLfGenerator.setVoltageControl`` for each unit (its reactive range, a unit not
    started, a target P outside the active range, a regulated bus out of reach, an
    implausible target), then ``LfNetworkLoaderImpl.createVoltageControls`` for each bus:
    with ``svc_voltage_monitoring``, an SVC whose automaton is in standby is a monitor; two
    monitors on a bus both regulate instead, and a monitor sharing its bus with a regulating
    unit is switched off (``monitor_with_regulator``); then the units left on a bus,
    monitors included, are checked for consistency (``disableInconsistentVoltageControls``).

    Returns a ``pandas.DataFrame`` indexed by unit id with ``kind`` (see ``UNIT_KINDS``), a
    boolean column per reason (``DISCARD_REASONS``), ``discarded`` (any of them) and
    ``monitor``. A unit not regulating in the first place reads False everywhere.
    ``gen`` is the generators' frame (by default :func:`generators`), ``bus_nominal_v`` maps a
    bus-view id to its nominal voltage (kV), read from the network when not given.
    """
    gen = generators(network) if gen is None else gen
    units = _units(network, gen)
    buses = network.get_buses(attributes=["voltage_level_id", "synchronous_component"])
    if bus_nominal_v is None:
        vls = network.get_voltage_levels(attributes=["nominal_v"])
        bus_nominal_v = buses["voltage_level_id"].map(vls["nominal_v"])
    reg = units["regulating"].to_numpy(bool)
    reg_bus = units["regulated_bus_id"]
    sync = buses["synchronous_component"]

    res = pd.DataFrame(index=units.index)
    res["kind"] = units["kind"]
    # AbstractLfGenerator.setVoltageControl, in its order
    if params.reactive_limits:
        res["reactive_range_too_small"] = reg & (units["range_q"] < _MIN_REACTIVE_RANGE_MVAR).to_numpy(bool)
    else:
        res["reactive_range_too_small"] = False
    forced = units["forced"].to_numpy(bool)
    if params.fictitious_voltage_control_forced:
        forced = forced | units["fictitious"].to_numpy(bool)
    not_started = ((units["target_p"].abs() < _NOT_STARTED_P_TOL_MW) & (units["min_p"] > _NOT_STARTED_P_TOL_MW)).to_numpy(bool)
    res["not_started"] = reg & params.zero_mw_target_not_started & not_started & ~forced
    outside = ((units["target_p"] < units["min_p"]) | (units["target_p"] > units["max_p"])).to_numpy(bool)
    res["outside_active_limits"] = (reg & outside & params.use_active_limits
                                    & params.disable_voltage_control_outside_active_limits)
    other_component = (reg_bus.map(sync) != units["bus_id"].map(sync)).to_numpy(bool)
    res["regulated_bus_unresolved"] = reg & ((reg_bus == "").to_numpy(bool) | other_component)
    reg_nominal_v = reg_bus.map(bus_nominal_v).to_numpy(float)
    implausible = generator_target_v_implausible(units, reg_nominal_v, params).to_numpy(bool)
    res["target_v_implausible"] = reg & implausible
    kept = reg & ~res[list(DISCARD_REASONS[:5])].any(axis=1).to_numpy(bool)

    # LfNetworkLoaderImpl.createVoltageControls: the monitors of each bus
    monitor = kept & params.svc_voltage_monitoring & units["standby"].to_numpy(bool)
    bus = units["bus_id"]
    n_monitors = pd.Series(monitor, index=units.index).groupby(bus).transform("sum").to_numpy()
    n_regulators = pd.Series(kept & ~monitor, index=units.index).groupby(bus).transform("sum").to_numpy()
    monitor = monitor & (n_monitors == 1)          # several on a bus: they all regulate
    res["monitor_with_regulator"] = monitor & (n_regulators > 0)
    monitor = monitor & (n_regulators == 0)
    kept = kept & ~res["monitor_with_regulator"].to_numpy(bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        target_v_pu = pd.Series(units["target_v"].to_numpy(float) / reg_nominal_v, index=units.index)
    res["inconsistent_controls"] = generator_inconsistent_controls(
        units, pd.Series(kept, index=units.index), reg_bus, target_v_pu, params).to_numpy(bool)
    res["discarded"] = res[list(DISCARD_REASONS)].any(axis=1)
    res["monitor"] = monitor & ~res["inconsistent_controls"].to_numpy(bool)
    return res


def generator_target_q(network, gen=None, params=OlfLoadingParameters()):
    """What a non-regulating generator injects (MVar, generator convention): its target Q,
    0 when unset, clamped into its reactive limits at its target P with
    ``forceTargetQInReactiveLimits`` (``LfGeneratorImpl.getTargetQ``). A
    ``pandas.Series`` indexed like ``gen``."""
    gen = generators(network) if gen is None else gen
    target_q = gen["target_q"].fillna(0.)
    if not (params.reactive_limits and params.force_target_q_in_reactive_limits):
        return target_q
    qmin, qmax = generator_limits_at_target_p(network, gen, params)
    # a NaN limit clamps nothing
    return target_q.clip(lower=qmin.fillna(-np.inf), upper=qmax.fillna(np.inf))


# ---------------------------------------------------------------------------------------
# participation in the distributed slack (PROPORTIONAL_TO_GENERATION_P_MAX)
# ---------------------------------------------------------------------------------------

def target_p_range(min_p, max_p, min_target_p, max_target_p):
    """``(min_target_p, max_target_p)`` with OpenLoadFlow's defaults (``min_p`` / ``max_p``)
    where the ``activePowerControl`` extension does not set them (NaN)."""
    min_p = np.asarray(min_p, dtype=float)
    max_p = np.asarray(max_p, dtype=float)
    min_target_p = np.asarray(min_target_p, dtype=float)
    max_target_p = np.asarray(max_target_p, dtype=float)
    return (np.where(np.isfinite(min_target_p), min_target_p, min_p),
            np.where(np.isfinite(max_target_p), max_target_p, max_p))


def participation_weight(target_p, min_p, max_p, participate, droop, min_target_p, max_target_p,
                         params=OlfLoadingParameters()):
    """OpenLoadFlow's distributed-slack key of each unit (``PROPORTIONAL_TO_GENERATION_P_MAX``),
    0 for a unit that does not take part, generators and batteries alike
    (``ActivePowerControlHelper``, ``AbstractLfGenerator.checkActivePowerControl`` and
    ``GenerationActivePowerDistributionStep.isParticipating``):

    * the ``activePowerControl`` extension's ``participate`` (True without the extension);
    * a target of zero (below ``_ZERO_P_TOL``, in MW here) does not take part;
    * nor a max P above ``plausible_active_power_limit``;
    * with ``use_active_limits``, nor a target outside ``[min_target_p, max_target_p]`` or a
      range narrower than ``_ZERO_P_TOL``;
    * the key is ``max_p / droop``, the droop being the extension's or 4 when it sets none
      (NaN); a droop or a key of 0 does not take part.

    Every input is an array aligned on the units, in MW and in the generator convention; a
    NaN ``min_target_p`` / ``max_target_p`` means ``min_p`` / ``max_p``. The connectivity
    and the synchronous component are the caller's to check."""
    target_p = np.asarray(target_p, dtype=float)
    min_p = np.asarray(min_p, dtype=float)
    max_p = np.asarray(max_p, dtype=float)
    droop = np.asarray(droop, dtype=float)
    min_tp, max_tp = target_p_range(min_p, max_p, min_target_p, max_target_p)
    droop_used = np.where(np.isfinite(droop), droop, _OLF_DEFAULT_DROOP)
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = np.where(droop_used != 0., max_p / droop_used, 0.)
        ok = np.asarray(participate, dtype=bool) & (max_p <= params.plausible_active_power_limit)
        if params.zero_mw_target_not_started:
            ok &= np.abs(target_p) >= _ZERO_P_TOL
        if params.use_active_limits:
            ok &= (target_p <= max_tp) & (target_p >= min_tp) & ((max_tp - min_tp) >= _ZERO_P_TOL)
        ok &= (droop_used != 0.) & np.isfinite(weight) & (weight != 0.)
    return np.where(ok, weight, 0.)


def generator_active_power_control(network, gen_index):
    """The ``activePowerControl`` extension of the generators ``gen_index``: a
    ``pandas.DataFrame`` with ``participate`` (True for a unit without the extension),
    ``droop``, ``min_target_p`` and ``max_target_p`` (NaN where unset)."""
    res = pd.DataFrame({"participate": True, "droop": np.nan, "min_target_p": np.nan,
                        "max_target_p": np.nan}, index=gen_index)
    try:
        apc = network.get_extensions("activePowerControl")
    except Exception:  # noqa: BLE001 - extension unknown to this pypowsybl
        return res
    if apc is None or not len(apc):
        return res
    apc = apc.reindex(gen_index)
    listed = apc.index.isin(apc.dropna(how="all").index)
    if "participate" in apc.columns:
        res.loc[listed, "participate"] = apc.loc[listed, "participate"].fillna(True).astype(bool)
    for col in ("droop", "min_target_p", "max_target_p"):
        if col in apc.columns:
            res[col] = apc[col].astype(float)
    return res


def generator_participation_weight(network, gen, params=OlfLoadingParameters()):
    """:func:`participation_weight` of each generator of ``gen`` (which carries
    ``target_p``, ``min_p`` and ``max_p``), as a ``pandas.Series`` indexed like ``gen``."""
    apc = generator_active_power_control(network, gen.index)
    weight = participation_weight(gen["target_p"].to_numpy(float), gen["min_p"].to_numpy(float),
                                  gen["max_p"].to_numpy(float), apc["participate"].to_numpy(bool),
                                  apc["droop"].to_numpy(float), apc["min_target_p"].to_numpy(float),
                                  apc["max_target_p"].to_numpy(float), params)
    return pd.Series(weight, index=gen.index)


# ---------------------------------------------------------------------------------------
# branches
# ---------------------------------------------------------------------------------------

def branch_on_same_bus(branches):
    """``LfNetworkLoaderImpl.addBranch``: a branch whose two ends are connected to the same
    bus is discarded. ``branches`` is a frame of pypowsybl's ``get_lines`` /
    ``get_2_windings_transformers`` (``bus1_id``, ``bus2_id``, ``connected1``,
    ``connected2``), read on the network's own buses -- before lightsim2grid fuses any of
    them: OpenLoadFlow keeps a branch whose ends are only joined by a zero-impedance one.
    A numpy boolean array aligned on ``branches``."""
    bus1 = branches["bus1_id"].astype(object)
    bus2 = branches["bus2_id"].astype(object)
    both_connected = (branches["connected1"].fillna(False).astype(bool) &
                      branches["connected2"].fillna(False).astype(bool))
    return (both_connected & bus1.notna() & (bus1 != "") & (bus1 == bus2)).to_numpy(bool)
