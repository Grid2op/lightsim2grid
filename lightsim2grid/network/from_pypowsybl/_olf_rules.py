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
    past its ends by its end segment, as OpenLoadFlow does with ``extrapolateReactiveLimits``.

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
    qmin[rows] = qmin_pt[i1] + (qmin_pt[i2] - qmin_pt[i1]) * (p_el - p1) / (p2 - p1)
    qmax[rows] = qmax_pt[i1] + (qmax_pt[i2] - qmax_pt[i1]) * (p_el - p1) / (p2 - p1)
    return qmin, qmax


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


def generator_max_reactive_range(network, gen_index):
    """OpenLoadFlow's ``reactiveRangeCheckMode`` = MAX: the widest ``max_q - min_q`` across
    the whole active-power range of the capability curve of a CURVE-kind generator, simply
    ``max_q - min_q`` for a MIN_MAX (fixed box) one. A ``pandas.Series`` of ranges (MVar)
    indexed like ``gen_index``, NaN for an id not found."""
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
    read from the network when not given.
    """
    gen = generators(network) if gen is None else gen
    regulating = gen["voltage_regulator_on"].fillna(False).astype(bool) & gen["connected"].astype(bool)
    if bus_nominal_v is None:
        buses = network.get_buses(attributes=["voltage_level_id"])
        vls = network.get_voltage_levels(attributes=["nominal_v"])
        bus_nominal_v = buses["voltage_level_id"].map(vls["nominal_v"])
    reg_bus = generator_regulated_bus(gen)
    reg_nominal_v = reg_bus.map(bus_nominal_v).to_numpy(float)

    res = pd.DataFrame(index=gen.index)
    # AbstractLfGenerator.setVoltageControl, in its order: consistency (range, started),
    # then the regulated terminal, then the plausibility of the target
    res["reactive_range_too_small"] = regulating & generator_reactive_range_too_small(network, gen, params)
    res["not_started"] = regulating & generator_not_started(gen, params)
    res["regulated_bus_unresolved"] = regulating & (reg_bus == "")
    res["target_v_implausible"] = regulating & generator_target_v_implausible(gen, reg_nominal_v, params)
    kept = regulating & ~res.any(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        target_v_pu = pd.Series(gen["target_v"].to_numpy(float) / reg_nominal_v, index=gen.index)
    res["inconsistent_controls"] = generator_inconsistent_controls(gen, kept, reg_bus, target_v_pu, params)
    res["discarded"] = res.any(axis=1)
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
