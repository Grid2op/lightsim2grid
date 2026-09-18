# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The distributed-slack participation of the batteries, as OpenLoadFlow decides it.

OpenLoadFlow distributes the slack on batteries with the same rule as on generators
(``LfBatteryImpl`` runs ``checkActivePowerControl`` on its target P), under the loop-free
``PROPORTIONAL_TO_GENERATION_P_MAX`` balance type:

* the ``activePowerControl`` extension gives ``participate`` (default true), ``droop``
  (default, and for a NaN, :data:`_OLF_DEFAULT_DROOP`) and ``min_target_p`` /
  ``max_target_p`` (default ``min_p`` / ``max_p``); ``participation_factor`` is not used
  by this balance type;
* a battery is excluded when ``|target_p| < 1e-4`` MW (generatorsWithZeroMwTargetAreNotStarted),
  ``max_p > 10000`` MW (plausibleActivePowerLimit), ``target_p`` lies outside
  ``[min_target_p, max_target_p]``, that range is degenerate, or ``droop == 0``; a
  *charging* battery (``target_p < 0``) does participate;
* its key is ``max_p / droop`` (``max_p`` itself, not ``max_target_p``).

Measured against OpenLoadFlow on a real RTE snapshot, battery by battery (shares x12.4,
x1, x0.144 and 0 for the ``droop = 0`` / ``participate = false`` ones).

pypowsybl up to 1.16.1 does not expose that extension on a battery
(``get_extensions("activePowerControl")`` only lists generators; fixed by pypowsybl PR
#1276, not released yet), so :func:`battery_active_power_control` can fall back to reading
it off an XIIDM export of the network.
"""

import io
import xml.etree.ElementTree as ET

import numpy as np
from packaging import version

from ._aux_common import PYPOWSYBL_VER
# OpenLoadFlow's `POWER_EPSILON_SI` and default `plausibleActivePowerLimit`, the ones the
# bake uses for the generators (_olf_bake imports this module lazily, for that reason)
from ._olf_bake import _ZERO_P_TOL, _MAX_PLAUSIBLE_ACTIVE_POWER_MW

# OpenLoadFlow's hardcoded fallback droop (`AbstractLfGenerator.DEFAULT_DROOP`, "why not"),
# for a unit whose `activePowerControl` extension does not set its own
_OLF_DEFAULT_DROOP = 4.0

# the last pypowsybl release whose `activePowerControl` extension ignores batteries
_PYPOWSYBL_NO_BATTERY_APC = version.parse("1.16.1")

BATTERY_APC_SOURCES = ("auto", "extension", "default")


def olf_participation_weight(target_p, min_p, max_p, participate, droop, min_target_p, max_target_p):
    """OpenLoadFlow's distributed-slack key of each unit (``PROPORTIONAL_TO_GENERATION_P_MAX``),
    0 for a unit ``checkActivePowerControl`` excludes. Every input is an array aligned on the
    units, in MW and in the generator convention; a NaN ``droop`` means the default one and
    a NaN ``min_target_p`` / ``max_target_p`` means ``min_p`` / ``max_p``. The
    connectivity and main-component filters are the caller's."""
    target_p = np.asarray(target_p, dtype=float)
    min_p = np.asarray(min_p, dtype=float)
    max_p = np.asarray(max_p, dtype=float)
    droop = np.asarray(droop, dtype=float)
    min_tp, max_tp = olf_target_p_range(min_p, max_p, min_target_p, max_target_p)
    droop_used = np.where(np.isfinite(droop), droop, _OLF_DEFAULT_DROOP)
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = np.where(droop_used != 0., max_p / droop_used, 0.)
        ok = (np.asarray(participate, dtype=bool)
              & (np.abs(target_p) >= _ZERO_P_TOL)
              & (max_p <= _MAX_PLAUSIBLE_ACTIVE_POWER_MW)
              & (target_p <= max_tp) & (target_p >= min_tp)
              & ((max_tp - min_tp) >= _ZERO_P_TOL)
              & (droop_used != 0.)
              & np.isfinite(weight) & (weight > 0.))
    return np.where(ok, weight, 0.)


def olf_target_p_range(min_p, max_p, min_target_p, max_target_p):
    """``(min_target_p, max_target_p)`` with OpenLoadFlow's defaults (``min_p`` / ``max_p``)
    where the extension does not set them (NaN)."""
    min_p = np.asarray(min_p, dtype=float)
    max_p = np.asarray(max_p, dtype=float)
    min_target_p = np.asarray(min_target_p, dtype=float)
    max_target_p = np.asarray(max_target_p, dtype=float)
    return (np.where(np.isfinite(min_target_p), min_target_p, min_p),
            np.where(np.isfinite(max_target_p), max_target_p, max_p))


def _pypowsybl_exposes_battery_apc():
    return PYPOWSYBL_VER > _PYPOWSYBL_NO_BATTERY_APC


def _apc_from_extension(net, batt_ids):
    """The rows of the `activePowerControl` extension that are batteries, as a dict
    ``{battery_id: (participate, droop, min_target_p, max_target_p)}``."""
    try:
        apc = net.get_extensions("activePowerControl")
    except Exception:
        return {}
    if apc is None or not len(apc):
        return {}
    common = apc.index.intersection(batt_ids)
    res = {}
    for batt_id in common:
        row = apc.loc[batt_id]
        res[batt_id] = (bool(row["participate"]) if "participate" in apc.columns else True,
                        float(row["droop"]) if "droop" in apc.columns else np.nan,
                        float(row["min_target_p"]) if "min_target_p" in apc.columns else np.nan,
                        float(row["max_target_p"]) if "max_target_p" in apc.columns else np.nan)
    return res


def _apc_from_xiidm(net, batt_ids):
    """Same as :func:`_apc_from_extension`, read off an XIIDM export of ``net``: what a
    pypowsybl that does not expose the extension on batteries leaves as the only way.
    The export is the size of the network file (~200 MB for a 7k-bus grid)."""
    wanted = set(str(el) for el in batt_ids)
    res = {}
    xml_text = net.save_to_string("XIIDM")
    current_id = None
    for event, elem in ET.iterparse(io.StringIO(xml_text), events=("start", "end")):
        tag = elem.tag.rsplit("}", 1)[-1]
        if event == "start":
            if tag == "extension":
                current_id = elem.get("id")
            continue
        if tag == "activePowerControl" and current_id in wanted:
            def _float(name):
                val = elem.get(name)
                return float(val) if val is not None else np.nan
            res[current_id] = (elem.get("participate", "true").lower() == "true",
                               _float("droop"),
                               _float("minTargetP"),
                               _float("maxTargetP"))
        elif tag == "extension":
            current_id = None
        elem.clear()
    return res


def battery_active_power_control(net, df_batt, source="auto"):
    """``(participate, droop, min_target_p, max_target_p)`` of each battery of ``df_batt``
    (aligned on it), from its ``activePowerControl`` extension; OpenLoadFlow's defaults --
    participating, NaN droop / target range -- for a battery without one.

    ``source``:

    * ``"auto"`` (default): the extension as pypowsybl lists it when it knows batteries do
      carry one, else read off an XIIDM export of the network (pypowsybl <= 1.16.1);
    * ``"extension"``: only what pypowsybl lists, never exporting the network (with
      pypowsybl <= 1.16.1 every battery then gets the defaults);
    * ``"default"``: OpenLoadFlow's defaults for every battery.
    """
    if source not in BATTERY_APC_SOURCES:
        raise RuntimeError(f"Unknown source {source!r} for the batteries' active power control, "
                           f"use one of {BATTERY_APC_SOURCES}.")
    n = len(df_batt)
    participate = np.ones(n, dtype=bool)
    droop = np.full(n, np.nan)
    min_target_p = np.full(n, np.nan)
    max_target_p = np.full(n, np.nan)
    if n == 0 or source == "default":
        return participate, droop, min_target_p, max_target_p

    if source == "extension" or _pypowsybl_exposes_battery_apc():
        rows = _apc_from_extension(net, df_batt.index)
    else:
        rows = _apc_from_xiidm(net, df_batt.index)

    for pos, batt_id in enumerate(df_batt.index):
        row = rows.get(batt_id, rows.get(str(batt_id)))
        if row is None:
            continue
        participate[pos], droop[pos], min_target_p[pos], max_target_p[pos] = row
    return participate, droop, min_target_p, max_target_p
