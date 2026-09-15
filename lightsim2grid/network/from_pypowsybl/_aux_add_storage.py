# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings

import numpy as np

from ._aux_common import _aux_get_bus


def _aux_battery_voltage_regulation(net, df_batt, batt_bus, batt_disco, voltage_levels):
    """The voltage side of the batteries: ``(voltage_regulator_on, target_vm_pu)``
    read from the IIDM ``voltageRegulation`` extension (OpenLoadFlow runs such a
    battery as a PV bus), aligned with ``df_batt``. Only LOCAL regulation is
    supported by lightsim2grid's storage units: a battery whose extension points
    at another element is left as a PQ injection, with a warning."""
    n = len(df_batt)
    vreg = np.zeros(n, dtype=bool)
    target_vm = np.ones(n, dtype=float)
    try:
        ext = net.get_extensions("voltageRegulation")
    except Exception:
        return vreg, target_vm
    if ext is None or not len(ext) or "voltage_regulator_on" not in ext.columns:
        return vreg, target_vm
    ext = ext.reindex(df_batt.index)
    on = ext["voltage_regulator_on"].fillna(False).to_numpy(bool) & ~np.asarray(batt_disco, dtype=bool)
    if not on.any():
        return vreg, target_vm
    if "regulated_element_id" in ext.columns:
        reg = ext["regulated_element_id"].fillna("").astype(str).to_numpy()
        remote = on & (reg != "") & (reg != df_batt.index.to_numpy().astype(str))
        if remote.any():
            warnings.warn("Batteries regulating the voltage of another element are not supported "
                          "(lightsim2grid storage units only regulate their own bus): "
                          f"{list(df_batt.index[remote])} are converted as PQ injections.")
            on &= ~remote
    nominal_v = voltage_levels.loc[df_batt["voltage_level_id"].values, "nominal_v"].to_numpy(float)
    tv = ext["target_v"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        pu = tv / nominal_v
    ok = on & np.isfinite(pu) & (pu > 0.)
    vreg[ok] = True
    target_vm[ok] = pu[ok]
    return vreg, target_vm


def _aux_battery_q_limits(df_batt):
    """``(min_q, max_q)`` of the batteries in MVAr (generator convention), the
    capability curve at the target P when the battery has one, the fixed box
    otherwise; NaN / absurd values become the float32 "unbounded" sentinels the
    generators use (see `_aux_add_generators.py`)."""
    def col(name, fallback):
        if name in df_batt.columns:
            s = df_batt[name]
            return s.where(s.notna(), df_batt[fallback] if fallback in df_batt.columns else np.nan).to_numpy(float)
        return df_batt[fallback].to_numpy(float) if fallback in df_batt.columns else np.full(len(df_batt), np.nan)
    min_q = col("min_q_at_target_p", "min_q")
    max_q = col("max_q_at_target_p", "max_q")
    min_float_value = np.finfo(np.float32).min * 1e-4 + 1.
    max_float_value = np.finfo(np.float32).max * 1e-4 + 1.
    swapped = np.isfinite(min_q) & np.isfinite(max_q) & (min_q > max_q)
    if swapped.any():
        min_q[swapped], max_q[swapped] = max_q[swapped], min_q[swapped].copy()
    min_q = np.where(np.isfinite(min_q) & (min_q >= min_float_value), min_q, min_float_value)
    max_q = np.where(np.isfinite(max_q) & (np.abs(max_q) <= max_float_value), max_q, max_float_value)
    return min_q.astype(np.float32).astype(float), max_q.astype(np.float32).astype(float)


def _aux_add_storage(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every storage unit (IIDM battery) of ``net`` to ``model``. IIDM gives
    the battery setpoints in the *generator* convention (positive target_p =
    power produced / injected) while lightsim2grid stores storage as PQ in the
    *load* convention (positive = power drawn from the grid, *ie* charging), same
    as pandapower and grid2op. We negate to convert, and sanitize NaN (IIDM
    allows an unset target_q). A battery whose IIDM ``voltageRegulation``
    extension is on regulates the voltage of its own bus (a PV bus, as
    OpenLoadFlow runs it; see :func:`_aux_battery_voltage_regulation`). Returns
    ``(df_batt, batt_sub)``, used by the final substation-id bookkeeping and
    ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_batt = net.get_batteries().sort_index()
    else:
        df_batt = net.get_batteries()
    batt_bus, batt_disco, batt_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "storage", df_batt)
    batt_p = df_batt["target_p"].values.astype(float)
    batt_q = df_batt["target_q"].values.astype(float)
    batt_p = np.where(np.isfinite(batt_p), batt_p, 0.)
    batt_q = np.where(np.isfinite(batt_q), batt_q, 0.)
    vreg, target_vm = _aux_battery_voltage_regulation(net, df_batt, batt_bus, batt_disco, voltage_levels)
    min_q, max_q = _aux_battery_q_limits(df_batt)
    model.init_storages_full(-batt_p,  # IIDM generator convention -> lightsim2grid load convention
                             -batt_q,
                             [bool(el) for el in vreg],
                             target_vm,
                             min_q,
                             max_q,
                             batt_bus
                             )
    for batt_id, disco in enumerate(batt_disco):
        if disco:
           model.deactivate_storage(batt_id)
    model.set_storage_names(df_batt.index)

    return df_batt, batt_sub
