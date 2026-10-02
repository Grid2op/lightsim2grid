# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings

import numpy as np
import pandas as pd

from ._aux_common import _aux_get_bus


def _hvdc_pmax_per_direction(net, hvdc_ids, max_p_mw):
    """Maximum active power of each hvdc line, from side 1 to side 2 and the other way
    (MW), as OpenLoadFlow sees them: the ``hvdcOperatorActivePowerRange`` extension when
    the line carries it (``opr_from_cs1_to_cs2`` / ``opr_from_cs2_to_cs1``), the line's
    own ``max_p`` in both directions otherwise. It is the limit OLF's
    ``AcHvdcAcEmulationLimits`` outer loop saturates an angle-droop line at."""
    pmax_1to2 = np.array(max_p_mw, dtype=float)
    pmax_2to1 = np.array(max_p_mw, dtype=float)
    try:
        df_opr = net.get_extensions("hvdcOperatorActivePowerRange")
    except Exception:
        # extension tables may be unavailable on (very) old pypowsybl versions
        df_opr = None
    if df_opr is None or not df_opr.shape[0]:
        return pmax_1to2, pmax_2to1
    df_opr = df_opr.reindex(pd.Index(hvdc_ids))
    for pmax, col in ((pmax_1to2, "opr_from_cs1_to_cs2"), (pmax_2to1, "opr_from_cs2_to_cs1")):
        opr = df_opr[col].to_numpy(dtype=float)
        has_opr = np.isfinite(opr)
        pmax[has_opr] = opr[has_opr]
    return pmax_1to2, pmax_2to1


def _aux_add_hvdc(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, can_be_pv=None,
                  ac_emulation_frozen=None, olf_vc=None):
    """Add every HVDC line of ``net`` (VSC / LCC converter stations, possibly
    carrying the angle-droop ("AC emulation") extension) to ``model``. Returns
    ``(df_dc, hvdc_sub_from_id, hvdc_sub_to_id)``, used by the final
    substation-id bookkeeping and ``return_sub_id`` in `initLSGrid.py`.

    The converter station ids ``can_be_pv`` holds (what `bake_outer_loops` returns for
    the VSC stations it froze at a reactive limit, or a boolean Series indexed by id)
    are flagged with ``LSGrid.set_hvdc_can_be_pv``, per line and side.

    The hvdc line ids ``ac_emulation_frozen`` holds (what `bake_outer_loops` froze at their
    AC-emulation limit) are flagged with ``LSGrid.set_hvdc_ac_emulation_frozen``, their droop
    parameters read off the extension although it is disabled.

    ``olf_vc`` (:func:`._olf_rules.voltage_controllers`, with ``olf_rules``) says which
    voltage-regulating stations OpenLoadFlow lets regulate: the others inject their target Q."""
    if sort_index:
        df_dc = net.get_hvdc_lines().sort_index()
    else:
        df_dc = net.get_hvdc_lines()
    # all_attributes: the capability curve at the target P (min_q_at_target_p...) is not a
    # default column
    df_vsc = net.get_vsc_converter_stations(all_attributes=True)
    df_lcc = net.get_lcc_converter_stations()
    # the vsc / lcc frames have different columns (target_v / power_factor...):
    # the concatenation puts NaN where an attribute does not exist for a type
    df_stations = pd.concat([df_vsc, df_lcc])
    nb_dc = df_dc.shape[0]
    _max_hvdc_mva = 1.0e7  # when pypowsybl exposes NaN limits

    df_station1 = df_stations.loc[df_dc["converter_station1_id"].values]
    df_station2 = df_stations.loc[df_dc["converter_station2_id"].values]
    hvdc_bus_from_id, hvdc_from_disco, hvdc_sub_from_id = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "hvdc (side 1)", df_station1)
    hvdc_bus_to_id, hvdc_to_disco, hvdc_sub_to_id = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "hvdc (side 2)", df_station2)

    def _aux_hvdc_station_data(df_side):
        # type: 0 = VSC, 1 = LCC (ConverterStationContainer convention)
        is_lcc = df_side.index.isin(df_lcc.index)
        types = np.where(is_lcc, 1, 0).astype(int)
        loss_factor = df_side["loss_factor"].values / 100.  # pypowsybl % -> fraction
        loss_factor = np.where(np.isfinite(loss_factor), loss_factor, 0.)
        vreg_on = df_side["voltage_regulator_on"].values.astype(bool) if nb_dc else np.zeros(0, dtype=bool)
        vreg_on = vreg_on & ~is_lcc  # lcc never regulates (NaN -> random bool otherwise)
        if olf_vc is not None:
            vreg_on = vreg_on & ~olf_vc["discarded"].reindex(df_side.index).fillna(False).to_numpy(bool)
        vl_kv = voltage_levels.loc[df_side["voltage_level_id"].values]["nominal_v"].values
        vset_pu = df_side["target_v"].values / vl_kv
        vset_pu = np.where(np.isfinite(vset_pu), vset_pu, 1.0)
        qset = df_side["target_q"].values
        qset = np.where(np.isfinite(qset), qset, 0.)
        # as for the generators (see `_aux_add_generators.py`): "min_q" / "max_q" are NaN
        # for a station whose reactive_limits_kind is CURVE, its limits are the curve at
        # its target P -- read alone, such a station was unlimited
        no_curve = pd.Series(np.nan, index=df_side.index)
        min_q = df_side.get("min_q_at_target_p", no_curve).fillna(df_side["min_q"]).to_numpy(float)
        max_q = df_side.get("max_q_at_target_p", no_curve).fillna(df_side["max_q"]).to_numpy(float)
        # malformed curve data can give min_q > max_q at the target P (as for the generators)
        swapped = np.isfinite(min_q) & np.isfinite(max_q) & (min_q > max_q)
        min_q[swapped], max_q[swapped] = max_q[swapped], min_q[swapped].copy()
        min_q = np.where(np.isfinite(min_q), min_q, -_max_hvdc_mva)
        max_q = np.where(np.isfinite(max_q), max_q, _max_hvdc_mva)
        power_factor = df_side["power_factor"].values
        power_factor = np.where(np.isfinite(power_factor), power_factor, 1.0)
        return types, loss_factor, vreg_on, vset_pu, qset, min_q, max_q, power_factor

    type1, lf1, vreg1, vset1, qset1, min_q1, max_q1, pf1 = _aux_hvdc_station_data(df_station1)
    type2, lf2, vreg2, vset2, qset2, min_q2, max_q2, pf2 = _aux_hvdc_station_data(df_station2)

    # 0 = side 1 rectifier, 1 = side 2 rectifier (HvdcLineContainer convention)
    converters_mode = np.where(df_dc["converters_mode"].values.astype(str) == "SIDE_1_RECTIFIER_SIDE_2_INVERTER", 0, 1).astype(int)
    p_setpoint_mw = df_dc["target_p"].values.astype(float).copy()
    if (~np.isfinite(p_setpoint_mw)).any():
        warnings.warn("Some non finite values are found for hvdc target_p, they have been replaced by 0.")
        p_setpoint_mw[~np.isfinite(p_setpoint_mw)] = 0.
    r_ohm = df_dc["r"].values.astype(float)
    nominal_v_kv = df_dc["nominal_v"].values.astype(float)
    max_p_mw = df_dc["max_p"].values.astype(float)
    max_p_mw = np.where(np.isfinite(max_p_mw), max_p_mw, _max_hvdc_mva)
    pmax_1to2_mw, pmax_2to1_mw = _hvdc_pmax_per_direction(net, df_dc.index, max_p_mw)

    # the angle-droop active power control ("AC emulation"), an IIDM extension
    droop_enabled = np.zeros(nb_dc, dtype=bool)
    droop_p0_mw = np.zeros(nb_dc)
    droop_mw_per_deg = np.zeros(nb_dc)
    try:
        df_droop = net.get_extensions("hvdcAngleDroopActivePowerControl")
    except Exception:
        # extension tables may be unavailable on (very) old pypowsybl versions
        df_droop = None
    frozen_ids = set(str(el) for el in ac_emulation_frozen) if ac_emulation_frozen is not None else set()
    unknown = frozen_ids.difference(df_dc.index)
    if unknown:
        raise ValueError(f"`hvdc_ac_emulation_frozen`: unknown hvdc line id(s) {sorted(unknown)[:10]}.")
    frozen = np.zeros(nb_dc, dtype=bool)
    if df_droop is not None and df_droop.shape[0]:
        for hvdc_pos, line_id in enumerate(df_dc.index):
            if line_id not in df_droop.index:
                continue
            is_frozen = line_id in frozen_ids
            if not bool(df_droop.loc[line_id, "enabled"]) and not is_frozen:
                continue
            # a line frozen at its limit keeps its droop parameters (droop off): the check of
            # its release reads them
            droop_enabled[hvdc_pos] = bool(df_droop.loc[line_id, "enabled"])
            frozen[hvdc_pos] = is_frozen and not droop_enabled[hvdc_pos]
            droop_p0_mw[hvdc_pos] = float(df_droop.loc[line_id, "p0"])
            droop_mw_per_deg[hvdc_pos] = float(df_droop.loc[line_id, "droop"])

    model.init_hvdc_lines(hvdc_bus_from_id.astype(np.int32),
                          hvdc_bus_to_id.astype(np.int32),
                          [int(el) for el in type1],
                          [int(el) for el in type2],
                          lf1, lf2,
                          [bool(el) for el in vreg1],
                          [bool(el) for el in vreg2],
                          vset1, vset2,
                          qset1, qset2,
                          min_q1, max_q1, min_q2, max_q2,
                          pf1, pf2,
                          [int(el) for el in converters_mode],
                          p_setpoint_mw,
                          r_ohm,
                          nominal_v_kv,
                          [bool(el) for el in droop_enabled],
                          droop_p0_mw,
                          droop_mw_per_deg,
                          pmax_1to2_mw,
                          pmax_2to1_mw,
                          )
    for hvdc_id, (is_or_disc, is_ex_disc, line_conn1, line_conn2) in enumerate(
            zip(hvdc_from_disco, hvdc_to_disco, df_dc["connected1"].values, df_dc["connected2"].values)):
        # a converter station with its own terminal open (eg its DC partner is
        # switched off, or its whole substation is dead) is NOT a dead branch: real
        # VSC stations (and OpenLoadFlow) keep the still-connected converter
        # injecting its scheduled P / regulating Q-V as a local device. Only fully
        # deactivate when BOTH stations are disconnected.
        or_disc = is_or_disc or (not line_conn1)
        ex_disc = is_ex_disc or (not line_conn2)
        if or_disc and ex_disc:
            model.deactivate_dcline(hvdc_id)
        elif or_disc:
            model.deactivate_dcline_side1(hvdc_id)
        elif ex_disc:
            model.deactivate_dcline_side2(hvdc_id)
    model.set_dcline_names(df_dc.index)
    if frozen.any():
        model.set_hvdc_ac_emulation_frozen([bool(el) for el in frozen])

    # the VSC stations an outer loop froze at a reactive limit: nothing in the powerflow
    # reads the flag, it only opens them to the physical check of their PQ -> PV release
    if can_be_pv is not None and nb_dc and not (isinstance(can_be_pv, np.ndarray) and can_be_pv.dtype == bool):
        if isinstance(can_be_pv, pd.Series):
            ids = pd.Index(can_be_pv.index[can_be_pv.astype(bool).to_numpy()].astype(str))
        else:
            ids = pd.Index([str(el) for el in can_be_pv])
        side_1 = df_dc["converter_station1_id"].isin(ids).to_numpy(bool)
        side_2 = df_dc["converter_station2_id"].isin(ids).to_numpy(bool)
        if side_1.any() or side_2.any():
            model.set_hvdc_can_be_pv(side_1, side_2)

    return df_dc, hvdc_sub_from_id, hvdc_sub_to_id
