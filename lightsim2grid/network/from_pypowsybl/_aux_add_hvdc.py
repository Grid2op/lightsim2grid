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


def _aux_add_hvdc(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every HVDC line of ``net`` (VSC / LCC converter stations, possibly
    carrying the angle-droop ("AC emulation") extension) to ``model``. Returns
    ``(df_dc, hvdc_sub_from_id, hvdc_sub_to_id)``, used by the final
    substation-id bookkeeping and ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_dc = net.get_hvdc_lines().sort_index()
    else:
        df_dc = net.get_hvdc_lines()
    df_vsc = net.get_vsc_converter_stations()
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
        vl_kv = voltage_levels.loc[df_side["voltage_level_id"].values]["nominal_v"].values
        vset_pu = df_side["target_v"].values / vl_kv
        vset_pu = np.where(np.isfinite(vset_pu), vset_pu, 1.0)
        qset = df_side["target_q"].values
        qset = np.where(np.isfinite(qset), qset, 0.)
        min_q = df_side["min_q"].values
        min_q = np.where(np.isfinite(min_q), min_q, -_max_hvdc_mva)
        max_q = df_side["max_q"].values
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

    # the angle-droop active power control ("AC emulation"), an IIDM extension
    droop_enabled = np.zeros(nb_dc, dtype=bool)
    droop_p0_mw = np.zeros(nb_dc)
    droop_mw_per_deg = np.zeros(nb_dc)
    try:
        df_droop = net.get_extensions("hvdcAngleDroopActivePowerControl")
    except Exception:
        # extension tables may be unavailable on (very) old pypowsybl versions
        df_droop = None
    if df_droop is not None and df_droop.shape[0]:
        for hvdc_pos, line_id in enumerate(df_dc.index):
            if line_id not in df_droop.index:
                continue
            if not bool(df_droop.loc[line_id, "enabled"]):
                continue
            droop_enabled[hvdc_pos] = True
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
                          max_p_mw,  # pmax 1 -> 2: IIDM has a single max_p (open-loadflow convention)
                          max_p_mw,  # pmax 2 -> 1
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

    return df_dc, hvdc_sub_from_id, hvdc_sub_to_id
