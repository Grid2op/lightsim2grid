# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import pandas as pd

from ._aux_common import (
    _aux_get_bus,
    _aux_current_limits,
    _aux_get_2wt_all_attrs,
    _aux_trafo_alpha,
    _aux_trafo_rho,
)


def _aux_phase_shift_rx_tables(trafo_index, net):
    """Per-transformer phase-shift -> series-impedance dependency for
    :meth:`LSGrid.set_trafo_shift_dependent_rx`.

    Returns two lists-of-lists aligned to ``trafo_index``: the phase-shift sample
    points ``alpha`` (in **radian**, matching lightsim2grid's ``shift_``) and the
    matching r/x correction (**in percent**), read from the pypowsybl
    phase-tap-changer steps. Inner lists are empty for transformers without a
    phase-tap-changer.

    pypowsybl exposes the transformer r / x at the *neutral* tap and only folds the
    tap into ``rho`` / ``alpha``; the per-step r/x deltas (percent) are dropped. They
    matter for phase-shifting transformers whose series impedance varies with the
    shift (phase-shifting transformers on real grids: tens of MW of through-flow).
    On those grids ``r% == x%`` for
    every step, so a single correction table is applied to both r and x."""
    n = len(trafo_index)
    alpha = [[] for _ in range(n)]
    corr = [[] for _ in range(n)]
    try:
        ptc = net.get_phase_tap_changers()
        steps = net.get_phase_tap_changer_steps()
    except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
        return alpha, corr
    if ptc.shape[0] == 0 or steps.shape[0] == 0:
        return alpha, corr
    have = set(ptc.index)
    for i, tid in enumerate(trafo_index):
        if tid not in have:
            continue
        try:
            st = steps.loc[tid]
        except KeyError:
            continue
        a = np.deg2rad(np.atleast_1d(st["alpha"].to_numpy(dtype=float)))
        x = np.atleast_1d(st["x"].to_numpy(dtype=float))  # r% == x% on these PSTs
        alpha[i] = [float(v) for v in a]
        corr[i] = [float(v) for v in x]
    return alpha, corr


_PHASE_REGULATION_MODES = {"CURRENT_LIMITER": "CURRENT_LIMITER", "ACTIVE_POWER_CONTROL": "ACTIVE_POWER"}
_SIDES = {"ONE": 1, "TWO": 2}


def _aux_bus_pu(bus_ids, bus_df, voltage_levels):
    """The grid bus id (-1 if unknown) and the nominal voltage (kV, NaN) of each bus-view id."""
    known = pd.Index(bus_ids).isin(bus_df.index)
    bus = np.full(len(bus_ids), -1, dtype=int)
    vn = np.full(len(bus_ids), np.nan)
    if known.any():
        rows = bus_df.loc[np.asarray(bus_ids)[known]]
        bus[known] = rows["bus_global_id"].to_numpy(int)
        vn[known] = voltage_levels.loc[rows["voltage_level_id"].values, "nominal_v"].to_numpy(float)
    return bus, vn


def _aux_tap_changers(model, net, trafo_index, bus_df, voltage_levels, olf_rules=None):
    """The ratio and phase tap changers of the 2-winding transformers ``trafo_index``: their
    step tables and positions (the pi model is then taken at the taps, as OpenLoadFlow takes
    it), and what they regulate. Legs of 3-winding transformers are not in ``trafo_index``.
    With ``olf_rules``, a ratio tap changer with an implausible target voltage does not regulate
    (OpenLoadFlow's VoltageControl.checkTargetV, as for the generators)."""
    pos = pd.Series(np.arange(len(trafo_index)), index=trafo_index)
    for phase in (False, True):
        try:
            if phase:
                changers = net.get_phase_tap_changers(all_attributes=True)
                steps = net.get_phase_tap_changer_steps(all_attributes=True)
            else:
                changers = net.get_ratio_tap_changers(all_attributes=True)
                steps = net.get_ratio_tap_changer_steps(all_attributes=True)
        except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
            continue
        changers = changers.loc[changers.index.isin(trafo_index)]
        if not len(changers):
            continue
        steps = steps.loc[steps.index.get_level_values(0).isin(changers.index)]
        for tid, st in steps.groupby(level=0, sort=False):
            row = changers.loc[tid]
            st = st.sort_index(level=1)
            args = [int(pos[tid]), int(row["low_tap"]), int(row["tap"]), st["rho"].to_numpy(float)]
            if phase:
                args.append(st["alpha"].to_numpy(float))
            args += [st[col].to_numpy(float) for col in ("r", "x", "g", "b")]
            if phase:
                model.set_trafo_phase_tap_changer(*args)
            else:
                model.set_trafo_ratio_tap_changer(*args)
        # what they regulate
        if phase:
            from lightsim2grid.lightsim2grid_cpp import RegulationMode
            for tid, row in changers.iterrows():
                side = _SIDES.get(row.get("regulated_side", ""), 0)
                if not side:
                    continue  # regulating something else than the transformer itself
                mode = getattr(RegulationMode, _PHASE_REGULATION_MODES.get(row["regulation_mode"], "FIXED"))
                value = float(row["regulation_value"]) if np.isfinite(row["regulation_value"]) else 0.
                deadband = float(row["target_deadband"]) if np.isfinite(row["target_deadband"]) else 0.
                model.set_trafo_phase_tap_regulation(int(pos[tid]), mode, bool(row["regulating"]), value, deadband,
                                                     side)
        else:
            bus, vn = _aux_bus_pu(changers["regulating_bus_id"].fillna("").to_numpy(object), bus_df, voltage_levels)
            regulating = changers["regulating"].to_numpy(bool)
            if "oltc" in changers:
                # OpenLoadFlow only regulates with a changer able to move on load
                regulating &= changers["oltc"].to_numpy(bool)
            target = changers["target_v"].to_numpy(float) / vn
            deadband = changers["target_deadband"].to_numpy(float) / vn
            if olf_rules is not None:
                implausible = ((vn > olf_rules.min_nominal_voltage_target_voltage_check) &
                               ((target < olf_rules.min_plausible_target_v) | (target > olf_rules.max_plausible_target_v)))
                regulating &= ~implausible
            for k, tid in enumerate(changers.index):
                model.set_trafo_ratio_tap_regulation(
                    int(pos[tid]), bool(regulating[k] and bus[k] >= 0),
                    float(target[k]) if np.isfinite(target[k]) else 0.,
                    float(deadband[k]) if np.isfinite(deadband[k]) else 0., int(bus[k]))


def _aux_add_trafos(model, net, net_pu, sort_index, voltage_levels, bus_df, first_bus_per_vl,
                    ol_current, keep_half_open_lines, fuse_zero_impedance_branches, fused_trafo_ids, olf_rules=None):
    """Add every 2-winding transformer of ``net`` to ``model``. ``ol_current``
    (``net.get_operational_limits()`` filtered to CURRENT, or ``None``) is shared
    with `_aux_add_lines.py`, computed once in `initLSGrid.py`. Returns
    ``(df_trafo, tor_sub, tex_sub)``, used by the final substation-id bookkeeping
    and ``return_sub_id`` in `initLSGrid.py`."""
    # I extract trafo with `all_attributes=True` so that I have access to the `rho`
    df_trafo_not_sorted = _aux_get_2wt_all_attrs(net)

    if sort_index:
        df_trafo = df_trafo_not_sorted.sort_index()
    else:
        df_trafo = df_trafo_not_sorted

    df_trafo_pu = _aux_get_2wt_all_attrs(net_pu)
    df_trafo_pu = df_trafo_pu.loc[df_trafo.index]
    ratio_tap_changer = net_pu.get_ratio_tap_changers()

    shift_ = _aux_trafo_alpha(df_trafo_pu, net)
    # tap is side 2 in IIDM
    is_tap_side1 = np.zeros(df_trafo.shape[0], dtype=bool)
    # neutral-tap impedance (the phase-shift -> r/x dependence of phase-shifting
    # transformers is handled by lightsim2grid as a function of the shift alpha, see
    # the model.set_trafo_shift_dependent_rx(...) call below).
    trafo_r = df_trafo_pu["r"].values
    trafo_x = df_trafo_pu["x"].values
    trafo_h = (df_trafo_pu["g"].values + 1j * df_trafo_pu["b"].values)

    # now get the ratio
    # in lightsim2grid (cpp)
    ratio = _aux_trafo_rho(df_trafo_pu, ratio_tap_changer)

    tor_bus, tor_disco, tor_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "trafo (side 1)", df_trafo, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    tex_bus, tex_disco, tex_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "trafo (side 2)", df_trafo, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_h,
                     ratio,
                     shift_,  # in degree !
                     is_tap_side1,
                     tor_bus,
                     tex_bus,
                     False,  # ignore_tap_side_for_phase_shift is False for pypowsybl
                     )
    for t_id, (is_or_disc, is_ex_disc) in enumerate(zip(tor_disco, tex_disco)):
        if is_or_disc and is_ex_disc:
            model.deactivate_trafo(t_id)
        elif is_or_disc:
            model.deactivate_trafo_side1(t_id) if keep_half_open_lines else model.deactivate_trafo(t_id)
        elif is_ex_disc:
            model.deactivate_trafo_side2(t_id) if keep_half_open_lines else model.deactivate_trafo(t_id)
        elif fuse_zero_impedance_branches and df_trafo.index[t_id] in fused_trafo_ids:
            # both terminal buses already fused into one node above
            model.deactivate_trafo(t_id)
    model.set_trafo_names(df_trafo.index)
    if "selected_limits_group_1" in df_trafo.columns:
        trafo_group_1 = df_trafo["selected_limits_group_1"]
        trafo_group_2 = df_trafo["selected_limits_group_2"]
    else:
        # not available on legacy pypowsybl
        trafo_group_1 = pd.Series(np.nan, index=df_trafo.index)
        trafo_group_2 = pd.Series(np.nan, index=df_trafo.index)
    trafo_limit_a1_ka, trafo_limit_a2_ka = _aux_current_limits(
        df_trafo.index, trafo_group_1, trafo_group_2, ol_current
    )
    model.set_trafo_current_limit_side1(trafo_limit_a1_ka)
    model.set_trafo_current_limit_side2(trafo_limit_a2_ka)
    # phase-shifting transformers: declare the (alpha -> r/x correction) dependency so
    # lightsim2grid keeps the series impedance right when the shift changes through
    # change_shift_trafo (between two taps: at a tap, the tap changer below decides).
    ps_alpha, ps_rx_corr = _aux_phase_shift_rx_tables(df_trafo.index, net)
    if any(len(a) for a in ps_alpha):
        model.set_trafo_shift_dependent_rx(True, ps_alpha, ps_rx_corr)
    # the tap changers: the pi model at the taps (r, x, g, b corrected by each step)
    _aux_tap_changers(model, net, df_trafo.index, bus_df, voltage_levels, olf_rules=olf_rules)

    return df_trafo, tor_sub, tex_sub
