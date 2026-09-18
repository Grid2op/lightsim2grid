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


def _aux_add_trafos(model, net, net_pu, sort_index, voltage_levels, bus_df, first_bus_per_vl,
                    ol_current, keep_half_open_lines, fuse_zero_impedance_branches, fused_trafo_ids):
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
    # lightsim2grid keeps the series impedance right when the shift changes, without any
    # "tap" concept (the per-step r/x deltas pypowsybl carries only in its tap steps).
    ps_alpha, ps_rx_corr = _aux_phase_shift_rx_tables(df_trafo.index, net)
    if any(len(a) for a in ps_alpha):
        model.set_trafo_shift_dependent_rx(True, ps_alpha, ps_rx_corr)

    return df_trafo, tor_sub, tex_sub
