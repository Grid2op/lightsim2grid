# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import pandas as pd

from ._aux_common import _aux_get_bus, _aux_current_limits
from ._my_const import _DANGLING_BOUNDARY_LINE_SUFFIX


def _aux_add_lines(model, net, net_pu, sort_index, voltage_levels, bus_df, first_bus_per_vl,
                   df_dl, ol_current, keep_half_open_lines, fuse_zero_impedance_branches, fused_line_ids):
    """Add every line of ``net`` to ``model``, plus the equivalent branch (local
    bus -> fictitious boundary bus) for each dangling line (``df_dl``, from
    `_aux_add_buses.py`). ``ol_current`` (``net.get_operational_limits()``
    filtered to CURRENT, or ``None``) is shared with `_aux_add_trafos.py`, computed
    once in `initLSGrid.py`. Returns ``(df_line, lor_sub, lex_sub)``, used by the
    final substation-id bookkeeping and ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_line = net.get_lines().sort_index()
    else:
        df_line = net.get_lines()
    try:
        line_limit_groups = net.get_lines(all_attributes=True)[["selected_limits_group_1", "selected_limits_group_2"]]
    except (TypeError, KeyError):
        # not available on legacy pypowsybl / grid with no limit group at all
        line_limit_groups = pd.DataFrame(
            {"selected_limits_group_1": np.nan, "selected_limits_group_2": np.nan}, index=df_line.index
        )

    df_line_pu = net_pu.get_lines().loc[df_line.index]
    if df_dl.shape[0] > 0:
        # equivalent branch (local bus -> fictitious boundary bus) for each
        # dangling line, see `_aux_dangling_lines_fictitious`. Shunt admittance
        # (g, b) sits entirely on the local (network) side, none on the
        # boundary side -- same convention already used for the single `h` of
        # a transformer (see `trafo_h` in `_aux_add_trafos.py`).
        df_dl_pu = net_pu.get_dangling_lines().loc[df_dl.index]
        line_ids_extra = [f"{dl_id}{_DANGLING_BOUNDARY_LINE_SUFFIX}" for dl_id in df_dl.index]
        df_line_extra = pd.DataFrame({
            "bus1_id": df_dl["bus_id"].values,
            "bus2_id": df_dl["boundary_bus_id"].values,
            "connected1": df_dl["connected"].values,
            "connected2": True,
            "voltage_level1_id": df_dl["voltage_level_id"].values,
            "voltage_level2_id": df_dl["boundary_vl_id"].values,
        }, index=line_ids_extra)
        df_line_extra_pu = pd.DataFrame({
            "r": df_dl_pu["r"].values,
            "x": df_dl_pu["x"].values,
            "g1": df_dl_pu["g"].values,
            "b1": df_dl_pu["b"].values,
            "g2": 0.,
            "b2": 0.,
        }, index=line_ids_extra)
        df_line = pd.concat([df_line, df_line_extra])
        df_line_pu = pd.concat([df_line_pu, df_line_extra_pu])
    line_r = df_line_pu["r"].values
    line_x = df_line_pu["x"].values
    line_h_or = (df_line_pu["g1"].values + 1j * df_line_pu["b1"].values)
    line_h_ex = (df_line_pu["g2"].values + 1j * df_line_pu["b2"].values)
    lor_bus, lor_disco, lor_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "line (side 1)", df_line, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    lex_bus, lex_disco, lex_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "line (side 2)", df_line, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_powerlines_full(line_r,
                               line_x,
                               line_h_or,
                               line_h_ex,
                               lor_bus,
                               lex_bus
                              )
    for line_id, (is_or_disc, is_ex_disc) in enumerate(zip(lor_disco, lex_disco)):
        if is_or_disc and is_ex_disc:
            model.deactivate_powerline(line_id)
        elif is_or_disc:
            model.deactivate_powerline_side1(line_id) if keep_half_open_lines else model.deactivate_powerline(line_id)
        elif is_ex_disc:
            model.deactivate_powerline_side2(line_id) if keep_half_open_lines else model.deactivate_powerline(line_id)
        elif fuse_zero_impedance_branches and df_line.index[line_id] in fused_line_ids:
            # both terminal buses already fused into one node above: this line
            # would otherwise contribute a 1/Z admittance (Inf for an exact zero)
            model.deactivate_powerline(line_id)
    model.set_line_names(df_line.index)
    line_limit_a1_ka, line_limit_a2_ka = _aux_current_limits(
        df_line.index,
        line_limit_groups["selected_limits_group_1"],
        line_limit_groups["selected_limits_group_2"],
        ol_current,
    )
    model.set_line_current_limit_side1(line_limit_a1_ka)
    model.set_line_current_limit_side2(line_limit_a2_ka)

    return df_line, lor_sub, lex_sub
