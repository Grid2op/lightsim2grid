# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import pandas as pd

from ._aux_common import _aux_get_bus


def _aux_shunt_sections(model, net, df_shunt, shunt_kv, bus_df, voltage_levels):
    """The sections of every shunt of ``df_shunt`` (cumulative p / q per section count, the
    same convention as ``init_shunt``) and the voltage they regulate."""
    try:
        linear = net.get_linear_shunt_compensator_sections()
        non_linear = net.get_non_linear_shunt_compensator_sections()
    except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
        return
    kv2 = pd.Series(shunt_kv ** 2, index=df_shunt.index)
    count = df_shunt["section_count"].to_numpy(int)
    for k, sid in enumerate(df_shunt.index):
        if sid in linear.index:
            row = linear.loc[sid]
            n = np.arange(1, int(row["max_section_count"]) + 1)
            g, b = n * float(row["g_per_section"]), n * float(row["b_per_section"])
        elif sid in non_linear.index.get_level_values(0):
            sec = non_linear.loc[sid].sort_index()
            g, b = sec["g"].to_numpy(float), sec["b"].to_numpy(float)
        else:
            continue
        model.set_shunt_sections(k, int(count[k]), list(g * kv2[sid]), list(-b * kv2[sid]))
    if "voltage_regulation_on" in df_shunt:
        from ._aux_add_trafos import _aux_bus_pu
        bus, vn = _aux_bus_pu(df_shunt["regulating_bus_id"].fillna("").to_numpy(object), bus_df, voltage_levels)
        target = df_shunt["target_v"].to_numpy(float) / vn
        deadband = df_shunt["target_deadband"].to_numpy(float) / vn
        for k, on in enumerate(df_shunt["voltage_regulation_on"].to_numpy(bool)):
            if not (on or np.isfinite(target[k])):
                continue
            model.set_shunt_section_regulation(
                k, bool(on and bus[k] >= 0), float(target[k]) if np.isfinite(target[k]) else 0.,
                float(deadband[k]) if np.isfinite(deadband[k]) else 0., int(bus[k]))


def _aux_add_shunts(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every shunt compensator of ``net`` to ``model``. Returns
    ``(df_shunt, sh_sub)``, used by the final substation-id bookkeeping and
    ``return_sub_id`` in `initLSGrid.py`."""
    try:
        df_shunt = net.get_shunt_compensators(all_attributes=True)
    except TypeError:  # legacy pypowsybl
        df_shunt = net.get_shunt_compensators()
    if sort_index:
        df_shunt = df_shunt.sort_index()

    sh_bus, sh_disco, sh_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "shunts", df_shunt)
    shunt_kv = voltage_levels.loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(df_shunt["g"].values * shunt_kv**2,
                     -df_shunt["b"].values * shunt_kv**2,
                     sh_bus
                    )
    for shunt_id, disco in enumerate(sh_disco):
        if disco:
           model.deactivate_shunt(shunt_id)
    model.set_shunt_names(df_shunt.index)
    if "section_count" in df_shunt:
        _aux_shunt_sections(model, net, df_shunt, shunt_kv, bus_df, voltage_levels)

    return df_shunt, sh_sub
