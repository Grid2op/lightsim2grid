# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import pandas as pd

from ._aux_common import _aux_get_bus


def _aux_add_loads(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl, df_dl):
    """Add every load of ``net`` to ``model``, plus the equivalent constant-power
    load at each dangling line's boundary bus (``df_dl``, from
    `_aux_add_buses.py`, empty unless ``convert_dangling_lines=True``). Returns
    ``(df_load, load_sub)``, used by the final substation-id bookkeeping and
    ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_load = net.get_loads().sort_index()
    else:
        df_load = net.get_loads()
    if df_dl.shape[0] > 0:
        # equivalent constant-power load at each dangling line's boundary bus
        # (see `_aux_dangling_lines_fictitious`); p0/q0 are already in MW/MVAr,
        # same raw convention as `net.get_loads()`.
        df_load_extra = pd.DataFrame({
            "p0": df_dl["p0"].values,
            "q0": df_dl["q0"].values,
            "bus_id": df_dl["boundary_bus_id"].values,
            "connected": True,
            "voltage_level_id": df_dl["boundary_vl_id"].values,
        }, index=[f"{dl_id}@boundary_load" for dl_id in df_dl.index])
        df_load = pd.concat([df_load, df_load_extra])
    load_bus, load_disco, load_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "load", df_load)
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     load_bus
                     )
    for load_id, is_disco in enumerate(load_disco):
        if is_disco:
            model.deactivate_load(load_id)
    model.set_load_names(df_load.index)

    return df_load, load_sub
