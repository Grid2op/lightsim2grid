# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from ._aux_common import _aux_get_bus


def _aux_add_shunts(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every shunt compensator of ``net`` to ``model``. Returns
    ``(df_shunt, sh_sub)``, used by the final substation-id bookkeeping and
    ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_shunt = net.get_shunt_compensators().sort_index()
    else:
        df_shunt = net.get_shunt_compensators()

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

    return df_shunt, sh_sub
