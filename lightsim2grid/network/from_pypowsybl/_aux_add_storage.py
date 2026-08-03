# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._aux_common import _aux_get_bus


def _aux_add_storage(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every storage unit (IIDM battery) of ``net`` to ``model``. IIDM gives
    the battery setpoints in the *generator* convention (positive target_p =
    power produced / injected) while lightsim2grid stores storage as PQ in the
    *load* convention (positive = power drawn from the grid, *ie* charging), same
    as pandapower and grid2op. We negate to convert, and sanitize NaN (IIDM
    allows an unset target_q). Returns ``(df_batt, batt_sub)``, used by the final
    substation-id bookkeeping and ``return_sub_id`` in `initLSGrid.py`."""
    if sort_index:
        df_batt = net.get_batteries().sort_index()
    else:
        df_batt = net.get_batteries()
    batt_bus, batt_disco, batt_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "storage", df_batt)
    batt_p = df_batt["target_p"].values.astype(float)
    batt_q = df_batt["target_q"].values.astype(float)
    batt_p = np.where(np.isfinite(batt_p), batt_p, 0.)
    batt_q = np.where(np.isfinite(batt_q), batt_q, 0.)
    model.init_storages(-batt_p,  # IIDM generator convention -> lightsim2grid load convention
                        -batt_q,
                        batt_bus
                        )
    for batt_id, disco in enumerate(batt_disco):
        if disco:
           model.deactivate_storage(batt_id)
    model.set_storage_names(df_batt.index)

    return df_batt, batt_sub
