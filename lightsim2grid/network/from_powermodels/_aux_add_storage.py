# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._bus_remap import pm_bus_to_ls


def _aux_add_storage(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the storage units of `network["storage"]` into the lightsim2grid "model".

    Note
    ----
    PowerModels' `"ps"`/`"qs"` are "Active/Reactive power withdrawn" (same load-like
    sign convention as `"pd"`/`"qd"`), matching lightsim2grid's own `init_storages`
    convention exactly (`LoadContainer`'s doc: "positive storage: the unit is
    charging, power is taken from the grid"), so no sign conversion is needed --
    unlike shunts. Energy/rating/efficiency/inverter-impedance fields have no
    equivalent in lightsim2grid's (loss-free, rating-free) storage model, matching
    how `from_pandapower`'s own storage loader already ignores pandapower's own
    storage ratings.

    Parameters
    ----------
    model
    network: dict
        The PowerModels network data dictionary
    pm_to_ls: dict
        PowerModels bus number (`"bus_i"`) -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`bus_type == 4`) buses

    """
    storage = network.get("storage", {})
    if not storage:
        return

    storage_keys = sorted(storage, key=int)
    storage_bus = pm_bus_to_ls(np.array([int(storage[k]["storage_bus"]) for k in storage_keys]), pm_to_ls)
    ps = np.array([storage[k]["ps"] for k in storage_keys])
    qs = np.array([storage[k].get("qs", 0.0) for k in storage_keys])
    model.init_storages(ps, qs, storage_bus)

    status = np.array([storage[k].get("status", 1) for k in storage_keys]) != 0
    if isolated_ls_bus.size:
        status &= ~np.isin(storage_bus, isolated_ls_bus)
    for st_id, is_ok in enumerate(status):
        if not is_ok:
            model.deactivate_storage(st_id)
