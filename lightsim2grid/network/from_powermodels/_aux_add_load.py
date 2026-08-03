# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._bus_remap import pm_bus_to_ls


def _aux_add_load(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the loads of `network["load"]` into the lightsim2grid "model".

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
    load = network.get("load", {})
    if not load:
        return

    load_keys = sorted(load, key=int)
    load_bus = pm_bus_to_ls(np.array([int(load[k]["load_bus"]) for k in load_keys]), pm_to_ls)
    pd = np.array([load[k]["pd"] for k in load_keys])
    qd = np.array([load[k]["qd"] for k in load_keys])
    model.init_loads(pd, qd, load_bus)

    status = np.array([load[k].get("status", 1) for k in load_keys]) != 0
    if isolated_ls_bus.size:
        status &= ~np.isin(load_bus, isolated_ls_bus)
    for load_id, is_ok in enumerate(status):
        if not is_ok:
            model.deactivate_load(load_id)
