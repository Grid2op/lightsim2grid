# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._my_const import BUS_I, PD, QD
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def _aux_add_load(model, bus, mp_to_ls, isolated_ls_bus):
    """
    Add the loads described by the `PD` / `QD` columns of `mpc.bus` into the
    lightsim2grid "model". Matpower has no separate load table: a load is only
    created for buses with a non zero `PD` or `QD`.

    Parameters
    ----------
    model
    bus: numpy array
        raw `mpc.bus` matrix
    mp_to_ls: dict
        matpower bus number -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`BUS_TYPE == 4`) buses

    """
    has_load = (bus[:, PD] != 0.) | (bus[:, QD] != 0.)
    if not np.any(has_load):
        return

    load_bus = mp_bus_to_ls(bus[has_load, BUS_I], mp_to_ls)
    model.init_loads(bus[has_load, PD], bus[has_load, QD], load_bus)

    if isolated_ls_bus.size:
        for load_id, is_isolated in enumerate(np.isin(load_bus, isolated_ls_bus)):
            if is_isolated:
                model.deactivate_load(load_id)
