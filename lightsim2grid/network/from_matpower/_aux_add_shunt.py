# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._my_const import BUS_I, GS, BS
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def _aux_add_shunt(model, bus, mp_to_ls, isolated_ls_bus):
    """
    Add the shunts described by the `GS` / `BS` columns of `mpc.bus` into the
    lightsim2grid "model". Matpower has no separate shunt table: a shunt is only
    created for buses with a non zero `GS` or `BS`.

    Note
    ----
    Matpower's `GS` is "shunt conductance (MW **demanded** at V = 1.0 p.u.)" (load
    convention, positive = consumed) while `BS` is "shunt susceptance (MVAr
    **injected** at V = 1.0 p.u.)" (generator convention, positive = supplied) -- the
    opposite sign convention from lightsim2grid/pandapower's `q_mvar` (load
    convention, positive = consumed). `BS` must therefore be negated. This matches
    pandapower's own matpower/ppc importer (`pandapower.converter.pypower.from_ppc`),
    which does `q_mvar=-ppc["bus"][is_shunt, BS]`.

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
    has_shunt = (bus[:, GS] != 0.) | (bus[:, BS] != 0.)
    if not np.any(has_shunt):
        return

    shunt_bus = mp_bus_to_ls(bus[has_shunt, BUS_I], mp_to_ls)
    model.init_shunt(bus[has_shunt, GS], -bus[has_shunt, BS], shunt_bus)

    if isolated_ls_bus.size:
        for sh_id, is_isolated in enumerate(np.isin(shunt_bus, isolated_ls_bus)):
            if is_isolated:
                model.deactivate_shunt(sh_id)
