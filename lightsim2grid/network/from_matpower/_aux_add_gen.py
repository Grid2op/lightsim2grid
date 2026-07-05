# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._my_const import GEN_BUS, PG, QG, QMAX, QMIN, VG, GEN_STATUS
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def _aux_add_gen(model, gen, mp_to_ls, isolated_ls_bus):
    """
    Add the generators of `mpc.gen` into the lightsim2grid "model".

    Note
    ----
    Several rows of `mpc.gen` can share the same `GEN_BUS`: they are all kept as
    independent generators (no aggregation), which is the entire point of this
    native matpower loader -- lightsim2grid's `LSGrid` has no uniqueness constraint
    on a generator's bus id (`from_pypowsybl` already relies on this for grids with
    several generators on the same bus).

    Matpower has no column equivalent to lightsim2grid's `voltage_regulator_on`: every
    generator is assumed to regulate voltage (standard powerflow convention), which is
    the modeling assumption made here.

    Parameters
    ----------
    model
    gen: numpy array
        raw `mpc.gen` matrix
    mp_to_ls: dict
        matpower bus number -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`BUS_TYPE == 4`) buses

    """
    if gen.shape[0] == 0:
        return

    gen_bus = mp_bus_to_ls(gen[:, GEN_BUS], mp_to_ls)
    voltage_regulator_on = [True] * gen.shape[0]
    model.init_generators_full(gen[:, PG], gen[:, VG], gen[:, QG], voltage_regulator_on,
                               gen[:, QMIN], gen[:, QMAX], gen_bus)

    gen_status = gen[:, GEN_STATUS] > 0
    if isolated_ls_bus.size:
        gen_status &= ~np.isin(gen_bus, isolated_ls_bus)
    for gen_id, is_ok in enumerate(gen_status):
        if not is_ok:
            model.deactivate_gen(gen_id)
