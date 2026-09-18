# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._bus_remap import pm_bus_to_ls


def _aux_add_dc_line(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the HVDC lines of `network["dcline"]` into the lightsim2grid "model", using
    the same simplified (station-loss-free, resistance-free) dc line model already
    used by lightsim2grid's pandapower and matpower loaders (`model.init_dclines`).

    Note
    ----
    PowerModels.jl keeps `"pf"` at matpower's own value/sign when parsing a matpower
    case (only `"pt"`/`"qf"`/`"qt"` get sign-flipped by PowerModels.jl's parser, see
    `PowerModels.jl/src/io/matpower.jl`'s `_mp2pm_dcline!`), so the same sign
    convention as lightsim2grid's `from_matpower` loader applies: `"pf"` is the MW
    flow "from" -> "to" measured at the `"f_bus"` end (positive means the `"f_bus"`
    side sends power into the line), while `model.init_dclines` expects the power
    *received* at side 1 (the `"f_bus"` side) when positive, so `"pf"` must be negated.

    The loss model is `pt = pf - (loss0 + loss1 * pf)` with `loss1` a fraction of
    `pf`, while lightsim2grid/pandapower's `loss_percent` expresses that same
    coefficient as a percentage (`loss1 * 100`) and `loss_mw` is `loss0` directly.

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
    dcline = network.get("dcline", {})
    if not dcline:
        return

    dc_keys = sorted(dcline, key=int)
    f_bus = pm_bus_to_ls(np.array([int(dcline[k]["f_bus"]) for k in dc_keys]), pm_to_ls)
    t_bus = pm_bus_to_ls(np.array([int(dcline[k]["t_bus"]) for k in dc_keys]), pm_to_ls)

    p_mw = -np.array([dcline[k]["pf"] for k in dc_keys])
    loss_percent = 100. * np.array([dcline[k]["loss1"] for k in dc_keys])
    loss_mw = np.array([dcline[k]["loss0"] for k in dc_keys])
    vf = np.array([dcline[k].get("vf", 1.0) for k in dc_keys])
    vt = np.array([dcline[k].get("vt", 1.0) for k in dc_keys])
    qminf = np.array([dcline[k]["qminf"] for k in dc_keys])
    qmaxf = np.array([dcline[k]["qmaxf"] for k in dc_keys])
    qmint = np.array([dcline[k]["qmint"] for k in dc_keys])
    qmaxt = np.array([dcline[k]["qmaxt"] for k in dc_keys])

    model.init_dclines(f_bus, t_bus, p_mw, loss_percent, loss_mw, vf, vt,
                       qminf, qmaxf, qmint, qmaxt)

    status = np.array([dcline[k].get("br_status", 1) for k in dc_keys]) != 0
    if isolated_ls_bus.size:
        status &= ~np.isin(f_bus, isolated_ls_bus)
        status &= ~np.isin(t_bus, isolated_ls_bus)
    for dc_id, is_ok in enumerate(status):
        if not is_ok:
            model.deactivate_dcline(dc_id)
