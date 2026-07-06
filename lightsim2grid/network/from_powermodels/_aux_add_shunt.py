# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._bus_remap import pm_bus_to_ls


def _aux_add_shunt(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the shunts of `network["shunt"]` into the lightsim2grid "model".

    Note
    ----
    PowerModels' shunt admittance is `Y = gs + j*bs` (per-unit, standard convention,
    see https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/ and
    PowerModels.jl's own `src/core/data.jl`: `p_delta += gs*vm^2; q_delta -= bs*vm^2`).
    lightsim2grid's `ShuntContainer` instead computes `Y_pu = (p_mw - j*q_mvar)/sn_mva`
    (see its own `// TODO check the sign here for p_mw, it is suspicious!` comment).
    Equating the two gives `p_mw = gs*sn_mva`, `q_mvar = -bs*sn_mva` -- the same sign
    flip lightsim2grid's `from_matpower` loader already applies to matpower's `BS`
    (verified independently, twice, against that loader's own docstring and code).

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
    shunt = network.get("shunt", {})
    if not shunt:
        return

    sn_mva = float(network["baseMVA"])
    shunt_keys = sorted(shunt, key=int)
    shunt_bus = pm_bus_to_ls(np.array([int(shunt[k]["shunt_bus"]) for k in shunt_keys]), pm_to_ls)
    p_mw = np.array([shunt[k]["gs"] for k in shunt_keys]) * sn_mva
    q_mvar = -np.array([shunt[k]["bs"] for k in shunt_keys]) * sn_mva
    model.init_shunt(p_mw, q_mvar, shunt_bus)

    status = np.array([shunt[k].get("status", 1) for k in shunt_keys]) != 0
    if isolated_ls_bus.size:
        status &= ~np.isin(shunt_bus, isolated_ls_bus)
    for sh_id, is_ok in enumerate(status):
        if not is_ok:
            model.deactivate_shunt(sh_id)
