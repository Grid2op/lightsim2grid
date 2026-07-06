# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._bus_remap import pm_bus_to_ls


def _aux_add_gen(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the generators of `network["gen"]` into the lightsim2grid "model". Several
    generators can share the same `"gen_bus"`: they are all kept as independent
    generators (no aggregation), matching lightsim2grid's own `GeneratorContainer`
    (no uniqueness constraint on a generator's bus id).

    PowerModels has no field equivalent to lightsim2grid's `voltage_regulator_on`:
    every generator is assumed to regulate voltage (standard powerflow convention).

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
    gen = network.get("gen", {})
    if not gen:
        return

    gen_keys = sorted(gen, key=int)
    gen_bus = pm_bus_to_ls(np.array([int(gen[k]["gen_bus"]) for k in gen_keys]), pm_to_ls)
    pg = np.array([gen[k]["pg"] for k in gen_keys])
    vg = np.array([gen[k]["vg"] for k in gen_keys])
    qg = np.array([gen[k].get("qg", 0.0) for k in gen_keys])
    qmin = np.array([gen[k]["qmin"] for k in gen_keys])
    qmax = np.array([gen[k]["qmax"] for k in gen_keys])
    voltage_regulator_on = [True] * len(gen_keys)
    model.init_generators_full(pg, vg, qg, voltage_regulator_on, qmin, qmax, gen_bus)

    gen_status = np.array([gen[k].get("gen_status", 1) for k in gen_keys]) > 0
    if isolated_ls_bus.size:
        gen_status &= ~np.isin(gen_bus, isolated_ls_bus)
    for gen_id, is_ok in enumerate(gen_status):
        if not is_ok:
            model.deactivate_gen(gen_id)
