# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np

from ._my_const import REF
from ._bus_remap import pm_bus_to_ls


def _aux_add_slack(model, network, pm_to_ls, isolated_ls_bus):
    """
    Resolve the slack bus(es) from `network["bus"][k]["bus_type"] == 3` (the
    PowerModels/matpower reference bus) and assign every in-service generator
    connected to it as a (possibly distributed) slack generator, with equal weight.

    Note
    ----
    Unlike pandapower (which has a separate `ext_grid` table to fall back on),
    PowerModels/matpower has no fallback: if the reference bus has no in-service
    generator connected to it, there is no slack source and a `RuntimeError` is raised.

    Several generators can be connected to the reference bus: they are all assigned as
    slack, each with weight 1.0 (equal-weight distributed slack), consistent with the
    `add_gen_slackbus` weighting convention already used by the pandapower loader.

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
    bus = network["bus"]
    ref_keys = [k for k, b in bus.items() if b.get("bus_type") == REF]
    if len(ref_keys) == 0:
        warnings.warn("No bus with `bus_type == 3` (reference bus) found. lightsim2grid "
                      "could not assign any slack bus, you will need to do it manually "
                      "with `model.add_gen_slackbus(...)`.")
        return
    if len(ref_keys) > 1:
        warnings.warn(f"{len(ref_keys)} buses found with `bus_type == 3` (reference bus), "
                      "PowerModels/matpower normally expects a single one. All of them will "
                      "be considered, with equal weights.")

    ref_bus_pm = np.array([int(bus[k]["bus_i"]) for k in ref_keys])
    ref_bus_ls = pm_bus_to_ls(ref_bus_pm, pm_to_ls)

    gen = network.get("gen", {})
    if not gen:
        raise RuntimeError("No generator found at all, so no slack bus can be assigned.")

    gen_keys = sorted(gen, key=int)
    gen_bus_ls = pm_bus_to_ls(np.array([int(gen[k]["gen_bus"]) for k in gen_keys]), pm_to_ls)
    gen_in_service = np.array([gen[k].get("gen_status", 1) for k in gen_keys]) > 0
    if isolated_ls_bus.size:
        gen_in_service = gen_in_service & ~np.isin(gen_bus_ls, isolated_ls_bus)

    is_slack_gen = np.isin(gen_bus_ls, ref_bus_ls) & gen_in_service
    slack_gen_ids = np.where(is_slack_gen)[0]
    if slack_gen_ids.size == 0:
        raise RuntimeError(f"Could not find any in-service generator connected to the reference "
                           f"bus(es) {ref_bus_pm.tolist()}. PowerModels/matpower has no separate "
                           "slack-bus fallback (unlike pandapower's `ext_grid`): lightsim2grid "
                           "needs at least one in-service generator at the reference bus to use "
                           "as slack.")

    for gen_id in slack_gen_ids:
        model.add_gen_slackbus(int(gen_id), 1.0)
