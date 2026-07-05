# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np

from ._my_const import BUS_I, BUS_TYPE, REF, GEN_BUS, GEN_STATUS
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def _aux_add_slack(model, bus, gen, mp_to_ls, isolated_ls_bus):
    """
    Resolve the slack bus(es) from `mpc.bus[:, BUS_TYPE] == 3` (the matpower reference
    bus) and assign every in-service generator connected to it as a (possibly
    distributed) slack generator, with equal weight.

    Note
    ----
    Unlike pandapower (which has a separate `ext_grid` table to fall back on), matpower
    has no fallback: if the reference bus has no in-service generator connected to it,
    there is no slack source and a `RuntimeError` is raised.

    Several generators can be connected to the reference bus (this is precisely the
    multi-generator-per-bus case this loader is designed to support): they are all
    assigned as slack, each with weight 1.0 (equal-weight distributed slack), consistent
    with the `add_gen_slackbus` weighting convention already used by the pandapower
    loader. Matpower has no column to derive a different weighting from.

    Parameters
    ----------
    model
    bus: numpy array
        raw `mpc.bus` matrix
    gen: numpy array
        raw `mpc.gen` matrix
    mp_to_ls: dict
        matpower bus number -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`BUS_TYPE == 4`) buses

    """
    ref_mask = bus[:, BUS_TYPE] == REF
    n_ref = int(ref_mask.sum())
    if n_ref == 0:
        warnings.warn("No bus with `BUS_TYPE == 3` (reference bus) found. lightsim2grid "
                      "could not assign any slack bus, you will need to do it manually "
                      "with `model.add_gen_slackbus(...)`.")
        return
    if n_ref > 1:
        warnings.warn(f"{n_ref} buses found with `BUS_TYPE == 3` (reference bus), matpower "
                      f"normally expects a single one. All of them will be considered.")

    ref_bus_mp = bus[ref_mask, BUS_I]
    ref_bus_ls = mp_bus_to_ls(ref_bus_mp, mp_to_ls)

    if gen.shape[0] == 0:
        raise RuntimeError("No generator found at all, so no slack bus can be assigned.")

    gen_bus_ls = mp_bus_to_ls(gen[:, GEN_BUS], mp_to_ls)
    gen_in_service = gen[:, GEN_STATUS] > 0
    if isolated_ls_bus.size:
        gen_in_service = gen_in_service & ~np.isin(gen_bus_ls, isolated_ls_bus)

    is_slack_gen = np.isin(gen_bus_ls, ref_bus_ls) & gen_in_service
    slack_gen_ids = np.where(is_slack_gen)[0]
    if slack_gen_ids.size == 0:
        raise RuntimeError(f"Could not find any in-service generator connected to the reference "
                           f"bus(es) {ref_bus_mp.tolist()}. Matpower has no separate slack-bus "
                           f"fallback (unlike pandapower's `ext_grid`): lightsim2grid needs at "
                           f"least one in-service generator at the reference bus to use as slack.")

    for gen_id in slack_gen_ids:
        model.add_gen_slackbus(int(gen_id), 1.0)
