# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np


def mp_bus_to_ls(mp_bus_id, mp_to_ls_converter):
    """Map matpower bus numbers (arbitrary, not necessarily contiguous / 0-based)
    to lightsim2grid contiguous 0-based bus ids, given the dict built by
    `build_bus_remap`.
    """
    return np.array([mp_to_ls_converter[mp_id] for mp_id in mp_bus_id], dtype=int)


def build_bus_remap(bus_i):
    """
    Build the mapping from matpower bus number (`bus[:, BUS_I]`, arbitrary integers,
    not necessarily contiguous / 0-based / sorted) to a contiguous 0-based lightsim2grid
    bus id.

    Parameters
    ----------
    bus_i: numpy array
        The `BUS_I` column of `mpc.bus` (one entry per bus, in the original file order)

    Returns
    -------
    mp_to_ls: dict
        mapping from matpower bus number to lightsim2grid bus id
    ls_to_orig: numpy array
        for each lightsim2grid bus id, the original matpower bus number (useful to
        keep on the model, mirroring `model._ls_to_orig` in the pandapower loader)
    """
    bus_i = np.asarray(bus_i).astype(int)
    if np.unique(bus_i).shape[0] != bus_i.shape[0]:
        raise RuntimeError("Duplicated bus numbers found in `mpc.bus[:, BUS_I]`, "
                           "lightsim2grid cannot handle this.")
    # buses keep the order they appear in mpc.bus (no need to sort: matpower bus
    # numbers are opaque identifiers, not meaningful for ordering)
    ls_to_orig = bus_i
    mp_to_ls = {mp_id: ls_id for ls_id, mp_id in enumerate(bus_i)}
    return mp_to_ls, ls_to_orig
