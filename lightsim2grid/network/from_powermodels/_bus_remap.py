# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np


def pm_bus_to_ls(pm_bus_id, pm_to_ls_converter):
    """Map PowerModels bus numbers (``"bus_i"``, arbitrary integers, not necessarily
    contiguous / 0-based) to lightsim2grid contiguous 0-based bus ids, given the dict
    built by `build_bus_remap`.
    """
    return np.array([pm_to_ls_converter[pm_id] for pm_id in pm_bus_id], dtype=int)


def build_bus_remap(bus_keys, network):
    """
    Build the mapping from PowerModels bus number (``network["bus"][k]["bus_i"]``,
    arbitrary integers, not necessarily contiguous / 0-based / sorted) to a contiguous
    0-based lightsim2grid bus id.

    Parameters
    ----------
    bus_keys: list of str
        ``network["bus"]`` dict keys, in the order lightsim2grid buses should be laid
        out (see `_aux_check_legit`/`initLSGrid.init`: ``sorted(network["bus"], key=int)``)
    network: dict
        The PowerModels network data dictionary

    Returns
    -------
    pm_to_ls: dict
        mapping from PowerModels bus number (``bus_i``) to lightsim2grid bus id
    ls_to_orig: numpy array
        for each lightsim2grid bus id, the original PowerModels bus number (useful to
        keep on the model, mirroring `model._ls_to_orig` in the pandapower loader)
    """
    bus_i = np.array([int(network["bus"][k]["bus_i"]) for k in bus_keys], dtype=int)
    if np.unique(bus_i).shape[0] != bus_i.shape[0]:
        raise RuntimeError('Duplicated "bus_i" values found in this PowerModels network, '
                           "lightsim2grid cannot handle this.")
    ls_to_orig = bus_i
    pm_to_ls = {pm_id: ls_id for ls_id, pm_id in enumerate(bus_i)}
    return pm_to_ls, ls_to_orig
