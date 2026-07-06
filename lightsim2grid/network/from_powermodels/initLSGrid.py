# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Natively parse a PowerModels.jl network data dictionary (see
https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/) and initialize a
LSGrid c++ object from it directly (bus, branch, gen, load, shunt, storage and dcline),
without ever going through a raw MATPOWER matrix or a pandapower/pypowsybl network.

This is the shared, feature-complete engine other loaders build on: `from_matpower`
converts its raw `bus`/`gen`/`branch`/`dcline` matrices into a PowerModels-style dict
and delegates here, and `from_powermodels.init_from_pfdelta` (PFΔ dataset rows, see
arXiv:2510.22048) unwraps its row wrapper and delegates here directly, since PFΔ rows
already are PowerModels dicts.
"""

from typing import Optional
import numpy as np

from ...lightsim2grid_cpp import LSGrid
from ._bus_remap import build_bus_remap, pm_bus_to_ls
from ._aux_check_legit import _aux_check_legit
from ._aux_add_branch import _aux_add_branch, classify_branches
from ._aux_add_load import _aux_add_load
from ._aux_add_shunt import _aux_add_shunt
from ._aux_add_storage import _aux_add_storage
from ._aux_add_gen import _aux_add_gen
from ._aux_add_slack import _aux_add_slack
from ._aux_add_dc_line import _aux_add_dc_line
from ._my_const import NONE


def init(network: dict,
         n_busbar_per_sub: Optional[int] = None,  # max number of buses allowed per substation / voltage level
         ) -> LSGrid:
    """
    Convert a PowerModels.jl network data dictionary into a LSGrid.

    Parameters
    ----------
    network: dict
        A PowerModels network data dictionary (the top-level dict with `"bus"` /
        `"branch"` / `"gen"` / ... keys -- *not* a full PFΔ dataset row, which wraps
        this dict under a `"network"` key; use `init_from_pfdelta` for that).
    n_busbar_per_sub:
        There is always exactly one substation / voltage level per PowerModels bus
        (PowerModels has no notion of several busbar sections within a bus, so this is
        not configurable). This parameter only controls how many buses / busbar
        sections lightsim2grid allocates *per substation*, which is useful if you
        intend to perform grid2op-like topology actions on the resulting grid
        afterwards. Defaults to 1 (no extra busbar section). Any extra busbar section
        is deactivated, since nothing in the base network is ever connected to it.

    Returns
    -------
    model: :class:`lightsim2grid.network.LSGrid`
        The initialized network

    """
    _aux_check_legit(network)

    bus_keys = sorted(network["bus"], key=int)
    pm_to_ls, ls_to_orig = build_bus_remap(bus_keys, network)

    line_keys, trafo_keys = classify_branches(network)
    nb_line, nb_trafo = len(line_keys), len(trafo_keys)

    model = LSGrid()
    model.set_sn_mva(float(network["baseMVA"]))

    n_sub = len(bus_keys)
    if n_busbar_per_sub is None:
        n_busbar_per_sub = 1
    try:
        tmp = int(n_busbar_per_sub)
    except ValueError as exc_:
        raise RuntimeError("Impossible to convert n_busbar_per_sub to int") from exc_
    if tmp != n_busbar_per_sub:
        raise RuntimeError(f"n_busbar_per_sub should be a int, you provided {tmp} which cannot safely be converted to an int.")
    n_busbar_per_sub = tmp
    if n_busbar_per_sub <= 0:
        raise RuntimeError(f"You need to provide a grid with at least 1 busbar per "
                           f"substation / voltage level, provided n_busbar_per_sub={n_busbar_per_sub}")

    # lightsim2grid lays out global bus ids busbar-section-major: ids [0, n_sub) are
    # busbar section 1 of every substation (one per PowerModels bus, in `bus_keys`
    # order), ids [n_sub, 2*n_sub) are section 2, etc.
    base_kv = np.array([network["bus"][k].get("base_kv", 1.0) for k in bus_keys])
    vn_kv = np.tile(base_kv, n_busbar_per_sub)
    model.init_bus(n_sub, n_busbar_per_sub, vn_kv, nb_line, nb_trafo)
    if n_busbar_per_sub > 1:
        # extra busbar sections have no corresponding PowerModels bus at all
        ls_to_orig = np.concatenate((ls_to_orig,
                                     np.full(n_sub * (n_busbar_per_sub - 1), -1, dtype=int)))
    model._ls_to_orig = ls_to_orig

    if n_busbar_per_sub > 1:
        # nothing is connected to these extra busbar sections in the base case, so
        # deactivate them until a topology action moves something there
        for ls_bus_id in range(n_sub, n_sub * n_busbar_per_sub):
            model.deactivate_bus(ls_bus_id)

    isolated_keys = [k for k in bus_keys if network["bus"][k].get("bus_type") == NONE]
    if isolated_keys:
        isolated_bus_i = np.array([int(network["bus"][k]["bus_i"]) for k in isolated_keys])
        isolated_ls_bus = pm_bus_to_ls(isolated_bus_i, pm_to_ls)
        for ls_bus_id in isolated_ls_bus:
            model.deactivate_bus(int(ls_bus_id))
    else:
        isolated_ls_bus = np.array([], dtype=int)

    # init the powerlines and transformers
    _aux_add_branch(model, network, pm_to_ls, isolated_ls_bus)

    # init the shunts
    _aux_add_shunt(model, network, pm_to_ls, isolated_ls_bus)

    # init the loads
    _aux_add_load(model, network, pm_to_ls, isolated_ls_bus)

    # init the storage units
    _aux_add_storage(model, network, pm_to_ls, isolated_ls_bus)

    # init the generators
    _aux_add_gen(model, network, pm_to_ls, isolated_ls_bus)

    # deal with the slack bus(es)
    _aux_add_slack(model, network, pm_to_ls, isolated_ls_bus)

    # init the HVDC lines, if any
    _aux_add_dc_line(model, network, pm_to_ls, isolated_ls_bus)

    return model
