# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Natively parse a MATPOWER case (".m" file, ".mat" file, or already-parsed mpc dict /
object) and initialize a LSGrid c++ object from it, without ever going through a
pandapower or pypowsybl network.
"""

from typing import Optional, Union
import os
import numpy as np

from ...lightsim2grid_cpp import LSGrid
from ._parse_matpower_source import load_matpower_data
from ._mp_bus_to_ls_bus import build_bus_remap, mp_bus_to_ls
from ._aux_check_legit import _aux_check_legit
from ._aux_add_branch import get_branch_split, _aux_add_branch
from ._aux_add_load import _aux_add_load
from ._aux_add_shunt import _aux_add_shunt
from ._aux_add_gen import _aux_add_gen
from ._aux_add_slack import _aux_add_slack
from ._my_const import BUS_I, BUS_TYPE, BASE_KV, NONE


def init(source: Union[str, "os.PathLike", dict],
         n_sub: Optional[int] = None,  # number of voltage levels
         n_busbar_per_sub: Optional[int] = None,  # max number of buses allowed per substation / voltage level
         ) -> LSGrid:
    """
    Convert a MATPOWER case into a LSGrid.

    Unlike `lightsim2grid.network.init_from_pandapower`, this never constructs a
    pandapower network: it reads MATPOWER's raw `bus` / `gen` / `branch` matrices
    directly and initializes the `LSGrid` from them. In particular, several
    generators connected to the same bus are all kept as independent generators
    (no aggregation).

    This can fail to convert the grid and still not throw any error, use with care
    (for example, you can run a powerflow after this conversion and compare the
    results against another tool, e.g. pandapower or pypowsybl, on the same case).

    Cases for which conversion is not possible include, but are not limited to:

    - `mpc.gencost`, `mpc.areas` and `mpc.dcline` are ignored (not needed for a powerflow)
    - matpower's `.m` files require the optional `matpowercaseframes` package,
      `.mat` files require the optional `scipy` package.

    Parameters
    ----------
    source:
        Either a path (str or `os.PathLike`) to a ".m" or ".mat" matpower case file,
        or an already parsed matpower case: a dict with "bus" / "gen" / "branch" /
        "baseMVA" keys (e.g. as returned by `pypower`'s `caseN()` functions), or any
        object exposing these as attributes (e.g. a `matpowercaseframes.CaseFrames`
        instance built beforehand).
    n_sub:
        number of voltage levels / substations. If not provided, defaults to one
        substation per matpower bus (`n_busbar_per_sub` is then forced to 1).
    n_busbar_per_sub:
        max number of buses allowed per substation / voltage level.

    Returns
    -------
    model: :class:`lightsim2grid.network.LSGrid`
        The initialized network

    """
    bus, gen, branch, baseMVA = load_matpower_data(source)
    _aux_check_legit(bus, gen, branch)

    mp_to_ls, ls_to_orig = build_bus_remap(bus[:, BUS_I])

    is_trafo = get_branch_split(branch)
    nb_line = int((~is_trafo).sum())
    nb_trafo = int(is_trafo.sum())

    model = LSGrid()
    model.set_sn_mva(baseMVA)

    if n_sub is None:
        n_sub = bus.shape[0]
        if n_busbar_per_sub is not None and n_busbar_per_sub != 1:
            raise RuntimeError(f"If n_sub is None, n_busbar_per_sub must be None (or 1), found {n_busbar_per_sub}.")
        n_busbar_per_sub = 1

    try:
        tmp = int(n_sub)
    except ValueError as exc_:
        raise RuntimeError("Impossible to convert n_sub to int") from exc_
    if tmp != n_sub:
        raise RuntimeError(f"n_sub should be a int, you provided {tmp} which cannot safely be converted to an int.")
    n_sub = tmp
    if n_sub <= 0:
        raise RuntimeError(f"You need to provide a grid with at least 1 substation / voltage level, provided n_sub={n_sub}")

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

    model.init_bus(n_sub, n_busbar_per_sub, bus[:, BASE_KV], nb_line, nb_trafo)
    model._ls_to_orig = ls_to_orig

    isolated_mask = bus[:, BUS_TYPE] == NONE
    if np.any(isolated_mask):
        isolated_ls_bus = mp_bus_to_ls(bus[isolated_mask, BUS_I], mp_to_ls)
        for ls_bus_id in isolated_ls_bus:
            model.deactivate_bus(int(ls_bus_id))
    else:
        isolated_ls_bus = np.array([], dtype=int)

    # init the powerlines and transformers
    _aux_add_branch(model, branch, is_trafo, mp_to_ls, isolated_ls_bus)

    # init the shunts
    _aux_add_shunt(model, bus, mp_to_ls, isolated_ls_bus)

    # init the loads
    _aux_add_load(model, bus, mp_to_ls, isolated_ls_bus)

    # init the generators
    _aux_add_gen(model, gen, mp_to_ls, isolated_ls_bus)

    # deal with the slack bus(es)
    _aux_add_slack(model, bus, gen, mp_to_ls, isolated_ls_bus)

    return model
