# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Build a LSGrid from a single PFΔ (arXiv:2510.22048) dataset row: a PowerModels.jl-format
dict (MATPOWER-derived, per-unit) with a solved AC power-flow state attached.

This loader translates the row's PowerModels ``network`` dict into a raw-MATPOWER-shaped
dict and delegates the actual `LSGrid` construction to
`lightsim2grid.network.from_matpower`, which already implements bus-id remapping,
line/transformer splitting, generator/load/shunt conversion and (possibly distributed)
slack-bus handling.
"""

from typing import Optional, Union
import os
import json

from ...lightsim2grid_cpp import LSGrid
from ..from_matpower import init as _init_from_matpower
from ._pf_delta_to_mpc import network_to_mpc


def init(row: Union[dict, str, "os.PathLike"],
         n_sub: Optional[int] = None,  # number of voltage levels
         n_busbar_per_sub: Optional[int] = None,  # max number of buses allowed per substation / voltage level
         ) -> LSGrid:
    """
    Convert a PFΔ dataset row into a LSGrid.

    Parameters
    ----------
    row:
        Either a PFΔ row already parsed into a dict (with a top-level ``"network"``
        key), or a path (str or `os.PathLike`) to a ``.json`` file containing that same
        structure.
    n_sub:
        number of voltage levels / substations. If not provided, defaults to one
        substation per PFΔ bus (`n_busbar_per_sub` is then forced to 1).
    n_busbar_per_sub:
        max number of buses allowed per substation / voltage level.

    Returns
    -------
    model: :class:`lightsim2grid.network.LSGrid`
        The initialized network

    """
    if isinstance(row, (str, os.PathLike)):
        with open(row, "r") as f:
            row = json.load(f)

    if "network" not in row:
        raise RuntimeError('Expected a PFΔ row (a dict with a top-level "network" key), '
                           f"got a dict with keys {sorted(row.keys())}.")

    mpc = network_to_mpc(row["network"])
    return _init_from_matpower(mpc, n_sub=n_sub, n_busbar_per_sub=n_busbar_per_sub)
