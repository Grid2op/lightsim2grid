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
    n_busbar_per_sub:
        There is always exactly one substation / voltage level per PFΔ bus (PFΔ has no
        notion of several busbar sections within a bus, so this is not configurable).
        This parameter only controls how many buses / busbar sections lightsim2grid
        allocates *per substation*, which is useful if you intend to perform grid2op-like
        topology actions on the resulting grid afterwards. Defaults to 1 (no extra busbar
        section). Any extra busbar section is deactivated, since nothing in the base
        PFΔ row is ever connected to it. Passed through directly to
        `lightsim2grid.network.init_from_matpower`.

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
    return _init_from_matpower(mpc, n_busbar_per_sub=n_busbar_per_sub)
