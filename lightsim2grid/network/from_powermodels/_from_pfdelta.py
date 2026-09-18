# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Build a LSGrid from a single PFΔ (arXiv:2510.22048) dataset row: a dict wrapping a
PowerModels.jl network data dictionary (under a `"network"` key) alongside its solved
AC power-flow state (under a `"solution"` key). The network data itself is plain
PowerModels.jl, so this is a thin wrapper around `lightsim2grid.network.from_powermodels.init`.
"""

from typing import Optional, Union
import os
import json

from ...lightsim2grid_cpp import LSGrid
from .initLSGrid import init as _init_from_powermodels


def init_from_pfdelta(row: Union[dict, str, "os.PathLike"],
                      n_busbar_per_sub: Optional[int] = None,
                      ) -> LSGrid:
    """
    Convert a PFΔ dataset row into a LSGrid.

    Parameters
    ----------
    row:
        Either a PFΔ row already parsed into a dict (with a top-level `"network"`
        key), or a path (str or `os.PathLike`) to a `.json` file containing that same
        structure.
    n_busbar_per_sub:
        Passed through directly to `lightsim2grid.network.init_from_powermodels`.

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

    return _init_from_powermodels(row["network"], n_busbar_per_sub=n_busbar_per_sub)
