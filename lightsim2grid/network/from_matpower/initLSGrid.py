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

MATPOWER's raw matrix format is a strict subset of PowerModels.jl's network data
dictionary (no separate load/shunt/storage/dcline tables to begin with, a single
combined line-charging value, a "TAP == 0" sentinel for plain lines, angles in
degrees, ...), so this loader's only real job is to translate the raw matrices into
that richer dict and hand off to the shared, feature-complete
`lightsim2grid.network.from_powermodels` engine -- which is also what
`init_from_pfdelta` (PFΔ dataset rows) uses directly, since those are already
PowerModels dicts.
"""

from typing import Optional, Union
import os

from ...lightsim2grid_cpp import LSGrid
from ..from_powermodels import init as _init_from_powermodels
from ._parse_matpower_source import load_matpower_data
from ._aux_check_legit import _aux_check_legit
from ._mpc_to_powermodels import mpc_to_powermodels


def init(source: Union[str, "os.PathLike", dict],
         n_busbar_per_sub: Optional[int] = None,  # max number of buses allowed per substation / voltage level
         ) -> LSGrid:
    """
    Convert a MATPOWER case into a LSGrid.

    Unlike `lightsim2grid.network.init_from_pandapower`, this never constructs a
    pandapower network: it reads MATPOWER's raw `bus` / `gen` / `branch` / `dcline`
    matrices directly and initializes the `LSGrid` from them. In particular, several
    generators connected to the same bus are all kept as independent generators
    (no aggregation).

    This can fail to convert the grid and still not throw any error, use with care
    (for example, you can run a powerflow after this conversion and compare the
    results against another tool, e.g. pandapower or pypowsybl, on the same case).

    Cases for which conversion is not possible include, but are not limited to:

    - `mpc.gencost` and `mpc.areas` are ignored (not needed for a powerflow)
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
    n_busbar_per_sub:
        There is always exactly one substation / voltage level per matpower bus (matpower
        has no notion of several busbar sections within a bus, so this is not configurable).
        This parameter only controls how many buses / busbar sections lightsim2grid
        allocates *per substation*, which is useful if you intend to perform grid2op-like
        topology actions on the resulting grid afterwards. Defaults to 1 (no extra busbar
        section). Any extra busbar section is deactivated, since nothing in the base
        matpower case is ever connected to it.

    Returns
    -------
    model: :class:`lightsim2grid.network.LSGrid`
        The initialized network

    """
    bus, gen, branch, dcline, baseMVA = load_matpower_data(source)
    _aux_check_legit(bus, gen, branch, dcline)

    network = mpc_to_powermodels(bus, gen, branch, dcline, baseMVA)
    return _init_from_powermodels(network, n_busbar_per_sub=n_busbar_per_sub)
