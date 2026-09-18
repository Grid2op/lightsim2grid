# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings

from ._my_const import NONE


def _aux_check_legit(network):
    """
    Check that this PowerModels network data dictionary can be handled by lightsim2grid.

    Parameters
    ----------
    network: dict
        The PowerModels network data dictionary

    """
    if network.get("multinetwork"):
        raise RuntimeError("multinetwork PowerModels dicts are not supported by this loader.")

    n_isolated = sum(1 for bus in network["bus"].values() if bus.get("bus_type") == NONE)
    if n_isolated:
        warnings.warn(f"{n_isolated} bus(es) have `bus_type == {NONE}` (isolated). They will be "
                      "deactivated, along with any load, shunt, generator, storage, dcline or "
                      "branch side still connected to them.")
