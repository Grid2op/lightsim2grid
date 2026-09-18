# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# Deprecated: use lightsim2grid.network instead.
# This shim re-exports everything from lightsim2grid.network for backward compatibility
# and will be removed in a future release.

import warnings as _warnings
_warnings.warn(
    "lightsim2grid.gridmodel is deprecated and will be removed in a future release. "
    "Use lightsim2grid.network instead.",
    DeprecationWarning,
    stacklevel=2,
)

from lightsim2grid.network import *  # noqa: E402, F401, F403
from lightsim2grid.network import __all__  # noqa: E402, F401

# Backward-compat alias: GridModel is now LSGrid
GridModel = LSGrid  # noqa: F405
