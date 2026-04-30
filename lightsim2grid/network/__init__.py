# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["LSGrid"]

from .. import lightsim2grid_cpp as lightsim2grid_cpp
from ..lightsim2grid_cpp import LSGrid # type: ignore

try:
    from lightsim2grid.network.from_pandapower import init as init_from_pandapower  # noqa
    __all__.append("init_from_pandapower")
except ImportError:
    # pandapower is not installed
    pass

try:
    from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl  # noqa
    __all__.append("init_from_pypowsybl")
except ImportError:
    # pypowsybl is not installed
    pass

try:
    from lightsim2grid.network.compare_lsgrid import compare_lsgrid  # noqa
    __all__.append("compare_lsgrid")
except ImportError:
    pass

