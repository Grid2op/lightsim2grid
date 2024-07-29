# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["GridModel"]

from lightsim2grid_cpp import GridModel

try:
    from lightsim2grid.gridmodel.from_pandapower import init as init_from_pandapower
    __all__.append("init_from_pandapower")
except ImportError:
    # pandapower is not installed
    pass

try:
    from lightsim2grid.gridmodel.from_pypowsybl import init as init_from_pypowsybl
    __all__.append("init_from_pypowsybl")
except ImportError:
    # pandapower is not installed
    pass

