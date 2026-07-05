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
    from lightsim2grid.network.from_pypowsybl import bake_outer_loops  # noqa
    from lightsim2grid.network.from_pypowsybl import get_pypowsybl_loopfree_parameters  # noqa
    from lightsim2grid.network.from_pypowsybl import get_pypowsybl_loopfree_distributed_slack_parameters  # noqa
    from lightsim2grid.network.from_pypowsybl import compare_baked, ComparisonResult  # noqa
    __all__.append("init_from_pypowsybl")
    __all__.append("bake_outer_loops")
    __all__.append("get_pypowsybl_loopfree_parameters")
    __all__.append("get_pypowsybl_loopfree_distributed_slack_parameters")
    __all__.append("compare_baked")
    __all__.append("ComparisonResult")
except ImportError:
    # pypowsybl is not installed
    pass

try:
    from lightsim2grid.network.from_matpower import init as init_from_matpower  # noqa
    __all__.append("init_from_matpower")
except ImportError:
    # should not happen: from_matpower has no hard dependency beyond numpy,
    # matpowercaseframes/scipy are only imported lazily when actually reading
    # a ".m"/".mat" file
    pass

from lightsim2grid.network.from_pf_delta import init as init_from_pf_delta  # noqa
__all__.append("init_from_pf_delta")

try:
    from lightsim2grid.network.compare_lsgrid import compare_lsgrid  # noqa
    __all__.append("compare_lsgrid")
except ImportError:
    pass

