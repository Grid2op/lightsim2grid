# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["SolverType", "SparseLUSolver", "GaussSeidelSolver", "GaussSeidelSynchSolver", "DCSolver"]

try:
    from lightsim2grid_cpp import KLUSolver
    __all__.append("KLUSolver")
except Exception as exc_:
    # KLU is not available
    pass

try:
    from lightsim2grid_cpp import NICSLUSolver
    __all__.append("NICSLUSolver")
except Exception as exc_:
    # NICSLU is not available
    pass

from lightsim2grid_cpp import SolverType
from lightsim2grid_cpp import SparseLUSolver
from lightsim2grid_cpp import GaussSeidelSolver
from lightsim2grid_cpp import GaussSeidelSynchSolver
from lightsim2grid_cpp import DCSolver
