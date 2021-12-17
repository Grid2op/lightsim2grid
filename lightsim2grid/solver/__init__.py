# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# TODO add the other solver now !
__all__ = ["SolverType",
           "ErrorType",
           "AnySolver",  
           "GaussSeidelSolver",
           "GaussSeidelSynchSolver",
           "SparseLUSolver",
           "SparseLUSolverSingleSlack",
           "DCSolver"]

from lightsim2grid_cpp import SolverType
from lightsim2grid_cpp import ErrorType
from lightsim2grid_cpp import AnySolver

from lightsim2grid_cpp import GaussSeidelSolver
from lightsim2grid_cpp import GaussSeidelSynchSolver
from lightsim2grid_cpp import SparseLUSolver
from lightsim2grid_cpp import SparseLUSolverSingleSlack
from lightsim2grid_cpp import DCSolver

try:
    from lightsim2grid_cpp import KLUSolver
    from lightsim2grid_cpp import KLUSolverSingleSlack
    from lightsim2grid_cpp import KLUDCSolver
    __all__.append("KLUSolver")
    __all__.append("KLUSolverSingleSlack")
    __all__.append("KLUDCSolver")
except Exception as exc_:
    # KLU is not available
    pass

try:
    from lightsim2grid_cpp import NICSLUSolver
    from lightsim2grid_cpp import NICSLUSolverSingleSlack
    from lightsim2grid_cpp import NICSLUDCSolver
    __all__.append("NICSLUSolver")
    __all__.append("NICSLUSolverSingleSlack")
    __all__.append("NICSLUDCSolver")
except Exception as exc_:
    # NICSLU is not available
    pass
