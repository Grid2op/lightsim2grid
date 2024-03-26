# Copyright (c) 2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["klu_solver_available",
           "nicslu_solver_available",
           "cktso_solver_available",
           "compiled_march_native",
           "compiled_o3_optim",
           "version",
           "lightsim2grid_lib"]

from lightsim2grid_cpp import klu_solver_available
from lightsim2grid_cpp import nicslu_solver_available
from lightsim2grid_cpp import cktso_solver_available
from lightsim2grid_cpp import compiled_march_native
from lightsim2grid_cpp import compiled_o3_optim
from lightsim2grid_cpp import version
from lightsim2grid_cpp import __file__ as lightsim2grid_lib

try:
    from lightsim2grid_cpp import nicslu_lib
    __all__.append("nicslu_lib")
except ImportError :
    # NICSLU linear solver is not available
    pass

try:
    from lightsim2grid_cpp import cktso_lib
    __all__.append("cktso_lib")
except ImportError :
    # CKTSO linear solver is not available
    pass
