# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["SolverType",
           "ErrorType",
           "AnySolver",  
           "GaussSeidelSolver",
           "GaussSeidelSynchSolver",
           "SparseLUSolver",
           "SparseLUSolverSingleSlack",
           "DCSolver",
           "FDPF_XB_SparseLUSolver",
           "FDPF_BX_SparseLUSolver"]

from lightsim2grid_cpp import SolverType
from lightsim2grid_cpp import ErrorType
from lightsim2grid_cpp import AnySolver
               
from lightsim2grid_cpp import GaussSeidelSolver  # SolverType.GaussSeidel
from lightsim2grid_cpp import GaussSeidelSynchSolver  # SolverType.GaussSeidelSynch
from lightsim2grid_cpp import SparseLUSolver  # SolverType.SparseLU
from lightsim2grid_cpp import SparseLUSolverSingleSlack  # SolverType.SparseLUSingleSlack
from lightsim2grid_cpp import DCSolver  # SolverType.DC
from lightsim2grid_cpp import FDPF_XB_SparseLUSolver  # SolverType.FDPF_XB_SparseLU
from lightsim2grid_cpp import FDPF_BX_SparseLUSolver  # SolverType.FDPF_BX_SparseLU

try:
    from lightsim2grid_cpp import KLUSolver  # SolverType.KLU
    from lightsim2grid_cpp import KLUSolverSingleSlack  # SolverType.KLUSingleSlack
    from lightsim2grid_cpp import KLUDCSolver  # SolverType.KLUDC
    from lightsim2grid_cpp import FDPF_XB_KLUSolver  # SolverType.FDPF_XB_KLU
    from lightsim2grid_cpp import FDPF_BX_KLUSolver  # SolverType.FDPF_BX_KLU
    __all__.append("KLUSolver")
    __all__.append("KLUSolverSingleSlack")
    __all__.append("KLUDCSolver")
    __all__.append("FDPF_XB_KLUSolver")
    __all__.append("FDPF_BX_KLUSolver")
except Exception as exc_:
    # KLU is not available
    pass

try:
    from lightsim2grid_cpp import NICSLUSolver  # SolverType.NICSLU
    from lightsim2grid_cpp import NICSLUSolverSingleSlack  # SolverType.NICSLUSingleSlack
    from lightsim2grid_cpp import NICSLUDCSolver  # SolverType.NICSLUDC
    from lightsim2grid_cpp import FDPF_XB_NICSLUSolver  # SolverType.FDPF_XB_NICSLU
    from lightsim2grid_cpp import FDPF_BX_NICSLUSolver  # SolverType.FDPF_BX_NICSLU
    __all__.append("NICSLUSolver")
    __all__.append("NICSLUSolverSingleSlack")
    __all__.append("NICSLUDCSolver")
    __all__.append("FDPF_XB_NICSLUSolver")
    __all__.append("FDPF_BX_NICSLUSolver")
except Exception as exc_:
    # NICSLU is not available
    pass

try:
    from lightsim2grid_cpp import CKTSOSolver  # SolverType.CKTSO
    from lightsim2grid_cpp import CKTSOSolverSingleSlack  # SolverType.CKTSOSingleSlack
    from lightsim2grid_cpp import CKTSODCSolver  # SolverType.CKTSODC
    from lightsim2grid_cpp import FDPF_XB_CKTSOSolver  # SolverType.FDPF_XB_CKTSO
    from lightsim2grid_cpp import FDPF_BX_CKTSOSolver  # SolverType.FDPF_BX_CKTSO
    __all__.append("CKTSOSolver")
    __all__.append("CKTSOSolverSingleSlack")
    __all__.append("CKTSODCSolver")
    __all__.append("FDPF_XB_CKTSOSolver")
    __all__.append("FDPF_BX_CKTSOSolver")
except Exception as exc_:
    # NICSLU is not available
    pass
