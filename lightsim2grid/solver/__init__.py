# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# Deprecated: use lightsim2grid.algorithm instead.
# This shim re-exports everything from lightsim2grid.algorithm for backward compatibility
# and will be removed in a future release.

from enum import Enum

import warnings as _warnings
_warnings.warn(
    "lightsim2grid.solver is deprecated and will be removed in a future release. "
    "Use lightsim2grid.algorithm instead.",
    DeprecationWarning,
    stacklevel=2,
)

# import lightsim2grid
from lightsim2grid.algorithm import *  # noqa: E402, F401, F403
from lightsim2grid.algorithm import __all__  # noqa: E402, F401


class SolverType(Enum):
    """For backward compatibility, please use 
    `from lightsim2grid.algorithm import AlgorithmType` instead now.
    """
    GaussSeidel = AlgorithmType.GaussSeidel  # noqa: F405
    GaussSeidelSynch = AlgorithmType.GaussSeidelSynch  # noqa: F405
    
    SparseLU = AlgorithmType.NR_SparseLU  # noqa: F405
    SparseLUSingleSlack = AlgorithmType.NRSing_SparseLU  # noqa: F405
    DC = AlgorithmType.DC_SparseLU  # noqa: F405
    FDPF_XB_SparseLU = AlgorithmType.FDPF_XB_SparseLU  # noqa: F405
    FDPF_BX_SparseLU = AlgorithmType.FDPF_BX_SparseLU  # noqa: F405
    
    KLU = AlgorithmType.NR_KLU  # noqa: F405
    KLUSingleSlack = AlgorithmType.NRSing_KLU  # noqa: F405
    KLUDC = AlgorithmType.DC_KLU  # noqa: F405
    FDPF_XB_KLU = AlgorithmType.FDPF_XB_KLU  # noqa: F405
    FDPF_BX_KLU = AlgorithmType.FDPF_BX_KLU  # noqa: F405
    
    NICSLU = AlgorithmType.NR_NICSLU  # noqa: F405
    NICSLUSingleSlack = AlgorithmType.NR_NICSLU  # noqa: F405
    NICSLUDC = AlgorithmType.NR_NICSLU  # noqa: F405
    FDPF_XB_NICSLU = AlgorithmType.NR_NICSLU  # noqa: F405
    FDPF_BX_NICSLU = AlgorithmType.NR_NICSLU  # noqa: F405
    
    CKTSO = AlgorithmType.NR_CKTSO  # noqa: F405
    CKTSOSingleSlack = AlgorithmType.NR_CKTSO  # noqa: F405
    CKTSODC = AlgorithmType.NR_CKTSO  # noqa: F405
    FDPF_XB_CKTSO = AlgorithmType.NR_CKTSO  # noqa: F405
    FDPF_BX_CKTSO = AlgorithmType.NR_CKTSO  # noqa: F405
    
# also for backward compatibility, some alias of the old "solver" 
# names
GaussSeidelSolver = GaussSeidelAlgo  # noqa: F405
GaussSeidelSynchSolver = GaussSeidelSynchAlgo  # noqa: F405

SparseLUSolver = NR_SparseLU  # noqa: F405
SparseLUSolverSingleSlack = NRSing_SparseLU  # noqa: F405
DCSolver = DC_SparseLU  # noqa: F405
FDPF_XB_SparseLUSolver = FDPF_XB_SparseLU  # noqa: F405
FDPF_BX_SparseLUSolver = FDPF_BX_SparseLU  # noqa: F405


try:
   KLUSolver = NR_KLU  # noqa: F405
   KLUSolverSingleSlack = NRSing_KLU  # noqa: F405 
   KLUDCSolver = DC_KLU  # noqa: F405 
   FDPF_XB_KLUSolver = FDPF_XB_KLU  # noqa: F405 
   FDPF_BX_KLUSolver = FDPF_BX_KLU  # noqa: F405 
except Exception as exc_:  # noqa: F841
    # KLU is not available
    pass


try:
    NICSLUSolver = NR_NICSLU  # noqa: F405
    NICSLUSolverSingleSlack = NRSing_NICSLU  # noqa: F405
    NICSLUDCSolver = DC_NICSLU  # noqa: F405
    FDPF_XB_NICSLUSolver = FDPF_XB_NICSLU  # noqa: F405
    FDPF_BX_NICSLUSolver =  FDPF_BX_NICSLU # noqa: F405
except Exception as exc_:  # noqa: F841
    # NICSLU is not available
    pass


try:
    CKTSOSolver = NR_CKTSO  # noqa: F405
    CKTSOSolverSingleSlack = NRSing_CKTSO  # noqa: F405
    CKTSODCSolver = DC_CKTSO  # noqa: F405
    FDPF_XB_CKTSOSolver = FDPF_XB_CKTSO  # noqa: F405
    FDPF_BX_CKTSOSolver = FDPF_XB_CKTSO  # noqa: F405
except Exception as exc_:  # noqa: F841
    # NICSLU is not available
    pass
