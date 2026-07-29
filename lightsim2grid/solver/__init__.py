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

# defined in a private module (not here) so that other, internal, non-deprecated code
# (LightSimBackend's solver_type / SolverType back-compat bridging) can use SolverType
# without importing this module and triggering the warning above.
from lightsim2grid._solver_type import SolverType  # noqa: E402, F401

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
