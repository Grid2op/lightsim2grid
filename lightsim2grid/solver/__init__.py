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

from ..lightsim2grid_cpp import SolverType # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import ErrorType  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import AnySolver  # pyright: ignore[reportMissingImports]
               
from ..lightsim2grid_cpp import GaussSeidelSolver  # SolverType.GaussSeidel # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import GaussSeidelSynchSolver  # SolverType.GaussSeidelSynch  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import SparseLUSolver  # SolverType.SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import SparseLUSolverSingleSlack  # SolverType.SparseLUSingleSlack  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import DCSolver  # SolverType.DC  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import FDPF_XB_SparseLUSolver  # SolverType.FDPF_XB_SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import FDPF_BX_SparseLUSolver  # SolverType.FDPF_BX_SparseLU  # pyright: ignore[reportMissingImports]

try:
    from ..lightsim2grid_cpp import KLUSolver  # SolverType.KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import KLUSolverSingleSlack  # SolverType.KLUSingleSlack  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import KLUDCSolver  # SolverType.KLUDC  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_KLUSolver  # SolverType.FDPF_XB_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_KLUSolver  # SolverType.FDPF_BX_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("KLUSolver")
    __all__.append("KLUSolverSingleSlack")
    __all__.append("KLUDCSolver")
    __all__.append("FDPF_XB_KLUSolver")
    __all__.append("FDPF_BX_KLUSolver")
except Exception as exc_:  # noqa: F841
    # KLU is not available
    pass

try:
    from ..lightsim2grid_cpp import NICSLUSolver  # SolverType.NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import NICSLUSolverSingleSlack  # SolverType.NICSLUSingleSlack  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import NICSLUDCSolver  # SolverType.NICSLUDC  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_NICSLUSolver  # SolverType.FDPF_XB_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_NICSLUSolver  # SolverType.FDPF_BX_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("NICSLUSolver")
    __all__.append("NICSLUSolverSingleSlack")
    __all__.append("NICSLUDCSolver")
    __all__.append("FDPF_XB_NICSLUSolver")
    __all__.append("FDPF_BX_NICSLUSolver")
except Exception as exc_:  # noqa: F841
    # NICSLU is not available
    pass

try:
    from ..lightsim2grid_cpp import CKTSOSolver  # SolverType.CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import CKTSOSolverSingleSlack  # SolverType.CKTSOSingleSlack  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import CKTSODCSolver  # SolverType.CKTSODC  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_CKTSOSolver  # SolverType.FDPF_XB_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_CKTSOSolver  # SolverType.FDPF_BX_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("CKTSOSolver")
    __all__.append("CKTSOSolverSingleSlack")
    __all__.append("CKTSODCSolver")
    __all__.append("FDPF_XB_CKTSOSolver")
    __all__.append("FDPF_BX_CKTSOSolver")
except Exception as exc_:  # noqa: F841
    # NICSLU is not available
    pass
