# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["AlgorithmType",
           "ErrorType",
           "AlgorithmSelector",
           "GaussSeidelAlgo",
           "GaussSeidelSynchAlgo",
           "NR_SparseLU",
           "NRSing_SparseLU",
           "DC_SparseLU",
           "FDPF_XB_SparseLU",
           "FDPF_BX_SparseLU"]

from ..lightsim2grid_cpp import AlgorithmType  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import ErrorType  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import AlgorithmSelector  # pyright: ignore[reportMissingImports]

from ..lightsim2grid_cpp import GaussSeidelAlgo  # AlgorithmType.GaussSeidel  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import GaussSeidelSynchAlgo  # AlgorithmType.GaussSeidelSynch  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import NR_SparseLU  # AlgorithmType.NR_SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import NRSing_SparseLU  # AlgorithmType.NRSing_SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import DC_SparseLU  # AlgorithmType.DC_SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import FDPF_XB_SparseLU  # AlgorithmType.FDPF_XB_SparseLU  # pyright: ignore[reportMissingImports]
from ..lightsim2grid_cpp import FDPF_BX_SparseLU  # AlgorithmType.FDPF_BX_SparseLU  # pyright: ignore[reportMissingImports]

try:
    from ..lightsim2grid_cpp import NR_KLU  # AlgorithmType.NR_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import NRSing_KLU  # AlgorithmType.NRSing_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import DC_KLU  # AlgorithmType.DC_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_KLU  # AlgorithmType.FDPF_XB_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_KLU  # AlgorithmType.FDPF_BX_KLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("NR_KLU")
    __all__.append("NRSing_KLU")
    __all__.append("DC_KLU")
    __all__.append("FDPF_XB_KLU")
    __all__.append("FDPF_BX_KLU")
except Exception as exc_:  # noqa: F841
    # KLU is not available
    pass

try:
    from ..lightsim2grid_cpp import NR_NICSLU  # AlgorithmType.NR_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import NRSing_NICSLU  # AlgorithmType.NRSing_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import DC_NICSLU  # AlgorithmType.DC_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_NICSLU  # AlgorithmType.FDPF_XB_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_NICSLU  # AlgorithmType.FDPF_BX_NICSLU  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("NR_NICSLU")
    __all__.append("NRSing_NICSLU")
    __all__.append("DC_NICSLU")
    __all__.append("FDPF_XB_NICSLU")
    __all__.append("FDPF_BX_NICSLU")
except Exception as exc_:  # noqa: F841
    # NICSLU is not available
    pass

try:
    from ..lightsim2grid_cpp import NR_CKTSO  # AlgorithmType.NR_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import NRSing_CKTSO  # AlgorithmType.NRSing_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import DC_CKTSO  # AlgorithmType.DC_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_XB_CKTSO  # AlgorithmType.FDPF_XB_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    from ..lightsim2grid_cpp import FDPF_BX_CKTSO  # AlgorithmType.FDPF_BX_CKTSO  # pyright: ignore[reportMissingImports]  # noqa: F401
    __all__.append("NR_CKTSO")
    __all__.append("NRSing_CKTSO")
    __all__.append("DC_CKTSO")
    __all__.append("FDPF_XB_CKTSO")
    __all__.append("FDPF_BX_CKTSO")
except Exception as exc_:  # noqa: F841
    # CKTSO is not available
    pass
