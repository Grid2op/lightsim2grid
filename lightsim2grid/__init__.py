# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__version__ = "0.13.1"

__all__ = [
    "newtonpf",
    "SolverType",
    "ErrorType",
    "solver",
    "compilation_options",
    "load_solver_plugin"]

import ctypes as _ctypes
import os as _os
import sys as _sys

# Load the C++ extension with RTLD_GLOBAL so its symbols (BaseAlgo, SolverRegistry,
# etc.) are visible to solver plugins loaded later via load_solver_plugin().
# This is the standard pattern for extension modules that support dlopen plugins
# (used by PyTorch, JAX, and others for the same reason).
if hasattr(_sys, "getdlopenflags"):
    # Unix only: load with RTLD_GLOBAL so symbols are visible to dlopen'd plugins.
    _old_dlopen_flags = _sys.getdlopenflags()
    _sys.setdlopenflags(_old_dlopen_flags | _os.RTLD_GLOBAL)
    try:
        from lightsim2grid.solver import SolverType
        from lightsim2grid.solver import ErrorType
    finally:
        _sys.setdlopenflags(_old_dlopen_flags)
else:
    # Windows: no dlopen flags; just import normally.
    from lightsim2grid.solver import SolverType
    from lightsim2grid.solver import ErrorType


def load_solver_plugin(path: str) -> None:
    """Load a shared library containing a lightsim2grid solver plugin.

    The library must contain at least one static ``SolverRegistrar`` object
    in an anonymous namespace (see ``examples/external_solver/`` for a minimal
    example).  Its constructor fires when the library is loaded, which
    registers the new solver into the C++ ``SolverRegistry`` singleton.

    After this call the new solver name is usable via::

        grid.change_solver("MySolverName")

    and will appear in ``grid.available_solver_names()``.

    Parameters
    ----------
    path:
        Absolute or relative path to the ``.so`` / ``.dll`` file.
    """
    _ctypes.CDLL(path, mode=_ctypes.RTLD_GLOBAL)

try:
    from lightsim2grid.lightSimBackend import LightSimBackend  # noqa: F401
    __all__.append("LightSimBackend")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the Backend will not be available
    pass

try:
    from lightsim2grid.physical_law_checker import PhysicalLawChecker  # noqa: F401
    __all__.append("PhysicalLawChecker")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the PhysicalLawChecker will not be available
    pass

try:
    from lightsim2grid.timeSerie import TimeSerie  # noqa: F401
    __all__.append("TimeSerie")  
    __all__.append("timeSerie")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the TimeSeries module will not be available
    pass

try:
    from lightsim2grid.contingencyAnalysis import ContingencyAnalysis  # noqa: F401
    __all__.append("contingencyAnalysis")
    __all__.append("ContingencyAnalysis")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the SecurtiyAnalysis module will not be available
    pass
    
try:
    from lightsim2grid.rewards import N1ContingencyReward  # noqa: F401
    __all__.append("rewards")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the SecurtiyAnalysis module will not be available
    pass
