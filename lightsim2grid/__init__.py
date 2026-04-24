# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__version__ = "0.13.2.dev0"

__all__ = [
    "newtonpf",
    "SolverType",
    "ErrorType",
    "solver",
    "compilation_options",
    "load_solver_plugin",
    "get_include",
    "get_cmake_dir"]

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
    # Windows (Python >= 3.8): the DLL search path no longer includes the
    # package directory automatically.  Add it so that lightsim2grid_core.dll
    # (a dependency of lightsim2grid_cpp.pyd) is found before the import.
    if hasattr(_os, "add_dll_directory"):
        _os.add_dll_directory(_os.path.dirname(__file__))
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


def _installed_pkg_dir() -> str:
    """Return the *installed* package directory (site-packages/lightsim2grid/).

    For editable installs __file__ points to the source tree, but the compiled
    extension and installed data (headers, cmake files) live in site-packages.
    We anchor to the .so extension module, which is always in site-packages.
    """
    import importlib.util as _ilu
    spec = _ilu.find_spec("lightsim2grid.lightsim2grid_cpp")
    if spec and spec.origin:
        return _os.path.dirname(spec.origin)
    return _os.path.dirname(_os.path.abspath(__file__))


def get_include() -> str:
    """Return the path to the lightsim2grid_core public C++ headers.

    Mirrors pybind11.get_include(). Useful for building C++ extensions that
    depend on lightsim2grid_core without going through CMake find_package.
    """
    return _os.path.join(_installed_pkg_dir(), "include")


def get_cmake_dir() -> str:
    """Return the path to the lightsim2grid_core CMake config directory.

    Use this to locate the package when calling cmake::

        cmake -Dlightsim2grid_core_DIR=$(python -c "import lightsim2grid; print(lightsim2grid.get_cmake_dir())")

    Or from Python::

        import subprocess, lightsim2grid
        subprocess.run(["cmake", f"-Dlightsim2grid_core_DIR={lightsim2grid.get_cmake_dir()}", ...])
    """
    p = _os.path.join(_installed_pkg_dir(), "share", "cmake", "lightsim2grid_core")
    if not _os.path.exists(p):
        raise ImportError(
            "lightsim2grid CMake files not found. "
            "Reinstall the package: pip install lightsim2grid"
        )
    return p
