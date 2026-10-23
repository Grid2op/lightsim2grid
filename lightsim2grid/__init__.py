# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__version__ = "1.3.0.dev1"

__all__ = [
    "newtonpf",
    "AlgorithmType",
    "ErrorType",
    "algorithm",
    "solver",
    "compilation_options",
    "load_algorithm_plugin",
    "find_algorithm_plugin",
    "get_include",
    "get_cmake_dir"]

import os as _os
import sys as _sys

# Load the C++ extension with RTLD_GLOBAL so its symbols (BaseAlgo, AlgorithmRegistry,
# etc.) are visible to solver plugins loaded later via load_solver_plugin().
# This is the standard pattern for extension modules that support dlopen plugins
# (used by PyTorch, JAX, and others for the same reason).
if hasattr(_sys, "getdlopenflags"):
    # Unix only: load with RTLD_GLOBAL so symbols are visible to dlopen'd plugins.
    _old_dlopen_flags = _sys.getdlopenflags()
    _sys.setdlopenflags(_old_dlopen_flags | _os.RTLD_GLOBAL)
    try:
        from lightsim2grid.algorithm import AlgorithmType
        from lightsim2grid.algorithm import ErrorType
    finally:
        _sys.setdlopenflags(_old_dlopen_flags)
else:
    # Windows (Python >= 3.8): the DLL search path no longer includes the
    # package directory automatically.  Add it so that lightsim2grid_core.dll
    # (a dependency of lightsim2grid_cpp.pyd) is found before the import.
    if hasattr(_os, "add_dll_directory"):
        _os.add_dll_directory(_os.path.dirname(__file__))
    from lightsim2grid.algorithm import AlgorithmType
    from lightsim2grid.algorithm import ErrorType


def load_algorithm_plugin(path: str) -> None:
    """Load a shared library containing a lightsim2grid algorithm plugin.

    The library must export a ``ls2g_register_plugin`` entry point, which the
    ``LS2G_PLUGIN_ENTRY`` macro generates from the plugin's registration
    function (see ``examples/external_algorithm/`` for a minimal example).  The
    library is loaded, the entry point is looked up and called, and the new
    algorithm(s) it declares are registered into the C++ ``AlgorithmRegistry``
    singleton.

    Registration happens through a function the loader calls explicitly, *not*
    inside a static constructor during ``dlopen``, so any failure -- an ABI
    mismatch, a solver name that is already registered (which includes loading
    the same plugin twice), a library that is not a plugin, a missing file --
    is raised here as a normal, catchable Python exception instead of aborting
    the interpreter.

    After this call the new solver name is usable via::

        grid.change_algorithm("MySolverName")

    and will appear in ``grid.available_algorithm_names()``.

    Parameters
    ----------
    path:
        Absolute or relative path to the ``.so`` / ``.dll`` file, or the bare
        name of a plugin (``"rebuild_outer_algorithm"``, ``"rebuild_outer"``),
        looked up by :func:`find_algorithm_plugin`.

    Raises
    ------
    RuntimeError
        If the library cannot be loaded, does not export the plugin entry
        point, or its registration is refused (ABI mismatch, duplicate name),
        or if a bare name is not found.
    """
    from lightsim2grid.lightsim2grid_cpp import _load_algorithm_plugin
    path = str(path)
    if _is_plugin_name(path) and not _os.path.isfile(path):
        path = find_algorithm_plugin(path)
    _load_algorithm_plugin(path)


_PLUGIN_PATH_ENV = "LS2G_PLUGIN_PATH"
_PLUGIN_SUFFIXES = (".so", ".dll", ".dylib")


def _is_plugin_name(path: str) -> bool:
    """A bare name: no directory part and no shared library extension."""
    seps = [s for s in (_os.sep, _os.altsep, "/") if s]
    return not any(s in path for s in seps) and _os.path.splitext(path)[1] not in _PLUGIN_SUFFIXES


def _plugin_file_names(name: str):
    if _sys.platform == "win32":
        return [f"{name}.dll", f"lib{name}.dll"]
    # the example plugins are lib<name>.so on macOS too (see their CMakeLists)
    return [f"lib{name}.so", f"{name}.so", f"lib{name}.dylib", f"{name}.dylib"]


def _plugin_libraries_in(directory: str):
    """The shared libraries of a build folder (and of its Release sub-folder, MSVC)."""
    res = []
    for folder in (directory, _os.path.join(directory, "Release")):
        if _os.path.isdir(folder):
            res += sorted(_os.path.join(folder, f) for f in _os.listdir(folder)
                          if _os.path.splitext(f)[1] in _PLUGIN_SUFFIXES
                          and _os.path.isfile(_os.path.join(folder, f)))
    return res


def _checkout_examples_dir():
    """examples/ of the source checkout, when lightsim2grid is an editable install of it."""
    candidate = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), "examples")
    return candidate if _os.path.isdir(candidate) else None


def find_algorithm_plugin(name: str) -> str:
    """Return the path of a built algorithm plugin, given its name.

    ``name`` is the plugin's library name (``"dummy_solver"`` for
    ``libdummy_solver.so`` / ``dummy_solver.dll``) or, for the plugins of the
    ``examples/`` folder, the name of their folder; the ``_algorithm`` suffix
    can be left out (``"rebuild_outer"``). The first match is returned, looking in:

    1. the directories listed in the ``LS2G_PLUGIN_PATH`` environment variable
       (separated by ``os.pathsep``, like ``PATH``), and their ``Release`` sub-folder;
    2. the ``examples/<plugin>/build`` folders of the source checkout, when
       lightsim2grid is an editable install (``pip install -e .``) of it. A
       wheel install has no ``examples/`` folder: use ``LS2G_PLUGIN_PATH``.

    ``examples/build_plugins.py`` builds the example plugins where the second
    rule finds them.

    Raises
    ------
    RuntimeError
        If no plugin of that name is found.
    """
    names = [name] if name.endswith("_algorithm") else [name, f"{name}_algorithm"]
    file_names = [f for n in names for f in _plugin_file_names(n)]

    searched = []
    for directory in _os.environ.get(_PLUGIN_PATH_ENV, "").split(_os.pathsep):
        if not directory:
            continue
        for folder in (directory, _os.path.join(directory, "Release")):
            searched.append(folder)
            for file_name in file_names:
                cand = _os.path.join(folder, file_name)
                if _os.path.isfile(cand):
                    return _os.path.abspath(cand)

    examples = _checkout_examples_dir()
    if examples is not None:
        plugin_dirs = sorted(d for d in _os.listdir(examples)
                             if _os.path.isfile(_os.path.join(examples, d, "CMakeLists.txt")))
        # by folder name: the one library its build holds (examples/external_algorithm -> dummy_solver)
        for n in names:
            if n in plugin_dirs:
                build = _os.path.join(examples, n, "build")
                searched.append(build)
                libs = _plugin_libraries_in(build)
                if len(libs) == 1:
                    return _os.path.abspath(libs[0])
        # by library name, in every example's build
        for d in plugin_dirs:
            build = _os.path.join(examples, d, "build")
            for lib in _plugin_libraries_in(build):
                if _os.path.basename(lib) in file_names:
                    return _os.path.abspath(lib)
            searched.append(build)

    where = "\n  ".join(dict.fromkeys(searched)) or "(nowhere: LS2G_PLUGIN_PATH is not set and " \
                                                     "lightsim2grid is not an editable install)"
    raise RuntimeError(
        f"Algorithm plugin '{name}' not found (looked for {', '.join(file_names)}) in:\n  {where}\n"
        f"Build it (python examples/build_plugins.py {name}), add its folder to {_PLUGIN_PATH_ENV}, "
        f"or pass the full path of the library.")

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
    from lightsim2grid.injectionSweep import InjectionSweep  # noqa: F401
    __all__.append("InjectionSweep")
    __all__.append("injectionSweep")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the InjectionSweep module will not be available
    pass

try:
    from lightsim2grid.contingencyAnalysis import ContingencyAnalysis  # noqa: F401
    __all__.append("contingencyAnalysis")
    __all__.append("ContingencyAnalysis")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the SecurtiyAnalysis module will not be available
    pass

try:
    from lightsim2grid.scenarioSweep import ScenarioSweep  # noqa: F401
    __all__.append("ScenarioSweep")
    __all__.append("scenarioSweep")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the ScenarioSweep module will not be available
    pass

# the continuation powerflow works on an LSGrid directly, it needs no grid2op
from lightsim2grid.continuationPowerflow import ContinuationPowerFlow, run_cpf  # noqa: F401
__all__.append("ContinuationPowerFlow")
__all__.append("run_cpf")
__all__.append("continuationPowerflow")

try:
    from lightsim2grid.lightEnv import LightEnv  # noqa: F401
    __all__.append("LightEnv")
    __all__.append("lightEnv")
except ImportError as exc_:  # noqa: F841
    # grid2op is not installed, the LightEnv module will not be available
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
