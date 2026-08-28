#!/usr/bin/env python3
"""Smoke-test for the DummyExternal solver plugin.

Build the plugin first:
    cd examples/external_algorithm
    mkdir build && cd build
    cmake ..
    make  (or: cmake --build . --config Release on Windows)

The plugin's CMakeLists.txt automatically matches whatever -march=native /
-O3 the installed lightsim2grid_core was built with (see
examples/cmake/MatchLightsim2gridBuildFlags.cmake and
docs/solver_plugin.rst) -- this is required, not just a performance nicety:
lightsim2grid_core and this plugin are separate shared libraries, and a
mismatched -march=native changes Eigen's alignment assumptions, so an Eigen
object allocated in one and freed in the other silently corrupts the heap.
Nothing to do here unless you're hand-rolling your own build system instead
of the provided CMakeLists.txt.

Then run:
    python test_plugin.py
"""
import platform
import pathlib

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType


def find_plugin():
    build = pathlib.Path(__file__).parent / "build"
    if platform.system() == "Windows":
        candidates = [
            build / "Release" / "dummy_solver.dll",
            build / "dummy_solver.dll",
        ]
    else:
        candidates = [build / "libdummy_solver.so"]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(
        f"Plugin not found (tried {[str(c) for c in candidates]}). "
        "Build it first (see CMakeLists.txt)."
    )


# ------------------------------------------------------------------
# Load the plugin
# ------------------------------------------------------------------
plugin_path = find_plugin()
lightsim2grid.load_algorithm_plugin(plugin_path)
print("Plugin loaded successfully.")

# ------------------------------------------------------------------
# Verify registration
# ------------------------------------------------------------------
gm = LSGrid()
names = gm.available_algorithm_names()
assert "DummyExternal" in names, f"DummyExternal not in {names}"
print(f"Registered solvers: {sorted(names)}")

# ------------------------------------------------------------------
# Change to the plugin solver
# ------------------------------------------------------------------
gm.change_algorithm("DummyExternal")
assert gm.get_algo_type() == AlgorithmType.Custom, \
    f"Expected AlgorithmType.Custom, got {gm.get_algo_type()}"
print("change_algorithm('DummyExternal') OK — solver type is Custom as expected.")

print("All checks passed.")
