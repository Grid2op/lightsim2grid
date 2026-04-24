#!/usr/bin/env python3
"""Smoke-test for the DummyExternal solver plugin.

Build the plugin first:
    cd examples/external_solver
    mkdir build && cd build
    cmake ..
    make  (or: cmake --build . --config Release on Windows)

Then run:
    python test_plugin.py
"""
import platform
import pathlib

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import GridModel, SolverType


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
lightsim2grid.load_solver_plugin(plugin_path)
print("Plugin loaded successfully.")

# ------------------------------------------------------------------
# Verify registration
# ------------------------------------------------------------------
gm = GridModel()
names = gm.available_solver_names()
assert "DummyExternal" in names, f"DummyExternal not in {names}"
print(f"Registered solvers: {sorted(names)}")

# ------------------------------------------------------------------
# Change to the plugin solver
# ------------------------------------------------------------------
gm.change_solver("DummyExternal")
assert gm.get_solver_type() == SolverType.Custom, \
    f"Expected SolverType.Custom, got {gm.get_solver_type()}"
print("change_solver('DummyExternal') OK — solver type is Custom as expected.")

print("All checks passed.")
