#!/usr/bin/env python3
"""Smoke-test for the DummyExternal solver plugin.

Build the plugin first:
    cd examples/external_solver
    mkdir build && cd build
    cmake ..
    make

Then run:
    python test_plugin.py
"""
import os
import sys

# Make sure the installed lightsim2grid package is importable.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import GridModel, SolverType

# ------------------------------------------------------------------
# Load the plugin
# ------------------------------------------------------------------
plugin_path = os.path.join(os.path.dirname(__file__), "build", "libdummy_solver.so")
if not os.path.exists(plugin_path):
    print(f"Plugin not found at {plugin_path}. Build it first (see CMakeLists.txt).")
    sys.exit(1)

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
