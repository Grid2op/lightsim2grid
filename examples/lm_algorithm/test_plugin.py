#!/usr/bin/env python3
"""Smoke-test for the LMNRAlgo (Levenberg-Marquardt Newton-Raphson) plugin.

Build the plugin first:
    cd examples/lm_algorithm
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

This only checks that the plugin loads and registers correctly (same scope
as examples/external_algorithm/test_plugin.py and
lightsim2grid/tests/test_solver_registry.py::TestPluginLoading) -- it does
NOT run a real power flow (that needs a fully configured grid: buses, lines,
loads, generators). To actually exercise the Levenberg-Marquardt iteration on
a real case, build a LightSimBackend on a grid2op environment and call
grid._grid.change_algorithm("NR_LM_KLU") (or "NR_LM_SparseLU") before
runpf().
"""
import platform
import pathlib

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType


def find_plugin():
    build = pathlib.Path(__file__).parent / "build"
    if platform.system() == "Windows":
        candidates = [
            build / "Release" / "lm_algorithm.dll",
            build / "lm_algorithm.dll",
        ]
    else:
        candidates = [build / "liblm_algorithm.so"]
    for p in candidates:
        if p.exists():
            return str(p)
    raise FileNotFoundError(
        f"Plugin not found (tried {[str(c) for c in candidates]}). "
        "Build it first (see CMakeLists.txt)."
    )


def _make_grid():
    """Minimal LSGrid, enough to register/switch solvers (not to run a real pf)."""
    gm = LSGrid()
    gm.set_sn_mva(100.0)
    gm.set_init_vm_pu(1.0)
    return gm


# ------------------------------------------------------------------
# Load the plugin
# ------------------------------------------------------------------
plugin_path = find_plugin()
lightsim2grid.load_algorithm_plugin(plugin_path)
print("Plugin loaded successfully.")

# ------------------------------------------------------------------
# Verify registration
# ------------------------------------------------------------------
gm = _make_grid()
names = gm.available_solver_names()
assert "NR_LM_SparseLU" in names, f"NR_LM_SparseLU not in {names}"
print(f"Registered solvers: {sorted(names)}")

# ------------------------------------------------------------------
# Change to the plugin solver(s)
# ------------------------------------------------------------------
gm.change_algorithm("NR_LM_SparseLU")
assert gm.get_algo_type() == AlgorithmType.Custom, \
    f"Expected AlgorithmType.Custom, got {gm.get_algo_type()}"
print("change_algorithm('NR_LM_SparseLU') OK — solver type is Custom as expected.")

if "NR_LM_KLU" in names:
    gm2 = _make_grid()
    gm2.change_algorithm("NR_LM_KLU")
    assert gm2.get_algo_type() == AlgorithmType.Custom, \
        f"Expected AlgorithmType.Custom, got {gm2.get_algo_type()}"
    print("change_algorithm('NR_LM_KLU') OK — solver type is Custom as expected.")
else:
    print("NR_LM_KLU not registered (KLU not available in this build) — skipped.")

print("All checks passed.")
