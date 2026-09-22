#!/usr/bin/env python3
"""Smoke-test for the NRAlgoDistSlack (distributed slack as an outer loop) plugin.

Build the plugin first:
    cd examples/dist_slack_algorithm
    mkdir build && cd build
    cmake ..
    make  (or: cmake --build . --config Release on Windows)

Then run:
    python test_plugin.py

This only checks that the plugin loads and registers correctly (same scope as
examples/external_algorithm/test_plugin.py and examples/lm_algorithm/test_plugin.py)
-- it does NOT run a real power flow. To exercise it on a real case, load a grid
and call grid.change_algorithm("NRDistSlack_KLU") (or "NRDistSlack_SparseLU")
before ac_pf().
"""
import platform
import pathlib

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType


def find_plugin():
    build = pathlib.Path(__file__).parent / "build"
    if platform.system() == "Windows":
        candidates = [
            build / "Release" / "dist_slack_algorithm.dll",
            build / "dist_slack_algorithm.dll",
        ]
    else:
        candidates = [build / "libdist_slack_algorithm.so"]
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
names = gm.available_algorithm_names()
assert "NRDistSlack_SparseLU" in names, f"NRDistSlack_SparseLU not in {names}"
print(f"Registered solvers: {sorted(names)}")

# ------------------------------------------------------------------
# Change to the plugin solver(s)
# ------------------------------------------------------------------
for name in ("NRDistSlack_SparseLU", "NRDistSlack_KLU"):
    if name not in names:
        print(f"{name} not registered (KLU not available in this build) — skipped.")
        continue
    grid = _make_grid()
    grid.change_algorithm(name)
    assert grid.get_algo_type() == AlgorithmType.Custom, \
        f"Expected AlgorithmType.Custom, got {grid.get_algo_type()}"
    print(f"change_algorithm('{name}') OK — solver type is Custom as expected.")

print("All checks passed.")
