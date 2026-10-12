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

This checks that the plugin loads and registers correctly (same scope as
examples/external_algorithm/test_plugin.py and
lightsim2grid/tests/test_solver_registry.py::TestPluginLoading), and that a
batch row's voltage-control set-point reaches it (a small hand-built grid). To
exercise the Levenberg-Marquardt iteration on a real case, build a
LightSimBackend on a grid2op environment and call
grid._grid.change_algorithm("NR_LM_KLU") (or "NR_LM_SparseLU") before runpf().
"""
import platform
import pathlib

import numpy as np

import lightsim2grid
from lightsim2grid.lightsim2grid_cpp import LSGrid, AlgorithmType, InjectionSweepCPP


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
names = gm.available_algorithm_names()
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

# ------------------------------------------------------------------
# A batch row's set-point on a voltage-control group
# ------------------------------------------------------------------
# A 4-bus feeder 0-1-2-3 whose generator 1 (on bus 1) holds bus 3 remotely: a
# voltage-control group. A batch row changing its set-point (modify_gen_v) hands it
# to the solver through BaseAlgo::set_voltage_control_v_set; a plugin that does not
# forward it solves the row at the grid's own set-point, and reports it converged.
GEN_ID, REG_BUS = 1, 3


def _remote_vc_grid(algo_name):
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.0)
    grid.init_bus(4, 1, np.full(4, 138.), 0, 0)
    grid.init_powerlines(np.full(3, 0.01), np.full(3, 0.1), np.zeros(3, dtype=complex),
                         np.array([0, 1, 2]), np.array([1, 2, 3]))
    grid.init_loads(np.array([80.]), np.array([20.]), np.array([3]))
    grid.init_generators(np.array([0., 10.]), np.array([1.02, 1.01]),
                         np.full(2, -1e3), np.full(2, 1e3), np.array([0, 1]))
    grid.add_gen_slackbus(0, 1.)
    grid.set_gen_regulated_bus(GEN_ID, REG_BUS)
    grid.change_algorithm(algo_name)
    return grid


for name in ("NR_LM_SparseLU", "NR_LM_KLU"):
    if name not in names:
        continue
    grid = _remote_vc_grid(name)
    gen_v = np.tile([g.target_vm_pu for g in grid.get_generators()], (2, 1))
    gen_v[0, GEN_ID] -= 0.02
    sweep = InjectionSweepCPP(grid)
    sweep.modify_gen_p(np.tile(np.asarray(grid.get_gen_target_p(), dtype=float), (2, 1)))
    sweep.modify_gen_v(gen_v)
    sweep.compute(np.ones(grid.total_bus(), dtype=complex), 30, 1e-10)
    assert sweep.get_status() == 1, f"{name}: the batch did not converge"
    vm = np.abs(sweep.get_voltages()[:, REG_BUS])
    np.testing.assert_allclose(vm, gen_v[:, GEN_ID], rtol=0., atol=1e-8,
                               err_msg=f"{name}: bus {REG_BUS} is not held at each row's set-point")
    print(f"{name}: each batch row holds bus {REG_BUS} at its own set-point.")

print("All checks passed.")
