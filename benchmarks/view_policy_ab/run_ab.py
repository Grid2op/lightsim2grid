# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""A/B of the return value policy of LSGrid's result getters: ``reference`` (before 1.3.0)
against ``reference_internal`` (each array keeps the grid alive, since), ``copy`` for scale. See
README.md.

    AB_VIEWS_DIR=build_ab python run_ab.py time   # wall time per call, interleaved rounds
    AB_VIEWS_DIR=build_ab python run_ab.py ir     # instructions per call, under callgrind
    python run_ab.py _calls ...                   # internal: the loop callgrind measures

``AB_VIEWS_DIR``: the CMake build directory holding the ``ab_views`` module.
"""

import gc
import os
import random
import re
import statistics
import subprocess
import sys
import tempfile
import timeit

# before numpy loads: its BLAS worker threads spin after start-up, which callgrind counts (and
# a timing run competes with); nothing here uses BLAS
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.abspath(os.environ.get("AB_VIEWS_DIR", "build_ab")))

GETTERS = ["get_gen_target_p", "get_gen_res", "get_line_res1", "get_line_res1_full"]
POLICIES = ["reference", "reference_internal", "copy"]


def _check_lifetimes(grid):
    """what the policies are about: the reference views do not keep the grid alive"""
    for getter in GETTERS:
        ref = getattr(grid, getter + "_reference")()
        internal = getattr(grid, getter + "_reference_internal")()
        copy = getattr(grid, getter + "_copy")()
        first = (lambda res: res[0] if isinstance(res, tuple) else res)
        assert first(ref).base is None and not first(ref).flags.owndata
        assert first(internal).base is grid
        assert first(copy).flags.owndata


def _lightsim2grid_grid():
    """a solved grid of lightsim2grid's own binding (LTO, as shipped), None without it: its
    get_gen_res is the row of the policy the installed version ships, as users call it"""
    try:
        import numpy as np
        from lightsim2grid.lightsim2grid_cpp import LSGrid
    except ImportError:
        return None
    grid = LSGrid()
    grid.set_sn_mva(100.)
    grid.set_init_vm_pu(1.)
    grid.init_bus(3, 1, np.full(3, 138.), 0, 0)
    grid.init_powerlines(np.full(3, .01), np.full(3, .1), np.zeros(3, complex), np.array([0, 1, 0]), np.array([1, 2, 2]))
    grid.init_loads(np.array([10., 10.]), np.array([3., 3.]), np.array([1, 2]))
    grid.init_generators(np.full(1, 20.), np.full(1, 1.02), np.full(1, -1000.), np.full(1, 1000.), np.array([0]))
    grid.add_gen_slackbus(0, 1.)
    assert grid.ac_pf(np.ones(3, dtype=complex), 30, 1e-8).size == 3
    return grid


def time_mode(sizes=(20, 2000), number=200_000, rounds=15, seed=0):
    import ab_views
    rng = random.Random(seed)
    print(f"python {sys.version.split()[0]}, {number} calls per measure, {rounds} interleaved rounds")
    shipped = _lightsim2grid_grid()
    for n_bus in sizes:
        grid = ab_views.make_grid(n_bus)
        _check_lifetimes(grid)
        candidates = [("nop", grid.nop)]
        if shipped is not None and n_bus == sizes[0]:
            candidates.append(("shipped", shipped.get_gen_res))
        for getter in GETTERS:
            for policy in POLICIES:
                candidates.append((f"{getter}_{policy}", getattr(grid, f"{getter}_{policy}")))
        samples = {name: [] for name, _ in candidates}
        for _ in range(rounds):
            order = list(candidates)
            rng.shuffle(order)  # no variant always runs first, or last
            for name, fun in order:
                gc.collect()
                samples[name].append(timeit.Timer(fun).timeit(number) / number * 1e9)
        print(f"\n{n_bus} buses ({len(grid.get_gen_res_reference()[0])} generators, "
              f"{len(grid.get_line_res1_reference()[0])} lines): ns per call, median [min]")
        nop = statistics.median(samples["nop"])
        print(f"  {'nop (dispatch only)':<28}{nop:9.1f} [{min(samples['nop']):7.1f}]")
        if "shipped" in samples:
            print(f"  {'lightsim2grid get_gen_res':<28}{statistics.median(samples['shipped']):9.1f} "
                  f"[{min(samples['shipped']):7.1f}]   (the installed binding)")
        for getter in GETTERS:
            med = {p: statistics.median(samples[f"{getter}_{p}"]) for p in POLICIES}
            mins = {p: min(samples[f"{getter}_{p}"]) for p in POLICIES}
            delta = med["reference_internal"] - med["reference"]
            print(f"  {getter:<28}" + "".join(f"{p}: {med[p]:8.1f} [{mins[p]:7.1f}]   " for p in POLICIES)
                  + f"internal - reference: {delta:+6.1f} ns ({100. * delta / med['reference']:+5.1f} %)")


def _calls(getter_name, n_calls, n_bus):
    """the loop callgrind counts: build the grid, then call the getter n_calls times"""
    import ab_views
    grid = ab_views.make_grid(int(n_bus))
    fun = getattr(grid, getter_name)
    for _ in range(int(n_calls)):
        fun()


def _instructions(getter_name, n_calls, n_bus):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "callgrind.out")
        # a fixed hash seed: two runs then differ by their number of calls only
        subprocess.run(["valgrind", "--tool=callgrind", f"--callgrind-out-file={out}",
                        sys.executable, os.path.abspath(__file__), "_calls", getter_name,
                        str(n_calls), str(n_bus)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=dict(os.environ, PYTHONHASHSEED="0"))
        with open(out) as f:
            match = re.search(r"^summary:\s+(\d+)", f.read(), re.MULTILINE)
        return int(match.group(1))


def ir_mode(n_bus=20, n1=2_000, n2=22_000):
    """instructions per call: (Ir with n2 calls - Ir with n1 calls) / (n2 - n1), so that the
    interpreter start, the imports and the grid build cancel out (callgrind counts are
    deterministic: one run of each is enough)"""
    print(f"instructions per call (callgrind Ir), {n_bus} buses, ({n2} - {n1}) calls")
    per_call = {}
    for name in ["nop"] + [f"{g}_{p}" for g in GETTERS for p in POLICIES]:
        per_call[name] = (_instructions(name, n2, n_bus) - _instructions(name, n1, n_bus)) / (n2 - n1)
    print(f"  {'nop (dispatch only)':<28}{per_call['nop']:9.0f}")
    for getter in GETTERS:
        ref = per_call[f"{getter}_reference"]
        internal = per_call[f"{getter}_reference_internal"]
        print(f"  {getter:<28}" + "".join(f"{p}: {per_call[f'{getter}_{p}']:8.0f}   " for p in POLICIES)
              + f"internal - reference: {internal - ref:+6.0f} ({100. * (internal - ref) / ref:+5.1f} %)")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "time"
    if mode == "time":
        time_mode()
    elif mode == "ir":
        ir_mode()
    elif mode == "_calls":
        _calls(*sys.argv[2:5])
    else:
        raise SystemExit(__doc__)
