#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Step 1 of the C++-side powerflow profile: MATPOWER case -> LSGrid -> ".lsb".

Converting a MATPOWER ".m" file goes through pandas / matpowercaseframes and
costs orders of magnitude more than the powerflow it is there to set up, so it
is done ONCE, here, and never inside anything that gets timed or profiled.
benchmark_binary_powerflow.py then reloads the binaries -- which is pure C++ --
and runs the powerflows.

Needs the optional ``matpowercaseframes`` package (".m" files); ".mat" files
need ``scipy`` instead. Get the cases from MATPOWER, e.g.

    for c in case118 case1354pegase case2869pegase case9241pegase; do
        curl -sSfLO https://raw.githubusercontent.com/MATPOWER/matpower/master/data/$c.m
    done

Usage:  benchmark_matpower_to_binary.py [case_dir] [out_dir]
"""

import os
import sys
import time

from lightsim2grid.network import init_from_matpower
from lightsim2grid.network import LSGrid

DEFAULT_CASES = ["case118", "case1354pegase", "case2869pegase", "case9241pegase"]


def convert(case_dir, out_dir, cases=None):
    os.makedirs(out_dir, exist_ok=True)
    done = []
    for name in (cases if cases is not None else DEFAULT_CASES):
        src = os.path.join(case_dir, name + ".m")
        if not os.path.exists(src):
            print(f"{name:16s} SKIPPED (no {src})")
            continue

        t0 = time.perf_counter()
        grid = init_from_matpower(src)
        t_parse = time.perf_counter() - t0

        dst = os.path.join(out_dir, name + ".lsb")
        t0 = time.perf_counter()
        grid.save_binary(dst)
        t_save = time.perf_counter() - t0

        # a binary nothing can read back is not a benchmark input
        LSGrid.load_binary(dst)

        print(f"{name:16s} buses={grid.total_bus():6d} lines={len(grid.get_lines()):6d} "
              f"trafos={len(grid.get_trafos()):5d} gens={len(grid.get_generators()):5d} "
              f"loads={len(grid.get_loads()):5d} | matpower {t_parse:6.2f}s  "
              f"save {t_save * 1e3:6.1f}ms  {os.path.getsize(dst) / 1024:8.1f} KiB")
        done.append(dst)
    return done


if __name__ == "__main__":
    convert(sys.argv[1] if len(sys.argv) > 1 else "cases",
            sys.argv[2] if len(sys.argv) > 2 else "grids")
