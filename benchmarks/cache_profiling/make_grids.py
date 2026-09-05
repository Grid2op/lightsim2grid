#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Dump a few pandapower reference grids as lightsim2grid binary (``.lsb``) files, so
that the C++ profiling driver (``profile_cached_pf.cpp``) can load them without
needing python, pandapower or the conversion code in the profile.

    python make_grids.py [output_dir]
"""

import os
import sys

import pandapower.networks as pn

from lightsim2grid.network import init_from_pandapower


# (file stem, pandapower factory).  Small + large, as asked: one IEEE case that
# fits in cache and one large pegase case that does not.
GRIDS = [
    ("case30", pn.case30),
    ("case118", pn.case118),
    ("case1354pegase", pn.case1354pegase),
    ("case9241pegase", pn.case9241pegase),
]


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for name, factory in GRIDS:
        net = factory()
        grid = init_from_pandapower(net)
        path = os.path.join(out_dir, f"{name}.lsb")
        grid.save_binary(path)
        print(f"{name}: {grid.total_bus()} buses, "
              f"{len(net.line)} lines, {len(net.trafo)} trafos, "
              f"{len(net.load)} loads -> {path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "grids")
