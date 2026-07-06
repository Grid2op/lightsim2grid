# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from ._my_const import MIN_BUS_COLS, MIN_GEN_COLS, MIN_BRANCH_COLS, MIN_DCLINE_COLS


def _aux_check_legit(bus, gen, branch, dcline):
    """
    Check that the raw matpower matrices have enough columns to be interpreted at all.

    Semantic checks that don't depend on the raw-array representation (isolated buses,
    "weird" tap/shift combinations) are done once the matrices are converted to a
    PowerModels-style dict, so they also cover the `from_powermodels`/
    `init_from_pfdelta` entry points -- see
    `lightsim2grid.network.from_matpower._mpc_to_powermodels` (the `TAP == 0` / `SHIFT
    != 0` warning) and `lightsim2grid.network.from_powermodels._aux_check_legit` (the
    isolated-bus warning).

    Parameters
    ----------
    bus, gen, branch, dcline: numpy arrays
        The raw matpower matrices, as returned by `_parse_matpower_source.load_matpower_data`

    """
    if bus.shape[1] < MIN_BUS_COLS:
        raise RuntimeError(f"`mpc.bus` should have at least {MIN_BUS_COLS} columns, found {bus.shape[1]}.")
    if gen.shape[0] and gen.shape[1] < MIN_GEN_COLS:
        raise RuntimeError(f"`mpc.gen` should have at least {MIN_GEN_COLS} columns, found {gen.shape[1]}.")
    if branch.shape[0] and branch.shape[1] < MIN_BRANCH_COLS:
        raise RuntimeError(f"`mpc.branch` should have at least {MIN_BRANCH_COLS} columns, found {branch.shape[1]}.")
    if dcline.shape[0] and dcline.shape[1] < MIN_DCLINE_COLS:
        raise RuntimeError(f"`mpc.dcline` should have at least {MIN_DCLINE_COLS} columns, found {dcline.shape[1]}.")
