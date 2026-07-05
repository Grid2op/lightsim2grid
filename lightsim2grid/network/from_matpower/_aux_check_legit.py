# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np

from ._my_const import (MIN_BUS_COLS, MIN_GEN_COLS, MIN_BRANCH_COLS, MIN_DCLINE_COLS,
                        TAP, SHIFT, BUS_TYPE, NONE)


def _aux_check_legit(bus, gen, branch, dcline):
    """
    Check that the raw matpower matrices can be handled by lightsim2grid.

    Parameters
    ----------
    bus, gen, branch, dcline: numpy arrays
        The raw matpower matrices, as returned by `_parse_matpower_source.load_matpower_data`

    Returns
    -------

    """
    if bus.shape[1] < MIN_BUS_COLS:
        raise RuntimeError(f"`mpc.bus` should have at least {MIN_BUS_COLS} columns, found {bus.shape[1]}.")
    if gen.shape[0] and gen.shape[1] < MIN_GEN_COLS:
        raise RuntimeError(f"`mpc.gen` should have at least {MIN_GEN_COLS} columns, found {gen.shape[1]}.")
    if branch.shape[0] and branch.shape[1] < MIN_BRANCH_COLS:
        raise RuntimeError(f"`mpc.branch` should have at least {MIN_BRANCH_COLS} columns, found {branch.shape[1]}.")
    if dcline.shape[0] and dcline.shape[1] < MIN_DCLINE_COLS:
        raise RuntimeError(f"`mpc.dcline` should have at least {MIN_DCLINE_COLS} columns, found {dcline.shape[1]}.")

    if branch.shape[0]:
        is_line = branch[:, TAP] == 0.
        weird_shift = is_line & (branch[:, SHIFT] != 0.)
        if np.any(weird_shift):
            warnings.warn(f"{weird_shift.sum()} branch(es) have `TAP == 0` (so are treated as a plain "
                          f"powerline by lightsim2grid) but a non-zero `SHIFT`, which is not physically "
                          f"meaningful for a plain line. The `SHIFT` will be ignored for these branches.")

    if bus.shape[0] and np.any(bus[:, BUS_TYPE] == NONE):
        n_isolated = int((bus[:, BUS_TYPE] == NONE).sum())
        warnings.warn(f"{n_isolated} bus(es) have `BUS_TYPE == {NONE}` (isolated). They will be deactivated, "
                      f"along with any load, shunt, generator or branch side still connected to them.")
