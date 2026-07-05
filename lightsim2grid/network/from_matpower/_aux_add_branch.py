# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._my_const import F_BUS, T_BUS, BR_R, BR_X, BR_B, TAP, SHIFT, BR_STATUS
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def get_branch_split(branch):
    """A branch row is a transformer if and only if its `TAP` column is non zero,
    this is matpower's own convention (see e.g. `makeYbus.m`). `TAP == 0` means
    "plain line, ratio 1.0".
    """
    return branch[:, TAP] != 0.


def _aux_add_branch(model, branch, is_trafo, mp_to_ls, isolated_ls_bus):
    """
    Add the lines and transformers of `mpc.branch` into the lightsim2grid "model".

    Parameters
    ----------
    model
    branch: numpy array
        raw `mpc.branch` matrix
    is_trafo: numpy bool array
        for each row of `branch`, whether it should be treated as a transformer
        (see `get_branch_split`)
    mp_to_ls: dict
        matpower bus number -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`BUS_TYPE == 4`) buses

    """
    is_line = ~is_trafo
    from_bus = mp_bus_to_ls(branch[:, F_BUS], mp_to_ls)
    to_bus = mp_bus_to_ls(branch[:, T_BUS], mp_to_ls)

    # powerlines
    line_r = branch[is_line, BR_R]
    line_x = branch[is_line, BR_X]
    line_h = 1j * branch[is_line, BR_B]
    model.init_powerlines(line_r, line_x, line_h, from_bus[is_line], to_bus[is_line])

    line_status = branch[is_line, BR_STATUS] != 0
    if isolated_ls_bus.size:
        line_status &= ~np.isin(from_bus[is_line], isolated_ls_bus)
        line_status &= ~np.isin(to_bus[is_line], isolated_ls_bus)
    for line_id, is_ok in enumerate(line_status):
        if not is_ok:
            model.deactivate_powerline(line_id)

    # transformers
    n_trafo = int(is_trafo.sum())
    trafo_r = branch[is_trafo, BR_R]
    trafo_x = branch[is_trafo, BR_X]
    trafo_b = 1j * branch[is_trafo, BR_B]
    trafo_ratio = branch[is_trafo, TAP]
    trafo_shift_degree = branch[is_trafo, SHIFT]
    trafo_tap_hv = [True] * n_trafo
    model.init_trafo(trafo_r, trafo_x, trafo_b, trafo_ratio, trafo_shift_degree, trafo_tap_hv,
                     from_bus[is_trafo], to_bus[is_trafo], True)

    trafo_status = branch[is_trafo, BR_STATUS] != 0
    if isolated_ls_bus.size:
        trafo_status &= ~np.isin(from_bus[is_trafo], isolated_ls_bus)
        trafo_status &= ~np.isin(to_bus[is_trafo], isolated_ls_bus)
    for trafo_id, is_ok in enumerate(trafo_status):
        if not is_ok:
            model.deactivate_trafo(trafo_id)
