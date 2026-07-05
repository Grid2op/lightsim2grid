# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._my_const import (DC_F_BUS, DC_T_BUS, DC_BR_STATUS, DC_PF, DC_VF, DC_VT,
                        DC_QMINF, DC_QMAXF, DC_QMINT, DC_QMAXT, DC_LOSS0, DC_LOSS1)
from ._mp_bus_to_ls_bus import mp_bus_to_ls


def _aux_add_dc_line(model, dcline, mp_to_ls, isolated_ls_bus):
    """
    Add the HVDC lines described by `mpc.dcline` into the lightsim2grid "model", using
    the same simplified (station-loss-free, resistance-free) dc line model already used
    by lightsim2grid's pandapower loader (`model.init_dclines`).

    Note
    ----
    Matpower's `PF` is the MW flow "from" -> "to" measured at the `F_BUS` end (positive
    means the `F_BUS` side sends power into the line), which is the same convention as
    pandapower's dcline `p_mw`. `model.init_dclines` instead expects the power *received*
    at side 1 (the `F_BUS` side) when positive, so `PF` must be negated, exactly as
    lightsim2grid's pandapower loader already does for pandapower's `p_mw`.

    Matpower's loss model is `PT = PF - (LOSS0 + LOSS1 * PF)` with `LOSS1` a fraction of
    `PF`, while lightsim2grid/pandapower's `loss_percent` expresses that same coefficient
    as a percentage (`LOSS1 * 100`) and `loss_mw` is `LOSS0` directly.

    `PMIN` / `PMAX` (limits on `PF`) have no equivalent in `model.init_dclines` (this
    simplified model does not enforce them), consistent with how the pandapower loader
    already ignores pandapower's own dcline power limits.

    Parameters
    ----------
    model
    dcline: numpy array
        raw `mpc.dcline` matrix (can have 0 rows: most matpower cases have no HVDC line)
    mp_to_ls: dict
        matpower bus number -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`BUS_TYPE == 4`) buses

    """
    if dcline.shape[0] == 0:
        return

    from_bus = mp_bus_to_ls(dcline[:, DC_F_BUS], mp_to_ls)
    to_bus = mp_bus_to_ls(dcline[:, DC_T_BUS], mp_to_ls)

    p_mw = -dcline[:, DC_PF]
    loss_percent = 100. * dcline[:, DC_LOSS1]
    loss_mw = dcline[:, DC_LOSS0]

    model.init_dclines(from_bus, to_bus, p_mw, loss_percent, loss_mw,
                       dcline[:, DC_VF], dcline[:, DC_VT],
                       dcline[:, DC_QMINF], dcline[:, DC_QMAXF],
                       dcline[:, DC_QMINT], dcline[:, DC_QMAXT])

    dc_status = dcline[:, DC_BR_STATUS] != 0
    if isolated_ls_bus.size:
        dc_status &= ~np.isin(from_bus, isolated_ls_bus)
        dc_status &= ~np.isin(to_bus, isolated_ls_bus)
    for dc_id, is_ok in enumerate(dc_status):
        if not is_ok:
            model.deactivate_dcline(dc_id)
