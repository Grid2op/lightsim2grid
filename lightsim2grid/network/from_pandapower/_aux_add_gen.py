# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
from ._pp_bus_to_ls_bus import pp_bus_to_ls


def _aux_add_gen(model, pp_net, pp_to_ls):
    """
    Add the generators of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    model
    pp_net

    """
    if "parallel" in pp_net.gen and np.any(pp_net.gen["parallel"].to_numpy() != 1):
        raise RuntimeError("Cannot handle 'parallel' gen columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.line[\"parallel\"] != 1 it is not handled by lightsim yet.")
    model.init_generators(pp_net.gen["p_mw"].to_numpy(),
                          pp_net.gen["vm_pu"].to_numpy(),
                          pp_net.gen["min_q_mvar"].to_numpy(),
                          pp_net.gen["max_q_mvar"].to_numpy(),
                          pp_bus_to_ls(pp_net.gen["bus"].to_numpy(), pp_to_ls)
                          )
    for gen_id, is_connected in enumerate(pp_net.gen["in_service"].to_numpy()):
        if not is_connected:
            # generator is deactivated
            model.deactivate_gen(gen_id)


def _aux_add_gen_p_limits(model, pp_net):
    """
    Thread pandapower's optional ``min_p_mw`` / ``max_p_mw`` generator columns into the
    lightsim2grid model.

    They are not used by the powerflow: they are what says whether the active power a
    distributed slack ended up asking of a machine is one it could actually deliver (see
    the batch algorithms' ``compute_physical_violations``). A net that does not carry
    them, or that carries only NaN, leaves the model without any -- exactly as if this
    had never been called.

    Called AFTER ``_aux_add_slack``, which appends one generator per ext_grid when no
    generator stands on the slack bus: those have no pandapower row, hence no limit, and
    are padded with NaN so the two vectors still match the container.
    """
    nb_gen = len(model.get_generators())
    nb_pp_gen = pp_net.gen.shape[0]
    if nb_gen == 0:
        return

    def _col(name):
        if name not in pp_net.gen:
            return None
        vals = pp_net.gen[name].to_numpy().astype(np.float64)
        return vals if np.any(np.isfinite(vals)) else None

    min_p = _col("min_p_mw")
    max_p = _col("max_p_mw")
    if min_p is None and max_p is None:
        return

    def _pad(vals):
        out = np.full(nb_gen, np.nan, dtype=np.float64)
        if vals is not None:
            out[:nb_pp_gen] = vals[:nb_pp_gen]
        return out

    model.set_gen_p_limits(_pad(min_p), _pad(max_p))
