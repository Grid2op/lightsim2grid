# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
from ._pp_bus_to_ls_bus import pp_bus_to_ls


def _aux_add_shunt(model, pp_net, pp_to_ls):
    """
    Add the shunts of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    model
    pp_net

    Returns
    -------

    """
    if "parallel" in pp_net.shunt and np.any(pp_net.shunt["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' sgen columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.sgen[\"parallel\"] != 1 it is not handled by lightsim yet.")

    if pp_net.shunt.shape[0] == 0:
        # nothing to do if no shunt
        return
    
    model.init_shunt(pp_net.shunt["p_mw"].values,
                     pp_net.shunt["q_mvar"].values,
                     pp_bus_to_ls(pp_net.shunt["bus"].values, pp_to_ls)
                     )
    for sh_id, is_connected in enumerate(pp_net.shunt["in_service"].values):
        if not is_connected:
            # shunt is deactivated
            model.deactivate_shunt(sh_id)
