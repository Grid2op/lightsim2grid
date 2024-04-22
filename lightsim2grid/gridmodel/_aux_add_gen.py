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
    if "parallel" in pp_net.gen and np.any(pp_net.gen["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' gen columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.line[\"parallel\"] != 1 it is not handled by lightsim yet.")

    model.init_generators(pp_net.gen["p_mw"].values,
                          pp_net.gen["vm_pu"].values,
                          pp_net.gen["min_q_mvar"].values,
                          pp_net.gen["max_q_mvar"].values,
                          pp_bus_to_ls(pp_net.gen["bus"].values, pp_to_ls)
                          )
    for gen_id, is_connected in enumerate(pp_net.gen["in_service"].values):
        if not is_connected:
            # generator is deactivated
            model.deactivate_gen(gen_id)
