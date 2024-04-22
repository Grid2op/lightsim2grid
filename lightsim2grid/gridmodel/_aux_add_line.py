# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ._pp_bus_to_ls_bus import pp_bus_to_ls

def _aux_add_line(converter, model, pp_net, pp_to_ls=None):
    """
    Add the lines of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    converter
    model
    pp_net

    Returns
    -------

    """
    if "parallel" in pp_net.line and np.any(pp_net.line["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' lines columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.line[\"parallel\"] != 1 it is not handled by lightsim yet.")

    #### find the right powerline parameters
    line_r, line_x, line_h = \
        converter.get_line_param(
            pp_net.line["r_ohm_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["x_ohm_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["c_nf_per_km"].values * pp_net.line["length_km"].values,
            pp_net.line["g_us_per_km"].values * pp_net.line["length_km"].values,
            pp_net.bus.loc[pp_net.line["from_bus"]]["vn_kv"], 
            pp_net.bus.loc[pp_net.line["to_bus"]]["vn_kv"], 
            )

    ### add them to the grid
    model.init_powerlines(line_r, line_x, line_h,
                          pp_bus_to_ls(pp_net.line["from_bus"].values, pp_to_ls),
                          pp_bus_to_ls(pp_net.line["to_bus"].values, pp_to_ls)
                          )
    for line_id, is_connected in enumerate(pp_net.line["in_service"].values):
        if not is_connected:
            # powerline is deactivated
            model.deactivate_powerline(line_id)
