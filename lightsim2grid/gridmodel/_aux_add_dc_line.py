# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings

import numpy as np
from ._aux_add_sgen import SOME_KIND_OF_INF_FOR_PMIN_PMAX
from ._pp_bus_to_ls_bus import pp_bus_to_ls


def _aux_add_dc_line(model, pp_net, pp_to_ls):
    """
    Add the transformers of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    converter
    model
    pp_net

    Returns
    -------

    """
    if "parallel" in pp_net.dcline and np.any(pp_net.dcline["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' dcline columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.dcline[\"parallel\"] != 1 it is not handled by lightsim yet.")

    if pp_net.dcline.shape[0] == 0:
        # nothing to do if no dc line
        return
    
    branch_from_id = pp_bus_to_ls(pp_net.dcline["from_bus"].values, pp_to_ls)
    branch_to_id = pp_bus_to_ls(pp_net.dcline["to_bus"].values, pp_to_ls)
    p_mw = -pp_net.dcline["p_mw"].values
    if np.any(~np.isfinite(p_mw)):
        warnings.warn("Some non finite values are found for p_mw, they have been replaced by 0.")
        p_mw[~np.isfinite(p_mw)] = 0.
        
    loss_percent = pp_net.dcline["loss_percent"].values
    loss_mw = pp_net.dcline["loss_mw"].values
    vm_or_pu = pp_net.dcline["vm_from_pu"].values
    vm_ex_pu = pp_net.dcline["vm_to_pu"].values
    
    min_q_or = pp_net.dcline["min_q_from_mvar"].values
    if np.any(~np.isfinite(min_q_or)):
        min_q_or[~np.isfinite(min_q_or)] = -SOME_KIND_OF_INF_FOR_PMIN_PMAX
        
    max_q_or = pp_net.dcline["max_q_from_mvar"].values
    if np.any(~np.isfinite(max_q_or)):
        max_q_or[~np.isfinite(max_q_or)] = SOME_KIND_OF_INF_FOR_PMIN_PMAX
        
    min_q_ex = pp_net.dcline["min_q_to_mvar"].values
    if np.any(~np.isfinite(min_q_ex)):
        min_q_ex[~np.isfinite(min_q_ex)] = -SOME_KIND_OF_INF_FOR_PMIN_PMAX
        
    max_q_ex = pp_net.dcline["max_q_to_mvar"].values
    if np.any(~np.isfinite(max_q_ex)):
        max_q_ex[~np.isfinite(max_q_ex)] = +SOME_KIND_OF_INF_FOR_PMIN_PMAX
                          
    model.init_dclines(branch_from_id.astype(np.int32),
                       branch_to_id.astype(np.int32),
                       p_mw,
                       loss_percent,
                       loss_mw,
                       vm_or_pu,
                       vm_ex_pu,
                       min_q_or,
                       max_q_or,
                       min_q_ex,
                       max_q_ex
                       )
    
    for dcl_id, is_connected in enumerate(pp_net.dcline["in_service"].values):
        if not is_connected:
            # trafo is deactivated
            model.deactivate_dcline(dcl_id)
