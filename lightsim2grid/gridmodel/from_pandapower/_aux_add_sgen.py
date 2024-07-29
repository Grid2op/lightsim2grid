# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
from ._pp_bus_to_ls_bus import pp_bus_to_ls

SOME_KIND_OF_INF_FOR_PMIN_PMAX = 99999.


def _aux_add_sgen(model, pp_net, pp_to_ls):
    """
    Add the static generators  (=PQ generators) of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    model
    pp_net

    Returns
    -------

    """
    if "parallel" in pp_net.sgen and np.any(pp_net.sgen["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' sgen columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.sgen[\"parallel\"] != 1 it is not handled by lightsim yet.")

    if pp_net.sgen.shape[0] == 0:
        # nothing to do if no static generators
        return

    if "min_p_mw" in pp_net.sgen:
        min_p_mw = pp_net.sgen["min_p_mw"].values
    else:
        min_p_mw = np.zeros(pp_net.sgen.shape[0]) - SOME_KIND_OF_INF_FOR_PMIN_PMAX
    if np.any(~np.isfinite(min_p_mw)):
        warnings.warn("There were some Nan in the pp_net.sgen[\"min_p_mw\"], they have been replaced by 0")
    min_p_mw[~np.isfinite(min_p_mw)] = 0.

    if "max_p_mw" in pp_net.sgen:
        max_p_mw = pp_net.sgen["max_p_mw"].values
    else:
        max_p_mw = np.zeros(pp_net.sgen.shape[0]) + SOME_KIND_OF_INF_FOR_PMIN_PMAX
    if np.any(~np.isfinite(max_p_mw)):
        warnings.warn("There were some Nan in the pp_net.sgen[\"max_p_mw\"], they have been replaced by 0")
    max_p_mw[~np.isfinite(max_p_mw)] = 0.

    if "min_q_mvar" in pp_net.sgen:
        min_q_mvar = pp_net.sgen["min_q_mvar"].values
    else:
        min_q_mvar = np.zeros(pp_net.sgen.shape[0]) - SOME_KIND_OF_INF_FOR_PMIN_PMAX
    if np.any(~np.isfinite(min_q_mvar)):
        warnings.warn("There were some Nan in the pp_net.sgen[\"min_q_mvar\"], they have been replaced by 0")
    min_q_mvar[~np.isfinite(min_q_mvar)] = 0.

    if "max_q_mvar" in pp_net.sgen:
        max_q_mvar = pp_net.sgen["max_q_mvar"].values
    else:
        max_q_mvar = np.zeros(pp_net.sgen.shape[0]) + SOME_KIND_OF_INF_FOR_PMIN_PMAX
    if np.any(~np.isfinite(max_q_mvar)):
        warnings.warn("There were some Nan in the pp_net.sgen[\"max_q_mvar\"], they have been replaced by 0")
    max_q_mvar[~np.isfinite(max_q_mvar)] = 0.

    ratio = 1.0
    if "scaling" in pp_net.sgen:
        ratio = pp_net.sgen["scaling"].values
        ratio[~np.isfinite(ratio)] = 1.0
        
    model.init_sgens(pp_net.sgen["p_mw"].values * ratio,
                     pp_net.sgen["q_mvar"].values * ratio,
                     min_p_mw,
                     max_p_mw,
                     min_q_mvar,
                     max_q_mvar,
                     pp_bus_to_ls(pp_net.sgen["bus"].values, pp_to_ls)
                     )
    for sgen_id, is_connected in enumerate(pp_net.sgen["in_service"].values):
        if not is_connected:
            # load is deactivated
            model.deactivate_sgen(sgen_id)
