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


def _aux_add_trafo(converter, model, pp_net, pp_to_ls):
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
    if "parallel" in pp_net.trafo and np.any(pp_net.trafo["parallel"].values != 1):
        raise RuntimeError("Cannot handle 'parallel' trafo columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.trafo[\"parallel\"] != 1 it is not handled by lightsim yet.")

    # fix the missing values
    tap_neutral = 1.0 * pp_net.trafo["tap_neutral"].values
    if np.any(~np.isfinite(tap_neutral)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_neutral\"], they have been replaced by 0")
    tap_neutral[~np.isfinite(tap_neutral)] = 0.

    if np.any(tap_neutral != 0.):
        raise RuntimeError("lightsim converter supposes that tap_neutral is 0 for the transformers")

    tap_step_pct = 1.0 * pp_net.trafo["tap_step_percent"].values
    if np.any(~np.isfinite(tap_step_pct)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_percent\"], they have been replaced by 0")
    tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

    tap_pos = 1.0 * pp_net.trafo["tap_pos"].values
    if np.any(~np.isfinite(tap_pos)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_pos\"], they have been replaced by 0")
    tap_pos[~np.isfinite(tap_pos)] = 0.

    shift_ = 1.0 * pp_net.trafo["shift_degree"].values
    if np.any(~np.isfinite(tap_pos)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"shift_degree\"], they have been replaced by 0")
    shift_[~np.isfinite(shift_)] = 0.

    is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
    if np.any(~np.isfinite(is_tap_hv_side)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_side\"], they have been replaced by \"hv\"")
    is_tap_hv_side[~np.isfinite(is_tap_hv_side)] = True

    if np.any(pp_net.trafo["tap_phase_shifter"].values):
        raise RuntimeError("ideal phase shifter are not modeled. Please remove all trafo with "
                           "pp_net.trafo[\"tap_phase_shifter\"] set to True.")

    tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
    if np.any(~np.isfinite(tap_angles_)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_degree\"], they have been replaced by 0")
    tap_angles_[~np.isfinite(tap_angles_)] = 0.
    tap_angles_ = np.deg2rad(tap_angles_)

    # compute physical parameters
    trafo_r, trafo_x, trafo_b = \
        converter.get_trafo_param(tap_step_pct,
                                  tap_pos,
                                  tap_angles_,  # in radian !
                                  is_tap_hv_side,
                                  pp_net.bus.loc[pp_net.trafo["hv_bus"]]["vn_kv"],
                                  pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"],
                                  pp_net.trafo["vk_percent"].values,
                                  pp_net.trafo["vkr_percent"].values,
                                  pp_net.trafo["sn_mva"].values,
                                  pp_net.trafo["pfe_kw"].values,
                                  pp_net.trafo["i0_percent"].values,
                                  )

    # initialize the grid
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_b,
                     tap_step_pct,
                     tap_pos,
                     shift_,
                     is_tap_hv_side,
                     pp_bus_to_ls(pp_net.trafo["hv_bus"].values, pp_to_ls),
                     pp_bus_to_ls(pp_net.trafo["lv_bus"].values, pp_to_ls))

    for tr_id, is_connected in enumerate(pp_net.trafo["in_service"].values):
        if not is_connected:
            # trafo is deactivated
            model.deactivate_trafo(tr_id)
