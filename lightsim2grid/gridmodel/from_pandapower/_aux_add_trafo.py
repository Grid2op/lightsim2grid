# Copyright (c) 2020-2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
import pandapower as pp
from packaging import version

from ._pp_bus_to_ls_bus import pp_bus_to_ls
from ._my_const import _MIN_PP_VERSION_ADV_GRID_MODEL

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

    if "tap_phase_shifter" in pp_net.trafo:
        if np.any(pp_net.trafo["tap_phase_shifter"].values):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all trafos with "
                               "pp_net.trafo[\"tap_phase_shifter\"] set to True.")
    elif "tap_changer_type" in pp_net.trafo:
        if np.any(pp_net.trafo["tap_changer_type"].values == "Ideal"):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 2-winding trafos "
                               "with \"tap_changer_type\" set to \"Ideal\".")
    elif "tap_changer_type" in pp_net.trafo3w:
        if np.any(pp_net.trafo3w["tap_changer_type"].values == "Ideal"):
            raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 3-winding trafos "
                               "with \"tap_changer_type\" set to \"Ideal\".")

    tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
    if np.any(~np.isfinite(tap_angles_)):
        warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_degree\"], they have been replaced by 0")
    tap_angles_[~np.isfinite(tap_angles_)] = 0.
    tap_angles_ = np.deg2rad(tap_angles_)

    trafo_model_is_t = True
    if "_options" in pp_net and "trafo_model" in pp_net._options:
        trafo_model_is_t = pp_net._options["trafo_model"] == "t"
    # compute physical parameters
    if version.parse(pp.__version__) >= _MIN_PP_VERSION_ADV_GRID_MODEL:
        trafo_r, trafo_x, trafo_b = \
            converter.get_trafo_param_pp3(tap_step_pct,
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
                                          trafo_model_is_t
                                          )
    else:
        if not trafo_model_is_t:
            raise RuntimeError("Cannot convert a transformer with model 'pi' to LightSim2grid (using pandapower < 3)")
        trafo_r, trafo_x, trafo_b = \
            converter.get_trafo_param_pp2(tap_step_pct,
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
    model.init_trafo_pandapower(trafo_r,
                                trafo_x,
                                trafo_b,
                                tap_step_pct,
                                tap_pos,
                                shift_,
                                is_tap_hv_side,
                                pp_bus_to_ls(pp_net.trafo["hv_bus"].values, pp_to_ls),
                                pp_bus_to_ls(pp_net.trafo["lv_bus"].values, pp_to_ls))
    # pp_net._ppc["internal"]["Bbus"][64, 67]
    # pp_net._ppc["internal"]["branch"][180]
    # trafo_r[7]
    # trafo_x[7]
    # trafo_b[7]
    # import pdb
    # pdb.set_trace()

    for tr_id, is_connected in enumerate(pp_net.trafo["in_service"].values):
        if not is_connected:
            # trafo is deactivated
            model.deactivate_trafo(tr_id)
