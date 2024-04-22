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


def _aux_add_slack(model, pp_net, pp_to_ls):
    """
    add the slack bus(es) to the lightsim2grid "model" based on the information in the pandapower network.

    Notes
    -----
    For now only single slack bus are supported in lightsim2grid. A warning will be issue if the pandapower
    network has a distributed slack bus.

    Parameters
    ----------
    model
    pp_net

    Returns
    -------

    """
    # TODO handle that better maybe, and warn only one slack bus is implemented
    slack_coeff = None
    if np.any(pp_net.gen["slack"].values):
        # most favorable cases
        # if np.sum(pp_net.gen["slack"].values) >= 2:
        #     # TODO SLACK remove this warning ! (will be done at the end when more tests will be done)
        #     warnings.warn("LightSim cannot handle multiple slack bus at the moment. Only the first "
        #                   "slack bus of pandapower will be used.")
        slack_gen_ids = np.where(pp_net.gen["slack"].values)[0]
        if "slack_weight" in pp_net.gen:
            slack_coeff = pp_net.gen["slack_weight"].values[slack_gen_ids]
        for slack_id in slack_gen_ids:
            model.change_v_gen(slack_id, pp_net.gen["vm_pu"].iloc[slack_id])

        if slack_coeff is not None:
            has_nan = np.any(~np.isfinite(slack_coeff))
            if has_nan:
                warnings.warn("Some slack coefficients were Nans. We set them all to 1.0")
                slack_coeff[:] = 1.0
            
        # in this case the ext grid is not taken into account, i raise a warning if
        # there is one
        slack_bus_ids = pp_bus_to_ls(pp_net.ext_grid["bus"].values, pp_to_ls)
        if pp_net.ext_grid.shape[0] >= 1:
            warnings.warn("LightSim will not consider the pandapower \"ext_grid\" as there "
                          "are already generators tagged as slack bus")
    else:
        # there is no slack bus in the generator of the pp grid
        warnings.warn("LightSim has not found any generators tagged as \"slack bus\" in the pandapower network."
                      "I will attempt to add some from the ext_grid.")
        # first i try to see if a generator is connected to a slack bus
        # TODO SLACK: deactivate warnings (will be done at the end when more tests will be done)
        slack_bus_ids = pp_bus_to_ls(pp_net.ext_grid["bus"].values, pp_to_ls)
        # if pp_net.ext_grid.shape[0] >= 2:
        #     warnings.warn("LightSim cannot handle multiple slack bus at the moment. Only the first "
        #                   "slack bus of pandapower will be used.")

        if np.all(np.isin(slack_bus_ids, pp_bus_to_ls(pp_net.gen["bus"].values, pp_to_ls))):
            # all slack buses have a generator connected to them
            # so i assume it was just a computation artifact, and assign these generators as slack buses
            slack_gen_ids = np.isin(pp_bus_to_ls(pp_net.gen["bus"].values, pp_to_ls), slack_bus_ids)  # id of generators connected to slack bus
            slack_gen_ids = np.where(slack_gen_ids)[0]  # keep only the id of the generators
            if "slack_weight" in pp_net.gen:
                slack_coeff = pp_net.gen["slack_weight"].values[slack_gen_ids]
                
            has_nan = np.any(~np.isfinite(slack_coeff))
            if has_nan:
                warnings.warn("We found some Nans in the slack coefficients. "
                              "We set them all to the same value !")
                slack_coeff[:] = 1.0
        else:
            # at least one slack bus has no generator connected to it
            # so I assume i need to add as many generators as number of slack bus
            nb_slack = len(slack_bus_ids)
            if "slack_weight" in pp_net.ext_grid:
                slack_coeff = 1.0 * pp_net.ext_grid["slack_weight"].values
            else:
                slack_coeff = np.ones(nb_slack)
            slack_coeff_norm = slack_coeff / slack_coeff.sum()
            
            has_nan = np.any(~np.isfinite(slack_coeff_norm))
            if has_nan:
                warnings.warn("We found some Nans in the slack coefficients \"slack_coeff_norm\", "
                              " (probably because the slack weights sum to 0.0 initially)"
                              "We set them all to the same value !")
                slack_coeff_norm[:] = 1.0 / slack_coeff_norm.shape[0]
                slack_coeff[:] = 1.0
                
            slack_gen_ids = np.arange(nb_slack) + pp_net.gen.shape[0]
            slack_contrib = -1.0 * (np.sum(pp_net.gen["p_mw"]) - np.sum(pp_net.load["p_mw"]) ) * slack_coeff_norm
            vm_pu = 1.0 * pp_net.ext_grid["vm_pu"].values
            gen_p = np.concatenate((pp_net.gen["p_mw"].values, slack_contrib))
            gen_v = np.concatenate((pp_net.gen["vm_pu"].values, vm_pu))
            gen_bus = np.concatenate((pp_bus_to_ls(pp_net.gen["bus"].values, pp_to_ls), slack_bus_ids))
            gen_min_q = np.concatenate((pp_net.gen["min_q_mvar"].values, [-999999. for _ in range(nb_slack)]))
            gen_max_q = np.concatenate((pp_net.gen["max_q_mvar"].values, [+99999. for _ in range(nb_slack)]))
            model.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus)

    # handle the possible distributed slack bus
    if slack_coeff is None:
        slack_coeff = np.ones(len(slack_gen_ids))
    if np.sum(slack_coeff) == 0. or np.any(slack_coeff < 0.):
        warnings.warn("We found either some slack coefficient to be < 0. or they were all 0."
                      "We set them all to 1.0 to avoid such issues")
        slack_coeff[:] = 1.

    for gid, slack_gen_id in enumerate(slack_gen_ids):
        model.add_gen_slackbus(slack_gen_id, slack_coeff[gid])
