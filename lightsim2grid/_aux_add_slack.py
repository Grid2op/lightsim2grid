# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
import warnings
import numpy as np


def _aux_add_slack(model, pp_net):
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
    if np.any(pp_net.gen["slack"].values):
        if np.sum(pp_net.gen["slack"].values) >= 2:
            warnings.warn("LightSim cannot handle multiple slack bus at the moment. Only the first "
                          "slack bus of pandapower will be used.")
        slack_gen_id = np.where(pp_net.gen["slack"].values)[0]
        model.change_v_gen(slack_gen_id, pp_net.gen["vm_pu"][slack_gen_id])
    else:
        # there is no slack bus in the generator of the pp grid

        # first i try to see if a generator is connected to a slack bus
        slack_bus_id = pp_net.ext_grid["bus"].values[0]
        if pp_net.ext_grid.shape[0] >= 2:
            warnings.warn("LightSim cannot handle multiple slack bus at the moment. Only the first "
                          "slack bus of pandapower will be used.")

        if np.any(pp_net.gen["bus"].values == slack_bus_id):
            slack_gen_id = np.where(pp_net.gen["bus"].values == slack_bus_id)[0]
        else:
            # no gen is connected to a slack bus, so i create one.
            gen_p = np.concatenate((pp_net.gen["p_mw"].values, [np.sum(pp_net.load["p_mw"]) - np.sum(pp_net.gen["p_mw"])]))
            gen_v = np.concatenate((pp_net.gen["vm_pu"].values, [pp_net.ext_grid["vm_pu"].values[0]]))
            gen_bus = np.concatenate((pp_net.gen["bus"].values, [slack_bus_id]))
            gen_min_q = np.concatenate((pp_net.gen["min_q_mvar"].values, [-999999.]))
            gen_max_q = np.concatenate((pp_net.gen["max_q_mvar"].values, [+99999.]))
            model.init_generators(gen_p, gen_v, gen_min_q, gen_max_q, gen_bus)
            slack_gen_id = pp_net.gen["bus"].shape[0]

    model.add_gen_slackbus(slack_gen_id)

