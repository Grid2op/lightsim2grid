# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


def _aux_add_shunt(model, pp_net):
    """
    Add the shunts of the pp_net into the lightsim2grid "model"

    Parameters
    ----------
    model
    pp_net

    Returns
    -------

    """
    model.init_shunt(pp_net.shunt["p_mw"].values,
                     pp_net.shunt["q_mvar"].values,
                     pp_net.shunt["bus"].values
                     )
    for sh_id, is_connected in enumerate(pp_net.shunt["in_service"].values):
        if not is_connected:
            # shunt is deactivated
            model.deactivate_shunt(sh_id)
