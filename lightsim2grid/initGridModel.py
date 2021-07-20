# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Use the pandapower converter to properly initialized a GridModel c++ object.
"""

import numpy as np
import warnings
from lightsim2grid_cpp import GridModel, PandaPowerConverter, SolverType
from lightsim2grid._aux_add_sgen import _aux_add_sgen
from lightsim2grid._aux_add_load import _aux_add_load
from lightsim2grid._aux_add_trafo import _aux_add_trafo
from lightsim2grid._aux_add_line import _aux_add_line
from lightsim2grid._aux_add_gen import _aux_add_gen
from lightsim2grid._aux_add_shunt import _aux_add_shunt
from lightsim2grid._aux_check_legit import _aux_check_legit
from lightsim2grid._aux_add_slack import _aux_add_slack
from lightsim2grid._aux_add_storage import _aux_add_storage


def init(pp_net):
    """
    Convert a pandapower network as input into a GridModel.

    This does not throw any error at the moment when the conversion is not possible.

    Cases for which conversion is not possible include, but are not limited to:

    - the pandapower grid has 3 winding transformers
    - the pandapower grid has xwards
    - the pandapower grid any parrallel "elements" (at least one of the column "parrallel" is not 1)
    - some `g_us_per_km` for some lines are not zero
    - some `p_mw` for some shunts are not zero
    - some `tap_step_degre` are non zero for some trafo
    - no "ext_grid" is reported on the initial grid

    if you really need any of the above, please submit a github issue and we will work on their support.

    This conversion has been extensively studied for the case118() of pandapower.networks and should work
    really well for this grid. Actually, this grid is used for testing the GridModel class.

    Parameters
    ----------
    pp_net: :class:`pandapower.grid`
        The initial pandapower network you want to convert

    Returns
    -------
    model: :class:`GridModel`
        The initialize gridmodel

    """
    # check for things not supported and raise if needed
    _aux_check_legit(pp_net)

    # initialize and use converters
    converter = PandaPowerConverter()
    converter.set_sn_mva(pp_net.sn_mva)  # TODO raise an error if not set !
    converter.set_f_hz(pp_net.f_hz)

    # set up the data model accordingly
    model = GridModel()
    if "_options" in pp_net:
        if "init_vm_pu" in pp_net["_options"]:
            model.set_init_vm_pu(pp_net["_options"]["init_vm_pu"])

    model.set_sn_mva(pp_net.sn_mva)

    tmp_bus_ind = np.argsort(pp_net.bus.index)
    model.init_bus(pp_net.bus.iloc[tmp_bus_ind]["vn_kv"].values,
                   pp_net.line.shape[0],
                   pp_net.trafo.shape[0])

    # init the powerlines
    _aux_add_line(converter, model, pp_net)

    # init the shunts
    _aux_add_shunt(model, pp_net)

    # handle the trafos
    _aux_add_trafo(converter, model, pp_net)

    # handle loads
    _aux_add_load(model, pp_net)

    # handle static generators (PQ generator)
    _aux_add_sgen(model, pp_net)

    # handle generators
    _aux_add_gen(model, pp_net)

    # handle storage units
    _aux_add_storage(model, pp_net)

    # deal with slack bus
    _aux_add_slack(model, pp_net)

    return model
