# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Use the pandapower converter to properly initialize a GridModel c++ object.
"""

__all__ = ["init", "GridModel"]

import numpy as np
from numbers import Number
import warnings
from lightsim2grid_cpp import GridModel, PandaPowerConverter
from lightsim2grid.gridmodel._aux_add_sgen import _aux_add_sgen
from lightsim2grid.gridmodel._aux_add_load import _aux_add_load
from lightsim2grid.gridmodel._aux_add_trafo import _aux_add_trafo
from lightsim2grid.gridmodel._aux_add_line import _aux_add_line
from lightsim2grid.gridmodel._aux_add_gen import _aux_add_gen
from lightsim2grid.gridmodel._aux_add_shunt import _aux_add_shunt
from lightsim2grid.gridmodel._aux_check_legit import _aux_check_legit
from lightsim2grid.gridmodel._aux_add_slack import _aux_add_slack
from lightsim2grid.gridmodel._aux_add_storage import _aux_add_storage


def init(pp_net):
    """
    Convert a pandapower network as input into a GridModel.

    This can fail to convert the grid and still not throw any error, use with care (for example, you can run a powerflow
    after this conversion, run a powerflow with pandapower, and compare the results to make sure they match !)

    Cases for which conversion is not possible include, but are not limited to:

    - the pandapower grid has 3 winding transformers
    - the pandapower grid has xwards
    - the pandapower grid has dcline
    - the pandapower grid has switch, motor, assymetric loads, etc.
    - the pandapower grid any parrallel "elements" (at least one of the column "parrallel" is not 1)
    - the bus indexes in pandapower do not start at 0 or are not contiguous (you can check `pp_net.bus.index`)
    - some `g_us_per_km` for some lines are not zero ? TODO not sure if that is still the case !
    - some `p_mw` for some shunts are not zero ? TODO not sure if that is still the case !

    if you really need any of the above, please submit a github issue and we will work on their support.

    This conversion has been extensively studied for the case118() of pandapower.networks and should work
    really well for this grid. Actually, this grid is used for testing the GridModel class.

    Parameters
    ----------
    pp_net: :class:`pandapower.grid`
        The initial pandapower network you want to convert

    Returns
    -------
    model: :class:`lightsim2grid.gridmodel.GridModel`
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
            tmp_ = pp_net["_options"]["init_vm_pu"]
            if isinstance(tmp_, Number):
                model.set_init_vm_pu(float(tmp_))
    model.set_sn_mva(pp_net.sn_mva)

    tmp_bus_ind = np.argsort(pp_net.bus.index)
    model.init_bus(pp_net.bus.iloc[tmp_bus_ind]["vn_kv"].values,
                   pp_net.line.shape[0],
                   pp_net.trafo.shape[0])

    # deactivate in lightsim the deactivated bus in pandapower
    for bus_id in range(pp_net.bus.shape[0]):
        if not pp_net.bus["in_service"].values[bus_id]:
            model.deactivate_bus(bus_id)

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
