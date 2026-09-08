# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Use the pandapower converter to properly initialize a LSGrid c++ object.
"""

from typing import Optional
import numpy as np
from numbers import Number

import pandapower
from ...lightsim2grid_cpp import LSGrid, PandaPowerConverter
from ._aux_add_sgen import _aux_add_sgen
from ._aux_add_load import _aux_add_load
from ._aux_add_trafo import _aux_add_trafo
from ._aux_add_line import _aux_add_line
from ._aux_add_gen import _aux_add_gen
from ._aux_add_shunt import _aux_add_shunt
from ._aux_check_legit import _aux_check_legit
from ._aux_add_slack import _aux_add_slack
from ._aux_add_storage import _aux_add_storage
from ._aux_add_dc_line import _aux_add_dc_line
from ._my_const import ALLOWED_PP_ORIG_FILE
from ._pp_bus_to_ls_bus import pp_bus_to_ls


def init(pp_net: "pandapower.auxiliary.pandapowerNet",
         n_sub: Optional[int]=None,  # number of voltage levels
         n_busbar_per_sub: Optional[int]=None,  # max number of buses allowed per substation / voltage level
         pp_orig_file : ALLOWED_PP_ORIG_FILE = "pandapower_v2",
         init_subid: bool = True,
         ) -> LSGrid:
    """
    Convert a pandapower network as input into a LSGrid.

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
    really well for this grid. Actually, this grid is used for testing the LSGrid class.

    Parameters
    ----------
    pp_net: :class:`pandapower.auxiliary.pandapowerNet`
        The initial pandapower network you want to convert
        
    pp_orig_file: 
        Pandapower change the formula they used internally to compute the "equations" parameters
        of the transformers between pandapower 2.xx and 3.xx.
        
        If you are using a recent (=> 3.xx) version of pandapower, you can pass use the
        ad-hoc trafo converter of lightsim2grid. 
        
        For grid2op environment, we recommed **NOT** to use it if the environment has been released
        before 2026 as the case files came from pandapower 2 (so it's better to use the pandapower 2 
        converter).

    init_subid:
        Whether to tell the resulting `LSGrid` which substation / voltage level each element
        belongs to (`set_gen_to_subid` and friends). Default ``True``, and the same rule as
        every other loader: one substation per pandapower bus, everything on busbar section
        1, and `n_busbar_per_sub` sections allocated per substation.

        :class:`lightsim2grid.lightSimBackend.LightSimBackend` passes ``False`` here, and
        only here. Its pandapower path is built on grid2op's own pandapower backend, which
        has been the one deciding what a substation is for as long as it has existed; the
        `LSGrid` gets its substation ids from there, and this flag is what keeps that true.
        Nothing else should need it.

        .. versionadded:: 1.0.1

    Returns
    -------
    model: :class:`lightsim2grid.network.LSGrid`
        The initialize network

    """
    if pp_orig_file not in ALLOWED_PP_ORIG_FILE.__args__:
        raise RuntimeError(f"pp_orig_file argument should be one of {sorted(ALLOWED_PP_ORIG_FILE.__args__)}")
    
    # check for things not supported and raise if needed
    _aux_check_legit(pp_net)

    # initialize and use converters
    converter = PandaPowerConverter()
    converter.set_sn_mva(pp_net.sn_mva)  # TODO raise an error if not set !
    converter.set_f_hz(pp_net.f_hz)

    # set up the data model accordingly
    model = LSGrid()
    if "_options" in pp_net:
        if "init_vm_pu" in pp_net["_options"]:
            tmp_ = pp_net["_options"]["init_vm_pu"]
            if isinstance(tmp_, Number):
                model.set_init_vm_pu(float(tmp_))
    model.set_sn_mva(pp_net.sn_mva)
    if n_sub is None:
        n_sub = pp_net.bus.shape[0]
        if n_busbar_per_sub is not None and n_busbar_per_sub != 1:
            raise RuntimeError(f"If n_sub is None, n_busbar_per_sub must be None (or 1), found {n_busbar_per_sub}.")
        n_busbar_per_sub = 1
    # input data check
    try:
        tmp = int(n_sub)
    except ValueError as exc_:
        raise RuntimeError("Impossible to convert n_sub to int") from exc_
    if tmp != n_sub:
        raise RuntimeError(f"n_sub should be a int, you provided {tmp} which cannot safely be converted to an int.")
    n_sub = tmp
    if n_sub <= 0:
        raise RuntimeError(f"You need to provide a grid with at least 1 substation / voltage level, provided n_sub={n_sub}")
    
    try:
        tmp = int(n_busbar_per_sub)
    except ValueError as exc_:
        raise RuntimeError("Impossible to convert n_busbar_per_sub to int") from exc_
    if tmp != n_busbar_per_sub:
        raise RuntimeError(f"n_busbar_per_sub should be a int, you provided {tmp} which cannot safely be converted to an int.")
    n_busbar_per_sub = tmp
    if n_busbar_per_sub <= 0:
        raise RuntimeError(f"You need to provide a grid with at least 1 busbar per "
                           f"substation / voltage level, provided n_busbar_per_sub={n_busbar_per_sub}")
    
    tmp_bus_ind = np.argsort(pp_net.bus.index)
    model.init_bus(n_sub,
                   n_busbar_per_sub,
                   pp_net.bus.iloc[tmp_bus_ind]["vn_kv"].to_numpy(),
                   pp_net.line.shape[0],
                   pp_net.trafo.shape[0])
    if np.any(np.sort(pp_net.bus.index) != np.arange(pp_net.bus.shape[0])):
        model._ls_to_orig = 1 * pp_net.bus.index.to_numpy().astype(int)
        pp_to_ls = {pp_bus: ls_bus for pp_bus, ls_bus in zip(pp_net.bus.index, tmp_bus_ind)}
    else:
        pp_to_ls = None
    # deactivate in lightsim the deactivated bus in pandapower
    for bus_id in range(pp_net.bus.shape[0]):
        if not pp_net.bus["in_service"].to_numpy()[bus_id]:
            if pp_to_ls is None:
                pp_bus_id = bus_id
            else:
                pp_bus_id = pp_to_ls[bus_id]
            model.deactivate_bus(pp_bus_id)

    # init the powerlines
    _aux_add_line(converter, model, pp_net, pp_to_ls)

    # init the shunts
    _aux_add_shunt(model, pp_net, pp_to_ls)

    # handle the trafos
    _aux_add_trafo(converter, model, pp_net, pp_to_ls, pp_orig_file)

    # handle loads
    _aux_add_load(model, pp_net, pp_to_ls)

    # handle static generators (PQ generator)
    _aux_add_sgen(model, pp_net, pp_to_ls)

    # handle generators
    _aux_add_gen(model, pp_net, pp_to_ls)

    # handle storage units
    _aux_add_storage(model, pp_net, pp_to_ls)

    # handle dc line
    _aux_add_dc_line(model, pp_net, pp_to_ls)

    # deal with slack bus
    added_gen_bus = _aux_add_slack(model, pp_net, pp_to_ls, pp_orig_file)

    if init_subid:
        # tell the LSGrid which substation / voltage level each element belongs to.
        #
        # One substation per pandapower bus, everything on busbar section 1: lightsim2grid
        # lays its global bus ids out busbar-section-major, so section 1 is [0, n_sub) and
        # section k is [k * n_sub, (k+1) * n_sub) -- an element's substation is therefore
        # its bus modulo n_sub. The modulo is not decoration: `LightSimBackend` hands this
        # loader a pandapower net grid2op has already widened to n_sub * n_busbar_per_sub
        # buses, so a bus id there can genuinely be past the first section.
        #
        # The bus each element was BUILT on, and not `el.bus_id`, because an element
        # pandapower declares out of service reads back as `bus_id == -1` while the
        # substation it belongs to is a property of the grid, not of its status.
        _aux_set_subid(model, pp_net, pp_to_ls, n_sub, added_gen_bus)

    # make sure the grid we just built is internally consistent (bus / substation
    # / topology-vector indices in range, no NaN/Inf in the physical inputs)
    model.check_grid()

    return model


def _aux_set_subid(model, pp_net, pp_to_ls, n_sub, added_gen_bus):
    """Assign every element's substation id, see the call site above."""
    def sub_of(pp_bus_ids):
        """pandapower bus ids -> substation ids"""
        if len(pp_bus_ids) == 0:
            return np.array([], dtype=int)
        return np.asarray(pp_bus_to_ls(np.asarray(pp_bus_ids), pp_to_ls), dtype=int) % n_sub

    # `_aux_add_slack` appends one generator per ext_grid when no generator stands on
    # the slack bus, so the generator container is not always as long as `pp_net.gen`.
    # The buses it returns are already lightsim2grid ones, hence the separate modulo.
    gen_sub = np.concatenate((sub_of(pp_net.gen["bus"].to_numpy()),
                              np.asarray(added_gen_bus, dtype=int) % n_sub))
    model.set_gen_to_subid(gen_sub.astype(int))
    model.set_load_to_subid(sub_of(pp_net.load["bus"].to_numpy()))
    model.set_storage_to_subid(sub_of(pp_net.storage["bus"].to_numpy()))
    model.set_shunt_to_subid(sub_of(pp_net.shunt["bus"].to_numpy()))
    model.set_line_to_sub1_id(sub_of(pp_net.line["from_bus"].to_numpy()))
    model.set_line_to_sub2_id(sub_of(pp_net.line["to_bus"].to_numpy()))
    model.set_trafo_to_sub1_id(sub_of(pp_net.trafo["hv_bus"].to_numpy()))
    model.set_trafo_to_sub2_id(sub_of(pp_net.trafo["lv_bus"].to_numpy()))
