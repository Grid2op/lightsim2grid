# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import copy
import numpy as np
import pandas as pd
import pypowsybl as pypo
from typing import Dict, Iterable, Optional, Union
from packaging import version
from lightsim2grid_cpp import GridModel


from ._aux_handle_slack import handle_slack_iterable, handle_slack_one_el

PP_BUG_RATIO_TAP_CHANGER = version.parse("1.9")
PYPOWSYBL_VER = version.parse(pypo.__version__)


def _aux_get_bus(vl_df, bus_df, df, conn_key="connected", bus_key="bus_id", vl_key="voltage_level_id"):
    if df.shape[0] == 0:
        # no element of this type so no problem
        return np.zeros(0, dtype=int), np.ones(0, dtype=bool), np.zeros(0, dtype=int)
    # retrieve which elements are disconnected 
    mask_disco = ~df[conn_key]
    
    # retrieve the bus where the element are
    tmp_bus_id = df[bus_key].copy()
    # element disconnected are, by default assigned to first bus of their substation
    tmp_disco = vl_df.loc[df.loc[mask_disco, vl_key], "vl_id"].values
    bus_el_disco = bus_df.iloc[tmp_disco]
    tmp_bus_id[mask_disco] = bus_el_disco.index
    bus_id = bus_df.loc[tmp_bus_id.values]["bus_global_id"].values
    # deactivate the element not on the main component
    # wrong_component = bus_df.loc[tmp_bus_id.values]["connected_component"].values != 0
    # mask_disco[wrong_component] = True
    # assign bus -1 to disconnected elements
    bus_id[mask_disco] = -1
    
    sub_id = bus_df.loc[tmp_bus_id.values]["glop_sub_id"].values
    return bus_id, mask_disco.values, sub_id


def init(net : pypo.network.Network,
         gen_slack_id: Union[int, str, Iterable[str], Dict[str, float]] = None,
         slack_bus_id: int = None,
         sn_mva : float = 100.,  # only used if not present in the grid
         sort_index : bool =True, 
         f_hz : float = 50.,  # unused
         net_pu : Optional[pypo.network.Network] = None,
         only_main_component : bool =True,
         return_sub_id: bool=False,
         n_busbar_per_sub: Optional[int]=None,  # new in 0.9.1
         buses_for_sub:Optional[bool]=None,  # new in 0.9.1
         init_vm_pu:float=1.06,
         ) -> GridModel:
    """
    This function is available under the `init_from_pypowsybl` in lightsim2grid
    

    .. code-block:: python
    
        from lightsim2grid.gridmodel import init_from_pypowsybl
        
    .. warning::
        It is not available if the `pypowsybl` python package is not installed.
    
    :param net: The pypowsybl network
    :type net: pypo.network.Network
    
    :param gen_slack_id: The id of the generator that should be used as the slack
                         (either it's given by id (int) or by name (str))
    :type gen_slack_id: Union[int, str]
    
    :param slack_bus_id: If you don't provide a generator ID as a slack bus, you can
            provide a bus id (int). We do not recommend setting the slack this way.            
    :type slack_bus_id: int
    
    :param sn_mva: The nominal apparent power used when converting the grid to 
                   per unit. It is only used if the pypowsybl grid 
                   has no `_nominal_apparent_power` attribute. 
                   **Advanced usage**.            
    :type sn_mva: float
    
    :param sort_index: Whether you want to sort the indexes of all the 
                       pypowsybl tables (*eg* get_loads() or *get_buses()*) or not.
                       Sorting the grid tables is preferable if you want to be 
                       "future proof" and don't want to depend on pandas version
                       (same order is guaranteed). Not sorting the grid will give
                       easier comparison of results with pypowsybl.          
    :type sn_mva: bool
    
    :param f_hz: Not used currently (frequency of the grid)
    :type net_pu: float
    
    :param net_pu: If you have already converted the grid in "per unit" then
                   you can pass it as the `net_pu` argument. Otherwise this 
                   function will do it.
                   **Advanced usage**.
    :type net_pu: Optional[pypo.network.Network]
    
    :param only_main_component: If this is True, then only the main component (*ie*
                                the one containing the slack bus) will be used. All
                                equipments not part of this component will be
                                deactivated (switched-off). **NB** currently
                                lightsim2grid will diverge if the grid is not connected,
                                this option might then "hide" some equipements from
                                the grid (silently) but you have higher chances of
                                convergence.
    :type only_main_component: bool
    
    :param return_sub_id: **Advanced usage**. If you want to retrieve the id of the
                          equipments as "tables". Used only for `LightSimBackend`
    :type return_sub_id: bool
                          
    :param n_busbar_per_sub: Currently, lightsim2grid works well with a constant 
                             number of independant buses that can be made at each 
                             substations. It can be infered from the grid or 
                             set with this attribute. We recommend to leave it 
                             to `None` (which corresponds to the "infer it from 
                             the grid" behaviour) in most cases.
    :type n_busbar_per_sub: Opional[int]
    
    :param buses_for_sub: Whether the lightsim2grid substation will correspond to buses
                          of the pypowsybl grid (if buses_for_sub is `True`).
                          Alternatively, if buses_for_sub is `False`, the
                          lightsim2grid susbtation will correspond to
                          pypowsybl voltage level (read from net.get_voltage_levels()).
                          buses_for_sub==`True` is a "legacy" behaviour.
    :type buses_for_sub: bool
    
    :param init_vm_pu: The voltage magnitude with which the init vector of AC powerflow
                       will be set.
    :type init_vm_pu: float
    
    :return: The properly initialized gridmodel.
    :rtype: :class:`GridModel`
    """
    model = GridModel()
    if hasattr(net, "_nominal_apparent_power"):
        sn_mva_used = getattr(net, "_nominal_apparent_power")
    else:
        sn_mva_used = float(sn_mva)
    model.set_sn_mva(sn_mva_used)
    model.set_init_vm_pu(float(init_vm_pu))
    
    if gen_slack_id is not None and slack_bus_id is not None:
        raise RuntimeError("Impossible to intialize a grid with both gen_slack_id and slack_bus_id")
    
    # assign unique id to the buses
    bus_df_orig = net.get_buses()
    if sort_index:
        bus_df = bus_df_orig.sort_index().copy()
    else:
        bus_df = bus_df_orig.copy()
    bus_df["orig_id"] = np.arange(bus_df.shape[0])
    
    if sort_index:
        voltage_levels = net.get_voltage_levels().sort_index()
    else:
        voltage_levels = net.get_voltage_levels()
    
    all_buses_vn_kv = voltage_levels.loc[bus_df["voltage_level_id"].values]["nominal_v"].values
    sub_unique = None
    sub_unique_id = None
    nb_bus_per_vl = bus_df[["voltage_level_id", "name"]].groupby("voltage_level_id").count()
    
    if buses_for_sub is not None and buses_for_sub:
        # I am in a compatibility mode,
        # the "substation" in lightsim2grid will be read
        # from the buses in the original grid (and not from the
        # voltage levels)
        if n_busbar_per_sub is None:
            # setting automatically n_busbar_per_sub
            # to 1
            # TODO logger here
            n_busbar_per_sub = 1
            
        all_buses_vn_kv = voltage_levels.loc[bus_df["voltage_level_id"], "nominal_v"].values
        if n_busbar_per_sub > 1:
            all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
        n_sub_ls = bus_df.shape[0]
        ls_to_orig = np.zeros(all_buses_vn_kv.shape[0], dtype=int) - 1
        ls_to_orig[:n_sub_ls] = np.arange(n_sub_ls)
        n_busbar_per_sub_ls = n_busbar_per_sub
        bus_df["bus_global_id"] = np.arange(n_sub_ls)
        bus_df["glop_sub_id"] = np.arange(n_sub_ls)  # np.concatenate([np.arange(n_sub_ls) for _ in range(n_busbar_per_sub)])
        sub_names = bus_df.index.values.astype(str)
        voltage_levels["vl_id"] = bus_df[["voltage_level_id", "bus_global_id"]].groupby("voltage_level_id").min()
    else:        
        # the "substation" in lightsim2grid
        voltage_levels["nb_bus_per_vl"] = nb_bus_per_vl["name"]        
        bus_df["name"] = [[el] for el in bus_df.index]
        voltage_levels["bus_names"] = bus_df[["name", "voltage_level_id"]].groupby("voltage_level_id").sum()   

        bus_df["local_id"] = [voltage_levels.loc[el, "bus_names"].index(id_) + 1 
                            for id_, el in zip(bus_df.index,
                                                bus_df["voltage_level_id"].values)]
        n_vl = voltage_levels.shape[0]
        voltage_levels["vl_id"] = np.arange(n_vl)
        bus_df["bus_global_id"] = [(loc_id - 1) * n_vl + voltage_levels.loc[vl, "vl_id"]
                                   for loc_id, vl in zip(
                                       bus_df["local_id"],
                                       bus_df["voltage_level_id"]
                                   )]
        nb_bus_per_vl_in_grid = nb_bus_per_vl.values.max()
        if n_busbar_per_sub is None:
            # setting automatically n_busbar_per_sub
            # to the value read from the grid
            # TODO logger here
            n_busbar_per_sub = int(nb_bus_per_vl_in_grid)
        elif n_busbar_per_sub < nb_bus_per_vl_in_grid:
            raise RuntimeError(f"The input pypowsybl grid counts some voltage levels "
                               f"with {nb_bus_per_vl_in_grid} independant buses, "
                               f"which is not compatible with the n_busbar_per_sub={n_busbar_per_sub} "
                               "given as input.")
        all_buses_vn_kv = voltage_levels["nominal_v"].values
        if n_busbar_per_sub > 1:
            all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
        n_sub_ls = voltage_levels.shape[0]
        voltage_levels["glop_sub_id"] = np.arange(voltage_levels.shape[0])
        n_busbar_per_sub_ls = n_busbar_per_sub
        ls_to_orig = np.zeros(all_buses_vn_kv.shape[0], dtype=int) - 1
        ls_to_orig[bus_df["bus_global_id"].values] = np.arange(bus_df.shape[0])
        sub_names = voltage_levels.index.values.astype(str)
        bus_df["glop_sub_id"] = voltage_levels.loc[bus_df["voltage_level_id"].values, "glop_sub_id"].values
        
    # all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
    model.init_bus(n_sub_ls,
                   n_busbar_per_sub_ls,
                   all_buses_vn_kv,
                   0, 0  # unused
                   )
    model._ls_to_orig = ls_to_orig
    model._max_nb_bus_per_sub = n_busbar_per_sub_ls
    model.init_substation_names(sub_names)
        
    # do the generators
    if sort_index:
        df_gen = net.get_generators().sort_index()
    else:
        df_gen = net.get_generators()
        
    # to handle encoding in 32 bits and overflow when "splitting" the Q values among 
    min_float_value = np.finfo(np.float32).min * 1e-4 + 1.
    max_float_value = np.finfo(np.float32).max * 1e-4 + 1.
    min_q_aux = 1. * df_gen["min_q"].values
    too_small = min_q_aux < min_float_value
    min_q_aux[too_small] = min_float_value
    min_q = min_q_aux.astype(np.float32)
    
    max_q_aux = 1. * df_gen["max_q"].values
    too_big = np.abs(max_q_aux) > max_float_value
    max_q_aux[too_big] = np.sign(max_q_aux[too_big]) * max_float_value
    max_q = max_q_aux.astype(np.float32)
    min_q[~np.isfinite(min_q)] = min_float_value
    max_q[~np.isfinite(max_q)] = max_float_value
    gen_bus, gen_disco, gen_sub = _aux_get_bus(voltage_levels, bus_df, df_gen)

    # dirty fix for when regulating elements are not the same
    bus_reg = copy.deepcopy(df_gen["regulated_element_id"].values)
    # for oldest pypowsybl version, we could have "" there
    bus_reg = np.where(bus_reg == "", df_gen.index, bus_reg)
    vl_reg = copy.deepcopy(df_gen["voltage_level_id"].values)
    mask_ref_bbs = bus_reg != df_gen.index
    
    bbs_df = net.get_busbar_sections().copy()
    if not (np.isin(bus_reg[mask_ref_bbs], bbs_df.index)).all():
        raise RuntimeError("At least some generator are in 'remote control' mode "
                           "and does not control a busbar section, this is not supported "
                           "at the moment.")
    vl_reg[mask_ref_bbs] = bbs_df.loc[bus_reg[mask_ref_bbs], "voltage_level_id"].values
    model.init_generators_full(df_gen["target_p"].values,
                            #    df_gen["target_v"].values / voltage_levels.loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
                               df_gen["target_v"].values / voltage_levels.loc[vl_reg]["nominal_v"].values,
                               df_gen["target_q"].values,
                               df_gen["voltage_regulator_on"].values,
                               min_q,
                               max_q,
                               gen_bus
                               )
    for gen_id, is_disco in enumerate(gen_disco):
        if is_disco:
            model.deactivate_gen(gen_id)
    model.set_gen_names(df_gen.index)   
    
    # for loads
    if sort_index:
        df_load = net.get_loads().sort_index()
    else:
        df_load = net.get_loads()
    load_bus, load_disco, load_sub = _aux_get_bus(voltage_levels, bus_df, df_load)
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     load_bus
                     )
    for load_id, is_disco in enumerate(load_disco):
        if is_disco:
            model.deactivate_load(load_id)
    model.set_load_names(df_load.index)
    
    # for lines
    if sort_index:
        df_line = net.get_lines().sort_index()
    else:
        df_line = net.get_lines()
        
    # per unit
    if net_pu is None:
        if hasattr(net, "per_unit"):
            net_pu = copy.deepcopy(net)
            net_pu.per_unit = True
        else:
            # legacy pypowsybl mode: this did not exist
            from pypowsybl.network import PerUnitView
            net_pu = PerUnitView(net)
            warnings.warn("The `PerUnitView` (python side) is less efficient and less "
                          "tested that the equivalent java class. Please upgrade pypowsybl version")
    df_line_pu = net_pu.get_lines().loc[df_line.index]
    line_r = df_line_pu["r"].values
    line_x = df_line_pu["x"].values
    line_h_or = (df_line_pu["g1"].values + 1j * df_line_pu["b1"].values)
    line_h_ex = (df_line_pu["g2"].values + 1j * df_line_pu["b2"].values)
    lor_bus, lor_disco, lor_sub = _aux_get_bus(voltage_levels, bus_df, df_line, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    lex_bus, lex_disco, lex_sub = _aux_get_bus(voltage_levels, bus_df, df_line, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_powerlines_full(line_r,
                               line_x,
                               line_h_or,
                               line_h_ex,
                               lor_bus,
                               lex_bus
                              )
    for line_id, (is_or_disc, is_ex_disc) in enumerate(zip(lor_disco, lex_disco)):
        if is_or_disc or is_ex_disc:
            model.deactivate_powerline(line_id)
    model.set_line_names(df_line.index) 
    
    # for trafo
    # I extract trafo with `all_attributes=True` so that I have access to the `rho`
    try:
        df_trafo_not_sorted = net.get_2_windings_transformers(all_attributes=True)
    except TypeError:
        # not available in legacy pypowsybl version
        df_trafo_not_sorted = net.get_2_windings_transformers()
        
    if sort_index:
        df_trafo = df_trafo_not_sorted.sort_index()
    else:
        df_trafo = df_trafo_not_sorted
    
    try :
        df_trafo_pu = net_pu.get_2_windings_transformers(all_attributes=True)
    except TypeError:
        df_trafo_pu = net_pu.get_2_windings_transformers()
    df_trafo_pu = df_trafo_pu.loc[df_trafo.index]
    ratio_tap_changer = net_pu.get_ratio_tap_changers()
    
    if 'alpha' in df_trafo_pu:
        shift_ = np.rad2deg(df_trafo_pu['alpha'].values)  # given in radian by pypowsybl
    else:
        if net.get_phase_tap_changers().shape[0] > 0:
            raise RuntimeError("Phase tap changer are not handled by the pypowsybl converter "
                               "when not accessible using the 'alpha' columns "
                               "of the net (once per unit). Please upgrade pypowsybl."
                               "NB: phase tap change are handled by lightsim2grid)")
        shift_ = np.zeros(df_trafo.shape[0])
    is_tap_hv_side = np.zeros(df_trafo.shape[0], dtype=bool)  # TODO    
    trafo_r = df_trafo_pu["r"].values
    trafo_x = df_trafo_pu["x"].values
    trafo_h = (df_trafo_pu["g"].values + 1j * df_trafo_pu["b"].values)
    
    # now get the ratio    
    # in lightsim2grid (cpp)
    if "rho" in df_trafo_pu:
        ratio = 1. * df_trafo_pu["rho"].values
    else:
        # in powsybl (https://javadoc.io/doc/com.powsybl/powsybl-core/latest/com/powsybl/iidm/network/TwoWindingsTransformer.html)
        #  rho = transfo.getRatedU2() / transfo.getRatedU1()
        # * (transfo.getRatioTapChanger() != null ? transfo.getRatioTapChanger().getCurrentStep().getRho() : 1);
        # * (transfo.getPhaseTapChanger() != null ? transfo.getPhaseTapChanger().getCurrentStep().getRho() : 1);

        ratio = 1. * (df_trafo_pu["rated_u2"].values / df_trafo_pu["rated_u1"].values)
        has_r_tap_changer = np.isin(df_trafo_pu.index, ratio_tap_changer.index)
    
        if PYPOWSYBL_VER <= PP_BUG_RATIO_TAP_CHANGER:
            # bug in per unit view in both python and java
            ratio[has_r_tap_changer] = 1. * ratio_tap_changer.loc[df_trafo_pu.loc[has_r_tap_changer].index, "rho"].values

    tor_bus, tor_disco, tor_sub = _aux_get_bus(voltage_levels, bus_df, df_trafo, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    tex_bus, tex_disco, tex_sub = _aux_get_bus(voltage_levels, bus_df, df_trafo, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_h,
                     ratio,
                     shift_,  # in degree !
                     is_tap_hv_side,
                     tor_bus, # TODO do I need to change hv / lv
                     tex_bus)    
    for t_id, (is_or_disc, is_ex_disc) in enumerate(zip(tor_disco, tex_disco)):
        if is_or_disc or is_ex_disc:
            model.deactivate_trafo(t_id)
    model.set_trafo_names(df_trafo.index)
    
    # for shunt
    if sort_index:
        df_shunt = net.get_shunt_compensators().sort_index()
    else:
        df_shunt = net.get_shunt_compensators()
        
    sh_bus, sh_disco, sh_sub = _aux_get_bus(voltage_levels, bus_df, df_shunt)    
    shunt_kv = voltage_levels.loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(df_shunt["g"].values * shunt_kv**2,
                     -df_shunt["b"].values * shunt_kv**2,
                     sh_bus
                    )
    for shunt_id, disco in enumerate(sh_disco):
        if disco:
           model.deactivate_shunt(shunt_id) 
    model.set_shunt_names(df_shunt.index)
           
    # for hvdc (TODO not tested yet)
    if sort_index:
        df_dc = net.get_hvdc_lines().sort_index()
        df_sations = net.get_vsc_converter_stations().sort_index()
    else:
        df_dc = net.get_hvdc_lines()
        df_sations = net.get_vsc_converter_stations()
    hvdc_bus_from_id, hvdc_from_disco, hvdc_sub_from_id  = _aux_get_bus(voltage_levels, bus_df, df_sations.loc[df_dc["converter_station1_id"].values]) 
    hvdc_bus_to_id, hvdc_to_disco, hvdc_sub_to_id = _aux_get_bus(voltage_levels, bus_df, df_sations.loc[df_dc["converter_station2_id"].values]) 
    loss_percent = np.zeros(df_dc.shape[0])  # TODO 
    loss_mw = np.zeros(df_dc.shape[0])  # TODO
    model.init_dclines(hvdc_bus_from_id,
                       hvdc_bus_to_id,
                       df_dc["target_p"].values,
                       loss_percent,
                       loss_mw,
                       voltage_levels.loc[df_sations.loc[df_dc["converter_station1_id"].values]["voltage_level_id"].values]["nominal_v"].values,
                       voltage_levels.loc[df_sations.loc[df_dc["converter_station2_id"].values]["voltage_level_id"].values]["nominal_v"].values,
                       df_sations.loc[df_dc["converter_station1_id"].values]["min_q"].values,
                       df_sations.loc[df_dc["converter_station1_id"].values]["max_q"].values,
                       df_sations.loc[df_dc["converter_station2_id"].values]["min_q"].values,
                       df_sations.loc[df_dc["converter_station2_id"].values]["max_q"].values
                       )
    # TODO will probably not work !
    for hvdc_id, (is_or_disc, is_ex_disc) in enumerate(zip(hvdc_from_disco, hvdc_to_disco)):
        if is_or_disc or is_ex_disc:
            model.deactivate_hvdc(hvdc_id)
    model.set_dcline_names(df_sations.index)
                
    # storage units  (TODO not tested yet)
    if sort_index:
        df_batt = net.get_batteries().sort_index()
    else:
        df_batt = net.get_batteries()
    batt_bus, batt_disco, batt_sub = _aux_get_bus(voltage_levels, bus_df, df_batt)
    model.init_storages(df_batt["target_p"].values,
                        df_batt["target_q"].values,
                        batt_bus
                        )
    for batt_id, disco in enumerate(batt_disco):
        if disco:
           model.deactivate_storage(batt_id) 
    model.set_storage_names(df_batt.index)

    # TODO dist slack
    if gen_slack_id is None and slack_bus_id is None:
        # if nothing is given, by default I assign a slack bus to a bus where a lot of lines are connected
        # quite central in the grid
        bus_id, gen_id = model.assign_slack_to_most_connected()
        gen_slack_ids_int = [gen_id]
    elif gen_slack_id is not None:
        if slack_bus_id is not None:
            raise RuntimeError("You provided both gen_slack_id and slack_bus_id "
                               "which is not possible.")
        if isinstance(gen_slack_id, (str, int, np.int32, np.int64, np.str_, tuple)):
            single_slack = True
            fun_slack = handle_slack_one_el
        else:
            single_slack = False
            fun_slack = handle_slack_iterable
        gen_slack_ids_int, gen_slack_weights = fun_slack(df_gen, gen_slack_id)
        if single_slack:
            if gen_slack_weights is None:
                raise RuntimeError(f"The slack {gen_slack_id} is disconnected.")
            gen_slack_ids_int = [gen_slack_ids_int]
            gen_slack_weights_fixed = [1.]
        else:
            gen_slack_weights_fixed = np.asarray([el if el is not None else np.nan for el in gen_slack_weights])
            mask_finite = np.isfinite(gen_slack_weights_fixed)
            if not mask_finite.any():
                raise RuntimeError(f"No connected generators match the slack {gen_slack_id}")
            gen_slack_weights_fixed[mask_finite] /= gen_slack_weights_fixed[mask_finite].sum()

        for gen_slack_id_int, gen_slack_weight in zip(gen_slack_ids_int, gen_slack_weights_fixed):
            if np.isfinite(gen_slack_weight):
                model.add_gen_slackbus(gen_slack_id_int, gen_slack_weight)
    elif slack_bus_id is not None:
        gen_bus = np.array([el.bus_id for el in model.get_generators()])
        gen_is_conn_slack = gen_bus == model._orig_to_ls[slack_bus_id]
        nb_conn = gen_is_conn_slack.sum()
        if nb_conn == 0:
            raise RuntimeError(f"There is no generator connected to bus {slack_bus_id}. It cannot be the slack")
        gen_slack_ids_int = []
        for gen_id, is_slack in enumerate(gen_is_conn_slack):
            if is_slack:
                gen_slack_ids_int.append(gen_id)
                model.add_gen_slackbus(gen_id, 1. / nb_conn)    
    else:
        raise RuntimeError("You need to provide at least one slack with `gen_slack_id` or `slack_bus_id`") 
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
        
    # and now deactivate all elements and nodes not in the main component
    if only_main_component:
        model.consider_only_main_component()
    else:
        # automatically disconnect non connected buses
        # (this is automatically done by consider_only_main_component)
        model.init_bus_status()  
        
    gen_sub = pd.DataFrame(index=df_gen.index, data={"sub_id": gen_sub})
    gen_sub["desired_slack"] = False
    gen_sub.loc[gen_sub.index[gen_slack_ids_int], "desired_slack"] = True
    load_sub = pd.DataFrame(index=df_load.index, data={"sub_id": load_sub})
    lor_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lor_sub})
    lex_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lex_sub})
    tor_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tor_sub})
    tex_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tex_sub})
    batt_sub = pd.DataFrame(index=df_batt.index, data={"sub_id": batt_sub})
    sh_sub = pd.DataFrame(index=df_shunt.index, data={"sub_id": sh_sub})
    hvdc_sub_from_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_from_id})
    hvdc_sub_to_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_to_id})
    
    # set the substation ID to which each object belong
    model.set_gen_to_subid(gen_sub["sub_id"].values)
    model.set_load_to_subid(load_sub["sub_id"].values)
    model.set_storage_to_subid(batt_sub["sub_id"].values)
    model.set_shunt_to_subid(sh_sub["sub_id"].values)
    model.set_line_or_to_subid(lor_sub["sub_id"].values)
    model.set_line_ex_to_subid(lex_sub["sub_id"].values)
    model.set_trafo_hv_to_subid(tor_sub["sub_id"].values)
    model.set_trafo_lv_to_subid(tex_sub["sub_id"].values)
    if not return_sub_id:
        return model
    else:
        return model, (gen_sub, load_sub, (lor_sub, tor_sub), (lex_sub, tex_sub), batt_sub, sh_sub, hvdc_sub_from_id, hvdc_sub_to_id)
