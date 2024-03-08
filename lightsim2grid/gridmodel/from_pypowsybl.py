# Copyright (c) 2023-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
import pandas as pd
import pypowsybl as pypo

from lightsim2grid_cpp import GridModel


def _aux_get_bus(bus_df, df, conn_key="connected", bus_key="bus_id"):
    if df.shape[0] == 0:
        # no element of this type so no problem
        return np.zeros(0, dtype=int), np.ones(0, dtype=bool)
    # retrieve which elements are disconnected 
    mask_disco = ~df[conn_key]
    if mask_disco.all():
        raise RuntimeError("All element of the same type are disconnected, the init will not work.")
    first_el_co = np.where(~mask_disco.values)[0][0]
    # retrieve the bus where the element are
    tmp_bus_id = df[bus_key].copy()
    tmp_bus_id[mask_disco] = df.iloc[first_el_co][bus_key]  # assign a "random" bus to disco element
    bus_id = bus_df.loc[tmp_bus_id.values]["bus_id"].values
    # deactivate the element not on the main component
    # wrong_component = bus_df.loc[tmp_bus_id.values]["connected_component"].values != 0
    # mask_disco[wrong_component] = True
    # assign bus -1 to disconnected elements
    bus_id[mask_disco] = -1
    return bus_id, mask_disco.values


def init(net : pypo.network,
         gen_slack_id: int = None,
         slack_bus_id: int = None,
         sn_mva = 100.,
         sort_index=True,
         f_hz = 50.,  # unused
         only_main_component=True,
         return_sub_id=False):
    model = GridModel()
    # model.set_f_hz(f_hz)
    
    # for substation
    # network.get_voltage_levels()["substation_id"]
    # network.get_substations()
    # network.get_busbar_sections()
    
    if gen_slack_id is not None and slack_bus_id is not None:
        raise RuntimeError("Impossible to intialize a grid with both gen_slack_id and slack_bus_id")
    
    # assign unique id to the buses
    bus_df_orig = net.get_buses()
    if sort_index:
        bus_df = bus_df_orig.sort_index()
    else:
        bus_df = bus_df_orig
    bus_df["bus_id"] = np.arange(bus_df.shape[0])
    bus_df_orig["bus_id"] = bus_df.loc[bus_df_orig.index]["bus_id"]
    voltage_levels = net.get_voltage_levels()
    model.set_sn_mva(sn_mva)
    model.set_init_vm_pu(1.06)
    model.init_bus(voltage_levels.loc[bus_df["voltage_level_id"].values]["nominal_v"].values,
                   0, 0  # unused
                   )
    model._orig_to_ls = 1 * bus_df_orig["bus_id"].values
    
    # do the generators
    if sort_index:
        df_gen = net.get_generators().sort_index()
    else:
        df_gen = net.get_generators()
    # to handle encoding in 32 bits and overflow when "splitting" the Q values among 
    min_q = df_gen["min_q"].values.astype(np.float32)
    max_q = df_gen["max_q"].values.astype(np.float32)
    min_q[~np.isfinite(min_q)] = np.finfo(np.float32).min * 0.5 + 1.
    max_q[~np.isfinite(max_q)] = np.finfo(np.float32).max * 0.5 - 1.
    gen_bus, gen_disco = _aux_get_bus(bus_df, df_gen)
    model.init_generators_full(df_gen["target_p"].values,
                               df_gen["target_v"].values / voltage_levels.loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
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
    load_bus, load_disco = _aux_get_bus(bus_df, df_load)
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
    branch_from_kv = voltage_levels.loc[df_line["voltage_level1_id"].values]["nominal_v"].values
    branch_to_kv = voltage_levels.loc[df_line["voltage_level2_id"].values]["nominal_v"].values
    
    # only valid for lines with same voltages at both side...
    # branch_from_pu = branch_from_kv * branch_from_kv / sn_mva
    # line_r = df_line["r"].values / branch_from_pu
    # line_x = df_line["x"].values / branch_from_pu  
    # line_h_or = (1j*df_line["g1"].values + df_line["b1"].values) * branch_from_pu
    # line_h_ex = (1j*df_line["g2"].values + df_line["b2"].values) * branch_from_pu
    # real per unit conversion 
    # see https://github.com/powsybl/pypowsybl/issues/642
    # see https://github.com/powsybl/powsybl-core/blob/266442cbbd84f630acf786018618eaa3d496c6ba/ieee-cdf/ieee-cdf-converter/src/main/java/com/powsybl/ieeecdf/converter/IeeeCdfImporter.java#L347
    # for right formula
    v1 = branch_from_kv
    v2 = branch_to_kv
    line_r = sn_mva *  df_line["r"].values / v1 / v2
    line_x = sn_mva *  df_line["x"].values / v1 / v2
    tmp_ = np.reciprocal(df_line["r"].values + 1j*df_line["x"].values)
    b1 = df_line["b1"].values * v1*v1/sn_mva + (v1-v2)*tmp_.imag*v1/sn_mva
    b2 = df_line["b2"].values * v2*v2/sn_mva + (v2-v1)*tmp_.imag*v2/sn_mva
    g1 = df_line["g1"].values * v1*v1/sn_mva + (v1-v2)*tmp_.real*v1/sn_mva
    g2 = df_line["g2"].values * v2*v2/sn_mva + (v2-v1)*tmp_.real*v2/sn_mva
    line_h_or = (b1 + 1j * g1)
    line_h_ex = (b2 + 1j * g2)
    lor_bus, lor_disco = _aux_get_bus(bus_df, df_line, conn_key="connected1", bus_key="bus1_id")
    lex_bus, lex_disco = _aux_get_bus(bus_df, df_line, conn_key="connected2", bus_key="bus2_id")
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
    if sort_index:
        df_trafo = net.get_2_windings_transformers().sort_index()
    else:
        df_trafo = net.get_2_windings_transformers()
    # TODO net.get_ratio_tap_changers()
    # TODO net.get_phase_tap_changers()
    shift_ = np.zeros(df_trafo.shape[0])
    tap_pos = 1.0 * shift_
    is_tap_hv_side = np.ones(df_trafo.shape[0], dtype=bool)  # TODO
    
    # per unit
    trafo_from_kv = voltage_levels.loc[df_trafo["voltage_level1_id"].values]["nominal_v"].values
    trafo_to_kv = voltage_levels.loc[df_trafo["voltage_level2_id"].values]["nominal_v"].values
    trafo_to_pu = trafo_to_kv * trafo_to_kv / sn_mva
    # tap
    tap_step_pct = (df_trafo["rated_u1"] / trafo_from_kv - 1.) * 100.
    has_tap = tap_step_pct != 0.
    tap_pos[has_tap] += 1
    tap_step_pct[~has_tap] = 1.0  # or any other values...
    tor_bus, tor_disco = _aux_get_bus(bus_df, df_trafo, conn_key="connected1", bus_key="bus1_id")
    tex_bus, tex_disco = _aux_get_bus(bus_df, df_trafo, conn_key="connected2", bus_key="bus2_id")
    model.init_trafo(df_trafo["r"].values / trafo_to_pu,
                     df_trafo["x"].values / trafo_to_pu,
                     2.*(1j*df_trafo["g"].values + df_trafo["b"].values) * trafo_to_pu,
                     tap_step_pct,
                     tap_pos,
                     shift_,
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
        
    sh_bus, sh_disco = _aux_get_bus(bus_df, df_shunt)    
    shunt_kv = voltage_levels.loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(-df_shunt["g"].values * shunt_kv**2,
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
    # bus_from_id = df_sations.loc[df_dc["converter_station1_id"].values]["bus_id"].values
    # bus_to_id = df_sations.loc[df_dc["converter_station2_id"].values]["bus_id"].values
    hvdc_bus_from_id, hvdc_from_disco = _aux_get_bus(bus_df, df_sations.loc[df_dc["converter_station1_id"].values]) 
    hvdc_bus_to_id, hvdc_to_disco = _aux_get_bus(bus_df, df_sations.loc[df_dc["converter_station2_id"].values]) 
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
    batt_bus, batt_disco = _aux_get_bus(bus_df, df_batt)
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
    elif gen_slack_id is not None:
        if slack_bus_id is not None:
            raise RuntimeError(f"You provided both gen_slack_id and slack_bus_id which is not possible.")
        model.add_gen_slackbus(gen_slack_id, 1.)
    elif slack_bus_id is not None:
        gen_bus = np.array([el.bus_id for el in model.get_generators()])
        gen_is_conn_slack = gen_bus == model._orig_to_ls[slack_bus_id]
        nb_conn = gen_is_conn_slack.sum()
        if nb_conn == 0:
            raise RuntimeError(f"There is no generator connected to bus {slack_bus_id}. It cannot be the slack")
        for gen_id, is_slack in enumerate(gen_is_conn_slack):
            if is_slack:
                model.add_gen_slackbus(gen_id, 1. / nb_conn)    
    
    # TODO
    # sgen => regular gen (from net.get_generators()) with voltage_regulator off TODO 
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
    if net.get_phase_tap_changers().shape[0] > 0:
        warnings.warn("There are tap changers in the iidm grid which are not taken "
                      "into account in the lightsim2grid at the moment. "
                      "NB: lightsim2grid gridmodel can handle tap changer, it is just not "
                      "handled by the 'from_pypowsybl` function at the moment.")
        
    # and now deactivate all elements and nodes not in the main component
    if only_main_component:
        model.consider_only_main_component()
    if not return_sub_id:
        # for backward compatibility
        return model
    else:
        # voltage_level_id is kind of what I call "substation" in grid2op
        vl_unique = bus_df["voltage_level_id"].unique()
        sub_df = pd.DataFrame(index=np.sort(vl_unique), data={"sub_id": np.arange(vl_unique.size)})
        buses_sub_id = pd.merge(left=bus_df, right=sub_df, how="left", left_on="voltage_level_id", right_index=True)[["bus_id", "sub_id"]]
        gen_sub = pd.merge(left=df_gen, right=sub_df, how="left", left_on="voltage_level_id", right_index=True)[["sub_id"]]
        load_sub = pd.merge(left=df_load, right=sub_df, how="left", left_on="voltage_level_id", right_index=True)[["sub_id"]]
        lor_sub = pd.merge(left=df_line, right=sub_df, how="left", left_on="voltage_level1_id", right_index=True)[["sub_id"]]
        lex_sub = pd.merge(left=df_line, right=sub_df, how="left", left_on="voltage_level2_id", right_index=True)[["sub_id"]]
        tor_sub = pd.merge(left=df_trafo, right=sub_df, how="left", left_on="voltage_level1_id", right_index=True)[["sub_id"]]
        tex_sub = pd.merge(left=df_trafo, right=sub_df, how="left", left_on="voltage_level2_id", right_index=True)[["sub_id"]]
        batt_sub = pd.merge(left=df_batt, right=sub_df, how="left", left_on="voltage_level_id", right_index=True)[["sub_id"]]
        sh_sub = pd.merge(left=df_shunt, right=sub_df, how="left", left_on="voltage_level_id", right_index=True)[["sub_id"]]
        hvdc_vl_info = pd.DataFrame(index=df_dc.index,
                                    data={"voltage_level1_id": df_sations.loc[df_dc["converter_station1_id"].values]["voltage_level_id"].values,
                                          "voltage_level2_id": df_sations.loc[df_dc["converter_station2_id"].values]["voltage_level_id"].values
                                          })
        hvdc_sub_from_id = pd.merge(left=hvdc_vl_info, right=sub_df, how="left", left_on="voltage_level1_id", right_index=True)[["sub_id"]]
        hvdc_sub_to_id = pd.merge(left=hvdc_vl_info, right=sub_df, how="left", left_on="voltage_level2_id", right_index=True)[["sub_id"]]
        return model, (buses_sub_id, gen_sub, load_sub, (lor_sub, tor_sub), (lex_sub, tex_sub), batt_sub, sh_sub, hvdc_sub_from_id, hvdc_sub_to_id)
