# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy
import numpy as np
import pypowsybl as pypo
# import pypowsybl.loadflow as lf

from lightsim2grid_cpp import GridModel


def init(net : pypo.network,
         gen_slack_id: int = None,
         slack_bus_id: int = None,
         sn_mva = 100.,
         sort_index=True,
         f_hz = 50.):
    model = GridModel()
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
    model._ls_to_orig = 1 * bus_df_orig["bus_id"].values
    voltage_levels = net.get_voltage_levels()
    model.set_sn_mva(sn_mva)
    model.set_init_vm_pu(1.06)
    model.init_bus(voltage_levels.loc[bus_df["voltage_level_id"].values]["nominal_v"].values,
                   0, 0  # unused
                   )
        
    # do the generators
    if sort_index:
        df_gen = net.get_generators().sort_index()
    else:
        df_gen = net.get_generators()
    # to handle encoding in 32 bits and overflow when "splitting" the Q values among 
    min_q = df_gen["min_q"].values.astype(np.float32)
    max_q = df_gen["max_q"].values.astype(np.float32)
    min_q[~np.isfinite(min_q)] = np.finfo(np.float32).min / 2. + 1.
    max_q[~np.isfinite(max_q)] = np.finfo(np.float32).max / 2. - 1.
    model.init_generators(df_gen["target_p"].values,
                          df_gen["target_v"].values / voltage_levels.loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
                          min_q,
                          max_q,
                          1 * bus_df.loc[df_gen["bus_id"].values]["bus_id"].values
                          )
    # TODO dist slack
    if gen_slack_id is not None:
        model.add_gen_slackbus(gen_slack_id, 1.)
    elif slack_bus_id is not None:
        gen_bus = np.array([el.bus_id for el in model.get_generators()])
        gen_is_conn_slack = gen_bus == model._ls_to_orig[slack_bus_id]
        nb_conn = gen_is_conn_slack.sum()
        if nb_conn == 0:
            raise RuntimeError(f"There is no generator connected to bus {slack_bus_id}. It cannot be the slack")
        for gen_id, is_slack in enumerate(gen_is_conn_slack):
            if is_slack:
                model.add_gen_slackbus(gen_id, 1. / nb_conn)   
    else:
        model.add_gen_slackbus(0, 1.)
        
    # for loads
    if sort_index:
        df_load = net.get_loads().sort_index()
    else:
        df_load = net.get_loads()
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     1 * bus_df.loc[df_load["bus_id"].values]["bus_id"].values
                     )
    
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
    model.init_powerlines_full(line_r,
                               line_x,
                               line_h_or,
                               line_h_ex,
                               1 * bus_df.loc[df_line["bus1_id"].values]["bus_id"].values,
                               1 * bus_df.loc[df_line["bus2_id"].values]["bus_id"].values
                              )
            
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
    model.init_trafo(df_trafo["r"].values / trafo_to_pu,
                     df_trafo["x"].values / trafo_to_pu,
                     2.*(1j*df_trafo["g"].values + df_trafo["b"].values) * trafo_to_pu,
                     tap_step_pct,
                     tap_pos,
                     shift_,
                     is_tap_hv_side,
                     1 * bus_df.loc[df_trafo["bus1_id"].values]["bus_id"].values, # TODO do I need to change hv / lv
                     1 * bus_df.loc[df_trafo["bus2_id"].values]["bus_id"].values)    
    
    # for shunt
    if sort_index:
        df_shunt = net.get_shunt_compensators().sort_index()
    else:
        df_shunt = net.get_shunt_compensators()
        
    is_on = copy.deepcopy(df_shunt["connected"])
    if (~is_on).any():
        df_shunt["connected"] = True
        net.update_shunt_compensators(df_shunt[["connected"]])
        if sort_index:
            df_shunt = net.get_shunt_compensators().sort_index()
        else:
            df_shunt = net.get_shunt_compensators()
        df_shunt["connected"] = is_on
        net.update_shunt_compensators(df_shunt[["connected"]])
        
    shunt_kv = voltage_levels.loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(-df_shunt["g"].values * shunt_kv**2,
                     -df_shunt["b"].values * shunt_kv**2,
                     1 * bus_df.loc[df_shunt["bus_id"].values]["bus_id"].values
                    )
    for shunt_id, conn in enumerate(is_on):
        if not conn:
           model.deactivate_shunt(shunt_id) 
           
    # for hvdc (TODO not tested yet)
    df_dc = net.get_hvdc_lines().sort_index()
    df_sations = net.get_vsc_converter_stations().sort_index()
    bus_from_id = df_sations.loc[df_dc["converter_station1_id"].values]["bus_id"].values
    bus_to_id = df_sations.loc[df_dc["converter_station2_id"].values]["bus_id"].values
    loss_percent = np.zeros(df_dc.shape[0])  # TODO 
    loss_mw = np.zeros(df_dc.shape[0])  # TODO
    model.init_dclines(bus_df.loc[bus_from_id]["bus_id"].values,
                       bus_df.loc[bus_to_id]["bus_id"].values,
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
    
    # storage units  (TODO not tested yet)
    if sort_index:
        df_batt = net.get_batteries().sort_index()
    else:
        df_batt = net.get_batteries()
    model.init_storages(df_batt["target_p"].values,
                        df_batt["target_q"].values,
                        1 * bus_df.loc[df_batt["bus_id"].values]["bus_id"].values
                        )
    
    # TODO
    # sgen => regular gen (from net.get_generators()) with voltage_regulator off TODO 
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
    if net.get_phase_tap_changers().shape[0] > 0:
        raise RuntimeError("Impossible currently to init a grid with tap changers at the moment.")
    return model
