# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import pypowsybl as pypo
from lightsim2grid_cpp import GridModel


def init(net : pypo.network, gen_slack_id: int = None):
    model = GridModel()
    # for substation
    # network.get_voltage_levels()["substation_id"]
    # network.get_substations()
    # network.get_busbar_sections()
    
    
    # initialize and use converters
    sn_mva_ = 100.  # TODO read from net
    f_hz = 50.   # TODO read from net
    
    # assign unique id to the buses
    bus_df = net.get_buses().copy()
    bus_df["bus_id"] = np.arange(bus_df.shape[0])
    model.set_sn_mva(sn_mva_)
    model.set_init_vm_pu(1.06)
    model.init_bus(net.get_voltage_levels().loc[bus_df["voltage_level_id"].values]["nominal_v"].values,
                   0, 0  # unused
                   )
        
    # do the generators
    df_gen = net.get_generators()
    model.init_generators(df_gen["target_p"].values,
                          df_gen["target_v"].values / net.get_voltage_levels().loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
                          df_gen["min_q"].values,
                          df_gen["max_q"].values,
                          1 * bus_df.loc[df_gen["bus_id"].values]["bus_id"].values
                          )
    # TODO dist slack
    if gen_slack_id is None:
        model.add_gen_slackbus(0, 1.)
    else:
        model.add_gen_slackbus(gen_slack_id, 1.)
    
    # for loads
    df_load = net.get_loads()
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     1 * bus_df.loc[df_load["bus_id"].values]["bus_id"].values
                     )
    
    # for lines
    df_line = net.get_lines()
    # TODO add g1 / b1 and g2 / b2 in lightsim2grid
    # line_h = (1j*df_line["g1"].values + df_line["b1"].values + 1j*df_line["g2"].values + df_line["b2"].values)
    # per unit
    branch_from_kv = net.get_voltage_levels().loc[df_line["voltage_level1_id"].values]["nominal_v"].values
    branch_to_kv = net.get_voltage_levels().loc[df_line["voltage_level2_id"].values]["nominal_v"].values
    
    # only valid for lines with same voltages at both side...
    # branch_from_pu = branch_from_kv * branch_from_kv / sn_mva_
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
    line_r = sn_mva_ *  df_line["r"].values / v1 / v2
    line_x = sn_mva_ *  df_line["x"].values / v1 / v2
    tmp_ = np.reciprocal(df_line["r"].values + 1j*df_line["x"].values)
    b1 = df_line["b1"].values * v1*v1/sn_mva_ + (v1-v2)*tmp_.imag*v1/sn_mva_
    b2 = df_line["b2"].values * v2*v2/sn_mva_ + (v2-v1)*tmp_.imag*v2/sn_mva_
    g1 = df_line["g1"].values * v1*v1/sn_mva_ + (v1-v2)*tmp_.real*v1/sn_mva_
    g2 = df_line["g2"].values * v2*v2/sn_mva_ + (v2-v1)*tmp_.real*v2/sn_mva_
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
    df_trafo = net.get_2_windings_transformers()
    # TODO net.get_ratio_tap_changers()
    # TODO net.get_phase_tap_changers()
    shift_ = np.zeros(df_trafo.shape[0])
    tap_pos = 1.0 * shift_
    is_tap_hv_side = np.ones(df_trafo.shape[0], dtype=bool)  # TODO
    
    # per unit
    trafo_from_kv = net.get_voltage_levels().loc[df_trafo["voltage_level1_id"].values]["nominal_v"].values
    trafo_to_kv = net.get_voltage_levels().loc[df_trafo["voltage_level2_id"].values]["nominal_v"].values
    trafo_to_pu = trafo_to_kv * trafo_to_kv / sn_mva_
    # tap
    tap_step_pct = (df_trafo["rated_u1"] / trafo_from_kv - 1.) * 100.
    tap_pos += 1
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
    df_shunt = net.get_shunt_compensators()
    shunt_kv = net.get_voltage_levels().loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(-df_shunt["g"].values * shunt_kv**2,
                     -df_shunt["b"].values * shunt_kv**2,
                     1 * bus_df.loc[df_shunt["bus_id"].values]["bus_id"].values
                    )
    
    # for hvdc (TODO not tested yet)
    df_dc = net.get_hvdc_lines()
    df_sations = net.get_vsc_converter_stations()
    bus_from_id = df_sations.loc[df_dc["converter_station1_id"].values]["bus_id"].values
    bus_to_id = df_sations.loc[df_dc["converter_station2_id"].values]["bus_id"].values
    loss_percent = np.zeros(df_dc.shape[0])  # TODO 
    loss_mw = np.zeros(df_dc.shape[0])  # TODO
    model.init_dclines(bus_df.loc[bus_from_id]["bus_id"].values,
                       bus_df.loc[bus_to_id]["bus_id"].values,
                       df_dc["target_p"].values,
                       loss_percent,
                       loss_mw,
                       net.get_voltage_levels().loc[df_sations.loc[df_dc["converter_station1_id"].values]["voltage_level_id"].values]["nominal_v"].values,
                       net.get_voltage_levels().loc[df_sations.loc[df_dc["converter_station2_id"].values]["voltage_level_id"].values]["nominal_v"].values,
                       df_sations.loc[df_dc["converter_station1_id"].values]["min_q"].values,
                       df_sations.loc[df_dc["converter_station1_id"].values]["max_q"].values,
                       df_sations.loc[df_dc["converter_station2_id"].values]["min_q"].values,
                       df_sations.loc[df_dc["converter_station2_id"].values]["max_q"].values
                       )
    
    # storage units  (TODO not tested yet)
    df_batt = net.get_batteries()
    model.init_storages(df_batt["target_p"].values,
                        df_batt["target_q"].values,
                        1 * bus_df.loc[df_batt["bus_id"].values]["bus_id"].values
                        )
    
    # TODO
    # sgen
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
    return model
