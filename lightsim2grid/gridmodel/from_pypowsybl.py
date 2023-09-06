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

def init(net : pypo.network):
    model = GridModel()
    
    # initialize and use converters
    sn_mva_ = 100.  # TODO read from net
    f_hz = 50.   # TODO read from net
    
    # assign unique id to the buses
    bus_df = net.get_buses().copy()
    bus_df["bus_id"] = np.arange(bus_df.shape[0])
    model.init_bus(net.get_voltage_levels().loc[bus_df["voltage_level_id"].values]["nominal_v"].values,
                   0, 0  # unused
                   )
        
    # do the generators
    df_gen = net.get_generators()
    model.init_generators(df_gen["target_p"].values,
                          df_gen["target_v"].values / net.get_voltage_levels().loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
                          df_gen["min_q"].values,
                          df_gen["max_q"].values,
                          bus_df.loc[df_gen["bus_id"].values]["bus_id"].values
                          )
    # TODO slack
    model.add_gen_slackbus(0, 1.)
    
    # for loads
    df_load = net.get_loads()
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     bus_df.loc[df_load["bus_id"].values]["bus_id"].values
                     )
    
    # for lines
    df_line = net.get_lines()
    # TODO add g1 / b1 and g2 / b2 in lightsim2grid
    line_h = (1j*df_line["g1"].values + df_line["b1"].values + 1j*df_line["g2"].values + df_line["b2"].values)
    # per unit
    branch_from_kv = net.get_voltage_levels().loc[df_line["voltage_level1_id"].values]["nominal_v"].values
    branch_from_pu = branch_from_kv * branch_from_kv / sn_mva_
    model.init_powerlines(df_line["r"].values / branch_from_pu,
                          df_line["x"].values / branch_from_pu,
                          line_h * branch_from_pu,
                          bus_df.loc[df_line["bus1_id"].values]["bus_id"].values,
                          bus_df.loc[df_line["bus2_id"].values]["bus_id"].values
                         )
            
    # for trafo
    df_trafo = net.get_2_windings_transformers()
    # TODO net.get_ratio_tap_changers()
    # TODO net.get_phase_tap_changers()
    shift_ = np.zeros(df_trafo.shape[0])
    tap_step_pct = shift_
    tap_pos = shift_
    is_tap_hv_side = np.ones(df_trafo.shape[0], dtype=bool)  # TODO
    
    # per unit
    trafo_from_kv = net.get_voltage_levels().loc[df_trafo["voltage_level1_id"].values]["nominal_v"].values
    trafo_to_kv = net.get_voltage_levels().loc[df_trafo["voltage_level2_id"].values]["nominal_v"].values
    trafo_from_pu = trafo_from_kv * trafo_to_kv / sn_mva_
    model.init_trafo(df_trafo["r"].values / trafo_from_pu,
                     df_trafo["x"].values / trafo_from_pu,
                     2.*(df_trafo["g"].values + 1j*df_trafo["b"].values) * trafo_from_pu,
                     tap_step_pct,
                     tap_pos,
                     shift_,
                     is_tap_hv_side,
                     bus_df.loc[df_trafo["bus1_id"].values]["bus_id"].values, # TODO do I need to change hv / lv
                     bus_df.loc[df_trafo["bus2_id"].values]["bus_id"].values)    
    
    # for shunt
    df_shunt = net.get_shunt_compensators()
    model.init_shunt(df_shunt["g"].values,
                     df_shunt["b"].values,
                     bus_df.loc[df_shunt["bus_id"].values]["bus_id"].values
                    )
    
    # TODO
    # dcline
    # sgen
    # storage
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
    return model
