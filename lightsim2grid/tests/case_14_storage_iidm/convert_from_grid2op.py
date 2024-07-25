# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import pandapower as pp
import pypowsybl as pypo
import pandas as pd
import numpy as np
import pypowsybl.network as pypo_n
import grid2op

from lightsim2grid.gridmodel import init_from_pandapower
from lightsim2grid.gridmodel import init_from_pypowsybl



def get_voltage_level_id(el, nb_dig_bus):
    return f'vl_{f"{el}".zfill(nb_dig_bus)}'


def get_bus_id(el, nb_dig_bus):
    return f'bus_{f"{el}".zfill(nb_dig_bus)}'

env = grid2op.make("educ_case14_storage", test=True)
pp_grid = env.backend._grid
pypo_grid = pypo_n.create_empty()  # network_id='educ_case14_storage'

# add substation NO BECAUSE pypowsybl call "site" substations
for el in range(env.n_sub):
    pypo_grid.create_substations(id=f"sub_{el}", name=f"sub_{el}")

# pypowsybl._pypowsybl.PyPowsyblError: Could not create transformer 3_6_0: both voltage ids must be on the same substation
corresp_sub = {6 : 3, 8 : 3, 5: 4, 7 : 3}

nb_dig_bus = len(str(env.n_sub))
# then add voltage levels
for el in range(env.n_sub):
    nominal_v = pp_grid.bus.iloc[el]["vn_kv"]
    pypo_grid.create_voltage_levels(id=get_voltage_level_id(el, nb_dig_bus),
                                    substation_id=f"sub_{el}" if not el in corresp_sub else f"sub_{corresp_sub[el]}",
                                    name=f"vl_{el}",
                                    topology_kind="BUS_BREAKER",
                                    nominal_v=nominal_v,
                                    low_voltage_limit= 0.90 * nominal_v,
                                    high_voltage_limit= 1.1 * nominal_v)

# add the buses (not modified grid)
for el in range(env.n_sub):
    pypo_grid.create_buses(id=get_bus_id(el, nb_dig_bus),
                           voltage_level_id=get_voltage_level_id(el, nb_dig_bus))

# add loads
nb_dig = len(str(env.n_load))
for el in range(env.n_load):
    el_pp = pp_grid.load.iloc[el]
    pypo_grid.create_loads(id=f'load_{f"{el}".zfill(nb_dig)}',
                           name=f"load_{el_pp['bus']}_{el}",
                           voltage_level_id=get_voltage_level_id(el_pp["bus"], nb_dig_bus),
                           bus_id=get_bus_id(el_pp["bus"], nb_dig_bus),
                           p0=el_pp["p_mw"],
                           q0=el_pp["q_mvar"],
                           )

# add generators
nb_dig = len(str(env.n_gen))
for el in range(env.n_gen):
    el_pp = pp_grid.gen.iloc[el]
    tg_v_pu = el_pp["vm_pu"]
    tg_v_kv = tg_v_pu * pp_grid.bus.loc[el_pp['bus']]["vn_kv"]
    pypo_grid.create_generators(id=f'gen_{f"{el}".zfill(nb_dig)}',
                                voltage_level_id=get_voltage_level_id(el_pp['bus'], nb_dig_bus),
                                bus_id=get_bus_id(el_pp["bus"], nb_dig_bus),
                                max_p=9999999.,
                                min_p=0.,
                                target_p=el_pp["p_mw"],
                                target_v=tg_v_kv,
                                voltage_regulator_on=True,
                                )
    
# add storage units
nb_dig = len(str(env.n_storage))
for el in range(env.n_storage):
    el_pp = pp_grid.storage.iloc[el]
    pypo_grid.create_batteries(id=f'storage_{f"{el}".zfill(nb_dig)}',
                               name=f"storage_{el_pp['bus']}_{el}",
                               voltage_level_id=get_voltage_level_id(el_pp['bus'], nb_dig_bus),
                               bus_id=get_bus_id(el_pp["bus"], nb_dig_bus),
                               min_p=-9999999.,
                               max_p=9999999.,
                               target_p=el_pp["p_mw"],
                               target_q=0.,
                               )
    
# add shunts
nb_dig = len(str(env.n_shunt))
shunt_df = pd.DataFrame({"id": [],
                         "voltage_level_id": [],
                         "bus_id" : [],
                         "name": [],
                         "model_type": [],
                         "section_count": []
                         })
model_df = pd.DataFrame({'id': [],
                         "g_per_section": [],
                         "b_per_section": [],
                         "max_section_count": []})
for el in range(env.n_shunt):
    el_pp = pp_grid.shunt.iloc[el]
    shunt_kv = el_pp["vn_kv"]
    id_ = f'shunt_{f"{el}".zfill(nb_dig)}'
    sh_dict = {"id": id_,
               "voltage_level_id": get_voltage_level_id(el_pp['bus'], nb_dig_bus),
               "bus_id": get_bus_id(el_pp["bus"], nb_dig_bus),
               "name": f"shunt_{el_pp['bus']}_{el}",
               "model_type": "LINEAR",
               "section_count": 1
               }
    mod_dict = {"id": id_,
                "g_per_section": -el_pp["p_mw"] / shunt_kv**2,
                "b_per_section": -el_pp["q_mvar"] / shunt_kv**2,
                "max_section_count": 1}
    shunt_df = pd.concat((shunt_df, pd.DataFrame({k: [v] for k, v in sh_dict.items()})))
    model_df = pd.concat((model_df, pd.DataFrame({k: [v] for k, v in mod_dict.items()})))
shunt_df["id"] = shunt_df["id"].astype(str)
shunt_df.set_index("id", inplace=True)
shunt_df["section_count"] = shunt_df["section_count"].astype(int)
model_df["id"] = model_df["id"].astype(str)
model_df.set_index("id", inplace=True)
model_df["max_section_count"] = model_df["max_section_count"].astype(int)
pypo_grid.create_shunt_compensators(shunt_df, model_df)

# powerlines
f_hz = 50.
nb_dig = len(str(pp_grid.line.shape[0]))
for el in range(pp_grid.line.shape[0]):
    el_pp = pp_grid.line.iloc[el]
    id_ = f'line_{f"{el}".zfill(nb_dig)}'
    nm_ = f"{el_pp['from_bus']}_{el_pp['to_bus']}_{el}"
    length_ = el_pp["length_km"]
    line_b = el_pp["c_nf_per_km"] * length_ * 2.0 * f_hz * np.pi * 1e-9
    pypo_grid.create_lines(id=id_,
                           name=nm_,
                           voltage_level1_id=get_voltage_level_id(el_pp['from_bus'], nb_dig_bus),
                           bus1_id=get_bus_id(el_pp["from_bus"], nb_dig_bus),
                           voltage_level2_id=get_voltage_level_id(el_pp['to_bus'], nb_dig_bus),
                           bus2_id=get_bus_id(el_pp["to_bus"], nb_dig_bus),
                           r=el_pp["r_ohm_per_km"] * length_,
                           x=el_pp["x_ohm_per_km"] * length_,
                           g1=0.5 * el_pp["g_us_per_km"] * length_ * 1e-6,  # TODO
                           g2=0.5 * el_pp["g_us_per_km"] * length_ * 1e-6,  # TODO
                           b1=0.5 * line_b,  # TODO
                           b2=0.5 * line_b,  # TODO
                           )

# transformers
nb_dig = len(str(pp_grid.trafo.shape[0]))
from pandapower.build_branch import _calc_tap_from_dataframe, _wye_delta, _calc_r_x_y_from_dataframe
net = {"_options": {"calculate_voltage_angles": True, "mode": "pf", "trafo_model": "t"}}
ppc = {}
sequence = 1
vn_trafo_hv, vn_trafo_lv, shift = _calc_tap_from_dataframe(net, pp_grid.trafo)
vn_lv = pp_grid.trafo["vn_lv_kv"].values / pp_grid.bus.loc[pp_grid.trafo["lv_bus"]]["vn_kv"].values
r, x, y = _calc_r_x_y_from_dataframe(pp_grid, pp_grid.trafo, vn_trafo_lv, vn_lv, ppc, sequence=sequence)
sn_mva = pp_grid.sn_mva
for el in range(pp_grid.trafo.shape[0]):
    el_pp = pp_grid.trafo.iloc[el]
    id_ = f'trafo_{f"{el}".zfill(nb_dig)}'
    nm_ = f"{el_pp['hv_bus']}_{el_pp['lv_bus']}_{el}"
    assert el_pp["tap_side"] == "hv" or el_pp["tap_side"] is None or not el_pp["tap_side"], f'{el_pp["tap_side"]} is not hv or None'
    # TODO I suppose here that tap is hv side !!!
    # otherwise reverse the tap here
    vn_kv_2 = pp_grid.bus.loc[el_pp["lv_bus"]]["vn_kv"]  
    vn_kv_1 = pp_grid.bus.loc[el_pp["hv_bus"]]["vn_kv"]
    if np.isfinite(el_pp["tap_step_percent"]):
        vn_kv_1 *= (1. + el_pp["tap_pos"] * el_pp["tap_step_percent"] / 100.)
    
    # trafo r and trafo x
    vn_lv_kv = el_pp["vn_lv_kv"]
    # parallel = 1
    # sn_mva = 1.
    # sn_trafo_mva = 1.
    # vk_percent = el_pp["vk_percent"]
    # vkr_percent = el_pp["vkr_percent"]
    # vn_lv = el_pp["vn_lv_kv"]
    # sn = 1.
    # tap_lv = np.square(vn_trafo_lv[el] / vn_lv) * sn_mva
    # parallel = 1

    # # r and x
    # z_sc = vk_percent / 100. / sn_trafo_mva * tap_lv
    # trafo_r = vkr_percent / 100. / sn_trafo_mva * tap_lv
    # trafo_x = np.sign(z_sc) * np.sqrt((z_sc ** 2 - trafo_r ** 2).astype(float))
    
    # # y which is needed for b and g
    # baseR = np.square(vn_lv) / sn_mva
    # pfe = el_pp["pfe_kw"] * 1e-3
    # vnl_squared = vn_lv_kv ** 2
    # b_real = pfe / vnl_squared * baseR
    # i0 = el_pp["i0_percent"]
    # b_img = (i0 / 100. * sn) ** 2 - pfe ** 2
    # if b_img < 0:
    #     b_img = 0.
    # b_img = np.sqrt(b_img) * baseR / vnl_squared
    # y = - b_real * 1j - b_img * np.sign(i0)
    
    # # then convert star to pi
    # trafo_r, trafo_x, y =_wye_delta(np.array([trafo_r]), np.array([trafo_x]), np.array([y]))
    trafo_to_pu = 1. / sn_mva
    trafo_r = r[el]
    trafo_x = x[el]
    
    pypo_grid.create_2_windings_transformers(id=id_,
                                             name=nm_,
                                             voltage_level1_id=get_voltage_level_id(el_pp['hv_bus'], nb_dig_bus),
                                             bus1_id=get_bus_id(el_pp["hv_bus"], nb_dig_bus),
                                             voltage_level2_id=get_voltage_level_id(el_pp['lv_bus'], nb_dig_bus),
                                             bus2_id=get_bus_id(el_pp["lv_bus"], nb_dig_bus),
                                             r=trafo_r * trafo_to_pu,
                                             x=trafo_x * trafo_to_pu,
                                             g=y[el].real / trafo_to_pu,
                                             b=y[el].imag / trafo_to_pu,
                                             rated_u1=vn_kv_1,
                                             rated_u2=vn_kv_2
                                             )


# check that grid are equals
ls_grid_pp = init_from_pandapower(pp_grid)
ls_grid_pypo = init_from_pypowsybl(pypo_grid,
                                   gen_slack_id=np.where(pp_grid.gen["slack"])[0],
                                   sn_mva=1.)

# check the elements are consistent
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_buses(), ls_grid_pypo.get_buses())):
    assert el_pp == el_pypo, f"error for {i}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_loads(), ls_grid_pypo.get_loads())):
    assert el_pp.bus_id == el_pypo.bus_id, f"error for {i}"
    assert el_pp.target_p_mw == el_pypo.target_p_mw, f"error for {i}"
    assert el_pp.target_q_mvar == el_pypo.target_q_mvar, f"error for {i}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_generators(), ls_grid_pypo.get_generators())):
    assert el_pp.bus_id == el_pypo.bus_id, f"error for {i}"
    assert el_pp.target_p_mw == el_pypo.target_p_mw, f"error for {i}"
    assert el_pp.target_vm_pu == el_pypo.target_vm_pu, f"error for {i}"
    assert el_pp.is_slack == el_pypo.is_slack, f"error for {i}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_storages(), ls_grid_pypo.get_storages())):
    assert el_pp.bus_id == el_pypo.bus_id, f"error for {i}"
    assert el_pp.target_p_mw == el_pypo.target_p_mw, f"error for {i}"
    assert el_pp.target_q_mvar == el_pypo.target_q_mvar, f"error for {i}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_shunts(), ls_grid_pypo.get_shunts())):
    assert el_pp.bus_id == el_pypo.bus_id, f"error for {i}"
    assert el_pp.target_p_mw == el_pypo.target_p_mw, f"error for {i}"
    assert el_pp.target_q_mvar == el_pypo.target_q_mvar, f"error for {i}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_lines(), ls_grid_pypo.get_lines())):
    assert np.allclose(el_pp.bus_or_id, el_pypo.bus_or_id), f"error for {i}"
    assert np.allclose(el_pp.bus_ex_id, el_pypo.bus_ex_id), f"error for {i}"
    assert np.allclose(el_pp.r_pu, el_pypo.r_pu), f"error for {i}: {el_pp.r_pu} vs {el_pypo.r_pu}"
    assert np.allclose(el_pp.x_pu, el_pypo.x_pu), f"error for {i}: {el_pp.x_pu} vs {el_pypo.x_pu}"
    assert np.allclose(el_pp.h_pu, el_pypo.h_pu), f"error for {i}: {el_pp.h_pu} vs {el_pypo.h_pu}"
    assert np.allclose(el_pp.h_or_pu, el_pypo.h_or_pu), f"error for {i}: {el_pp.h_or_pu} vs {el_pypo.h_or_pu}"
    assert np.allclose(el_pp.h_ex_pu, el_pypo.h_ex_pu), f"error for {i}: {el_pp.h_ex_pu} vs {el_pypo.h_ex_pu}"
    
for i, (el_pp, el_pypo) in enumerate(zip(ls_grid_pp.get_trafos(), ls_grid_pypo.get_trafos())):    
    assert np.allclose(el_pp.bus_lv_id, el_pypo.bus_lv_id), f"error for {i}"
    assert np.allclose(el_pp.bus_hv_id, el_pypo.bus_hv_id), f"error for {i}"
    assert np.allclose(el_pp.ratio, el_pypo.ratio), f"error for {i}: {el_pp.ratio} vs {el_pypo.ratio}"
    assert np.allclose(el_pp.shift_rad, el_pypo.shift_rad), f"error for {i}: {el_pp.shift_rad} vs {el_pypo.shift_rad}"
    assert np.allclose(el_pp.r_pu, el_pypo.r_pu), f"error for {i}: {el_pp.r_pu} vs {el_pypo.r_pu}"
    assert np.allclose(el_pp.x_pu, el_pypo.x_pu), f"error for {i}: {el_pp.x_pu} vs {el_pypo.x_pu}"
    assert np.allclose(el_pp.h_pu, el_pypo.h_pu), f"error for {i}: {el_pp.h_pu} vs {el_pypo.h_pu}"
    
# check powerflow is the same
V_init_pp = np.zeros(len(ls_grid_pp.get_buses()), dtype=complex) * 1.06
V_init_pypo = np.zeros(len(ls_grid_pypo.get_buses()), dtype=complex) * 1.06
V_pp = ls_grid_pp.ac_pf(1. * V_init_pp, 10, 1e-8)
V_pypo = ls_grid_pypo.ac_pf(1. * V_init_pypo, 10, 1e-8)
assert np.allclose(V_pp[:env.n_sub], V_pypo[:env.n_sub])
V_pp = ls_grid_pp.dc_pf(1. * V_init_pp, 10, 1e-8)
V_pypo = ls_grid_pypo.dc_pf(1. * V_init_pypo, 10, 1e-8)
assert np.allclose(V_pp[:env.n_sub], V_pypo[:env.n_sub])

# TODO dc_lines, static generators

# and now save the grid
pypo_grid.dump("grid.xiidm")
