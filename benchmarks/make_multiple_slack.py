import pandapower as pp
import numpy as np
import pandapower.networks as pn
import pdb
import copy
import tempfile
from lightsim2grid.initGridModel import init


# LF PARAMETERS
max_it = 10
tol = 1e-8

# retrieve the case14 and remove the "ext_grid" => put a generator as slack bus instead
case = pn.case14()
pp.runpp(case)  # forced to do it to retrieve the power of the slack bus
# adding a slack bus
var_pp = ["bus", 'in_service', "name", "max_p_mw", "min_p_mw", "max_q_mvar", "min_q_mvar", "slack_weight"]
if "slack_weight" not in case.ext_grid:
    var_pp = var_pp[:-1]

slack_bus_gen_id_ppc = 0  # checked manually
pp.create_gen(case,
              **case.ext_grid[var_pp].iloc[0],
              p_mw=case._ppc['gen'][slack_bus_gen_id_ppc, 1],
              slack=True)
# "deactivate" the "ext_grid" 
case.ext_grid["in_service"][0] = False
pp.runpp(case)

# now create a copy of it, by removing the ext_grid completely (to be sure)
net = pp.create_empty_network("case14_custom", sn_mva=case.sn_mva, f_hz=case.f_hz)
# create bus
for i in range(case.bus.shape[0]):
    pp.create_bus(net, **case.bus.iloc[i])

# create lines
var_line = ["from_bus", "to_bus", "length_km", "r_ohm_per_km", "x_ohm_per_km", "c_nf_per_km", "max_i_ka", "name",
            "index", "type", "geodata", "in_service", "df", "parallel", "g_us_per_km",
            "max_loading_percent", "alpha", "temperature_degree_celsius"]
var_line = [el for el in var_line if el in case.line]
for i in range(case.line.shape[0]):
    pp.create_line_from_parameters(net, **case.line[var_line].iloc[i])

# create trafos
var_trafo = ["hv_bus", "lv_bus", "sn_mva", "vn_hv_kv", "vn_lv_kv", "vkr_percent", "vk_percent", "pfe_kw", "i0_percent", "shift_degree",
             "tap_side", "tap_neutral", "tap_max", "tap_min", "tap_step_percent", "tap_step_degree", 
             "tap_pos", "tap_phase_shifter", "in_service", "name", "index", "max_loading_percent",
             "parallel", "df"]
var_trafo = [el for el in var_trafo if el in case.trafo]
for i in range(case.trafo.shape[0]):
    pp.create_transformer_from_parameters(net, **case.trafo[var_trafo].iloc[i])

# create shunts
var_shunt = ["bus", "q_mvar", "p_mw", "vn_kv", "step", "max_step", "name",
             "in_service", "index"]
var_shunt = [el for el in var_shunt if el in case.shunt]
for i in range(case.shunt.shape[0]):
    pp.create_shunt(net, **case.shunt[var_shunt].iloc[i])

# create loads
var_load = ["bus", "p_mw", "q_mvar", "const_z_percent", "const_i_percent", "sn_mva", "name",
            "scaling", "index", "in_service", "type", "max_p_mw", "min_p_mw", "max_q_mvar",
            "min_q_mvar", "controllable"]
var_load = [el for el in var_load if el in case.load]
for i in range(case.load.shape[0]):
    pp.create_load(net, **case.load[var_load].iloc[i])

# create gens
var_gen = ["bus", "p_mw", "vm_pu", "sn_mva", "name", "index", "max_q_mvar", "min_q_mvar", "min_p_mw",
            "max_p_mw", "scaling", "type", "slack", "controllable", "vn_kv", "xdss_pu", "rdss_pu",
            "cos_phi", "in_service"]
var_gen = [el for el in var_gen if el in case.gen]
for i in range(case.gen.shape[0]):
    pp.create_gen(net, **case.gen[var_gen].iloc[i])

id_ref_slack = net.gen.shape[0]-1  # initial generator added as the slack bus added
net.gen["min_p_mw"][[id_ref_slack]] = 0.
net.gen["max_p_mw"][[id_ref_slack]] = 300.
assert (net.gen["min_p_mw"] >= 0.).all()
assert (net.gen["min_p_mw"] <= 400.).all()

# check the powerflow match
print("Setpoint values")
print(net.gen["p_mw"])
if True:
    print()
    print()
    print()
    print()
    print("With only one slack (case 2) ")
    pp.runpp(case)
    pp.runpp(net)

    assert ((case.res_gen - net.res_gen) <= 1e-6).all().all()
    assert ((case.res_load - net.res_load) <= 1e-6).all().all()
    assert ((case.res_shunt - net.res_shunt) <= 1e-6).all().all()
    assert ((case.res_line - net.res_line) <= 1e-6).all().all()
    assert ((case.res_trafo - net.res_trafo) <= 1e-6).all().all()
    assert ((case.res_bus.iloc[:net.bus.shape[0]] - net.res_bus) <= 1e-6).all().all()

    gen_without_single_slack = copy.deepcopy(net.res_gen)
    print(gen_without_single_slack["p_mw"])

    ls_grid_single = init(net)
    nb_bus_total = 14
    ls_grid_single.deactivate_result_computation()
    V = np.ones(nb_bus_total, dtype=np.complex_) * ls_grid_single.get_init_vm_pu()
    Vdc = ls_grid_single.dc_pf(copy.deepcopy(V), max_it, tol)
    ls_grid_single.reactivate_result_computation()
    V = ls_grid_single.ac_pf(Vdc, max_it, tol)
    pdb.set_trace()
    assert np.all(np.abs([el.res_p_or_mw for el in ls_grid_single.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)

print()
print()
print()
print()
net.gen["slack_weight"][[id_ref_slack]] = 1.0
pp.runpp(net, init="dc", distributed_slack=True)
pdb.set_trace()
# now activate more slack bus
net.gen["slack"][[1]] = True
if "slack_weight" in net.gen:
    net.gen["slack_weight"][[1, id_ref_slack]] = 0.5
pp.rundcpp(net)

if False:
    # run the powerflow without distributed_slack=True
    print("Multiple slack, but nothing change in runpp (case 3)")
    pp.runpp(net, init="dc")
    gen_without_dist_slack = copy.deepcopy(net.res_gen)
    print(gen_without_dist_slack["p_mw"])

if False:
    ls_grid = init(net)
    nb_bus_total = 14
    ls_grid.deactivate_result_computation()
    V = np.ones(nb_bus_total, dtype=np.complex_) * ls_grid.get_init_vm_pu()
    Vdc = ls_grid.dc_pf(copy.deepcopy(V), max_it, tol)
    ls_grid.reactivate_result_computation()
    V = ls_grid.ac_pf(Vdc, max_it, tol)
    pdb.set_trace()


# run the powerflow with distributed_slack=True
print("Multiple slack, and distributed_slack=True (case 4)")
pp.runpp(net, init="dc", distributed_slack=True)
print(net.res_gen["p_mw"])
pdb.set_trace()

