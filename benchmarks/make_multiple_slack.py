from numpy.lib.npyio import load
import pandapower as pp
import numpy as np
import pandapower.networks as pn
import pdb
import copy
import tempfile
from lightsim2grid.initGridModel import init

# perf
# net = pn.case300()
# ls_grid_single = init(net)
# ls_grid_single.deactivate_result_computation()
# V = np.ones(net.bus.shape[0], dtype=np.complex_)
# # Vdc = ls_grid_single.dc_pf(copy.deepcopy(V), max_it, tol)
# ls_grid_single.reactivate_result_computation()
# V = ls_grid_single.ac_pf(V, 10, 1e-8)
# import sys
# sys.exit()

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
net = pp.create_empty_network("case14_custom", sn_mva=1.0 * case.sn_mva, f_hz= 1.0 * case.f_hz)
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
net.gen["slack_weight"][[id_ref_slack]] = 1.0

if False:  # test pass now
    print()
    print()
    print()
    print()
    print("With only one slack (case 2) ")
    pp.runpp(case, init="flat")
    print("\n")
    print("Real testcase")
    pp.runpp(net, init="flat")

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
    V = np.ones(nb_bus_total, dtype=np.complex_)
    # Vdc = ls_grid_single.dc_pf(copy.deepcopy(V), max_it, tol)
    ls_grid_single.reactivate_result_computation()
    V = ls_grid_single.ac_pf(V, max_it, tol)
    real_J = np.load("J_dist_slack_one_slack_firstIter.npy")  # from pandapower
    real_J_ls = np.load(file="J_ref_single_slack_ls_firstIter.npy")  # from ls single slack
    my_J = ls_grid_single.get_J()
    my_first_row = my_J[0].todense()
    ref_first_row = real_J[0]
    ref_first_row_ls = real_J_ls[0]
    Ybus_me = ls_grid_single.get_Ybus()
    Ybus_ref = net._ppc["internal"]["Ybus"]
    assert np.abs((Ybus_me - Ybus_ref).todense()).max() <= 1e-6, "wrong Ybus"
    # n_pv = 5
    # n_pq = 8
    pdb.set_trace()
    assert np.all(np.abs([el.res_p_or_mw for el in ls_grid_single.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)

print()
print()
print()
print()
print("Run PF with distributed slack (but still one slack)")
pp.runpp(net, init="flat", distributed_slack=True)
# now activate more slack bus
net.gen["slack"][[1]] = True
if "slack_weight" in net.gen:
    net.gen["slack_weight"][[1, id_ref_slack]] = 0.5
pp.rundcpp(net)

if False:
    # run the powerflow without distributed_slack=True
    print("Multiple slack, but nothing change in runpp (case 3)")
    pp.runpp(net, init="flat")
    gen_without_dist_slack = copy.deepcopy(net.res_gen)
    print(gen_without_dist_slack["p_mw"])

# run the powerflow with distributed_slack=True (same weights)
print("Multiple slack, and distributed_slack=True (case 4)")
pp.runpp(net, init="flat", distributed_slack=True)
print(net.res_gen["p_mw"])

ls_grid = init(net)
nb_bus_total = 14
# ls_grid.deactivate_result_computation()
V = np.ones(nb_bus_total, dtype=np.complex_)  # * ls_grid.get_init_vm_pu()
# Vdc = ls_grid.dc_pf(copy.deepcopy(V), max_it, tol)
# ls_grid.reactivate_result_computation()
V = ls_grid.ac_pf(V, max_it, tol)
my_ref = np.where(np.angle(V) == 0.)[0][0]
V_pp = net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  net.res_bus["va_degree"].values)
V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
assert np.abs(V_pp - V).max() <= 1e-6
assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)
assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - net.res_gen["q_mvar"].values) <= 1e-6)
assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - net.res_gen["p_mw"].values) <= 1e-6)

load_p = ls_grid.get_loads_res()[0]
gen_p = ls_grid.get_gen_res()[0]
p_or = ls_grid.get_lineor_res()[0]
p_ex = ls_grid.get_lineex_res()[0]
p_hv = ls_grid.get_trafohv_res()[0]
p_lv = ls_grid.get_trafolv_res()[0]
Sbus = ls_grid.get_Sbus()
Ybus = ls_grid.get_Ybus()
mis = V * (Ybus * V).conjugate() - Sbus
mis_pp = V_pp * (Ybus * V_pp).conjugate() - Sbus

# bus id 2
# gen 1
# load 1
# l2 ex
# l5 or
assert abs(-gen_p[1] + load_p[1] + p_ex[2] + p_or[5]) <= 1e-6

# bus id 0
# l0 or
# l1 or
# gen 4
assert abs(-gen_p[4] + p_or[0] + p_or[1]) <= 1e-6

# bus id 3
# l3 => to
# l5 => to
# l6 => from
# load 2
# t0 hv
# t1 hv

assert abs(load_p[2] + p_ex[3] + p_ex[5] + p_or[6] + p_hv[0] + p_hv[1]) <= 1e-6

# bus 5
# l7 or
# l8 or
# l9 or
# g2
# l4
# t2 lv
assert abs(p_or[7] + p_or[8] + p_or[9] - gen_p[2] + load_p[4] + p_lv[2]) <= 1e-6

mis = ls_grid.check_solution(V, False)
assert np.abs(mis).max() <= 1e-6, "error for lighsim2grid"
mis_pp = ls_grid.check_solution(V_pp, False)
assert np.abs(mis_pp).max() <= 1e-6, "error for pandapower"


### re run the powerflow but with different weights
net.gen["slack"][[1]] = True
if "slack_weight" in net.gen:
    net.gen["slack_weight"][[1, id_ref_slack]] = 0.7, 0.3
pp.rundcpp(net)

pp.runpp(net, init="flat", distributed_slack=True)
ls_grid2 = init(net)
V = np.ones(nb_bus_total, dtype=np.complex_) 
V = ls_grid2.ac_pf(V, max_it, tol)
my_ref = np.where(np.angle(V) == 0.)[0][0]
V_pp = net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  net.res_bus["va_degree"].values)
V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
assert np.abs(V_pp - V).max() <= 1e-6
assert np.all(np.abs([el.res_p_or_mw for el in ls_grid2.get_lines()] - net.res_line["p_from_mw"].values) <= 1e-6)

assert np.all(np.abs([el.res_q_mvar for el in ls_grid2.get_generators()] - net.res_gen["q_mvar"].values) <= 1e-6)
assert np.all(np.abs([el.res_p_mw for el in ls_grid2.get_generators()] - net.res_gen["p_mw"].values) <= 1e-6)

# TODO SLACK: hardest test 1: multiple slack gen at the same bus
# TODO SLACK: hardest test 2: at a bus, have both slack gen and non slack gen
# TODO SLACK: measure the performance of all this