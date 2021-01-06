import pandapower as pp
from lightsim2grid import LightSimBackend
import numpy as np
import copy
import scipy
import pandapower
import pdb
import sys
import warnings
try:
    from lightsim2grid.solver import KLUSolver
    ClassSolver = KLUSolver
except ImportError as exc_:
    from lightsim2grid.solver import SparseLUSolver
    ClassSolver = SparseLUSolver

import pandapower.networks as pn

# To create the json file
# case = pn.test_case6495rte()
# case_name = "test_case6495rte.json"
# pp.to_json(pn_net, case_name)

case_name = "test_case6495rte.json"
nb_sub = 6495
# case_name = "test_case2848rte.json"
# nb_sub = 2848
# case_name = "test_case118.json"
# nb_sub = 118
# case_name = "test_case14.json"
# nb_sub = 14
tol = 1e-4  # results are equal if they match up to tol

good_working_case = pn.case118()
good_working_case = pn.case14()

#### test if it works
print(f"Basic check for {case_name}")
real_init_file = pp.from_json(case_name)
backend = LightSimBackend()
backend.load_grid(case_name)
pp_net = backend.init_pp_backend._grid
# first i deactivate all slack bus in pp that are connected but not handled in ls
pp_net.ext_grid["in_service"].loc[:] = False
pp_net.ext_grid["in_service"].iloc[0] = True
conv = backend.runpf()
conv_pp = backend.init_pp_backend.runpf()

if not conv_pp:
    print("Error: pandapower do not converge, impossible to perform the necessary checks")
    sys.exit()
por_pp, qor_pp, vor_pp, aor_pp = copy.deepcopy(backend.init_pp_backend.lines_or_info())
pex_pp, qex_pp, vex_pp, aex_pp = copy.deepcopy(backend.init_pp_backend.lines_ex_info())

print("I- Check for divergence and equality of flows")
if conv and conv_pp:
    print("\t Info: Both converged")
    por_ls, qor_ls, vor_ls, aor_ls = backend.lines_or_info()
    test_ok = True
    max_mis = np.max(np.abs(por_ls - por_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: por do not match, maximum absolute error is {max_mis:.5f} MW")

    max_mis = np.max(np.abs(qor_ls - qor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr")
    max_mis = np.max(np.abs(vor_ls - vor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: vor do not match, maximum absolute error is {max_mis:.5f} kV")
    max_mis = np.max(np.abs(aor_ls - aor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: aor do not match, maximum absolute error is {max_mis:.5f} A")
    if test_ok:
        print(f"Info: Everything ok for the tolerance {tol} for {case_name}")
        # sys.exit()
elif not conv and conv_pp:
    print("\t Error: LS do not converge, but pandapower does")


#### test if it's the solver that cannot handle the jacobian
## in this case, i simply take the Ybus, Sbus and all from pandapower, fill it to lightsim
## and theck if results matches
print(f"II - Check for possible solver issues")
pp.runpp(backend.init_pp_backend._grid, v_debug=True)
v_tmp = backend.init_pp_backend._grid.res_bus["vm_pu"].values[:nb_sub] + 0j
v_tmp *= np.exp(1j * np.pi/180. * backend.init_pp_backend._grid.res_bus["va_degree"].values[:nb_sub])
v_tmp = np.concatenate((v_tmp, v_tmp))
backend._grid.ac_pf(v_tmp, 1000, 1e-5)

Y_pp = backend.init_pp_backend._grid._ppc["internal"]["Ybus"]
Sbus = backend.init_pp_backend._grid._ppc["internal"]["Sbus"]
pv_ = backend.init_pp_backend._grid._ppc["internal"]["pv"]
pq_ = backend.init_pp_backend._grid._ppc["internal"]["pq"]
max_iter = 10
tol_this = 1e-8
All_Vms = backend.init_pp_backend._grid._ppc["internal"]["Vm_it"]
AllVas = backend.init_pp_backend._grid._ppc["internal"]["Va_it"]

for index_V in range(All_Vms.shape[1]-1, -1, -1):
    # i check from easiest to hardest, so from the last iteartion of pandapower to the first iteration of pandapower
    # take the same V as pandapower
    V_init = All_Vms[:, index_V] * (np.cos(AllVas[:, index_V]) + 1j * np.sin(AllVas[:, index_V]))
    # V_init *= np.exp(1j * AllVas[:, 0])
    V_init_ref = copy.deepcopy(V_init)
    solver = ClassSolver()
    solver.solve(scipy.sparse.csc_matrix(Y_pp), V_init, Sbus, pv_, pq_, max_iter, tol_this)
    time_for_nr = solver.get_timers()[3]
    print(f"\t Info: Time to perform the NR for {case_name}: {1000.*time_for_nr:.2f}ms")
    error_va = np.abs(solver.get_Va() - np.angle(backend.init_pp_backend._grid._ppc["internal"]["V"]))
    test_ok = True
    if np.max(error_va) > tol:
        test_ok = False
        print(f"\t Error: VA do not match for iteration {index_V}, maximum absolute error "
              f"is {np.max(error_va):.5f} rad")
    error_vm = np.abs(np.abs(solver.get_Vm() - np.abs(backend.init_pp_backend._grid._ppc["internal"]["V"])))
    if np.max(error_vm) > tol:
        test_ok = False
        print(f"\t Error: VM do not match for iteration {index_V}, maximum absolute error "
              f"is {np.max(error_vm):.5f} pu")
    solver.reset()
    if test_ok:
        print(f"\t Info: solver ok if initialized with pandapower data from iteration "
              f"{index_V} / {int(All_Vms.shape[1]-1)}")


### Now i check for solver issue
print('III - Check the data conversion')
pp_vect_converter = backend.init_pp_backend._grid._pd2ppc_lookups["bus"][:nb_sub]
### first check the Sbus vector
print("1) Checking Sbus conversion")
pp_net = backend.init_pp_backend._grid
Sbus_pp = backend.init_pp_backend._grid._ppc["internal"]["Sbus"]
Sbus_pp_right_order = Sbus_pp[pp_vect_converter]
Sbus_me = backend._grid.get_Sbus()
test_ok = True
error_p = np.abs(np.real(Sbus_me) - np.real(Sbus_pp_right_order))
if np.max(error_p) > tol:
    test_ok = False
    print(f"\t Error: P do not match for Sbus, maximum absolute error is {np.max(error_p):.5f} MW")
    print(f"\t Error: significative difference for bus index (lightsim): {np.where(error_p > tol)[0]}")

error_q = np.abs(np.imag(Sbus_me) - np.imag(Sbus_pp_right_order))
if np.max(error_q) > tol:
    test_ok = False
    print(f"\t Error: Q do not match for Sbus, maximum absolute error is {np.max(error_q):.5f} MVAr")
    print(f"\t Error: significative difference for bus index (lightsim): {np.where(error_q > tol)[0]}")

if test_ok:
    print("\t Info: Sbus is correct")
else:
    print("\t Error: there is a problem with Sbus, entering debug")
    errors = np.where(error_p > tol)[0]
    np.sum(Sbus_me[errors])
    np.sum(Sbus_pp_right_order[errors])

    for bus_id_error in np.where(error_p > tol)[0]:
        bus_pp_orig = np.where(backend.init_pp_backend._grid._pd2ppc_lookups["bus"][:nb_sub] == bus_id_error)[0][0]
        gen_this_bus = pp_net.gen.loc[pp_net.gen["bus"] == bus_id_error]
        load_this_bus = pp_net.load.loc[pp_net.load["bus"] == bus_id_error]
        sgen_this_bus = pp_net.sgen.loc[pp_net.sgen["bus"] == bus_id_error]


### then check the Ybus matrix
print("2)  Checking Ybus conversion")
Y_me = backend._grid.get_Ybus()
Y_pp = backend.init_pp_backend._grid._ppc["internal"]["Ybus"]
Y_pp_right_order = Y_pp[pp_vect_converter.reshape(nb_sub, 1), pp_vect_converter.reshape(1, nb_sub)]

test_ok = True
error_p = np.abs(np.real(Y_me) - np.real(Y_pp_right_order))
if np.max(error_p) > tol:
    test_ok = False
    print(f"\t Error: P do not match for Ybus, maximum absolute error is {np.max(error_p):.5f}")
    diffs = error_p > tol
    diffs = diffs.todense()
    ind_row, ind_cols = np.where(diffs)
    rest = f"... {len(ind_row)} in total" if len(ind_row) >= 10 else ""
    print(f"\t Error: significative difference for bus index (lightsim): {ind_row[:10], ind_cols[:10]} "
          f"{rest}")
    where_max = np.where(error_p.todense() >= 100)

error_q = np.abs(np.imag(Y_me) - np.imag(Y_pp_right_order))
if np.max(error_q) > tol:
    test_ok = False
    print(f"\t Error: Q do not match for Ybus, maximum absolute error is {np.max(error_q):.5f}")
    diffs = error_q > tol
    diffs = diffs.todense()
    ind_row, ind_cols = np.where(diffs)
    rest = f"... {len(ind_row)} in total" if len(ind_row) >= 10 else ""
    print(f"\t Error: significative difference for bus index (lightsim): {ind_row[:10], ind_cols[:10]} "
          f"{rest}")
    where_max = np.where(error_p.todense() >= 100)

if test_ok:
    print("\t Info: Ybus is correct")
else:
    where_errors = np.where(diffs)
    nb_error = where_errors[0].shape[0]
    error_diag_coeff = where_errors[0] == where_errors[1]
    error_q_d = error_q.todense()
    indx_max = np.where(error_q_d == np.max(error_q_d))

    # test trafo
    # for issue_id in range(where_errors[0].shape[0]):
    #     bus_id_or_error = where_errors[0][issue_id]
    #     bus_id_ex_error = where_errors[1][issue_id]
    #     trafo_bug = pp_net.trafo.loc[(pp_net.trafo["hv_bus"] == bus_id_or_error) & (pp_net.trafo["lv_bus"] == bus_id_ex_error)]
    #     if trafo_bug.shape[0] == 1:
    #         print("hv id_or")
    #         import pdb
    #         pdb.set_trace()
    #     trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_id_or_error) & (pp_net.trafo["hv_bus"] == bus_id_ex_error)]
    #     if trafo_bug.shape[0] == 1:
    #         print("hv id_ex")
    #         import pdb
    #         pdb.set_trace()

    # error of side of tap
    if False:
        # 'test_case6495rte.json'
        bus_lv = 30
        bus_hv = 32
        trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_lv) & (pp_net.trafo["hv_bus"] == bus_hv)]
        me_val_lv_hv = Y_me[bus_lv, bus_hv]
        pp_val_lv_hv = Y_pp_right_order[bus_lv, bus_hv]
        me_val_hv_lv = Y_me[bus_hv, bus_lv]
        pp_val_hv_lv = Y_pp_right_order[bus_hv, bus_lv]
        trafo_bug["tap_step_percent"]
        np.sqrt(np.real(me_val_hv_lv / pp_val_hv_lv))
        pdb.set_trace()

    # error of side of tap
    if False:
        # don't know which case...
        bus_lv = 25
        bus_hv = 26
        trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_lv) & (pp_net.trafo["hv_bus"] == bus_hv)]
        me_val = Y_me[bus_lv, bus_hv]
        pp_val = Y_pp_right_order[bus_lv, bus_hv]
        pdb.set_trace()

        # error for non diagonal value of phase shifter
    if False:
        bus_lv = 151
        bus_hv = 155
        trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_lv) & (pp_net.trafo["hv_bus"] == bus_hv)]
        me_val = Y_me[bus_lv, bus_hv]
        pp_val = Y_pp_right_order[bus_lv, bus_hv]

    # error for non diagonal value when tap_pos is +1 (usually it was -1)
    if False:
        # case_name = "test_case6495rte.json"
        bus_lv = 1259
        bus_hv = 5941
        trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_lv) & (pp_net.trafo["hv_bus"] == bus_hv)]
        me_val_lv_hv = Y_me[bus_lv, bus_hv]
        pp_val_lv_hv = Y_pp_right_order[bus_lv, bus_hv]
        me_val_hv_lv = Y_me[bus_hv, bus_lv]
        pp_val_hv_lv = Y_pp_right_order[bus_hv, bus_lv]

    # same bug other example tap_pos is +1
    if False:
        # case_name = "test_case6495rte.json"
        bus_lv = 6187
        bus_hv = 5991
        trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_lv) & (pp_net.trafo["hv_bus"] == bus_hv)]
        trafo_id = np.where((pp_net.trafo["lv_bus"].values == bus_lv) & (pp_net.trafo["hv_bus"].values == bus_hv))
        me_val_lv_hv = Y_me[bus_lv, bus_hv]
        pp_val_lv_hv = Y_pp_right_order[bus_lv, bus_hv]
        me_val_hv_lv = Y_me[bus_hv, bus_lv]
        pp_val_hv_lv = Y_pp_right_order[bus_hv, bus_lv]
        np.sqrt(np.real(me_val_hv_lv / pp_val_hv_lv))  # is exactly trafo_bug["tap_step_percent"]

    # test trafo
    # threshold = 1000.
    # for issue_id in range(where_errors[0].shape[0]):
    #     bus_id_or_error = where_errors[0][issue_id]
    #     trafo_bug = pp_net.trafo.loc[(pp_net.trafo["hv_bus"] == bus_id_or_error)]
    #     if (trafo_bug.shape[0] == 1) and (error_q[bus_id_or_error, bus_id_or_error] > threshold):
    #         print("hv bus")
    #         pdb.set_trace()
    #     trafo_bug = pp_net.trafo.loc[(pp_net.trafo["lv_bus"] == bus_id_or_error)]
    #     if (trafo_bug.shape[0] == 1) and (error_q[bus_id_or_error, bus_id_or_error] > threshold):
    #         print("lv bus")
    #         pdb.set_trace()

    if False:
        # case_name = "test_case6495rte.json"
        bus_hv = 4164
        shunt_bug = pp_net.shunt.loc[pp_net.shunt["bus"] == bus_hv]
        trafo_bug = pp_net.trafo.loc[pp_net.trafo["hv_bus"] == bus_hv]
        id_trafo = np.where((pp_net.trafo["hv_bus"].values == bus_hv))
        trafo_bug
        tap_step = trafo_bug["tap_step_percent"].values[0]
        ratio = 1 + (trafo_bug["tap_pos"] - trafo_bug["tap_neutral"]) * trafo_bug["tap_step_percent"] / 100.
        bus_lv = int(trafo_bug["lv_bus"])
        me_val_lv_hv = Y_me[bus_lv, bus_hv]
        pp_val_lv_hv = Y_pp_right_order[bus_lv, bus_hv]
        me_val_hv_lv = Y_me[bus_hv, bus_lv]
        pp_val_hv_lv = Y_pp_right_order[bus_hv, bus_lv]
        me_val_hv_hv = Y_me[bus_hv, bus_hv]
        me_val_lv_lv = Y_me[bus_lv, bus_lv]
        pp_val_hv_hv = Y_pp_right_order[bus_hv, bus_hv]
        pp_val_lv_lv = Y_pp_right_order[bus_lv, bus_lv]
        np.abs(me_val_lv_lv / pp_val_lv_lv)
        np.abs(me_val_hv_hv / pp_val_hv_hv)
        # see pandapower.build_branch._calc_tap_from_dataframe

    if False:
        # case_name = "test_case6495rte.json"
        bus_lv = 6195
        trafo_bug = pp_net.trafo.loc[pp_net.trafo["lv_bus"] == bus_lv]
        shunt_bug = pp_net.shunt.loc[pp_net.shunt["bus"] == bus_lv]
        bus_hv = int(trafo_bug["hv_bus"])
        me_val_hv_hv = Y_me[bus_hv, bus_hv]
        me_val_lv_lv = Y_me[bus_lv, bus_lv]
        pp_val_hv_hv = Y_pp_right_order[bus_hv, bus_hv]
        pp_val_lv_lv = Y_pp_right_order[bus_lv, bus_lv]
        np.abs(me_val_lv_lv / pp_val_lv_lv)
        np.abs(me_val_hv_hv / pp_val_hv_hv)

        # TODO add the first two lines (the function _calc_tap_from_dataframe) to the converter !
        from pandapower.build_branch import _calc_tap_from_dataframe
        vn_trafo_hv, vn_trafo_lv, shift = _calc_tap_from_dataframe(pp_net, pp_net.trafo)
        tap_angles = np.deg2rad(trafo_bug["tap_step_degree"].values)
        if ~np.isfinite(tap_angles):
            tap_angles = 0.
        vn = trafo_bug["vn_lv_kv"].values
        u1 = vn
        tap_steps = (trafo_bug["tap_pos"].values - trafo_bug["tap_neutral"].values) * trafo_bug["tap_step_percent"].values / 100.
        du = u1 * tap_steps
        np.sqrt((u1 + du * np.cos(tap_angles)) ** 2 + (du * np.sin(tap_angles)) ** 2)
        pdb.set_trace()

    if False:
        if case_name == 'test_case2848rte.json':
            bus_lv = 97
            trafo_bug = pp_net.trafo.loc[pp_net.trafo["lv_bus"] == bus_lv]
            shunt_bug = pp_net.shunt.loc[pp_net.shunt["bus"] == bus_lv]
            bus_hv = int(trafo_bug["hv_bus"].values)
            trafo_id = np.where(pp_net.trafo["lv_bus"].values == bus_lv)[0]
            ratio = 1 + (trafo_bug["tap_pos"] - trafo_bug["tap_neutral"]) * trafo_bug["tap_step_percent"] / 100.

            me_val_hv_hv = Y_me[bus_hv, bus_hv]
            me_val_lv_lv = Y_me[bus_lv, bus_lv]
            pp_val_hv_hv = Y_pp_right_order[bus_hv, bus_hv]
            pp_val_lv_lv = Y_pp_right_order[bus_lv, bus_lv]

            me_val_hv_lv = Y_me[bus_hv, bus_lv]
            me_val_lv_hv = Y_me[bus_lv, bus_hv]
            pp_val_hv_lv = Y_pp_right_order[bus_hv, bus_lv]
            pp_val_lv_hv = Y_pp_right_order[bus_lv, bus_hv]

            np.abs(me_val_lv_lv / pp_val_lv_lv)
            np.abs(me_val_hv_hv / pp_val_hv_hv)

            # TODO add the first two lines (the function _calc_tap_from_dataframe) to the converter !
            from pandapower.build_branch import _calc_tap_from_dataframe
            vn_trafo_hv, vn_trafo_lv, shift = _calc_tap_from_dataframe(pp_net, pp_net.trafo)
            tap_angles = np.deg2rad(trafo_bug["tap_step_degree"].values)
            if ~np.isfinite(tap_angles):
                tap_angles = 0.
            vn = trafo_bug["vn_lv_kv"].values
            u1 = vn
            tap_steps = (trafo_bug["tap_pos"].values - trafo_bug["tap_neutral"].values) * trafo_bug["tap_step_percent"].values / 100.
            du = u1 * tap_steps
            v_modif = np.sqrt((u1 + du * np.cos(tap_angles)) ** 2 + (du * np.sin(tap_angles)) ** 2)

            from pandapower.build_branch import _calc_branch_values_from_trafo_df, get_trafo_values
            from pandapower.build_branch import _calc_nominal_ratio_from_dataframe, _calc_r_x_y_from_dataframe
            from pandapower.build_branch import _calc_tap_from_dataframe, BASE_KV, _calc_r_x_from_dataframe

            ppc = copy.deepcopy(pp_net._ppc)
            bus_lookup = pp_net["_pd2ppc_lookups"]["bus"]
            trafo_df = pp_net["trafo"]
            lv_bus = get_trafo_values(trafo_df, "lv_bus")
            vn_lv = ppc["bus"][bus_lookup[lv_bus], BASE_KV]
            vn_trafo_hv, vn_trafo_lv, shift_pp = _calc_tap_from_dataframe(pp_net, trafo_df)
            ratio_pp = _calc_nominal_ratio_from_dataframe(ppc, trafo_df, vn_trafo_hv, vn_trafo_lv, bus_lookup)
            r_pp, x_pp, b_pp = _calc_r_x_y_from_dataframe(pp_net, trafo_df, vn_trafo_lv, vn_lv, pp_net.sn_mva)

            y_trafo = 1. / (r_pp[trafo_id] + 1j * x_pp[trafo_id])
            ratio_pp[trafo_id]
            y_trafo

            v_modif / vn
            vn / v_modif
            # pdb.set_trace()

    if case_name == "test_case2848rte.json":
        np.where((error_q == np.max(error_q)).todense())
        bus_bug = 685
        shunt_bug = pp_net.shunt.loc[pp_net.shunt["bus"] == bus_bug]
        trafo_hv_bug = pp_net.trafo.loc[pp_net.trafo["hv_bus"] == bus_bug]
        trafo_lv_bug = pp_net.trafo.loc[pp_net.trafo["lv_bus"] == bus_bug]
        line_from_bug = pp_net.line.loc[pp_net.line["from_bus"] == bus_bug]
        line_to_bug = pp_net.line.loc[pp_net.line["to_bus"] == bus_bug]
        print("\t\tA trafo is connected to a bus where there is a mismatch in the Ybus: "
              f"{int(np.where(pp_net.trafo['lv_bus'].values == bus_bug)[0])}")

print("IV - Check for the initialization (dc powerflow)")
print("1) check that the results are same for dc lightsim and dc pandapower")
Vinit = np.ones(backend.nb_bus_total, dtype=np.complex_) * pp_net["_options"]["init_vm_pu"]
backend._grid.deactivate_result_computation()
Vdc = backend._grid.dc_pf(Vinit, max_iter, tol_this)
backend._grid.reactivate_result_computation()
test_ok = True

Ydc_me = backend._grid.get_Ybus()
Sdc_me = backend._grid.get_Sbus()
if np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) >= 100.*tol:
    test_ok = False
    print(f"\t Error for the DC approximation: resulting voltages are different "
          f"{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f}pu")
    highest_bug = np.argmax(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub]))

elif np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) >= tol:
    test_ok = False
    print("\t Warning: maximum difference after DC approximation is "
          f"{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f} which is higher than "
          f"the tolerance (this is just a warning because we noticed this could happen even if the results "
          f"match perfectly. Probably some conversion issue with complex number and radian / degree.")
if test_ok:
    print("\t Info: Same result for the dc approximation (same V).")
else:
    pass

print("2) check that the Sbus vector is same for PP and lightisim in DC")
# get the Sbus for pandapower
from pandapower.pd2ppc import _pd2ppc
from pandapower.pf.run_newton_raphson_pf import  _get_pf_variables_from_ppci
from pandapower.pypower.idx_brch import F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS
from pandapower.pypower.idx_bus import VA, GS
from pandapower.pypower.makeBdc import makeBdc
from pandapower.pypower.makeSbus import makeSbus

pp_net._pd2ppc_lookups = {"bus": np.array([], dtype=int), "ext_grid": np.array([], dtype=int),
                          "gen": np.array([], dtype=int), "branch": np.array([], dtype=int)}
# convert pandapower net to ppc
ppc, ppci = _pd2ppc(pp_net)
baseMVA, bus, gen, branch, ref, pv, pq, on, gbus, _, refgen = _get_pf_variables_from_ppci(ppci)
Va0 = bus[:, VA] * (np.pi / 180.)
B, Bf, Pbusinj, Pfinj = makeBdc(bus, branch)
from pandapower.pypower.idx_brch import F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS

Pbus = makeSbus(baseMVA, bus, gen) - Pbusinj - bus[:, GS] / baseMVA
Pbus_pp_ro = Pbus[pp_vect_converter]
error_p = np.abs(np.real(Sdc_me) - np.real(Pbus_pp_ro))
test_ok = True

#### pandapower DC algo (yet another one)
Va = copy.deepcopy(Va0)
pvpq = np.r_[pv, pq]
pvpq_matrix = B[pvpq.T, :].tocsc()[:, pvpq]
ref_matrix = np.transpose(Pbus[pvpq] - B[pvpq.T, :].tocsc()[:, ref] * Va0[ref])
Va[pvpq] = np.real(scipy.sparse.linalg.spsolve(pvpq_matrix, ref_matrix))
####

if np.max(error_p) > tol:
    test_ok = False
    print(f"\t Error: P do not match for Sbus (dc), maximum absolute error is {np.max(error_p):.5f} MW")
    print(f"\t Error: significative difference for bus index (lightsim): {np.where(error_p > tol)[0]}")

error_q = np.abs(np.imag(Sdc_me) - np.imag(Pbus_pp_ro))
if np.max(error_q) > tol:
    test_ok = False
    print(f"\t Error: Q do not match for Sbus (dc), maximum absolute error is {np.max(error_q):.5f} MVAr")
    print(f"\t Error: significative difference for bus index (lightsim): {np.where(error_q > tol)[0]}")
if test_ok:
    print("\t Info: Sbus (dc) is the same")

print("3) check that the Ybus matrix is same for PP and lightisim in DC")
pandapower.rundcpp(pp_net)
Ydc_pp = backend.init_pp_backend._grid._ppc["internal"]["Bbus"]
Ydc_pp_right_order = Ydc_pp[pp_vect_converter.reshape(nb_sub, 1), pp_vect_converter.reshape(1, nb_sub)]
test_ok = True
error_p = np.abs(np.real(Ydc_me) - np.real(Ydc_pp_right_order))
if np.max(error_p) > tol:
    test_ok = False
    print(f"\t Error: P do not match for Ybus (dc mode), maximum absolute error is {np.max(error_p):.5f}")
    diffs_p_dc = error_p > tol
    diffs_p_dc = diffs_p_dc.todense()
    ind_row, ind_cols = np.where(diffs_p_dc)
    rest = f"... {len(ind_row)} in total" if len(ind_row) >= 10 else ""
    print(f"\t Error: significative difference for bus index (lightsim): {ind_row[:10], ind_cols[:10]} "
          f"{rest}")
    where_max = np.where(error_p.todense() >= 100)

error_q = np.abs(np.imag(Ydc_me) - np.imag(Ydc_pp_right_order))
if np.max(error_q) > tol:
    test_ok = False
    print(f"\t Error: Q do not match for Ybus (dc mdoe), maximum absolute error is {np.max(error_q):.5f}")
    diffs_q_dc = error_q > tol
    diffs_q_dc = diffs.todense()
    ind_row, ind_cols = np.where(diffs_q_dc)
    rest = f"... {len(ind_row)} in total" if len(ind_row) >= 10 else ""
    print(f"\t Error: significative difference for bus index (lightsim): {ind_row[:10], ind_cols[:10]} "
          f"{rest}")
    where_max = np.where(error_p.todense() >= 100)

if test_ok:
    print("\t Info: Ybus (dc mode) is correct")
else:
    if case_name == "test_case2848rte":
        where_errors = np.where(diffs_p_dc)
        nb_error = where_errors[0].shape[0]
        error_diag_coeff = where_errors[0] == where_errors[1]
        error_p_d = error_p.todense()
        indx_max = np.where(error_p_d == np.max(error_p_d))

        bus_lv = 24
        trafo_bug = pp_net.trafo.loc[pp_net.trafo["lv_bus"] == bus_lv]
        bus_hv = int(trafo_bug["hv_bus"].values)
        trafo_id = np.where(pp_net.trafo["lv_bus"].values == bus_lv)[0]
        ratio = 1 + (trafo_bug["tap_pos"] - trafo_bug["tap_neutral"]) * trafo_bug["tap_step_percent"] / 100.

        me_val_hv_hv = Ydc_me[bus_hv, bus_hv]
        me_val_lv_lv = Ydc_me[bus_lv, bus_lv]
        pp_val_hv_hv = Ydc_pp_right_order[bus_hv, bus_hv]
        pp_val_lv_lv = Ydc_pp_right_order[bus_lv, bus_lv]

        me_val_hv_lv = Ydc_me[bus_hv, bus_lv]
        me_val_lv_hv = Ydc_me[bus_lv, bus_hv]
        pp_val_hv_lv = Ydc_pp_right_order[bus_hv, bus_lv]
        pp_val_lv_hv = Ydc_pp_right_order[bus_lv, bus_hv]

        from pandapower.build_branch import _calc_branch_values_from_trafo_df, get_trafo_values
        from pandapower.build_branch import _calc_nominal_ratio_from_dataframe, _calc_r_x_y_from_dataframe
        from pandapower.build_branch import _calc_tap_from_dataframe, BASE_KV, _calc_r_x_from_dataframe

        ppc = copy.deepcopy(pp_net._ppc)
        bus_lookup = pp_net["_pd2ppc_lookups"]["bus"]
        trafo_df = pp_net["trafo"]
        lv_bus = get_trafo_values(trafo_df, "lv_bus")
        vn_lv = ppc["bus"][bus_lookup[lv_bus], BASE_KV]
        vn_trafo_hv, vn_trafo_lv, shift_pp = _calc_tap_from_dataframe(pp_net, trafo_df)
        ratio_pp = _calc_nominal_ratio_from_dataframe(ppc, trafo_df, vn_trafo_hv, vn_trafo_lv, bus_lookup)
        r_pp, x_pp, b_pp = _calc_r_x_y_from_dataframe(pp_net, trafo_df, vn_trafo_lv, vn_lv, pp_net.sn_mva)
        pdb.set_trace()

# np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub]))  # 0.061585492793420286

print("3) check that lightsim ac pf init with pp dc pf give same results (than pp)")
Vinit = np.ones(backend.nb_bus_total, dtype=np.complex_) * pp_net["_options"]["init_vm_pu"]
Vinit[:nb_sub] = V_init_ref[pp_vect_converter]
conv = backend._grid.ac_pf(Vinit, max_iter, tol_this)

if conv.shape[0] == 0:
    print("\t Error: the lightsim diverge when initialized with pandapower Vinit_dc")
    test_ok = False
else:
    lpor, lqor, lvor, laor = backend._grid.get_lineor_res()
    tpor, tqor, tvor, taor = backend._grid.get_trafohv_res()
    tpex, tqex, tvex, taex = backend._grid.get_trafolv_res()
    nb_trafo = tpor.shape[0]
    nb_powerline = lpor.shape[0]
    p_or_me2 = np.concatenate((lpor, tpor))
    q_or_me2 = np.concatenate((lqor, tqor))
    v_or_me2 = np.concatenate((lvor, tvor))
    a_or_me2 = 1000. * np.concatenate((laor, taor))
    test_ok = True
    # pdb.set_trace()
    max_mis = np.max(np.abs(p_or_me2 - por_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: por do not match, maximum absolute error is {max_mis:.5f} MW")
    max_mis = np.max(np.abs(q_or_me2 - qor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr")
    max_mis = np.max(np.abs(v_or_me2 - vor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: vor do not match, maximum absolute error is {max_mis:.5f} kV")
    max_mis = np.max(np.abs(a_or_me2 - aor_pp))
    if max_mis > tol:
        test_ok = False
        print(f"\t Error: aor do not match, maximum absolute error is {max_mis:.5f} A")
    if test_ok:
        print(f"\t flows ok for {case_name} after initialized with DC approximation from PP")

    #for test_case2848rte
    if case_name == 'test_case2848rte.json':
        max_id_bug = np.argmax(np.abs(p_or_me2 - por_pp))
        id_trafo = max_id_bug - nb_powerline
        if id_trafo >= 0:
            # id_trafo = 322
            np.where(backend.init_pp_backend._grid.trafo.columns == "shift_degree")
            backend.init_pp_backend._grid.trafo.values[id_trafo, 9]
            angle_error = np.angle((p_or_me2 - por_pp)[max_id_bug] + (q_or_me2 - qor_pp)[max_id_bug] * 1j, deg=True)
            angle_me = np.angle(p_or_me2[max_id_bug] + q_or_me2[max_id_bug] * 1j, deg=True)
            angle_pp = np.angle(por_pp[max_id_bug] + qor_pp[max_id_bug] * 1j, deg=True)
            pdb.set_trace()

if test_ok:
    print("\t Info: All good for DC approximation")
else:
    # for case 'test_case2848rte.json'
    if case_name == 'test_case2848rte.json':
        max_mis_p = np.max(np.abs(p_or_me2 - por_pp))
        id_trafo_bug = np.where(np.abs(p_or_me2 - por_pp) >= max_mis_p -1.0)[0] - nb_powerline
        if np.all(id_trafo_bug >= 0):
            pp_hv, pp_lv = pp_net.res_trafo.iloc[id_trafo_bug][["p_hv_mw", "p_lv_mw"]].iloc[0]
            me_lv = tpex[id_trafo_bug]
            me_hv = tpor[id_trafo_bug]
            pdb.set_trace()

print("V - Check trafo proper conversion to r,x, b")
test_ok = True
from lightsim2grid_cpp import GridModel, PandaPowerConverter, SolverType
from pandapower.build_branch import _calc_branch_values_from_trafo_df, get_trafo_values
from pandapower.build_branch import _calc_nominal_ratio_from_dataframe, _calc_r_x_y_from_dataframe
from pandapower.build_branch import _calc_tap_from_dataframe, BASE_KV, _calc_r_x_from_dataframe

# my values
# trafo_df = pp_net["trafo"]
# vn_trafo_hv, vn_trafo_lv, shift = _calc_tap_from_dataframe(pp_net, trafo_df)

converter = PandaPowerConverter()
converter.set_sn_mva(pp_net.sn_mva)  # TODO raise an error if not set !
converter.set_f_hz(pp_net.f_hz)

# fix the missing values
tap_neutral = 1.0 * pp_net.trafo["tap_neutral"].values
if np.any(~np.isfinite(tap_neutral)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"tap_neutral\"], they have been replaced by 0")
tap_neutral[~np.isfinite(tap_neutral)] = 0.

if np.any(tap_neutral != 0.):
    raise RuntimeError("lightsim converter supposes that tap_neutral is 0 for the transformers")

tap_step_pct = 1.0 * pp_net.trafo["tap_step_percent"].values
if np.any(~np.isfinite(tap_step_pct)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_percent\"], they have been replaced by 0")
tap_step_pct[~np.isfinite(tap_step_pct)] = 0.

tap_pos = 1.0 * pp_net.trafo["tap_pos"].values
if np.any(~np.isfinite(tap_pos)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"tap_pos\"], they have been replaced by 0")
tap_pos[~np.isfinite(tap_pos)] = 0.

shift_ = 1.0 * pp_net.trafo["shift_degree"].values
if np.any(~np.isfinite(tap_pos)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"shift_degree\"], they have been replaced by 0")
shift_[~np.isfinite(shift_)] = 0.

is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
if np.any(~np.isfinite(is_tap_hv_side)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"tap_side\"], they have been replaced by \"hv\"")
is_tap_hv_side[~np.isfinite(is_tap_hv_side)] = True

if np.any(pp_net.trafo["tap_phase_shifter"].values):
    raise RuntimeError("ideal phase shifter are not modeled. Please remove all trafo with "
                       "pp_net.trafo[\"tap_phase_shifter\"] set to True.")

tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
if np.any(~np.isfinite(tap_pos)):
    warnings.warn("There were some Nan in the pp_net.trafo[\"tap_step_degree\"], they have been replaced by 0")
tap_angles_[~np.isfinite(tap_angles_)] = 0.
tap_angles_ = np.deg2rad(tap_angles_)
trafo_r, trafo_x, trafo_b = \
    converter.get_trafo_param(tap_step_pct,
                              tap_pos,
                              tap_angles_,  # in radian !
                              is_tap_hv_side,
                              pp_net.bus.loc[pp_net.trafo["hv_bus"]]["vn_kv"],
                              pp_net.bus.loc[pp_net.trafo["lv_bus"]]["vn_kv"],
                              pp_net.trafo["vk_percent"].values,
                              pp_net.trafo["vkr_percent"].values,
                              pp_net.trafo["sn_mva"].values,
                              pp_net.trafo["pfe_kw"].values,
                              pp_net.trafo["i0_percent"].values,
                              )

# todo remove
tap_steps = tap_step_pct * (tap_pos - tap_neutral) / 100.

u_hv = 1.0 * pp_net.trafo["vn_hv_kv"].values
u_lv = 1.0 * pp_net.trafo["vn_lv_kv"].values
du_hv = pp_net.trafo["vn_hv_kv"].values * tap_steps
du_lv = pp_net.trafo["vn_lv_kv"].values * tap_steps
vn_trafo_hv_me = 1.0 * pp_net.trafo["vn_hv_kv"].values
vn_trafo_hv_me[is_tap_hv_side] = np.sqrt((u_hv[is_tap_hv_side] + du_hv[is_tap_hv_side] * np.cos(tap_angles_[is_tap_hv_side])) ** 2 +
                                         (du_hv[is_tap_hv_side] * np.sin(tap_angles_[is_tap_hv_side])) ** 2)

vn_trafo_lv_me = 1.0 * pp_net.trafo["vn_lv_kv"].values
vn_trafo_lv_me[~is_tap_hv_side] = np.sqrt((u_lv[~is_tap_hv_side] + du_lv[~is_tap_hv_side] * np.cos(tap_angles_[~is_tap_hv_side])) ** 2 +
                                          (du_lv[~is_tap_hv_side] * np.sin(tap_angles_[~is_tap_hv_side])) ** 2)
# todo end remove

# pp values
ppc = copy.deepcopy(pp_net._ppc)
bus_lookup = pp_net["_pd2ppc_lookups"]["bus"]
trafo_df = pp_net["trafo"]
lv_bus = get_trafo_values(trafo_df, "lv_bus")
vn_lv = ppc["bus"][bus_lookup[lv_bus], BASE_KV]
vn_trafo_hv, vn_trafo_lv, shift_pp = _calc_tap_from_dataframe(pp_net, trafo_df)
ratio = _calc_nominal_ratio_from_dataframe(ppc, trafo_df, vn_trafo_hv, vn_trafo_lv, bus_lookup)
r_t, x_t, b_t = _calc_r_x_y_from_dataframe(pp_net, trafo_df, vn_trafo_lv, vn_lv, pp_net.sn_mva)

# todo remove
matches_hv = vn_trafo_hv_me == vn_trafo_hv
matches_lv = vn_trafo_lv_me == vn_trafo_lv
import pdb
pdb.set_trace()
# todo end remove

# check where there are mismatch if any
val_r_pp = r_t
val_r_me = trafo_r
all_equals_r = np.abs(val_r_pp - val_r_me) <= tol
if not np.all(all_equals_r):
    test_ok = False
    print(f"\t Error: some trafo resistance are not equal, max error: {np.max(np.abs(val_r_pp - val_r_me)):.5f}")

val_x_pp = x_t
val_x_me = trafo_x
all_equals_x = np.abs(val_x_pp - val_x_me) <= tol
if not np.all(all_equals_x):
    test_ok = False
    print(f"\t Error: some trafo x are not equal, max error: {np.max(np.abs(val_x_pp - val_x_me)):.5f}")

val_ib_pp = np.imag(b_t)
val_ib_me = np.imag(trafo_b)
all_equals_imag_b = np.abs(val_ib_pp - val_ib_me) <= tol
if not np.all(all_equals_imag_b):
    test_ok = False
    print(f"\t Error: some trafo (imag) b are not equal, max error: {np.max(np.abs(val_ib_pp - val_ib_me)):.5f}")

val_reb_pp = np.real(b_t)
val_reb_me = np.real(trafo_b)
all_equals_real_b = np.abs(val_reb_pp - val_reb_me) <= tol
if not np.all(all_equals_real_b):
    test_ok = False
    print(f"\t Error: some trafo (real) b are not equal, max error: {np.max(np.abs(val_reb_pp - val_reb_me)):.5f}")

if test_ok:
    print("\t Info: all good for the transformers parameters")
vals_problem = np.real(b_t)[~all_equals_real_b]
bus_vals_problem = lv_bus[~all_equals_real_b]

if False:
    trafo_id_ls = 70
    trafo_bug_ybus = pp_net.trafo.iloc[trafo_id_ls]
    id_trafo_pp = int(trafo_bug_ybus.name)
    val_reb_pp[trafo_id_ls]
    val_reb_me[trafo_id_ls]
print("All checks are over. I hope the debugging information helps (-: ")

#########################"
#garbage

# to clear a pandapower grid
# from pandapower.pf.run_newton_raphson_pf import _run_dc_pf, _get_pf_variables_from_ppci
# from pandapower.pd2ppc import _pd2ppc, _calc_pq_elements_and_add_on_ppc, _ppc2ppci
# from pandapower.pypower.idx_bus import VM, VA
# pp_net.res_bus.drop(axis=0, index=np.arange(pp_net.res_bus.shape[0]), inplace=True)
# pp_net.res_line.drop(axis=0, index=np.arange(pp_net.res_line.shape[0]), inplace=True)
# pp_net.res_trafo.drop(axis=0, index=np.arange(pp_net.res_trafo.shape[0]), inplace=True)
# pp_net.res_ext_grid.drop(axis=0, index=np.arange(pp_net.res_ext_grid.shape[0]), inplace=True)
# pp_net.res_load.drop(axis=0, index=np.arange(pp_net.res_load.shape[0]), inplace=True)
# pp_net.res_shunt.drop(axis=0, index=np.arange(pp_net.res_shunt.shape[0]), inplace=True)
# pp_net.res_shunt.drop(axis=0, index=np.arange(pp_net.res_shunt.shape[0]), inplace=True)
# pp_net._ppc = None
# pandapower.rundcpp(pp_net)
# ppc, ppci = _pd2ppc(pp_net)
# ppci = _run_dc_pf(ppci)
# baseMVA, bus, gen, branch, ref, pv, pq, _, _, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
# V0_pp = copy.deepcopy(V0[ppci._pd2ppc_lookups["bus"]])
# V0_pp[np.isnan(V0_pp)] = pp_net["_options"]["init_vm_pu"]