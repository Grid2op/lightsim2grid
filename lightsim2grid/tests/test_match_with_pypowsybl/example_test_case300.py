import json
import os
import time

import numpy as np

import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf

from lightsim2grid.gridmodel import init_from_pypowsybl
from utils_for_slack import get_same_slack, get_pypowsybl_parameters


def debug_with_olf_physical_params(
    el,
    olf_par,
    tol,
    is_line,
    row_ids=None,
    col_ids=None,
    val_ids=None,
    mat_id=None):
    olf_r = olf_par["r"]
    olf_x = olf_par["x"]
    if "b1" in olf_par:
        olf_b1 = olf_par["b1"]
    else:
        olf_b1 = 0.
    if "b2"in olf_par:
        olf_b2 = olf_par["b2"]
    else:
        olf_b2 = 0.
    if "g1" in olf_par:
        olf_g1 = olf_par["g1"]
    else:
        olf_g1 = 0.
    if "g2"in olf_par:
        olf_g2 = olf_par["g2"]
    else:
        olf_g2 = 0.
    olf_bus_f = int(olf_par["num1"])
    olf_bus_t = int(olf_par["num2"])
    
    if is_line:
        lf_bus_f = int(el.bus_or_id)
        lf_bus_t = int(el.bus_or_id)
    else:
        lf_bus_f = int(el.bus_hv_id)
        lf_bus_t = int(el.bus_lv_id)
        
    lf_r = el.r_pu
    lf_x = el.x_pu
    if is_line:
        lf_b1 = np.imag(el.h_or_pu)
        lf_b2 = np.imag(el.h_ex_pu) 
        lf_g1 = np.real(el.h_or_pu)
        lf_g2 = np.real(el.h_ex_pu)
    else:
        lf_b1 = np.imag(el.h_pu) * 0.5
        lf_b2 = np.imag(el.h_pu) * 0.5
        lf_g1 = np.real(el.h_pu) * 0.5
        lf_g2 = np.real(el.h_pu) * 0.5
        if "a1" in olf_par:
            import pdb
            pdb.set_trace()
    
    ok_ = True
    tmp_ = {}
    if abs(olf_r - lf_r) >= tol:
        tmp_["r_ls"] = float(lf_r)
        tmp_["r_olf"] = float(olf_r)
        ok_ = False
        print(f"{el.name} r: {olf_r} vs {lf_r}")
    if abs(olf_x - lf_x) >= tol:
        tmp_["x_ls"] = float(lf_x)
        tmp_["x_olf"] = float(olf_x)
        ok_ = False
        print(f"{el.name} x: {olf_x} vs {lf_x}")
    if abs(olf_b1 - lf_b1) >= tol:
        tmp_["b1_ls"] = float(lf_b1)
        tmp_["b1_olf"] = float(olf_b1)
        ok_ = False
        print(f"{el.name} b1: {olf_b1} vs {lf_b1}")
    if abs(olf_b2 - lf_b2) >= tol:
        tmp_["b2_ls"] = float(lf_b2)
        tmp_["b2_olf"] = float(olf_b2)
        ok_ = False
        print(f"{el.name} b2: {olf_b2} vs {lf_b2}")
    if abs(olf_g1 - lf_g1) >= tol:
        tmp_["g1_ls"] = float(lf_g1)
        tmp_["g1_olf"] = float(olf_g1)
        ok_ = False
        print(f"{el.name} g1: {olf_g1} vs {lf_g1}")
    if abs(olf_g2 - lf_g2) >= tol:
        tmp_["g2_ls"] = float(lf_g2)
        tmp_["g2_olf"] = float(olf_g2)
        ok_ = False
        print(f"{el.name} g2: {olf_g2} vs {lf_g2}")   

         
    # fill sp mat
    if row_ids is not None:
        ys = 1. / (olf_r + 1j * olf_x)
        h_or = (olf_b1 + 1j * olf_g1)
        h_ex = (olf_b2 + 1j * olf_g2)
        # yac_ff_(i) = (ys + h_or)
        # yac_tt_(i) = (ys + h_ex)
        # yac_tf_(i) = -ys
        # yac_ft_(i) = -ys
        tau = olf_par.get('r1', 1.)
        # real_type cos_theta = std::cos(theta_shift);
        # real_type sin_theta = std::sin(theta_shift);
        # eitheta_shift = {cos_theta, sin_theta};
        # emitheta_shift = {cos_theta, -sin_theta};
        eitheta_shift = 1.  # TODO phase shifter
        emitheta_shift = 1.
        
        id_ff = mat_id
        id_tt = mat_id + 1
        id_tf = mat_id + 2
        id_ft = mat_id + 3
        
        val_ids[id_ff] = (ys + h_or) / (tau * tau)
        row_ids[id_ff] = olf_bus_f
        col_ids[id_ff] = olf_bus_f
        
        val_ids[id_tt] = (ys + h_ex)
        row_ids[id_tt] = olf_bus_t
        col_ids[id_tt] = olf_bus_t
        
        val_ids[id_tf] = -ys / tau * emitheta_shift
        row_ids[id_tf] = olf_bus_t
        col_ids[id_tf] = olf_bus_f
        
        val_ids[id_ft] = -ys / tau * eitheta_shift
        row_ids[id_ft] = olf_bus_f
        col_ids[id_ft] = olf_bus_t
    return ok_, tmp_


if __name__ == "__main__":
    debug_dir = "/tmp"

    case_name = "ieee300"
    
    slack_pypowysbl, slack_ls = get_same_slack(case_name)
    
    pypow_grid = getattr(pypow.network, f"create_{case_name}")()
    res = pypow_lf.run_ac(pypow_grid)
    ls_grid = init_from_pypowsybl(
        pypow_grid,
        slack_bus_id=slack_ls,
        sort_index=False,
        buses_for_sub=True,
        n_busbar_per_sub=1
        )
    
    pypowsybl_parameters = get_pypowsybl_parameters(slack_pypowysbl)
    if debug_dir is not None:
        if not os.path.exists(debug_dir):
            os.mkdir(debug_dir)
        pypowsybl_parameters.provider_parameters["debugDir"] = str(debug_dir)  # os.path.abspath(os.path.join(".", "olf_debug"))
        pypowsybl_parameters.provider_parameters["reportedFeatures"] = "NEWTON_RAPHSON_LOAD_FLOW"
    
    pypow_lf.run_ac(pypow_grid, parameters=pypowsybl_parameters)
    
    res_ls = ls_grid.ac_pf(
        np.ones(ls_grid.get_bus_vn_kv().shape[0], dtype=complex),
        10,
        1e-6)
    
    if debug_dir is not None:
        fn_li = [el for el in sorted(os.listdir(debug_dir)) if el.endswith("json")]
        if len(fn_li) == 0:
            raise RuntimeError(f"OLF powerflow has not been run with debug dir {debug_dir}")
        fn = fn_li[-1]
        with open(os.path.join(debug_dir, fn), "r", encoding="utf-8") as f:
            physical_params_olf = json.load(f)
        
        # check branches
        olf_branch_params = physical_params_olf["branches"]
        olf_branch_dict = {el["id"]: el for el in olf_branch_params}
        tol = 1e-9
        row_ids = np.zeros(4*len(olf_branch_params), dtype=int)
        col_ids = np.zeros(4*len(olf_branch_params), dtype=int)
        val_ids = np.zeros(4*len(olf_branch_params), dtype=complex)
        mat_id = 0
        # check lines
        for el in ls_grid.get_lines():
            if not el.connected:
                continue
            if el.name not in olf_branch_dict:
                print(f"ERROR {el.name} (line) not found in olf branch dictionary")
                continue
            olf_par = olf_branch_dict[el.name]
            debug_with_olf_physical_params(
                el, olf_par, tol, True,
                row_ids, col_ids, val_ids,
                mat_id)
            mat_id += 4
        # check trafo
        for el in ls_grid.get_trafos():
            if not el.connected:
                continue
            if el.name not in olf_branch_dict:
                print(f"ERROR {el.name} (trafo) not found in olf branch dictionary")
                continue
            olf_par = olf_branch_dict[el.name]
            debug_with_olf_physical_params(
                el, olf_par, tol, False,
                row_ids, col_ids, val_ids, mat_id)
            mat_id += 4
        from scipy.sparse import coo_array

        Ybus_csr_olf = coo_array((val_ids, (row_ids, col_ids))).tocsr()
        # TODO missing shunt in OLF this way !!!
        # shunts at bus : [ 0,  3,  7, 25, 34, 62, 75]
        Ybus_csr_ls = ls_grid.get_Ybus()