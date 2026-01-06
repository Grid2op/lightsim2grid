# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import tempfile
import pandapower as pp
import pandapower.networks as pn
import os
import copy
import numpy as np
import scipy
import warnings

from lightsim2grid import LightSimBackend
import lightsim2grid
from lightsim2grid_cpp import PandaPowerConverter
from pandapower.build_branch import _calc_branch_values_from_trafo_df
try:
    from lightsim2grid.solver import KLUSolver
    ClassSolver = KLUSolver
except ImportError as exc_:
    from lightsim2grid.solver import SparseLUSolver
    ClassSolver = SparseLUSolver

from global_var_tests import (
    MAX_PP2_DATAREADER,
    CURRENT_PP_VERSION,
)

    
TIMER_INFO = False  # do i print information regarding computation time


def create_transformers(n_pdp):
    # bus = n_pdp.bus[['vn_kv']]
    # trafo_and_bus = n_pdp.trafo.merge(bus.rename(columns=lambda x: x + '_lv_bus'), left_on='lv_bus', right_index=True, how='left')
    trafo_and_bus = n_pdp.trafo
    n_tap = np.where(~np.isnan(trafo_and_bus['tap_pos']) & ~np.isnan(trafo_and_bus['tap_neutral']) & ~np.isnan(trafo_and_bus['tap_step_percent']),
                     1.0 + (trafo_and_bus['tap_pos'] - trafo_and_bus['tap_neutral']) * trafo_and_bus['tap_step_percent'] / 100.0, 1.0)
    rated_u1 = np.where(trafo_and_bus['tap_side'] == "hv", trafo_and_bus['vn_hv_kv'] * n_tap, trafo_and_bus['vn_hv_kv'])
    rated_u2 = np.where(trafo_and_bus['tap_side'] == "lv", trafo_and_bus['vn_lv_kv'] * n_tap, trafo_and_bus['vn_lv_kv'])
    c = n_pdp.sn_mva / n_pdp.trafo['sn_mva']
    rk = trafo_and_bus['vkr_percent'] / 100 * c
    zk = trafo_and_bus['vk_percent'] / 100 * c
    xk = np.sqrt(zk ** 2 - rk ** 2)
    ym = trafo_and_bus['i0_percent'] / 100
    gm = trafo_and_bus['pfe_kw'] / (trafo_and_bus['sn_mva'] * 1000) / c
    bm = -1 * np.sign(ym) * np.sqrt(ym ** 2 - gm ** 2)
    
    zb_tr = (rated_u2 ** 2) / n_pdp.sn_mva
    r = rk * zb_tr / trafo_and_bus['parallel']
    x = xk * zb_tr / trafo_and_bus['parallel']
    g = gm / zb_tr * trafo_and_bus['parallel']
    b = bm / zb_tr * trafo_and_bus['parallel']
    return r, x, g, b
    
    
class TestACPF(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 3e-5  # results are equal if they match up to tol
        self.raise_error = True

    def test_case14(self):
        case = pn.case14()
        self._aux_test(case)

    def test_case14_with_phaseshift(self):
        case = pn.case14()        
        hv_bus=0
        lv_bus=2
        pp.create_transformer_from_parameters(case,
                                              hv_bus=hv_bus,
                                              lv_bus=lv_bus,
                                            #   sn_mva=1184.0,   # case RTE
                                              sn_mva=9900.0,
                                              vn_hv_kv=case.bus.iloc[hv_bus]["vn_kv"],
                                              vn_lv_kv=case.bus.iloc[lv_bus]["vn_kv"],
                                            #   i0_percent=-0.05152,   # case RTE
                                              i0_percent=0.0,
                                            #   vk_percent=0.404445,  # case RTE
                                              vk_percent=2070.288000,
                                            #   vkr_percent=0.049728, # case RTE
                                              vkr_percent=0.0,
                                              shift_degree=-10.0,  # case RTE
                                            #   shift_degree=0.,
                                              pfe_kw=0.
                                              )
        # self.tol = 2e-3
        case.name = "case14_2"
        self._aux_test(case)
        
    def test_case39(self):
        case = pn.case39()
        self.tol = 3.1e-5  # results are equal if they match up to tol
        self._aux_test(case)
        
    def test_case57(self):
        case = pn.case57()
        self.tol = 3.1e-5  # results are equal if they match up to tol
        self._aux_test(case)

    def test_case118(self):        
        case = pn.case118()
        self.tol = 3.1e-5  # results are equal if they match up to tol
        self._aux_test(case)

    def test_case300(self):
        case = pn.case300()
        self.tol = 3.1e-5  # results are equal if they match up to tol
        self._aux_test(case)

    def test_case1888rte(self):
        case = pn.case1888rte()
        self.tol = 3.1e-5  # results are equal if they match up to tol
        self._aux_test(case)

    # def test_case9241pegase(self):
    #     # too large for CI
    #     case = pn.case9241pegase()
    #     self.tol = 7e-5  # results are equal if they match up to tol
    #     self._aux_test(case)

    # def test_case2848rte(self):
    #     # too large for CI
    #     case = pn.case2848rte()
    #     self.tol = 4e-5  # results are equal if they match up to tol
    #     self._aux_test(case)

    # def test_case6470rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6470rte()
    #     self._aux_test(case)

    # def test_case6495rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6495rte()
    #     self._aux_test(case)

    # def test_case6515rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6515rte()
    #     self._aux_test(case)

    def test_case_illinois200(self):
        case = pn.case_illinois200()
        self._aux_test(case)

    def _assert_or_print(self, boolean: bool, msg: str):
        if not boolean:
            if not self.raise_error:
                print(msg)
            else:
                raise AssertionError(msg)
            
    def _aux_test(self, pn_net):
        n_busbar = 2
        with tempfile.TemporaryDirectory() as path:
            case_name = os.path.join(path, "this_case.json")
            pp.to_json(pn_net, case_name)

            # real_init_file = pp.from_json(case_name)
            LightSimBackend._clear_grid_dependant_class_attributes()
            LightSimBackend.set_env_name(type(self).__name__ + case_name)
            LightSimBackend.set_n_busbar_per_sub(n_busbar)
            if MAX_PP2_DATAREADER < CURRENT_PP_VERSION:
                loader_kwargs = {"pp_orig_file": "pandapower_v3"}
            else:
                loader_kwargs = {"pp_orig_file": "pandapower_v2"}
            backend = LightSimBackend(loader_kwargs=loader_kwargs)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                backend.load_grid(case_name)
            backend.assert_grid_correct()

        nb_sub = backend.n_sub
        pp_net = backend._init_pp_backend._grid
        # first i deactivate all slack bus in pp that are connected but not handled in ls
        pp_net.ext_grid["in_service"] = False
        pp_net.ext_grid.loc[pp_net.ext_grid.index[0], "in_service"] = True
        conv_pp, exc_pp = backend._init_pp_backend.runpf()
        conv, exc_ = backend.runpf()
        
        # backend._debug_Vdc
        # backend._init_pp_backend._grid._ppc["internal"]["Va_it"][0]
        # np.abs(backend._grid.check_solution(backend._init_pp_backend._grid["_ppc"]["internal"]["V"], False)).max()
        self._assert_or_print(conv_pp, "Error: pandapower do not converge, impossible to perform the necessary checks")
        self._assert_or_print(conv, "Error: lightsim do not converge")

        por_pp, qor_pp, vor_pp, aor_pp = copy.deepcopy(backend._init_pp_backend.lines_or_info())
        pex_pp, qex_pp, vex_pp, aex_pp = copy.deepcopy(backend._init_pp_backend.lines_ex_info())

        # I- Check for divergence and equality of flows"
        por_ls, qor_ls, vor_ls, aor_ls = backend.lines_or_info()
        max_mis = np.max(np.abs(por_ls - por_pp))
        self._assert_or_print(max_mis <= self.tol, f"Error: por do not match, maximum absolute error is {max_mis:.2e} MW")
        max_mis = np.max(np.abs(qor_ls - qor_pp))
        self._assert_or_print(max_mis <= self.tol, f"Error: qor do not match, maximum absolute error is {max_mis:.2e} MVAr")
        max_mis = np.max(np.abs(vor_ls - vor_pp))
        self._assert_or_print(max_mis <= self.tol, f"Error: vor do not match, maximum absolute error is {max_mis:.2e} kV")
        max_mis = np.max(np.abs(aor_ls - aor_pp) * 0.001)
        self._assert_or_print(max_mis <= self.tol, f"Error: aor do not match, maximum absolute error is {max_mis:.2e} kA")

        # "II - Check for possible solver issues"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(backend._init_pp_backend._grid, v_debug=True, lightsim2grid=False)
        v_tmp = backend._init_pp_backend._grid.res_bus["vm_pu"].values[:nb_sub] + 0j
        v_tmp *= np.exp(1j * np.pi / 180. * backend._init_pp_backend._grid.res_bus["va_degree"].values[:nb_sub])
        v_tmp = np.concatenate((v_tmp, v_tmp))

        V = backend._grid.ac_pf(v_tmp, 10, 1e-5)
        assert V.shape[0], "? lightsim diverge when initialized with pp final voltage ?"
        backend._grid.tell_solver_need_reset()

        Y_pp = backend._init_pp_backend._grid._ppc["internal"]["Ybus"]
        Sbus = backend._init_pp_backend._grid._ppc["internal"]["Sbus"]
        pv_ = backend._init_pp_backend._grid._ppc["internal"]["pv"]
        pq_ = backend._init_pp_backend._grid._ppc["internal"]["pq"]
        max_iter = 10
        tol_this = 1e-8
        All_Vms = backend._init_pp_backend._grid._ppc["internal"]["Vm_it"]
        AllVas = backend._init_pp_backend._grid._ppc["internal"]["Va_it"]

        for index_V in range(All_Vms.shape[1] - 1, -1, -1):
            nb_iter = All_Vms.shape[1] - 1
            # i check from easiest to hardest, so from the last iteartion of pandapower to the first iteration of pandapower
            # take the same V as pandapower
            V_init = All_Vms[:, index_V] * (np.cos(AllVas[:, index_V]) + 1j * np.sin(AllVas[:, index_V]))
            # V_init *= np.exp(1j * AllVas[:, 0])
            V_init_ref = copy.deepcopy(V_init)
            solver = ClassSolver()
            
            # extract the slack bus
            ref = set(np.arange(V_init.shape[0])) - set(pv_) - set(pq_)
            ref = np.array(list(ref))
            # build the slack weights
            slack_weights = np.zeros(V_init.shape[0])
            slack_weights[ref] = 1.0 / ref.shape[0]

            solver.solve(scipy.sparse.csc_matrix(Y_pp), V_init, Sbus,
                         ref, slack_weights, pv_, pq_, max_iter, tol_this)
            time_for_nr = solver.get_timers()[3]
            if TIMER_INFO:
                print(f"\t Info: Time to perform {nb_iter - index_V} NR iterations for a grid with {nb_sub} "
                      f"buses: {1000. * time_for_nr:.2f}ms")
            error_va = np.abs(solver.get_Va() - np.angle(backend._init_pp_backend._grid._ppc["internal"]["V"]))
            self._assert_or_print(np.max(error_va) <= self.tol, f"Error: VA do not match for iteration {index_V}, maximum absolute " \
                                                 f"error is {np.max(error_va):.5f} rad")

            error_vm = np.abs(np.abs(solver.get_Vm() - np.abs(backend._init_pp_backend._grid._ppc["internal"]["V"])))
            self._assert_or_print(np.max(error_vm) <= self.tol, f"\t Error: VM do not match for iteration {index_V}, maximum absolute " \
                                                 f"error  is {np.max(error_vm):.5f} pu")
            solver.reset()

        # if TIMER_INFO:
        #     print("")

        # 'III - Check the data conversion'
        pp_vect_converter = backend._init_pp_backend._grid._pd2ppc_lookups["bus"][:nb_sub]
        pp_net = backend._init_pp_backend._grid

        # 1) Checking Sbus conversion
        Sbus_pp = backend._init_pp_backend._grid._ppc["internal"]["Sbus"]
        Sbus_pp_right_order = Sbus_pp[pp_vect_converter]
        Sbus_me = backend._grid.get_Sbus_solver()

        # slack bus is not the same
        all_but_slack = np.array(list(set(pv_).union(set(pq_))))
        error_p = np.abs(np.real(Sbus_me[all_but_slack]) - np.real(Sbus_pp_right_order[all_but_slack]))
        assert np.max(error_p) <= self.tol, f"\t Error: P do not match for Sbus, maximum absolute error is " \
                                            f"{np.max(error_p):.5f} MW, \t Error: significative difference for bus " \
                                            f"index (lightsim): {(error_p > self.tol).nonzero()[0]}"

        error_q = np.abs(np.imag(Sbus_me[all_but_slack]) - np.imag(Sbus_pp_right_order[all_but_slack]))
        assert np.max(error_q) <= self.tol, f"\t Error: Q do not match for Sbus, maximum absolute error is " \
                                            f"{np.max(error_q):.5f} MVAr, \t Error: significative difference for bus " \
                                            f"index (lightsim): {(error_q > self.tol).nonzero()[0]}"

        # 2)  Checking Ybus conversion"
        Y_me = backend._grid.get_Ybus_solver()
        Y_pp = backend._init_pp_backend._grid._ppc["internal"]["Ybus"]
        Y_pp_right_order = Y_pp[pp_vect_converter.reshape(nb_sub, 1), pp_vect_converter.reshape(1, nb_sub)]
        error_p = np.abs(np.real(Y_me) - np.real(Y_pp_right_order))
        # arr_ = error_p.toarray()
        # (arr_ >= 1).nonzero()
        self._assert_or_print(np.max(error_p) <= self.tol, f"Error: P do not match for Ybus, maximum absolute error " \
                                            f"is {np.max(error_p):.5f}")

        error_q = np.abs(np.imag(Y_me) - np.imag(Y_pp_right_order))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: Q do not match for Ybus, maximum absolute error is " \
                                            f"{np.max(error_q):.5f}")

        # "IV - Check for the initialization (dc powerflow)"
        # 1) check that the results are same for dc lightsim and dc pandapower
        if isinstance(pp_net["_options"]["init_vm_pu"], str):
            mult_ = 1.0
        else:
            mult_ = pp_net["_options"]["init_vm_pu"]
        Vinit = np.ones(backend.nb_bus_total, dtype=complex) * mult_
        backend._grid.deactivate_result_computation()
        Vdc = backend._grid.dc_pf(Vinit, max_iter, tol_this)
        backend._grid.reactivate_result_computation()
        backend._grid.tell_solver_need_reset()
        Ydc_me = copy.deepcopy(backend._grid.get_dcYbus_solver())
        Sdc_me = copy.deepcopy(backend._grid.get_dcSbus_solver())
        self._assert_or_print(np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) <= 100.*self.tol,\
            f"\t Error for the DC approximation: resulting voltages are different " \
            f"{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f}pu")

        if np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) >= self.tol:
            warnings.warn("\t Warning: maximum difference after DC approximation is "
                          "{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f} which is higher than "
                          "the tolerance (this is just a warning because we noticed this could happen even if the "
                          "results match perfectly. Probably some conversion issue with complex number and "
                          "radian / degree.")
        # "2) check that the Sbus vector is same for PP and lightisim in DC"
        from pandapower.pd2ppc import _pd2ppc
        from pandapower.pf.run_newton_raphson_pf import _get_pf_variables_from_ppci
        from pandapower.pypower.idx_bus import VA, GS
        from pandapower.pypower.makeBdc import makeBdc
        from pandapower.pypower.makeSbus import makeSbus

        pp_net._pd2ppc_lookups = {"bus": np.array([], dtype=int), "ext_grid": np.array([], dtype=int),
                                  "gen": np.array([], dtype=int), "branch": np.array([], dtype=int)}
        # convert pandapower net to ppc
        ppc, ppci = _pd2ppc(pp_net)
        try:
            baseMVA, bus, gen, branch, ref, pv, pq, on, gbus, _, refgen = _get_pf_variables_from_ppci(ppci)
        except ValueError:
            try:
                # change in pandapower 2.12
                baseMVA, bus, gen, branch, svc, tcsc, ref, pv, pq, on, gbus, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
            except ValueError:
                # change in pandapower 2.14
                try:
                    baseMVA, bus, gen, branch, svc, tcsc, ssc, ref, pv, pq, on, gbus, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
                except ValueError:
                    # change in pandapower 3.
                    baseMVA, bus, gen, branch, svc, tcsc, ssc, vsc, ref, pv, pq, *_, V0, ref_gens = _get_pf_variables_from_ppci(ppci, True)
        
        Va0 = bus[:, VA] * (np.pi / 180.)
        try:
            B, Bf, Pbusinj, Pfinj = makeBdc(bus, branch)
        except ValueError:
            # change in pandapower 2.12
            B, Bf, Pbusinj, Pfinj, Cft = makeBdc(bus, branch)
            
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

        self._assert_or_print(np.max(error_p) <= self.tol, f"\t Error: P do not match for Sbus (dc), maximum absolute error is " \
                                            f"{np.max(error_p):.5f} MW, \nError: significative difference for bus " \
                                            f"index (lightsim): {(error_p > self.tol).nonzero()[0]}")

        error_q = np.abs(np.imag(Sdc_me) - np.imag(Pbus_pp_ro))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: Q do not match for Sbus (dc), maximum absolute error is " \
                                            f"{np.max(error_q):.5f} MVAr, \n\t Error: significative difference for " \
                                            f"bus index (lightsim): {(error_q > self.tol).nonzero()[0]}")

        # "3) check that the Ybus matrix is same for PP and lightisim in DC"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(pp_net)
        Ydc_pp = backend._init_pp_backend._grid._ppc["internal"]["Bbus"]
        Ydc_pp_right_order = Ydc_pp[pp_vect_converter.reshape(nb_sub, 1), pp_vect_converter.reshape(1, nb_sub)]
        if False:
            # The way i handle DC PF is different than pandapower
            # for example, matrices are not the same (i put the reactive of PP as "active" and don't 
            # add the "active")
            # so these tests do not really make sense here
            error_p = np.abs(np.real(Ydc_me) - np.real(Ydc_pp_right_order))
            assert np.max(error_p) <= self.tol, f"Error: P do not match for Ybus (dc mode), maximum absolute error " \
                                                f"is {np.max(error_p):.5f}"
            error_q = np.abs(np.imag(Ydc_me) - np.imag(Ydc_pp_right_order))
            assert np.max(error_q) <= self.tol, f"\t Error: Q do not match for Ybus (dc mdoe), maximum absolute error " \
                                                f"is {np.max(error_q):.5f}"

        # "3) check that lightsim ac pf init with pp dc pf give same results (than pp)"
        if isinstance(pp_net["_options"]["init_vm_pu"], str):
            mult_coeff = 1.0
        else:
            mult_coeff = pp_net["_options"]["init_vm_pu"]
        Vinit = np.ones(backend.nb_bus_total, dtype=complex) * mult_coeff
        Vinit[:nb_sub] = V_init_ref[pp_vect_converter]
        conv = backend._grid.ac_pf(Vinit, max_iter, tol_this)
        self._assert_or_print(conv.shape[0] > 0, "\t Error: the lightsim diverge when initialized with pandapower Vinit_dc")
        lpor, lqor, lvor, laor = backend._grid.get_line_res1()
        tpor, tqor, tvor, taor = backend._grid.get_trafo_res1()
        tpex, tqex, tvex, taex = backend._grid.get_trafo_res2()
        nb_trafo = tpor.shape[0]
        nb_powerline = lpor.shape[0]
        p_or_me2 = np.concatenate((lpor, tpor))
        q_or_me2 = np.concatenate((lqor, tqor))
        v_or_me2 = np.concatenate((lvor, tvor))
        a_or_me2 = 1000. * np.concatenate((laor, taor))
        test_ok = True

        max_mis = np.max(np.abs(p_or_me2 - por_pp))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: por do not match, maximum absolute error is {max_mis:.5f} MW")
        max_mis = np.max(np.abs(q_or_me2 - qor_pp))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr")
        max_mis = np.max(np.abs(v_or_me2 - vor_pp))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: vor do not match, maximum absolute error is {max_mis:.5f} kV")
        max_mis = np.max(np.abs(a_or_me2 - aor_pp))
        self._assert_or_print(np.max(error_q) <= self.tol, f"\t Error: aor do not match, maximum absolute error is {max_mis:.5f} A")

        # "V - Check trafo proper conversion to r,x, b"
        
        # my trafo parameters
        converter = PandaPowerConverter()
        converter.set_sn_mva(pp_net.sn_mva)
        converter.set_f_hz(pp_net.f_hz)
        tap_neutral = 1.0 * pp_net.trafo["tap_neutral"].values
        tap_neutral[~np.isfinite(tap_neutral)] = 0.
        if (np.abs(tap_neutral) > 1e-6).any():
            raise RuntimeError("lightsim converter supposes that tap_neutral is 0. for the transformers")
        tap_step_pct = 1.0 * pp_net.trafo["tap_step_percent"].values
        tap_step_pct[~np.isfinite(tap_step_pct)] = 0.
        tap_pos = 1.0 * pp_net.trafo["tap_pos"].values
        tap_pos[~np.isfinite(tap_pos)] = 0.
        shift_ = 1.0 * pp_net.trafo["shift_degree"].values
        shift_[~np.isfinite(shift_)] = 0.
        is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
        is_tap_hv_side[~np.isfinite(is_tap_hv_side)] = True
        # if np.any(pp_net.trafo["tap_phase_shifter"].values):
            # raise RuntimeError("ideal phase shifter are not modeled. Please remove all trafo with "
                            #    "pp_net.trafo[\"tap_phase_shifter\"] set to True.")
        tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
        tap_angles_[~np.isfinite(tap_angles_)] = 0.
        tap_angles_ = np.deg2rad(tap_angles_)
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            trafo_r, trafo_x, trafo_b = \
                converter.get_trafo_param_pp2(tap_step_pct,
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
        else:
            trafo_r, trafo_x, trafo_b = \
                converter.get_trafo_param_pp3(tap_step_pct,
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
                                              True  # trafo_model_is_t
                                              )
            
        # pandapower trafo parameters
        garbage_pp = copy.deepcopy(pp_net)
        if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
            r_pu, x_pu, g_pu, b_pu, g_asym, b_asym, ratio, shift_xx = _calc_branch_values_from_trafo_df(
                    garbage_pp, garbage_pp._ppc)
        else:
            r_pu, x_pu, y_tr_pp, ratio, shift_xx = _calc_branch_values_from_trafo_df(
                    garbage_pp, garbage_pp._ppc)
            g_pu = -y_tr_pp.imag # inverted (conjugate)
            b_pu = y_tr_pp.real  # inverted (conjugate)
            
        
        assert np.allclose(trafo_r, r_pu, atol=self.tol, rtol=self.tol)
        assert np.allclose(trafo_x, x_pu, atol=self.tol, rtol=self.tol)
        assert np.allclose(trafo_b.real, g_pu, atol=self.tol, rtol=self.tol)
        assert np.allclose(trafo_b.imag, b_pu, atol=self.tol, rtol=self.tol)


if __name__ == '__main__':
    unittest.main()
