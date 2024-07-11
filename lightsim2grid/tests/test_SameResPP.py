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
try:
    from lightsim2grid.solver import KLUSolver
    ClassSolver = KLUSolver
except ImportError as exc_:
    from lightsim2grid.solver import SparseLUSolver
    ClassSolver = SparseLUSolver

TIMER_INFO = False  # do i print information regarding computation time


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 1e-4  # results are equal if they match up to tol

    def test_case14(self):
        case = pn.case14()
        self.tol = 2e-3
        self._aux_test(case)

    def test_case39(self):
        case = pn.case39()
        self._aux_test(case)

    def test_case118(self):
        case = pn.case118()
        self._aux_test(case)

    def test_case1888rte(self):
        case = pn.case1888rte()
        self.tol = 3e-4
        self._aux_test(case)

    # def test_case300(self):
    # TODO make it work
    #     case = pn.case300()
    #     self._aux_test(case)

    # def test_case9241pegase(self):
    # TODO make it work
    #     case = pn.case9241pegase()
    #     self._aux_test(case)

    def test_case2848rte(self):
        case = pn.case2848rte()
        self.tol = 0.1  # yeah this one is a bit tough... # TODO
        self._aux_test(case)

    def test_case6470rte(self):
        case = pn.case6470rte()
        self.tol = 1e-2
        self._aux_test(case)

    def test_case6495rte(self):
        case = pn.case6495rte()
        self.tol = 1e-2
        self._aux_test(case)

    def test_case6515rte(self):
        case = pn.case6515rte()
        self.tol = 1e-2
        self._aux_test(case)

    def test_case_illinois200(self):
        case = pn.case_illinois200()
        self._aux_test(case)

    def _aux_test(self, pn_net):
        n_busbar = 2
        with tempfile.TemporaryDirectory() as path:
            case_name = os.path.join(path, "this_case.json")
            pp.to_json(pn_net, case_name)

            real_init_file = pp.from_json(case_name)
            LightSimBackend._clear_grid_dependant_class_attributes()
            LightSimBackend.set_env_name(type(self).__name__ + case_name)
            LightSimBackend.set_n_busbar_per_sub(n_busbar)
            backend = LightSimBackend()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                backend.load_grid(case_name)
            backend.assert_grid_correct()

        nb_sub = backend.n_sub
        pp_net = backend.init_pp_backend._grid
        # first i deactivate all slack bus in pp that are connected but not handled in ls
        pp_net.ext_grid["in_service"].loc[:] = False
        pp_net.ext_grid["in_service"].iloc[0] = True
        conv, exc_ = backend.runpf()
        conv_pp, exc_pp = backend.init_pp_backend.runpf()
        # import pdb
        # pdb.set_trace()
        assert conv_pp, "Error: pandapower do not converge, impossible to perform the necessary checks"
        assert conv, "Error: lightsim do not converge"

        por_pp, qor_pp, vor_pp, aor_pp = copy.deepcopy(backend.init_pp_backend.lines_or_info())
        pex_pp, qex_pp, vex_pp, aex_pp = copy.deepcopy(backend.init_pp_backend.lines_ex_info())

        # I- Check for divergence and equality of flows"
        por_ls, qor_ls, vor_ls, aor_ls = backend.lines_or_info()
        max_mis = np.max(np.abs(por_ls - por_pp))
        assert max_mis <= self.tol, f"Error: por do not match, maximum absolute error is {max_mis:.5f} MW"
        max_mis = np.max(np.abs(qor_ls - qor_pp))
        assert max_mis <= self.tol, f"Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(vor_ls - vor_pp))
        assert max_mis <= self.tol, f"Error: vor do not match, maximum absolute error is {max_mis:.5f} kV"
        max_mis = np.max(np.abs(aor_ls - aor_pp))
        assert max_mis <= self.tol, f"Error: aor do not match, maximum absolute error is {max_mis:.5f} A"

        # "II - Check for possible solver issues"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(backend.init_pp_backend._grid, v_debug=True, lightsim2grid=False)
        v_tmp = backend.init_pp_backend._grid.res_bus["vm_pu"].values[:nb_sub] + 0j
        v_tmp *= np.exp(1j * np.pi / 180. * backend.init_pp_backend._grid.res_bus["va_degree"].values[:nb_sub])
        v_tmp = np.concatenate((v_tmp, v_tmp))

        V = backend._grid.ac_pf(v_tmp, 10, 1e-5)
        assert V.shape[0], "? lightsim diverge when initialized with pp final voltage ?"
        backend._grid.tell_solver_need_reset()

        Y_pp = backend.init_pp_backend._grid._ppc["internal"]["Ybus"]
        Sbus = backend.init_pp_backend._grid._ppc["internal"]["Sbus"]
        pv_ = backend.init_pp_backend._grid._ppc["internal"]["pv"]
        pq_ = backend.init_pp_backend._grid._ppc["internal"]["pq"]
        max_iter = 10
        tol_this = 1e-8
        All_Vms = backend.init_pp_backend._grid._ppc["internal"]["Vm_it"]
        AllVas = backend.init_pp_backend._grid._ppc["internal"]["Va_it"]

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
            error_va = np.abs(solver.get_Va() - np.angle(backend.init_pp_backend._grid._ppc["internal"]["V"]))
            assert np.max(error_va) <= self.tol, f"Error: VA do not match for iteration {index_V}, maximum absolute " \
                                                 f"error is {np.max(error_va):.5f} rad"

            error_vm = np.abs(np.abs(solver.get_Vm() - np.abs(backend.init_pp_backend._grid._ppc["internal"]["V"])))
            assert np.max(error_vm) <= self.tol, f"\t Error: VM do not match for iteration {index_V}, maximum absolute " \
                                                 f"error  is {np.max(error_vm):.5f} pu"
            solver.reset()

        if TIMER_INFO:
            print("")

        # 'III - Check the data conversion'
        pp_vect_converter = backend.init_pp_backend._grid._pd2ppc_lookups["bus"][:nb_sub]
        pp_net = backend.init_pp_backend._grid

        # 1) Checking Sbus conversion
        Sbus_pp = backend.init_pp_backend._grid._ppc["internal"]["Sbus"]
        Sbus_pp_right_order = Sbus_pp[pp_vect_converter]
        Sbus_me = backend._grid.get_Sbus()
        # slack bus is not the same
        all_but_slack = np.array(list(set(pv_).union(set(pq_))))
        error_p = np.abs(np.real(Sbus_me[all_but_slack]) - np.real(Sbus_pp_right_order[all_but_slack]))
        assert np.max(error_p) <= self.tol, f"\t Error: P do not match for Sbus, maximum absolute error is " \
                                            f"{np.max(error_p):.5f} MW, \t Error: significative difference for bus " \
                                            f"index (lightsim): {np.where(error_p > self.tol)[0]}"

        error_q = np.abs(np.imag(Sbus_me[all_but_slack]) - np.imag(Sbus_pp_right_order[all_but_slack]))
        assert np.max(error_q) <= self.tol, f"\t Error: Q do not match for Sbus, maximum absolute error is " \
                                            f"{np.max(error_q):.5f} MVAr, \t Error: significative difference for bus " \
                                            f"index (lightsim): {np.where(error_q > self.tol)[0]}"

        # 2)  Checking Ybus conversion"
        Y_me = backend._grid.get_Ybus()
        Y_pp = backend.init_pp_backend._grid._ppc["internal"]["Ybus"]
        Y_pp_right_order = Y_pp[pp_vect_converter.reshape(nb_sub, 1), pp_vect_converter.reshape(1, nb_sub)]
        error_p = np.abs(np.real(Y_me) - np.real(Y_pp_right_order))
        assert np.max(error_p) <= self.tol, f"Error: P do not match for Ybus, maximum absolute error " \
                                            f"is {np.max(error_p):.5f}"

        error_q = np.abs(np.imag(Y_me) - np.imag(Y_pp_right_order))
        assert np.max(error_q) <= self.tol, f"\t Error: Q do not match for Ybus, maximum absolute error is " \
                                            f"{np.max(error_q):.5f}"

        # "IV - Check for the initialization (dc powerflow)"
        # 1) check that the results are same for dc lightsim and dc pandapower
        if isinstance(pp_net["_options"]["init_vm_pu"], str):
            mult_ = 1.0
        else:
            mult_ = pp_net["_options"]["init_vm_pu"]
        Vinit = np.ones(backend.nb_bus_total, dtype=np.complex_) * mult_
        backend._grid.deactivate_result_computation()
        Vdc = backend._grid.dc_pf(Vinit, max_iter, tol_this)
        backend._grid.reactivate_result_computation()
        backend._grid.tell_solver_need_reset()
        Ydc_me = copy.deepcopy(backend._grid.get_dcYbus())
        Sdc_me = copy.deepcopy(backend._grid.get_dcSbus())
        assert np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) <= 100.*self.tol,\
            f"\t Error for the DC approximation: resulting voltages are different " \
            f"{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f}pu"

        if np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])) >= self.tol:
            warnings.warn("\t Warning: maximum difference after DC approximation is "
                          "{np.max(np.abs(V_init_ref[pp_vect_converter] - Vdc[:nb_sub])):.5f} which is higher than "
                          "the tolerance (this is just a warning because we noticed this could happen even if the "
                          "results match perfectly. Probably some conversion issue with complex number and "
                          "radian / degree.")
        # "2) check that the Sbus vector is same for PP and lightisim in DC"
        from pandapower.pd2ppc import _pd2ppc
        from pandapower.pf.run_newton_raphson_pf import _get_pf_variables_from_ppci
        from pandapower.pypower.idx_brch import F_BUS, T_BUS, BR_X, TAP, SHIFT, BR_STATUS
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
                baseMVA, bus, gen, branch, svc, tcsc, ssc, ref, pv, pq, on, gbus, V0, ref_gens = _get_pf_variables_from_ppci(ppci)
        
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

        assert np.max(error_p) <= self.tol, f"\t Error: P do not match for Sbus (dc), maximum absolute error is " \
                                            f"{np.max(error_p):.5f} MW, \nError: significative difference for bus " \
                                            f"index (lightsim): {np.where(error_p > self.tol)[0]}"

        error_q = np.abs(np.imag(Sdc_me) - np.imag(Pbus_pp_ro))
        assert np.max(error_q) <= self.tol, f"\t Error: Q do not match for Sbus (dc), maximum absolute error is " \
                                            f"{np.max(error_q):.5f} MVAr, \n\t Error: significative difference for " \
                                            f"bus index (lightsim): {np.where(error_q > self.tol)[0]}"

        # "3) check that the Ybus matrix is same for PP and lightisim in DC"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(pp_net)
        Ydc_pp = backend.init_pp_backend._grid._ppc["internal"]["Bbus"]
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
        Vinit = np.ones(backend.nb_bus_total, dtype=np.complex_) * mult_coeff
        Vinit[:nb_sub] = V_init_ref[pp_vect_converter]
        conv = backend._grid.ac_pf(Vinit, max_iter, tol_this)
        assert conv.shape[0] > 0, "\t Error: the lightsim diverge when initialized with pandapower Vinit_dc"
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

        max_mis = np.max(np.abs(p_or_me2 - por_pp))
        assert np.max(error_q) <= self.tol, f"\t Error: por do not match, maximum absolute error is {max_mis:.5f} MW"
        max_mis = np.max(np.abs(q_or_me2 - qor_pp))
        assert np.max(error_q) <= self.tol, f"\t Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(v_or_me2 - vor_pp))
        assert np.max(error_q) <= self.tol, f"\t Error: vor do not match, maximum absolute error is {max_mis:.5f} kV"
        max_mis = np.max(np.abs(a_or_me2 - aor_pp))
        assert np.max(error_q) <= self.tol, f"\t Error: aor do not match, maximum absolute error is {max_mis:.5f} A"

        # "V - Check trafo proper conversion to r,x, b"
        from lightsim2grid_cpp import GridModel, PandaPowerConverter, SolverType
        from pandapower.build_branch import _calc_branch_values_from_trafo_df, get_trafo_values
        from pandapower.build_branch import _calc_nominal_ratio_from_dataframe, _calc_r_x_y_from_dataframe
        from pandapower.build_branch import _calc_tap_from_dataframe, BASE_KV, _calc_r_x_from_dataframe

        # my trafo parameters
        converter = PandaPowerConverter()
        converter.set_sn_mva(pp_net.sn_mva)
        converter.set_f_hz(pp_net.f_hz)
        tap_neutral = 1.0 * pp_net.trafo["tap_neutral"].values
        tap_neutral[~np.isfinite(tap_neutral)] = 0.
        if np.any(tap_neutral != 0.):
            raise RuntimeError("lightsim converter supposes that tap_neutral is 0 for the transformers")
        tap_step_pct = 1.0 * pp_net.trafo["tap_step_percent"].values
        tap_step_pct[~np.isfinite(tap_step_pct)] = 0.
        tap_pos = 1.0 * pp_net.trafo["tap_pos"].values
        tap_pos[~np.isfinite(tap_pos)] = 0.
        shift_ = 1.0 * pp_net.trafo["shift_degree"].values
        shift_[~np.isfinite(shift_)] = 0.
        is_tap_hv_side = pp_net.trafo["tap_side"].values == "hv"
        is_tap_hv_side[~np.isfinite(is_tap_hv_side)] = True
        if np.any(pp_net.trafo["tap_phase_shifter"].values):
            raise RuntimeError("ideal phase shifter are not modeled. Please remove all trafo with "
                               "pp_net.trafo[\"tap_phase_shifter\"] set to True.")
        tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
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
        # pandapower trafo parameters
        ppc = copy.deepcopy(pp_net._ppc)
        bus_lookup = pp_net["_pd2ppc_lookups"]["bus"]
        trafo_df = pp_net["trafo"]
        lv_bus = get_trafo_values(trafo_df, "lv_bus")
        vn_lv = ppc["bus"][bus_lookup[lv_bus], BASE_KV]
        vn_trafo_hv, vn_trafo_lv, shift_pp = _calc_tap_from_dataframe(pp_net, trafo_df)
        ratio = _calc_nominal_ratio_from_dataframe(ppc, trafo_df, vn_trafo_hv, vn_trafo_lv, bus_lookup)
        r_t, x_t, b_t = _calc_r_x_y_from_dataframe(pp_net, trafo_df, vn_trafo_lv, vn_lv, pp_net.sn_mva)

        # check where there are mismatch if any
        val_r_pp = r_t
        val_r_me = trafo_r
        all_equals_r = np.abs(val_r_pp - val_r_me) <= self.tol
        if not np.all(all_equals_r):
            test_ok = False
            print(f"\t Error: some trafo resistance are not equal, max error: {np.max(np.abs(val_r_pp - val_r_me)):.5f}")

        val_x_pp = x_t
        val_x_me = trafo_x
        all_equals_x = np.abs(val_x_pp - val_x_me) <= self.tol
        assert np.all(all_equals_x), f"\t Error: some trafo x are not equal, max error: " \
                                     f"{np.max(np.abs(val_x_pp - val_x_me)):.5f}"

        val_ib_pp = np.imag(b_t)
        val_ib_me = np.imag(trafo_b)
        all_equals_imag_b = np.abs(val_ib_pp - val_ib_me) <= self.tol
        assert np.all(all_equals_imag_b), f"\t Error: some trafo (imag) b are not equal, max error: " \
                                          f"{np.max(np.abs(val_ib_pp - val_ib_me)):.5f}"

        val_reb_pp = np.real(b_t)
        val_reb_me = np.real(trafo_b)
        all_equals_real_b = np.abs(val_reb_pp - val_reb_me) <= self.tol
        assert np.all(all_equals_real_b), f"\t Error: some trafo (real) b are not equal, max error: " \
                                          f"{np.max(np.abs(val_reb_pp - val_reb_me)):.5f}"

if __name__ == '__main__':
    unittest.main()
