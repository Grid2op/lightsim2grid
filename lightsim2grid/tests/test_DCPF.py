# Copyright (c) 2020-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from re import M
import time
import unittest
import tempfile
import pandapower as pp
import pandapower.networks as pn
import os
import copy
import numpy as np
import pdb
import warnings


from lightsim2grid import LightSimBackend
try:
    from lightsim2grid.solver import KLUSolver
    ClassSolver = KLUSolver
except ImportError as exc_:
    from lightsim2grid.solver import SparseLUSolver
    ClassSolver = SparseLUSolver


from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION

TIMER_INFO = False  # do i print information regarding computation time


class TestDCPF(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 3e-5  # results are equal if they match up to tol
        self.tol_kcl = 1e-5 # for P = C
        # if CURRENT_PP_VERSION > MAX_PP2_DATAREADER:
        #     self.skipTest("Test not correct: pp changed the way it computed trafo params")

    def test_case14(self):
        case = pn.case14()
        case.name = "case14"
        # self.tol = 2e-3
        self._aux_test(case)

    def test_case14_with_phaseshift(self):
        self.skipTest("pypowsybl and pandapower does not align on DC computation and phase shift... Lightsim2grid choses pypowsybl")
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
        case.name = "case39"
        self.tol = 4e-5
        self._aux_test(case)

    def test_case118(self):
        case = pn.case118()
        case.name = "case118"
        self.tol = 4e-5
        self._aux_test(case)

    # def test_case1888rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case1888rte()
    #     case.name = "case1888rte"
    #     self.tol_kcl = 2e-3
    #     self._aux_test(case)

    def test_case300(self):
        case = pn.case300()
        self._aux_test(case)

    # def test_case2848rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case2848rte()
    #     case.name = "case2848rte"
    #     self.tol = 4e-5
    #     self.tol_kcl = 7e-4
    #     self._aux_test(case)

    # def test_case6470rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6470rte()
    #     case.name = "case6470rte"
    #     self.tol = 2e-4
    #     self.tol_kcl = 7e-4
    #     self._aux_test(case)

    # def test_case6495rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6495rte()
    #     case.name = "case6495rte"
    #     self.tol = 2e-4
    #     self.tol_kcl = 5e-3
    #     self._aux_test(case)

    # def test_case6515rte(self):
    #     # does not work probably with None in converters
    #     case = pn.case6515rte()
    #     case.name = "case6515rte"
    #     self.tol = 2e-4
    #     self.tol_kcl = 9e-3
    #     self._aux_test(case)

    def test_case_illinois200(self):
        case = pn.case_illinois200()
        case.name = "case_illinois200"
        self.tol_kcl = 3e-4
        self._aux_test(case)

    # def test_case9241pegase(self):
    #     # too large for CI
    #     case = pn.case9241pegase()
    #     self.tol = 4e-5
    #     self.tol_kcl = 1.6e-2
    #     self._aux_test(case)
    
    def _aux_make_grid(self, pp_net):
        with tempfile.TemporaryDirectory() as path:
            case_name = os.path.join(path, "this_case.json")
            pp.to_json(pp_net, case_name)

            # real_init_file = pp.from_json(case_name)
            if MAX_PP2_DATAREADER < CURRENT_PP_VERSION:
                loader_kwargs = {"pp_orig_file": "pandapower_v3"}
            else:
                loader_kwargs = {"pp_orig_file": "pandapower_v2"}
            backend = LightSimBackend(loader_kwargs=loader_kwargs)
            
            type(backend)._clear_grid_dependant_class_attributes()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                type(backend).env_name = pp_net.name if pp_net.name != "" else "_test"
                backend.load_grid(case_name)
                backend.assert_grid_correct()
                # backend._init_pp_backend.assert_grid_correct()
                
        # first i deactivate all slack bus in pp that are connected but not handled in ls
        backend._init_pp_backend._grid.ext_grid["in_service"] = False
        backend._init_pp_backend._grid.ext_grid.loc[0, "in_service"] = True
        backend._init_pp_backend._grid.gen["slack"] = False
        return backend

    def _aux_test(self, pn_net):
        backend = self._aux_make_grid(pn_net)
        nb_sub = backend.n_sub
        pp_net = backend._init_pp_backend._grid
        conv, exc_ = backend.runpf(is_dc=True)
        # print(f"{backend.n_sub}: {backend._timer_solver = }")
        # PTDF = backend._grid.get_ptdf_solver()
        # print(f"{backend.n_sub}: {backend._grid.get_dc_solver().get_timers_ptdf_lodf()[0] = }")
        # nb_powerflow = 1_000
        # Sbuses = np.vstack([backend._grid.get_dcSbus_solver().real for _ in range(nb_powerflow)]).T.copy()
        # beg_ = time.perf_counter()
        # PTDF @ Sbuses
        # end_ = time.perf_counter()
        # LODF = backend._grid.get_lodf()
        # print(f"{backend.n_sub}: {backend._grid.get_dc_solver().get_timers_ptdf_lodf()[1] = }")
        # print(f"time {nb_powerflow} powerflows: {end_-beg_:.2e}")
        conv_pp, exc_pp = backend._init_pp_backend.runpf(is_dc=True)
        assert conv_pp, "Error: pandapower do not converge, impossible to perform the necessary checks"
        assert conv, f"Error: lightsim do not converge with error: {exc_}"
        
        # check input data
        Ybus =  backend._grid.get_dcYbus_solver()
        Bbus = pp_net._ppc["internal"]["Bbus"]
        assert np.abs(Ybus - Bbus).max() <= self.tol, f"max error for Ybus {np.abs(Ybus - Bbus).max()}"
        
        from pandapower.pf.run_dc_pf import (makeSbus, GS)
        baseMVA = pp_net._ppc["internal"]['baseMVA']
        bus = pp_net._ppc["internal"]['bus']
        gen =  pp_net._ppc["internal"]['gen']
        Pbusinj =  pp_net._ppc["internal"]['Pbusinj']
        Sbus = backend._grid.get_dcSbus_solver()
        Pbus = makeSbus(baseMVA, bus, gen) - Pbusinj - bus[:, GS] / baseMVA
        slack_bus = pp_net.ext_grid.iloc[0]["bus"]
        mask_not_slack = np.ones(Sbus.shape[0], dtype=bool)
        mask_not_slack[slack_bus] = False
        # python -m unittest test_DCPF.TestDCPF.test_case6495rte
        # test_case6495rte
        # diff > 1e-5
        # (np.abs(Sbus.real - Pbus.real) >= 1e-5).nonzero()
        # array([6077, 6161, 6305, 6306, 6307, 6308])
        # 6077 is the slack so ok !
        # large_error = np.array([6161, 6305, 6306, 6307, 6308], dtype=int)
        # Pbusinj[large_error] => all 0
        # (bus[:, GS] / baseMVA)[large_error] => all 0
        # Sbus_tmp = makeSbus(baseMVA, bus, gen)
        # pp_net.trafo["lv_bus"] == 6160
        # pp_net.trafo.loc[pp_net.trafo["lv_bus"] == 6160]
        assert np.abs(Sbus[mask_not_slack].real - Pbus[mask_not_slack].real).max() <= self.tol, f"max error in Sbus {np.abs(Sbus[mask_not_slack].real - Pbus[mask_not_slack].real).max()}"

        por_pp, qor_pp, vor_pp, aor_pp = copy.deepcopy(backend._init_pp_backend.lines_or_info())
        pex_pp, qex_pp, vex_pp, aex_pp = copy.deepcopy(backend._init_pp_backend.lines_ex_info())
        load_p_pp, load_q_pp, load_v_pp = copy.deepcopy(backend._init_pp_backend.loads_info())
        gen_p_pp, gen_q_pp, gen_v_pp = copy.deepcopy(backend._init_pp_backend.generators_info())
        sh_p_pp, sh_q_pp, sh_v_pp, *_ = copy.deepcopy(backend._init_pp_backend.shunt_info())
        sgen_p_pp = copy.deepcopy(backend._init_pp_backend._grid.res_sgen["p_mw"].values)
        init_gen_p = copy.deepcopy(backend._init_pp_backend._grid.gen["p_mw"].values)
        init_load_p = copy.deepcopy(backend._init_pp_backend._grid.load["p_mw"].values)
        init_sgen_p = copy.deepcopy(backend._init_pp_backend._grid.sgen["p_mw"].values)

        # I- Check for divergence and equality of flows"
        por_ls, qor_ls, vor_ls, aor_ls = backend.lines_or_info()
        big_err_lid = np.where(np.abs(por_ls - por_pp) > 5000)[0]
        backend.line_ex_to_subid[big_err_lid]
        psub_ls, qsub_ls, pbus_ls, qbus_ls, diff_v_bus_ls = backend.check_kirchhoff()
        # below it does not work due to a bug fixed in dev_1.8.2 (after 1.8.2.dev4)
        # psub_pp, qsub_pp, pbus_pp, qbus_pp, diff_v_bus_pp = backend._init_pp_backend.check_kirchoff()
        
        # check voltages
        Va_pp = pp_net.res_bus["va_degree"].values[:nb_sub]
        Va_ls = np.rad2deg(np.angle(backend.V[:nb_sub]))
        assert np.abs(Va_pp - Va_ls).max() <= self.tol, f"max error for voltages {np.abs(Va_pp - Va_ls).max()}"
        
        line_or_theta_pp, line_ex_theta_pp, *_ = backend._init_pp_backend.get_theta()
        line_or_theta_ls, line_ex_theta_ls, *_ = backend.get_theta()
        assert np.all(np.abs(line_or_theta_ls - line_or_theta_pp) <= self.tol), "error in voltage angles (theta_or)"
        assert np.all(np.abs(line_ex_theta_ls - line_ex_theta_pp) <= self.tol), "error in voltage angles (theta_ex)"
        
        max_mis = np.max(np.abs(por_ls - por_pp))
        # max_error_id = np.argmax(np.abs(por_ls - por_pp))
        # nb_line = pp_net.line.shape[0]
        # trafo_id = max_error_id - nb_line
        # por_ls[max_error_id]
        # np.abs(por_ls - por_pp)[25 + nb_line]
        # backend._grid.get_trafos()[trafo_id]
        # pp_net.trafo.iloc[trafo_id]
        # check the Ybus for DC
        # check the voltage angles

        # check flows
        assert max_mis <= self.tol, f"Error: por do not match, maximum absolute error is {max_mis:.5f} MW"
        max_mis = np.max(np.abs(qor_ls - qor_pp))
        assert max_mis <= self.tol, f"Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(vor_ls - vor_pp))
        assert max_mis <= self.tol, f"Error: vor do not match, maximum absolute error is {max_mis:.5f} kV"
        max_mis = np.max(np.abs(aor_ls - aor_pp))
        assert max_mis <= self.tol * 1e3, f"Error: aor do not match, maximum absolute error is {max_mis:.5f} kA"
        
        load_p, load_q, load_v = backend.loads_info()
        max_mis = np.max(np.abs(load_p - load_p_pp))
        assert max_mis <= self.tol, f"Error: load_p do not match, maximum absolute error is {max_mis:.5f} MW"
        # PP does not set "load_q" to 0. in DC
        # max_mis = np.max(np.abs(load_q - load_q_pp))
        # assert max_mis <= self.tol, f"Error: load_q do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(load_v - load_v_pp))
        assert max_mis <= self.tol, f"Error: load_v do not match, maximum absolute error is {max_mis:.5f} kV"

        gen_p, gen_q, gen_v = backend.generators_info()
        sgen_p, sgen_q, sgen_v = backend._grid.get_sgens_res()
        # pandapower is not correct on dc...
        # max_mis = np.max(np.abs(gen_p - gen_p_pp))
        # assert max_mis <= self.tol, f"Error: gen_p do not match, maximum absolute error is {max_mis:.5f} MW"
        
        # np.sum(gen_p_pp) + np.sum(sgen_p_pp) - np.sum(load_p_pp)
        # pandapower also does weird things in dc for gen_q... lightsim2grid puts everything at 0.
        # max_mis = np.max(np.abs(gen_q - gen_q_pp))
        # assert max_mis <= self.tol, f"Error: gen_q do not match, maximum absolute error is {max_mis:.5f} MVAr"
        assert np.max(np.abs(gen_q)) <= self.tol
        if sgen_q.size:
            assert np.max(np.abs(sgen_q)) <= self.tol

        max_mis = np.max(np.abs(gen_v - gen_v_pp))
        assert max_mis <= self.tol, f"Error: gen_v do not match, maximum absolute error is {max_mis:.5f} kV"

        sh_p, sh_q, sh_v, *_ = backend.shunt_info()
        if sh_p.size > 0:
            max_mis = np.max(np.abs(sh_p - sh_p_pp))
            assert max_mis <= self.tol, f"Error: sh_p do not match, maximum absolute error is {max_mis:.5f} MW"
        mismatch = np.sum(gen_p) + np.sum(sgen_p) - np.sum(load_p) - sh_p.sum()
        assert abs(mismatch) <= self.tol_kcl, f"max sum P = sum C mismatch: {mismatch} MW"
        # max_mis = np.max(np.abs(sh_q - sh_q_pp))
        # assert max_mis <= self.tol, f"Error: sh_q do not match, maximum absolute error is {max_mis:.5f} MVAr"
        # again pandapower does weird stuff in dc...

        # assert np.max(np.abs(sh_q)) <= self.tol
        # max_mis = np.max(np.abs(sh_v - sh_v_pp))
        # assert max_mis <= self.tol, f"Error: sh_v do not match, maximum absolute error is {max_mis:.5f} kV"
        # again, pandapower put nan for the voltages...


class TestDCPF_LODF(TestDCPF):
    """test all the `get_xxx` and `get_xxx_solver` can be accessed
    without crash.
    
    Also tests the LODF formula
    """
    def _aux_aux_test_accessors(self, gridmodel):
        id_dc_solver_to_me = np.array(gridmodel.id_dc_solver_to_me())
        PTDF = gridmodel.get_ptdf()        
        assert PTDF.shape == (len(gridmodel.get_lines()) + len(gridmodel.get_trafos()), gridmodel.total_bus())        
        PTDF_solver = gridmodel.get_ptdf_solver()
        assert PTDF_solver.shape == (len(gridmodel.get_lines()) + len(gridmodel.get_trafos()), gridmodel.nb_connected_bus()) 
        assert (PTDF[:, id_dc_solver_to_me] == PTDF_solver).all()      
        with self.assertRaises(RuntimeError):
            Ybus = gridmodel.get_Ybus()              
        Ybus_solver = gridmodel.get_Ybus_solver()
        assert Ybus_solver.shape == (0, 0)
        dcYbus = gridmodel.get_dcYbus()     
        assert dcYbus.shape == (gridmodel.total_bus(), gridmodel.total_bus())    
        dcYbus_solver = gridmodel.get_dcYbus_solver()
        assert dcYbus_solver.shape == (gridmodel.nb_connected_bus(), gridmodel.nb_connected_bus())   
        assert (dcYbus[id_dc_solver_to_me.reshape(-1,1), id_dc_solver_to_me.reshape(1,-1)] != dcYbus_solver).nnz == 0    
        with self.assertRaises(RuntimeError):
            Sbus = gridmodel.get_Sbus()                
        Sbus_solver = gridmodel.get_Sbus_solver()
        assert Sbus_solver.shape == (0, )
        dcSbus = gridmodel.get_dcSbus()      
        assert dcSbus.shape == (gridmodel.total_bus(), )            
        dcSbus_solver = gridmodel.get_dcSbus_solver()
        assert dcSbus_solver.shape == (gridmodel.nb_connected_bus(), )   
        assert (dcSbus[id_dc_solver_to_me] == dcSbus_solver).all()    
        pv = gridmodel.get_pv()                  
        pv_solver = gridmodel.get_pv_solver()
        assert len(pv_solver) == len(np.unique([el.bus_id for el in gridmodel.get_generators() if not el.is_slack]))
        assert len(pv) == len(pv_solver)
        assert (np.sort(pv) == np.sort(pv_solver)).all()
        assert (id_dc_solver_to_me[pv_solver] == pv).all()
        pq = gridmodel.get_pq()                  
        pq_solver = gridmodel.get_pq_solver()
        assert (np.sort(pq) == np.sort(pq_solver)).all()
        assert (id_dc_solver_to_me[pq_solver] == pq).all()
        with self.assertRaises(RuntimeError):
            slack_ids = gridmodel.get_slack_ids()           
        slack_ids_solver = gridmodel.get_slack_ids_solver()
        assert slack_ids_solver.shape == (0, )
        slack_ids_dc = gridmodel.get_slack_ids_dc()    
        this_slack = np.sort(np.unique([el.bus_id for el in gridmodel.get_generators() if el.is_slack]))
        assert (np.sort(slack_ids_dc) == this_slack).all()   
        slack_ids_dc_solver = gridmodel.get_slack_ids_dc_solver()     
        assert (id_dc_solver_to_me[slack_ids_dc_solver] == this_slack).all() 
        slack_weights = gridmodel.get_slack_weights()    
        assert slack_weights.shape == (gridmodel.total_bus(), )      
        slack_weights_solver = gridmodel.get_slack_weights_solver()
        assert slack_weights_solver.shape == (gridmodel.nb_connected_bus(), )    
        assert (slack_weights[id_dc_solver_to_me] == slack_weights_solver).all()
        Bf = gridmodel.get_Bf()    
        assert Bf.shape == (len(gridmodel.get_lines()) + len(gridmodel.get_trafos()), gridmodel.total_bus())
        Bf_solver = gridmodel.get_Bf_solver()
        assert Bf_solver.shape == (len(gridmodel.get_lines()) + len(gridmodel.get_trafos()), gridmodel.nb_connected_bus())
        
    def _aux_test(self, pn_net):
        backend = self._aux_make_grid(pn_net)
        nb_sub = backend.n_sub
        pp_net = backend._init_pp_backend._grid
        conv, exc_ = backend.runpf(is_dc=True)
        assert conv
        gridmodel = backend._grid
        
        # test I can access all that without crash
        self._aux_aux_test_accessors(gridmodel)
        # test the lodf formula
        LODF_mat = gridmodel.get_lodf()      
        lor_p, *_ = gridmodel.get_line_res1()
        tor_p, *_ = gridmodel.get_trafo_res1()
        init_powerflow = np.concatenate((lor_p, tor_p))
        prng = np.random.default_rng(0)
        nb_powerline = len(gridmodel.get_lines())
        total_nb = nb_powerline + len(gridmodel.get_trafos())
        for i, l_id in enumerate(prng.integers(0, total_nb, size=10)):
            if pp_net.name == "case1888rte" and l_id == 103:
                # does not really work, LODF diverges I don't really know why
                # maybe something to do with (abs(diag_coeff - 1.) > tol_equal_float) cpp side
                # see BaseDCAlgo.tpp
                continue
            por_lodf = init_powerflow + LODF_mat[:, l_id] * init_powerflow[l_id]
            gridmodel_cpy = gridmodel.copy()
            if l_id >= nb_powerline:
                gridmodel_cpy.deactivate_trafo(l_id - nb_powerline)
            else:
                gridmodel_cpy.deactivate_powerline(l_id)
            Vinit = np.ones(gridmodel_cpy.total_bus(), dtype=complex)
            Vdc = gridmodel_cpy.dc_pf(Vinit, 1, 1e-8)
            if Vdc.shape == (0, ):
                # divergence, so it should be Nan as per LODF
                assert (~np.isfinite(por_lodf)).all(), f"error for line / trafo {l_id} (iter {i})"
                continue
            else:
                # no divergence, so it should be finite as per LODF
                assert np.isfinite(LODF_mat[:, l_id]).all(), f"error for line / trafo {l_id} : divergence for LODF but not for regular LF(iter {i})"
                
            # convergence, flows should match
            lor_p_tmp, *_ = gridmodel_cpy.get_line_res1()
            tor_p_tmp, *_ = gridmodel_cpy.get_trafo_res1()
            real_val = np.concatenate((lor_p_tmp, tor_p_tmp))
            assert (np.abs(por_lodf - real_val) <= 1e-6).all(), f"error for line / trafo {l_id} (iter {i}): {por_lodf - real_val}"
        
        
if __name__ == "__main__":
    unittest.main()
