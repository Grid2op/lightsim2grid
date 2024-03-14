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
from pandapower.std_types import parameter_from_std_type
import scipy
import pdb
import warnings

from lightsim2grid import LightSimBackend
try:
    from lightsim2grid.solver import KLUSolver
    ClassSolver = KLUSolver
except ImportError as exc_:
    from lightsim2grid.solver import SparseLUSolver
    ClassSolver = SparseLUSolver

import pdb
TIMER_INFO = False  # do i print information regarding computation time


class TestDCPF(unittest.TestCase):
    def setUp(self) -> None:
        self.tol = 1e-4  # results are equal if they match up to tol
        self.tol_big = 0.01  # for P = C

    def test_case14(self):
        case = pn.case14()
        self.tol = 2e-3
        self._aux_test(case)

    def test_case14_with_phaseshift(self):
        case = pn.case14()
        # case.trafo.iloc[[0]]["shift_degree"] = -5.0
        
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
        self.tol = 2e-3
        self._aux_test(case)

    def test_case39(self):
        case = pn.case39()
        self.tol = 3e-4
        self._aux_test(case)

    def test_case118(self):
        case = pn.case118()
        self._aux_test(case)

    def test_case1888rte(self):
        case = pn.case1888rte()
        self.tol = 3e-4
        self._aux_test(case)

    # def test_case300(self):
    #     # issue with Bbus / Ybus matrix
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
        self.tol_big = 0.1  # for P = C
        self.tol = 1e-2
        self._aux_test(case)

    def test_case6495rte(self):
        case = pn.case6495rte()
        self.tol = 1e-2
        self._aux_test(case)

    def test_case6515rte(self):
        case = pn.case6515rte()
        self.tol_big = 0.1  # for P = C
        self.tol = 1e-2
        self._aux_test(case)

    def test_case_illinois200(self):
        case = pn.case_illinois200()
        self.tol = 3e-4
        self._aux_test(case)

    def _aux_test(self, pn_net):
        with tempfile.TemporaryDirectory() as path:
            case_name = os.path.join(path, "this_case.json")
            pp.to_json(pn_net, case_name)

            real_init_file = pp.from_json(case_name)
            backend = LightSimBackend()
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                backend.load_grid(case_name)
                backend.assert_grid_correct()
                # backend.init_pp_backend.assert_grid_correct()
                
        nb_sub = backend.n_sub
        pp_net = backend.init_pp_backend._grid
        # first i deactivate all slack bus in pp that are connected but not handled in ls
        pp_net.ext_grid["in_service"].loc[:] = False
        pp_net.ext_grid["in_service"].iloc[0] = True
        conv, exc_ = backend.runpf(is_dc=True)
        conv_pp, exc_pp = backend.init_pp_backend.runpf(is_dc=True)

        assert conv_pp, "Error: pandapower do not converge, impossible to perform the necessary checks"
        assert conv, f"Error: lightsim do not converge with error: {exc_}"

        por_pp, qor_pp, vor_pp, aor_pp = copy.deepcopy(backend.init_pp_backend.lines_or_info())
        pex_pp, qex_pp, vex_pp, aex_pp = copy.deepcopy(backend.init_pp_backend.lines_ex_info())
        load_p_pp, load_q_pp, load_v_pp = copy.deepcopy(backend.init_pp_backend.loads_info())
        gen_p_pp, gen_q_pp, gen_v_pp = copy.deepcopy(backend.init_pp_backend.generators_info())
        sh_p_pp, sh_q_pp, sh_v_pp, *_ = copy.deepcopy(backend.init_pp_backend.shunt_info())
        sgen_p_pp = copy.deepcopy(backend.init_pp_backend._grid.res_sgen["p_mw"].values)
        init_gen_p = copy.deepcopy(backend.init_pp_backend._grid.gen["p_mw"].values)
        init_load_p = copy.deepcopy(backend.init_pp_backend._grid.load["p_mw"].values)
        init_sgen_p = copy.deepcopy(backend.init_pp_backend._grid.sgen["p_mw"].values)

        # I- Check for divergence and equality of flows"
        por_ls, qor_ls, vor_ls, aor_ls = backend.lines_or_info()
        big_err_lid = np.where(np.abs(por_ls - por_pp) > 5000)[0]
        backend.line_ex_to_subid[big_err_lid]
        psub_ls, qsub_ls, pbus_ls, qbus_ls, diff_v_bus_ls = backend.check_kirchoff()
        # below it does not work due to a bug fixed in dev_1.8.2 (after 1.8.2.dev4)
        # psub_pp, qsub_pp, pbus_pp, qbus_pp, diff_v_bus_pp = backend.init_pp_backend.check_kirchoff()
        
        # check voltages
        line_or_theta_pp, line_ex_theta_pp, *_ = backend.init_pp_backend.get_theta()
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
        Ybus =  backend._grid.get_dcYbus()
        Bbus = pp_net._ppc["internal"]["Bbus"]
        assert np.abs(Ybus - Bbus).max() <= self.tol
        # check the voltage angles
        Va_pp = pp_net.res_bus["va_degree"].values[:nb_sub]
        Va_ls = np.rad2deg(np.angle(backend.V[:nb_sub]))

        # check flows
        assert np.abs(Va_pp - Va_ls).max() <= self.tol
        assert max_mis <= self.tol, f"Error: por do not match, maximum absolute error is {max_mis:.5f} MW"
        max_mis = np.max(np.abs(qor_ls - qor_pp))
        assert max_mis <= self.tol, f"Error: qor do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(vor_ls - vor_pp))
        assert max_mis <= self.tol, f"Error: vor do not match, maximum absolute error is {max_mis:.5f} kV"
        max_mis = np.max(np.abs(aor_ls - aor_pp))
        assert max_mis <= 10. * self.tol, f"Error: aor do not match, maximum absolute error is {max_mis:.5f} A"
        
        load_p, load_q, load_v = backend.loads_info()
        max_mis = np.max(np.abs(load_p - load_p_pp))
        assert max_mis <= self.tol, f"Error: load_p do not match, maximum absolute error is {max_mis:.5f} MW"
        max_mis = np.max(np.abs(load_q - load_q_pp))
        assert max_mis <= self.tol, f"Error: load_q do not match, maximum absolute error is {max_mis:.5f} MVAr"
        max_mis = np.max(np.abs(load_v - load_v_pp))
        assert max_mis <= self.tol, f"Error: load_v do not match, maximum absolute error is {max_mis:.5f} kV"

        gen_p, gen_q, gen_v = backend.generators_info()
        sgen_p, sgen_q, sgen_v = backend._grid.get_sgens_res()
        # pandapower is not correct on dc...
        # max_mis = np.max(np.abs(gen_p - gen_p_pp))
        # assert max_mis <= self.tol, f"Error: gen_p do not match, maximum absolute error is {max_mis:.5f} MW"
        assert abs(np.sum(gen_p) + np.sum(sgen_p) - np.sum(load_p)) <= self.tol_big
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
        if sh_p.size:
            max_mis = np.max(np.abs(sh_p - sh_p_pp))
            assert max_mis <= self.tol, f"Error: sh_p do not match, maximum absolute error is {max_mis:.5f} MW"
        # max_mis = np.max(np.abs(sh_q - sh_q_pp))
        # assert max_mis <= self.tol, f"Error: sh_q do not match, maximum absolute error is {max_mis:.5f} MVAr"
        # again pandapower does weird stuff in dc...

        # assert np.max(np.abs(sh_q)) <= self.tol
        # max_mis = np.max(np.abs(sh_v - sh_v_pp))
        # assert max_mis <= self.tol, f"Error: sh_v do not match, maximum absolute error is {max_mis:.5f} kV"
        # again, pandapower put nan for the voltages...
