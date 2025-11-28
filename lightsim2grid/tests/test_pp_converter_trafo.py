# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy
import warnings
import numpy as np
from lightsim2grid_cpp import PandaPowerConverter
import pandapower as pp
import pandapower.networks as pn
from pandapower.build_branch import (
    _calc_branch_values_from_trafo_df,
    _wye_delta
    )
import unittest
from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


class _TestTrafoConverter:
    def get_pp_case():
        return pn.case118()
    
    def setUp(self):
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            self.skipTest(f"This tests only pandapower3 data reader, which is not usable with the current pandapower version used {CURRENT_PP_VERSION}.")
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.pp_net = pn.case118()
        self.converter = PandaPowerConverter()
        self.converter.set_sn_mva(self.pp_net.sn_mva)
        self.converter.set_f_hz(self.pp_net.f_hz)
    
    def get_tap(self, pp_net):
        # fix the missing values
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

        if "tap_phase_shifter" in pp_net.trafo:
            if np.any(pp_net.trafo["tap_phase_shifter"].values):
                raise RuntimeError("Ideal phase shifters are not modeled. Please remove all trafos with "
                                "pp_net.trafo[\"tap_phase_shifter\"] set to True.")
        elif "tap_changer_type" in pp_net.trafo:
            if np.any(pp_net.trafo["tap_changer_type"].values == "Ideal"):
                raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 2-winding trafos "
                                "with \"tap_changer_type\" set to \"Ideal\".")
        elif "tap_changer_type" in pp_net.trafo3w:
            if np.any(pp_net.trafo3w["tap_changer_type"].values == "Ideal"):
                raise RuntimeError("Ideal phase shifters are not modeled. Please remove all 3-winding trafos "
                                "with \"tap_changer_type\" set to \"Ideal\".")

        tap_angles_ = 1.0 * pp_net.trafo["tap_step_degree"].values
        tap_angles_[~np.isfinite(tap_angles_)] = 0.
        tap_angles_ = np.deg2rad(tap_angles_)
        return tap_step_pct, tap_pos, tap_angles_, is_tap_hv_side
    
    def test_model_in_pi(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.pp_net, lightsim2grid=False, numba=False, trafo_model="pi")
            garbage_pp = copy.deepcopy(self.pp_net)
            r_pu, x_pu, g_pu, b_pu, g_asym, b_asym, ratio, shift_xx = _calc_branch_values_from_trafo_df(
                    garbage_pp, garbage_pp._ppc)
        
        tap_step_pct, tap_pos, tap_angles_, is_tap_hv_side = self.get_tap(self.pp_net)
        trafo_r, trafo_x, trafo_b = \
            self.converter.get_trafo_param_pp3(tap_step_pct,
                                               tap_pos,
                                               tap_angles_,  # in radian !
                                               is_tap_hv_side,
                                               self.pp_net.bus.loc[self.pp_net.trafo["hv_bus"]]["vn_kv"],
                                               self.pp_net.bus.loc[self.pp_net.trafo["lv_bus"]]["vn_kv"],
                                               self.pp_net.trafo["vk_percent"].values,
                                               self.pp_net.trafo["vkr_percent"].values,
                                               self.pp_net.trafo["sn_mva"].values,
                                               self.pp_net.trafo["pfe_kw"].values,
                                               self.pp_net.trafo["i0_percent"].values,
                                               False  # trafo_model_is_t
                                               )
        assert np.allclose(trafo_r, r_pu)
        assert np.allclose(trafo_x, x_pu)
        assert np.allclose(trafo_b.real, g_pu)
        assert np.allclose(trafo_b.imag, b_pu)
        
    def test_model_in_t(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(self.pp_net, lightsim2grid=False, numba=False, trafo_model="t")
            garbage_pp = copy.deepcopy(self.pp_net)
            r_pu, x_pu, g_pu, b_pu, g_asym, b_asym, ratio, shift_xx = _calc_branch_values_from_trafo_df(
                    garbage_pp, garbage_pp._ppc)
        
        tap_step_pct, tap_pos, tap_angles_, is_tap_hv_side = self.get_tap(self.pp_net)
        
        ###############################
        #### pi data from lightsim2grid
        #### wye delta from pandapower
        trafo_r_pi, trafo_x_pi, trafo_b_pi = \
            self.converter.get_trafo_param_pp3(tap_step_pct,
                                               tap_pos,
                                               tap_angles_,  # in radian !
                                               is_tap_hv_side,
                                               self.pp_net.bus.loc[self.pp_net.trafo["hv_bus"]]["vn_kv"],
                                               self.pp_net.bus.loc[self.pp_net.trafo["lv_bus"]]["vn_kv"],
                                               self.pp_net.trafo["vk_percent"].values,
                                               self.pp_net.trafo["vkr_percent"].values,
                                               self.pp_net.trafo["sn_mva"].values,
                                               self.pp_net.trafo["pfe_kw"].values,
                                               self.pp_net.trafo["i0_percent"].values,
                                               False  # trafo_model_is_t
                                               )
        r_ratio = np.ones(trafo_r_pi.shape[0]) * 0.5
        x_ratio = np.ones(trafo_r_pi.shape[0]) * 0.5
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            r_pi_pp, x_pi_pp, g_pi_pp, b_pi_pp, g_asym, b_asym = _wye_delta(
                trafo_r_pi.copy(),
                trafo_x_pi.copy(),
                trafo_b_pi.real.copy(),
                trafo_b_pi.imag.copy(),
                r_ratio,
                x_ratio)
        assert np.allclose(r_pi_pp, r_pu)
        assert np.allclose(x_pi_pp, x_pu)
        assert np.allclose(g_pi_pp, g_pu)
        assert np.allclose(b_pi_pp, b_pu)
        #####################################
        
        trafo_r, trafo_x, trafo_b = \
            self.converter.get_trafo_param_pp3(tap_step_pct,
                                               tap_pos,
                                               tap_angles_,  # in radian !
                                               is_tap_hv_side,
                                               self.pp_net.bus.loc[self.pp_net.trafo["hv_bus"]]["vn_kv"],
                                               self.pp_net.bus.loc[self.pp_net.trafo["lv_bus"]]["vn_kv"],
                                               self.pp_net.trafo["vk_percent"].values,
                                               self.pp_net.trafo["vkr_percent"].values,
                                               self.pp_net.trafo["sn_mva"].values,
                                               self.pp_net.trafo["pfe_kw"].values,
                                               self.pp_net.trafo["i0_percent"].values,
                                               True  # trafo_model_is_t
                                               )
        assert np.allclose(trafo_r, r_pu)
        assert np.allclose(trafo_x, x_pu)
        assert np.allclose(trafo_b.real, g_pu)
        assert np.allclose(trafo_b.imag, b_pu)
        
        
class TestTrafoConverter_Case9(_TestTrafoConverter, unittest.TestCase):
    def get_pp_case():
        return pn.case9()
    
class TestTrafoConverter_Case14(_TestTrafoConverter, unittest.TestCase):
    def get_pp_case():
        return pn.case14()
    
class TestTrafoConverter_Case57(_TestTrafoConverter, unittest.TestCase):
    def get_pp_case():
        return pn.case57()
    
class TestTrafoConverter_Case118(_TestTrafoConverter, unittest.TestCase):
    def get_pp_case():
        return pn.case118()
    
class TestTrafoConverter_Case300(_TestTrafoConverter, unittest.TestCase):
    def get_pp_case():
        return pn.case300()
    
if __name__ == "__main__":
    unittest.main()