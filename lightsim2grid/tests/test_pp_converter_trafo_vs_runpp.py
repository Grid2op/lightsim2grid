# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""issue #274: the pandapower transformer converter against ``pp.runpp``, one transformer
variant at a time, on the smallest grid that has one (ext_grid -- trafo -- line -- load)."""

import copy
import unittest
import warnings

import numpy as np
import pandapower as pp

from lightsim2grid.network import init_from_pandapower
from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION

# both solvers converge far below it: what is left is the model
TOL_V = 1e-8


def make_net(**trafo_kw):
    net = pp.create_empty_network()
    b_hv = pp.create_bus(net, vn_kv=110.)
    b_lv = pp.create_bus(net, vn_kv=20.)
    b_ld = pp.create_bus(net, vn_kv=20.)
    pp.create_ext_grid(net, b_hv, vm_pu=1.02)
    par = dict(sn_mva=40., vn_hv_kv=110., vn_lv_kv=20., vk_percent=10., vkr_percent=0.3, pfe_kw=30.,
               i0_percent=0.1, shift_degree=0., tap_side="hv", tap_neutral=0, tap_min=-9, tap_max=9,
               tap_pos=3, tap_step_percent=1.5, tap_step_degree=0., tap_changer_type="Ratio")
    par.update(trafo_kw)
    pp.create_transformer_from_parameters(net, b_hv, b_lv, **par)
    pp.create_line_from_parameters(net, b_lv, b_ld, 2., 0.1, 0.3, 10., 0.4)
    pp.create_load(net, b_ld, p_mw=25., q_mvar=8.)
    return net


def pp_voltages(net, trafo_model="t"):
    net = copy.deepcopy(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp(net, calculate_voltage_angles=True, init="flat", trafo_model=trafo_model,
                 tolerance_mva=1e-11, lightsim2grid=False, numba=False)
    return net.res_bus["vm_pu"].to_numpy() * np.exp(1j * np.deg2rad(net.res_bus["va_degree"].to_numpy()))


def ls_model(net, trafo_model="t"):
    net = copy.deepcopy(net)
    net._options = dict(trafo_model=trafo_model)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pandapower(net, pp_orig_file="pandapower_v3")


def ls_voltages(model, n_bus):
    V = model.ac_pf(np.full(model.total_bus(), 1. + 0j), 20, 1e-11)
    assert V.shape[0] > 0, "the lightsim2grid powerflow diverged"
    return V[:n_bus]


class _Base(unittest.TestCase):
    def setUp(self):
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            self.skipTest("compares with the pandapower 3 transformer model")

    def assert_same_as_pp(self, net, trafo_model="t"):
        V_pp = pp_voltages(net, trafo_model)
        V_ls = ls_voltages(ls_model(net, trafo_model), V_pp.shape[0])
        np.testing.assert_allclose(V_ls, V_pp, rtol=0., atol=TOL_V)


class TestTrafoVsRunpp(_Base):
    def test_ratio_tap(self):
        """the plain case: a ratio tap changer off neutral, with iron losses and a magnetising
        current (the magnetising susceptance used pfe instead of pfe**2, issue #274)"""
        self.assert_same_as_pp(make_net())

    def test_ratio_tap_lv(self):
        """issue #274: with the tap on the low voltage side, the magnetising admittance was
        multiplied by (vn_trafo_lv / vn_lv)**2 where pandapower divides by it"""
        self.assert_same_as_pp(make_net(tap_side="lv"))

    def test_ratio_tap_neutral(self):
        self.assert_same_as_pp(make_net(tap_neutral=2, tap_pos=5))

    def test_shift_degree(self):
        self.assert_same_as_pp(make_net(shift_degree=30.))
        self.assert_same_as_pp(make_net(shift_degree=30., tap_side="lv"))

    def test_asymmetrical_tap(self):
        """issue #274: a tap with both a step in % and a step angle moves the voltage along a
        complex step; pandapower puts the angle of the result in the phase shift, which was
        lost here (the angle only entered the impedance base)"""
        self.assert_same_as_pp(make_net(tap_step_degree=30.))
        self.assert_same_as_pp(make_net(tap_step_degree=30., tap_side="lv"))
        self.assert_same_as_pp(make_net(tap_step_degree=30., shift_degree=30., tap_neutral=1, tap_pos=-4))

    def test_off_nominal_rated_voltages(self):
        """issue #274: a transformer rated for other voltages than its buses' (vn_hv_kv / vn_lv_kv
        != the buses' vn_kv) has an off-nominal ratio, and pandapower moves the tap from the rated
        voltage; the bus voltages were used for both here, which dropped that ratio"""
        self.assert_same_as_pp(make_net(vn_hv_kv=115.))
        self.assert_same_as_pp(make_net(vn_lv_kv=21.))
        self.assert_same_as_pp(make_net(vn_lv_kv=21., tap_side="lv"))
        self.assert_same_as_pp(make_net(vn_hv_kv=107., vn_lv_kv=21., tap_step_degree=30.))

    def test_tap_ignored_by_pandapower(self):
        """issue #274: pandapower 3 applies a tap only for a "Ratio", "Symmetrical" (or "Ideal")
        tap_changer_type and a tap_side of "hv" or "lv"; without them the transformer is at
        ratio 1 whatever tap_pos says, where this converter applied the tap anyway"""
        self.assert_same_as_pp(make_net(tap_changer_type=None))
        self.assert_same_as_pp(make_net(tap_side=None))

    def test_symmetrical_tap(self):
        """issue #274: pandapower 3 computes a "Symmetrical" tap changer as a "Ratio" one"""
        self.assert_same_as_pp(make_net(tap_changer_type="Symmetrical"))
        self.assert_same_as_pp(make_net(tap_changer_type="Symmetrical", tap_step_degree=30.))


if __name__ == "__main__":
    unittest.main()
