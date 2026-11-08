# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""issue #274: the input checks of the pandapower transformer converter (`_aux_add_trafo`)."""

import unittest
import warnings
from unittest import mock

import numpy as np
import pandapower as pp

from lightsim2grid.network import init_from_pandapower
from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION


def _small_net():
    """ext_grid -- trafo (with a tap off neutral) -- load"""
    net = pp.create_empty_network()
    b_hv = pp.create_bus(net, vn_kv=110.)
    b_lv = pp.create_bus(net, vn_kv=20.)
    pp.create_ext_grid(net, b_hv)
    pp.create_transformer_from_parameters(net, b_hv, b_lv, sn_mva=40., vn_hv_kv=110., vn_lv_kv=20.,
                                          vk_percent=10., vkr_percent=0.3, pfe_kw=20., i0_percent=0.05,
                                          shift_degree=0., tap_side="hv", tap_neutral=0, tap_min=-5,
                                          tap_max=5, tap_pos=2, tap_step_percent=2.5)
    pp.create_load(net, b_lv, p_mw=20., q_mvar=5.)
    return net


def _pp_orig_file():
    return "pandapower_v2" if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER else "pandapower_v3"


def _convert(net):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        model = init_from_pandapower(net, pp_orig_file=_pp_orig_file())
    return model, [str(el.message) for el in w]


def _solve(model):
    V0 = np.full(model.total_bus(), 1.04, dtype=complex)
    V = model.ac_pf(V0, 10, 1e-8)
    assert V.shape[0] > 0, "the powerflow diverged"
    return V


class TestTrafoConverterChecks(unittest.TestCase):
    def test_leakage_ratio_warning_only_when_unsupported(self):
        """issue #274: the 'leakage_resistance_ratio_hv != 0.5 is not supported' warning
        fired when the value was 0.5 (the supported one) and never otherwise."""
        net = _small_net()
        net.trafo["leakage_resistance_ratio_hv"] = 0.5
        _, msgs = _convert(net)
        self.assertFalse(any("leakage_resistance_ratio_hv" in m for m in msgs), msgs)

        net = _small_net()
        net.trafo["leakage_resistance_ratio_hv"] = 0.3
        _, msgs = _convert(net)
        self.assertTrue(any("leakage_resistance_ratio_hv" in m for m in msgs), msgs)

    def test_nan_shift_degree_warns(self):
        """issue #274: the NaN check of shift_degree once tested tap_pos instead, so a NaN
        shift with a finite tap became 0 without a warning (already fixed on this branch;
        kept as a guard)."""
        net = _small_net()
        net.trafo["shift_degree"] = np.nan
        _, msgs = _convert(net)
        self.assertTrue(any("shift_degree" in m for m in msgs), msgs)

    def test_missing_tap_side_is_hv_and_warns(self):
        """issue #274: a missing tap_side became the low voltage side (`== "hv"` is False for
        None), and the NaN check ran on that boolean array so it never warned. It must warn
        and give the same transformer as tap_side="hv", as the warning says."""
        ref = _small_net()
        model_ref, _ = _convert(ref)
        V_ref = _solve(model_ref)

        net = _small_net()
        net.trafo["tap_side"] = None
        model, msgs = _convert(net)
        self.assertTrue(any("tap_side" in m for m in msgs), msgs)
        V = _solve(model)
        np.testing.assert_allclose(V, V_ref, rtol=0., atol=1e-10)

    def test_ideal_trafo3w_is_checked(self):
        """issue #274: the "Ideal" check of the 3-winding transformers sat in an `elif` behind
        the 2-winding one, which pandapower 3 always takes, so it could not be reached. The
        3-winding transformers are rejected earlier by `_aux_check_legit`, which is bypassed
        here to reach the transformer converter itself."""
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            self.skipTest("tap_changer_type only exists with pandapower 3")
        net = _small_net()
        b_mv = pp.create_bus(net, vn_kv=20.)
        b_lv = pp.create_bus(net, vn_kv=10.)
        pp.create_transformer3w_from_parameters(net, 0, b_mv, b_lv, vn_hv_kv=110., vn_mv_kv=20., vn_lv_kv=10.,
                                                sn_hv_mva=40., sn_mv_mva=20., sn_lv_mva=20.,
                                                vk_hv_percent=10., vk_mv_percent=10., vk_lv_percent=10.,
                                                vkr_hv_percent=0.3, vkr_mv_percent=0.3, vkr_lv_percent=0.3,
                                                pfe_kw=20., i0_percent=0.05)
        net.trafo3w["tap_changer_type"] = "Ideal"
        with mock.patch("lightsim2grid.network.from_pandapower.initLSGrid._aux_check_legit"):
            with self.assertRaisesRegex(RuntimeError, "3-winding"):
                _convert(net)


if __name__ == "__main__":
    unittest.main()
