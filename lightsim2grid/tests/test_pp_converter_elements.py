# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The pandapower converter against ``pp.runpp``, for the elements other than the transformers
(those are in test_pp_converter_trafo_vs_runpp.py): one variant at a time, on a small meshed
grid. See docs/from_pandapower.rst for what is converted and how."""

import copy
import unittest
import warnings

import numpy as np
import pandapower as pp

from lightsim2grid.network import init_from_pandapower
from global_var_tests import MAX_PP2_DATAREADER, CURRENT_PP_VERSION

# both solvers converge far below it: what is left is the model
TOL_V = 1e-8


def make_net():
    """ext_grid, a generator, a load, a static generator and a storage unit on a 4-bus loop"""
    net = pp.create_empty_network()
    b = [pp.create_bus(net, vn_kv=110.) for _ in range(4)]
    pp.create_ext_grid(net, b[0], vm_pu=1.02)
    for i, j in ((0, 1), (1, 2), (2, 3), (0, 3)):
        pp.create_line_from_parameters(net, b[i], b[j], 10., 0.1, 0.4, 10., 0.6)
    pp.create_load(net, b[2], p_mw=40., q_mvar=10.)
    pp.create_gen(net, b[1], p_mw=20., vm_pu=1.01)
    pp.create_sgen(net, b[3], p_mw=5., q_mvar=1., min_q_mvar=-10., max_q_mvar=10.)
    pp.create_storage(net, b[3], p_mw=3., q_mvar=0.5, max_e_mwh=10.)
    return net


def pp_voltages(net):
    net = copy.deepcopy(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp.runpp(net, init="flat", tolerance_mva=1e-11, lightsim2grid=False, numba=False)
    return net.res_bus["vm_pu"].to_numpy() * np.exp(1j * np.deg2rad(net.res_bus["va_degree"].to_numpy()))


def ls_model(net):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return init_from_pandapower(copy.deepcopy(net), pp_orig_file="pandapower_v3")


def ls_voltages(model, n_bus):
    V = model.ac_pf(np.full(model.total_bus(), 1. + 0j), 30, 1e-11)
    assert V.shape[0] > 0, "the lightsim2grid powerflow diverged"
    return V[:n_bus]


class TestElementsVsRunpp(unittest.TestCase):
    def setUp(self):
        if CURRENT_PP_VERSION <= MAX_PP2_DATAREADER:
            self.skipTest("compares with pandapower 3")

    def assert_same_as_pp(self, net):
        V_pp = pp_voltages(net)
        np.testing.assert_allclose(ls_voltages(ls_model(net), V_pp.shape[0]), V_pp, rtol=0., atol=TOL_V)

    def test_reference(self):
        self.assert_same_as_pp(make_net())

    def test_load_and_sgen_scaling(self):
        net = make_net()
        net.load["scaling"] = 0.5
        net.sgen["scaling"] = 0.5
        self.assert_same_as_pp(net)

    def test_gen_scaling(self):
        """pandapower multiplies a generator's p_mw by its `scaling`; the converter did not read it"""
        net = make_net()
        net.gen["scaling"] = 0.5
        self.assert_same_as_pp(net)

    def test_storage_scaling(self):
        """pandapower multiplies a storage unit's p_mw and q_mvar by its `scaling`; the converter
        did not read it"""
        net = make_net()
        net.storage["scaling"] = 0.5
        self.assert_same_as_pp(net)

    def test_ext_grid_out_of_service(self):
        """an ext_grid out of service is not a slack for pandapower; the converter made it one"""
        net = make_net()
        pp.create_ext_grid(net, 2, vm_pu=1.0, in_service=False)
        self.assert_same_as_pp(net)

    def test_no_slack_refused(self):
        """with its only ext_grid out of service and no slack generator, the grid has no slack"""
        net = make_net()
        net.ext_grid["in_service"] = False
        with self.assertRaisesRegex(RuntimeError, "No slack"):
            ls_model(net)




if __name__ == "__main__":
    unittest.main()
