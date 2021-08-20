# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import unittest
import copy
import numpy as np
from scipy import sparse
from lightsim2grid.initGridModel import init
import pandapower.networks as pn
import pandapower as pp
import warnings
from test_GridModel import BaseTests
import pdb


class MakeACTestsDisco(BaseTests, unittest.TestCase):
    def setUp(self):
        self.net = pn.case118()
        self.last_real_bus = self.net.bus.shape[0]
        pp.create_bus(self.net, vn_kv=self.net.bus["vn_kv"][0])
        self.net.bus["in_service"][self.last_real_bus] = False
        self.net_ref = copy.deepcopy(self.net)
        self.net_datamodel = copy.deepcopy(self.net)
        self.n_bus = self.net.bus.shape[0]

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.model = init(self.net)

        self.model.deactivate_bus(self.last_real_bus)

        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 1e-5  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)

    def run_me_pf(self, V0):
        return self.model.ac_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.runpp(net, init="flat")

    def do_i_skip(self, test_nm):
        return


class MakeDCTestsDisco(BaseTests, unittest.TestCase):
    def setUp(self):
        self.net = pn.case118()
        self.last_real_bus = self.net.bus.shape[0]
        pp.create_bus(self.net, vn_kv=self.net.bus["vn_kv"][0])
        self.net.bus["in_service"][self.last_real_bus] = False
        self.net_ref = copy.deepcopy(self.net)
        self.net_datamodel = copy.deepcopy(self.net)
        self.n_bus = self.net.bus.shape[0]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.model = init(self.net)
        self.model.deactivate_bus(self.last_real_bus)

        self.max_it = 10
        self.tol = 1e-8  # tolerance for the solver
        self.tol_test = 1e-5  # tolerance for the test (2 matrices are equal if the l_1 of their difference is less than this)

    def run_me_pf(self, V0):
        return self.model.dc_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pp.rundcpp(net, init="flat")

    def do_i_skip(self, test_nm):
        pass

    def check_res(self, Vfinal, net):
        assert Vfinal.shape[0] > 0, "powerflow diverged !"
        tmp_bus_ind = np.argsort(net.bus.index)
        tmp_bus_ind = tmp_bus_ind[tmp_bus_ind != self.last_real_bus]
        va_deg = net.res_bus["va_degree"].values
        self.assert_equal(np.angle(Vfinal)[:self.last_real_bus], va_deg[tmp_bus_ind] / 180. * np.pi)


if __name__ == "__main__":
    unittest.main()
