# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import warnings
import pandapower as pp
import numpy as np
import pandapower.networks as pn
import grid2op
import pdb
import copy
import tempfile
import unittest
from lightsim2grid import LightSimBackend
from lightsim2grid.initGridModel import init
from grid2op.Chronics import GridStateFromFileWithForecastsWithoutMaintenance as GridStateFromFile

class TestMultipleL2RPN(unittest.TestCase):
    def setUp(self) -> None:
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 118

    def test_neurips_track2(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            # "l2rpn_neurips_2020_track2", test=True  ## TODO SLACK: make a failing test with test=True
            env = grid2op.make("l2rpn_neurips_2020_track2_small",
                               data_feeding_kwargs={"gridvalueClass": GridStateFromFile})  

        li_envs = list(env.keys())

        for el in li_envs:
            self.pp_net = env[el].backend._grid
            self.idx_slack = np.where(self.pp_net.gen["slack"].values)[0]
            self.pp_net.gen["slack_weight"][self.idx_slack] = 1.
            pp.rundcpp(self.pp_net)  # to forget the result
            pp.runpp(self.pp_net, distributed_slack=True, init_vm_pu="flat", init_va_degree="flat")
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                ls_grid = init(self.pp_net)
            ls_grid.tell_topo_changed()
            V = np.ones(2 * self.nb_bus_total, dtype=np.complex_)
            V = ls_grid.ac_pf(V, self.max_it, self.tol)
            self.check_results(V[:self.nb_bus_total], ls_grid, self.pp_net)

    def check_results(self, V_ls, ls_grid, pp_net):
        assert len(V_ls), "lightsim diverged !"
        my_ref = np.where(np.angle(V_ls) == 0.)[0][0]
        V_pp = pp_net.res_bus["vm_pu"].values * np.exp(1j*np.pi / 180. *  pp_net.res_bus["va_degree"].values)
        V_pp *= np.exp(-1j * np.angle(V_pp)[my_ref])
        V_pp = V_pp[:self.nb_bus_total]
        assert np.abs(V_pp - V_ls).max() <= 1e-6, "wrong voltages"
        assert np.all(np.abs([el.res_p_or_mw for el in ls_grid.get_lines()] - pp_net.res_line["p_from_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_a_or_ka for el in ls_grid.get_lines()] - pp_net.res_line["i_from_ka"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_hv_mw for el in ls_grid.get_trafos()] - pp_net.res_trafo["p_hv_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_p_mw for el in ls_grid.get_generators()] - pp_net.res_gen["p_mw"].values) <= 1e-6)
        assert np.all(np.abs([el.res_q_mvar for el in ls_grid.get_generators()] - pp_net.res_gen["q_mvar"].values) <= 1e-6)
        

class Test118LightsimBackend(unittest.TestCase):
    def setUp(self) -> None:
        self.max_it = 10
        self.tol = 1e-8
        self.nb_bus_total = 118

    def test_make_and_pf(self):
        env_name_input = "l2rpn_neurips_2020_track2_small"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env_ls = grid2op.make(env_name_input, backend=LightSimBackend(),
                                       data_feeding_kwargs={"gridvalueClass": GridStateFromFile})


if __name__ == "__main__":
    unittest.main()