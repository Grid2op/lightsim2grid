# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings

import numpy as np
import pandapower as pp
from lightsim2grid.newtonpf import newtonpf_new
import pdb


class TestIssueJacobian(unittest.TestCase):
    """
    see conversation at https://github.com/e2nIEE/pandapower/pull/1455
    for a description of the bug
    """
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            self.net = pp.from_json("./dist_slack_test.json")
            pp.runpp(self.net, lightsim2grid=False, numba=True, distributed_slack=True)
    
    def test_jacobian_dist_slack(self):    
        """that's the grid described in the issue 
        https://github.com/e2nIEE/pandapower/pull/1455
        """
        options = {"max_iteration": 10, "tolerance_mva": 1e-8, "distributed_slack": True}
        V, converged, iterations, J, Vm_it, Va_it = newtonpf_new(self.net._ppc["internal"]["Ybus"], 
                                                                 self.net._ppc["internal"]["Sbus"],
                                                                 np.array([1.+0.j, 1.+0.j, 1.+0.j]),
                                                                 self.net._ppc["internal"]["ref"],
                                                                 self.net._ppc["internal"]["pv"],
                                                                 self.net._ppc["internal"]["pq"],
                                                                 self.net._ppc["internal"],
                                                                 options)
        assert converged
        assert iterations == 4 
        assert np.max(np.abs(V - self.net._ppc["internal"]["V"])) <= 1e-6
        if self.net._ppc["internal"]["J"].shape == (4,4):
            # with earlier pandapower version, distributed slack are
            # not taken into account, so this test does not work
            # because pandapower Jacobian has not the same shape.
            assert np.max(np.abs(J - self.net._ppc["internal"]["J"])) <= 1e-6
    
    def test_jacobian_single_slack(self):    
        """that's the grid described in the issue 
        https://github.com/e2nIEE/pandapower/pull/1455
        """
        pp.runpp(self.net, lightsim2grid=False, numba=True, distributed_slack=False)
        options = {"max_iteration": 10, "tolerance_mva": 1e-8, "distributed_slack": False}
        V, converged, iterations, J, Vm_it, Va_it = newtonpf_new(self.net._ppc["internal"]["Ybus"], 
                                                                 self.net._ppc["internal"]["Sbus"],
                                                                 np.array([1.+0.j, 1.+0.j, 1.+0.j]),
                                                                 self.net._ppc["internal"]["ref"],
                                                                 self.net._ppc["internal"]["pv"],
                                                                 self.net._ppc["internal"]["pq"],
                                                                 self.net._ppc["internal"],
                                                                 options)
        assert converged
        assert iterations == 4 
        assert np.max(np.abs(V - self.net._ppc["internal"]["V"])) <= 1e-6
        if self.net._ppc["internal"]["J"].shape == (3,3):
            # with earlier pandapower version, distributed slack are
            # not taken into account, so this test does not work
            # because pandapower Jacobian has not the same shape.
            assert np.max(np.abs(J - self.net._ppc["internal"]["J"])) <= 1e-6


if __name__ == "__main__":
    unittest.main()
