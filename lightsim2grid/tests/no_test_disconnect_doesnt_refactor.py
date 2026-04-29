# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import unittest
import warnings
import numpy as np
import grid2op
from grid2op.Action import CompleteAction

from lightsim2grid import LightSimBackend
from lightsim2grid.algorithm import AlgorithmType


class TestRefactorNotUsed(unittest.TestCase):
    def _aux_setup_grid(self):
        self.need_dc = True  # is it worth it to run DC powerflow ?
        self.can_dist_slack = True
        self.gridmodel.change_algorithm(AlgorithmType.NR_SparseLU)
        self.gridmodel.change_algorithm(AlgorithmType.DC_SparseLU)
        
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("educ_case14_storage",
                                    test=True,
                                    action_class=CompleteAction,
                                    backend=LightSimBackend())
        self.gridmodel = self.env.backend._grid
        self.iter = 10
        self._aux_setup_grid()
        self.v_init = 0.0 * self.env.backend.V + 1.04  # just to have a vector with the right dimension
        self.tol_solver = 1e-8  # solver
        self.tol_equal = 1e-10  #  for comparing with and without the "smarter solver" things, and make sure everything is really equal!
            
    def test_disco_line_doesnt_refactor(self):
        # perform init powerflow
        self.gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
        
        el_id = 0
        self.gridmodel.deactivate_powerline(el_id)
        print("now now now")
        self.gridmodel.ac_pf(self.v_init, self.iter, self.tol_solver)
        print("PF done")
        t_fx, t_solve, t_refactor, t_init, t_check, t_sbus, t_fillJ, t_va, t_preproc_, t_total = self.gridmodel.get_solver().get_timers_jacobian()
        print(f"{t_refactor = }")
        print(f"{t_init = }")