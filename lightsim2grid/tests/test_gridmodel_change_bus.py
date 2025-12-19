# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import numpy as np

import grid2op
from grid2op.Action import CompleteAction

from lightsim2grid import LightSimBackend

import warnings
import pdb


class BaseTests(unittest.TestCase):
    """
    Make sure the gridmodel can correctly handle changing the bus
    status of the elements (loads, generators, lines, ...) when
    input as "local bus".
    
    """
    def n_busbar_per_sub(self):
        return 2
    
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(
                "educ_case14_storage",
                test=True,
                n_busbar=self.n_busbar_per_sub(),
                _add_to_name=f"{type(self).__name__}",
                backend=LightSimBackend(),
                action_class=CompleteAction
            )
        self.gridmodel = self.env.backend._grid
        
    def tearDown(self):
        self.env.close()
        return super().tearDown()
    
    # TODO change also twosides
    
    def _aux_test_oneside_works(
        self,
        el_to_subid,
        keys="loads_id",
        fun="get_loads",
        ):
        n_sub = type(self.env).n_sub
        for el_id in range(len(getattr(self.gridmodel, fun)())):
            sub_id = el_to_subid[el_id]
            for bus_id in range(1, self.n_busbar_per_sub() + 1):
                bk_act = type(self.env.backend).my_bk_act_class()
                bk_act += self.env.action_space({"set_bus": {keys: [(el_id, bus_id)]}})
                th_bus_id = (bus_id - 1) * n_sub + sub_id
                self.gridmodel.update_topo(
                    bk_act.current_topo.changed,
                    bk_act.current_topo.values
                    )
                el_gridmodel = getattr(self.gridmodel, fun)()[el_id]
                assert el_gridmodel.bus_id == th_bus_id, f"error for {bus_id} : {el_gridmodel.bus_id} vs {th_bus_id}"
                
    def test_change_bus_load_works(self):
        self._aux_test_oneside_works(
            type(self.env).load_to_subid,
            "loads_id",
            "get_loads")
        
    def test_change_bus_gen_works(self):
        self._aux_test_oneside_works(
            type(self.env).gen_to_subid,
            "generators_id",
            "get_generators")
        
    def test_change_bus_storage_works(self):
        self._aux_test_oneside_works(
            type(self.env).storage_to_subid,
            "storages_id",
            "get_storages")
    
    # TODO change also twosides !
    def _aux_test_change_bus_oneside_fails(
        self,
        el_pos_topo_vect,
        fun="get_loads"
        ):
        for el_id in range(len(getattr(self.gridmodel, fun)())):
            pos_tpv = el_pos_topo_vect[el_id]
            changed = np.zeros(type(self.env).dim_topo, dtype=bool)
            values = np.zeros(type(self.env).dim_topo, dtype=int)
            
            changed[pos_tpv] = True
            values[pos_tpv] = -2
            with self.assertRaises(IndexError):
                self.gridmodel.update_topo(changed, values)
            values[pos_tpv] = self.n_busbar_per_sub() + 1
            with self.assertRaises(IndexError):
                self.gridmodel.update_topo(changed, values)
        
    def test_change_bus_load_fails(self):
        self._aux_test_change_bus_oneside_fails(
            type(self.env).load_pos_topo_vect,
            "get_loads")
        
    def test_change_bus_gen_fails(self):
        self._aux_test_change_bus_oneside_fails(
            type(self.env).gen_pos_topo_vect,
            "get_generators")
        
    def test_change_bus_storage_fails(self):
        self._aux_test_change_bus_oneside_fails(
            type(self.env).storage_pos_topo_vect,
            "get_storages")
        
    
if __name__ == "__main__":
    unittest.main()
