# Copyright (c) 2023, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import unittest
import warnings
import os
from lightsim2grid import LightSimBackend
import grid2op


class BackendTester(unittest.TestCase):
    """issue is still not replicated and these tests pass"""
    def setUp(self) -> None:
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.path = os.path.join(dir_path, "case_14_iidm")
        self.file_name = "grid.xiidm"

    def _aux_prep_backend(self, backend):
        backend.set_env_name("case_14_iidm")
        backend.load_grid(self.path, self.file_name)
        backend.load_storage_data(self.path)
        backend.load_redispacthing_data(self.path)
        backend.assert_grid_correct()  
    
    def _aux_get_loader_kwargs(self):
        return {"use_buses_for_sub": True, "double_bus_per_sub": True}
    
    def test_load_grid(self):
        backend = LightSimBackend(loader_method="pypowsybl",
                                  loader_kwargs=self._aux_get_loader_kwargs())
        self._aux_prep_backend(backend)
        cls = type(backend)
        assert cls.n_line == 20, f"wrong number of line {cls.n_line} vs 20"
        assert cls.n_load == 11, f"wrong number of load {cls.n_load} vs 11"
        assert cls.n_gen == 5, f"wrong number of gen {cls.n_gen} vs 5"
        assert cls.n_sub == 14, f"wrong number of gen {cls.n_sub} vs 14"
        assert cls.n_storage == 0, f"wrong number of storage { cls.n_storage} vs 0"
        assert cls.n_shunt == 1, f"wrong number of shunt { cls.n_shunt} vs 1"
        assert cls.shunts_data_available
        
        assert (backend._LightSimBackend__init_topo_vect == 1).all()
        assert backend._LightSimBackend__nb_powerline == 17
        assert backend._LightSimBackend__nb_bus_before == 14
        assert backend.nb_bus_total == 28
        
    def test_runpf(self):
        backend = LightSimBackend(loader_method="pypowsybl", loader_kwargs=self._aux_get_loader_kwargs())
        self._aux_prep_backend(backend)
        # AC powerflow
        conv, exc_ = backend.runpf()
        assert conv
        # DC powerflow
        conv, exc_ = backend.runpf(is_dc=True)
        assert conv
                
# TODO env tester
if __name__ == "__main__":
    unittest.main()
        