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

try:
    from grid2op._create_test_suite import create_test_suite
    CAN_DO_TEST_SUITE = True
except ImportError as exc_:
    CAN_DO_TEST_SUITE = False


def _aux_get_loader_kwargs():
    return {"use_buses_for_sub": True, "double_bus_per_sub": True, "gen_slack_id": 6}
    
    
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
    
    def test_load_grid(self):
        backend = LightSimBackend(loader_method="pypowsybl",
                                  loader_kwargs=_aux_get_loader_kwargs())
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
        backend = LightSimBackend(loader_method="pypowsybl", loader_kwargs=_aux_get_loader_kwargs())
        self._aux_prep_backend(backend)
        # AC powerflow
        conv, exc_ = backend.runpf()
        assert conv
        # DC powerflow
        conv, exc_ = backend.runpf(is_dc=True)
        assert conv

if CAN_DO_TEST_SUITE:
    dir_path = os.path.dirname(os.path.realpath(__file__))
    path_case_14_storage_iidm = os.path.join(dir_path, "case_14_storage_iidm")
    from grid2op.tests.aaa_test_backend_interface import AAATestBackendAPI
    
    class TestBackendAPI_PyPoBk(AAATestBackendAPI, unittest.TestCase):
        def get_path(self):
            return path_case_14_storage_iidm
        
        def get_casefile(self):
            return "grid.xiidm"
        
        def make_backend(self, detailed_infos_for_cascading_failures=False):
            return  LightSimBackend(loader_method="pypowsybl",
                                    loader_kwargs=_aux_get_loader_kwargs(),
                                    detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures)
                            
    # # add test of grid2op for the backend based on pypowsybl
    # def this_make_backend(self, detailed_infos_for_cascading_failures=False):
    #     return LightSimBackend(
    #             loader_method="pypowsybl",
    #             loader_kwargs=_aux_get_loader_kwargs(),
    #             detailed_infos_for_cascading_failures=detailed_infos_for_cascading_failures
    #         )
    # add_name_cls = "test_LightSimBackend_pypowsybl"

    # def get_path_test_api(self):
    #     return path
    
    # def get_casefile(self):
    #     return "grid.xiidm"
    
    # res = create_test_suite(make_backend_fun=this_make_backend,
    #                         add_name_cls=add_name_cls,
    #                         add_to_module=__name__,
    #                         extended_test=False,  # for now keep `extended_test=False` until all problems are solved
    #                         get_paths={"AAATestBackendAPI": get_path_test_api},
    #                         get_casefiles={"AAATestBackendAPI": get_casefile}
    #                         )         
# TODO env tester
if __name__ == "__main__":
    unittest.main()
        