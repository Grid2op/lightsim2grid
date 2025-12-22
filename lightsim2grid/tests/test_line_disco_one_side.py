# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import pickle
import tempfile
import unittest
import warnings
import numpy as np
import pypowsybl
import pypowsybl.network as pp_network
import pypowsybl as pp
import pypowsybl.loadflow as pp_lf

import grid2op
from lightsim2grid import LightSimBackend

from test_GridModel_pypowsybl import BaseTests, MakeACTests

from test_match_with_pypowsybl.utils_for_slack import (
    get_pypowsybl_parameters,
    get_same_slack
)

class BaseDiscoOneSide:      
    """Test the information are correctly computed on the gridmodel side in all 4 cases"""  
    def synch_status_both_side(self):
        return True
    def ignore_status_global(self):
        return False
    
    def setUp(self):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make("l2rpn_case14_sandbox", test=True, backend=LightSimBackend())
        self.model = self.env.backend._grid
        self.model.set_ignore_status_global(self.ignore_status_global())
        self.model.set_synch_status_both_side(self.synch_status_both_side())
    
    def test_pickle(self):
        """test pickle does not lose the info"""
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "test_pickle.pickle"), "wb") as f:
                pickle.dump(self.model, f)
            with open(os.path.join(tmpdir, "test_pickle.pickle"), "rb") as f:
                backend_1 = pickle.load(f)
        assert self.model.get_ignore_status_global() == backend_1.get_ignore_status_global()
        assert self.model.get_synch_status_both_side() == backend_1.get_synch_status_both_side()
    
    def test_copy(self):
        """test copy preserve the flags"""
        backend_1 = self.model.copy()
        assert self.model.get_ignore_status_global() == backend_1.get_ignore_status_global()
        assert self.model.get_synch_status_both_side() == backend_1.get_synch_status_both_side()
        
    def test_gridmodel_line_global(self):
        """test disconnecting globally the line has the correct effect"""
        el_id = 0
        self.model.deactivate_powerline(el_id)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        else:
            assert not self.model.get_lines()[el_id].connected_global
        assert not self.model.get_lines()[el_id].connected_1
        assert not self.model.get_lines()[el_id].connected_2
        
        self.model.reactivate_powerline(el_id)
        assert self.model.get_lines()[el_id].connected_global
        assert self.model.get_lines()[el_id].connected_1
        assert self.model.get_lines()[el_id].connected_2
        
    def test_gridmodel_line_side1_topo(self):
        """test disconnecting the line side1 when updated from topology"""
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id]
        change[el_tp] = True
        new_values[el_tp] = -1
        
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert not self.model.get_lines()[el_id].connected_1
        
        if self.synch_status_both_side():
            if not self.ignore_status_global():
                assert not self.model.get_lines()[el_id].connected_global
            assert not self.model.get_lines()[el_id].connected_2
        else:
            assert self.model.get_lines()[el_id].connected_2
        
        # now reconnects it
        new_values[el_tp] = 1
        self.model.update_topo(change, new_values)
        if self.ignore_status_global():
            assert self.model.get_lines()[el_id].connected_global
        assert self.model.get_lines()[el_id].connected_1
        
        if self.synch_status_both_side():
            assert self.model.get_lines()[el_id].connected_global
            assert self.model.get_lines()[el_id].connected_2
        else:
            assert self.model.get_lines()[el_id].connected_2
            
    # def test_gridmodel_line_side1_changebus(self):
    #     """test disconnecting the line side1, when updated from change_bus"""
    #     # TODO lightsim2grid forbid this atm
    #     el_id = 0        
    #     self.model.change_bus_powerline_or(el_id, -1)
    #     if self.ignore_status_global():
    #         assert self.model.get_lines()[el_id].connected_global
    #     assert not self.model.get_lines()[el_id].connected_1
        
    #     if self.synch_status_both_side():
    #         if not self.ignore_status_global():
    #             assert not self.model.get_lines()[el_id].connected_global
    #         assert not self.model.get_lines()[el_id].connected_2
    #     else:
    #         assert self.model.get_lines()[el_id].connected_2
        
    #     # now reconnects it
    #     self.model.change_bus_powerline_or(el_id, self.model.get_lines()[el_id].sub_1_id)
    #     if self.ignore_status_global():
    #         assert self.model.get_lines()[el_id].connected_global
    #     assert self.model.get_lines()[el_id].connected_1
        
    #     if self.synch_status_both_side():
    #         assert self.model.get_lines()[el_id].connected_global
    #         assert self.model.get_lines()[el_id].connected_2
    #     else:
    #         assert self.model.get_lines()[el_id].connected_2
            
    def test_gridmodel_line_both_sides_topo(self):
        """test disconnecting both sides of the line"""
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp1 = self.env.line_or_pos_topo_vect[el_id]
        el_tp2 = self.env.line_ex_pos_topo_vect[el_id]
        change[el_tp1] = True
        new_values[el_tp1] = -1
        change[el_tp2] = True
        new_values[el_tp2] = -1
        
        self.model.update_topo(change, new_values)
        assert not self.model.get_lines()[el_id].connected_1
        assert not self.model.get_lines()[el_id].connected_2
        if self.ignore_status_global():
            # flag not modified
            assert self.model.get_lines()[el_id].connected_global
        else:
            # automatic disconnection: both sides are disconnected
            assert not self.model.get_lines()[el_id].connected_global
        
        # now reconnects it
        new_values[el_tp1] = 1
        new_values[el_tp2] = 1
        self.model.update_topo(change, new_values)
        assert self.model.get_lines()[el_id].connected_1
        assert self.model.get_lines()[el_id].connected_2
        # both sides are reconnected, so this should reconnect this automatically
        # if it was disconnected
        assert self.model.get_lines()[el_id].connected_global
        
    
class GridModelDiscoOneSideTF(BaseDiscoOneSide, unittest.TestCase):
    pass


class GridModelDiscoOneSideTT(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return True
    def ignore_status_global(self):
        return True
    
    
class GridModelDiscoOneSideFT(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return False
    def ignore_status_global(self):
        return True
    
    
class GridModelDiscoOneSideFF(BaseDiscoOneSide, unittest.TestCase):
    def synch_status_both_side(self):
        return False
    def ignore_status_global(self):
        return False
    
    
class TestPFOk(unittest.TestCase):    
    """test powerflows are correctly working when powerlines are disconnected"""
    def setUp(self, grid_nm="ieee118"):
        BaseTests.setUp(self, grid_nm)
        
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.path = os.path.join(dir_path, "case_14_iidm")
        
        self.env = grid2op.make(
            "blank", 
            test=True,
            grid_path=self.path,
            n_busbar=1,
            backend = LightSimBackend(
                loader_method="pypowsybl",
                gen_slack_id=[el.id for el in self.model.get_generators() if el.is_slack],
                loader_kwargs={
                    "grid": self.net_datamodel,
                    "use_buses_for_sub": True,
                    "sort_index": False,
                    "use_grid2op_default_names": False,
                }))
        self.model = self.env.backend._grid
        self.model.set_ignore_status_global(True)
        self.model.set_synch_status_both_side(False)
    
    def make_v0(self, net):
        V0 = np.full(self.model.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        return V0
    
    def assert_equal(self, tmp, ref, error="", tol=None):
        return BaseTests.assert_equal(self, tmp, ref, error, tol)
        
    def run_me_pf(self, V0):
        self.dc = False
        return self.model.ac_pf(V0, self.max_it, self.tol)

    def run_ref_pf(self, net):
        MakeACTests.run_ref_pf(self, net)
    
    def test_gridmodel_line_side1(self):
        """test disconnecting the line side1"""
        
        # update lightsim2grid
        el_id = 0
        change = np.zeros(self.env.dim_topo, dtype=bool)
        new_values = np.zeros(self.env.dim_topo, dtype=int)
        el_tp = self.env.line_or_pos_topo_vect[el_id]
        change[el_tp] = True
        new_values[el_tp] = -1
        self.model.update_topo(change, new_values)
        
        
        from lightsim2grid.gridmodel import init_from_pypowsybl
        model2 = init_from_pypowsybl(
            self.net_datamodel,
            slack_bus_id=self.ls_slack,
            buses_for_sub=True,
            sort_index=False,
            n_busbar_per_sub=2)
        model2.tell_solver_need_reset()
        model2.change_bus_powerline_or(el_id, model2.get_lines()[el_id].sub_1_id + len(model2.get_substations()))
        V0 = np.full(model2.total_bus(),
                     fill_value=1.0,
                     dtype=complex)
        Vres = model2.ac_pf(V0, 10, 1e-6)
            
        # update pypowsybl
        self.net_ref.update_lines(
            # id=self.net_ref.get_lines().index[el_id],
            id=self.model.get_lines()[el_id].name,
            connected1=False,
        )
        
        # compute powerflows
        Vfinal = BaseTests._run_both_pf(self, self.net_ref)
        
        # sanity check (for Vres, should match)
        model_init = self.model
        self.model = model2
        BaseTests.check_res(self, Vres[:118], self.net_ref)
        
        # compare Vres[:118] and Vfinal
        self.model = model_init
        assert np.abs(Vres[:self.model.total_bus()] - Vfinal).max() <= self.tol, f"{np.abs(Vres[:self.model.total_bus()] - Vfinal).max()}"
        BaseTests.check_res(self, Vfinal, self.net_ref)

# TODO remove the BaseTests and MakeACTests
# TODO DC powerflow
# TODO trafo with ratio
# TODO trafo with alpha (phase shift)
# TODO FDPF powerflow too
