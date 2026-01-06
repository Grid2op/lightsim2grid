# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import pypowsybl as pp
import unittest

from lightsim2grid.gridmodel import init_from_pypowsybl


# TODO test when an element is disconnected too
class _AuxCorrectGraph:    
    """Test that the graph, when reading from pypowsybl is correct lightsim2grid side"""
    def use_buses_for_sub(self):
        return True
    
    def sort_index(self):
        return True
    
    def get_pypo_grid_name(self):
        return "ieee14"
    
    def fetch_network_ref(self):
        self.network_ref = getattr(pp.network, f"create_{self.pypo_grid_name}")()
    
    def get_n_busbar_per_sub(self):
        return None
    
    def setUp(self) -> None:
        self.pypo_grid_name = self.get_pypo_grid_name()
        self.fetch_network_ref()
        
        # init lightsim2grid model
        self.gridmodel, self.el_ids = init_from_pypowsybl(self.network_ref,
                                                          gen_slack_id=0,  # we don't really care here
                                                          sort_index=self.sort_index(),
                                                          return_sub_id=True,
                                                          buses_for_sub=self.use_buses_for_sub(),
                                                          n_busbar_per_sub=self.get_n_busbar_per_sub())
        
        # use some data
        if self.use_buses_for_sub():
            self.nb_bus_total = self.network_ref.get_buses().shape[0]
        else:
            self.nb_bus_total = self.gridmodel.get_bus_vn_kv().shape[0]
        return super().setUp()
    
    def get_sub_vn_kv(self, df, name):
        if self.use_buses_for_sub():
            return self.network_ref.get_voltage_levels().loc[df.loc[name, "voltage_level_id"], "nominal_v"]
        else:
            return df.loc[name, "nominal_v"]
        
    def test_correct_substations(self):
        if self.use_buses_for_sub():
            df = self.network_ref.get_buses()
        else:
            df = self.network_ref.get_voltage_levels()
        assert len(self.gridmodel.get_substations()) == df.shape[0]
        
        for sub in self.gridmodel.get_substations():
            nm = sub.name
            assert nm != ""  # name should be set
            assert nm in df.index  # name should be in the dataframe
            
            th_vn_kv = self.get_sub_vn_kv(df, nm)
            assert sub.vn_kv == th_vn_kv
        if not self.sort_index():
            assert (self.gridmodel.get_substation_names() == df.index).all()
        else:
            assert (self.gridmodel.get_substation_names() == df.index.sort_values()).all()
            
    def test_correct_lines(self):
        df_els = self.network_ref.get_lines()
        vls = self.network_ref.get_voltage_levels()
        
        for line in self.gridmodel.get_lines():
            assert line.name in df_els.index
            sub1 = self.gridmodel.get_substations()[line.sub1_id]
            sub2 = self.gridmodel.get_substations()[line.sub2_id]
            assert sub1.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level1_id"], "nominal_v"]
            assert sub2.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level2_id"], "nominal_v"]
            if self.use_buses_for_sub():
                if line.connected1:
                    assert sub1.name == df_els.loc[line.name, "bus1_id"]
                else:
                    pass # TODO
                if line.connected2:
                    assert sub2.name == df_els.loc[line.name, "bus2_id"]
                else:
                    pass # TODO
            else:
                assert sub1.name == df_els.loc[line.name, "voltage_level1_id"]
                assert sub2.name == df_els.loc[line.name, "voltage_level2_id"]
                
        if not self.sort_index():
            assert ([el.name for el in self.gridmodel.get_lines()] == df_els.index).all()
        else:
            assert ([el.name for el in self.gridmodel.get_lines()] == df_els.index.sort_values()).all()
            
    def test_correct_trafos(self):
        df_els = self.network_ref.get_2_windings_transformers()
        vls = self.network_ref.get_voltage_levels()
        
        for line in self.gridmodel.get_trafos():
            assert line.name in df_els.index
            sub1 = self.gridmodel.get_substations()[line.sub1_id]
            sub2 = self.gridmodel.get_substations()[line.sub2_id]
            assert sub1.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level1_id"], "nominal_v"]
            assert sub2.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level2_id"], "nominal_v"]
            if self.use_buses_for_sub():
                if line.connected1:
                    assert sub1.name == df_els.loc[line.name, "bus1_id"]
                else:
                    pass # TODO
                if line.connected2:
                    assert sub2.name == df_els.loc[line.name, "bus2_id"]
                else:
                    pass # TODO
            else:
                assert sub1.name == df_els.loc[line.name, "voltage_level1_id"]
                assert sub2.name == df_els.loc[line.name, "voltage_level2_id"]
            
        if not self.sort_index():
            assert ([el.name for el in self.gridmodel.get_trafos()] == df_els.index).all()
        else:
            assert ([el.name for el in self.gridmodel.get_trafos()] == df_els.index.sort_values()).all()
            
    def test_correct_loads(self):
        df_els = self.network_ref.get_loads()
        vls = self.network_ref.get_voltage_levels()
        
        for line in self.gridmodel.get_loads():
            assert line.name in df_els.index
            sub = self.gridmodel.get_substations()[line.sub_id]
            assert sub.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level_id"], "nominal_v"]
            if self.use_buses_for_sub():
                if line.connected:
                    assert sub.name == df_els.loc[line.name, "bus_id"]
                else:
                    pass  # TODO
            else:
                assert sub.name == df_els.loc[line.name, "voltage_level_id"]
            
        if not self.sort_index():
            assert ([el.name for el in self.gridmodel.get_loads()] == df_els.index).all()
        else:
            assert ([el.name for el in self.gridmodel.get_loads()] == df_els.index.sort_values()).all()
            
    def test_correct_gens(self):
        df_els = self.network_ref.get_generators()
        vls = self.network_ref.get_voltage_levels()
        
        for line in self.gridmodel.get_generators():
            assert line.name in df_els.index
            sub = self.gridmodel.get_substations()[line.sub_id]
            assert sub.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level_id"], "nominal_v"]
            if self.use_buses_for_sub():
                if line.connected:
                    assert sub.name == df_els.loc[line.name, "bus_id"]
                else:
                    pass  # TODO
            else:
                assert sub.name == df_els.loc[line.name, "voltage_level_id"]
            
        if not self.sort_index():
            assert ([el.name for el in self.gridmodel.get_generators()] == df_els.index).all()
        else:
            assert ([el.name for el in self.gridmodel.get_generators()] == df_els.index.sort_values()).all()
            
    def test_correct_shunts(self):
        df_els = self.network_ref.get_shunt_compensators()
        vls = self.network_ref.get_voltage_levels()
        
        for line in self.gridmodel.get_shunts():
            assert line.name in df_els.index
            sub = self.gridmodel.get_substations()[line.sub_id]
            assert sub.vn_kv == vls.loc[df_els.loc[line.name, "voltage_level_id"], "nominal_v"]
            if self.use_buses_for_sub():
                assert sub.name == df_els.loc[line.name, "bus_id"]
            else:
                assert sub.name == df_els.loc[line.name, "voltage_level_id"]
            
        if not self.sort_index():
            assert ([el.name for el in self.gridmodel.get_shunts()] == df_els.index).all()
        else:
            assert ([el.name for el in self.gridmodel.get_shunts()] == df_els.index.sort_values()).all()
            
        
class Test_CorrectGraph_TT(_AuxCorrectGraph, unittest.TestCase):
    pass


class Test_CorrectGraph_TF(_AuxCorrectGraph, unittest.TestCase):
    def sort_index(self):
        return False
    
    
class Test_CorrectGraph_FT(_AuxCorrectGraph, unittest.TestCase):
    def use_buses_for_sub(self):
        return False


class Test_CorrectGraph_FF(_AuxCorrectGraph, unittest.TestCase):
    def use_buses_for_sub(self):
        return False
    
    def sort_index(self):
        return False
        
        
class Test_CorrectGraph_TT_300(_AuxCorrectGraph, unittest.TestCase):
    def get_pypo_grid_name(self):
        return "ieee300"


class Test_CorrectGraph_TF_300(_AuxCorrectGraph, unittest.TestCase):
    def get_pypo_grid_name(self):
        return "ieee300"
    
    def sort_index(self):
        return False
    
    
class Test_CorrectGraph_FT_300(_AuxCorrectGraph, unittest.TestCase):
    def get_pypo_grid_name(self):
        return "ieee300"
    
    def use_buses_for_sub(self):
        return False


class Test_CorrectGraph_FF_300(_AuxCorrectGraph, unittest.TestCase):
    
    def get_pypo_grid_name(self):
        return "ieee300"
    
    def use_buses_for_sub(self):
        return False
    
    def sort_index(self):
        return False
    
        
class Test_CorrectGraph_TT_disco(_AuxCorrectGraph, unittest.TestCase):
    def fetch_network_ref(self):
        self.network_ref = getattr(pp.network, f"create_{self.pypo_grid_name}")()
        # disconnect a line
        df = self.network_ref.get_lines()
        self.network_ref.update_lines(id=df.index[0], connected1=False, connected2=False)


class Test_CorrectGraph_TF_disco(_AuxCorrectGraph, unittest.TestCase):
    def fetch_network_ref(self):
        self.network_ref = getattr(pp.network, f"create_{self.pypo_grid_name}")()
        # disconnect a trafo
        df = self.network_ref.get_2_windings_transformers()
        self.network_ref.update_2_windings_transformers(id=df.index[0], connected1=False, connected2=False)
    
    def sort_index(self):
        return False
    
    
class Test_CorrectGraph_FT_disco(_AuxCorrectGraph, unittest.TestCase):
    def fetch_network_ref(self):
        self.network_ref = getattr(pp.network, f"create_{self.pypo_grid_name}")()
        # disconnect a load
        df = self.network_ref.get_loads()
        self.network_ref.update_loads(id=df.index[0], connected=False)
    
    def use_buses_for_sub(self):
        return False


class Test_CorrectGraph_FF_disco(_AuxCorrectGraph, unittest.TestCase):
    def get_n_busbar_per_sub(self):
        return 3
    
    def use_buses_for_sub(self):
        return False
    
    def sort_index(self):
        return False
        
        
class Test_CorrectGraph_TT_mult_busbar(_AuxCorrectGraph, unittest.TestCase):
    def get_n_busbar_per_sub(self):
        return 3
            
    def test_line_change_bus(self):
        el_id = 0
        el = self.gridmodel.get_lines()[el_id]
        n_sub = len(self.gridmodel.get_substations())
        self.gridmodel.change_bus1_powerline(el_id, el.sub1_id + n_sub)
        self.test_correct_lines()
        
    def test_trafo_change_bus(self):
        el_id = 0
        el = self.gridmodel.get_trafos()[el_id]
        n_sub = len(self.gridmodel.get_substations())
        self.gridmodel.change_bus1_trafo(el_id, el.sub1_id + n_sub)
        self.test_correct_trafos()
        
    def test_load_change_bus(self):
        el_id = 0
        el = self.gridmodel.get_loads()[el_id]
        n_sub = len(self.gridmodel.get_substations())
        self.gridmodel.change_bus_load(el_id, el.sub_id + n_sub)
        self.test_correct_loads()
        
    def test_gen_change_bus(self):
        el_id = 0
        el = self.gridmodel.get_generators()[el_id]
        n_sub = len(self.gridmodel.get_substations())
        self.gridmodel.change_bus_gen(el_id, el.sub_id + n_sub)
        self.test_correct_gens()
        
    def test_shunt_change_bus(self):
        el_id = 0
        el = self.gridmodel.get_shunts()[el_id]
        n_sub = len(self.gridmodel.get_substations())
        self.gridmodel.change_bus_shunt(el_id, el.sub_id + n_sub)
        self.test_correct_shunts()


class Test_CorrectGraph_TF_mult_busbar(Test_CorrectGraph_TT_mult_busbar):
    def get_n_busbar_per_sub(self):
        return 3
    
    def sort_index(self):
        return False
    
    
class Test_CorrectGraph_FT_mult_busbar(Test_CorrectGraph_TT_mult_busbar):
    def get_n_busbar_per_sub(self):
        return 3
    
    def use_buses_for_sub(self):
        return False


class Test_CorrectGraph_FF_mult_busbar(Test_CorrectGraph_TT_mult_busbar):
    def get_n_busbar_per_sub(self):
        return 3
    
    def use_buses_for_sub(self):
        return False
    
    def sort_index(self):
        return False


if __name__ == "__main__":
    unittest.main()
