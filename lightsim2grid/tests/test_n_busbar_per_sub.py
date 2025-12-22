# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.
import unittest
import warnings
import numpy as np

import grid2op
from grid2op.Action import CompleteAction

from lightsim2grid import LightSimBackend
 
 
class TestLightSimBackend_3busbars(unittest.TestCase):
    def get_nb_bus(self):
        return 3
    
    def get_env_nm(self):
        return "educ_case14_storage"
    
    def get_backend_kwargs(self):
        return dict()
    
    def setUp(self) -> None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = grid2op.make(self.get_env_nm(),
                                    backend=LightSimBackend(**self.get_backend_kwargs()),
                                    action_class=CompleteAction,
                                    test=True,
                                    n_busbar=self.get_nb_bus(),
                                    _add_to_name=type(self).__name__ + f'_{self.get_nb_bus()}')
        self.list_loc_bus = [-1] + list(range(1, type(self.env).n_busbar_per_sub + 1))
        return super().setUp()
    
    def tearDown(self) -> None:
        self.env.close()
        return super().tearDown()
    
    def test_right_bus_made(self):
        assert len(self.env.backend._grid.get_bus_vn_kv()) == self.get_nb_bus() * type(self.env).n_sub
        assert (~np.array(self.env.backend._grid.get_bus_status())[type(self.env).n_sub:]).all()
        assert (np.array(self.env.backend._grid.get_bus_status())[:type(self.env).n_sub]).all()

    @staticmethod
    def _aux_find_sub(env, obj_col):
        """find a sub with 4 elements, the type of elements and at least 2 lines"""
        cls = type(env)
        res = None
        for sub_id in range(cls.n_sub):
            this_sub_mask = cls.grid_objects_types[:,cls.SUB_COL] == sub_id
            this_sub = cls.grid_objects_types[this_sub_mask, :]
            if this_sub.shape[0] <= 3:
                # not enough element
                continue
            if (this_sub[:, obj_col] == -1).all():
                # no load
                continue
            if ((this_sub[:, cls.LOR_COL] != -1) | (this_sub[:, cls.LEX_COL] != -1)).sum() <= 1:
                # only 1 line
                continue
            el_id = this_sub[this_sub[:, obj_col] != -1, obj_col][0]
            if (this_sub[:, cls.LOR_COL] != -1).any():
                line_or_id = this_sub[this_sub[:, cls.LOR_COL] != -1, cls.LOR_COL][0]
                line_ex_id = None
            else:
                line_or_id = None
                line_ex_id = this_sub[this_sub[:, cls.LEX_COL] != -1, cls.LEX_COL][0]
            res = (sub_id, el_id, line_or_id, line_ex_id)
            break
        return res
    
    @staticmethod
    def _aux_find_sub_shunt(env):
        """find a sub with 4 elements, the type of elements and at least 2 lines"""
        cls = type(env)
        res = None
        for el_id in range(cls.n_shunt):
            sub_id = cls.shunt_to_subid[el_id]
            this_sub_mask = cls.grid_objects_types[:,cls.SUB_COL] == sub_id
            this_sub = cls.grid_objects_types[this_sub_mask, :]
            if this_sub.shape[0] <= 3:
                # not enough element
                continue
            if ((this_sub[:, cls.LOR_COL] != -1) | (this_sub[:, cls.LEX_COL] != -1)).sum() <= 1:
                # only 1 line
                continue
            if (this_sub[:, cls.LOR_COL] != -1).any():
                line_or_id = this_sub[this_sub[:, cls.LOR_COL] != -1, cls.LOR_COL][0]
                line_ex_id = None
            else:
                line_or_id = None
                line_ex_id = this_sub[this_sub[:, cls.LEX_COL] != -1, cls.LEX_COL][0]
            res = (sub_id, el_id, line_or_id, line_ex_id)
            break
        return res
        
    def test_move_load(self):
        cls = type(self.env)            
        res = self._aux_find_sub(self.env, cls.LOA_COL)
        if res is None:
            raise RuntimeError("Cannot carry the test 'test_move_load' as "
                               "there are no suitable subastation in your grid.")
        (sub_id, el_id, line_or_id, line_ex_id) = res
        for new_bus in self.list_loc_bus:
            if line_or_id is not None:
                act = self.env.action_space({"set_bus": {"loads_id": [(el_id, new_bus)], "lines_or_id": [(line_or_id, new_bus)]}})
            else:
                act = self.env.action_space({"set_bus": {"loads_id": [(el_id, new_bus)], "lines_ex_id": [(line_ex_id, new_bus)]}})
            bk_act = self.env._backend_action_class()
            bk_act += act
            self.env.backend.apply_action(bk_act)
            global_bus = sub_id + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_loads()[el_id].bus_id == global_bus, f"error for new_bus {new_bus}: {self.env.backend._grid.get_loads()[el_id].bus_id} vs {global_bus}"
                if line_or_id is not None:
                    assert self.env.backend._grid.get_lines()[line_or_id].bus_1_id == global_bus
                else:
                    assert self.env.backend._grid.get_lines()[line_ex_id].bus_2_id == global_bus
                assert self.env.backend._grid.get_bus_status()[global_bus]
            else:
                assert not self.env.backend._grid.get_loads()[el_id].connected
                if line_or_id is not None:
                    assert not self.env.backend._grid.get_lines()[line_or_id].connected_global
                else:
                    assert not self.env.backend._grid.get_lines()[line_ex_id].connected_global
            topo_vect = 1 * self.env.backend.get_topo_vect()
            assert topo_vect[cls.load_pos_topo_vect[el_id]] == new_bus, f"{topo_vect[cls.load_pos_topo_vect[el_id]]} vs {new_bus}"
        
    def test_move_gen(self):
        cls = type(self.env)            
        res = self._aux_find_sub(self.env, cls.GEN_COL)
        if res is None:
            raise RuntimeError("Cannot carry the test 'test_move_gen' as "
                               "there are no suitable subastation in your grid (>= 4 elements with 2 lines and a gen).")
        (sub_id, el_id, line_or_id, line_ex_id) = res
        for new_bus in self.list_loc_bus:
            if line_or_id is not None:
                act = self.env.action_space({"set_bus": {"generators_id": [(el_id, new_bus)], "lines_or_id": [(line_or_id, new_bus)]}})
            else:
                act = self.env.action_space({"set_bus": {"generators_id": [(el_id, new_bus)], "lines_ex_id": [(line_ex_id, new_bus)]}})
            bk_act = self.env._backend_action_class()
            bk_act += act
            self.env.backend.apply_action(bk_act)
            global_bus = sub_id + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_generators()[el_id].bus_id == global_bus, f"error for new_bus {new_bus}: {self.env.backend._grid.get_generators()[el_id].bus_id} vs {global_bus}"
                if line_or_id is not None:
                    assert self.env.backend._grid.get_lines()[line_or_id].bus_1_id == global_bus
                else:
                    assert self.env.backend._grid.get_lines()[line_ex_id].bus_2_id == global_bus
                assert self.env.backend._grid.get_bus_status()[global_bus]
            else:
                assert not self.env.backend._grid.get_generators()[el_id].connected
                if line_or_id is not None:
                    line_id = line_or_id
                else:
                    line_id = line_ex_id
                assert not self.env.backend._grid.get_lines()[line_id].connected_global
            topo_vect = 1 * self.env.backend.get_topo_vect()
            assert topo_vect[cls.gen_pos_topo_vect[el_id]] == new_bus, f"{topo_vect[cls.gen_pos_topo_vect[el_id]]} vs {new_bus}"
        
    def test_move_storage(self):
        cls = type(self.env)            
        res = self._aux_find_sub(self.env, cls.STORAGE_COL)
        if res is None:
            raise RuntimeError("Cannot carry the test 'test_move_storage' as "
                               "there are no suitable subastation in your grid.")
        (sub_id, el_id, line_or_id, line_ex_id) = res
        for new_bus in self.list_loc_bus:
            if line_or_id is not None:
                act = self.env.action_space({"set_bus": {"storages_id": [(el_id, new_bus)], "lines_or_id": [(line_or_id, new_bus)]}})
            else:
                act = self.env.action_space({"set_bus": {"storages_id": [(el_id, new_bus)], "lines_ex_id": [(line_ex_id, new_bus)]}})
            bk_act = self.env._backend_action_class()
            bk_act += act
            self.env.backend.apply_action(bk_act)
            global_bus = sub_id + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_storages()[el_id].bus_id == global_bus, f"error for new_bus {new_bus}: {self.env.backend._grid.get_sotrages()[el_id].bus_id} vs {global_bus}"
                if line_or_id is not None:
                    assert self.env.backend._grid.get_lines()[line_or_id].bus_1_id == global_bus
                else:
                    assert self.env.backend._grid.get_lines()[line_ex_id].bus_2_id == global_bus
                assert self.env.backend._grid.get_bus_status()[global_bus]
            else:
                assert not self.env.backend._grid.get_storages()[el_id].connected
                if line_or_id is not None:
                    assert not self.env.backend._grid.get_lines()[line_or_id].connected_global
                else:
                    assert not self.env.backend._grid.get_lines()[line_ex_id].connected_global
            topo_vect = 1 * self.env.backend.get_topo_vect()
            assert topo_vect[cls.storage_pos_topo_vect[el_id]] == new_bus, f"{topo_vect[cls.storage_pos_topo_vect[el_id]]} vs {new_bus}"
    
    def test_move_line_or(self):
        cls = type(self.env)            
        line_id = 0
        bk_act = self.env._backend_action_class()   # TODO: this should be init with the state of the gridmodel (especially topology)
        for new_bus in self.list_loc_bus:
            act = self.env.action_space({"set_bus": {"lines_or_id": [(line_id, new_bus)]}})
            bk_act += act
            self.env.backend.apply_action(bk_act)
            global_bus = cls.line_or_to_subid[line_id] + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_lines()[line_id].bus_1_id == global_bus
                assert self.env.backend._grid.get_bus_status()[global_bus]
            else:
                assert not self.env.backend._grid.get_lines()[line_id].connected_global
            topo_vect = 1 * self.env.backend.get_topo_vect()
            assert topo_vect[cls.line_or_pos_topo_vect[line_id]] == new_bus, f"{topo_vect[cls.line_or_pos_topo_vect[line_id]]} vs {new_bus}"
                
    def test_move_line_ex(self):
        cls = type(self.env)            
        line_id = 0
        bk_act = self.env._backend_action_class()   # TODO: this should be init with the state of the gridmodel (especially topology)
        for new_bus in self.list_loc_bus:
            act = self.env.action_space({"set_bus": {"lines_ex_id": [(line_id, new_bus)]}})
            bk_act += act
            self.env.backend.apply_action(bk_act)
            global_bus = cls.line_ex_to_subid[line_id] + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_lines()[line_id].bus_2_id == global_bus
                assert self.env.backend._grid.get_lines()[line_id].connected_global
                assert self.env.backend._grid.get_lines()[line_id].connected_1
                assert self.env.backend._grid.get_lines()[line_id].connected_2
                assert self.env.backend._grid.get_bus_status()[global_bus], f"for new_bus = {new_bus} and line_id={line_id}"
            else:
                assert not self.env.backend._grid.get_lines()[line_id].connected_global
                assert not self.env.backend._grid.get_lines()[line_id].connected_1
                assert not self.env.backend._grid.get_lines()[line_id].connected_2
            topo_vect = 1 * self.env.backend.get_topo_vect()
            assert topo_vect[cls.line_ex_pos_topo_vect[line_id]] == new_bus, f"{topo_vect[cls.line_ex_pos_topo_vect[line_id]]} vs {new_bus}"
            
    def test_move_shunt(self):
        cls = type(self.env)            
        res = self._aux_find_sub_shunt(self.env)
        if res is None:
            raise RuntimeError("Cannot carry the test 'test_move_shunt' as "
                               "there are no suitable subastation in your grid (>= 4 elements with a shunt).")
        (sub_id, el_id, line_or_id, line_ex_id) = res
        assert self.env.backend._grid.get_shunts()[el_id].sub_id == sub_id, f"{self.env.backend._grid.get_shunts()[el_id].sub_id} vs {sub_id}"
        bk_act = self.env._backend_action_class()  # TODO: this should be init with the state of the gridmodel (especially topology)
        for new_bus in self.list_loc_bus:
            if line_or_id is not None:
                act = self.env.action_space({"shunt": {"set_bus": [(el_id, new_bus)]}, "set_bus": {"lines_or_id": [(line_or_id, new_bus)]}})
            else:
                act = self.env.action_space({"shunt": {"set_bus": [(el_id, new_bus)]}, "set_bus": {"lines_ex_id": [(line_ex_id, new_bus)]}})
            bk_act += act
            self.env.backend.apply_action(bk_act)
            bk_act.reset()
            self.env.backend._set_shunt_info()
            sh_p, sh_q, sh_v, sh_bus = self.env.backend.shunt_info()
            assert sh_bus[el_id] == new_bus, f"{sh_bus[el_id]} vs {new_bus}"
            global_bus = sub_id + (new_bus -1) * cls.n_sub 
            if new_bus >= 1:
                assert self.env.backend._grid.get_shunts()[el_id].bus_id == global_bus, f"error for new_bus {new_bus}: {self.env.backend._grid.get_shunts()[el_id].bus_id} vs {global_bus}"
                if line_or_id is not None:
                    assert self.env.backend._grid.get_lines()[line_or_id].bus_1_id == global_bus
                else:
                    assert self.env.backend._grid.get_lines()[line_ex_id].bus_2_id == global_bus
                assert self.env.backend._grid.get_bus_status()[global_bus]
            else:
                assert not self.env.backend._grid.get_shunts()[el_id].connected
                if line_or_id is not None:
                    assert not self.env.backend._grid.get_lines()[line_or_id].connected_global
                else:
                    assert not self.env.backend._grid.get_lines()[line_ex_id].connected_global
    
    def test_check_kirchoff(self):
        cls = type(self.env)            
        res = self._aux_find_sub(self.env, cls.LOA_COL)
        if res is None:
            raise RuntimeError("Cannot carry the test 'test_move_load' as "
                               "there are no suitable subastation in your grid.")
        (sub_id, el_id, line_or_id, line_ex_id) = res
        for new_bus in self.list_loc_bus:
            if new_bus <= -1:
                continue
            if line_or_id is not None:
                act = self.env.action_space({"set_bus": {"loads_id": [(el_id, new_bus)], "lines_or_id": [(line_or_id, new_bus)]}})
            else:
                act = self.env.action_space({"set_bus": {"loads_id": [(el_id, new_bus)], "lines_ex_id": [(line_ex_id, new_bus)]}})
            bk_act = self.env._backend_action_class()
            bk_act += act
            self.env.backend.apply_action(bk_act)
            conv, maybe_exc = self.env.backend.runpf()
            assert conv, f"error : {maybe_exc}"
            p_subs, q_subs, p_bus, q_bus, diff_v_bus = self.env.backend.check_kirchhoff()
            # assert laws are met
            assert np.abs(p_subs).max() <= 1e-5, f"error for busbar {new_bus}: {np.abs(p_subs).max():.2e}"
            assert np.abs(q_subs).max() <= 1e-5, f"error for busbar {new_bus}: {np.abs(q_subs).max():.2e}"
            assert np.abs(p_bus).max() <= 1e-5, f"error for busbar {new_bus}: {np.abs(p_bus).max():.2e}"
            assert np.abs(q_bus).max() <= 1e-5, f"error for busbar {new_bus}: {np.abs(q_bus).max():.2e}"
            assert np.abs(diff_v_bus).max() <= 1e-5, f"error for busbar {new_bus}: {np.abs(diff_v_bus).max():.2e}"


class TestLightSimBackend_1busbar(TestLightSimBackend_3busbars):
    def get_nb_bus(self):
        return 1
    

class TestLightSimBackend_3busbars_iidm(TestLightSimBackend_3busbars):    
    def get_env_nm(self):
        return "./case_14_storage_iidm"
    
    def get_gen_slack_id(self):
        return 5
    
    def get_backend_kwargs(self):
        return dict(loader_method="pypowsybl",
                    gen_slack_id=self.get_gen_slack_id(),
                    loader_kwargs={"use_buses_for_sub": True,
                                   "double_bus_per_sub": True}
                    )
    
    
if __name__ == "__main__":
    unittest.main()
