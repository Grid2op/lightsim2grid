from grid2op import make
from grid2op.Agent import DoNothingAgent
from grid2op.Parameters import Parameters
from grid2op.Rules import AlwaysLegal
from grid2op.Chronics import ChangeNothing
from lightsim2grid.LightSimBackend import LightSimBackend
import numpy as np
import unittest
import warnings
import pdb


class TestDN(unittest.TestCase):
    def setUp(self):
        self.param = Parameters()
        self.param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        self.max_ts = 100
        self.tol = 1e-5
        self.env_name = "case5_example"
        backend = LightSimBackend()
        env_name = self.env_name
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            self.env = make(self.env_name, param=self.param, backend=backend,
                            gamerules_class=AlwaysLegal,
                            chronics_class=ChangeNothing)

    def tearDown(self):
        self.env.close()

    def test_discoline(self):
        return
        # disconnect every powerline one by one, check that it is properly disconnected
        # then reconnect it, and check it is reconneccted
        env = self.env
        backend = env.backend
        for l_id in range(backend.n_line):
            if l_id == 2:
                # this makes the powerflow diverge when disconnected for case5_example
                continue
            action = env.action_space.disconnect_powerline(l_id)
            obs, reward, done, info = env.step(action)
            assert not done, "divergence after disconnection of powerline {}".format(l_id)
            pos_or = backend.line_or_pos_topo_vect[l_id]
            pos_ex = backend.line_ex_pos_topo_vect[l_id]
            all_rest = np.full(env.action_space.dim_topo, fill_value=True, dtype=np.bool)
            all_rest[pos_or] = all_rest[pos_ex] = False
            assert backend.topo_vect[pos_or] == -1
            assert backend.topo_vect[pos_ex] == -1
            assert np.all(backend.topo_vect[all_rest] == 1)
            assert obs.a_or[l_id] == 0
            assert obs.a_ex[l_id] == 0
            assert not obs.line_status[l_id], "line is not disconnected for {}".format(l_id)
            assert np.sum(obs.line_status) == env.n_line -1

            action = env.action_space.reconnect_powerline(l_id, bus_or=1, bus_ex=1)
            obs, reward, done, info = env.step(action)
            assert not done, "divergence after reconnection of powerline {}".format(l_id)
            all_rest = np.full(env.action_space.dim_topo, fill_value=True, dtype=np.bool)
            assert backend.topo_vect[pos_or] == 1
            assert backend.topo_vect[pos_ex] == 1
            assert np.all(backend.topo_vect[all_rest] == 1)
            assert obs.line_status[l_id]
            assert np.sum(obs.line_status) == env.n_line

    def test_topo_action(self):
        return
        env = self.env
        backend = env.backend
        action_restore = env.action_space({"set_bus": {"substations_id": [(2, np.array([1,1,1,1]))]}})

        l_1 = 4
        l_2 = 6
        new_topo = np.array([1, 2, 1, 2])
        action = env.action_space({"set_bus": {"substations_id": [(2, new_topo)]}})
        pos_1 = backend.line_ex_pos_topo_vect[l_1]
        pos_2 = backend.line_or_pos_topo_vect[l_2]
        obs, reward, done, info = env.step(action)
        assert backend.topo_vect[pos_1] == 2
        assert backend.topo_vect[pos_2] == 2
        assert np.sum(obs.line_status) == env.n_line
        assert not done
        assert np.all(backend.nb_obj_per_bus == [3, 3, 2, 3, 2, 0, 0, 2, 0, 0])

    def test_disco_line_6(self):
        env = self.env
        backend = env.backend
        action_space = env.action_space
        l_id = 6
        action = env.action_space.disconnect_powerline(l_id)
        obs, reward, done, info = env.step(action)
        assert not done, "divergence after disconnection of powerline {}".format(l_id)

        action_reco = action_space({"set_bus": np.ones(action_space.dim_topo, dtype=np.int),
                                    "set_line_status": np.ones(action_space.n_line, dtype=np.int)})
        obs, reward, done, info = env.step(action_reco)
        assert not done, "divergence after disconnection of powerline {}".format(l_id)
        pos_or = backend.line_or_pos_topo_vect[l_id]
        pos_ex = backend.line_ex_pos_topo_vect[l_id]
        all_rest = np.full(env.action_space.dim_topo, fill_value=True, dtype=np.bool)
        assert backend.topo_vect[pos_or] == 1
        assert backend.topo_vect[pos_ex] == 1
        assert np.all(backend.topo_vect[all_rest] == 1)
        assert obs.line_status[l_id]
        assert np.sum(obs.line_status) == env.n_line


if __name__ == "__main__":
    unittest.main()
