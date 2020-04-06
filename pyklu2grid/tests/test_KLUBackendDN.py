from grid2op import make
from grid2op.Agent import DoNothingAgent
from grid2op.Parameters import Parameters
from grid2op.Rules import AlwaysLegal
from pyklu2grid.PyKLUBackend import PyKLUBackend
import numpy as np
import unittest
from abc import ABC, abstractmethod
import warnings
import pdb


class TestDN(ABC):
    def setUp(self):
        self.param = Parameters()
        self.param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        self.max_ts = 100
        self.tol = 1e-5
        self.agent_class = DoNothingAgent

    @abstractmethod
    def _get_env_name(self):
        pass

    def _run_env(self, env):
        aor = np.zeros((self.max_ts, env.n_line))
        agent = self.agent_class(action_space=env.action_space)
        obs = env.get_obs()
        done = False
        reward = env.reward_range[0]
        nb_ts = 0
        while not done:
            act = agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            aor[nb_ts, :] = obs.a_or
            nb_ts += 1
            if nb_ts >= self.max_ts:
                break
        return nb_ts, aor

    def test_do_nothing(self):
        backend = PyKLUBackend()
        env_name = self._get_env_name()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with make(env_name, param=self.param, backend=backend,  gamerules_class=AlwaysLegal) as env:
                nb_ts_klu, aor_klu = self._run_env(env)
            with make(env_name, param=self.param, gamerules_class=AlwaysLegal) as env:
                nb_ts_pp, aor_pp = self._run_env(env)

        assert nb_ts_klu == nb_ts_pp, "not same number of timesteps for {}".format(env_name)
        assert np.max(np.abs(aor_klu - aor_pp)) <= self.tol, "l inf different for {}".format(env_name)
        assert np.mean(np.abs(aor_klu - aor_pp)) <= self.tol, "l1 different for {}".format(env_name)


class Testcase5(TestDN, unittest.TestCase):
    def _get_env_name(self):
        return "case5_example"


class Testcase14(TestDN, unittest.TestCase):
    def _get_env_name(self):
        return "case14_realistic"


if __name__ == "__main__":
    unittest.main()

