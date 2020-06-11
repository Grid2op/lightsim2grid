from grid2op import make
from grid2op.Agent import AgentWithConverter
from grid2op.Parameters import Parameters
from grid2op.Rules import AlwaysLegal
from grid2op.Converter import IdToAct
from lightsim2grid.LightSimBackend import LightSimBackend
import numpy as np
import unittest
from abc import ABC, abstractmethod
import warnings
import pdb

MAX_TS = 1000  # does not work on ieee118 (pp diverge at 455s)
MAX_TS = 100


class TestAgent(AgentWithConverter):
    def __init__(self, action_space, env_name, action_space_converter=IdToAct, **kwargs_converter):
        AgentWithConverter.__init__(self, action_space, action_space_converter=action_space_converter, **kwargs_converter)
        self.action_space.all_actions = []

        # do nothing
        all_actions_tmp = [action_space()]

        # powerline switch: disconnection
        for i in range(action_space.n_line):
            if env_name == "rte_case14_realistic":
                if i == 18:
                    continue
            elif env_name == "rte_case5_example":
                pass
            elif env_name == "rte_case118_example":
                if i == 6:
                    continue
                if i == 26:
                    continue
                if i == 72:
                    continue
                if i == 73:
                    continue
                if i == 80:
                    continue
                if i == 129:
                    continue
                if i == 140:
                    continue
                if i == 176:
                    continue
                if i == 177:
                    continue
            elif env_name == "l2rpn_wcci_2020":
                if i == 2:
                    continue
            all_actions_tmp.append(action_space.disconnect_powerline(line_id=i))

        # other type of actions
        all_actions_tmp += action_space.get_all_unitary_topologies_set(action_space)
        # self.action_space.all_actions += action_space.get_all_unitary_redispatch(action_space)

        if env_name == "rte_case14_realistic":
            # remove action that makes the powerflow diverge
            breaking_acts = [action_space({"set_bus": {"lines_or_id": [(7,2), (8,1), (9,1)],
                                                       "lines_ex_id": [(17,2)],
                                                       "generators_id": [(2,2)],
                                                       "loads_id": [(4,1)]}}),
                             action_space({"set_bus": {"lines_or_id": [(10, 2), (11, 1), (19,2)],
                                                       "lines_ex_id": [(16, 2)],
                                                       "loads_id": [(5, 1)]}}),
                             action_space({"set_bus": {"lines_or_id": [(5, 1)],
                                                       "lines_ex_id": [(2, 2)],
                                                       "generators_id": [(1, 2)],
                                                       "loads_id": [(1, 1)]}}),
                             action_space({"set_bus": {"lines_or_id": [(6, 2), (15, 2), (16, 1)],
                                                       "lines_ex_id": [(3, 2), (5, 2)],
                                                       "loads_id": [(2, 1)]}}),
                            action_space({"set_bus": {"lines_or_id": [(18, 1)],
                                                      "lines_ex_id": [(15, 2), (19, 2)],
                                                      }})
            ]
        elif env_name == "rte_case118_example":
            breaking_acts = [action_space({"set_bus": {"lines_or_id": [(100, 2), (129, 1), (173, 2)],
                                                       # "lines_ex_id": [(17,2)],
                                                       "generators_id": [(2, 2)],
                                                       "loads_id": [(6, 1)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(100, 2), (129, 1), (173, 2)],
                                                       # "lines_ex_id": [(17,2)],
                                                       "generators_id": [(2, 2)],
                                                       "loads_id": [(6, 2)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(100, 2), (129, 1), (173, 2)],
                                                       # "lines_ex_id": [(17,2)],
                                                       "generators_id": [(2, 1)],
                                                       "loads_id": [(6, 1)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(140, 1)],
                                                       "lines_ex_id": [(129, 2)],
                                                       # "generators_id": [(2, 1)],
                                                       # "loads_id": [(6, 1)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(57, 2), (80, 1), (83, 2)],
                                                       "lines_ex_id": [(2, 2), (13, 2), (24, 2), (35, 2)],
                                                       "generators_id": [(6, 2)],
                                                       "loads_id": [(8, 2)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(57, 2), (80, 1), (83, 2)],
                                                       "lines_ex_id": [(2, 2), (13, 2), (24, 2), (35, 2)],
                                                       "generators_id": [(6, 2)],
                                                       "loads_id": [(8, 1)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(57, 2), (80, 1), (83, 2)],
                                                       "lines_ex_id": [(2, 2), (13, 2), (24, 2), (35, 2)],
                                                       "generators_id": [(6, 1)],
                                                       "loads_id": [(8, 2)]
                                                       }}),
                             action_space({"set_bus": {"lines_or_id": [(57, 2), (80, 1), (83, 2)],
                                                       "lines_ex_id": [(2, 2), (13, 2), (24, 2), (35, 2)],
                                                       "generators_id": [(6, 1)],
                                                       "loads_id": [(8, 1)]
                                                       }}),
            ]

        elif env_name == "l2rpn_wcci_2020":
            breaking_acts = [action_space({"set_bus": {"lines_or_id": [(5, 2), (6, 2)],
                                                       "lines_ex_id": [(1, 2), (2, 1), (4, 2), (55, 2)],
                                                       }}),
            ]
        else:
            breaking_acts = [action_space({"set_bus": {"lines_or_id": [(0,2), (1,2), (2,2), (3,1)],
                                                       "generators_id": [(0,1)],
                                                       "loads_id": [(0,1)]}}),
                             ]

        # filter out actions that break everything
        all_actions = []
        for el in all_actions_tmp:
            if not el in breaking_acts:
                all_actions.append(el)

        # set the action to the action space
        self.action_space.all_actions = all_actions

        # add the action "reset everything to 1 bus"
        self.action_space.all_actions.append(action_space({"set_bus": np.ones(action_space.dim_topo, dtype=np.int),
                                                           "set_line_status": np.ones(action_space.n_line, dtype=np.int)}))
        self.nb_act_done = 0
        self.act_this = True

    def my_act(self, transformed_obs, reward, done=False):
        if self.act_this:
            res = self.nb_act_done
            self.nb_act_done += 1
            self.nb_act_done %= len(self.action_space.all_actions)
            self.act_this = False
        else:
            res = -1
            self.act_this = True
        return res


class TestAgentAllMove(ABC):
    def setUp(self):
        self.param = Parameters()
        self.param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})
        self.max_ts = MAX_TS
        self.tol = 1e-9
        self.tol_q = 2e-5
        self.agent_class = TestAgent

    @abstractmethod
    def _get_env_name(self):
        pass

    def _run_env(self, env):
        aor = np.zeros((self.max_ts, env.n_line))
        gen_p = np.zeros((self.max_ts, env.n_gen))
        gen_q = np.zeros((self.max_ts, env.n_gen))
        agent = self.agent_class(action_space=env.action_space, env_name=self._get_env_name())
        obs = env.get_obs()
        done = False
        reward = env.reward_range[0]
        nb_ts = 0
        while not done:
            act = agent.act(obs, reward, done)
            obs, reward, done, info = env.step(act)
            aor[nb_ts, :] = obs.a_or
            gen_p[nb_ts, :] = obs.prod_p
            gen_q[nb_ts, :] = obs.prod_q
            nb_ts += 1
            if nb_ts >= self.max_ts:
                break
        return nb_ts, aor, gen_p, gen_q

    def test_do_nothing(self):
        backend = LightSimBackend()
        env_name = self._get_env_name()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            with make(env_name, test=True,
                      param=self.param,
                      backend=backend,
                      gamerules_class=AlwaysLegal,
                      data_feeding_kwargs={"chunk_size": 128, "max_iter": MAX_TS}) as env:
                nb_ts_klu, aor_klu, gen_p_klu, gen_q_klu = self._run_env(env)
            with make(env_name,
                      test=True,
                      param=self.param,
                      gamerules_class=AlwaysLegal,
                      data_feeding_kwargs={"chunk_size": 128, "max_iter": MAX_TS}) as env:
                nb_ts_pp, aor_pp, gen_p_pp, gen_q_pp = self._run_env(env)

        assert nb_ts_klu == nb_ts_pp, "not same number of timesteps for {}: lightsim: {}, pp: {}" \
                                      "".format(env_name, nb_ts_klu, nb_ts_pp)
        assert np.max(np.abs(aor_klu - aor_pp)) <= self.tol, "aor l inf different for {}: {}" \
                                                             "".format(env_name, np.max(np.abs(aor_klu - aor_pp)))
        assert np.mean(np.abs(aor_klu - aor_pp)) <= self.tol, "aor l1 different for {}: {}" \
                                                              "".format(env_name, np.mean(np.abs(aor_klu - aor_pp)))
        assert np.max(np.abs(gen_p_klu - gen_p_pp)) <= self.tol, "gen_p l inf different for {}: {}" \
                                                                 "".format(env_name, np.max(np.abs(gen_p_klu - gen_p_pp)))
        assert np.mean(np.abs(gen_p_klu - gen_p_pp)) <= self.tol, "gen_p l1 different for {}: {}" \
                                                                  "".format(env_name, np.mean(np.abs(gen_p_klu - gen_p_pp)))

        # a slightly different algorithm compare to pandapower (see DataGen.cpp, set_q)
        assert np.max(np.abs(gen_q_klu - gen_q_pp)) <= self.tol_q, "l inf different for {} gen_q: {}" \
                                                                   "".format(env_name, np.max(np.abs(gen_q_klu - gen_q_pp)))
        assert np.mean(np.abs(gen_q_klu - gen_q_pp)) <= self.tol_q, "l1 different for {} gen_q".format(env_name)


class Testcase5(TestAgentAllMove, unittest.TestCase):
    def _get_env_name(self):
        return "rte_case5_example"


class Testcase14(TestAgentAllMove, unittest.TestCase):
    def _get_env_name(self):
        return "rte_case14_realistic"


# takes a looooong time
class Testcase118(TestAgentAllMove, unittest.TestCase):
    def _get_env_name(self):
        return "rte_case118_example"


# requires additional data to be downloaded
# class TestcaseSandbox(TestDN, unittest.TestCase):
#    def _get_env_name(self):
#        return "l2rpn_case14_sandbox"


if __name__ == "__main__":
    unittest.main()

