from grid2op import make
from grid2op.Agent import AgentWithConverter, DoNothingAgent
from grid2op.Parameters import Parameters
from grid2op.Rules import AlwaysLegal
from grid2op.Converter import IdToAct
from pyklu.PyKLUBackend import PyKLUBackend
import time
from tqdm import tqdm
import numpy as np
import sys
import pdb


ENV_NAME = "case5_example"
max_ts = 1000 #000  # env_me.chronics_handler.max_timestep()


class TestAgent(AgentWithConverter):
    def __init__(self, action_space, action_space_converter=IdToAct, **kwargs_converter):
        AgentWithConverter.__init__(self, action_space, action_space_converter=action_space_converter, **kwargs_converter)
        self.action_space.all_actions = []

        # do nothing
        # powerline switch: disconnection
        for i in range(action_space.n_line):
            self.action_space.all_actions.append(action_space.disconnect_powerline(line_id=i))

        # other type of actions
        # self.action_space.all_actions += action_space.get_all_unitary_topologies_set(action_space)
        # self.action_space.all_actions += action_space.get_all_unitary_redispatch(action_space)

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

# TestAgent = DoNothingAgent

backend = PyKLUBackend()
param = Parameters()
param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

env_me = make(ENV_NAME, backend=backend, param=param, gamerules_class=AlwaysLegal)
aor_me = np.zeros((max_ts, env_me.n_line))
aor_pp = np.zeros((max_ts, env_me.n_line))
agent = TestAgent(action_space=env_me.action_space)
obs = env_me.get_obs()
done = False
reward = env_me.reward_range[0]
nb_ts = 0
prev_act = None
with tqdm(total=max_ts) as pbar:
    while not done:
        act = agent.act(obs, reward, done)
        # print("new action at time step {}".format(nb_ts))
        # print("{}".format(act))
        obs, reward, done, info = env_me.step(act)
        aor_me[nb_ts, :] = obs.a_or
        nb_ts += 1
        pbar.update(1)
        if np.sum(obs.line_status) < obs.n_line - 1 * (nb_ts % 2 == 1):
            pdb.set_trace()
        if nb_ts >= max_ts:
            break
        if done:
            print("Action performed: {}".format(act))
            print("Error message: {}".format(info))
            # pdb.set_trace()
        prev_act = act

print("My Backend {} time steps".format(nb_ts))
print("\tTime apply act: {:.2f}ms".format(1000.*env_me._time_apply_act/nb_ts))
print("\tTime powerflow: {:.2f}ms".format(1000.*env_me._time_powerflow/nb_ts))
print("\tTime extract observation: {:.2f}ms".format(1000.*env_me._time_extract_obs/nb_ts))

env_pp = make(ENV_NAME, param=param, gamerules_class=AlwaysLegal)
agent = TestAgent(action_space=env_pp.action_space)
obs = env_pp.get_obs()
done = False
reward = env_pp.reward_range[0]
nb_ts = 0
prev_act = None
with tqdm(total=max_ts) as pbar:
    while not done:
        act = agent.act(obs, reward, done)
        obs, reward, done, info = env_pp.step(act)
        aor_pp[nb_ts, :] = obs.a_or
        nb_ts += 1
        pbar.update(1)
        if np.sum(obs.line_status) < obs.n_line - 1 * (nb_ts % 2 == 1):
            pdb.set_trace()
        if nb_ts >= max_ts:
            break
        if done:
            pdb.set_trace()
            print("Action performed: {}".format(act))
            print("Error message: {}".format(info))
        prev_act = act

print("Pandapower Backend {} time steps".format(nb_ts))
print("\tTime apply act: {:.2f}ms".format(1000.*env_pp._time_apply_act/nb_ts))
print("\tTime powerflow: {:.2f}ms".format(1000.*env_pp._time_powerflow/nb_ts))
print("\tTime extract observation: {:.2f}ms".format(1000.*env_pp._time_extract_obs/nb_ts))

print("Absolute value of the difference: {}".format(np.max(np.abs(aor_me - aor_pp))))
pdb.set_trace()
