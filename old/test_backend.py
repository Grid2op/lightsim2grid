from grid2op import make
from grid2op.Agent import DoNothingAgent
from grid2op.Parameters import Parameters
from pyklu.PyKLUBackend import PyKLUBackend
import time
from tqdm import tqdm
import numpy as np
import sys

backend = PyKLUBackend()
param = Parameters()
param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

env_me = make(backend=backend, param=param)
# sys.exit()
agent = DoNothingAgent(action_space=env_me.action_space)
obs = env_me.get_obs()
done = False
reward = env_me.reward_range[0]
nb_ts = 0
max_ts = 1000  # env_me.chronics_handler.max_timestep()
aor_me = np.zeros((max_ts, env_me.n_line))
aor_pp = np.zeros((max_ts, env_me.n_line))
with tqdm(total=max_ts) as pbar:
    while not done:
        nb_ts += 1
        act = agent.act(obs, reward, done)
        obs, reward, done, info = env_me.step(act)
        aor_me[nb_ts-1, :] = obs.a_or
        pbar.update(1)
        if nb_ts >= max_ts:
            break
print("My Backend")
print("\tTime apply act: {:.2f}ms".format(1000.*env_me._time_apply_act/nb_ts))
print("\tTime powerflow: {:.2f}ms".format(1000.*env_me._time_powerflow/nb_ts))
print("\tTime extract observation: {:.2f}ms".format(1000.*env_me._time_extract_obs/nb_ts))

env_pp = make(param=param)
agent = DoNothingAgent(action_space=env_me.action_space)
obs = env_pp.get_obs()
done = False
reward = env_pp.reward_range[0]
nb_ts = 0
with tqdm(total=max_ts) as pbar:
    while not done:
        nb_ts += 1
        act = agent.act(obs, reward, done)
        obs, reward, done, info = env_pp.step(act)
        aor_pp[nb_ts-1, :] = obs.a_or
        pbar.update(1)
        if nb_ts >= max_ts:
            break
print("Pandapower Backend")
print("\tTime apply act: {:.2f}ms".format(1000.*env_pp._time_apply_act/nb_ts))
print("\tTime powerflow: {:.2f}ms".format(1000.*env_pp._time_powerflow/nb_ts))
print("\tTime extract observation: {:.2f}ms".format(1000.*env_pp._time_extract_obs/nb_ts))

print("Absolute value of the difference: {}".format(np.max(np.abs(aor_me - aor_pp))))