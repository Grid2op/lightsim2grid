#!/bin/python3
import grid2op
from lightsim2grid import LightSimBackend
import warnings

# usage:
# perf record ./test_profile.py
# perf report
env_name = "l2rpn_neurips_2020_track2_small"
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    env_tmp = grid2op.make(env_name)

    param = env_tmp.parameters
    param.NO_OVERFLOW_DISCONNECTION = True
    env = grid2op.make(env_name, backend=LightSimBackend(), param=param)


def make_steps(env, nb=1000):
    for i in range(nb):
        _ = env.step(env.action_space())


if __name__ == "__main__":
    make_steps(env, nb=1000)
