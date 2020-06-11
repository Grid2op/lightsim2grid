# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import numpy as np
import os

from grid2op import make
from grid2op.Agent import DoNothingAgent
from grid2op.Chronics import GridStateFromFile
from grid2op.Parameters import Parameters
from lightsim2grid.LightSimBackend import LightSimBackend
from utils_benchmark import print_res, run_env, str2bool
import pdb

MAX_TS = 1000
ENV_NAME = "rte_case14_realistic"


def main(max_ts, ENV_NAME, test=True):
    backend = LightSimBackend()
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

    env_klu = make(ENV_NAME, backend=backend, param=param, test=test,
                   data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
    agent = DoNothingAgent(action_space=env_klu.action_space)
    nb_ts_klu, time_klu, aor_klu, gen_p_klu, gen_q_klu = run_env(env_klu, max_ts, agent, chron_id=0)

    env_pp = make(ENV_NAME, param=param, test=test,
                   data_feeding_kwargs={"gridvalueClass": GridStateFromFile})
    agent = DoNothingAgent(action_space=env_pp.action_space)
    nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent, chron_id=0)

    print_res(env_klu, env_pp,
              nb_ts_klu, nb_ts_pp,
              time_klu, time_pp,
              aor_klu, aor_pp,
              gen_p_klu, gen_p_pp,
              gen_q_klu, gen_q_pp
              )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark pyKLU and Pandapower Backend for a "do nothing" agent')
    parser.add_argument('--name', default=ENV_NAME, type=str,
                        help='Environment name to be used for the benchmark.')
    parser.add_argument('--number', type=int, default=MAX_TS,
                        help='Maximum number of time steps for which the benchamark will be run.')
    parser.add_argument('--no_test', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='Do not use test environment for the benchmark (default False: use test environment)')

    args = parser.parse_args()

    max_ts = int(args.number)
    name = str(args.name)
    test_env = not args.no_test
    main(max_ts, name, test_env)
