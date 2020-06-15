# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid a implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from grid2op import make
from grid2op.Agent import AgentWithConverter, DoNothingAgent
from grid2op.Parameters import Parameters
from grid2op.Converter import IdToAct
from grid2op.Rules import AlwaysLegal
from grid2op.Chronics import GridStateFromFile

from lightsim2grid.LightSimBackend import LightSimBackend

from utils_benchmark import print_res, run_env, str2bool

ENV_NAME = "rte_case5_example"
ENV_NAME = "rte_case14_realistic"
MAX_TS = 1000


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
                             action_space({"set_bus": {"lines_or_id": [(100, 2), (129, 1), (173, 2)],
                                                       # "lines_ex_id": [(2, 2), (13, 2), (24, 2), (35, 2)],
                                                       "generators_id": [(2, 1)],
                                                       "loads_id": [(6, 2)]
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


def main(max_ts, name, test=True):
    backend = LightSimBackend()
    param = Parameters()
    param.init_from_dict({"NO_OVERFLOW_DISCONNECTION": True})

    env_klu = make(name, backend=backend, test=test,
                   param=param, gamerules_class=AlwaysLegal,
                   data_feeding_kwargs={"chunk_size": 128, "max_iter": max_ts, "gridvalueClass": GridStateFromFile})
    agent = TestAgent(action_space=env_klu.action_space, env_name=name)
    nb_ts_klu, time_klu, aor_klu, gen_p_klu, gen_q_klu = run_env(env_klu, max_ts, agent,
                                                                 chron_id=0, keep_forecast=False)

    env_pp = make(name, param=param, gamerules_class=AlwaysLegal, test=test,
                   data_feeding_kwargs={"chunk_size": 128, "max_iter": max_ts, "gridvalueClass": GridStateFromFile})
    agent = TestAgent(action_space=env_pp.action_space, env_name=name)
    nb_ts_pp, time_pp, aor_pp, gen_p_pp, gen_q_pp = run_env(env_pp, max_ts, agent,
                                                            chron_id=0, keep_forecast=False)

    print_res(env_klu, env_pp,
              nb_ts_klu, nb_ts_pp,
              time_klu, time_pp,
              aor_klu, aor_pp,
              gen_p_klu, gen_p_pp,
              gen_q_klu, gen_q_pp
              )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark pyKLU and Pandapower Backend for an agent that takes every '
                                                 'topological action possible')
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
