# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from grid2op.Reward import BaseReward
from grid2op.Environment import Environment
from grid2op.Backend import PandaPowerBackend
from grid2op.Action._backendAction import _BackendAction

from lightsim2grid import LightSimBackend, ContingencyAnalysis
from lightsim2grid.compilation_options import klu_solver_available
from lightsim2grid.solver import SolverType


class N1ContingencyReward(BaseReward):
    """
    This class implements a reward that leverage the :class:`lightsim2grid.ContingencyAnalysis`
    to compute the number of unsafe contingency at any given time.

    Examples
    --------

    This can be used as:

    .. code-block:: python

        import grid2op
        from lightsim2grid.rewards import N1ContingencyReward
        l_ids = [0, 1, 7]
        env = grid2op.make("l2rpn_case14_sandbox",
                           reward_class=N1ContingencyReward(l_ids=l_ids)
                          )
        obs = env.reset()
        obs, reward, *_ = env.step(env.action_space())
        print(f"reward: {reward:.3f}")
        
    """

    def __init__(self,
                 l_ids=None,
                 threshold_margin=1.,
                 dc=False,
                 normalize=False,
                 logger=None,):
        BaseReward.__init__(self, logger=logger)
        self._backend = None
        self._backend_action = None
        self._l_ids = None
        self._dc = dc
        self._normalize = normalize
        if l_ids is not None:
            self._l_ids = [int(el) for el in l_ids]
        self._threshold_margin = float(threshold_margin)
        if klu_solver_available:
            if self._dc:
                self._solver_type = SolverType.KLUDC
            else:
                self._solver_type = SolverType.KLU
        else:
            if self._dc:
                self._solver_type = SolverType.DC
            else:
                self._solver_type = SolverType.SparseLU
            
    def initialize(self, env: Environment):
        if not isinstance(env.backend, (PandaPowerBackend, LightSimBackend)):
            raise RuntimeError("Impossible to use the `N1ContingencyReward` with "
                               "a environment with a backend that is not "
                               "``PandaPowerBackend` nor `LightSimBackend`."
                               )
        if isinstance(env.backend, LightSimBackend):
            self._backend = env.backend.copy()
        elif isinstance(env.backend, PandaPowerBackend):
            from lightsim2grid.gridmodel import init
            self._backend = init(env.backend)
        else:
            raise NotImplementedError()
        
        bk_act_cls = _BackendAction.init_grid(env.backend)
        self._backend_action = bk_act_cls()
        if self._l_ids is None:
            self._l_ids = list(range(type(env).n_line))
            
        self.reward_min = 0.
        self.reward_max = type(env).n_line if not self._normalize else 1.
        self._contingecy_analyzer = ContingencyAnalysis(self._backend)
        self._contingecy_analyzer.add_multiple_contingencies(self._l_ids)

    def __call__(self, action, env, has_error, is_done, is_illegal, is_ambiguous):
        if is_done:
            return self.reward_min
        self._backend_action.reset()
        act = env.backend.get_action_to_set()
        th_lim = env.get_thermal_limit()
        th_lim[th_lim <= 1] = 1  # assign 1 for the thermal limit
        this_n1 = copy.deepcopy(act)
        self._backend_action += this_n1
            
        self._backend.apply_action(self._backend_action)
        self._backend._disconnect_line(self.l_id)
        div_exc_ = None
        try:
            # TODO there is a bug in lightsimbackend that make it crash instead of diverging
            conv, div_exc_ = self._backend.runpf()
        except Exception as exc_:
            conv = False
            div_exc_ = exc_

        if conv:
            flow = self._backend.get_line_flow()
            res = (flow / th_lim).max()
        else:
            self.logger.info(f"Divergence of the backend at step {env.nb_time_step} for N1Reward with error `{div_exc_}`")
            res = self.reward_min
        return res

    def close(self):
        self._backend.close()
        del self._backend
        self._backend = None
