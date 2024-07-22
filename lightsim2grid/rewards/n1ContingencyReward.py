# Copyright (c) 2020-2024, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import time
import numpy as np

import grid2op
from grid2op.Reward import BaseReward
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
                 logger=None,
                 tol=1e-8,
                 nb_iter=10):
        BaseReward.__init__(self, logger=logger)
        self._backend : LightSimBackend = None
        self._backend_action = None
        self._l_ids = None
        self._dc : bool = dc
        self._normalize : bool = normalize
        if l_ids is not None:
            self._l_ids = [int(el) for el in l_ids]
        self._threshold_margin :float = float(threshold_margin)
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
        self._backend_ls = False
        self._tol = tol
        self._nb_iter = nb_iter
        self._timer_call = 0.
        self._timer_pre_proc = 0.
        self._timer_compute = 0.
        self._timer_post_proc = 0.
            
    def initialize(self, env: "grid2op.Environment.Environment"):
        from grid2op.Environment import BaseEnv
        from grid2op.Backend import PandaPowerBackend  # lazy import because grid2op -> pandapower-> lightsim2grid -> grid2op
        if not isinstance(env, BaseEnv):
            raise RuntimeError("You can only initialize this reward with a "
                               "proper grid2op environment (`BaseEnv`)")
             
        if not isinstance(env.backend, (PandaPowerBackend, LightSimBackend)):
            raise RuntimeError("Impossible to use the `N1ContingencyReward` with "
                               "an environment with a backend that is not "
                               "``PandaPowerBackend` nor `LightSimBackend`."
                               )
        if isinstance(env.backend, LightSimBackend):
            self._backend : LightSimBackend = env.backend.copy()
            self._backend_ls : bool  = True
        elif isinstance(env.backend, PandaPowerBackend):
            self._backend = LightSimBackend.init_grid(type(env.backend))()
            self._backend.init_from_loaded_pandapower(env.backend)
            self._backend.is_loaded = True
        else:
            raise NotImplementedError()
        
        self._backend.set_solver_type(self._solver_type)
        conv, exc_ = self._backend.runpf()
        if not conv:
            raise RuntimeError(f"The reward N1ContingencyReward diverge with error {exc_}")
        bk_act_cls = _BackendAction.init_grid(type(env.backend))
        self._backend_action = bk_act_cls()
        if self._l_ids is None:
            self._l_ids = list(range(type(env).n_line))
        
        if len(self._l_ids) == 0:
            raise RuntimeError("Impossible to use the N1ContingencyReward "
                               "without any contingencies !")
        self.reward_min = 0.
        self.reward_max = len(self._l_ids) if not self._normalize else 1.
        # self._contingecy_analyzer = ContingencyAnalysis(self._backend)
        # self._contingecy_analyzer.add_multiple_contingencies(self._l_ids)

    def __call__(self, action, env, has_error, is_done, is_illegal, is_ambiguous):
        if is_done:
            return self.reward_min
        
        beg = time.perf_counter()
        # retrieve the state of the grid
        self._backend_action.reset()
        act = env.backend.get_action_to_set()
        th_lim_a = 1. * env.get_thermal_limit()
        th_lim_a[th_lim_a <= 1.] = 1.  # assign 1 for the thermal limit
        
        # apply it to the backend
        self._backend_action += act
        self._backend.apply_action(self._backend_action)
        conv, exc_ = self._backend.runpf()
        if not conv:
            self.logger.warn("Cannot set the backend of the `N1ContingencyReward` => divergence")
            return self.reward_min
        
        # synch the contingency analyzer
        contingecy_analyzer = ContingencyAnalysis(self._backend)
        contingecy_analyzer.computer.change_solver(self._solver_type)
        contingecy_analyzer.add_multiple_contingencies(*self._l_ids)
        now_ = time.perf_counter()
        self._timer_pre_proc += now_ - beg
        tmp = contingecy_analyzer.get_flows()
        self.logger.info(f"{contingecy_analyzer.computer.nb_solved()} converging contingencies")
        now_2 = time.perf_counter()
        self._timer_compute += now_2 - now_
        if self._dc:
            # In DC is study p, but take into account q in the limits
            tmp_res = np.abs(tmp[0])  # this is Por
            # now transform the limits in A in MW
            por, qor, vor, aor = env.backend.lines_or_info()
            p_sq = (1e-3 * th_lim_a)**2 * 3. * vor**2 - qor**2
            p_sq[p_sq <= 0.] = 0.
            limits = np.sqrt(p_sq)
        else:
            tmp_res = 1. * tmp[1]
            limits = th_lim_a
        res = ((tmp_res > self._threshold_margin * limits) | (~np.isfinite(tmp_res))).any(axis=1)  # whether one powerline is above its limit, per cont
        res |=  (np.abs(tmp_res) <= self._tol).all(axis=1)  # other type of divergence: all 0.
        res = res.sum()  # count total of n-1 unsafe 
        res = len(self._l_ids) - res  # reward = things to maximise
        if self._normalize:
            res /= len(self._l_ids)
        now_3 = time.perf_counter()
        self._timer_post_proc += now_3 - now_2
        self._timer_call += time.perf_counter() - beg
        return res

    def reset(self, env: "grid2op.Environment.BaseEnv") -> None:
        self._timer_call = 0.
        self._timer_pre_proc = 0.
        self._timer_compute = 0.
        self._timer_post_proc = 0.
        return super().reset(env)
    
    def close(self):
        if self._backend is not None:
            self._backend.close()
        del self._backend
        self._backend = None
