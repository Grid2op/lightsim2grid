# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
import numpy as np
import warnings

import grid2op
from lightsim2grid.LightSimBackend import LightSimBackend
from grid2op.Chronics import Multifolder, GridStateFromFile
from lightsim2grid_cpp import Computers


class TimeSerie:
    """
    This helper class, that only works with grid2op when using a LightSimBackend allows to compute
    the flows (at the origin side of the powerline / transformers). It is roughly equivalent to the
    grid2op code:

    .. code-block:: python

        import grid2op
        import numpy as np
        from grid2op.Parameters import Parameters
        from lightsim2grid.LightSimBackend import LightSimBackend

        env_name = ...
        param = Parameters()
        param.NO_OVERFLOW_DISCONNECTION = True

        env = grid2op.make(env_name, param=param, backend=LightSimBackend())

        done = False
        obs = env.reset()
        nb_step = obs.max_step
        Vs = np.zeros((nb_step, 2 * env.n_sub), dtype=complex)
        As = np.zeros((nb_step, env.n_line), dtype=float)
        while not done:
            obs, reward, done, info = env.step(env.action_space())
            Vs[i, :env.n_sub] = env.backend.V
            As[i] = obs.a_or

    Compare to the previous code, it avoid all grid2op code (coded in python) and can be roughly 3-5 times 
    faster.
    It allows also to use python threading, as the c++ computation can be done in different thread.

    Examples
    ----------

    It can be used as:

    .. code-block:: python

        from lightsim2grid import TimeSerie
        import grid2op
        from lightsim2grid.LightSimBackend import LightSimBackend

        env_name = ...
        env = grid2op.make(env_name, param=param, backend=LightSimBackend())

        time_series = TimeSerie(env)
        Vs = time_series.compute_V(scenario_id=..., seed=...)
        As = time_series.compute_A()

    """
    def __init__(self, grid2op_env):
        if not isinstance(grid2op_env.backend, LightSimBackend):
            raise RuntimeError("This class only works with LightSimBackend")
        self.grid2op_env = grid2op_env.copy()
        self.computer = Computers(self.grid2op_env.backend._grid)
        self.prod_p = None
        self.load_p = None
        self.load_q = None
        self.__computed = False
    
    def get_injections(self, scenario_id=None, seed=None):
        """
        This function allows to retrieve the injection of the given scenario, for the given seed
        from the grid2op internal environment.
        """
        if scenario_id is not None:
            self.grid2op_env.set_id(scenario_id)
        if seed is not None:
            self.grid2op_env.seed(seed)
        self.grid2op_env.reset()
        self.__computed = False
        return self._extract_inj()

    def compute_V(self, scenario_id=None, seed=None, ignore_errors=False):
        """
        This function allows to retrieve the complex voltage at each bus of the grid for each step.
        """
        prod_p, load_p, load_q = self.get_injections()
        v_init = self.grid2op_env.backend.V
        status = self.computer.compute_Vs(prod_p,
                                          np.zeros((prod_p.shape[0], 0)),  # no static generators for now !
                                          load_p,
                                          load_q,
                                          v_init,
                                          self.grid2op_env.backend.max_it,
                                          self.grid2op_env.backend.tol)
        if status != 1 and not ignore_errors:
            # raise an error if the powerflow diverged
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        elif status != 1:
            # only raise a warning in this case
            warnings.warn(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        
        Vs = self.computer.get_voltages()
        self.__computed = True
        return Vs

    def compute_A(self):
        """
        This function returns the current flows (in Amps, A) at the origin / high voltage side
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        ampss = self.computer.compute_flows()
        return 1000. * ampss

    def _extract_inj(self):
        data_loader = None
        if isinstance(self.grid2op_env.chronics_handler.real_data, Multifolder):
            data_loader = self.grid2op_env.chronics_handler.real_data.data
        else:
            data_loader = self.grid2op_env.chronics_handler.data

        if not isinstance(data_loader, GridStateFromFile):
            raise RuntimeError("This function only work with chronics coming from files at the moment")
        
        self.prod_p = 1.0 * data_loader.prod_p
        self.load_p = 1.0 * data_loader.load_p
        self.load_q = 1.0 * data_loader.load_q
        return self.prod_p, self.load_p, self.load_q
