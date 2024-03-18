# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["TimeSerieCPP", "TimeSerie",
           # deprecated
           "Computers"]

import numpy as np
import warnings

from grid2op.Chronics import Multifolder, GridStateFromFile

from lightsim2grid.lightSimBackend import LightSimBackend
from lightsim2grid.solver import SolverType
from lightsim2grid_cpp import TimeSeriesCPP

# deprecated
Computers = TimeSeriesCPP


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

    Compare to the previous code, it avoid all grid2op code and can be more than 15 times 
    faster (on the case 118).
    
    It also allows to use python threading module, as the c++ computation can be done in different python threads (the GIL is not locked
    during the c++ computation).

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
        res_p, res_a, res_v = time_series.get_flows(scenario_id=..., seed=...)

    """
    def __init__(self, grid2op_env):
        from grid2op.Environment import Environment  # otherwise i got issues...
        if not isinstance(grid2op_env.backend, LightSimBackend):
            raise RuntimeError("This class only works with LightSimBackend")
        if not isinstance(grid2op_env, Environment):
            raise RuntimeError("Please an environment of class \"Environment\", "
                               "and not \"MultimixEnv\" or \"BaseMultiProcessEnv\"")
        self.grid2op_env = grid2op_env.copy()
        self.computer = TimeSeriesCPP(self.grid2op_env.backend._grid)
        self.prod_p = None
        self.load_p = None
        self.load_q = None
        self.__computed = False
        
        self.available_solvers = self.computer.available_solvers()
        if SolverType.KLU in self.available_solvers:
            # use the faster KLU if available
            self.computer.change_solver(SolverType.KLU)
    
    def get_injections(self, scenario_id=None, seed=None):
        """
        This function allows to retrieve the injection of the given scenario, for the given seed
        from the grid2op internal environment.
        """
        if scenario_id is not None:
            self.grid2op_env.set_id(scenario_id)
        if seed is not None:
            self.grid2op_env.seed(seed)
        obs = self.grid2op_env.reset()
        self.__computed = False
        return self._extract_inj()

    def compute_V_from_inj(self, prod_p, load_p, load_q, v_init=None, ignore_errors=False):
        """
        This function allows to compute the voltages, at each bus given a list of
        productions and loads.

        We do not recommend to use it directly, as the order of the load or generators might vary !
        """

        if len(prod_p.shape) != 2:
            raise RuntimeError("prod_p should be a matrix with rows representing time steps "
                               "and columns representing individual production.")
        if len(load_p.shape) != 2:
            raise RuntimeError("load_p should be a matrix with rows representing time steps "
                               "and columns representing individual loads.")
        if len(load_q.shape) != 2:
            raise RuntimeError("load_q should be a matrix with rows representing time steps "
                               "and columns representing individual loads.")
        if prod_p.shape[0] != load_p.shape[0] or prod_p.shape[0] != load_q.shape[0]:
            raise RuntimeError(f"prod_p, load_p and load_q should have the same number of "
                               f"rows. We found: prod_p.shape[0] = {prod_p.shape[0]}, load_p.shape[0] = {load_p.shape[0]} "
                               f"and load_q.shape[0] = {load_q.shape[0]}")
        if prod_p.shape[1] != self.grid2op_env.n_gen:
            raise RuntimeError(f"The number of generators on the grid {self.grid2op_env.n_gen} "
                               f"is different that the number of columns of the provided prod_p data: "
                               f"prod_p.shape[1] = {prod_p.shape[1]}")
        if load_p.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid {self.grid2op_env.n_load} "
                               f"is different that the number of columns of the provided load_p data: "
                               f"load_p.shape[1] = {load_p.shape[1]}")
        if load_q.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid {self.grid2op_env.n_load} "
                               f"is different that the number of columns of the provided load_q data: "
                               f"load_q.shape[1] = {load_q.shape[1]}")
        if v_init is None:
            v_init_comp = self.grid2op_env.backend.V
        else:
            v_init_comp = 1.0 * v_init  # make a copy !
        status = self.computer.compute_Vs(prod_p,
                                          np.zeros((prod_p.shape[0], 0)),  # no static generators for now !
                                          load_p,
                                          load_q,
                                          v_init_comp,
                                          self.grid2op_env.backend.max_it,
                                          self.grid2op_env.backend.tol)
        if status != 1 and not ignore_errors:
            # raise an error if the powerflow diverged
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        elif status != 1:
            # only raise a warning in this case
            warnings.warn(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        Vs = 1.0 * self.computer.get_voltages()  # If I don't copy, lazy eval may break stuff... 
        # eg test_time_series_dc.py does behave stochastically
        self.__computed = True
        return Vs
        
    def compute_V(self, scenario_id=None, seed=None, v_init=None, ignore_errors=False):
        """
        This function allows to retrieve the complex voltage at each bus of the grid for each step.

        .. warning:: Topology fixed = no maintenance, no attacks, etc.

            As the topology is fixed, this class does not allow to simulate the effect of maintenance or attacks !
        """
        prod_p, load_p, load_q = self.get_injections(scenario_id=scenario_id, seed=seed)
        Vs = self.compute_V_from_inj(prod_p, load_p, load_q, v_init, ignore_errors)
        return Vs

    def compute_A(self):
        """
        This function returns the current flows (in Amps, A) at the origin (for powerline) / high voltage (for transformer) 
        side
        
        It does not recompute the voltages at each buses, it uses the information get from `compute_V` and
        This is why you must call `compute_V(...)` first !
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        ampss = self.computer.compute_flows()
        return 1000. * ampss

    def compute_P(self):
        """
        This function returns the active power flows (in MW) at the origin (for powerline) / high voltage (for transformer) 
        side
        
        It does not recompute the voltages at each buses, it uses the information get from `compute_V` and
        This is why you must call `compute_V(...)` first !
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        mws = 1. * self.computer.compute_power_flows() # If I don't copy, lazy eval may break stuff... 
        # eg test_time_series_dc.py does behave stochastically
        return mws

    def get_flows(self, scenario_id=None, seed=None, v_init=None, ignore_errors=False):
        """
        Retrieve the flows for each step simulated.

        Each row of the resulting flow matrix will correspond to a step.

        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import TimeSerie
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            timeserie = TimeSerie(env)
            res_p, res_a, res_v = timeserie.get_flows(scenario_id, seed, v_init, ignore_errors)

            # in this results, then
            # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
            # you can retrieve it with `security_analysis.contingency_order[row_id]`
        """
        
        Vs = self.compute_V(scenario_id, seed, v_init, ignore_errors)
        amps = self.compute_A()
        Ps = self.compute_P()
        return Ps, amps, Vs
    
    def clear(self):
        """
        Clear everything, as if nothing has been computed
        """
        self.computer.clear()
        self.__computed = False
        
        self.prod_p = None
        self.load_p = None
        self.load_q = None
    
    def close(self):
        """permanently close the object"""
        self.grid2op_env.close()
        self.clear()
        self.computer.close()
         
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
