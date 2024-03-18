# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["ContingencyAnalysisCPP", "ContingencyAnalysis",
           # deprecated
           "SecurityAnalysisCPP", "SecurityAnalysis",
           ]

import copy
import numpy as np
from collections.abc import Iterable

from lightsim2grid.lightSimBackend import LightSimBackend
from lightsim2grid.solver import SolverType
from lightsim2grid_cpp import ContingencyAnalysisCPP


class ContingencyAnalysis(object):
    """
    This class allows to perform a "security analysis" from a given grid state.

    For now, you cannot change the grid state, and it only computes the security analysis with 
    current flows at origin of powerlines.
    
    Feel free to post a feature request if you want to extend it.

    This class is used in 4 phases:
    
    0) you create it from a grid2op environment (the grid topology will not be modified from this environment)
    1) you add some contingencies to simulate
    2) you start the simulation
    3) you read back the results
    
    
    Examples
    --------
    An example is given here
    
    .. code-block:: python

        import grid2op
        from lightsim2grid import SecurityAnalysis
        from lightsim2grid import LightSimBackend
        env_name = ...
        env = grid2op.make(env_name, backend=LightSimBackend())

        0) you create
        security_analysis = SecurityAnalysis(env)
        
        1) you add some contingencies to simulate
        security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
        
        2) you start the simulation (done automatically)
        3) you read back the results
        res_p, res_a, res_v = security_analysis.get_flows()

        # in this results, then
        # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
        # you can retrieve it with `security_analysis.contingency_order[row_id]`

    Notes
    ------

    Sometimes, the behaviour might differ from grid2op. For example, if simulating a contingency
    leads to a non connected grid, then this function will return "Nan" for the flows and 0. for
    the voltages.

    In grid2op, it would be, in this case, 0. for the flows and 0. for the voltages.

    """
    STR_TYPES = (str, np.str_)  # np.str deprecated in numpy 1.20 and earlier versions not supported anyway
        
    def __init__(self, grid2op_env):
        if not isinstance(grid2op_env.backend, LightSimBackend):
            raise RuntimeError("This class only works with LightSimBackend")
        self.grid2op_env = grid2op_env.copy()
        self.computer = ContingencyAnalysisCPP(self.grid2op_env.backend._grid)
        self._contingency_order = {}  # key: contingency (as tuple), value: order in which it is entered
        self._all_contingencies = []
        self.__computed = False
        self._vs = None
        self._ampss = None

        self.available_solvers = self.computer.available_solvers()
        if SolverType.KLU in self.available_solvers:
            # use the faster KLU if available
            self.computer.change_solver(SolverType.KLU)

    @property
    def all_contingencies(self):
        return copy.deepcopy(self._all_contingencies)

    @all_contingencies.setter
    def all_contingencies(self, val):
        raise RuntimeError("Impossible to add new topologies like this. Please use `add_single_contingency` "
                           "or `add_multiple_contingencies`.")

    def clear(self):
        """
        Clear the list of contingencies to simulate
        """
        self.computer.clear()
        self._contingency_order = {}
        self.__computed = False
        self._vs = None
        self._ampss = None
        self._all_contingencies = []

    def _single_cont_to_li_int(self, single_cont):
        li_disc = []
        if isinstance(single_cont, int):
            single_cont = [single_cont]

        for stuff in single_cont:
            if isinstance(stuff, type(self).STR_TYPES):
                stuff = np.where(self.grid2op_env.name_line == stuff)
                stuff = stuff[0]
                if stuff.size == 0:
                    # name is not found
                    raise RuntimeError(f"Impossible to find a powerline named \"{stuff}\" in the environment")
                stuff = int(stuff[0])
            else:
                stuff = int(stuff)
            li_disc.append(stuff)
        return li_disc

    def add_single_contingency(self, *args):
        """
        This function allows to add a single contingency specified by either the powerlines names
        (which should match env.name_line) or by their ID.

        The contingency added can be a "n-1" which will simulate a single powerline disconnection
        or a "n-k" which will simulate the disconnection of multiple powerlines.

        It does not accept any keword arguments.

        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import SecurityAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_anlysis = SecurityAnalysis(env)
            # the single (n-1) contingency "disconnect powerline 0" is added
            security_anlysis.add_single_contingency(0)

            # add the single (n-1) contingency "disconnect line 1
            security_anlysis.add_single_contingency(env.name_line[1])

            # add a single contingency that disconnect powerline 2 and 3 at the same time
            security_anlysis.add_single_contingency(env.name_line[2], 3)

        Notes
        -----
        If it raises an error for a given contingency, the object might be not properly initialized.
        In this case, we recommend you to clear it (using the `clear()` method and to attempt to 
        add contingencies again.)

        """
        li_disc = self._single_cont_to_li_int(args)
        li_disc_tup = tuple(li_disc)
        if li_disc_tup not in self._contingency_order:
            # this is really the first time this contingency is seen
            try:
                self.computer.add_nk(li_disc)
                my_id = len(self._contingency_order)
                self._contingency_order[li_disc_tup] = my_id
                self._all_contingencies.append(li_disc_tup)
            except Exception as exc_:
                raise RuntimeError(f"Impossible to add the contingency {args}. The most likely cause "
                                   f"is that you try to disconnect a powerline that is not present "
                                   f"on the grid") from exc_

    def add_multiple_contingencies(self, *args):
        """
        This function will add multiple contingencies at the same time.

        This code is equivalent to:

        .. code-block:: python

            for single_cont in args:
                self.add_single_contingency(single_cont)

        It does not accept any keword arguments.

        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import SecurityAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_anlysis = SecurityAnalysis(env)

            # add a single contingency that disconnect powerline 2 and 3 at the same time
            security_anlysis.add_single_contingency(env.name_line[2], 3)

            # add a multiple contingencies the first one disconnect powerline 2 and 
            # and the second one disconnect powerline 3
            security_anlysis.add_multiple_contingencies(env.name_line[2], 3)
        """        
        for single_cont in args:
            if isinstance(single_cont, Iterable) and not isinstance(single_cont, type(self).STR_TYPES):
                # this is a contingency consisting in cutting multiple powerlines
                self.add_single_contingency(*single_cont)
            else:
                # this is likely an int or a string representing a contingency
                self.add_single_contingency(single_cont)

    def add_all_n1_contingencies(self):
        """
        This method registers as the contingencies that will be computed all the contingencies that disconnects 1 powerline

        This is equivalent to:

        .. code-block:: python

            for single_cont_id in range(env.n_line):
                self.add_single_contingency(single_cont_id)
        """
        for single_cont_id in range(self.grid2op_env.n_line):
            self.add_single_contingency(single_cont_id)

    def get_flows(self, *args):
        """
        Retrieve the flows after each contingencies has been simulated.

        Each row of the resulting flow matrix will correspond to a contingency simulated in the arguments.

        You can require only the result on some contingencies with the `args` argument, but in each case, all the results will
        be computed. If you don't specify anything, the results will be returned for all contingencies (which we recommend to do)
        
        Examples
        --------

        .. code-block:: python

            import grid2op
            from lightsim2grid import SecurityAnalysis
            from lightsim2grid import LightSimBackend
            env_name = ...
            env = grid2op.make(env_name, backend=LightSimBackend())

            security_analysis = SecurityAnalysis(env)
            security_analysis.add_multiple_contingencies(...) # or security_analysis.add_single_contingency(...)
            res_p, res_a, res_v = security_analysis.get_flows()

            # in this results, then
            # res_a[row_id] will be the flows, on all powerline corresponding to the `row_id` contingency.
            # you can retrieve it with `security_analysis.contingency_order[row_id]`
        """
        
        all_defaults = self.computer.my_defaults()
        if len(args) == 0:
            # default: i consider all contingencies
            orders_ = np.zeros(len(all_defaults), dtype=int)
            for id_cpp, cont_ in enumerate(all_defaults):
                tup_ = tuple(cont_)
                orders_[self._contingency_order[tup_]] = id_cpp
        else:
            # a list of interesting contingencies has been provided
            orders_ = np.zeros(len(args), dtype=int)
            all_defaults = [tuple(cont) for cont in all_defaults]
            for id_me, cont_ in enumerate(args):
                cont_li = self._single_cont_to_li_int(cont_)
                tup_ = tuple(cont_li)
                if tup_ not in self._contingency_order:
                    raise RuntimeError(f"Contingency {cont_} is not simulated by this class. Have you called "
                                       f"`add_single_contingency` or `add_multiple_contingencies` ?")
                id_cpp = all_defaults.index(tup_)
                orders_[id_me] = id_cpp

        if not self.__computed:
            self.compute_V()
            self.compute_A()
            self.compute_P()
        
        return self._mws[orders_], self._ampss[orders_], self._vs[orders_]

    def compute_V(self):
        """
        This function allows to retrieve the complex voltage at each bus of the grid for each contingency.

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results

        """
        v_init = self.grid2op_env.backend.V
        self.computer.compute(v_init,
                              self.grid2op_env.backend.max_it,
                              self.grid2op_env.backend.tol)
        self._vs = self.computer.get_voltages()
        self.__computed = True
        return self._vs

    def compute_A(self):
        """
        This function returns the current flows (in Amps, A) at the origin / high voltage side

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results !

        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        self._ampss = 1e3 * self.computer.compute_flows()
        return self._ampss

    def compute_P(self):
        """
        This function returns the active power flows (in MW) at the origin / high voltage side

        .. warning:: Order of the results

            The order in which the results are returned is NOT necessarily the order in which the contingencies have
            been entered. Please use `get_flows()` method for easier reading back of the results !

        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute_V has been sucessfully called")
        self._mws = 1.0 * self.computer.compute_power_flows()
        return self._mws

    def close(self):
        """permanently close the object"""
        self.grid2op_env.close()
        self.clear()
        self.computer.close()
