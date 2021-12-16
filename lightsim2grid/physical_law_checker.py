# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import os
from lightsim2grid.lightSimBackend import LightSimBackend


class PhysicalLawChecker:
    """
    Utility tools to check that a given vector meets the KCL (Kirchhoff's Current Law)
    for a given grid (given as a grid2op observation).

    Notes
    ------
    Due to floating point precision between grid2op (float32) and lightsim2grid (float64) there might be a slight
    difference in the "mismatch" in the order of 1e-5 / 1e-6. If you want to check the kirchoffs law, then you won't be
    able to check with a tolerance less than 1e-5.

    .. warning::
        The grid2op environment is read from a grid.json file. Make sure to use an environment that can be
        loaded by lightsim2grid !

    Examples
    --------

    It can be used as:

    .. code-block:: python

        import grid2op
        import numpy as np
        from lightsim2grid import PhysicalLawChecker

        # create a grid2op environment
        env_name = "l2rpn_case14_sandbox"
        env = grid2op.make(env_name, ...)

        # create the checker
        checker = PhysicalLawChecker(env)

        # get an observation
        obs = env.reset()

        # retrieve somehow a complex voltage
        v = np.zeros(2*env.n_sub, dtype=complex)
        v[...] = ...  # put here the value of the complex voltage you want to get

        # check if it meets the KCL (Kirchhoff's Current Law)
        mismatch = checker.check_solution(v, obs)
        # mistmatch has same size as v and contains the (complex) current mismatch at each bus of the grid.


    """
    def __init__(self, grid2op_env):
        # lazy import to avoir circular import  TODO why ?
        from grid2op.Environment.Environment import Environment
        from grid2op.Environment.MultiMixEnv import MultiMixEnvironment
        if isinstance(grid2op_env, Environment):
            env_path = grid2op_env.get_path_env()
            gridobj_ref = type(grid2op_env)
        elif isinstance(grid2op_env, MultiMixEnvironment):
            env_path = grid2op_env.current_env.get_path_env()
            gridobj_ref = type(grid2op_env.current_env)
        else:
            raise RuntimeError("Only grid2op Environment and MultiMixEnvironment are supported at the moment.")
        grid_path = os.path.join(env_path, "grid.json")
        if not os.path.exists(grid_path):
            raise RuntimeError("Unable to locate the environment grid file. This feature will not work.")

        self.init_env = grid2op_env
        self._this_backend = LightSimBackend.init_grid(gridobj_ref)()
        self._this_backend.load_grid(grid_path)

    def check_solution(self, vcomplex, grid2op_obs, with_qlim=False):
        """
        This function checks that the complete `vcomplex` vector meets the kirchhoff current laws for a given grid.
        This vector is expected to be a complex vector representing the complex voltage at each bus of a given
        powergrid.

        The given grid is given in the form of a grid2op observation that should come from the same environment as the
        environment used to initialize the `PhysicalLawChecker`.

        Parameters
        ----------
        vcomplex:
            A numpy complex vector representing the complex voltage at each bus (should have the same size as the total number of buses in the 
            grid2op environment)

        grid2op_obs:
            A grid2op observation representing the state of the grid

        with_qlim:
            Whether to take into account the reactive limits for the generators

        Returns
        -------
        res:
            A numpy array (complex) that has the same shape as `vcomplex` and that computes the power mismatch
            at each bus of the grid.

        Notes
        -----
        The vector `vcomplex` needs to have as many components as the possible number of buses (typically `2 * env.n_sub`)
        for regular grid2op environments.

        If some buses are "empty" / "deactivated" (*ie* no elements are connected to it) then the given components of `vcomplex` is
        ignored: you can put whatever you want there it will have no impact on the final result.

        The order of the buses are expected to be given in the following order: First `n_sub` components of
        `vcomplex` (and `res`) represents bus 1 of the substation 1, 2, ..., n_sub and last `n_sub` components
        of `vcomplex` (and `res`) represents bus 2 of the same substation.
        """

        if vcomplex.shape[0] != 2 * self.init_env.n_sub:
            raise RuntimeError(f"The size of the input vector vcomplex needs to match the total number of possible "
                               f"buses on the grid (so basically 2 * env.n_sub). It is {vcomplex.shape[0]} while the "
                               f"total number of possible buses is {self.init_env.n_sub}. [NB: disconnected buses "
                               f"will not be taken into account anyway, you can put whatever you like there.]")
        self._this_backend.update_from_obs(grid2op_obs)
        res = self._this_backend._grid.check_solution(vcomplex, with_qlim)
        return res
