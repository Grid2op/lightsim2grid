# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["ScenarioSweepCPP"]

import numpy as np
import warnings


try:
    from lightsim2grid.lightSimBackend import LightSimBackend
    __all__.append("ScenarioSweep")
    GRID2OP_INSTALLED = True
except ImportError as exc_:  # noqa: F841
    # grid2Op is not installed
    GRID2OP_INSTALLED = False

from lightsim2grid.algorithm import AlgorithmType
from .lightsim2grid_cpp import ScenarioSweepCPP


class ScenarioSweep:
    """
    Batch powerflow that varies **both** the injection and a contingency (line / trafo
    disconnection) per simulation, independently, row-aligned: row `i` of every
    ``modify_*`` input is solved together with row `i` of ``set_contingency_lines`` /
    ``set_contingency_trafos``.

    Unlike :class:`lightsim2grid.timeSerie.TimeSerie` /
    :class:`lightsim2grid.injectionSweep.InjectionSweep` (a single bundled
    ``compute_V_from_inj`` call), this class uses a setter-based API: call as many of
    ``modify_gen_p`` / ``modify_sgen_p`` / ``modify_load_p`` / ``modify_load_q`` /
    ``set_contingency_lines`` / ``set_contingency_trafos`` as are relevant, then
    ``compute()``. Any axis you never set defaults to the grid's own state for every
    row (its own target injection, or "nothing disconnected"). The first setter you
    call fixes the number of simulations; every later setter is checked against it
    immediately.

    This is deliberately a *different* API from
    :class:`lightsim2grid.contingencyAnalysis.ContingencyAnalysis`'s ``add_n1`` /
    ``add_nk`` (a set of distinct scenarios applied to one shared base case): here
    each row pairs its own injection with its own contingency, so a dense,
    row-aligned mask is the natural shape -- the two APIs are not interchangeable.

    Examples
    ---------

    .. code-block:: python

        import numpy as np
        import grid2op
        from lightsim2grid import ScenarioSweep
        from lightsim2grid import LightSimBackend

        env_name = ...
        env = grid2op.make(env_name, backend=LightSimBackend())

        sweep = ScenarioSweep(env)
        sweep.modify_load_p(load_p)  # shape (n_simul, n_load)
        line_mask = np.zeros((load_p.shape[0], env.n_line), dtype=bool)
        line_mask[:, 3] = True  # disconnect line 3 for every simulation
        sweep.set_contingency_lines(line_mask)
        sweep.compute()
        Vs = sweep.get_voltages()
        Ps, amps = sweep.compute_P(), sweep.compute_A()
    """

    #: the c++ class this wrapper drives.
    _CPP_CLASS = ScenarioSweepCPP

    def __init__(self, grid2op_env):
        if not GRID2OP_INSTALLED:
            raise RuntimeError(f"Impossible to use the python wrapper `{type(self).__name__}` "
                               "when grid2op is not installed. Please fall back to the "
                               "c++ version (available in python) with:\n"
                               f"\tfrom lightsim2grid.scenarioSweep import {type(self)._CPP_CLASS.__name__}\n"
                               "and refer to the appropriate documentation.")

        from grid2op.Environment import Environment  # type: ignore # otherwise i got issues...
        if not isinstance(grid2op_env.backend, LightSimBackend):
            raise RuntimeError("This class only works with LightSimBackend")
        if not isinstance(grid2op_env, Environment):
            raise RuntimeError("Please an environment of class \"Environment\", "
                               "and not \"MultimixEnv\" or \"BaseMultiProcessEnv\"")
        self.grid2op_env = grid2op_env.copy()
        self.computer = type(self)._CPP_CLASS(self.grid2op_env.backend._grid)
        self.__computed = False

        self.available_default_algorithms = self.computer.available_default_algorithms()
        if AlgorithmType.NR_KLU in self.available_default_algorithms:
            # use the faster KLU if available
            self.computer.change_algorithm(AlgorithmType.NR_KLU)

    @property
    def nb_thread(self):
        """Number of OS threads used to compute the steps (default: ``1``).

        The steps are split into contiguous ranges, each solved by its own thread with its
        own solver (and its own admittance-matrix copy, since a contingency edits it), writing
        to disjoint rows of the result matrix: the results do **not** depend on the number of
        threads. Values ``< 1`` are clamped to ``1``.

        Must be set before ``compute()`` actually runs; it has no effect on a batch that has
        already been computed.
        """
        return self.computer.nb_thread

    @nb_thread.setter
    def nb_thread(self, val: int):
        if int(val) != val:
            raise ValueError("The `nb_thread` attribute must be an integer.")
        self.computer.nb_thread = int(val)

    @property
    def init_from_n_powerflow(self):
        """Whether to initialize the complex voltages of **each** simulation with the results
        of a "n" powerflow (a powerflow with no injection change and no contingency) instead of
        the vector given to ``compute``. Default: ``False``. Must be set before ``compute()``
        actually runs.
        """
        return self.computer.init_from_n_powerflow

    @init_from_n_powerflow.setter
    def init_from_n_powerflow(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `init_from_n_powerflow` attribute must be a boolean.")
        self.computer.init_from_n_powerflow = bool(val)

    def _check_2d(self, arr, name):
        arr = np.asarray(arr)
        if arr.ndim != 2:
            raise RuntimeError(f"{name} should be a matrix with rows representing simulations "
                               f"and columns representing individual elements.")
        return arr

    def modify_gen_p(self, gen_p):
        """Per-step active generator setpoints, shape ``(n_simul, n_gen)``."""
        gen_p = self._check_2d(gen_p, "gen_p")
        if gen_p.shape[1] != self.grid2op_env.n_gen:
            raise RuntimeError(f"The number of generators on the grid ({self.grid2op_env.n_gen}) "
                               f"differs from the number of columns of gen_p ({gen_p.shape[1]}).")
        self.computer.modify_gen_p(gen_p)
        self.__computed = False

    def modify_sgen_p(self, sgen_p):
        """Per-step active static generator setpoints, shape ``(n_simul, n_sgen)``."""
        sgen_p = self._check_2d(sgen_p, "sgen_p")
        self.computer.modify_sgen_p(sgen_p)
        self.__computed = False

    def modify_load_p(self, load_p):
        """Per-step active load setpoints, shape ``(n_simul, n_load)``."""
        load_p = self._check_2d(load_p, "load_p")
        if load_p.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid ({self.grid2op_env.n_load}) "
                               f"differs from the number of columns of load_p ({load_p.shape[1]}).")
        self.computer.modify_load_p(load_p)
        self.__computed = False

    def modify_load_q(self, load_q):
        """Per-step reactive load setpoints, shape ``(n_simul, n_load)``."""
        load_q = self._check_2d(load_q, "load_q")
        if load_q.shape[1] != self.grid2op_env.n_load:
            raise RuntimeError(f"The number of loads on the grid ({self.grid2op_env.n_load}) "
                               f"differs from the number of columns of load_q ({load_q.shape[1]}).")
        self.computer.modify_load_q(load_q)
        self.__computed = False

    def set_contingency_lines(self, mask):
        """Per-step powerline contingency mask, shape ``(n_simul, n_line)``, dtype bool.
        ``True`` means "deactivate this powerline for this simulation"."""
        mask = self._check_2d(mask, "mask").astype(dtype=bool, copy=False)
        n_line = len(self.grid2op_env.backend._grid.get_lines())
        if mask.shape[1] != n_line:
            raise RuntimeError(f"The number of powerlines on the grid ({n_line}) "
                               f"differs from the number of columns of the mask ({mask.shape[1]}).")
        self.computer.set_contingency_lines(mask)
        self.__computed = False

    def set_contingency_trafos(self, mask):
        """Per-step trafo contingency mask, shape ``(n_simul, n_trafo)``, dtype bool.
        See :func:`set_contingency_lines`."""
        mask = self._check_2d(mask, "mask").astype(dtype=bool, copy=False)
        n_trafo = len(self.grid2op_env.backend._grid.get_trafos())
        if mask.shape[1] != n_trafo:
            raise RuntimeError(f"The number of trafos on the grid ({n_trafo}) "
                               f"differs from the number of columns of the mask ({mask.shape[1]}).")
        self.computer.set_contingency_trafos(mask)
        self.__computed = False

    def compute(self, v_init=None, max_iter=None, tol=None, ignore_errors=False):
        """
        Run the batch: one powerflow per simulation, using whatever was set by the
        ``modify_*`` / ``set_contingency_*`` setters. Raises if nothing was ever set.

        ``max_iter`` / ``tol`` default to the backend's own values
        (``self.grid2op_env.backend.max_it`` / ``.tol``) when not given.
        """
        if v_init is None:
            v_init_comp = self.grid2op_env.backend.V
        else:
            v_init_comp = 1.0 * v_init  # make a copy !
        if max_iter is None:
            max_iter = self.grid2op_env.backend.max_it
        if tol is None:
            tol = self.grid2op_env.backend.tol
        self.computer.compute(v_init_comp, max_iter, tol)
        status = self.computer.get_status()
        if status != 1 and not ignore_errors:
            raise RuntimeError(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        elif status != 1:
            warnings.warn(f"Some error occurred, the powerflow has diverged after {self.computer.nb_solved()} step(s)")
        self.__computed = True

    def get_voltages(self):
        """Complex voltage, at each bus, for each simulation. Must be called after ``compute()``."""
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        return self.computer.get_voltages().copy()

    def compute_A(self):
        """
        Current flows (in Amps, A) at the origin (for powerline) / high voltage (for
        transformer) side, per simulation. Does not recompute the voltages; you must
        call ``compute()`` first.
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        ampss = self.computer.compute_flows()
        return 1000. * ampss

    def compute_P(self):
        """
        Active power flows (in MW) at the origin (for powerline) / high voltage (for
        transformer) side, per simulation. Does not recompute the voltages; you must
        call ``compute()`` first.
        """
        if not self.__computed:
            raise RuntimeError("This function can only be used if compute() has been successfully called")
        mws = 1. * self.computer.compute_power_flows()  # If I don't copy, lazy eval may break stuff...
        return mws

    def get_flows(self, v_init=None, max_iter=None, tol=None, ignore_errors=False):
        """
        Run ``compute()`` (using whatever was set by the setters) then retrieve the
        active power / current / voltage results for every simulation, in one call.

        Each row of the resulting flow matrix corresponds to a simulation.
        """
        self.compute(v_init, max_iter, tol, ignore_errors)
        amps = self.compute_A()
        Ps = self.compute_P()
        Vs = self.get_voltages()
        return Ps, amps, Vs

    def clear(self):
        """Clear everything, as if nothing had ever been set / computed."""
        self.computer.clear()
        self.__computed = False

    def close(self):
        """permanently close the object"""
        self.grid2op_env.close()
        self.clear()
        self.computer.close()
