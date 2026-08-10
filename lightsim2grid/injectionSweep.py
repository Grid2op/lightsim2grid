# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

__all__ = ["InjectionSweepCPP"]

from lightsim2grid.timeSerie import TimeSerie, GRID2OP_INSTALLED
from .lightsim2grid_cpp import InjectionSweepCPP

if GRID2OP_INSTALLED:
    __all__.append("InjectionSweep")


class InjectionSweep(TimeSerie):
    """
    Same computation as :class:`lightsim2grid.timeSerie.TimeSerie` -- a fixed grid topology,
    one powerflow per set of injections -- but every powerflow starts from the same voltage
    instead of from the result of the previous one.

    That makes the steps independent of one another: the result of a step does not depend on
    the steps computed before it, nor on the order in which they were given. Use this class
    when the "steps" are unrelated scenarios rather than consecutive instants of a time
    series. Two practical consequences:

    - the computation can be spread over several OS threads, see ``nb_thread`` below
      (:class:`lightsim2grid.timeSerie.TimeSerie` cannot: splitting a chained computation
      would make its results depend on how it was split);
    - a step that is far from its neighbours does not inherit a bad starting point from
      them -- but a step that IS close to its neighbours no longer benefits from their
      solution, so a genuine time series usually converges in fewer iterations with
      :class:`lightsim2grid.timeSerie.TimeSerie`.

    Examples
    ---------

    It is used exactly like :class:`lightsim2grid.timeSerie.TimeSerie`:

    .. code-block:: python

        import grid2op
        from lightsim2grid import InjectionSweep
        from lightsim2grid import LightSimBackend

        env_name = ...
        env = grid2op.make(env_name, backend=LightSimBackend())

        sweep = InjectionSweep(env)
        sweep.nb_thread = 4  # optional, the results do not depend on it
        res_p, res_a, res_v = sweep.get_flows(scenario_id=..., seed=...)

    """

    _CPP_CLASS = InjectionSweepCPP

    @property
    def nb_thread(self):
        """Number of OS threads used to compute the steps (default: ``1``).

        The steps are split into contiguous ranges, each solved by its own thread with its
        own solver, writing to disjoint rows of the result matrix: the results do **not**
        depend on the number of threads. Values ``< 1`` are clamped to ``1``.

        Must be set before the computation actually runs (eg before ``compute_V`` is
        called); it has no effect on a batch that has already been computed.
        """
        return self.computer.nb_thread

    @nb_thread.setter
    def nb_thread(self, val: int):
        if int(val) != val:
            raise ValueError("The `nb_thread` attribute must be an integer.")
        self.computer.nb_thread = int(val)

    @property
    def init_from_n_powerflow(self):
        """Whether to initialize the complex voltages of **each** step of the batch with the
        results of a "n" powerflow (a powerflow at the current state of the grid) instead of
        the vector given to ``compute_V``. Default: ``False``.

        Unlike :attr:`lightsim2grid.timeSerie.TimeSerie.init_from_n_powerflow`, this applies
        to every step and not only the first one -- here every step starts from that same
        voltage. Must be set before the computation actually runs.
        """
        return self.computer.init_from_n_powerflow

    @init_from_n_powerflow.setter
    def init_from_n_powerflow(self, val: bool):
        if bool(val) != val:
            raise ValueError("The `init_from_n_powerflow` attribute must be a boolean.")
        self.computer.init_from_n_powerflow = bool(val)
