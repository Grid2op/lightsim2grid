# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
A grid2op environment whose powergrid is a MATPOWER case (``grid.m``).

Build it with the matpower loader, the same way a pypowsybl / iidm environment is
built with the pypowsybl one -- the backend has to be given explicitly, because the
`backend` entry of this config is instantiated by grid2op with no argument and would
otherwise default to reading a `grid.json`:

.. code-block:: python

    import grid2op
    from lightsim2grid import LightSimBackend

    env = grid2op.make(path_to_this_directory,
                       backend=LightSimBackend(loader_method="matpower"))

MATPOWER carries a single operating point and no time series, so the environment
"plays" that operating point over and over (``ChangeNothing``) -- enough to exercise
the backend through a real environment, a Runner and a copy. Same reason for the
thermal limits below: MATPOWER's RATE_A column is a branch MVA rating, and 0
("unlimited") in this case, so they are declared here, which is where a grid2op
environment declares them anyway. They are the base-case flows with a 50% margin, so
the environment starts around rho = 0.67 rather than at rho = 0.
"""

from grid2op.Action import TopologyAndDispatchAction
from grid2op.Reward import L2RPNReward
from grid2op.Rules import DefaultRules
from grid2op.Chronics import ChangeNothing


config = {
    "backend": None,  # see the note above: pass `backend=LightSimBackend(loader_method="matpower")`
    "action_class": TopologyAndDispatchAction,
    "observation_class": None,
    "reward_class": L2RPNReward,
    "gamerules_class": DefaultRules,
    "chronics_class": ChangeNothing,
    "volagecontroler_class": None,
    "names_chronics_to_grid": None,
    # the 17 powerlines then the 3 transformers, in the order the matpower `branch`
    # table lists them (a plain "ratio == 0" branch is a powerline, a branch with a
    # ratio is a transformer, and lightsim2grid puts all the lines before all the
    # transformers)
    "thermal_limits": [
        979.0, 468.0, 461.0, 353.0, 262.0, 156.0, 407.0, 200.0, 201.0, 470.0,
        167.0, 251.0, 103.0, 44.0, 148.0, 1273.0, 709.0, 749.0, 402.0, 1099.0,
    ],
}
