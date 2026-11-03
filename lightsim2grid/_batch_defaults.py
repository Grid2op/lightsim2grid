# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The algorithm a batch wrapper (TimeSerie, InjectionSweep, ScenarioSweep,
ContingencyAnalysis) starts with."""

from lightsim2grid.algorithm import AlgorithmType


def use_default_batch_algorithm(computer, grid):
    """A batch built from `grid` inherits its AC algorithm. One that runs outer loops
    (NROuter_*) is kept: every row is then solved with the grid's outer loops, as the backend
    solves the grid. Any other is replaced by the fastest Newton-Raphson available (NR_KLU)."""
    if grid.get_algo().supports_outer_loops():
        return
    if AlgorithmType.NR_KLU in computer.available_default_algorithms():
        computer.change_algorithm(AlgorithmType.NR_KLU)
