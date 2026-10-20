# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Reverse-mode differentiation of powerflows, for pytorch: a batch
(:class:`BatchCPUPowerFlow`) or a single grid (:class:`GridPowerFlow`).

Not imported by ``lightsim2grid`` itself: torch is an optional dependency, so this
subpackage is reached explicitly.

    from lightsim2grid.differentiable import BatchCPUPowerFlow, GridPowerFlow
"""

__all__ = ["BatchCPUPowerFlow", "GridPowerFlow", "BranchFlows", "compute_branch_flows"]

from ._flows import BranchFlows, compute_branch_flows
from ._batch_cpu_power_flow import BatchCPUPowerFlow
from ._grid_power_flow import GridPowerFlow
