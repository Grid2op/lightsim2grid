# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Checks of the data the outer loops of the ``NROuter_*`` algorithms read on a grid."""

import warnings

__all__ = ["distributed_slack_is_idle", "warn_if_distributed_slack_idle"]


def _slack_bus(grid):
    """The bus the single-slack Newton of an ``NROuter_*`` algorithm keeps as its slack: the
    first (lowest id) bus of a connected slack participant, ``None`` without one."""
    buses = [el.bus_id for el in list(grid.get_generators()) + list(grid.get_storages())
             if el.is_slack and el.connected and el.bus_id >= 0]
    return min(buses) if buses else None


def distributed_slack_is_idle(grid) -> bool:
    """
    Whether the ``DistributedSlack`` outer loop is in the outer loops of ``grid`` but has
    nothing to act on: no connected unit flagged "can participate in the slack" stands off
    the slack bus, so the slack bus takes all of the active power imbalance. The loop is then
    dropped at run time, which says nothing on its own.

    This is what a single-slack grid looks like (one ``slack=True`` generator or one ext_grid
    of a pandapower grid): see ``init_from_pandapower(slack_weight_can_participate=True)``
    or ``LSGrid.set_gen_can_participate_slack`` to give it units to share on.
    """
    from lightsim2grid.algorithm import DistributedSlack
    if not any(isinstance(loop, DistributedSlack) for loop in grid.get_outer_loops()):
        return False
    slack_bus = _slack_bus(grid)
    if slack_bus is None:
        return False
    for el in list(grid.get_generators()) + list(grid.get_storages()):
        if el.connected and el.bus_id >= 0 and el.can_participate_slack_weight > 0. and el.bus_id != slack_bus:
            return False
    return True


def warn_if_distributed_slack_idle(grid) -> bool:
    """Warn when ``grid`` runs an ``NROuter_*`` algorithm and :func:`distributed_slack_is_idle`.
    Returns whether it warned."""
    if not str(grid.get_algo_name()).startswith("NROuter"):
        return False
    if not distributed_slack_is_idle(grid):
        return False
    warnings.warn("The DistributedSlack outer loop has no unit to share the active power imbalance on "
                  "(no unit that can participate in the slack stands off the slack bus): it will not run, "
                  "and the slack bus takes all of the imbalance. With a pandapower grid, see "
                  "`init_from_pandapower(..., slack_weight_can_participate=True)` (or "
                  "`LightSimBackend(loader_kwargs={'slack_weight_can_participate': True})`); otherwise "
                  "`LSGrid.set_gen_can_participate_slack`.",
                  UserWarning)
    return True
