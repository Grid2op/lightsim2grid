# Copyright (c) 2025-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Python side of the light environment: a :class:`LightEnv` whose actions can be given
as grid2op actions."""

__all__ = ["LightEnv", "TopoAction", "Protections", "ElementType", "topo_action_from_grid2op"]

from typing import Iterable

import numpy as np

from grid2op.Action import BaseAction

from lightsim2grid.lightsim2grid_cpp import (LightEnv as _LightEnvCPP,  # pyright: ignore[reportMissingImports]
                                             TopoAction,
                                             Protections,
                                             ElementType)


# what a grid2op action can carry that the light environment does not model
_UNSUPPORTED_MODIF = (
    ("_modif_change_bus", "change_bus"),
    ("_modif_change_status", "change_line_status"),
    ("_modif_redispatch", "redispatch"),
    ("_modif_storage", "set_storage"),
    ("_modif_curtailment", "curtail"),
    ("_modif_inj", "injection"),
    ("_modif_shunt", "shunt"),
    ("_modif_detach_load", "detach_load"),
    ("_modif_detach_gen", "detach_gen"),
    ("_modif_detach_storage", "detach_storage"),
)


def topo_action_from_grid2op(action: BaseAction) -> TopoAction:
    """Convert a grid2op action into a (not yet checked) :class:`TopoAction`.

    Only ``set_bus`` and ``set_line_status`` are modelled by the light environment. An
    action using anything else that affects the grid (``change_bus``, ``change_line_status``,
    ``redispatch``, ``curtail``, storage, injections, shunts, detachment) raises a
    ``ValueError``. Alarms and alerts are ignored, they do not affect the grid.

    The values are copied as they are: the check that every element exists and every busbar
    exists is done against the grid by :func:`TopoAction.check_validity`, which
    :func:`LightEnv.init_actions` calls.
    """
    if not isinstance(action, BaseAction):
        raise TypeError(f"topo_action_from_grid2op expects a grid2op action, you provided a {type(action)}")
    unsupported = [name for attr, name in _UNSUPPORTED_MODIF if getattr(action, attr, False)]
    if unsupported:
        raise ValueError("The light environment only supports set_bus and set_line_status actions, "
                         f"this action also modifies: {', '.join(unsupported)}")
    cls = type(action)
    res = TopoAction()
    if action._modif_set_bus:
        set_bus = action._set_topo_vect
        for el_type, pos_topo_vect in ((ElementType.load, cls.load_pos_topo_vect),
                                       (ElementType.gen, cls.gen_pos_topo_vect),
                                       (ElementType.line_or, cls.line_or_pos_topo_vect),
                                       (ElementType.line_ex, cls.line_ex_pos_topo_vect),
                                       (ElementType.storage, cls.storage_pos_topo_vect)):
            for el_id, pos in enumerate(pos_topo_vect):
                bus = int(set_bus[pos])
                if bus != 0:
                    res.add_element(el_type, el_id, bus)
    if action._modif_set_status:
        for line_id, status in enumerate(action._set_line_status):
            status = int(status)
            if status != 0:
                res.set_line_status(line_id, status)
    return res


class LightEnv(_LightEnvCPP):
    """The light environment, see the c++ documentation. On top of it,
    :func:`LightEnv.init_actions` also accepts grid2op actions."""

    def init_actions(self, actions: Iterable) -> None:
        """Register the actions the agent can take: ``step(i)`` plays ``actions[i]``.

        ``actions`` is a list of grid2op actions (``set_bus`` and ``set_line_status`` only,
        see :func:`topo_action_from_grid2op`) or of :class:`TopoAction`, mixing both is
        fine. Every action is checked against the initial grid: an element that does not
        exist, a busbar that does not exist (eg ``-2``, or ``3`` on a grid with 2 busbars per
        substation), an unsupported modification or a contradiction raises a ``ValueError``
        naming the action, and nothing is registered.
        """
        converted = []
        for i, act in enumerate(actions):
            if isinstance(act, TopoAction):
                converted.append(act)
            elif isinstance(act, BaseAction):
                try:
                    converted.append(topo_action_from_grid2op(act))
                except (ValueError, TypeError) as exc_:
                    raise ValueError(f"LightEnv.init_actions: action {i} is invalid: {exc_}") from exc_
            else:
                raise ValueError(f"LightEnv.init_actions: action {i} is invalid: expected a grid2op action "
                                 f"or a TopoAction, got a {type(act)}")
        super().init_actions(converted)
