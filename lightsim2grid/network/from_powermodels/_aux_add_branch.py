# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np

from ._bus_remap import pm_bus_to_ls


def is_transformer(branch: dict) -> bool:
    """A PowerModels branch is a transformer if its own `"transformer"` flag says so,
    or its `"tap"` / `"shift"` are away from the "plain line" neutral (`tap == 1.0`,
    `shift == 0.0`) -- PowerModels normalizes a matpower `TAP == 0` (its own "this is a
    plain line" sentinel) to `tap = 1.0`, so `tap != 0` (matpower's own test) does not
    apply here.
    """
    return bool(branch.get("transformer", False)) or branch.get("tap", 1.0) != 1.0 or branch.get("shift", 0.0) != 0.0


def classify_branches(network: dict):
    """Returns `(line_keys, trafo_keys)`, each a list of `network["branch"]` string
    keys, sorted by `int(key)` (this order is a documented, deterministic contract:
    it is what `_aux_add_branch` below feeds `init_powerlines_full`/`init_trafo` in,
    and callers -- e.g. a validation helper matching lightsim2grid's own solved
    results back to the original PowerModels branch keys -- can reproduce it
    independently from `network` alone).
    """
    line_keys, trafo_keys = [], []
    for k in sorted(network["branch"], key=int):
        (trafo_keys if is_transformer(network["branch"][k]) else line_keys).append(k)
    return line_keys, trafo_keys


def _aux_add_branch(model, network, pm_to_ls, isolated_ls_bus):
    """
    Add the lines and transformers of `network["branch"]` into the lightsim2grid "model".

    Note
    ----
    Unlike lightsim2grid's `from_matpower` (which has to route through matpower's raw
    `mpc.branch` matrix, whose single `BR_B` column cannot represent asymmetric line
    charging), powerlines here keep `g_fr`/`b_fr` and `g_to`/`b_to` exactly as given by
    PowerModels, independently per side (`LineContainer`'s two-sided `init` overload
    supports this directly). Transformers are still limited to a single, symmetrically
    split charging admittance -- `TrafoContainer` has no "asymmetric" overload -- so an
    asymmetric transformer's `g_fr`/`g_to` (or `b_fr`/`b_to`) is summed and re-split
    50/50, with a warning (this never loses information for MATPOWER-derived data,
    where PowerModels always sets `b_fr == b_to` and `g_fr == g_to == 0`).

    Parameters
    ----------
    model
    network: dict
        The PowerModels network data dictionary
    pm_to_ls: dict
        PowerModels bus number (`"bus_i"`) -> lightsim2grid bus id
    isolated_ls_bus: numpy array
        lightsim2grid bus ids of isolated (`bus_type == 4`) buses

    """
    branch = network["branch"]
    line_keys, trafo_keys = classify_branches(network)

    f_bus_line = pm_bus_to_ls(np.array([int(branch[k]["f_bus"]) for k in line_keys]), pm_to_ls)
    t_bus_line = pm_bus_to_ls(np.array([int(branch[k]["t_bus"]) for k in line_keys]), pm_to_ls)
    r = np.array([branch[k]["br_r"] for k in line_keys])
    x = np.array([branch[k]["br_x"] for k in line_keys])
    h_or = np.array([complex(branch[k].get("g_fr", 0.), branch[k].get("b_fr", 0.)) for k in line_keys])
    h_ex = np.array([complex(branch[k].get("g_to", 0.), branch[k].get("b_to", 0.)) for k in line_keys])
    model.init_powerlines_full(r, x, h_or, h_ex, f_bus_line, t_bus_line)

    line_status = np.array([branch[k].get("br_status", 1) for k in line_keys]) != 0
    if isolated_ls_bus.size:
        line_status &= ~np.isin(f_bus_line, isolated_ls_bus)
        line_status &= ~np.isin(t_bus_line, isolated_ls_bus)
    for line_id, is_ok in enumerate(line_status):
        if not is_ok:
            model.deactivate_powerline(line_id)

    f_bus_trafo = pm_bus_to_ls(np.array([int(branch[k]["f_bus"]) for k in trafo_keys]), pm_to_ls)
    t_bus_trafo = pm_bus_to_ls(np.array([int(branch[k]["t_bus"]) for k in trafo_keys]), pm_to_ls)
    trafo_r = np.array([branch[k]["br_r"] for k in trafo_keys])
    trafo_x = np.array([branch[k]["br_x"] for k in trafo_keys])
    asymmetric = [k for k in trafo_keys
                 if branch[k].get("g_fr", 0.) != branch[k].get("g_to", 0.)
                 or branch[k].get("b_fr", 0.) != branch[k].get("b_to", 0.)]
    if asymmetric:
        warnings.warn(f"{len(asymmetric)} transformer(s) have an asymmetric charging "
                      f"admittance (\"g_fr\"/\"b_fr\" != \"g_to\"/\"b_to\"), which "
                      "lightsim2grid's transformer model cannot represent (only a single, "
                      "symmetrically-split value): the total charging admittance is kept, "
                      "but re-split 50/50 between the two sides.")
    trafo_b = np.array([complex(branch[k].get("g_fr", 0.) + branch[k].get("g_to", 0.),
                                branch[k].get("b_fr", 0.) + branch[k].get("b_to", 0.))
                       for k in trafo_keys])
    trafo_ratio = np.array([branch[k].get("tap", 1.0) for k in trafo_keys])
    # PowerModels stores angles in radians; `init_trafo` expects degrees
    trafo_shift_degree = np.rad2deg([branch[k].get("shift", 0.0) for k in trafo_keys])
    trafo_tap_hv = [True] * len(trafo_keys)  # PowerModels always applies tap/shift on the "from" side
    model.init_trafo(trafo_r, trafo_x, trafo_b, trafo_ratio, trafo_shift_degree, trafo_tap_hv,
                     f_bus_trafo, t_bus_trafo, True)

    trafo_status = np.array([branch[k].get("br_status", 1) for k in trafo_keys]) != 0
    if isolated_ls_bus.size:
        trafo_status &= ~np.isin(f_bus_trafo, isolated_ls_bus)
        trafo_status &= ~np.isin(t_bus_trafo, isolated_ls_bus)
    for trafo_id, is_ok in enumerate(trafo_status):
        if not is_ok:
            model.deactivate_trafo(trafo_id)
