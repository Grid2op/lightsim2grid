# Copyright (c) 2020, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import numpy as np
from ._pp_bus_to_ls_bus import pp_bus_to_ls


def _aux_add_shunt(model, pp_net, pp_to_ls):
    """
    Add the shunts of the pp_net into the lightsim2grid "model"

    The power of a shunt is pandapower's: ``p_mw`` and ``q_mvar`` are per step, at ``vn_kv``,
    so a shunt with ``step`` steps on draws ``step * q_mvar * (vn_bus / vn_kv) ** 2`` at 1 pu.
    Its ``1`` .. ``max_step`` steps are its sections (``LSGrid.set_shunt_sections``), and a
    ``DiscreteShuntController`` says what they regulate (``set_shunt_section_regulation``):
    the data the ShuntVoltageControl outer loop reads.

    Parameters
    ----------
    model
    pp_net

    Returns
    -------

    """
    if "parallel" in pp_net.shunt and np.any(pp_net.shunt["parallel"].to_numpy() != 1):
        raise RuntimeError("Cannot handle 'parallel' sgen columns. Please duplicate the rows if that is the case. "
                           "Some pp_net.sgen[\"parallel\"] != 1 it is not handled by lightsim yet.")

    if pp_net.shunt.shape[0] == 0:
        # nothing to do if no shunt
        return

    nb_shunt = pp_net.shunt.shape[0]
    # pandapower scales p_mw / q_mvar (given at vn_kv) to the nominal voltage of the bus
    v_ratio = np.ones(nb_shunt)
    if "vn_kv" in pp_net.shunt:
        vn_kv = pp_net.shunt["vn_kv"].to_numpy(dtype=float)
        vn_bus = pp_net.bus.loc[pp_net.shunt["bus"], "vn_kv"].to_numpy(dtype=float)
        ok = np.isfinite(vn_kv) & (vn_kv > 0.)
        v_ratio[ok] = (vn_bus[ok] / vn_kv[ok]) ** 2
    step = np.ones(nb_shunt)
    if "step" in pp_net.shunt:
        step = pp_net.shunt["step"].to_numpy(dtype=float).copy()
        step[~np.isfinite(step)] = 1.
    has_table = np.zeros(nb_shunt, dtype=bool)
    if "step_dependency_table" in pp_net.shunt:
        has_table = pp_net.shunt["step_dependency_table"].fillna(False).to_numpy(dtype=bool)
        if has_table.any():
            warnings.warn("Some shunts have a step_dependency_table, which lightsim2grid does not read: their "
                          "power is taken as linear in the number of steps (q_mvar * step).")
    p_step = pp_net.shunt["p_mw"].to_numpy(dtype=float) * v_ratio
    q_step = pp_net.shunt["q_mvar"].to_numpy(dtype=float) * v_ratio

    model.init_shunt(p_step * step,
                     q_step * step,
                     pp_bus_to_ls(pp_net.shunt["bus"].to_numpy(), pp_to_ls)
                     )
    for sh_id, is_connected in enumerate(pp_net.shunt["in_service"].to_numpy()):
        if not is_connected:
            # shunt is deactivated
            model.deactivate_shunt(sh_id)

    with_sections = _aux_shunt_sections(model, pp_net, step, p_step, q_step, has_table)
    _aux_shunt_regulation(model, pp_net, pp_to_ls, with_sections)


def _aux_shunt_sections(model, pp_net, step, p_step, q_step, has_table):
    """One section per step, ``k`` sections drawing ``k`` times the power of one step (the
    convention of ``set_shunt_sections``: cumulative, at 1 pu), for every shunt with a finite
    integer ``max_step`` >= 1 and an integer ``step`` in ``[0, max_step]``. Returns the mask
    of the shunts given sections."""
    nb_shunt = pp_net.shunt.shape[0]
    with_sections = np.zeros(nb_shunt, dtype=bool)
    if "max_step" not in pp_net.shunt:
        return with_sections
    max_step = pp_net.shunt["max_step"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        ok = (np.isfinite(max_step) & (max_step >= 1.) & (max_step == np.round(max_step)) &
              (step == np.round(step)) & (step >= 0.) & (step <= max_step) & ~has_table)
    for sh_id in np.flatnonzero(ok):
        n = np.arange(1, int(max_step[sh_id]) + 1, dtype=float)
        model.set_shunt_sections(int(sh_id), int(step[sh_id]), list(n * p_step[sh_id]), list(n * q_step[sh_id]))
        with_sections[sh_id] = True
    return with_sections


def _aux_shunt_regulation(model, pp_net, pp_to_ls, with_sections):
    """What the sections regulate, read from pandapower's ``DiscreteShuntController``: the
    voltage of its ``controlled_bus`` at ``vm_set_pu``, with a deadband of twice its ``tol``."""
    if "controller" not in pp_net or pp_net.controller.shape[0] == 0:
        return
    try:
        from pandapower.control import DiscreteShuntController
    except ImportError:  # older pandapower: no shunt controller
        return
    skipped = []
    for obj, in_service in zip(pp_net.controller["object"], pp_net.controller["in_service"]):
        if not in_service or not isinstance(obj, DiscreteShuntController):
            continue
        label = obj.shunt_index
        if label not in pp_net.shunt.index:
            continue
        sh_id = int(pp_net.shunt.index.get_loc(label))
        if not with_sections[sh_id]:
            skipped.append(label)
            continue
        target = float(obj.vm_set_pu)
        deadband = 2. * float(obj.tol)
        ls_bus = int(pp_bus_to_ls(np.array([obj.controlled_bus]), pp_to_ls)[0])
        model.set_shunt_section_regulation(sh_id, bool(np.isfinite(target)),
                                           target if np.isfinite(target) else 0.,
                                           deadband if np.isfinite(deadband) else 0.,
                                           ls_bus)
    if skipped:
        warnings.warn(f"The controllers of the shunts {skipped} are ignored: these shunts have no sections "
                      "lightsim2grid can model (see `_aux_shunt_sections`).")
