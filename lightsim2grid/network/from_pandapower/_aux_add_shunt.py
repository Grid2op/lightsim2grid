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

    With ``step_dependency_table``, the power is instead the row of ``step`` in
    ``shunt_characteristic_table`` (same rated voltage), and ``p_mw`` / ``q_mvar`` are not read.

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
    p_step = pp_net.shunt["p_mw"].to_numpy(dtype=float) * v_ratio
    q_step = pp_net.shunt["q_mvar"].to_numpy(dtype=float) * v_ratio
    p_now = p_step * step
    q_now = q_step * step
    tables = _aux_step_tables(pp_net, v_ratio, step)
    for sh_id, table in tables.items():
        p_now[sh_id], q_now[sh_id] = table[step[sh_id]]

    model.init_shunt(p_now,
                     q_now,
                     pp_bus_to_ls(pp_net.shunt["bus"].to_numpy(), pp_to_ls)
                     )
    for sh_id, is_connected in enumerate(pp_net.shunt["in_service"].to_numpy()):
        if not is_connected:
            # shunt is deactivated
            model.deactivate_shunt(sh_id)

    with_sections = _aux_shunt_sections(model, pp_net, step, p_step, q_step, tables)
    _aux_shunt_regulation(model, pp_net, pp_to_ls, with_sections)


def _aux_step_tables(pp_net, v_ratio, step):
    """The power per step of the shunts with a ``step_dependency_table``, as pandapower reads
    it: ``{position of the shunt: {step: (p, q)}}``, at 1 pu, a NaN in the table being 0.

    Refuses what ``pp.runpp`` cannot read either, a shunt whose step has no row included
    (for a grid with a single shunt, ``pp.runpp`` then silently drops it)."""
    if "step_dependency_table" not in pp_net.shunt:
        return {}
    has_table = pp_net.shunt["step_dependency_table"].fillna(False).to_numpy(dtype=bool)
    if not has_table.any():
        return {}
    if "shunt_characteristic_table" not in pp_net:
        raise RuntimeError("Some shunts have a step_dependency_table but the network has no "
                           "shunt_characteristic_table.")
    all_rows = pp_net.shunt_characteristic_table
    labels = pp_net.shunt.index
    ids = pp_net.shunt["id_characteristic_table"]
    tables = {}
    for sh_id in np.flatnonzero(has_table):
        id_characteristic = ids.iloc[sh_id]
        if ids.isna().iloc[sh_id]:
            raise RuntimeError(f"Shunt {labels[sh_id]} has a step_dependency_table and no id_characteristic_table.")
        rows = all_rows.loc[all_rows["id_characteristic"] == id_characteristic]
        if rows["step"].duplicated().any():
            raise RuntimeError(f"The characteristic {id_characteristic} of shunt {labels[sh_id]} has several rows "
                               "for one step.")
        p_mw = np.nan_to_num(rows["p_mw"].to_numpy(dtype=float)) * v_ratio[sh_id]
        q_mvar = np.nan_to_num(rows["q_mvar"].to_numpy(dtype=float)) * v_ratio[sh_id]
        table = {int(k): (p, q) for k, p, q in zip(rows["step"].to_numpy(), p_mw, q_mvar)}
        if step[sh_id] not in table:
            raise RuntimeError(f"The characteristic {id_characteristic} of shunt {labels[sh_id]} has no row for its "
                               f"step {step[sh_id]}.")
        tables[int(sh_id)] = table
    return tables


def _table_sections(table, max_step):
    """The cumulative p and q of sections ``1`` .. ``n`` read from a step table, ``n`` the
    last step ``max_step`` allows and the table has. ``None`` if a section table cannot say
    it: a table with power at step 0 (no section on means no power), or a step without a row."""
    if any(power != 0. for power in table.get(0, (0., 0.))):
        return None
    last = min(max_step, max(table))
    if last < 1 or any(k not in table for k in range(1, last + 1)):
        return None
    return [table[k][0] for k in range(1, last + 1)], [table[k][1] for k in range(1, last + 1)]


def _aux_shunt_sections(model, pp_net, step, p_step, q_step, tables):
    """One section per step, for every shunt with a finite integer ``max_step`` >= 1 and an
    integer ``step`` in ``[0, max_step]``. ``k`` sections draw ``k`` times the power of one
    step, or the row ``k`` of its step table (the convention of ``set_shunt_sections``:
    cumulative, at 1 pu). Returns the mask of the shunts given sections."""
    nb_shunt = pp_net.shunt.shape[0]
    with_sections = np.zeros(nb_shunt, dtype=bool)
    max_step = np.full(nb_shunt, np.nan)
    if "max_step" in pp_net.shunt:
        max_step = pp_net.shunt["max_step"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        ok = (np.isfinite(max_step) & (max_step >= 1.) & (max_step == np.round(max_step)) &
              (step == np.round(step)) & (step >= 0.) & (step <= max_step))
    for sh_id in np.flatnonzero(ok):
        sh_id = int(sh_id)
        if sh_id in tables:
            sections = _table_sections(tables[sh_id], int(max_step[sh_id]))
            if sections is None or step[sh_id] > len(sections[0]):
                continue
            p_sections, q_sections = sections
        else:
            n = np.arange(1, int(max_step[sh_id]) + 1, dtype=float)
            p_sections, q_sections = list(n * p_step[sh_id]), list(n * q_step[sh_id])
        model.set_shunt_sections(sh_id, int(step[sh_id]), p_sections, q_sections)
        with_sections[sh_id] = True
    no_sections = pp_net.shunt.index[[sh_id for sh_id in tables if not with_sections[sh_id]]].tolist()
    if no_sections:
        warnings.warn(f"The shunts {no_sections} have no sections: their table cannot be read as cumulative "
                      "sections 1 .. max_step (power at step 0, a missing step, or a step above max_step). "
                      "Their power is right at their current step, but it cannot be moved.")
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
