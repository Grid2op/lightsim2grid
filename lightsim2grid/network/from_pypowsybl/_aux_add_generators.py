# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import copy

import numpy as np

from ._aux_common import _aux_get_bus, _aux_regulated_bus_view_ids


def _aux_add_generators(model, net, sort_index, voltage_levels, bus_df, first_bus_per_vl):
    """Add every generator of ``net`` to ``model``. Returns ``(df_gen, gen_sub)``:
    ``df_gen`` is reused by the slack-assignment phase (`_aux_add_slack.py`),
    ``gen_sub`` by the final substation-id bookkeeping in `initLSGrid.py`."""
    gen_attrs = [
        "connected", "max_p", "target_p", "target_v", "target_q", "p",
        "voltage_regulator_on", "regulated_element_id", "voltage_level_id", "bus_id",
        "min_q", "max_q", "min_q_at_target_p", "max_q_at_target_p",
    ]
    if sort_index:
        df_gen = net.get_generators(attributes=gen_attrs).sort_index()
    else:
        df_gen = net.get_generators(attributes=gen_attrs)

    # to handle encoding in 32 bits and overflow when "splitting" the Q values among
    min_float_value = np.finfo(np.float32).min * 1e-4 + 1.
    max_float_value = np.finfo(np.float32).max * 1e-4 + 1.
    # "min_q"/"max_q" (the flat MIN_MAX box) are NaN for a generator whose
    # reactive_limits_kind is CURVE -- pypowsybl only populates those through
    # "min_q_at_target_p"/"max_q_at_target_p" (the capability curve evaluated at the
    # generator's own target P, available even before any loadflow has been run,
    # unlike "min_q_at_p" which depends on a solved "p"). Without this, every
    # CURVE-kind generator silently got the "no limit" float32 sentinel below
    # regardless of its real reactive range.
    min_q_src = df_gen["min_q_at_target_p"].where(df_gen["min_q_at_target_p"].notna(), df_gen["min_q"])
    max_q_src = df_gen["max_q_at_target_p"].where(df_gen["max_q_at_target_p"].notna(), df_gen["max_q"])
    min_q_aux = 1. * min_q_src.values
    max_q_aux = 1. * max_q_src.values
    # malformed source curve data (eg a reactive capability curve point entered with
    # min_q/max_q swapped) can make the "at target p" interpolation yield min_q > max_q.
    # OpenLoadFlow tolerates this silently; lightsim2grid's GeneratorContainer::init
    # hard-rejects it (real case found on a real grid snapshot). Restore a valid
    # interval by sorting the pair instead of crashing -- this only ever affects the
    # (already tiny) width of the interval, never which generators get a reactive
    # constraint at all.
    swapped = min_q_aux > max_q_aux
    if swapped.any():
        min_q_aux[swapped], max_q_aux[swapped] = max_q_aux[swapped], min_q_aux[swapped].copy()
    too_small = min_q_aux < min_float_value
    min_q_aux[too_small] = min_float_value
    min_q = min_q_aux.astype(np.float32)

    too_big = np.abs(max_q_aux) > max_float_value
    max_q_aux[too_big] = np.sign(max_q_aux[too_big]) * max_float_value
    max_q = max_q_aux.astype(np.float32)
    min_q[~np.isfinite(min_q)] = min_float_value
    max_q[~np.isfinite(max_q)] = max_float_value
    gen_bus, gen_disco, gen_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "gen", df_gen)

    # remote voltage control: a generator may regulate the voltage of a *different*
    # bus, identified by `regulated_element_id` (its own id when controlling locally).
    # Resolve that element to the bus-view id of its terminal, then to the regulated
    # voltage level (used for the target_v -> pu conversion) and, further down, to the
    # lightsim2grid global bus id threaded into the C++ container.
    bus_reg = copy.deepcopy(df_gen["regulated_element_id"].values)
    # for oldest pypowsybl version, we could have "" there
    bus_reg = np.where(bus_reg == "", df_gen.index, bus_reg)
    vl_reg = copy.deepcopy(df_gen["voltage_level_id"].values)
    # a disconnected generator is deactivated below regardless of its (possibly
    # itself disconnected / unresolvable) regulated element, so exclude it here:
    # otherwise a disconnected generator remotely "regulating" a disconnected
    # busbar section (empty bus-view id) crashes the (moot, since deactivated)
    # bus resolution below.
    mask_remote_gen = (bus_reg != df_gen.index.values) & ~gen_disco
    gen_reg_bus_view = None
    if mask_remote_gen.any():
        remote_idx = np.nonzero(mask_remote_gen)[0]
        gen_reg_bus_view = _aux_regulated_bus_view_ids(net, bus_reg[mask_remote_gen])
        # a *connected* generator can still remotely regulate a busbar section that is
        # itself disconnected (found on a real grid snapshot: a de-energized
        # voltage level). pypowsybl's bus-view id for a disconnected element is '',
        # not NaN, and can't be resolved to any bus_df row. OLF converges fine on such
        # grids, so it must fall back to local voltage control in this situation;
        # mirror that instead of crashing.
        unresolved = gen_reg_bus_view == ""
        if unresolved.any():
            mask_remote_gen[remote_idx[unresolved]] = False
            gen_reg_bus_view = gen_reg_bus_view[~unresolved]
        vl_reg[mask_remote_gen] = bus_df.loc[gen_reg_bus_view, "voltage_level_id"].values
    model.init_generators_full(df_gen["target_p"].values,
                            #    df_gen["target_v"].values / voltage_levels.loc[df_gen["voltage_level_id"].values]["nominal_v"].values,
                               df_gen["target_v"].values / voltage_levels.loc[vl_reg]["nominal_v"].values,
                               df_gen["target_q"].values,
                               df_gen["voltage_regulator_on"].values,
                               min_q,
                               max_q,
                               gen_bus
                               )
    for gen_id, is_disco in enumerate(gen_disco):
        if is_disco:
            model.deactivate_gen(gen_id)
    model.set_gen_names(df_gen.index)

    # thread the regulated bus to the C++ generator container. Local generators keep
    # their own bus (already the C++ default), so a grid without any remote control
    # stays byte-identical to before this feature.
    # TODO: the regulated bus is resolved once, here, from the pypowsybl grid. If the
    # regulated element later changes bus inside lightsim2grid, this stays frozen and
    # the two grids desynchronise (see `_aux_regulated_bus_view_ids` warning).
    if mask_remote_gen.any():
        gen_reg_bus_global = bus_df.loc[gen_reg_bus_view, "bus_global_id"].values
        for gen_id, reg_bus in zip(np.nonzero(mask_remote_gen)[0], gen_reg_bus_global):
            model.set_gen_regulated_bus(int(gen_id), int(reg_bus))

    return df_gen, gen_sub
