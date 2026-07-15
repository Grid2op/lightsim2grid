# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Low-level helpers shared by two or more of the ``_aux_add_*.py`` element
converters (bus resolution, current limits, per-unit view, transformer ratio /
phase-shift, remote voltage-control bus resolution)."""

import copy
import warnings

import numpy as np
import pandas as pd
import pypowsybl as pypo
from packaging import version

PP_BUG_RATIO_TAP_CHANGER = version.parse("1.9")
PYPOWSYBL_VER = version.parse(pypo.__version__)


def _aux_get_bus(vl_df, bus_df, first_bus_per_vl, el_type, df, conn_key="connected", bus_key="bus_id", vl_key="voltage_level_id"):
    if df.shape[0] == 0:
        # no element of this type so no problem
        return np.zeros(0, dtype=int), np.ones(0, dtype=bool), np.zeros(0, dtype=int)
    # retrieve which elements are disconnected
    mask_disco = ~df[conn_key]
    if mask_disco.any() and first_bus_per_vl is None:
        raise RuntimeError(f"Some elements of type {el_type} are disconnected (no bus attached to them) and you "
                           "are in 'use pypowsybl bus as grid2op substation' mode. The converter cannot assign "
                           f"substation for {el_type}:\n {mask_disco[mask_disco].index.to_numpy()}\n, this is not supported. "
                           "Please upgrade the pypowsybl grid (*eg* grid.xiidm of the grid2op env) "
                           "and manually reconnect any non connected elements there.")
    # retrieve the bus where the element are
    tmp_bus_id = df[bus_key].copy()

    if mask_disco.any():
        # element disconnected are, by default assigned to first bus of their substation
        el_disco = vl_df.loc[df.loc[mask_disco, vl_key]].index
        tmp_bus_id.loc[mask_disco] = first_bus_per_vl.loc[el_disco, "first_bus_name"].values

    # assign bus id
    bus_id = bus_df.loc[tmp_bus_id.values]["bus_global_id"].values.copy()

    if mask_disco.any():
        # deactivate the element not on the main component
        # wrong_component = bus_df.loc[tmp_bus_id.values]["connected_component"].values != 0
        # mask_disco[wrong_component] = True
        # assign bus -1 to disconnected elements
        bus_id[mask_disco] = -1

    sub_id = bus_df.loc[tmp_bus_id.values]["glop_sub_id"].values
    return bus_id, mask_disco.values, sub_id


def _aux_current_limits(element_ids, group_1, group_2, operational_limits):
    """Per-side thermal (current) limit, in **kA**, for a set of branch elements
    (lines or transformers), aligned to ``element_ids``.

    ``group_1``/``group_2`` are the elements' ``selected_limits_group_1``/``_2``
    (pandas Series, from ``net.get_lines()``/``net.get_2_windings_transformers()``
    called with ``all_attributes=True``) -- a branch can carry several unused limit
    groups (eg summer/winter), only the *selected* one per side is relevant.
    ``operational_limits`` is ``net.get_operational_limits()`` already filtered to
    ``type == "CURRENT"``.

    For each (element, side), the minimum finite value across all durations
    (permanent + temporary) within the side's selected limit group is kept --
    pypowsybl's ``~1.8e308`` "no limit for this duration" placeholder is ignored.
    NaN for a given (element, side) if it carries no current limit (or its selected
    group has none)."""
    n = len(element_ids)
    limit_a1_ka = np.full(n, np.nan)
    limit_a2_ka = np.full(n, np.nan)
    if operational_limits is None or operational_limits.shape[0] == 0 or n == 0:
        return limit_a1_ka, limit_a2_ka
    # ignore the "no limit for this duration" placeholder
    operational_limits = operational_limits[operational_limits["value"] < 1e30]
    pos = pd.Series(np.arange(n), index=element_ids)
    for out, side_str, group in ((limit_a1_ka, "ONE", group_1), (limit_a2_ka, "TWO", group_2)):
        rows = operational_limits[operational_limits.index.get_level_values("side") == side_str]
        if rows.shape[0] == 0:
            continue
        rows = rows.reset_index()[["element_id", "group_name", "value"]]
        sel_group = group.reindex(element_ids).rename("sel_group").rename_axis("element_id").reset_index()
        rows = rows.merge(sel_group, on="element_id", how="inner")
        rows = rows[rows["group_name"] == rows["sel_group"]]
        if rows.shape[0] == 0:
            continue
        min_per_el = rows.groupby("element_id")["value"].min() / 1000.  # A -> kA
        idx = pos.reindex(min_per_el.index).dropna().astype(int)
        out[idx.values] = min_per_el.loc[idx.index].values
    return limit_a1_ka, limit_a2_ka


def _aux_get_operational_limits_current(net):
    """``net.get_operational_limits()`` filtered to ``type == "CURRENT"``, shared
    by the line and transformer converters (see :func:`_aux_current_limits`).
    ``None`` on legacy pypowsybl versions that don't expose it."""
    try:
        ol_current = net.get_operational_limits()
        return ol_current[ol_current.index.get_level_values("type") == "CURRENT"]
    except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
        return None


def _aux_ensure_net_pu(net, net_pu):
    """Return `net_pu` unchanged if already provided by the caller; otherwise
    build the per-unit view of `net` once. Factored out so the zero-impedance
    bus-fusion pre-pass (which needs per-unit r/x before the "normal" per-unit
    view is built) and the main line/trafo processing can share a single,
    lazily-built `net_pu` instead of constructing it twice."""
    if net_pu is not None:
        return net_pu
    if hasattr(net, "per_unit"):
        net_pu = copy.deepcopy(net)
        net_pu.per_unit = True
    else:
        # legacy pypowsybl mode: this did not exist
        from pypowsybl.network import PerUnitView
        net_pu = PerUnitView(net)
        warnings.warn("The `PerUnitView` (python side) is less efficient and less "
                      "tested that the equivalent java class. Please upgrade pypowsybl version")
    return net_pu


def _aux_get_2wt_all_attrs(net_like):
    """`get_2_windings_transformers(all_attributes=True)`, falling back to the
    plain call on legacy pypowsybl versions that don't support it."""
    try:
        return net_like.get_2_windings_transformers(all_attributes=True)
    except TypeError:
        return net_like.get_2_windings_transformers()


def _aux_trafo_alpha(df_trafo_pu, net):
    """Phase-shift angle, in **degree**, for every transformer in `df_trafo_pu`
    (a per-unit `get_2_windings_transformers()` table)."""
    if 'alpha' in df_trafo_pu:
        return np.rad2deg(df_trafo_pu['alpha'].values)  # given in radian by pypowsybl
    if net.get_phase_tap_changers().shape[0] > 0:
        raise RuntimeError("Phase tap changer are not handled by the pypowsybl converter "
                           "when not accessible using the 'alpha' columns "
                           "of the net (once per unit). Please upgrade pypowsybl."
                           "NB: phase tap change are handled by lightsim2grid)")
    return np.zeros(df_trafo_pu.shape[0])


def _aux_trafo_rho(df_trafo_pu, ratio_tap_changer):
    """Turns ratio for every transformer in `df_trafo_pu` (a per-unit
    `get_2_windings_transformers()` table)."""
    if "rho" in df_trafo_pu:
        return 1. * df_trafo_pu["rho"].values
    # in powsybl (https://javadoc.io/doc/com.powsybl/powsybl-core/latest/com/powsybl/iidm/network/TwoWindingsTransformer.html)
    #  rho = transfo.getRatedU2() / transfo.getRatedU1()
    # * (transfo.getRatioTapChanger() != null ? transfo.getRatioTapChanger().getCurrentStep().getRho() : 1);
    # * (transfo.getPhaseTapChanger() != null ? transfo.getPhaseTapChanger().getCurrentStep().getRho() : 1);
    ratio = 1. * (df_trafo_pu["rated_u2"].values / df_trafo_pu["rated_u1"].values)
    has_r_tap_changer = np.isin(df_trafo_pu.index, ratio_tap_changer.index)
    if PYPOWSYBL_VER <= PP_BUG_RATIO_TAP_CHANGER:
        # bug in per unit view in both python and java
        ratio[has_r_tap_changer] = 1. * ratio_tap_changer.loc[df_trafo_pu.loc[has_r_tap_changer].index, "rho"].values
    return ratio


def _aux_regulated_bus_view_ids(net, regulated_ids):
    """Resolve voltage-controller regulated elements to their terminal bus.

    A remote voltage controller (generator or SVC) points at a `regulated_element_id`
    which, depending on the grid topology, may be a busbar section (node/breaker) or
    any bus-connected equipment (load, generator, ...). OpenLoadFlow regulates the
    voltage of that element's terminal bus. This returns, for each id in
    `regulated_ids`, the bus-view bus id of its terminal (which is also the index of
    the `bus_df` table used elsewhere in this converter), as a numpy object array.

    .. warning::
        This mapping is resolved **once**, when the grid is imported. lightsim2grid
        then stores the regulated bus by its (fixed) global id. If, afterwards, the
        regulated element is moved to another bus *in lightsim2grid* (e.g. through a
        ``change_bus_*`` / topology change), the controller keeps regulating the bus
        resolved here: the two grids will desynchronise. Re-import the grid (or update
        the regulated bus manually with ``set_gen_regulated_bus``) if you need to
        follow such a topology change.

    .. todo::
        Track the regulated *element* (not the resolved bus) so that moving it to
        another bus inside lightsim2grid updates the regulated bus automatically.
    """
    lookup = {}
    for getter in ("get_busbar_sections", "get_loads", "get_generators",
                   "get_static_var_compensators", "get_shunt_compensators",
                   "get_batteries", "get_vsc_converter_stations",
                   "get_lcc_converter_stations"):
        try:
            df = getattr(net, getter)()
        except Exception:
            continue
        if df.shape[0] == 0 or "bus_id" not in df.columns:
            continue
        for el_id, bus_view_id in zip(df.index, df["bus_id"].values):
            lookup.setdefault(el_id, bus_view_id)
    missing = [rid for rid in regulated_ids if rid not in lookup]
    if missing:
        raise RuntimeError("Some voltage controllers regulate an element whose bus could "
                           "not be resolved (unknown or disconnected regulated element): "
                           f"{missing}. This is not supported at the moment.")
    return np.array([lookup[rid] for rid in regulated_ids], dtype=object)
