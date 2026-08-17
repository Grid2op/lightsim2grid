# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Bus / substation numbering: the first phase of `init_from_pypowsybl`, every
other element converter (`_aux_add_generators.py`, `_aux_add_lines.py`, ...)
depends on its output (`voltage_levels`, `bus_df`, `first_bus_per_vl`). Also
covers dangling-line fictitious boundary buses and zero-impedance branch fusion,
both of which are bus-numbering concerns resolved before any other element type
is read (see `_aux_add_buses`'s docstring)."""

import warnings

import numpy as np
import pandas as pd

from typing import Dict

from ...lightsim2grid_cpp import LSGrid  # type: ignore
from ._aux_common import (
    _aux_ensure_net_pu,
    _aux_get_2wt_all_attrs,
    _aux_trafo_alpha,
    _aux_trafo_rho,
)
from ._my_const import _DANGLING_BOUNDARY_LINE_SUFFIX


def _aux_dangling_lines_fictitious(net, sort_index):
    """One fictitious 1-bus "substation" per dangling line, at the far
    (boundary) end of its series impedance.

    Real full grids never have any dangling line -- they only appear when
    zooming into a sub-area with
    ``network.reduce_by_ids_and_depths(..., with_boundary_lines=True)``,
    which cuts the grid and represents everything beyond the cut as a
    ``DanglingLine`` (series r/x/g/b + a constant-power p0/q0 "load" at the
    boundary). Without this, the converter silently drops that
    boundary injection (it never reads ``get_dangling_lines()``), which can
    be tens to hundreds of MW and makes the reduced-grid comparison
    meaningless.

    IIDM convention: the shunt admittance (g, b) of a dangling line sits
    entirely on the connectable (network) side; the boundary side carries
    none (same convention already used for transformers' single ``h``, see
    ``trafo_h`` in `_aux_add_lines.py`).

    Returns ``(bus_extra, vl_extra, df_dl)``; ``bus_extra``/``vl_extra`` are
    ``None`` and ``df_dl`` is empty if the grid has no dangling line.
    ``df_dl`` carries two extra columns, ``boundary_bus_id`` /
    ``boundary_vl_id``, used later to attach the equivalent branch + load.

    Deliberately built from ``net`` (physical units), not ``net_pu``: every
    field actually read off ``df_dl`` downstream is unit-agnostic (ids,
    ``connected``) or physical on purpose (``p0``/``q0`` in MW/MVAr, matching
    every other load in this module). The per-unit-sensitive series
    impedance (r/x/g/b) is read separately, from ``net_pu.get_dangling_lines()``
    (see the ``df_dl_pu`` fetch in `_aux_add_lines.py`) -- ``net``/``df_dl`` is
    never used for those.
    """
    df_dl = net.get_dangling_lines()
    if sort_index:
        df_dl = df_dl.sort_index()
    if df_dl.shape[0] == 0:
        return None, None, df_dl
    df_dl = df_dl.copy()
    df_dl["boundary_bus_id"] = [f"{dl_id}@boundary_bus" for dl_id in df_dl.index]
    df_dl["boundary_vl_id"] = [f"{dl_id}@boundary_vl" for dl_id in df_dl.index]
    nominal_v = net.get_voltage_levels().loc[df_dl["voltage_level_id"].values, "nominal_v"].values
    bus_extra = pd.DataFrame({"voltage_level_id": df_dl["boundary_vl_id"].values, "name": ""},
                             index=pd.Index(df_dl["boundary_bus_id"].values, name=df_dl.index.name or "id"))
    vl_extra = pd.DataFrame({"nominal_v": nominal_v},
                            index=pd.Index(df_dl["boundary_vl_id"].values, name="id"))
    return bus_extra, vl_extra, df_dl


def dangling_line_boundary_bus(model: LSGrid) -> Dict[str, int]:
    """Dangling-line id -> lightsim2grid global bus id of its fictitious
    boundary bus, for a grid built with ``init(..., convert_dangling_lines=True)``
    (see :func:`_aux_dangling_lines_fictitious`). Empty if the grid has no
    dangling line, or was built with ``convert_dangling_lines=False``.

    Unlike every other bus, a dangling line's boundary bus has no counterpart in
    ``network.get_buses()`` so it is not reachable through ``model._ls_to_orig``.
    Rather than stashing extra state on ``model`` (a pybind11 object without
    ``dynamic_attr``, so arbitrary attributes cannot be set on it), this is
    reconstructed from the synthetic boundary branch's own (real, C++-backed)
    ``bus2_id`` -- once built, that branch is an ordinary lightsim2grid line.
    """
    suffix = _DANGLING_BOUNDARY_LINE_SUFFIX
    return {el.name[:-len(suffix)]: el.bus2_id
            for el in model.get_lines() if el.name.endswith(suffix)}


def _aux_add_buses(model, net, net_pu, sort_index, buses_for_sub, n_busbar_per_sub,
                   convert_dangling_lines, fuse_zero_impedance_branches, zero_impedance_threshold_pu):
    """Assign every pypowsybl bus a lightsim2grid (substation, local-bus) pair,
    call ``model.init_bus``/``set_substation_names``/``set_bus_voltage_limits``,
    and (optionally) fuse (near-)zero-impedance branches' terminal buses
    together. Every other element converter depends on this phase's
    ``voltage_levels``/``bus_df``/``first_bus_per_vl`` outputs to resolve its
    own elements' buses (see ``_aux_get_bus`` in `_aux_common.py`).

    Returns ``(voltage_levels, bus_df, first_bus_per_vl, df_dl, net_pu,
    fused_line_ids, fused_trafo_ids)``.
    """
    # assign unique id to the buses
    bus_df_orig = net.get_buses()
    if sort_index:
        bus_df = bus_df_orig.sort_index().copy()
    else:
        bus_df = bus_df_orig.copy()
    bus_df["orig_id"] = np.arange(bus_df.shape[0])

    if sort_index:
        voltage_levels = net.get_voltage_levels().sort_index()
    else:
        voltage_levels = net.get_voltage_levels()

    # dangling lines (only present when the grid was cut with
    # `reduce_by_ids_and_depths(..., with_boundary_lines=True)`, see
    # `_aux_dangling_lines_fictitious`) each get their own fictitious 1-bus
    # "substation" at their boundary end, appended *after* the real buses so
    # `orig_id` / row-position based mappings (eg `_ls_to_orig`, used by
    # `lightsim_bus_to_iidm`) keep pointing at real pypowsybl buses only.
    if convert_dangling_lines:
        bus_extra, vl_extra, df_dl = _aux_dangling_lines_fictitious(net, sort_index)
    else:
        raw_dl = net.get_dangling_lines()
        bus_extra, vl_extra, df_dl = None, None, raw_dl.iloc[0:0]
        if raw_dl.shape[0] > 0:
            warnings.warn(f"{raw_dl.shape[0]} dangling line(s) found in the grid (eg from "
                          "`network.reduce_by_ids_and_depths(..., with_boundary_lines=True)`) "
                          "but `convert_dangling_lines=False` (the default): their boundary "
                          "injection is ignored. Pass `convert_dangling_lines=True` to model "
                          "them as an equivalent branch + constant-power load instead.")
    if df_dl.shape[0] > 0:
        if buses_for_sub is not None and buses_for_sub:
            raise RuntimeError("Dangling lines (eg from `network.reduce_by_ids_and_depths("
                               "..., with_boundary_lines=True)`) are not supported "
                               "together with buses_for_sub=True (legacy mode). Use "
                               "buses_for_sub=False (or None).")
        bus_extra = bus_extra.copy()
        bus_extra["orig_id"] = np.arange(bus_df.shape[0], bus_df.shape[0] + bus_extra.shape[0])
        bus_df = pd.concat([bus_df, bus_extra])
        voltage_levels = pd.concat([voltage_levels, vl_extra])

    all_buses_vn_kv = voltage_levels.loc[bus_df["voltage_level_id"].values]["nominal_v"].values
    nb_bus_per_vl = bus_df[["voltage_level_id", "name"]].groupby("voltage_level_id").count()

    if buses_for_sub is not None and buses_for_sub:
        # I am in a compatibility mode,
        # the "substation" in lightsim2grid will be read
        # from the buses in the original grid (and not from the
        # voltage levels)
        if n_busbar_per_sub is None:
            # setting automatically n_busbar_per_sub
            # to 1
            # TODO logger here
            n_busbar_per_sub = 1

        all_buses_vn_kv = voltage_levels.loc[bus_df["voltage_level_id"], "nominal_v"].values
        all_buses_vmin_kv = voltage_levels.loc[bus_df["voltage_level_id"], "low_voltage_limit"].values
        all_buses_vmax_kv = voltage_levels.loc[bus_df["voltage_level_id"], "high_voltage_limit"].values
        if n_busbar_per_sub > 1:
            all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
            all_buses_vmin_kv = np.concatenate([all_buses_vmin_kv for _ in range(n_busbar_per_sub)])
            all_buses_vmax_kv = np.concatenate([all_buses_vmax_kv for _ in range(n_busbar_per_sub)])
        n_sub_ls = bus_df.shape[0]
        ls_to_orig = np.zeros(all_buses_vn_kv.shape[0], dtype=int) - 1
        ls_to_orig[:n_sub_ls] = np.arange(n_sub_ls)
        n_busbar_per_sub_ls = n_busbar_per_sub
        bus_df["bus_global_id"] = np.arange(n_sub_ls)
        bus_df["glop_sub_id"] = np.arange(n_sub_ls)  # np.concatenate([np.arange(n_sub_ls) for _ in range(n_busbar_per_sub)])
        sub_names = bus_df.index.values.astype(str)
        voltage_levels["vl_id"] = bus_df[["voltage_level_id", "bus_global_id"]].groupby("voltage_level_id").min()
        first_bus_per_vl = None
    else:
        # the "substation" in lightsim2grid
        voltage_levels["nb_bus_per_vl"] = nb_bus_per_vl["name"]
        bus_df["name"] = [[el] for el in bus_df.index]
        voltage_levels["bus_names"] = bus_df[["name", "voltage_level_id"]].groupby("voltage_level_id").sum()

        bus_df["local_id"] = [voltage_levels.loc[el, "bus_names"].index(id_) + 1
                            for id_, el in zip(bus_df.index,
                                                bus_df["voltage_level_id"].values)]
        n_vl = voltage_levels.shape[0]
        voltage_levels["vl_id"] = np.arange(n_vl)
        bus_df["bus_global_id"] = [(loc_id - 1) * n_vl + voltage_levels.loc[vl, "vl_id"]
                                   for loc_id, vl in zip(
                                       bus_df["local_id"],
                                       bus_df["voltage_level_id"]
                                   )]
        nb_bus_per_vl_in_grid = nb_bus_per_vl.values.max()
        if n_busbar_per_sub is None:
            # setting automatically n_busbar_per_sub
            # to the value read from the grid
            # TODO logger here
            n_busbar_per_sub = int(nb_bus_per_vl_in_grid)
        elif n_busbar_per_sub < nb_bus_per_vl_in_grid:
            raise RuntimeError(f"The input pypowsybl grid counts some voltage levels "
                               f"with {nb_bus_per_vl_in_grid} independant buses, "
                               f"which is not compatible with the n_busbar_per_sub={n_busbar_per_sub} "
                               "given as input.")
        all_buses_vn_kv = voltage_levels["nominal_v"].values
        all_buses_vmin_kv = voltage_levels["low_voltage_limit"].values
        all_buses_vmax_kv = voltage_levels["high_voltage_limit"].values
        if n_busbar_per_sub > 1:
            all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
            all_buses_vmin_kv = np.concatenate([all_buses_vmin_kv for _ in range(n_busbar_per_sub)])
            all_buses_vmax_kv = np.concatenate([all_buses_vmax_kv for _ in range(n_busbar_per_sub)])
        n_sub_ls = voltage_levels.shape[0]
        voltage_levels["glop_sub_id"] = np.arange(voltage_levels.shape[0])
        n_busbar_per_sub_ls = n_busbar_per_sub
        ls_to_orig = np.zeros(all_buses_vn_kv.shape[0], dtype=int) - 1
        ls_to_orig[bus_df["bus_global_id"].values] = np.arange(bus_df.shape[0])
        sub_names = voltage_levels.index.values.astype(str)
        bus_df["glop_sub_id"] = voltage_levels.loc[bus_df["voltage_level_id"].values, "glop_sub_id"].values

        # retrieve the "first bus of every substation"
        # this is used to connected reconnected element
        first_bus_per_vl = bus_df[["glop_sub_id", "voltage_level_id", "name"]].groupby("glop_sub_id").first().reset_index().set_index("voltage_level_id")
        first_bus_per_vl["first_bus_name"] = [el[0] for el in first_bus_per_vl["name"].values]
        first_bus_per_vl.drop(["glop_sub_id", "name"], axis=1, inplace=True)

        # make sure every voltage level as a "first bus"
        # I would raise an error in the case a voltage is fully disconnected but...
        check_grid_ok = voltage_levels.index.isin(first_bus_per_vl.index)
        if np.any(~check_grid_ok):
            warnings.warn("There are some voltage levels without any connected buses. "
                          f"Check voltage levels {voltage_levels[~check_grid_ok].index}")
            # name_added_bus = voltage_levels[~check_grid_ok].index
            nm_vl_without_bus = voltage_levels.loc[~check_grid_ok, ["vl_id", "glop_sub_id"]]
            to_add = voltage_levels.loc[~check_grid_ok, ["vl_id", "glop_sub_id"]].copy()
            to_add.index = to_add.index.astype(str) + "added_bus"
            to_add["name"] = [[el] for el in to_add.index.astype(str)]
            to_add["v_mag"] = np.nan
            to_add["v_angle"] = np.nan
            to_add["connected_component"] = 99999
            to_add["synchronous_component"] = 99999
            to_add["voltage_level_id"] = nm_vl_without_bus.index
            to_add["orig_id"] = 9999
            to_add["local_id"] = 1
            to_add["bus_global_id"] = to_add["glop_sub_id"]
            for added_bus_nm in to_add.index:
                bus_df.loc[added_bus_nm] = to_add.loc[added_bus_nm]
            for vl_bus_added in nm_vl_without_bus.index:
                first_bus_per_vl.loc[vl_bus_added, "first_bus_name"] = vl_bus_added + "added_bus"
    # import pdb
    # pdb.set_trace()
    # all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
    model.init_bus(n_sub_ls,
                   n_busbar_per_sub_ls,
                   all_buses_vn_kv,
                   0, 0  # unused
                   )
    model._ls_to_orig = ls_to_orig
    model._max_nb_bus_per_sub = n_busbar_per_sub_ls
    # relevant kwargs downstream consumers (eg a pypowsybl-shaped result view) need to
    # recover conversion-time settings that are otherwise plain Python arguments lost
    # after this function returns -- see LSGrid._init_kwargs.
    model._init_kwargs = {"sort_index": str(sort_index), "buses_for_sub": str(buses_for_sub)}
    model.set_substation_names(sub_names)
    model.set_bus_voltage_limits(all_buses_vmin_kv.astype(float), all_buses_vmax_kv.astype(float))

    # fuse the two terminal buses of (near-)zero-impedance lines / neutral-tap
    # transformers, mirroring OpenLoadFlow's `lowImpedanceBranchMode`. Must run
    # before any `_aux_get_bus` call (generators, just below) so every element
    # type transparently picks up the fused `bus_global_id` for free -- see
    # `_aux_get_bus`'s single choke point at `bus_df["bus_global_id"]`.
    fused_line_ids = set()
    fused_trafo_ids = set()
    if fuse_zero_impedance_branches:
        net_pu = _aux_ensure_net_pu(net, net_pu)
        parent = np.arange(all_buses_vn_kv.shape[0])

        def _uf_find(x):
            while parent[x] != x:
                x = parent[x]
            return x

        def _uf_union(a, b):
            ra, rb = _uf_find(a), _uf_find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)  # deterministic: smaller id wins

        # -- zero-impedance LINES --
        df_line_fuse = net.get_lines()
        if sort_index:
            df_line_fuse = df_line_fuse.sort_index()
        df_line_fuse_pu = net_pu.get_lines().loc[df_line_fuse.index]
        both_conn = (df_line_fuse["connected1"].values & df_line_fuse["connected2"].values)
        is_zero = ((np.abs(df_line_fuse_pu["r"].values) < zero_impedance_threshold_pu)
                   & (np.abs(df_line_fuse_pu["x"].values) < zero_impedance_threshold_pu))
        line_candidate = both_conn & is_zero
        if line_candidate.any():
            cand = df_line_fuse[line_candidate]
            vn1 = voltage_levels.loc[cand["voltage_level1_id"].values, "nominal_v"].values
            vn2 = voltage_levels.loc[cand["voltage_level2_id"].values, "nominal_v"].values
            bad = ~np.isclose(vn1, vn2)
            if bad.any():
                raise RuntimeError(
                    f"Zero-impedance line(s) {list(cand.index[bad])} connect two buses at "
                    "different nominal voltages: this is inconsistent grid data, refusing "
                    "to fuse them (fuse_zero_impedance_branches=True)."
                )
            gid1 = bus_df.loc[cand["bus1_id"].values, "bus_global_id"].values
            gid2 = bus_df.loc[cand["bus2_id"].values, "bus_global_id"].values
            for a, b in zip(gid1, gid2):
                _uf_union(int(a), int(b))
            fused_line_ids.update(cand.index)

        # -- zero-impedance, neutral-tap, SAME-NOMINAL-VOLTAGE TRANSFORMERS --
        # NB: pypowsybl's per-unit `rho` is the deviation from the transformer's OWN
        # rated_u1/rated_u2 ratio (tap-changer effect), not the absolute turns ratio:
        # a genuine step-down transformer (eg rated_u1=225kV/rated_u2=90kV) reports
        # rho=1 at neutral tap even though it performs a real voltage transformation
        # (implicit in the differing per-unit bus voltage bases on each side, not in
        # `rho` itself). So `rho~=1` alone only means "no ADDITIONAL tap deviation" --
        # it does NOT mean "no transformation at all". A transformer is only a true,
        # fusable ideal wire when it is ALSO between two buses of the same nominal
        # voltage (same requirement as for lines, just not an error here since most
        # real transformers legitimately span different nominal voltages).
        df_trafo_fuse = _aux_get_2wt_all_attrs(net)
        if sort_index:
            df_trafo_fuse = df_trafo_fuse.sort_index()
        if df_trafo_fuse.shape[0] > 0:
            df_trafo_fuse_pu = _aux_get_2wt_all_attrs(net_pu).loc[df_trafo_fuse.index]
            ratio_tap_changer_fuse = net_pu.get_ratio_tap_changers()
            alpha_fuse = _aux_trafo_alpha(df_trafo_fuse_pu, net)
            rho_fuse = _aux_trafo_rho(df_trafo_fuse_pu, ratio_tap_changer_fuse)
            both_conn_t = (df_trafo_fuse["connected1"].values & df_trafo_fuse["connected2"].values)
            is_zero_t = ((np.abs(df_trafo_fuse_pu["r"].values) < zero_impedance_threshold_pu)
                        & (np.abs(df_trafo_fuse_pu["x"].values) < zero_impedance_threshold_pu))
            is_ideal_t = (np.abs(rho_fuse - 1.) < 1e-9) & (np.abs(alpha_fuse) < 1e-9)
            vn1_t = voltage_levels.loc[df_trafo_fuse["voltage_level1_id"].values, "nominal_v"].values
            vn2_t = voltage_levels.loc[df_trafo_fuse["voltage_level2_id"].values, "nominal_v"].values
            same_vn_t = np.isclose(vn1_t, vn2_t)
            trafo_candidate = both_conn_t & is_zero_t & is_ideal_t & same_vn_t
            if trafo_candidate.any():
                cand_t = df_trafo_fuse[trafo_candidate]
                gid1 = bus_df.loc[cand_t["bus1_id"].values, "bus_global_id"].values
                gid2 = bus_df.loc[cand_t["bus2_id"].values, "bus_global_id"].values
                for a, b in zip(gid1, gid2):
                    _uf_union(int(a), int(b))
                fused_trafo_ids.update(cand_t.index)

        rep = np.array([_uf_find(i) for i in range(parent.shape[0])])
        bus_df["bus_global_id"] = rep[bus_df["bus_global_id"].values]
        # a fused-away ("loser") bus keeps its own row in the lightsim2grid model
        # (fixed-size bus containers) but loses every element to its representative,
        # so it ends up disconnected with no solved voltage of its own. Stash the
        # representative lookup so a downstream result view (eg
        # `LightsimResultNetwork`) can report that bus's voltage as its
        # representative's -- the two are electrically the same node -- instead of
        # the disconnected bus's stale/zero value. See `LSGrid._bus_fusion_rep`.
        model._bus_fusion_rep = rep.astype(int)
        if fused_line_ids or fused_trafo_ids:
            # ids of the fused-away branches themselves, so a downstream result view
            # (eg `LightsimResultNetwork`) can tell "deactivated because fused" apart
            # from "deactivated because genuinely out of service" and, where the
            # physics allow it, reconstruct the flow through it (see
            # `LightsimResultNetwork._reconstruct_fused_branches`). "\x1f" (ASCII unit
            # separator) is used as delimiter: pypowsybl element ids never contain it.
            kwargs = dict(model._init_kwargs)
            kwargs["fused_line_ids"] = "\x1f".join(sorted(str(i) for i in fused_line_ids))
            kwargs["fused_trafo_ids"] = "\x1f".join(sorted(str(i) for i in fused_trafo_ids))
            model._init_kwargs = kwargs

    return voltage_levels, bus_df, first_bus_per_vl, df_dl, net_pu, fused_line_ids, fused_trafo_ids
