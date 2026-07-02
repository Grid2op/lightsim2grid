# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings
import copy
from collections import deque
import numpy as np
import pandas as pd
import pypowsybl as pypo
from typing import Dict, Iterable, Optional, Union
from packaging import version
from ...lightsim2grid_cpp import LSGrid # type: ignore


from ._aux_handle_slack import handle_slack_iterable, handle_slack_one_el

PP_BUG_RATIO_TAP_CHANGER = version.parse("1.9")
PYPOWSYBL_VER = version.parse(pypo.__version__)

# suffix of the synthetic branch lightsim2grid line id for a dangling line's
# equivalent boundary branch (see `_aux_dangling_lines_fictitious`)
_DANGLING_BOUNDARY_LINE_SUFFIX = "@boundary_line"

# OpenLoadFlow's hardcoded fallback droop (in `AbstractLfGenerator.DEFAULT_DROOP`,
# "why not") used for every generator whose `activePowerControl` extension does not
# set its own droop, under the ``PROPORTIONAL_TO_GENERATION_P_MAX`` balance type --
# see `_default_distributed_slack`.
_OLF_DEFAULT_DROOP = 4.0


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


def _aux_dangling_lines_fictitious(net, sort_index):
    """One fictitious 1-bus "substation" per dangling line, at the far
    (boundary) end of its series impedance.

    Real full grids never have any dangling line -- they only appear when
    zooming into a sub-area with
    ``network.reduce_by_ids_and_depths(..., with_boundary_lines=True)``
    (*eg* in ``reduce_and_compare.py`` / ``validate_olf_vs_lightsim.py``),
    which cuts the grid and represents everything beyond the cut as a
    ``DanglingLine`` (series r/x/g/b + a constant-power p0/q0 "load" at the
    boundary). Without this, ``_from_pypowsybl`` silently drops that
    boundary injection (it never reads ``get_dangling_lines()``), which can
    be tens to hundreds of MW and makes the reduced-grid comparison
    meaningless.

    IIDM convention: the shunt admittance (g, b) of a dangling line sits
    entirely on the connectable (network) side; the boundary side carries
    none (same convention already used for transformers' single ``h``, see
    ``trafo_h`` below).

    Returns ``(bus_extra, vl_extra, df_dl)``; ``bus_extra``/``vl_extra`` are
    ``None`` and ``df_dl`` is empty if the grid has no dangling line.
    ``df_dl`` carries two extra columns, ``boundary_bus_id`` /
    ``boundary_vl_id``, used later to attach the equivalent branch + load.
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


def _aux_phase_shift_rx_tables(trafo_index, net):
    """Per-transformer phase-shift -> series-impedance dependency for
    :meth:`LSGrid.set_trafo_shift_dependent_rx`.

    Returns two lists-of-lists aligned to ``trafo_index``: the phase-shift sample
    points ``alpha`` (in **radian**, matching lightsim2grid's ``shift_``) and the
    matching r/x correction (**in percent**), read from the pypowsybl
    phase-tap-changer steps. Inner lists are empty for transformers without a
    phase-tap-changer.

    pypowsybl exposes the transformer r / x at the *neutral* tap and only folds the
    tap into ``rho`` / ``alpha``; the per-step r/x deltas (percent) are dropped. They
    matter for phase-shifting transformers whose series impedance varies with the
    shift (RTE PSTs: tens of MW of through-flow). On those grids ``r% == x%`` for
    every step, so a single correction table is applied to both r and x."""
    n = len(trafo_index)
    alpha = [[] for _ in range(n)]
    corr = [[] for _ in range(n)]
    try:
        ptc = net.get_phase_tap_changers()
        steps = net.get_phase_tap_changer_steps()
    except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
        return alpha, corr
    if ptc.shape[0] == 0 or steps.shape[0] == 0:
        return alpha, corr
    have = set(ptc.index)
    for i, tid in enumerate(trafo_index):
        if tid not in have:
            continue
        try:
            st = steps.loc[tid]
        except KeyError:
            continue
        a = np.deg2rad(np.atleast_1d(st["alpha"].to_numpy(dtype=float)))
        x = np.atleast_1d(st["x"].to_numpy(dtype=float))  # r% == x% on these PSTs
        alpha[i] = [float(v) for v in a]
        corr[i] = [float(v) for v in x]
    return alpha, corr


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


def _default_distributed_slack(net, df_gen):
    """Build the default distributed slack matching OpenLoadFlow's behaviour.

    Returns an *ordered* ``{generator_name: weight}`` dict (the reference generator
    first, so it becomes ``slack_ids[0]`` -- the angle datum) suitable for
    ``gen_slack_id``, or ``None`` if no participating generator is found (the caller
    then falls back to :meth:`LSGrid.assign_slack_to_most_connected`).

    The participating set and the sharing key mirror OLF's default distributed
    slack:

    * participants are the connected generators with ``max_p > 0`` that take part in
      active-power control (``participate`` flag of the ``activePowerControl``
      extension); if that extension is absent, every connected producing generator
      participates;
    * participants are restricted to the *main synchronous component* (OLF's
      ``ActivePowerDistribution.run`` only ever considers
      ``LfSynchronousNetwork.getBuses()``) -- **not** to any particular country:
      ``ActivePowerDistribution.filterParticipatingBuses`` has no country logic at
      all, and OLF's own ``slackBusCountryFilter`` is empty by default, so foreign
      border-equivalent generators (*eg* a cross-border interconnection's
      equivalent generator, which can carry a very large ``max_p``/weight) do
      participate. An earlier version of this function restricted participants to
      the country hosting the most buses; that was a guess, not something OLF
      actually does, and it was found to materially skew which generators absorb
      the slack mismatch on real multi-country RTE grids -- removed;
    * the per-generator weight matches OLF's default ``PROPORTIONAL_TO_GENERATION_P_MAX``
      balance type exactly (``GenerationActivePowerDistributionStep.getParticipationFactor``,
      ``MAX`` case): ``max_p / droop``, where ``droop`` is the ``activePowerControl``
      extension's own droop when it sets one (> 0), otherwise OLF's hardcoded
      ``DEFAULT_DROOP = 4`` for *every* generator regardless of the extension.
      **Not** the extension's ``participation_factor`` -- that key is only used
      under the (different, not reproduced here) ``PROPORTIONAL_TO_GENERATION_PARTICIPATION_FACTOR``
      balance type, even though real grids often set both fields together.

    The reference generator (angle datum, ``slack_ids[0]``) is a separate concern
    from the participant/weight logic above and IS restricted to the main country
    (the country hosting the most buses) when available: follow each candidate's
    step-up transformer to the highest-voltage bus it reaches, pick the bus
    carrying the most generator active power, then the generator injecting the
    most into it -- a lightsim2grid-only proxy for OLF's own
    ``slackBusSelectionMode=MOST_MESHED``. Without the country restriction here, a
    single cross-border aggregate generator (*eg* a "whole neighbouring country"
    equivalent injection, `max_p` in the tens of GW) reaches a high-voltage tie bus
    and dominates this heuristic, becoming the reference -- which behaves nothing
    like a real MOST_MESHED domestic bus and was empirically much worse (see the
    country-filter history above; this restriction previously covered the
    participant set too, until that was found to disagree with OLF and narrowed
    down to just the reference pick).
    """
    names = df_gen.index
    n = len(names)
    if n == 0:
        return None
    connected = df_gen["connected"].to_numpy(bool)
    max_p = df_gen["max_p"].to_numpy(float)
    target_p = df_gen["target_p"].to_numpy(float)
    gen_bus = df_gen["bus_id"].to_numpy()

    # participation set + sharing key from the activePowerControl extension
    try:
        apc = net.get_extensions("activePowerControl")
    except Exception:
        apc = None
    if apc is not None and len(apc) and "participate" in apc.columns:
        # a generator listed in the extension uses its flag; one absent from it
        # participates by default (OLF treats a missing extension as participating)
        participate = apc["participate"].reindex(names).fillna(True).to_numpy(bool)
        if "droop" in apc.columns:
            droop = apc["droop"].reindex(names).to_numpy(float)
        else:
            droop = np.full(n, np.nan)
    else:
        participate = np.ones(n, bool)  # no (or empty) extension -> everything participates
        droop = np.full(n, np.nan)

    # sharing key (PROPORTIONAL_TO_GENERATION_P_MAX): max_p / droop, OLF's own
    # DEFAULT_DROOP=4 when the extension does not set a droop of its own (pypowsybl
    # reports an unset droop as 0.0, not NaN, hence the `> 0.` guard rather than
    # `np.isfinite`).
    droop_used = np.where(np.isfinite(droop) & (droop > 0.), droop, _OLF_DEFAULT_DROOP)
    weight = max_p / droop_used

    # participants: connected, *started* (positive target P -- OLF does not
    # distribute on zero-MW generators), participating, with a usable positive weight.
    mask = connected & participate & (target_p > 0.) & np.isfinite(weight) & (weight > 0.)

    # main-component filter: a slack generator must sit in the main component,
    # otherwise lightsim2grid's `consider_only_main_component` deactivates its
    # (islanded) bus and the solver then trips on a disconnected slack. "Main" is
    # the largest synchronous component, which matches the line + transformer
    # connectivity lightsim2grid keeps (HVDC excluded) and mirrors OLF's own
    # `LfSynchronousNetwork.getBuses()` restriction. Real grids carry many small
    # boundary islands that satisfy the participation tests above.
    df_bus = net.get_buses(attributes=["synchronous_component", "voltage_level_id"])
    main_sync = df_bus["synchronous_component"].value_counts().idxmax()
    gen_sync = df_gen["bus_id"].map(df_bus["synchronous_component"]).to_numpy()
    mask = mask & (gen_sync == main_sync)

    if not mask.any():
        return None

    # reference (angle datum): follow the step-up transformer(s) to the
    # highest-voltage node, take the node with the most generator active power, then
    # the generator injecting the most into it.
    vls = net.get_voltage_levels()
    nomv = net.get_buses()["voltage_level_id"].map(vls["nominal_v"]).to_dict()
    tw = net.get_2_windings_transformers(attributes=["bus1_id", "bus2_id", "connected1", "connected2"])
    tc = tw[tw["connected1"] & tw["connected2"]]
    adj = {}
    for b1, b2 in zip(tc["bus1_id"], tc["bus2_id"]):
        adj.setdefault(b1, []).append(b2)
        adj.setdefault(b2, []).append(b1)

    def follow_gsu(bus, max_hops=4):
        seen = {bus}
        q = deque([(bus, 0)])
        best = bus
        while q:
            b, d = q.popleft()
            if nomv.get(b, 0.) > nomv.get(best, 0.):
                best = b
            if d < max_hops:
                for nb in adj.get(b, []):
                    if nb not in seen:
                        seen.add(nb)
                        q.append((nb, d + 1))
        return best

    # candidate pool for the *reference* (angle datum) only: restricted to the main
    # country when available. This is a lightsim2grid-only heuristic proxy for OLF's
    # own `slackBusSelectionMode=MOST_MESHED` reference pick (unrelated to country),
    # added because a cross-border tie's equivalent generator (*eg* a single "whole
    # neighbouring country" aggregate with a very large `max_p`) otherwise dominates
    # the "most generation at the highest voltage reached" heuristic below and gets
    # picked as reference, which is a poor MOST_MESHED proxy and was empirically much
    # worse than restricting to the domestic candidates. Unlike this reference pick,
    # the *participant* set (`mask` above) is intentionally left unrestricted, since
    # OLF's actual active-power distribution (`ActivePowerDistribution.
    # filterParticipatingBuses`) has no country logic at all.
    ref_mask = mask
    subs = net.get_substations()
    if "country" in subs.columns and subs["country"].notna().any():
        vl2sub = vls["substation_id"]
        bus_country = df_bus["voltage_level_id"].map(vl2sub).map(subs["country"])
        main_country = bus_country.value_counts().idxmax()
        gen_country = df_gen["voltage_level_id"].map(vl2sub).map(subs["country"]).to_numpy()
        if (mask & (gen_country == main_country)).any():
            ref_mask = mask & (gen_country == main_country)

    cand = np.where(mask)[0]
    ref_cand = np.where(ref_mask)[0]
    conn = {i: follow_gsu(gen_bus[i]) for i in ref_cand}
    cvolt = {i: nomv.get(conn[i], 0.) for i in ref_cand}
    vmax = max(cvolt.values())
    node_p = {}
    for i in ref_cand:
        if cvolt[i] == vmax:
            node_p[conn[i]] = node_p.get(conn[i], 0.) + target_p[i]
    ref_node = max(node_p, key=node_p.get)
    ref_i = max((i for i in ref_cand if cvolt[i] == vmax and conn[i] == ref_node),
                key=lambda i: target_p[i])

    order = [int(ref_i)] + [int(i) for i in cand if i != ref_i]
    return {names[i]: float(weight[i]) for i in order}


def init(net : pypo.network.Network,
         gen_slack_id: Union[int, str, Iterable[str], Dict[str, float]] = None,
         slack_bus_id: int = None,
         sn_mva : float = 100.,  # only used if not present in the grid
         sort_index : bool =True, 
         f_hz : float = 50.,  # unused
         net_pu : Optional[pypo.network.Network] = None,
         only_main_component : bool =True,
         return_sub_id: bool=False,
         n_busbar_per_sub: Optional[int]=None,  # new in 0.9.1
         buses_for_sub:Optional[bool]=None,  # new in 0.9.1
         init_vm_pu:float=1.06,
         keep_half_open_lines: bool=False,
         convert_dangling_lines: bool=False,
         ) -> LSGrid:
    """
    This function is available under the `init_from_pypowsybl` in lightsim2grid
    

    .. code-block:: python
    
        from lightsim2grid.network import init_from_pypowsybl
        
    .. warning::
        It is not available if the `pypowsybl` python package is not installed.
    
    :param net: The pypowsybl network
    :type net: pypo.network.Network
    
    :param gen_slack_id: The id of the generator that should be used as the slack
                         (either it's given by id (int) or by name (str))
    :type gen_slack_id: Union[int, str]
    
    :param slack_bus_id: If you don't provide a generator ID as a slack bus, you can
            provide a bus id (int). We do not recommend setting the slack this way.            
    :type slack_bus_id: int
    
    :param sn_mva: The nominal apparent power used when converting the grid to 
                   per unit. It is only used if the pypowsybl grid 
                   has no `_nominal_apparent_power` attribute. 
                   **Advanced usage**.            
    :type sn_mva: float
    
    :param sort_index: Whether you want to sort the indexes of all the 
                       pypowsybl tables (*eg* get_loads() or *get_buses()*) or not.
                       Sorting the grid tables is preferable if you want to be 
                       "future proof" and don't want to depend on pandas version
                       (same order is guaranteed). Not sorting the grid will give
                       easier comparison of results with pypowsybl.          
    :type sn_mva: bool
    
    :param f_hz: Not used currently (frequency of the grid)
    :type net_pu: float
    
    :param net_pu: If you have already converted the grid in "per unit" then
                   you can pass it as the `net_pu` argument. Otherwise this 
                   function will do it.
                   **Advanced usage**.
    :type net_pu: Optional[pypo.network.Network]
    
    :param only_main_component: If this is True, then only the main component (*ie*
                                the one containing the slack bus) will be used. All
                                equipments not part of this component will be
                                deactivated (switched-off). **NB** currently
                                lightsim2grid will diverge if the grid is not connected,
                                this option might then "hide" some equipements from
                                the grid (silently) but you have higher chances of
                                convergence.
    :type only_main_component: bool
    
    :param return_sub_id: **Advanced usage**. If you want to retrieve the id of the
                          equipments as "tables". Used only for `LightSimBackend`
    :type return_sub_id: bool
                          
    :param n_busbar_per_sub: Currently, lightsim2grid works well with a constant 
                             number of independant buses that can be made at each 
                             substations. It can be infered from the grid or 
                             set with this attribute. We recommend to leave it 
                             to `None` (which corresponds to the "infer it from 
                             the grid" behaviour) in most cases.
    :type n_busbar_per_sub: Opional[int]
    
    :param buses_for_sub: Whether the lightsim2grid substation will correspond to buses
                          of the pypowsybl grid (if buses_for_sub is `True`).
                          Alternatively, if buses_for_sub is `False`, the
                          lightsim2grid susbtation will correspond to
                          pypowsybl voltage level (read from net.get_voltage_levels()).
                          buses_for_sub==`True` is a "legacy" behaviour.
    :type buses_for_sub: bool
    
    :param init_vm_pu: The voltage magnitude with which the init vector of AC powerflow
                       will be set.
    :type init_vm_pu: float

    :param keep_half_open_lines: If True, a powerline or transformer connected on only
                       one terminal (``connected1 != connected2``, *eg* a dangling
                       boundary stub in a real grid) is modeled as "half-open": the
                       energized side is kept in the admittance matrix and the open end
                       is Kron-reduced out, instead of deactivating the whole branch.
                       This sets ``synch_status_both_side=False`` on the returned model,
                       so a later one-sided topology change is no longer mirrored to the
                       other side. Branches disconnected on *both* sides are still fully
                       deactivated. Default False (whole-branch deactivation, as before).
    :type keep_half_open_lines: bool

    :param convert_dangling_lines: If True, every IIDM ``DanglingLine`` (*eg*
                       produced by ``network.reduce_by_ids_and_depths(...,
                       with_boundary_lines=True)`` when zooming into a
                       sub-area) is converted to its equivalent branch +
                       constant-power load: a fictitious 1-bus "substation" at
                       the boundary end, a line carrying the dangling line's own
                       r/x/g/b (shunt entirely on the local side, matching
                       transformers' single ``h``), and a load consuming its
                       p0/q0. Without this (the default), dangling lines are
                       silently ignored -- fine for real full grids, which
                       never have any, but it drops a real boundary injection
                       (can be hundreds of MW) whenever they do appear. Off by
                       default to keep existing behaviour unchanged; the
                       reduce/validate debug scripts turn it on.
    :type convert_dangling_lines: bool

    :return: The properly initialized network.
    :rtype: :class:`LSGrid`
    """
    model = LSGrid()
    if hasattr(net, "_nominal_apparent_power"):
        sn_mva_used = getattr(net, "_nominal_apparent_power")
    else:
        sn_mva_used = float(sn_mva)
    model.set_sn_mva(sn_mva_used)
    model.set_init_vm_pu(float(init_vm_pu))
    if keep_half_open_lines:
        # allow branches connected on a single terminal (the open end is Kron-reduced
        # in the C++ model); a one-sided disconnection is no longer mirrored to the
        # other side.
        model.set_synch_status_both_side(False)

    if gen_slack_id is not None and slack_bus_id is not None:
        raise RuntimeError("Impossible to intialize a grid with both gen_slack_id and slack_bus_id")
    
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
            
    # all_buses_vn_kv = np.concatenate([all_buses_vn_kv for _ in range(n_busbar_per_sub)])
    model.init_bus(n_sub_ls,
                   n_busbar_per_sub_ls,
                   all_buses_vn_kv,
                   0, 0  # unused
                   )
    model._ls_to_orig = ls_to_orig
    model._max_nb_bus_per_sub = n_busbar_per_sub_ls
    model.set_substation_names(sub_names)
    model.set_bus_voltage_limits(all_buses_vmin_kv.astype(float), all_buses_vmax_kv.astype(float))
        
    # do the generators
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
    too_small = min_q_aux < min_float_value
    min_q_aux[too_small] = min_float_value
    min_q = min_q_aux.astype(np.float32)

    max_q_aux = 1. * max_q_src.values
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
        gen_reg_bus_view = _aux_regulated_bus_view_ids(net, bus_reg[mask_remote_gen])
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
    
    # for loads
    if sort_index:
        df_load = net.get_loads().sort_index()
    else:
        df_load = net.get_loads()
    if df_dl.shape[0] > 0:
        # equivalent constant-power load at each dangling line's boundary bus
        # (see `_aux_dangling_lines_fictitious`); p0/q0 are already in MW/MVAr,
        # same raw convention as `net.get_loads()`.
        df_load_extra = pd.DataFrame({
            "p0": df_dl["p0"].values,
            "q0": df_dl["q0"].values,
            "bus_id": df_dl["boundary_bus_id"].values,
            "connected": True,
            "voltage_level_id": df_dl["boundary_vl_id"].values,
        }, index=[f"{dl_id}@boundary_load" for dl_id in df_dl.index])
        df_load = pd.concat([df_load, df_load_extra])
    load_bus, load_disco, load_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "load", df_load)
    model.init_loads(df_load["p0"].values,
                     df_load["q0"].values,
                     load_bus
                     )
    for load_id, is_disco in enumerate(load_disco):
        if is_disco:
            model.deactivate_load(load_id)
    model.set_load_names(df_load.index)

    # thermal (current) limits for lines / trafos, see `_aux_current_limits`
    try:
        ol_current = net.get_operational_limits()
        ol_current = ol_current[ol_current.index.get_level_values("type") == "CURRENT"]
    except Exception:  # noqa: BLE001 - not available on legacy pypowsybl
        ol_current = None

    # for lines
    if sort_index:
        df_line = net.get_lines().sort_index()
    else:
        df_line = net.get_lines()
    try:
        line_limit_groups = net.get_lines(all_attributes=True)[["selected_limits_group_1", "selected_limits_group_2"]]
    except (TypeError, KeyError):
        # not available on legacy pypowsybl / grid with no limit group at all
        line_limit_groups = pd.DataFrame(
            {"selected_limits_group_1": np.nan, "selected_limits_group_2": np.nan}, index=df_line.index
        )

    # per unit
    if net_pu is None:
        if hasattr(net, "per_unit"):
            net_pu = copy.deepcopy(net)
            net_pu.per_unit = True
        else:
            # legacy pypowsybl mode: this did not exist
            from pypowsybl.network import PerUnitView
            net_pu = PerUnitView(net)
            warnings.warn("The `PerUnitView` (python side) is less efficient and less "
                          "tested that the equivalent java class. Please upgrade pypowsybl version")
    df_line_pu = net_pu.get_lines().loc[df_line.index]
    if df_dl.shape[0] > 0:
        # equivalent branch (local bus -> fictitious boundary bus) for each
        # dangling line, see `_aux_dangling_lines_fictitious`. Shunt admittance
        # (g, b) sits entirely on the local (network) side, none on the
        # boundary side -- same convention already used for the single `h` of
        # a transformer (see `trafo_h` below).
        df_dl_pu = net_pu.get_dangling_lines().loc[df_dl.index]
        line_ids_extra = [f"{dl_id}{_DANGLING_BOUNDARY_LINE_SUFFIX}" for dl_id in df_dl.index]
        df_line_extra = pd.DataFrame({
            "bus1_id": df_dl["bus_id"].values,
            "bus2_id": df_dl["boundary_bus_id"].values,
            "connected1": df_dl["connected"].values,
            "connected2": True,
            "voltage_level1_id": df_dl["voltage_level_id"].values,
            "voltage_level2_id": df_dl["boundary_vl_id"].values,
        }, index=line_ids_extra)
        df_line_extra_pu = pd.DataFrame({
            "r": df_dl_pu["r"].values,
            "x": df_dl_pu["x"].values,
            "g1": df_dl_pu["g"].values,
            "b1": df_dl_pu["b"].values,
            "g2": 0.,
            "b2": 0.,
        }, index=line_ids_extra)
        df_line = pd.concat([df_line, df_line_extra])
        df_line_pu = pd.concat([df_line_pu, df_line_extra_pu])
    line_r = df_line_pu["r"].values
    line_x = df_line_pu["x"].values
    line_h_or = (df_line_pu["g1"].values + 1j * df_line_pu["b1"].values)
    line_h_ex = (df_line_pu["g2"].values + 1j * df_line_pu["b2"].values)
    lor_bus, lor_disco, lor_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "line (side 1)", df_line, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    lex_bus, lex_disco, lex_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "line (side 2)", df_line, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_powerlines_full(line_r,
                               line_x,
                               line_h_or,
                               line_h_ex,
                               lor_bus,
                               lex_bus
                              )
    for line_id, (is_or_disc, is_ex_disc) in enumerate(zip(lor_disco, lex_disco)):
        if is_or_disc and is_ex_disc:
            model.deactivate_powerline(line_id)
        elif is_or_disc:
            model.deactivate_powerline_side1(line_id) if keep_half_open_lines else model.deactivate_powerline(line_id)
        elif is_ex_disc:
            model.deactivate_powerline_side2(line_id) if keep_half_open_lines else model.deactivate_powerline(line_id)
    model.set_line_names(df_line.index)
    line_limit_a1_ka, line_limit_a2_ka = _aux_current_limits(
        df_line.index,
        line_limit_groups["selected_limits_group_1"],
        line_limit_groups["selected_limits_group_2"],
        ol_current,
    )
    model.set_line_thermal_limit(line_limit_a1_ka, line_limit_a2_ka)

    # for trafo
    # I extract trafo with `all_attributes=True` so that I have access to the `rho`
    try:
        df_trafo_not_sorted = net.get_2_windings_transformers(all_attributes=True)
    except TypeError:
        # not available in legacy pypowsybl version
        df_trafo_not_sorted = net.get_2_windings_transformers()
        
    if sort_index:
        df_trafo = df_trafo_not_sorted.sort_index()
    else:
        df_trafo = df_trafo_not_sorted
    
    try :
        df_trafo_pu = net_pu.get_2_windings_transformers(all_attributes=True)
    except TypeError:
        df_trafo_pu = net_pu.get_2_windings_transformers()
    df_trafo_pu = df_trafo_pu.loc[df_trafo.index]
    ratio_tap_changer = net_pu.get_ratio_tap_changers()
    
    if 'alpha' in df_trafo_pu:
        shift_ = np.rad2deg(df_trafo_pu['alpha'].values)  # given in radian by pypowsybl
    else:
        if net.get_phase_tap_changers().shape[0] > 0:
            raise RuntimeError("Phase tap changer are not handled by the pypowsybl converter "
                               "when not accessible using the 'alpha' columns "
                               "of the net (once per unit). Please upgrade pypowsybl."
                               "NB: phase tap change are handled by lightsim2grid)")
        shift_ = np.zeros(df_trafo.shape[0])
    # tap is side 2 in IIDM
    is_tap_side1 = np.zeros(df_trafo.shape[0], dtype=bool)   
    # neutral-tap impedance (the phase-shift -> r/x dependence of phase-shifting
    # transformers is handled by lightsim2grid as a function of the shift alpha, see
    # the model.set_trafo_shift_dependent_rx(...) call below).
    trafo_r = df_trafo_pu["r"].values
    trafo_x = df_trafo_pu["x"].values
    trafo_h = (df_trafo_pu["g"].values + 1j * df_trafo_pu["b"].values)

    # now get the ratio
    # in lightsim2grid (cpp)
    if "rho" in df_trafo_pu:
        ratio = 1. * df_trafo_pu["rho"].values
    else:
        # in powsybl (https://javadoc.io/doc/com.powsybl/powsybl-core/latest/com/powsybl/iidm/network/TwoWindingsTransformer.html)
        #  rho = transfo.getRatedU2() / transfo.getRatedU1()
        # * (transfo.getRatioTapChanger() != null ? transfo.getRatioTapChanger().getCurrentStep().getRho() : 1);
        # * (transfo.getPhaseTapChanger() != null ? transfo.getPhaseTapChanger().getCurrentStep().getRho() : 1);

        ratio = 1. * (df_trafo_pu["rated_u2"].values / df_trafo_pu["rated_u1"].values)
        has_r_tap_changer = np.isin(df_trafo_pu.index, ratio_tap_changer.index)
    
        if PYPOWSYBL_VER <= PP_BUG_RATIO_TAP_CHANGER:
            # bug in per unit view in both python and java
            ratio[has_r_tap_changer] = 1. * ratio_tap_changer.loc[df_trafo_pu.loc[has_r_tap_changer].index, "rho"].values

    tor_bus, tor_disco, tor_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "trafo (side 1)", df_trafo, conn_key="connected1", bus_key="bus1_id", vl_key="voltage_level1_id")
    tex_bus, tex_disco, tex_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "trafo (side 2)", df_trafo, conn_key="connected2", bus_key="bus2_id", vl_key="voltage_level2_id")
    model.init_trafo(trafo_r,
                     trafo_x,
                     trafo_h,
                     ratio,
                     shift_,  # in degree !
                     is_tap_side1,
                     tor_bus,
                     tex_bus,
                     False,  # ignore_tap_side_for_phase_shift is False for pypowsybl
                     )
    for t_id, (is_or_disc, is_ex_disc) in enumerate(zip(tor_disco, tex_disco)):
        if is_or_disc and is_ex_disc:
            model.deactivate_trafo(t_id)
        elif is_or_disc:
            model.deactivate_trafo_side1(t_id) if keep_half_open_lines else model.deactivate_trafo(t_id)
        elif is_ex_disc:
            model.deactivate_trafo_side2(t_id) if keep_half_open_lines else model.deactivate_trafo(t_id)
    model.set_trafo_names(df_trafo.index)
    if "selected_limits_group_1" in df_trafo.columns:
        trafo_group_1 = df_trafo["selected_limits_group_1"]
        trafo_group_2 = df_trafo["selected_limits_group_2"]
    else:
        # not available on legacy pypowsybl
        trafo_group_1 = pd.Series(np.nan, index=df_trafo.index)
        trafo_group_2 = pd.Series(np.nan, index=df_trafo.index)
    trafo_limit_a1_ka, trafo_limit_a2_ka = _aux_current_limits(
        df_trafo.index, trafo_group_1, trafo_group_2, ol_current
    )
    model.set_trafo_thermal_limit(trafo_limit_a1_ka, trafo_limit_a2_ka)
    # phase-shifting transformers: declare the (alpha -> r/x correction) dependency so
    # lightsim2grid keeps the series impedance right when the shift changes, without any
    # "tap" concept (the per-step r/x deltas pypowsybl carries only in its tap steps).
    ps_alpha, ps_rx_corr = _aux_phase_shift_rx_tables(df_trafo.index, net)
    if any(len(a) for a in ps_alpha):
        model.set_trafo_shift_dependent_rx(True, ps_alpha, ps_rx_corr)

    # for shunt
    if sort_index:
        df_shunt = net.get_shunt_compensators().sort_index()
    else:
        df_shunt = net.get_shunt_compensators()
        
    sh_bus, sh_disco, sh_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "shunts", df_shunt)    
    shunt_kv = voltage_levels.loc[df_shunt["voltage_level_id"].values]["nominal_v"].values
    model.init_shunt(df_shunt["g"].values * shunt_kv**2,
                     -df_shunt["b"].values * shunt_kv**2,
                     sh_bus
                    )
    for shunt_id, disco in enumerate(sh_disco):
        if disco:
           model.deactivate_shunt(shunt_id)
    model.set_shunt_names(df_shunt.index)

    # for Static Var Compensators (SVC): VOLTAGE (local/remote, optional slope),
    # REACTIVE_POWER (fixed Q) or OFF, all solved through the bordered
    # VoltageControl NR extension. A grid with no SVC declares no controller and
    # stays byte-identical to before this feature.
    if sort_index:
        df_svc = net.get_static_var_compensators().sort_index()
    else:
        df_svc = net.get_static_var_compensators()
    nb_svc = df_svc.shape[0]
    svc_bus, svc_disco, svc_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "svc", df_svc)

    # SvcContainer.RegulationMode: OFF=0, VOLTAGE=1, REACTIVE_POWER=2
    OFF_MODE, VOLTAGE_MODE, REACTIVE_POWER_MODE = 0, 1, 2
    svc_mode = np.zeros(nb_svc, dtype=int)
    svc_reg_bus = svc_bus.copy()             # regulated bus (own bus unless remote)
    svc_reg_vn = np.ones(nb_svc)             # nominal v (kV) of the regulated bus
    svc_slope_pu = np.zeros(nb_svc)
    svc_target_q_inject = np.zeros(nb_svc)
    if nb_svc:
        mode_str = df_svc["regulation_mode"].values.astype(str)
        if "regulating" in df_svc.columns:
            regulating = df_svc["regulating"].values.astype(bool)
        else:
            # legacy pypowsybl: "OFF" is encoded directly in the regulation mode
            regulating = np.ones(nb_svc, dtype=bool)
        svc_mode[(mode_str == "VOLTAGE") & regulating] = VOLTAGE_MODE
        svc_mode[(mode_str == "REACTIVE_POWER") & regulating] = REACTIVE_POWER_MODE

        # resolve the regulated bus, mirroring the generator `regulated_element_id`
        # logic above (busbar section in node/breaker grids, or any bus-connected
        # element). Local SVCs keep their own bus.
        svc_vl = copy.deepcopy(df_svc["voltage_level_id"].values)
        if "regulated_element_id" in df_svc.columns:
            reg_id = copy.deepcopy(df_svc["regulated_element_id"].values)
            reg_id = np.where(reg_id == "", df_svc.index, reg_id)
            mask_svc_remote = reg_id != df_svc.index.values
            if mask_svc_remote.any():
                # TODO: resolved once at import; if the regulated element later changes
                # bus inside lightsim2grid this stays frozen and desynchronises from the
                # original grid (see `_aux_regulated_bus_view_ids` warning).
                svc_reg_bus_view = _aux_regulated_bus_view_ids(net, reg_id[mask_svc_remote])
                svc_reg_bus[mask_svc_remote] = bus_df.loc[svc_reg_bus_view, "bus_global_id"].values
                svc_vl[mask_svc_remote] = bus_df.loc[svc_reg_bus_view, "voltage_level_id"].values
        svc_reg_vn = voltage_levels.loc[svc_vl, "nominal_v"].values

        # IIDM gives the SVC reactive setpoint in the receptor (load) convention,
        # whereas lightsim2grid stamps Q with the generator-injection convention
        # (Phase 0 probe: SVC target_q=+30 absorbs 30 MVar) -> negate.
        target_q = df_svc["target_q"].values.astype(float)
        target_q = np.where(np.isfinite(target_q), target_q, 0.)
        svc_target_q_inject = -target_q

        # optional voltage/reactive-power slope ("droop"), in kV/MVar:
        #   s_pu = slope[kV/MVar] * sn_mva / vn_kv(regulated bus)   (Phase 0 probe #1)
        try:
            df_slope = net.get_extensions("voltagePerReactivePowerControl")
        except Exception:
            # extension tables may be unavailable on (very) old pypowsybl versions
            df_slope = None
        if df_slope is not None and df_slope.shape[0]:
            for svc_pos, svc_id in enumerate(df_svc.index):
                if svc_id in df_slope.index:
                    slope_kv_per_mvar = float(df_slope.loc[svc_id, "slope"])
                    svc_slope_pu[svc_pos] = slope_kv_per_mvar * sn_mva_used / svc_reg_vn[svc_pos]

    if nb_svc:
        target_v = df_svc["target_v"].values.astype(float)
        # target_v (kV) -> pu at the regulated bus; NaN (REACTIVE_POWER / OFF) -> 1 pu
        target_vm_pu = np.where(np.isfinite(target_v), target_v, svc_reg_vn) / svc_reg_vn
        # pypowsybl/IIDM gives b_min/b_max in SIEMENS (physical susceptance), while
        # lightsim2grid's SvcContainer expects them in per unit (base sn_mva, at the
        # regulated bus's nominal voltage -- same base as target_vm_pu/slope above).
        # Without this conversion the SVC's modeled Q range is smaller than its real one
        # by a factor of (nominal_v_kv)^2 / sn_mva (eg ~500x on a 225kV/100MVA bus),
        # making a perfectly healthy SVC collapse to a near-zero Q range: it saturates
        # (or, since check_solution never enforces SVC Q limits, silently hits a "hard"
        # voltage pin at its own target_vm_pu regardless of what Q that would truly take)
        # long before the real device would.
        b_min = df_svc["b_min"].values.astype(float) * (svc_reg_vn ** 2) / sn_mva_used
        b_max = df_svc["b_max"].values.astype(float) * (svc_reg_vn ** 2) / sn_mva_used
    else:
        target_vm_pu = np.zeros(0)
        b_min = np.zeros(0)
        b_max = np.zeros(0)

    model.init_svcs([int(m) for m in svc_mode],
                    target_vm_pu,
                    svc_target_q_inject,
                    svc_slope_pu,
                    b_min,
                    b_max,
                    svc_reg_bus.astype(np.int32),
                    svc_bus.astype(np.int32))
    for svc_id, disco in enumerate(svc_disco):
        if disco:
            model.deactivate_svc(svc_id)
    model.set_svc_names(df_svc.index)

    # for hvdc: vsc / lcc converter stations + hvdc lines, possibly carrying
    # the angle-droop ("AC emulation") extension
    if sort_index:
        df_dc = net.get_hvdc_lines().sort_index()
    else:
        df_dc = net.get_hvdc_lines()
    df_vsc = net.get_vsc_converter_stations()
    df_lcc = net.get_lcc_converter_stations()
    # the vsc / lcc frames have different columns (target_v / power_factor...):
    # the concatenation puts NaN where an attribute does not exist for a type
    df_stations = pd.concat([df_vsc, df_lcc])
    nb_dc = df_dc.shape[0]
    _max_hvdc_mva = 1.0e7  # when pypowsybl exposes NaN limits

    df_station1 = df_stations.loc[df_dc["converter_station1_id"].values]
    df_station2 = df_stations.loc[df_dc["converter_station2_id"].values]
    hvdc_bus_from_id, hvdc_from_disco, hvdc_sub_from_id = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "hvdc (side 1)", df_station1)
    hvdc_bus_to_id, hvdc_to_disco, hvdc_sub_to_id = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "hvdc (side 2)", df_station2)

    def _aux_hvdc_station_data(df_side):
        # type: 0 = VSC, 1 = LCC (ConverterStationContainer convention)
        is_lcc = df_side.index.isin(df_lcc.index)
        types = np.where(is_lcc, 1, 0).astype(int)
        loss_factor = df_side["loss_factor"].values / 100.  # pypowsybl % -> fraction
        loss_factor = np.where(np.isfinite(loss_factor), loss_factor, 0.)
        vreg_on = df_side["voltage_regulator_on"].values.astype(bool) if nb_dc else np.zeros(0, dtype=bool)
        vreg_on = vreg_on & ~is_lcc  # lcc never regulates (NaN -> random bool otherwise)
        vl_kv = voltage_levels.loc[df_side["voltage_level_id"].values]["nominal_v"].values
        vset_pu = df_side["target_v"].values / vl_kv
        vset_pu = np.where(np.isfinite(vset_pu), vset_pu, 1.0)
        qset = df_side["target_q"].values
        qset = np.where(np.isfinite(qset), qset, 0.)
        min_q = df_side["min_q"].values
        min_q = np.where(np.isfinite(min_q), min_q, -_max_hvdc_mva)
        max_q = df_side["max_q"].values
        max_q = np.where(np.isfinite(max_q), max_q, _max_hvdc_mva)
        power_factor = df_side["power_factor"].values
        power_factor = np.where(np.isfinite(power_factor), power_factor, 1.0)
        return types, loss_factor, vreg_on, vset_pu, qset, min_q, max_q, power_factor

    type1, lf1, vreg1, vset1, qset1, min_q1, max_q1, pf1 = _aux_hvdc_station_data(df_station1)
    type2, lf2, vreg2, vset2, qset2, min_q2, max_q2, pf2 = _aux_hvdc_station_data(df_station2)

    # 0 = side 1 rectifier, 1 = side 2 rectifier (HvdcLineContainer convention)
    converters_mode = np.where(df_dc["converters_mode"].values.astype(str) == "SIDE_1_RECTIFIER_SIDE_2_INVERTER", 0, 1).astype(int)
    p_setpoint_mw = df_dc["target_p"].values.astype(float).copy()
    if (~np.isfinite(p_setpoint_mw)).any():
        warnings.warn("Some non finite values are found for hvdc target_p, they have been replaced by 0.")
        p_setpoint_mw[~np.isfinite(p_setpoint_mw)] = 0.
    r_ohm = df_dc["r"].values.astype(float)
    nominal_v_kv = df_dc["nominal_v"].values.astype(float)
    max_p_mw = df_dc["max_p"].values.astype(float)
    max_p_mw = np.where(np.isfinite(max_p_mw), max_p_mw, _max_hvdc_mva)

    # the angle-droop active power control ("AC emulation"), an IIDM extension
    droop_enabled = np.zeros(nb_dc, dtype=bool)
    droop_p0_mw = np.zeros(nb_dc)
    droop_mw_per_deg = np.zeros(nb_dc)
    try:
        df_droop = net.get_extensions("hvdcAngleDroopActivePowerControl")
    except Exception:
        # extension tables may be unavailable on (very) old pypowsybl versions
        df_droop = None
    if df_droop is not None and df_droop.shape[0]:
        for hvdc_pos, line_id in enumerate(df_dc.index):
            if line_id not in df_droop.index:
                continue
            if not bool(df_droop.loc[line_id, "enabled"]):
                continue
            droop_enabled[hvdc_pos] = True
            droop_p0_mw[hvdc_pos] = float(df_droop.loc[line_id, "p0"])
            droop_mw_per_deg[hvdc_pos] = float(df_droop.loc[line_id, "droop"])

    model.init_hvdc_lines(hvdc_bus_from_id.astype(np.int32),
                          hvdc_bus_to_id.astype(np.int32),
                          [int(el) for el in type1],
                          [int(el) for el in type2],
                          lf1, lf2,
                          [bool(el) for el in vreg1],
                          [bool(el) for el in vreg2],
                          vset1, vset2,
                          qset1, qset2,
                          min_q1, max_q1, min_q2, max_q2,
                          pf1, pf2,
                          [int(el) for el in converters_mode],
                          p_setpoint_mw,
                          r_ohm,
                          nominal_v_kv,
                          [bool(el) for el in droop_enabled],
                          droop_p0_mw,
                          droop_mw_per_deg,
                          max_p_mw,  # pmax 1 -> 2: IIDM has a single max_p (open-loadflow convention)
                          max_p_mw,  # pmax 2 -> 1
                          )
    for hvdc_id, (is_or_disc, is_ex_disc, line_conn1, line_conn2) in enumerate(
            zip(hvdc_from_disco, hvdc_to_disco, df_dc["connected1"].values, df_dc["connected2"].values)):
        # a converter station with its own terminal open (eg its DC partner is
        # switched off, or its whole substation is dead) is NOT a dead branch: real
        # VSC stations (and OpenLoadFlow) keep the still-connected converter
        # injecting its scheduled P / regulating Q-V as a local device. Only fully
        # deactivate when BOTH stations are disconnected.
        or_disc = is_or_disc or (not line_conn1)
        ex_disc = is_ex_disc or (not line_conn2)
        if or_disc and ex_disc:
            model.deactivate_dcline(hvdc_id)
        elif or_disc:
            model.deactivate_dcline_side1(hvdc_id)
        elif ex_disc:
            model.deactivate_dcline_side2(hvdc_id)
    model.set_dcline_names(df_dc.index)
                
    # storage units (IIDM batteries). IIDM gives the battery setpoints in the
    # *generator* convention (positive target_p = power produced / injected) while
    # lightsim2grid stores storage as PQ in the *load* convention (positive = power
    # drawn from the grid, *ie* charging), same as pandapower and grid2op. We negate
    # to convert, and sanitize NaN (IIDM allows an unset target_q).
    if sort_index:
        df_batt = net.get_batteries().sort_index()
    else:
        df_batt = net.get_batteries()
    batt_bus, batt_disco, batt_sub = _aux_get_bus(voltage_levels, bus_df, first_bus_per_vl, "storage", df_batt)
    batt_p = df_batt["target_p"].values.astype(float)
    batt_q = df_batt["target_q"].values.astype(float)
    batt_p = np.where(np.isfinite(batt_p), batt_p, 0.)
    batt_q = np.where(np.isfinite(batt_q), batt_q, 0.)
    model.init_storages(-batt_p,  # IIDM generator convention -> lightsim2grid load convention
                        -batt_q,
                        batt_bus
                        )
    for batt_id, disco in enumerate(batt_disco):
        if disco:
           model.deactivate_storage(batt_id) 
    model.set_storage_names(df_batt.index)

    if gen_slack_id is None and slack_bus_id is None:
        # Default: reproduce OpenLoadFlow's distributed slack, sharing the
        # active-power mismatch over the participating generators (see
        # _default_distributed_slack). Returns None -- handled by the single
        # most-connected slack fallback below -- when no participating generator
        # is found.
        gen_slack_id = _default_distributed_slack(net, df_gen)

    if gen_slack_id is not None:
        if slack_bus_id is not None:
            raise RuntimeError("You provided both gen_slack_id and slack_bus_id "
                               "which is not possible.")
        if isinstance(gen_slack_id, (str, int, np.int32, np.int64, np.str_, tuple)):
            single_slack = True
            fun_slack = handle_slack_one_el
        else:
            single_slack = False
            fun_slack = handle_slack_iterable
        gen_slack_ids_int, gen_slack_weights = fun_slack(df_gen, gen_slack_id)
        if single_slack:
            if gen_slack_weights is None:
                raise RuntimeError(f"The slack {gen_slack_id} is disconnected.")
            gen_slack_ids_int = [gen_slack_ids_int]
            gen_slack_weights_fixed = [1.]
        else:
            gen_slack_weights_fixed = np.asarray([el if el is not None else np.nan for el in gen_slack_weights])
            mask_finite = np.isfinite(gen_slack_weights_fixed)
            if not mask_finite.any():
                raise RuntimeError(f"No connected generators match the slack {gen_slack_id}")
            gen_slack_weights_fixed[mask_finite] /= gen_slack_weights_fixed[mask_finite].sum()

        for gen_slack_id_int, gen_slack_weight in zip(gen_slack_ids_int, gen_slack_weights_fixed):
            if np.isfinite(gen_slack_weight):
                model.add_gen_slackbus(gen_slack_id_int, gen_slack_weight)
    elif slack_bus_id is not None:
        gen_bus = np.array([el.bus_id for el in model.get_generators()])
        gen_is_conn_slack = gen_bus == model._orig_to_ls[slack_bus_id]
        nb_conn = gen_is_conn_slack.sum()
        if nb_conn == 0:
            raise RuntimeError(f"There is no generator connected to bus {slack_bus_id}. It cannot be the slack")
        gen_slack_ids_int = []
        for gen_id, is_slack in enumerate(gen_is_conn_slack):
            if is_slack:
                gen_slack_ids_int.append(gen_id)
                model.add_gen_slackbus(gen_id, 1. / nb_conn)    
    else:
        # nothing provided and the default distributed slack found no participating
        # generator: fall back to a single slack on the most-connected generator bus
        bus_id, gen_id = model.assign_slack_to_most_connected()
        gen_slack_ids_int = [gen_id]
    
    # TODO checks
    # no 3windings trafo and other exotic stuff
        
    # and now deactivate all elements and nodes not in the main component
    if only_main_component:
        model.consider_only_main_component()
    else:
        # automatically disconnect non connected buses
        # (this is automatically done by consider_only_main_component)
        model.init_bus_status()  
        
    gen_sub = pd.DataFrame(index=df_gen.index, data={"sub_id": gen_sub})
    gen_sub["desired_slack"] = False
    gen_sub.loc[gen_sub.index[gen_slack_ids_int], "desired_slack"] = True
    load_sub = pd.DataFrame(index=df_load.index, data={"sub_id": load_sub})
    lor_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lor_sub})
    lex_sub = pd.DataFrame(index=df_line.index, data={"sub_id": lex_sub})
    tor_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tor_sub})
    tex_sub = pd.DataFrame(index=df_trafo.index, data={"sub_id": tex_sub})
    batt_sub = pd.DataFrame(index=df_batt.index, data={"sub_id": batt_sub})
    sh_sub = pd.DataFrame(index=df_shunt.index, data={"sub_id": sh_sub})
    hvdc_sub_from_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_from_id})
    hvdc_sub_to_id = pd.DataFrame(index=df_dc.index, data={"sub_id": hvdc_sub_to_id})
    
    # set the substation ID to which each object belong
    model.set_gen_to_subid(gen_sub["sub_id"].values)
    model.set_load_to_subid(load_sub["sub_id"].values)
    model.set_storage_to_subid(batt_sub["sub_id"].values)
    model.set_shunt_to_subid(sh_sub["sub_id"].values)
    model.set_line_to_sub1_id(lor_sub["sub_id"].values)
    model.set_line_to_sub2_id(lex_sub["sub_id"].values)
    model.set_trafo_to_sub1_id(tor_sub["sub_id"].values)
    model.set_trafo_to_sub2_id(tex_sub["sub_id"].values)
    if not return_sub_id:
        return model
    else:
        return model, (gen_sub, load_sub, (lor_sub, tor_sub), (lex_sub, tex_sub), batt_sub, sh_sub, hvdc_sub_from_id, hvdc_sub_to_id)
