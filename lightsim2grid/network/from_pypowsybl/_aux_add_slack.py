# Copyright (c) 2023-2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

from collections import deque

import numpy as np

from ._aux_handle_slack import handle_slack_iterable, handle_slack_one_el

# OpenLoadFlow's hardcoded fallback droop (in `AbstractLfGenerator.DEFAULT_DROOP`,
# "why not") used for every generator whose `activePowerControl` extension does not
# set its own droop, under the ``PROPORTIONAL_TO_GENERATION_P_MAX`` balance type --
# see `_default_distributed_slack`.
_OLF_DEFAULT_DROOP = 4.0


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
      the slack mismatch on real multi-country grids -- removed;
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


def _aux_add_slack(model, net, df_gen, gen_slack_id, slack_bus_id):
    """Resolve and assign the slack bus(es) of ``model``: an explicit
    ``gen_slack_id`` / ``slack_bus_id``, else OpenLoadFlow's default distributed
    slack (see :func:`_default_distributed_slack`), else a single slack on the
    most-connected generator bus. Returns ``gen_slack_ids_int``, the (0-based,
    ``df_gen``-indexed) list of slack generator ids -- used by `initLSGrid.py`
    to mark ``gen_sub["desired_slack"]``."""
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

    return gen_slack_ids_int
