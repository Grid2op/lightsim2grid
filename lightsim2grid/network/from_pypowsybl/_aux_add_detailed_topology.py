# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""The detailed topology (switches inside each substation) of
`init_from_pypowsybl`: read pypowsybl's node-breaker view and hand it to the
C++ model.

Two passes. `_aux_scan_detailed_topology` runs BEFORE the buses are numbered:
it reads every voltage level's nodes, busbar sections and switches (a
bus-breaker voltage level gets a synthetic node model, see below), which node
every element terminal stands on, and two things `_aux_add_buses` needs --
the capacity of each substation (the most buses any switch configuration can
make there) and the order in which the C++ projection will number the buses of
each voltage level, so that the bus-view buses get those local ids up front.
`_aux_add_detailed_topology` runs AFTER every element phase: it declares the
topology on the model, sets the node of every terminal, projects the switches
and checks that nothing moved -- the model's own reading of the switches must
put every element exactly where pypowsybl's bus view had put it.

Node model of a BUS_BREAKER voltage level (IIDM's own mapping of its
bus-breaker view onto a node-breaker one): one node per bus-breaker bus, each
carrying a synthetic busbar section named after the bus; one node per element
terminal, joined to its bus' node by a synthetic breaker (open iff the terminal
is disconnected); the voltage level's own switches join bus nodes.
"""

import warnings

import numpy as np
import pandas as pd

from ...lightsim2grid_cpp import SwitchKind  # type: ignore

_KIND_OF = {
    "BREAKER": int(SwitchKind.BREAKER),
    "DISCONNECTOR": int(SwitchKind.DISCONNECTOR),
    "LOAD_BREAK_SWITCH": int(SwitchKind.LOAD_BREAK_SWITCH),
}
_INTERNAL_CONNECTION = int(SwitchKind.INTERNAL_CONNECTION)
_BREAKER = int(SwitchKind.BREAKER)

# (key, pypowsybl getter, side suffix, is a branch terminal)
_TERMINAL_TABLES = [
    ("gen", "get_generators", "", False),
    ("load", "get_loads", "", False),
    ("shunt", "get_shunt_compensators", "", False),
    ("svc", "get_static_var_compensators", "", False),
    ("storage", "get_batteries", "", False),
    ("line1", "get_lines", "1", True),
    ("line2", "get_lines", "2", True),
    ("trafo1", "get_2_windings_transformers", "1", True),
    ("trafo2", "get_2_windings_transformers", "2", True),
    ("vsc", "get_vsc_converter_stations", "", True),
    ("lcc", "get_lcc_converter_stations", "", True),
]


def _aux_detailed_topology_wanted(net, detailed_topology, buses_for_sub,
                                  fuse_zero_impedance_branches, convert_dangling_lines):
    """Resolve `init(..., detailed_topology=...)` to a bool, refusing the
    combinations the first version does not support."""
    if detailed_topology is None or detailed_topology is False:
        return False
    if isinstance(detailed_topology, str):
        if detailed_topology != "auto":
            raise ValueError(f"init_from_pypowsybl: detailed_topology must be True, False or \"auto\", "
                             f"not {detailed_topology!r}")
        if buses_for_sub or fuse_zero_impedance_branches or convert_dangling_lines:
            return False
        if net.get_switches().shape[0] > 0:
            return True
        vls = net.get_voltage_levels(all_attributes=True)
        return bool((vls["topology_kind"] == "NODE_BREAKER").any())
    if detailed_topology is not True:
        raise ValueError(f"init_from_pypowsybl: detailed_topology must be True, False or \"auto\", "
                         f"not {detailed_topology!r}")
    if buses_for_sub:
        raise RuntimeError("init_from_pypowsybl: detailed_topology=True is not compatible with "
                           "buses_for_sub=True (legacy mode): a voltage level's switches span several "
                           "bus-view buses, so the substation must be the voltage level.")
    if fuse_zero_impedance_branches:
        raise RuntimeError("init_from_pypowsybl: detailed_topology=True is not compatible with "
                           "fuse_zero_impedance_branches=True (not supported yet).")
    if convert_dangling_lines:
        raise RuntimeError("init_from_pypowsybl: detailed_topology=True is not compatible with "
                           "convert_dangling_lines=True (not supported yet).")
    return True


class _DetailedTopologyScan:
    """What `_aux_scan_detailed_topology` read off the network. Every id here
    is what the C++ model will speak: substations in voltage-level order
    (`vl_pos`), nodes local to their voltage level, switches and busbar
    sections in their grid-wide order."""

    def __init__(self):
        self.vl_ids = None               # voltage level ids, in substation order
        self.nb_nodes = None             # one per voltage level
        self.bbs = None                  # DataFrame: vl_pos, node, name, bus (bus-view bus or "")
        self.switches = None             # DataFrame: vl_pos, node1, node2, kind, open, retained, name
        self.terminals = {}              # key -> DataFrame indexed by element id: vl_pos, node, bus
        self.bus_key = None              # DataFrame indexed by bus-view bus id: key_group, key_idx
        self.bound_per_vl = None         # Series indexed by voltage level id


def _aux_scan_detailed_topology(net, sort_index):
    """Read the node-breaker view of every voltage level (see the module docstring)."""
    scan = _DetailedTopologyScan()

    vls = net.get_voltage_levels(all_attributes=True)
    if sort_index:
        vls = vls.sort_index()
    scan.vl_ids = list(vls.index)
    n_vl = len(scan.vl_ids)
    vl_pos = pd.Series(np.arange(n_vl), index=vls.index)
    is_nb = vls["topology_kind"].astype(str) == "NODE_BREAKER"
    nb_vls = set(vls.index[is_nb.values])
    nb_nodes = np.zeros(n_vl, dtype=int)

    # ---- node-breaker voltage levels: node count and internal connections ----
    # (the one per-voltage-level call; everything else is a bulk table)
    ic_rows = []
    for vl in vls.index[is_nb.values]:
        topo = net.get_node_breaker_topology(vl)
        if topo.nodes.shape[0] > 0:
            nb_nodes[vl_pos[vl]] = int(topo.nodes.index.max()) + 1
        for n1, n2 in topo.internal_connections[["node1", "node2"]].itertuples(index=False):
            ic_rows.append((int(vl_pos[vl]), int(n1), int(n2)))

    # ---- bus-breaker voltage levels: one node per bus-breaker bus ----
    bbv = net.get_bus_breaker_view_buses()
    if sort_index:
        bbv = bbv.sort_index()
    bbv = bbv[~bbv["voltage_level_id"].isin(nb_vls)]
    bbv = bbv.assign(node=bbv.groupby("voltage_level_id").cumcount())
    bb_bus_node = bbv["node"]  # indexed by bus-breaker bus id
    if bbv.shape[0] > 0:
        counts = bbv.groupby("voltage_level_id").size()
        nb_nodes[vl_pos[counts.index].values] = counts.values
    # nodes allocated so far in each bus-breaker voltage level (terminals come next)
    next_node = pd.Series(nb_nodes.copy(), index=vls.index)

    # ---- the terminals: which node each element end stands on ----
    synthetic_sw = []  # (vl_pos, node1, node2, kind, open, retained, name, seq=2)
    branch_count = pd.Series(0, index=vls.index)
    feeder_count = pd.Series(0, index=vls.index)
    for key, getter, side, is_branch in _TERMINAL_TABLES:
        try:
            df = getattr(net, getter)(all_attributes=True)
        except TypeError:  # legacy pypowsybl
            df = getattr(net, getter)()
        if sort_index:
            df = df.sort_index()
        vl_col, node_col, bbb_col = f"voltage_level{side}_id", f"node{side}", f"bus_breaker_bus{side}_id"
        bus_col, conn_col = f"bus{side}_id", f"connected{side}"
        term = pd.DataFrame(index=df.index)
        term["vl_pos"] = vl_pos.reindex(df[vl_col].values).values
        term["bus"] = df[bus_col].astype(str).values if bus_col in df.columns else ""
        term["node"] = -1
        in_nb = df[vl_col].isin(nb_vls).values
        if node_col in df.columns:
            node = df[node_col].to_numpy(int)
            term.loc[in_nb & (node >= 0), "node"] = node[in_nb & (node >= 0)]
        # bus-breaker voltage levels: a fresh node per terminal, and a synthetic
        # breaker to the node of its bus-breaker bus (open iff disconnected)
        in_bb = ~in_nb
        if bbb_col in df.columns and in_bb.any():
            bbb = df[bbb_col].astype(str).values
            has_bus = in_bb & (bbb != "") & pd.Series(bbb).isin(bb_bus_node.index).values
            if has_bus.any():
                sub = df.loc[has_bus, [vl_col]].copy()
                sub["k"] = sub.groupby(vl_col).cumcount()
                base = next_node.reindex(sub[vl_col].values).values
                new_nodes = base + sub["k"].values
                term.loc[has_bus, "node"] = new_nodes
                added = sub.groupby(vl_col).size()
                next_node[added.index] += added.values
                bus_nodes = bb_bus_node.reindex(bbb[has_bus]).values
                connected = df.loc[has_bus, conn_col].to_numpy(bool) if conn_col in df.columns else np.ones(has_bus.sum(), bool)
                names = [f"{el}@{side}" if side else str(el) for el in df.index[has_bus]]
                for vp, n1, n2, on, nm in zip(term.loc[has_bus, "vl_pos"].values, new_nodes, bus_nodes, connected, names):
                    synthetic_sw.append((int(vp), int(n1), int(n2), _BREAKER, not bool(on), True, nm, 2))
        scan.terminals[key] = term
        described = term["node"].values >= 0
        vl_of = df[vl_col].values[described]
        if described.any():
            per_vl = pd.Series(1, index=vl_of).groupby(level=0).sum()
            feeder_count[per_vl.index] += per_vl.values
            if is_branch:
                branch_count[per_vl.index] += per_vl.values
    nb_nodes = next_node.values.astype(int)

    # ---- busbar sections, in grid-wide order ----
    bbs_nb = net.get_busbar_sections(all_attributes=True)
    if sort_index:
        bbs_nb = bbs_nb.sort_index()
    bbs_nb = bbs_nb[bbs_nb["voltage_level_id"].isin(nb_vls)]
    bbs_rows = pd.DataFrame({
        "vl_pos": vl_pos.reindex(bbs_nb["voltage_level_id"].values).values,
        "node": bbs_nb["node"].to_numpy(int),
        "name": bbs_nb.index.astype(str),
        "bus": bbs_nb["bus_id"].astype(str).values,
        "seq": 0,
    })
    bbs_synth = pd.DataFrame({
        "vl_pos": vl_pos.reindex(bbv["voltage_level_id"].values).values,
        "node": bbv["node"].to_numpy(int),
        "name": bbv.index.astype(str),
        "bus": bbv["bus_id"].astype(str).values,
        "seq": 1,
    })
    bbs = pd.concat([bbs_rows, bbs_synth], ignore_index=True)
    bbs["order"] = np.arange(bbs.shape[0])
    bbs = bbs.sort_values(["vl_pos", "seq", "order"], kind="mergesort").reset_index(drop=True)
    scan.bbs = bbs[["vl_pos", "node", "name", "bus"]]

    # ---- switches, in grid-wide order ----
    sw = net.get_switches(all_attributes=True)
    if sort_index:
        sw = sw.sort_index()
    rows = []
    sw_in_nb = sw["voltage_level_id"].isin(nb_vls).values
    for sw_id, r in sw[sw_in_nb].iterrows():
        kind = _KIND_OF.get(str(r["kind"]))
        if kind is None:
            raise RuntimeError(f"init_from_pypowsybl: switch '{sw_id}' has an unknown kind {r['kind']!r}.")
        rows.append((int(vl_pos[r["voltage_level_id"]]), int(r["node1"]), int(r["node2"]), kind,
                     bool(r["open"]), bool(r["retained"]), str(sw_id), 0))
    for vp, n1, n2 in ic_rows:
        rows.append((vp, n1, n2, _INTERNAL_CONNECTION, False, False, "", 1))
    skipped = []
    for sw_id, r in sw[~sw_in_nb].iterrows():
        b1, b2 = str(r["bus_breaker_bus1_id"]), str(r["bus_breaker_bus2_id"])
        if (b1 not in bb_bus_node.index) or (b2 not in bb_bus_node.index):
            skipped.append(sw_id)
            continue
        kind = _KIND_OF.get(str(r["kind"]))
        if kind is None:
            raise RuntimeError(f"init_from_pypowsybl: switch '{sw_id}' has an unknown kind {r['kind']!r}.")
        rows.append((int(vl_pos[r["voltage_level_id"]]), int(bb_bus_node[b1]), int(bb_bus_node[b2]), kind,
                     bool(r["open"]), bool(r["retained"]), str(sw_id), 0))
    if skipped:
        warnings.warn(f"init_from_pypowsybl: {len(skipped)} switch(es) of bus-breaker voltage levels "
                      f"join a bus that is not in the bus-breaker view and are ignored: {skipped[:5]}...")
    rows.extend(synthetic_sw)
    switches = pd.DataFrame(rows, columns=["vl_pos", "node1", "node2", "kind", "open", "retained", "name", "seq"])
    switches["order"] = np.arange(switches.shape[0])
    switches = switches.sort_values(["vl_pos", "seq", "order"], kind="mergesort").reset_index(drop=True)
    scan.switches = switches[["vl_pos", "node1", "node2", "kind", "open", "retained", "name"]]

    # a switch or a terminal may reference a node past the nodes table: the
    # count is whatever the highest reference needs
    highest = pd.Series(-1, index=np.arange(n_vl))
    for frame, cols in ((scan.switches, ["node1", "node2"]), (scan.bbs, ["node"])):
        if frame.shape[0] > 0:
            m = frame.groupby("vl_pos")[cols].max().max(axis=1)
            highest[m.index] = np.maximum(highest[m.index].values, m.values)
    for term in scan.terminals.values():
        t = term[term["node"] >= 0]
        if t.shape[0] > 0:
            m = t.groupby("vl_pos")["node"].max()
            highest[m.index] = np.maximum(highest[m.index].values, m.values)
    scan.nb_nodes = np.maximum(nb_nodes, highest.values + 1).astype(int)

    # ---- the projection order of the bus-view buses, and each substation's capacity ----
    # A component holding a busbar section is numbered first, by the section's
    # position; the others by the lowest node one of their terminals stands on
    # (exactly SubstationTopology::label's rule).
    with_bbs = scan.bbs[scan.bbs["bus"] != ""]
    key_bbs = with_bbs.reset_index().groupby("bus")["index"].min()
    term_all = pd.concat([t[(t["node"] >= 0) & (t["bus"] != "")][["bus", "node"]] for t in scan.terminals.values()])
    key_term = term_all.groupby("bus")["node"].min() if term_all.shape[0] > 0 else pd.Series(dtype=int)
    all_buses = key_bbs.index.union(key_term.index)
    bus_key = pd.DataFrame(index=all_buses)
    bus_key["key_group"] = 1
    bus_key["key_idx"] = key_term.reindex(all_buses).fillna(0).astype(int)
    bus_key.loc[key_bbs.index, "key_group"] = 0
    bus_key.loc[key_bbs.index, "key_idx"] = key_bbs.values
    scan.bus_key = bus_key

    # every bus holds a busbar section or a branch terminal, and a feeder
    n_bbs_per_vl = scan.bbs.groupby("vl_pos").size().reindex(np.arange(n_vl)).fillna(0).astype(int)
    bound = np.minimum(n_bbs_per_vl.values + branch_count.values, feeder_count.values)
    scan.bound_per_vl = pd.Series(np.maximum(bound, 1), index=vls.index)
    return scan


def _aux_terminal_state(model):
    """(connected, bus) of every terminal, to check a projection against."""
    res = {}
    res["load"] = [(el.connected, el.bus_id) for el in model.get_loads()]
    res["gen"] = [(el.connected, el.bus_id) for el in model.get_generators()]
    res["shunt"] = [(el.connected, el.bus_id) for el in model.get_shunts()]
    res["svc"] = [(el.connected, el.bus_id) for el in model.get_svcs()]
    res["storage"] = [(el.connected, el.bus_id) for el in model.get_storages()]
    res["line"] = [(el.connected1, el.bus1_id, el.connected2, el.bus2_id) for el in model.get_lines()]
    res["trafo"] = [(el.connected1, el.bus1_id, el.connected2, el.bus2_id) for el in model.get_trafos()]
    res["hvdc"] = [(el.connected1, el.bus1_id, el.connected2, el.bus2_id) for el in model.get_dclines()]
    return res


def _aux_add_detailed_topology(model, scan, df_gen, df_load, df_line, df_trafo, df_shunt, df_svc,
                               df_dc, df_batt, hvdc_sub_from_id, hvdc_sub_to_id):
    """Declare the scanned topology on `model`, set every terminal's node,
    project the switches and check that nothing moved."""
    sw = scan.switches
    bbs = scan.bbs
    model.init_detailed_topology(
        np.asarray(scan.nb_nodes, dtype=np.int32),
        bbs["vl_pos"].to_numpy(np.int32), bbs["node"].to_numpy(np.int32),
        sw["vl_pos"].to_numpy(np.int32), sw["node1"].to_numpy(np.int32), sw["node2"].to_numpy(np.int32),
        [int(k) for k in sw["kind"].values],
        [bool(o) for o in sw["open"].values],
        [bool(r) for r in sw["retained"].values])
    model.set_switch_names([str(n) for n in sw["name"].values])
    model.set_busbar_section_names([str(n) for n in bbs["name"].values])

    def nodes_of(key, index):
        return scan.terminals[key]["node"].reindex(index).fillna(-1).to_numpy(np.int32)

    model.set_gen_to_node_id(nodes_of("gen", df_gen.index))
    model.set_load_to_node_id(nodes_of("load", df_load.index))
    model.set_shunt_to_node_id(nodes_of("shunt", df_shunt.index))
    model.set_svc_to_node_id(nodes_of("svc", df_svc.index))
    model.set_storage_to_node_id(nodes_of("storage", df_batt.index))
    model.set_line_to_node1_id(nodes_of("line1", df_line.index))
    model.set_line_to_node2_id(nodes_of("line2", df_line.index))
    model.set_trafo_to_node1_id(nodes_of("trafo1", df_trafo.index))
    model.set_trafo_to_node2_id(nodes_of("trafo2", df_trafo.index))
    if df_dc.shape[0] > 0:
        stations = pd.concat([scan.terminals["vsc"], scan.terminals["lcc"]])["node"]
        model.set_dcline_to_sub1_id(np.asarray(hvdc_sub_from_id, dtype=np.int32))
        model.set_dcline_to_sub2_id(np.asarray(hvdc_sub_to_id, dtype=np.int32))
        model.set_dcline_to_node1_id(stations.reindex(df_dc["converter_station1_id"].values).fillna(-1).to_numpy(np.int32))
        model.set_dcline_to_node2_id(stations.reindex(df_dc["converter_station2_id"].values).fillna(-1).to_numpy(np.int32))

    # the model's own reading of the switches must put every element exactly
    # where pypowsybl's bus view put it: that is the check of the rule
    before = _aux_terminal_state(model)
    model.project_switches()
    after = _aux_terminal_state(model)
    diffs = []
    for key in before:
        for el_id, (b, a) in enumerate(zip(before[key], after[key])):
            if b != a:
                diffs.append(f"{key} {el_id}: bus view {b} but switches {a}")
    if diffs:
        shown = "\n  ".join(diffs[:10])
        raise RuntimeError(
            f"init_from_pypowsybl: the switch positions of the grid put {len(diffs)} element(s) "
            "on a different bus than pypowsybl's bus view (connected, bus id); the loader and "
            "pypowsybl disagree on which components are buses, please report it with the grid:\n  "
            f"{shown}")
