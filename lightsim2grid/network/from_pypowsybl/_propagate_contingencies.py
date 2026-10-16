# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Contingency propagation, the way powsybl does it on a node-breaker voltage level.

Tripping an element does not only disconnect that element: the protection opens the
first breaker on the way out of each of its terminals, and everything that hangs on the
same piece of the substation, without a breaker in between, goes down with it.

The rule is powsybl's own (`TrippingTopologyTraverser`): from the node of a tripped
terminal, walk the closed switches; a breaker stops the walk (and is the switch that
opens), a disconnector or a load-break switch does not, a busbar section is walked through,
and any other equipment met on the way is tripped and not walked through. An open switch
is not a way through. A bus-breaker voltage level has one breaker per terminal, so nothing
propagates there.
"""

import numpy as np
import pandas as pd

# element tables whose terminals sit on a node: (getter, side suffixes)
_NODE_TABLES = [
    ("get_generators", ("",)),
    ("get_loads", ("",)),
    ("get_shunt_compensators", ("",)),
    ("get_static_var_compensators", ("",)),
    ("get_batteries", ("",)),
    ("get_dangling_lines", ("",)),
    ("get_lines", ("1", "2")),
    ("get_2_windings_transformers", ("1", "2")),
    ("get_3_windings_transformers", ("1", "2", "3")),
    ("get_vsc_converter_stations", ("",)),
    ("get_lcc_converter_stations", ("",)),
    ("get_busbar_sections", ("",)),
]


class _PypowsyblTopology:
    """What the walk needs of a network, read from a few bulk tables (only the columns it needs).

    A node is an integer ``voltage level number * stride + node``, so that the whole network is
    one graph. Only the edges the walk can cross are kept: a closed switch that is not a breaker.
    The one thing no bulk table gives is the internal connections (a link between two nodes with
    no switch on it), which take one call per voltage level; see `internal_connections` of
    `propagate_contingencies`."""

    def __init__(self, net, internal_connections=None):
        self._net = net
        self._internal_connections = internal_connections
        vls = net.get_voltage_levels(attributes=["topology_kind"])
        self._vl_ids = list(vls.index[vls["topology_kind"].astype(str) == "NODE_BREAKER"])
        vl_index = pd.Index(self._vl_ids)

        self.terminals = {}        # element id -> [node]
        self.station_to_hvdc = {}  # converter station id -> hvdc line id
        self.hvdc_stations = {}    # hvdc line id -> [station ids]
        self.tie_dangling = {}     # tie line id -> [dangling line ids]
        self.dangling_to_tie = {}  # dangling line id -> tie line id
        self.known = set()         # every id the contingency may name

        # ---- the terminals: (element id, voltage level number, node, is a busbar section)
        parts = []
        for getter, sides in _NODE_TABLES:
            cols = [c for side in sides for c in (f"voltage_level{side}_id", f"node{side}")]
            if getter == "get_dangling_lines" and hasattr(net, "get_boundary_lines"):
                getter = "get_boundary_lines"  # the old name is deprecated
            df = getattr(net, getter)(attributes=cols)
            self.known.update(df.index)
            for side in sides:
                code = vl_index.get_indexer(df[f"voltage_level{side}_id"].values)
                keep = code >= 0
                parts.append((df.index.values[keep], code[keep].astype(np.int64),
                              df[f"node{side}"].values[keep].astype(np.int64),
                              getter == "get_busbar_sections"))
        sw = net.get_switches(attributes=["voltage_level_id", "node1", "node2", "kind", "open"])
        sw_code = vl_index.get_indexer(sw["voltage_level_id"].values)
        sw = sw[sw_code >= 0]
        sw_code = sw_code[sw_code >= 0].astype(np.int64)
        n1, n2 = sw["node1"].values.astype(np.int64), sw["node2"].values.astype(np.int64)

        highest = max([int(p[2].max()) for p in parts if p[2].size] + ([int(max(n1.max(), n2.max()))] if n1.size else [0]))
        self.stride = highest + 1
        ids = np.concatenate([p[0] for p in parts]) if parts else np.empty(0, dtype=object)
        keys = np.concatenate([p[1] * self.stride + p[2] for p in parts]) if parts else np.empty(0, np.int64)
        is_bbs = np.concatenate([np.full(p[0].shape[0], p[3]) for p in parts]) if parts else np.empty(0, bool)
        for el, key in zip(ids.tolist(), keys.tolist()):
            self.terminals.setdefault(el, []).append(key)
        self._equip = dict(zip(keys[~is_bbs].tolist(), ids[~is_bbs].tolist()))
        self._bbs = set(keys[is_bbs].tolist())

        # ---- the edges the walk can cross
        k1, k2 = sw_code * self.stride + n1, sw_code * self.stride + n2
        passes = (sw["kind"].astype(str).values != "BREAKER") & ~sw["open"].values.astype(bool)
        u, v = np.concatenate([k1[passes], k2[passes]]), np.concatenate([k2[passes], k1[passes]])
        order = np.argsort(u, kind="stable")
        u, self._nbr = u[order], v[order].tolist()
        first, count = np.unique(u, return_index=True, return_counts=True)[1:]
        self._span = dict(zip(np.unique(u).tolist(), zip(first.tolist(), (first + count).tolist())))
        # a node carrying equipment that no switch touches: only an internal connection can join it
        unswitched = ~np.isin(keys, np.concatenate([k1, k2]))
        self._suspect = set((keys[unswitched] // self.stride).tolist())
        self._extra = {}           # node -> [node], the internal connections
        self._ic_read = set()

        hvdc = net.get_hvdc_lines(attributes=["converter_station1_id", "converter_station2_id"])
        for line, s1, s2 in zip(hvdc.index, hvdc["converter_station1_id"], hvdc["converter_station2_id"]):
            self.known.add(line)
            self.hvdc_stations[line] = [s1, s2]
            self.station_to_hvdc[s1] = line
            self.station_to_hvdc[s2] = line
        ties = net.get_tie_lines(attributes=["dangling_line1_id", "dangling_line2_id"])
        for line, d1, d2 in zip(ties.index, ties["dangling_line1_id"], ties["dangling_line2_id"]):
            self.known.add(line)
            self.tie_dangling[line] = [d1, d2]
            self.dangling_to_tie[d1] = line
            self.dangling_to_tie[d2] = line

    def read_internal_connections(self, code):
        """Make the internal connections of one voltage level part of the graph (once)."""
        if code in self._ic_read:
            return
        self._ic_read.add(code)
        wanted = self._internal_connections
        if wanted is None:
            wanted = code in self._suspect
        if not wanted:
            return
        ic = self._net.get_node_breaker_topology(self._vl_ids[code]).internal_connections
        base = code * self.stride
        for a, b in zip(ic["node1"].values.tolist(), ic["node2"].values.tolist()):
            self._extra.setdefault(base + a, []).append(base + b)
            self._extra.setdefault(base + b, []).append(base + a)

    def walk(self, start, tripped):
        """Walk out of one terminal (see the module docstring); the equipment met goes to `tripped`."""
        self.read_internal_connections(start // self.stride)
        span, nbr, extra, equip, bbs = self._span, self._nbr, self._extra, self._equip, self._bbs
        seen = {start}
        stack = [start]
        while stack:
            node = stack.pop()
            lo_hi = span.get(node)
            others = nbr[lo_hi[0]:lo_hi[1]] if lo_hi is not None else []
            if node in extra:
                others = others + extra[node]
            for other in others:
                if other in seen:
                    continue
                seen.add(other)
                el = equip.get(other)
                if el is not None:
                    tripped.append(el)  # equipment is not a way through
                else:
                    stack.append(other)  # a busbar section, or a bare node


def _terminals_of(topo, name):
    """The nodes a contingency on `name` trips from (and the element ids it names)."""
    if name in topo.hvdc_stations:
        elements = list(topo.hvdc_stations[name])
    elif name in topo.tie_dangling:
        elements = list(topo.tie_dangling[name])
    else:
        elements = [name]
    return elements, [t for el in elements for t in topo.terminals.get(el, ())]


def propagate_contingencies(net, contingencies, internal_connections=None):
    """Propagate contingencies the way powsybl does on a node-breaker topology.

    Each contingency is a list of element names (lines, transformers, generators, loads,
    shunts, batteries, static var compensators, dangling and tie lines, converter stations,
    HVDC lines, busbar sections). The result has one list per contingency: the same names,
    in the same order, followed by the other equipment that goes down with them because it
    stands between the same breakers (see this module's docstring). A name appears once.

    Parameters
    ----------
    net:
        a pypowsybl network.
    contingencies:
        a list of lists of element names.
    internal_connections:
        pypowsybl has no bulk getter for the internal connections of a node-breaker voltage
        level (links between two nodes with no switch on them), only one call per voltage
        level, which is most of the running time on a large grid. ``None`` (default) makes that
        call only for the voltage levels walked where a node carrying equipment has no switch,
        the mark of an internal connection; ``True`` always makes it (exact, slower); ``False``
        never does.

    Returns
    -------
    list of lists of str.
        See :func:`group_equivalent_contingencies` to find the ones that need simulating once.

    Notes
    -----
    * Names come back as pypowsybl knows them; a line that is only tripped on one side (the
      other side stands behind a breaker of its own) is listed, as a name cannot say "one side".
    * Propagating onto one converter station lists its HVDC line, and onto one dangling line
      its tie line.
    * A name that is not an element of the network (a switch included) is an error, since a
      silent pass-through would hide a typo.
    """
    topo = _PypowsyblTopology(net, internal_connections)
    out = []
    for contingency in contingencies:
        if isinstance(contingency, str):
            raise TypeError("propagate_contingencies: each contingency must be a list of names, "
                            f"got the string {contingency!r}")
        result = list(dict.fromkeys(contingency))
        for name in result:
            if name not in topo.known:
                raise ValueError(f"propagate_contingencies: '{name}' is not a line, transformer, injection, "
                                 f"converter station, HVDC line or busbar section of this network.")
        tripped = []
        for name in result:
            elements, nodes = _terminals_of(topo, name)
            tripped.extend(elements)
            for node in nodes:
                topo.walk(node, tripped)
        # what the loss of a piece of equipment takes with it: a station its HVDC line (and the
        # other station), a dangling line its tie line (and the other dangling line)
        for el in list(tripped):
            for whole in (topo.station_to_hvdc.get(el), topo.dangling_to_tie.get(el)):
                if whole is not None:
                    tripped.append(whole)
                    tripped.extend(_terminals_of(topo, whole)[0])
        known = set(result)
        result.extend(el for el in dict.fromkeys(tripped) if el not in known)
        out.append(result)
    return out


def group_equivalent_contingencies(contingencies):
    """Find the contingencies that trip the same set of equipment, so each is simulated once.

    Meant for the output of :func:`propagate_contingencies`: tripping a line, or the load hanging
    on its feeder, ends in the same grid. Two contingencies are equivalent when they hold the same
    names, whatever the order.

    Returns
    -------
    unique:
        the first contingency of each distinct set, in order of first appearance.
    inverse:
        for each input contingency, the position in `unique` of its representative, so that
        ``results = [unique_results[i] for i in inverse]`` expands the simulated results back.
    """
    seen = {}
    unique = []
    inverse = []
    for contingency in contingencies:
        key = frozenset(contingency)
        pos = seen.get(key)
        if pos is None:
            pos = seen[key] = len(unique)
            unique.append(list(contingency))
        inverse.append(pos)
    return unique, inverse
