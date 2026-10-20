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
opens), a disconnector, a load-break switch or a fictitious breaker (powsybl never opens one)
does not, a busbar section is walked through,
and any other equipment met on the way is tripped. An open switch is not a way through.
The equipment tripped that way is tripped as a whole, so the walk starts again from each of
its other terminals, until nothing new is met: a line reached at one end of a breakerless
tee takes down what hangs on its far end too, as OpenLoadFlow does.

Behind each breaker the walk opens, what the remaining closed switches (breakers included)
still join is cut off from the walked part. With no busbar section, it is disconnected when it
holds a single piece of equipment (left hanging: a transformer whose bay hangs on a busbar that
a breakerless line takes down, or each feeder of a tripped busbar section), or only injections
(nothing feeds them any more), at that end only, so it is not walked any further. A part still
holding a busbar section, or several pieces of equipment one of which is a branch, stays: it is
still fed. A bus-breaker voltage level has one breaker per terminal, so nothing propagates there.

A branch cut off that way is opened at that end only: its other end stays connected, and the
branch keeps its charging from there (OpenLoadFlow opens the side of each terminal it
disconnects, see `ContingencyPropagator.opened_sides`). Equipment the walk meets is tripped as
a whole, at every end.

A branch alone behind an opened breaker is only disconnected when its terminal sits right behind
that breaker, with no other closed switch nor internal connection on that node (powsybl's
`NodeBreakerTraverser` then disconnects the terminal instead of opening the breaker). Further
behind (another closed switch, a fictitious breaker for instance, in between), powsybl opens the
breaker only and the branch stays in service, energized from its other end and feeding a bus that
holds nothing else: it is not tripped, and `ContingencyPropagator.dead_ends` lists it.

Two more of `NodeBreakerTraverser`'s rules decide which breakers open:

* on a node with a terminal (the tripped equipment's, or a busbar section's) whose closed
  switches are all breakers, powsybl opens none of them and disconnects the terminal alone:
  what is behind them stays joined through that node. Likewise when that node holds a single
  other switch, leading to a node whose other closed switches are all breakers: that switch is
  not crossed, and those breakers stay closed;
* a breaker leading to a node with no terminal whose other closed switches are all breakers is
  walked through, and the breakers of that node open instead.
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
    one graph. Only the edges the walk can cross are kept: a closed switch that is not a breaker
    (or a fictitious one); the closed breakers are kept apart.
    The one thing no bulk table gives is the internal connections (a link between two nodes with
    no switch on it), which take one call per voltage level; see `internal_connections` of
    `ContingencyPropagator`."""

    def __init__(self, net, internal_connections=None):
        self._net = net
        self._internal_connections = internal_connections
        vls = net.get_voltage_levels(attributes=["topology_kind"])
        self._vl_ids = list(vls.index[vls["topology_kind"].astype(str) == "NODE_BREAKER"])
        vl_index = pd.Index(self._vl_ids)

        self.terminals = {}        # element id -> [node]
        self.side_of = {}          # equipment node -> the side of the terminal on it (1, 2, 3; 0 one-sided)
        self.station_to_hvdc = {}  # converter station id -> hvdc line id
        self.hvdc_stations = {}    # hvdc line id -> [station ids]
        self.tie_dangling = {}     # tie line id -> [dangling line ids]
        self.dangling_to_tie = {}  # dangling line id -> tie line id
        self.known = set()         # every id the contingency may name
        self.branches = set()      # lines, transformers, and dangling lines of a tie line

        # ---- the terminals: (element id, voltage level number, node, is a busbar section)
        parts = []
        for getter, sides in _NODE_TABLES:
            cols = [c for side in sides for c in (f"voltage_level{side}_id", f"node{side}")]
            if getter == "get_dangling_lines" and hasattr(net, "get_boundary_lines"):
                getter = "get_boundary_lines"  # the old name is deprecated
            df = getattr(net, getter)(attributes=cols)
            self.known.update(df.index)
            if len(sides) > 1:
                self.branches.update(df.index)
            for side in sides:
                code = vl_index.get_indexer(df[f"voltage_level{side}_id"].values)
                keep = code >= 0
                parts.append((df.index.values[keep], code[keep].astype(np.int64),
                              df[f"node{side}"].values[keep].astype(np.int64),
                              getter == "get_busbar_sections", int(side) if side else 0))
        sw = net.get_switches(attributes=["voltage_level_id", "node1", "node2", "kind", "open", "fictitious"])
        sw_code = vl_index.get_indexer(sw["voltage_level_id"].values)
        sw = sw[sw_code >= 0]
        sw_code = sw_code[sw_code >= 0].astype(np.int64)
        n1, n2 = sw["node1"].values.astype(np.int64), sw["node2"].values.astype(np.int64)

        highest = max([int(p[2].max()) for p in parts if p[2].size] + ([int(max(n1.max(), n2.max()))] if n1.size else [0]))
        self.stride = highest + 1
        ids = np.concatenate([p[0] for p in parts]) if parts else np.empty(0, dtype=object)
        keys = np.concatenate([p[1] * self.stride + p[2] for p in parts]) if parts else np.empty(0, np.int64)
        is_bbs = np.concatenate([np.full(p[0].shape[0], p[3]) for p in parts]) if parts else np.empty(0, bool)
        sides = np.concatenate([np.full(p[0].shape[0], p[4]) for p in parts]) if parts else np.empty(0, np.int64)
        for el, key in zip(ids.tolist(), keys.tolist()):
            self.terminals.setdefault(el, []).append(key)
        self._equip = dict(zip(keys[~is_bbs].tolist(), ids[~is_bbs].tolist()))
        self.side_of = dict(zip(keys[~is_bbs].tolist(), sides[~is_bbs].tolist()))
        self._bbs = set(keys[is_bbs].tolist())

        # ---- the edges the walk can cross
        k1, k2 = sw_code * self.stride + n1, sw_code * self.stride + n2
        closed = ~sw["open"].values.astype(bool)
        # powsybl does not open a fictitious breaker: the walk goes through it
        breaker = (sw["kind"].astype(str).values == "BREAKER") & ~sw["fictitious"].values.astype(bool)
        passes = ~breaker & closed
        u, v = np.concatenate([k1[passes], k2[passes]]), np.concatenate([k2[passes], k1[passes]])
        order = np.argsort(u, kind="stable")
        u, self._nbr = u[order], v[order].tolist()
        first, count = np.unique(u, return_index=True, return_counts=True)[1:]
        self._span = dict(zip(np.unique(u).tolist(), zip(first.tolist(), (first + count).tolist())))
        # the closed breakers: where a walk stops (and what is behind them, see `cut_off`)
        self._breakers = {}        # node -> [node]
        brk = breaker & closed
        for a, b in zip(k1[brk].tolist(), k2[brk].tolist()):
            self._breakers.setdefault(a, []).append(b)
            self._breakers.setdefault(b, []).append(a)
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
            self.branches.update((d1, d2))

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

    def _has_terminal(self, node):
        return node in self._equip or node in self._bbs

    def _nb_passing(self, node):
        """The closed switches on `node` that are not breakers (a fictitious breaker is one)."""
        lo_hi = self._span.get(node)
        return 0 if lo_hi is None else lo_hi[1] - lo_hi[0]

    def _stops_at(self, node, nb_passing, dismissed):
        """powsybl's `traverserStopsAtOtherEdges`: but for the edge looked at (`nb_passing` 1 when it
        is a closed switch that is not a breaker, else 0), every closed switch on `node` is a
        breaker, and its internal connections (`dismissed` left out) lead only to nodes with no
        terminal where the same holds."""
        if self._nb_passing(node) != nb_passing:
            return False
        visited, stack = {dismissed, node}, [node]
        while stack:
            for other in self._extra.get(stack.pop(), ()):
                if other in visited:
                    continue
                visited.add(other)
                if self._has_terminal(other) or self._nb_passing(other):
                    return False
                stack.append(other)
        return True

    def walk(self, start, tripped, seen, behind_breakers, kept):
        """Walk out of one terminal (see the module docstring); the equipment met goes to `tripped`,
        the node on the far side of each breaker opened to `behind_breakers`, and the nodes left
        joined through a node with a terminal, where the walk stops without opening anything, to
        `kept` (one list per such node).
        `seen` holds the nodes already met by the walks of the same contingency: everything
        reachable from them is already found, so they are not walked again (`start` is, since an
        equipment node is met but not walked out of)."""
        self.read_internal_connections(start // self.stride)
        span, nbr, extra, equip, breakers = self._span, self._nbr, self._extra, self._equip, self._breakers
        seen.add(start)
        stack = [start]
        while stack:
            node = stack.pop()
            at_terminal = self._has_terminal(node)  # `start`, or a busbar section
            joined = []
            for other in breakers.get(node, ()):
                if other in seen:
                    continue
                if at_terminal and self._stops_at(node, 0, other):
                    joined.append(other)  # the terminal alone is disconnected: the breaker stays closed
                elif not self._has_terminal(other) and self._stops_at(other, 0, node):
                    seen.add(other)  # a junction of breakers only: walked through, its own breakers open
                    stack.append(other)
                else:
                    behind_breakers.append(other)
            lo_hi = span.get(node)
            others = [(other, 1) for other in (nbr[lo_hi[0]:lo_hi[1]] if lo_hi is not None else ())]
            others += [(other, 0) for other in extra.get(node, ())]
            for other, nb_passing in others:
                if other in seen:
                    continue
                if (at_terminal and not self._has_terminal(other) and self._stops_at(node, nb_passing, other)
                        and self._stops_at(other, nb_passing, node)):
                    joined.append(other)  # a junction of breakers right after the terminal: not walked
                    continue
                seen.add(other)
                el = equip.get(other)
                if el is not None:
                    tripped.append(el)  # equipment is not a way through
                else:
                    stack.append(other)  # a busbar section, or a bare node
            if joined:
                kept.append(joined)

    def cut_off(self, start, seen, visited, joined, kept=False):
        """The equipment left hanging behind an opened breaker (see the module docstring): the
        equipment `start` is joined to by the closed switches, breakers included, once the walked
        nodes (`seen`) are out, when no busbar section is joined and it is a single piece of
        equipment or only injections; else nothing. `visited` holds the nodes already looked at behind another breaker
        of the same contingency, `joined` the nodes still joined to each other through a walked node
        (see `walk`), and `kept` says that `start` is one of them (no breaker opened in front of it).
        Each piece comes with the side of its terminal there (0 for a one-sided element), the only
        end it is disconnected at.
        Also says whether that is a branch left at a dead end (not right behind an opened breaker,
        see the module docstring), which stays in service."""
        if start in seen or start in visited:
            return [], False
        self.read_internal_connections(start // self.stride)
        span, nbr, extra, equip, breakers = self._span, self._nbr, self._extra, self._equip, self._breakers
        visited.add(start)
        stack = [start]
        found = []
        live = False  # a busbar section stays energized, with whatever it still holds
        while stack:
            node = stack.pop()
            live = live or node in self._bbs
            el = equip.get(node)
            if el is not None:
                found.append((el, self.side_of.get(node, 0)))
            lo_hi = span.get(node)
            others = nbr[lo_hi[0]:lo_hi[1]] if lo_hi is not None else []
            others = others + breakers.get(node, []) + extra.get(node, []) + joined.get(node, [])
            for other in others:
                if other in seen or other in visited:
                    continue
                visited.add(other)
                stack.append(other)
        # (the whole part is walked before deciding, so that another breaker into it finds it visited)
        found = list(dict.fromkeys(found))
        if live or (len({el for el, _ in found}) > 1 and any(el in self.branches for el, _ in found)):
            return [], False  # still fed
        # a branch alone is disconnected when its terminal is right behind the breaker, the only
        # closed switch on its node (NodeBreakerTraverser.isEquivalentToStopAfterSwitch)
        lo_hi = span.get(start)
        nb_closed = (lo_hi[1] - lo_hi[0] if lo_hi is not None else 0) + len(breakers.get(start, ()))
        right_behind = not kept and equip.get(start) is not None and nb_closed == 1 and start not in extra
        return found, len(found) == 1 and found[0][0] in self.branches and not right_behind


def _terminals_of(topo, name):
    """The nodes a contingency on `name` trips from (and the element ids it names)."""
    if name in topo.hvdc_stations:
        elements = list(topo.hvdc_stations[name])
    elif name in topo.tie_dangling:
        elements = list(topo.tie_dangling[name])
    else:
        elements = [name]
    return elements, [t for el in elements for t in topo.terminals.get(el, ())]


class ContingencyPropagator:
    """Propagate contingencies the way powsybl does on a node-breaker topology.

    Built once from a pypowsybl network: it reads the topology (the switches, the terminals,
    the hvdc and tie lines) at construction and keeps it, so that contingencies can be added
    and propagated in several goes without reading the network again.

    Each contingency is a list of element names (lines, transformers, generators, loads,
    shunts, batteries, static var compensators, dangling and tie lines, converter stations,
    HVDC lines, busbar sections). Propagating it gives the same names, in the same order,
    followed by the other equipment that goes down with them because it stands between the
    same breakers (see this module's docstring). A name appears once.

    Examples
    --------

    .. code-block:: python

        from lightsim2grid.network import ContingencyPropagator

        propagator = ContingencyPropagator(net)
        propagator.add_contingencies([["LINE_AB"], ["LINE_BC"], ["GEN_B", "LOAD_C"]])
        propagated = propagator.propagate_contingencies()   # one list per contingency
        unique, inverse = propagator.group_contingencies()  # the ones to simulate

    Notes
    -----
    * Names come back as pypowsybl knows them; a branch that is only opened at one end (the
      other end stands behind a breaker of its own) is listed too, and :func:`opened_sides`
      says which end. A branch left in service at a dead end is not listed: see
      :func:`dead_ends`.
    * Propagating onto one converter station lists its HVDC line, and onto one dangling line
      its tie line.
    * A name that is not an element of the network (a switch included) is an error, raised by
      :func:`add_contingencies`, since a silent pass-through would hide a typo.
    * The topology is read at construction: build a new one after changing the network.
    """

    def __init__(self, net, internal_connections=None):
        """
        Parameters
        ----------
        net:
            a pypowsybl network.
        internal_connections:
            pypowsybl has no bulk getter for the internal connections of a node-breaker voltage
            level (links between two nodes with no switch on them), only one call per voltage
            level, which is most of the running time on a large grid. ``None`` (default) makes
            that call only for the voltage levels walked where a node carrying equipment has no
            switch, the mark of an internal connection; ``True`` always makes it (exact,
            slower); ``False`` never does.
        """
        self._topo = _PypowsyblTopology(net, internal_connections)
        self._contingencies = []
        self._propagated = []  # one entry per contingency, None until propagated
        self._sides = []       # aligned with _propagated: {branch: side opened}, see opened_sides
        self._dead_ends = []   # aligned with _propagated: {branch: side of the dead end}, see dead_ends

    @property
    def contingencies(self):
        """The contingencies added so far, as given (without repeated names)."""
        return [list(el) for el in self._contingencies]

    def add_contingencies(self, contingencies):
        """Add contingencies, each a list of element names.

        Raises ``TypeError`` for a contingency given as a single string, and ``ValueError`` for
        a name that is not an element of the network; nothing is added then.
        """
        checked = []
        for contingency in contingencies:
            if isinstance(contingency, str):
                raise TypeError("add_contingencies: each contingency must be a list of names, "
                                f"got the string {contingency!r}")
            names = list(dict.fromkeys(contingency))
            for name in names:
                if name not in self._topo.known:
                    raise ValueError(f"add_contingencies: '{name}' is not a line, transformer, injection, "
                                     f"converter station, HVDC line or busbar section of this network.")
            checked.append(names)
        self._contingencies.extend(checked)
        self._propagated.extend([None] * len(checked))
        self._sides.extend([None] * len(checked))
        self._dead_ends.extend([None] * len(checked))

    def clear(self):
        """Remove every contingency (the topology read at construction is kept)."""
        self._contingencies = []
        self._propagated = []
        self._sides = []
        self._dead_ends = []

    def propagate_contingencies(self):
        """The contingencies added, propagated: one list of names per contingency, in the order
        they were added. Each is propagated once, however many times this is called.

        Returns
        -------
        list of lists of str.
        """
        self._propagate_all()
        return [list(el) for el in self._propagated]

    def opened_sides(self):
        """The branches each propagated contingency opens at one end only: one dict per
        contingency, in the order they were added, ``{branch name: side opened}`` (1 or 2, 3 for
        the third leg of a three-winding transformer). Every other branch of
        :func:`propagate_contingencies` is opened at all its ends.

        A branch cut off behind a breaker the propagation opens (a feeder alone on its piece of
        the substation) is disconnected at that end only. OpenLoadFlow keeps it connected at its
        other end, where it still draws its charging current: simulate it opened at that side
        only (*eg* ``branch_ids_side1`` / ``branch_ids_side2`` of
        ``ContingencyAnalysisCPP.add_contingency``), not tripped. A tie line is always opened as
        a whole, as OpenLoadFlow does.

        Returns
        -------
        list of dicts.
        """
        self._propagate_all()
        return [dict(el) for el in self._sides]

    def dead_ends(self):
        """The branches each propagated contingency leaves in service at a dead end: one dict per
        contingency, in the order they were added, ``{branch name: side of the dead end}``.

        Such a branch stands alone behind a breaker the propagation opens, but not right behind
        it (another closed switch is in between). powsybl then opens the breaker only: the branch
        is not disconnected, it stays energized from its other end and feeds a bus that holds
        nothing else (OpenLoadFlow monitors that bus and the branch). It is not in
        :func:`propagate_contingencies`.

        A model built on the bus view (lightsim2grid's) has that bus merged with the rest of the
        substation, which may still be fed: simulate the branch opened at that side (as for
        :func:`opened_sides`). The flows at its other end are the same; only the voltage of the
        dead-end bus and the branch's own current are not seen.

        Returns
        -------
        list of dicts.
        """
        self._propagate_all()
        return [dict(el) for el in self._dead_ends]

    def _propagate_all(self):
        for i, contingency in enumerate(self._contingencies):
            if self._propagated[i] is None:
                self._propagated[i], self._sides[i], self._dead_ends[i] = self._propagate_one(contingency)

    def group_contingencies(self):
        """Group the propagated contingencies that trip the same set of equipment, so that
        each is simulated once (see :func:`group_equivalent_contingencies`).

        Returns
        -------
        unique:
            the first propagated contingency of each distinct set, in order of first appearance.
        inverse:
            for each contingency added, the position in `unique` of its representative.
        """
        # a branch at a dead end is simulated opened at that side: it tells the sets apart too
        sides = [{**opened, **dead} for opened, dead in zip(self.opened_sides(), self.dead_ends())]
        return group_equivalent_contingencies(self.propagate_contingencies(), sides)

    def _propagate_one(self, contingency):
        topo = self._topo
        result = list(contingency)
        tripped = []
        for name in result:
            tripped.extend(_terminals_of(topo, name)[0])
        # every element tripped is walked out of from all its terminals, the ones the walks trip
        # included (appended to `tripped` while it is read), until nothing new is met
        seen, done, behind_breakers, kept = set(), set(), [], []
        pos = 0
        while pos < len(tripped):
            el = tripped[pos]
            pos += 1
            if el in done:
                continue
            done.add(el)
            # what the loss of a piece of equipment takes with it: a station its HVDC line (and
            # the other station), a dangling line its tie line (and the other dangling line)
            for whole in (topo.station_to_hvdc.get(el), topo.dangling_to_tie.get(el)):
                if whole is not None:
                    tripped.append(whole)
                    tripped.extend(_terminals_of(topo, whole)[0])
            for node in topo.terminals.get(el, ()):
                topo.walk(node, tripped, seen, behind_breakers, kept)
        # what the opened breakers leave hanging, once the walked part is complete (disconnected at
        # that end only: not walked out of, but it still takes its HVDC or tie line with it), and
        # what a walk stopped in front of without opening a breaker (still joined through that node)
        whole_ids = set(tripped)
        cut_sides = {}  # element only cut off -> the sides it is cut at
        dead = {}       # branch left at a dead end -> the sides of its dead ends
        visited = set()
        joined = {}
        for group in kept:
            for node in group:
                joined.setdefault(node, []).extend(group)
        for node, was_kept in [(node, False) for node in behind_breakers] + [(node, True) for node in joined]:
            found, dead_end = topo.cut_off(node, seen, visited, joined, was_kept)
            for el, side in found:
                (dead if dead_end else cut_sides).setdefault(el, set()).add(side)
        # a branch at a dead end at one end only, and not reached in any other way, stays in
        # service; at several ends, or also cut off at another, nothing feeds it any more
        dead_ends = {}
        for el, el_sides in dead.items():
            if el in whole_ids:
                continue
            if len(el_sides) == 1 and el not in cut_sides:
                tie = topo.dangling_to_tie.get(el)
                if tie is None:
                    dead_ends[el] = next(iter(el_sides))
                else:
                    dead_ends[tie] = topo.tie_dangling[tie].index(el) + 1
                continue
            cut_sides.setdefault(el, set()).update(el_sides)
        for el in cut_sides:
            tripped.append(el)
            for whole in (topo.station_to_hvdc.get(el), topo.dangling_to_tie.get(el)):
                if whole is not None:
                    tripped.append(whole)
                    whole_ids.add(whole)
                    tripped.extend(_terminals_of(topo, whole)[0])
        known = set(result)
        result.extend(el for el in dict.fromkeys(tripped) if el not in known)
        # a branch cut off at one of its ends, and not tripped as a whole otherwise, stays
        # connected at the other ones (a dangling line goes with its tie line, a whole)
        sides = {}
        for el, el_sides in cut_sides.items():
            if el in whole_ids or el not in topo.branches or el in topo.dangling_to_tie:
                continue
            if len(el_sides) == 1:
                sides[el] = next(iter(el_sides))
        # (a tie line at a dead end at one end may still be tripped through its other one)
        listed = set(result)
        return result, sides, {el: side for el, side in dead_ends.items() if el not in listed}


def propagate_contingencies(net, contingencies, internal_connections=None):
    """Propagate contingencies in one call: a shortcut for

    .. code-block:: python

        propagator = ContingencyPropagator(net, internal_connections)
        propagator.add_contingencies(contingencies)
        propagator.propagate_contingencies()

    See :class:`ContingencyPropagator` for the parameters and the rule; build one directly to
    add contingencies in several goes, or to group them, without reading the network again.

    Returns
    -------
    list of lists of str.
    """
    propagator = ContingencyPropagator(net, internal_connections)
    propagator.add_contingencies(contingencies)
    return propagator.propagate_contingencies()


def group_equivalent_contingencies(contingencies, opened_sides=None):
    """Find the contingencies that trip the same set of equipment, so each is simulated once.

    Meant for the output of :func:`ContingencyPropagator.propagate_contingencies`: tripping a line, or the load hanging
    on its feeder, ends in the same grid. Two contingencies are equivalent when they hold the same
    names, whatever the order, and (`opened_sides`, one dict per contingency as
    :func:`ContingencyPropagator.opened_sides` gives, optional) open the same branches at the
    same single end.

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
    for pos, contingency in enumerate(contingencies):
        key = frozenset(contingency)
        if opened_sides is not None and opened_sides[pos]:
            key = (key, frozenset(opened_sides[pos].items()))
        pos = seen.get(key)
        if pos is None:
            pos = seen[key] = len(unique)
            unique.append(list(contingency))
        inverse.append(pos)
    return unique, inverse
