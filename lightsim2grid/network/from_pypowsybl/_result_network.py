# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Cast a converged :class:`LSGrid` (built by ``init_from_pypowsybl``) back into a
pypowsybl-``Network``-shaped result view: :class:`LightsimResultNetwork` exposes
``get_buses`` / ``get_lines`` / ``get_generators`` / ... methods that return
pandas DataFrames with the same index (pypowsybl element ids) and the same
result column names/units as the real ``pypowsybl.network.Network``, so
analysis code written against a solved pypowsybl network runs unmodified
against a solved lightsim2grid one.

Only supports grids built by :func:`_from_pypowsybl.init` (``init_from_pypowsybl``):
it relies on that function's invariants -- every non-bus element's ``.name`` set
verbatim to its pypowsybl id, and the grid's ``_init_kwargs``/``_orig_to_ls``
properties -- which do not hold for pandapower- or powermodels-built grids.

Sign conventions (lightsim2grid vs pypowsybl's post-solve ``p``/``q``, both in
the "load convention": positive = power flowing *into* the equipment):
loads, shunts and storage units already match (no flip needed: `_from_pypowsybl`
feeds storages a load-convention target_p/target_q, negating IIDM's own
generator-convention target_p on the way in). Generators and HVDC converter
stations use lightsim2grid's generation convention internally (positive =
production), the opposite of pypowsybl's post-solve columns, so those are
negated here -- confirmed against ``net.get_generators()["p"] == -prod_p``
(see ``lightsim2grid/tests/test_LSGrid_pypowsybl.py``) and, for HVDC converter
stations, by direct comparison against an OLF-solved
``create_four_substations_node_breaker_network()``. SVC is assumed to follow
the same generation convention as generators (it shares the same internal
voltage-control machinery and is fed target_q/target_v un-negated, exactly
like a generator, unlike storages) but this could not be independently
double-checked against a converged real grid.
"""

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import pypowsybl as pypo

from ...lightsim2grid_cpp import LSGrid  # type: ignore


class LightsimResultNetwork:
    """pypowsybl-``Network``-shaped view of a solved lightsim2grid ``LSGrid``.

    :param ls_grid: a grid built by ``init_from_pypowsybl(net, ...)`` and
        already solved (``ac_pf``/``dc_pf`` converged).
    :param net: the same pypowsybl network passed to that ``init()`` call.

    Every ``get_*`` method mirrors its ``pypo.network.Network`` namesake:
    it accepts an optional ``attributes`` list and returns a DataFrame
    indexed by the pypowsybl element id, built lazily on first call and
    cached afterwards.
    """

    def __init__(self, ls_grid: LSGrid, net: "pypo.network.Network"):
        self._grid = ls_grid
        self._net = net

        init_kwargs = ls_grid._init_kwargs
        self._sort_index = init_kwargs.get("sort_index", "True") == "True"
        if init_kwargs.get("buses_for_sub") == "True":
            raise NotImplementedError(
                "LightsimResultNetwork does not support a grid built with "
                "buses_for_sub=True (legacy 'substation = pypowsybl bus' mode)."
            )

        self._cache: Dict[str, pd.DataFrame] = {}
        self._bus_id_lookup: Optional[np.ndarray] = None  # see _ls_bus_to_pypo
        self._sub_names: Optional[List[str]] = None  # see _ls_sub_to_vl
        self._fused_ids_cache = None  # see _fused_ids

    def _bus_df(self) -> pd.DataFrame:
        net = self._net
        return net.get_buses().sort_index() if self._sort_index else net.get_buses()

    @staticmethod
    def _maybe_select(df: pd.DataFrame, attributes: Optional[List[str]]) -> pd.DataFrame:
        return df[attributes] if attributes is not None else df

    @staticmethod
    def _records_to_frame(records: List[dict], columns: List[str]) -> pd.DataFrame:
        # `columns=` must be passed explicitly: `pd.DataFrame(records)` alone
        # infers columns from the (possibly empty) records list, so an element
        # type with 0 instances (eg no battery in this grid) would otherwise
        # come back with no "id" column at all instead of an empty frame.
        return pd.DataFrame(records, columns=columns).set_index("id")

    def _ls_bus_to_pypo(self, ls_bus_id: int) -> Optional[str]:
        """lightsim2grid internal bus id (eg `GenInfo.bus_id`) -> pypowsybl bus id.

        Every ``*Info.bus_id``/``bus1_id``/``bus2_id`` field is lightsim2grid's own
        integer bus numbering, *not* the pypowsybl string id despite the shared
        name -- this reverses it through ``_ls_to_orig`` (opposite of
        ``_orig_to_ls``, used for :meth:`get_buses`) into ``_bus_df()``'s index.
        Returns ``None`` for a disconnected element with no bus (``ls_bus_id < 0``)
        or an ``_ls_to_orig`` slot with no counterpart in ``_bus_df()`` (eg a
        dangling-line boundary bus, see :meth:`_build_buses`).
        """
        if self._bus_id_lookup is None:
            bus_df = self._bus_df()
            n_bus = bus_df.shape[0]
            ls_to_orig = np.asarray(self._grid._ls_to_orig)
            lookup = np.full(ls_to_orig.shape[0], None, dtype=object)
            valid = (ls_to_orig >= 0) & (ls_to_orig < n_bus)
            lookup[valid] = bus_df.index.to_numpy()[ls_to_orig[valid]]
            self._bus_id_lookup = lookup
        if ls_bus_id < 0 or ls_bus_id >= self._bus_id_lookup.shape[0]:
            return None
        return self._bus_id_lookup[ls_bus_id]

    def _ls_sub_to_vl(self, sub_id: int) -> Optional[str]:
        """lightsim2grid internal substation id -> pypowsybl voltage_level id.

        Every ``*Info.voltage_level_id``/``voltage_level1_id``/``voltage_level2_id``
        field is bound to the same integer as ``sub_id`` (see
        ``src/bindings/python/binding_containers.cpp``), not the pypowsybl string
        id. Since the grid was built with ``buses_for_sub`` not ``True`` (enforced
        in ``__init__``), a lightsim2grid substation *is* a pypowsybl voltage
        level, and ``LSGrid.get_substations()[k].name`` is exactly the string
        `_from_pypowsybl.init()` set it to (``model.set_substation_names(...)``).
        """
        if self._sub_names is None:
            self._sub_names = [s.name for s in self._grid.get_substations()]
        if sub_id < 0 or sub_id >= len(self._sub_names):
            return None
        return self._sub_names[sub_id]

    # ------------------------------------------------------------------ #
    # buses
    # ------------------------------------------------------------------ #
    def get_buses(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "buses" not in self._cache:
            self._cache["buses"] = self._build_buses()
        return self._maybe_select(self._cache["buses"], attributes)

    def _build_buses(self) -> pd.DataFrame:
        grid = self._grid
        bus_df = self._bus_df()
        n_bus = bus_df.shape[0]

        # `_orig_to_ls` may be longer than `bus_df` (eg extra fictitious buses
        # appended for dangling-line boundaries when the grid was built with
        # `convert_dangling_lines=True`): those extra slots have no row in
        # `net.get_buses()`, so they are dropped here -- same guard as
        # `_olf_compare.py::lightsim_bus_to_iidm`.
        orig_to_ls = np.asarray(grid._orig_to_ls)[:n_bus]

        # a bus fused away by `fuse_zero_impedance_branches` (see `_from_pypowsybl.py`)
        # keeps its own row here (fixed-size bus containers) but lost every element to
        # its representative, so it is disconnected with no solved voltage of its own
        # -- redirect it to the representative bus, which is electrically the same node
        # and carries the real solved voltage. `_bus_fusion_rep` is empty for a grid
        # built without fusion (or not through `init_from_pypowsybl` at all), in which
        # case this is a no-op.
        fusion_rep = np.asarray(grid._bus_fusion_rep)
        orig_to_ls_v = fusion_rep[orig_to_ls] if fusion_rep.shape[0] else orig_to_ls

        v_cplx = np.asarray(grid.get_V())[orig_to_ls_v]
        vn_kv = np.asarray(grid.get_bus_vn_kv())[orig_to_ls]

        res = pd.DataFrame(index=bus_df.index)
        res["v_mag"] = np.abs(v_cplx) * vn_kv
        res["v_angle"] = np.degrees(np.angle(v_cplx))
        res["voltage_level_id"] = bus_df["voltage_level_id"].values

        # a uniform angle offset across every bus is just a difference of
        # reference-datum convention (lightsim2grid's slack bus needs not be
        # at the same angle-datum as whatever `net` last held), not physical:
        # remove it by aligning medians, mirroring the "offset-removed" angle
        # metric already used by `_olf_compare.py::compare_baked`.
        orig_angle = self._net.get_buses()["v_angle"]
        mask_orig = np.isfinite(orig_angle.to_numpy()) & (orig_angle.abs() > 1e-6)
        mask_new = np.isfinite(res["v_angle"].to_numpy()) & (res["v_angle"].abs() > 1e-6)
        if mask_orig.any() and mask_new.any():
            offset = res.loc[mask_new, "v_angle"].median() - orig_angle.loc[mask_orig].median()
            res["v_angle"] = res["v_angle"] - offset

        return res

    # ------------------------------------------------------------------ #
    # lines / transformers
    # ------------------------------------------------------------------ #
    def get_lines(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        self._ensure_lines_trafos()
        return self._maybe_select(self._cache["lines"], attributes)

    def get_2_windings_transformers(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        self._ensure_lines_trafos()
        return self._maybe_select(self._cache["trafos"], attributes)

    _TWO_SIDED_COLUMNS = ["id", "p1", "q1", "i1", "p2", "q2", "i2",
                          "bus1_id", "bus2_id", "connected1", "connected2",
                          "voltage_level1_id", "voltage_level2_id"]

    def _build_two_sided(self, container) -> pd.DataFrame:
        records = []
        for el in container:
            records.append({
                "id": el.name,
                "p1": el.res_p1_mw, "q1": el.res_q1_mvar, "i1": el.res_a1_ka * 1000.,
                "p2": el.res_p2_mw, "q2": el.res_q2_mvar, "i2": el.res_a2_ka * 1000.,
                "bus1_id": self._ls_bus_to_pypo(el.bus1_id), "bus2_id": self._ls_bus_to_pypo(el.bus2_id),
                "connected1": el.connected1, "connected2": el.connected2,
                "voltage_level1_id": self._ls_sub_to_vl(el.voltage_level1_id),
                "voltage_level2_id": self._ls_sub_to_vl(el.voltage_level2_id),
            })
        return self._records_to_frame(records, self._TWO_SIDED_COLUMNS)

    def _ensure_lines_trafos(self) -> None:
        """Build (and cache) `lines`/`trafos` together, then patch in the flow
        through any *fused* (near-zero-impedance) branch where recoverable --
        see :meth:`_reconstruct_fused_branches`. Built jointly (rather than each
        lazily on its own first call) because reconstructing a fused line may
        need an already-built *trafos* frame to sum up a shared bus's other
        elements, and vice versa; done as a raw build + a patch pass, rather than
        interleaving, to avoid the two triggering each other recursively.
        """
        if "lines" in self._cache and "trafos" in self._cache:
            return
        lines_df = self._build_two_sided(self._grid.get_lines())
        trafos_df = self._build_two_sided(self._grid.get_trafos())
        fused_line_ids, fused_trafo_ids = self._fused_ids()
        if fused_line_ids or fused_trafo_ids:
            self._reconstruct_fused_branches(lines_df, trafos_df, fused_line_ids, fused_trafo_ids)
        self._cache["lines"] = lines_df
        self._cache["trafos"] = trafos_df

    def _fused_ids(self):
        """(fused_line_ids, fused_trafo_ids): ids of the branches
        `fuse_zero_impedance_branches` fused away when the grid was built (see
        `_from_pypowsybl.py`), stashed by that function into `_init_kwargs` as
        `"\\x1f"`-joined strings. Empty frozensets for a grid built without
        fusion (or not through `init_from_pypowsybl` at all).
        """
        if self._fused_ids_cache is None:
            kwargs = self._grid._init_kwargs
            line_ids = frozenset(kwargs.get("fused_line_ids", "").split("\x1f")) - {""}
            trafo_ids = frozenset(kwargs.get("fused_trafo_ids", "").split("\x1f")) - {""}
            self._fused_ids_cache = (line_ids, trafo_ids)
        return self._fused_ids_cache

    def _bus_balance(self, bus_id, one_sided_topo, two_sided_topo, exclude):
        """Sum, in the common load-convention sign (positive = power flows
        *from* the bus *into* the element -- matching every `get_*` method here,
        see the module docstring), of every element attached to `bus_id`,
        excluding the branch id `exclude` on either of its sides. By Kirchhoff's
        current law this equals the negative of `exclude`'s own flow at this
        bus once every *other* element is accounted for -- see
        :meth:`_reconstruct_fused_branches`.

        ``one_sided_topo``/``two_sided_topo`` (built once by
        :meth:`_reconstruct_fused_branches`, not per bus) pair each element
        type's *true, pre-fusion* bus assignment -- read from ``self._net``'s
        own static topology -- with this view's own already-computed result
        frame. The true assignment must come from ``self._net``, NOT from this
        view's own ``bus_id``/``bus1_id``/``bus2_id`` columns: an element
        originally on a bus fused away by `fuse_zero_impedance_branches`
        reports the *fusion representative's* pypowsybl id there instead (see
        `_ls_bus_to_pypo`), which would silently hide it from this sum.
        """
        p = q = 0.0
        for orig_bus, df in one_sided_topo:
            common = df.index.intersection(orig_bus.index[orig_bus == bus_id])
            if len(common):
                p += df.loc[common, "p"].sum()
                q += df.loc[common, "q"].sum()
        for orig_bus1, orig_bus2, df in two_sided_topo:
            c1 = df.index.intersection(orig_bus1.index[(orig_bus1 == bus_id) & (orig_bus1.index != exclude)])
            c2 = df.index.intersection(orig_bus2.index[(orig_bus2 == bus_id) & (orig_bus2.index != exclude)])
            if len(c1):
                p += df.loc[c1, "p1"].sum()
                q += df.loc[c1, "q1"].sum()
            if len(c2):
                p += df.loc[c2, "p2"].sum()
                q += df.loc[c2, "q2"].sum()
        return p, q

    def _reconstruct_fused_branches(self, lines_df, trafos_df, fused_line_ids, fused_trafo_ids) -> None:
        """Recover the flow through a fused (near-)zero-impedance branch (see
        `_from_pypowsybl.py`'s `fuse_zero_impedance_branches`) wherever the
        physics make it unambiguous, and patch `lines_df`/`trafos_df` in place.

        A fused branch is deactivated in the lightsim2grid model (its two
        terminal buses were merged into one, see `_build_buses`), so
        `lines_df`/`trafos_df` initially carry it as disconnected with 0 flow --
        correct for the *reduced* network lightsim2grid actually solved, but not
        what a caller comparing against the original (unfused) topology expects.

        The true flow is recoverable by Kirchhoff's current law at an original
        endpoint bus that has EXACTLY ONE fused branch attached (a "leaf" of the
        fused sub-graph): every other element at that bus already has its
        correct flow (see :meth:`_bus_balance`), and their sum must be exactly
        zero, so the fused branch's own flow at that side is minus that sum.

        A bus where 2+ fused branches meet (an internal node of a longer fused
        chain/star, eg two zero-impedance lines in series) has no such single
        equation -- how the combined flow splits between those branches is
        genuinely indeterminate from the solved state alone, so it is left
        as-is (disconnected, 0 flow) rather than guessed at.

        Since these branches are (near-)zero-impedance by construction (the
        fusion precondition), a side reconstructed this way is mirrored,
        lossless, to the other side (p/q negated, current magnitude unchanged)
        when that other side is not itself independently reconstructable.
        """
        orig_lines = self._net.get_lines(
            attributes=["bus1_id", "bus2_id", "voltage_level1_id", "voltage_level2_id"]
        )
        orig_trafos = self._net.get_2_windings_transformers(
            attributes=["bus1_id", "bus2_id", "voltage_level1_id", "voltage_level2_id"]
        )

        # degree, in the fused sub-graph, of every original bus touched by ANY
        # fused branch (lines and trafos together)
        touches = []
        for ids, orig in ((fused_line_ids, orig_lines), (fused_trafo_ids, orig_trafos)):
            if not ids:
                continue
            sub = orig.loc[orig.index.intersection(list(ids))]
            touches.append(sub["bus1_id"])
            touches.append(sub["bus2_id"])
        if not touches:
            return
        degree = pd.concat(touches).value_counts()

        bus_v = self.get_buses(attributes=["v_mag"])["v_mag"]

        # true (pre-fusion) bus assignment of every element, straight from the
        # original network's own static topology -- built once here, not per bus,
        # see :meth:`_bus_balance` for why this must be `self._net`, not this
        # view's own bus_id/bus1_id/bus2_id columns.
        net = self._net
        one_sided_topo = [
            (net.get_loads(attributes=["bus_id"])["bus_id"], self.get_loads()),
            (net.get_generators(attributes=["bus_id"])["bus_id"], self.get_generators()),
            (net.get_shunt_compensators(attributes=["bus_id"])["bus_id"], self.get_shunt_compensators()),
            (net.get_static_var_compensators(attributes=["bus_id"])["bus_id"], self.get_static_var_compensators()),
            (net.get_batteries(attributes=["bus_id"])["bus_id"], self.get_batteries()),
            (net.get_vsc_converter_stations(attributes=["bus_id"])["bus_id"], self.get_vsc_converter_stations()),
            (net.get_lcc_converter_stations(attributes=["bus_id"])["bus_id"], self.get_lcc_converter_stations()),
        ]
        two_sided_topo = [
            (orig_lines["bus1_id"], orig_lines["bus2_id"], lines_df),
            (orig_trafos["bus1_id"], orig_trafos["bus2_id"], trafos_df),
        ]

        for ids, orig, df in ((fused_line_ids, orig_lines, lines_df),
                              (fused_trafo_ids, orig_trafos, trafos_df)):
            for br_id in ids:
                if br_id not in orig.index:
                    continue
                bus1, bus2 = orig.loc[br_id, "bus1_id"], orig.loc[br_id, "bus2_id"]
                vl1, vl2 = orig.loc[br_id, "voltage_level1_id"], orig.loc[br_id, "voltage_level2_id"]
                leaf1 = degree.get(bus1, 0) == 1
                leaf2 = degree.get(bus2, 0) == 1
                if not leaf1 and not leaf2:
                    continue  # genuinely indeterminate, leave as-is

                p1 = q1 = p2 = q2 = None
                if leaf1:
                    other_p, other_q = self._bus_balance(bus1, one_sided_topo, two_sided_topo, exclude=br_id)
                    p1, q1 = -other_p, -other_q
                if leaf2:
                    other_p, other_q = self._bus_balance(bus2, one_sided_topo, two_sided_topo, exclude=br_id)
                    p2, q2 = -other_p, -other_q
                if p1 is None:
                    p1, q1 = -p2, -q2
                if p2 is None:
                    p2, q2 = -p1, -q1

                v1 = bus_v.get(bus1, np.nan)
                v2 = bus_v.get(bus2, np.nan)
                i1 = np.hypot(p1, q1) / (np.sqrt(3) * v1) * 1000. if v1 else np.nan
                i2 = np.hypot(p2, q2) / (np.sqrt(3) * v2) * 1000. if v2 else np.nan

                df.loc[br_id, "p1"], df.loc[br_id, "q1"], df.loc[br_id, "i1"] = p1, q1, i1
                df.loc[br_id, "p2"], df.loc[br_id, "q2"], df.loc[br_id, "i2"] = p2, q2, i2
                df.loc[br_id, "bus1_id"], df.loc[br_id, "bus2_id"] = bus1, bus2
                df.loc[br_id, "voltage_level1_id"], df.loc[br_id, "voltage_level2_id"] = vl1, vl2
                df.loc[br_id, "connected1"], df.loc[br_id, "connected2"] = True, True

    # ------------------------------------------------------------------ #
    # generators / loads / shunts / svc / batteries: one-sided elements
    # ------------------------------------------------------------------ #
    def get_generators(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "generators" not in self._cache:
            self._cache["generators"] = self._build_one_sided(self._grid.get_generators(), flip_sign=True)
        return self._maybe_select(self._cache["generators"], attributes)

    def get_loads(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "loads" not in self._cache:
            self._cache["loads"] = self._build_one_sided(self._grid.get_loads(), flip_sign=False)
        return self._maybe_select(self._cache["loads"], attributes)

    def get_shunt_compensators(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "shunts" not in self._cache:
            self._cache["shunts"] = self._build_one_sided(self._grid.get_shunts(), flip_sign=False)
        return self._maybe_select(self._cache["shunts"], attributes)

    def get_static_var_compensators(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "svcs" not in self._cache:
            self._cache["svcs"] = self._build_one_sided(self._grid.get_svcs(), flip_sign=True)
        return self._maybe_select(self._cache["svcs"], attributes)

    def get_batteries(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "batteries" not in self._cache:
            self._cache["batteries"] = self._build_one_sided(self._grid.get_storages(), flip_sign=False)
        return self._maybe_select(self._cache["batteries"], attributes)

    _ONE_SIDED_COLUMNS = ["id", "p", "q", "bus_id", "connected", "voltage_level_id"]

    def _build_one_sided(self, container, flip_sign: bool) -> pd.DataFrame:
        sign = -1. if flip_sign else 1.
        records = []
        for el in container:
            records.append({
                "id": el.name,
                "p": sign * el.res_p_mw, "q": sign * el.res_q_mvar,
                "bus_id": self._ls_bus_to_pypo(el.bus_id), "connected": el.connected,
                "voltage_level_id": self._ls_sub_to_vl(el.voltage_level_id),
            })
        return self._records_to_frame(records, self._ONE_SIDED_COLUMNS)

    # ------------------------------------------------------------------ #
    # hvdc lines / converter stations
    # ------------------------------------------------------------------ #
    def get_hvdc_lines(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "hvdc_lines" not in self._cache:
            self._build_hvdc()
        return self._maybe_select(self._cache["hvdc_lines"], attributes)

    def get_vsc_converter_stations(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "vsc_stations" not in self._cache:
            self._build_hvdc()
        return self._maybe_select(self._cache["vsc_stations"], attributes)

    def get_lcc_converter_stations(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "lcc_stations" not in self._cache:
            self._build_hvdc()
        return self._maybe_select(self._cache["lcc_stations"], attributes)

    def _build_hvdc(self) -> None:
        # `HvdcLineInfo.station1`/`.station2` never carry the pypowsybl
        # `converter_station1_id`/`converter_station2_id` string (lightsim2grid
        # only stores them as numeric bus references), so they are recovered
        # here from the original network, keyed by the hvdc line's own `.name`
        # (== its pypowsybl id). Side 1 pairs with `converter_station1_id` and
        # side 2 with `converter_station2_id`: `_from_pypowsybl.init()` builds
        # its per-side station data from those two columns in that same order.
        orig_hvdc = self._net.get_hvdc_lines()

        line_records = []
        station_records = []
        for el in self._grid.get_dclines():
            st1_id, st2_id = orig_hvdc.loc[el.name, ["converter_station1_id", "converter_station2_id"]]
            line_records.append({
                "id": el.name,
                "p1": -el.station1.res_p_mw, "q1": -el.station1.res_q_mvar,
                "p2": -el.station2.res_p_mw, "q2": -el.station2.res_q_mvar,
                "connected1": el.connected1, "connected2": el.connected2,
                "converter_station1_id": st1_id, "converter_station2_id": st2_id,
            })
            for station_id, station in ((st1_id, el.station1), (st2_id, el.station2)):
                station_records.append({
                    "id": station_id,
                    "p": -station.res_p_mw, "q": -station.res_q_mvar,
                    "bus_id": self._ls_bus_to_pypo(station.bus_id), "connected": station.connected,
                    "voltage_level_id": self._ls_sub_to_vl(station.voltage_level_id),
                    # 0 = VSC, 1 = LCC (ConverterStationInfo::ConverterType)
                    "is_lcc": station.converter_type == 1,
                })

        hvdc_columns = ["id", "p1", "q1", "p2", "q2", "connected1", "connected2",
                        "converter_station1_id", "converter_station2_id"]
        station_columns = ["id", "p", "q", "bus_id", "connected", "voltage_level_id", "is_lcc"]
        self._cache["hvdc_lines"] = self._records_to_frame(line_records, hvdc_columns)
        stations = self._records_to_frame(station_records, station_columns)
        # `.astype(bool)`: an empty `station_records` (no hvdc line in this grid)
        # leaves "is_lcc" as an empty object-dtype column, and a non-bool-dtype
        # empty mask silently drops every column (not just rows) when used to
        # index an empty DataFrame.
        is_lcc = stations["is_lcc"].astype(bool)
        self._cache["vsc_stations"] = stations[~is_lcc].drop(columns="is_lcc")
        self._cache["lcc_stations"] = stations[is_lcc].drop(columns="is_lcc")
