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
        v_cplx = np.asarray(grid.get_V())[orig_to_ls]
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
        if "lines" not in self._cache:
            self._cache["lines"] = self._build_two_sided(self._grid.get_lines())
        return self._maybe_select(self._cache["lines"], attributes)

    def get_2_windings_transformers(self, attributes: Optional[List[str]] = None) -> pd.DataFrame:
        if "trafos" not in self._cache:
            self._cache["trafos"] = self._build_two_sided(self._grid.get_trafos())
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
