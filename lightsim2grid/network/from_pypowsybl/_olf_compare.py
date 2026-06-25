# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Thin helper to validate lightsim2grid against PowSyBl Open Load Flow (OLF) on
identical inputs.

The workflow:

1. Solve the network with OLF *with* outer loops.
2. :func:`bake_outer_loops` to freeze the converged outer-loop state into the
   network (see ``_olf_bake``).
3. Optionally apply the same topology change (e.g. a line outage) to both
   engines.
4. Solve the baked network loop-free in OLF (see
   :func:`get_pypowsybl_loopfree_parameters`) and in lightsim2grid.
5. Compare bus voltage magnitude (pu) and angle (deg).

Bus mapping
-----------
lightsim2grid uses its own internal solver-bus indexing, not the IIDM bus
order. The mapping is taken directly from ``grid._ls_to_orig``: solver bus
``i`` corresponds to pypowsybl original bus index ``_ls_to_orig[i]``, i.e.
``network.get_buses().index[_ls_to_orig[i]]``. Building the grid with
``sort_index=False`` makes ``_ls_to_orig`` the identity, so it aligns 1:1 with
``get_buses()``; ``buses_for_sub=False`` keeps the bus set simple. This covers
*every* solver bus, including injection-free junction buses, which a mapping
inferred from element names cannot.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pypowsybl as pp
import pypowsybl.loadflow as lf

from ._from_pypowsybl import init as init_from_pypowsybl
from ._olf_bake import bake_outer_loops
from ._olf_params import get_pypowsybl_loopfree_parameters


def iidm_bus_voltages(network) -> pd.DataFrame:
    """Per-bus |V| (pu) and angle (deg) from a solved IIDM network, indexed by
    bus id (in ``get_buses()`` order), with the kV->pu conversion applied."""
    buses = network.get_buses()
    vls = network.get_voltage_levels()[["nominal_v"]]
    buses = buses.join(vls, on="voltage_level_id")
    out = pd.DataFrame(index=buses.index)
    out["vm_pu"] = buses["v_mag"] / buses["nominal_v"]
    out["va_deg"] = buses["v_angle"]
    return out


def lightsim_bus_to_iidm(grid, network) -> dict:
    """Map lightsim2grid solver bus index -> IIDM bus id.

    Uses lightsim2grid's own ``_ls_to_orig`` array, which gives the index of
    each solver bus in pypowsybl's *original* bus order. Paired with
    ``network.get_buses().index`` (same original order) this pins every solver
    bus -- including injection-free junction buses -- to a concrete IIDM bus
    id, with no element-name parsing and no guessing.

    Requires the grid to have been built with ``init_from_pypowsybl(...,
    sort_index=False, buses_for_sub=False)`` so that the original order matches
    ``get_buses()`` one-to-one.
    """
    ls_to_orig = np.asarray(grid._ls_to_orig)
    orig_bus_ids = list(network.get_buses().index)
    return {
        ls: orig_bus_ids[orig]
        for ls, orig in enumerate(ls_to_orig)
        if orig < len(orig_bus_ids)
    }


@dataclass
class ComparisonResult:
    """Outcome of :func:`compare_baked`: how far lightsim2grid and OLF disagree.

    Attributes
    ----------
    max_dvm_pu : float
        Largest absolute voltage-magnitude difference, in per unit, over every
        bus common to both engines.
    max_dva_deg : float
        Largest absolute voltage-angle difference, in degrees (raw).
    max_dva_deg_offset_removed : float
        Same as ``max_dva_deg`` but with a uniform angle offset removed first.
        A constant offset on all buses is just a difference of reference-datum
        convention between the two engines, not a physical disagreement, so this
        is usually the meaningful angle metric.
    table : pandas.DataFrame
        Per-bus detail, indexed by IIDM bus id, with the OLF and lightsim2grid
        magnitudes / angles and their differences (columns ``olf_vm``,
        ``ls_vm``, ``olf_va``, ``ls_va``, ``dvm``, ``dva``).
    """

    max_dvm_pu: float
    max_dva_deg: float
    max_dva_deg_offset_removed: float
    table: pd.DataFrame

    def __repr__(self):
        return (
            f"ComparisonResult(max |dV| = {self.max_dvm_pu:.3e} pu, "
            f"max |dAngle| = {self.max_dva_deg:.3e} deg, "
            f"offset-removed = {self.max_dva_deg_offset_removed:.3e} deg)"
        )


def compare_baked(
    network_factory,
    slack_gen_id: str,
    line_outages=None,
    trafo_outages=None,
    olf_loop_params: "lf.Parameters | None" = None,
):
    """Bake, optionally apply outages, solve in both engines, and compare.

    Parameters
    ----------
    network_factory : callable returning a fresh pypowsybl Network
        Called twice (once per engine) so the two starts are identical.
        lightsim2grid mutates/consumes the network it is built from, so a fresh
        instance is needed for each side.
    slack_gen_id : str
        Generator id to use as the lightsim2grid slack.
    line_outages, trafo_outages : list of str, optional
        IIDM ids to disconnect identically in both engines after baking.
    olf_loop_params : pypowsybl.loadflow.Parameters, optional
        Parameters for the initial *with-loops* OLF solve. Defaults to
        distributed slack + reactive limits.

    Returns
    -------
    ComparisonResult
    """
    line_outages = line_outages or []
    trafo_outages = trafo_outages or []
    if olf_loop_params is None:
        # ``twt_split_shunt_admittance=True`` matches the transformer model used
        # by the loop-free OLF solve and by lightsim2grid, so the baked
        # setpoints are derived under the same model as the comparison solve.
        olf_loop_params = lf.Parameters(
            distributed_slack=True,
            use_reactive_limits=True,
            balance_type=pp.loadflow.BalanceType.PROPORTIONAL_TO_GENERATION_P_MAX,
            twt_split_shunt_admittance=True,
        )
    loop_free_params = get_pypowsybl_loopfree_parameters()

    # ---- OLF side: with-loops solve, bake, outage, loop-free solve ----
    n_olf = network_factory()
    lf.run_ac(n_olf, olf_loop_params)
    bake_outer_loops(n_olf)
    for lid in line_outages:
        n_olf.update_lines(id=lid, connected1=False, connected2=False)
    for tid in trafo_outages:
        n_olf.update_2_windings_transformers(id=tid, connected1=False, connected2=False)
    res = lf.run_ac(n_olf, loop_free_params)
    if res[0].status != pp.loadflow.ComponentStatus.CONVERGED:
        raise RuntimeError(f"OLF loop-free did not converge: {res[0].status}")
    olf_v = iidm_bus_voltages(n_olf)  # indexed by IIDM bus id

    # ---- lightsim2grid side: bake, outage, solve ----
    n_ls = network_factory()
    lf.run_ac(n_ls, olf_loop_params)
    bake_outer_loops(n_ls)
    grid = init_from_pypowsybl(
        n_ls,
        gen_slack_id=slack_gen_id,
        sort_index=False,
        buses_for_sub=False,
    )
    solver_to_iidm = lightsim_bus_to_iidm(grid, n_ls)
    line_ids = list(n_ls.get_lines().index)
    for lid in line_outages:
        grid.deactivate_powerline(line_ids.index(lid))
    trafo_ids = list(n_ls.get_2_windings_transformers().index)
    for tid in trafo_outages:
        grid.deactivate_trafo(trafo_ids.index(tid))

    V = grid.ac_pf(np.full(grid.total_bus(), 1.06 + 0j), 20, 1e-10)
    if len(V) == 0:
        raise RuntimeError("lightsim2grid did not converge")

    n_bus = len(V)
    ls_v = pd.DataFrame(
        {
            "vm_pu": [abs(V[s]) for s in range(n_bus)],
            "va_deg": [np.degrees(np.angle(V[s])) for s in range(n_bus)],
        },
        index=[solver_to_iidm.get(s) for s in range(n_bus)],
    )
    ls_v = ls_v[ls_v.index.notna()]  # drop any unmapped junction bus

    # ---- compare on the intersection of mapped IIDM bus ids ----
    common = olf_v.index.intersection(ls_v.index)
    table = pd.DataFrame(
        {
            "olf_vm": olf_v.loc[common, "vm_pu"],
            "ls_vm": ls_v.loc[common, "vm_pu"],
            "olf_va": olf_v.loc[common, "va_deg"],
            "ls_va": ls_v.loc[common, "va_deg"],
        }
    )
    table["dvm"] = (table["olf_vm"] - table["ls_vm"]).abs()
    table["dva"] = (table["olf_va"] - table["ls_va"]).abs()
    # A uniform angle offset is just a reference-datum convention difference;
    # report both raw and offset-removed.
    offset = (table["olf_va"] - table["ls_va"]).median()
    dva_rel = ((table["olf_va"] - table["ls_va"]) - offset).abs()

    return ComparisonResult(
        max_dvm_pu=float(table["dvm"].max()),
        max_dva_deg=float(table["dva"].max()),
        max_dva_deg_offset_removed=float(dva_rel.max()),
        table=table,
    )
