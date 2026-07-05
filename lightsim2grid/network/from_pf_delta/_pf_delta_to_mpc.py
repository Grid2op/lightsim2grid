# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import numpy as np

from ..from_matpower._my_const import (
    BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN,
    GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,
    F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS,
    DC_F_BUS, DC_T_BUS, DC_BR_STATUS, DC_PF, DC_VF, DC_VT,
    DC_QMINF, DC_QMAXF, DC_QMINT, DC_QMAXT, DC_LOSS0, DC_LOSS1,
    MIN_BUS_COLS, MIN_GEN_COLS, MIN_BRANCH_COLS, MIN_DCLINE_COLS,
)


def _is_transformer(branch):
    """A PowerModels branch is a transformer if either its own `transformer` flag says
    so, or its `tap`/`shift` are away from the "plain line" neutral (`tap == 1.0`,
    `shift == 0.0`) -- PowerModels normalizes a matpower `TAP == 0` (its own "this is a
    plain line" sentinel) to `tap = 1.0`, so this loader cannot just test `tap != 0`
    like `from_matpower.get_branch_split` does on a raw mpc branch matrix.
    """
    return bool(branch.get("transformer", False)) or branch.get("tap", 1.0) != 1.0 or branch.get("shift", 0.0) != 0.0


def network_to_mpc(network: dict) -> dict:
    """
    Convert a PFΔ / PowerModels.jl ``network`` dict into a raw-MATPOWER-shaped dict
    (``{"bus", "gen", "branch"}`` numpy arrays in matpower's own column layout, plus
    ``"baseMVA"``) consumable by :func:`lightsim2grid.network.from_matpower.init`.

    Row order is ``sorted(network[<elem>], key=int)`` for bus / gen / branch: this is a
    deterministic, documented contract (not just an implementation detail), relied upon
    by tests that need to map a built `LSGrid`'s element positions back to the original
    PFΔ row's bus / gen / branch keys.
    """
    if network.get("storage"):
        raise RuntimeError('"storage" is present in this PFΔ row but is not supported '
                           "by this loader (from_matpower has no storage container either).")
    if network.get("multinetwork"):
        raise RuntimeError("multinetwork PowerModels dicts are not supported by this loader.")

    base_mva = float(network["baseMVA"])

    bus_keys = sorted(network["bus"], key=int)
    gen_keys = sorted(network["gen"], key=int)
    branch_keys = sorted(network["branch"], key=int)

    # raw matpower embeds PD / QD / GS / BS directly on the bus row; PowerModels splits
    # them into separate "load" / "shunt" dicts, so they need to be aggregated back
    # (only active elements contribute -- matpower has no separate "load/shunt present
    # but off" column, so an inactive one is equivalent to 0 demand/admittance here)
    pd = {int(bus["bus_i"]): 0.0 for bus in network["bus"].values()}
    qd, gs, bs = dict(pd), dict(pd), dict(pd)
    for load in network.get("load", {}).values():
        if load.get("status", 1) == 0:
            continue
        b = int(load["load_bus"])
        pd[b] += load["pd"]
        qd[b] += load["qd"]
    for shunt in network.get("shunt", {}).values():
        if shunt.get("status", 1) == 0:
            continue
        b = int(shunt["shunt_bus"])
        gs[b] += shunt["gs"] * base_mva
        bs[b] += shunt["bs"] * base_mva

    bus = np.zeros((len(bus_keys), MIN_BUS_COLS))
    for i, k in enumerate(bus_keys):
        b = network["bus"][k]
        bus_i = int(b["bus_i"])
        bus[i, BUS_I] = bus_i
        bus[i, BUS_TYPE] = b["bus_type"]
        bus[i, PD] = pd[bus_i]
        bus[i, QD] = qd[bus_i]
        bus[i, GS] = gs[bus_i]
        bus[i, BS] = bs[bus_i]
        bus[i, BUS_AREA] = 1.0
        bus[i, VM] = 1.0
        bus[i, VA] = 0.0
        bus[i, BASE_KV] = b.get("base_kv", 1.0)
        bus[i, ZONE] = 1.0
        bus[i, VMAX] = b["vmax"]
        bus[i, VMIN] = b["vmin"]

    gen = np.zeros((len(gen_keys), MIN_GEN_COLS))
    for i, k in enumerate(gen_keys):
        g = network["gen"][k]
        gen[i, GEN_BUS] = int(g["gen_bus"])
        gen[i, PG] = g["pg"]
        gen[i, QG] = g["qg"]
        gen[i, QMAX] = g["qmax"]
        gen[i, QMIN] = g["qmin"]
        gen[i, VG] = g["vg"]
        gen[i, MBASE] = g.get("mbase", base_mva)
        gen[i, GEN_STATUS] = g.get("gen_status", 1)
        gen[i, PMAX] = g["pmax"]
        gen[i, PMIN] = g["pmin"]

    branch = np.zeros((len(branch_keys), MIN_BRANCH_COLS))
    for i, k in enumerate(branch_keys):
        br = network["branch"][k]
        if br.get("g_fr", 0.0) != 0.0 or br.get("g_to", 0.0) != 0.0:
            raise RuntimeError(f'branch "{k}" has a non-zero "g_fr"/"g_to" (line '
                               f"conductance); the raw MATPOWER format this loader "
                               f"delegates to cannot represent that (only susceptance).")
        is_trafo = _is_transformer(br)
        branch[i, F_BUS] = int(br["f_bus"])
        branch[i, T_BUS] = int(br["t_bus"])
        branch[i, BR_R] = br["br_r"]
        branch[i, BR_X] = br["br_x"]
        branch[i, BR_B] = br.get("b_fr", 0.0) + br.get("b_to", 0.0)
        branch[i, RATE_A] = br.get("rate_a", 0.0)
        branch[i, RATE_B] = br.get("rate_b", 0.0)
        branch[i, RATE_C] = br.get("rate_c", 0.0)
        # matpower's own "this is a plain line" sentinel is TAP == 0, which is what
        # `from_matpower.get_branch_split` tests for; PowerModels never uses 0 (it
        # normalizes a matpower TAP==0 to tap=1.0), so it has to be re-encoded here
        branch[i, TAP] = br.get("tap", 1.0) if is_trafo else 0.0
        # PowerModels stores shift/angles in radians (converted from matpower's
        # degrees at parse time); matpower's raw SHIFT column is in degrees
        branch[i, SHIFT] = np.rad2deg(br.get("shift", 0.0))
        branch[i, BR_STATUS] = br.get("br_status", 1)

    mpc = {"bus": bus, "gen": gen, "branch": branch, "baseMVA": base_mva}

    if network.get("dcline"):
        dcline_keys = sorted(network["dcline"], key=int)
        dcline = np.zeros((len(dcline_keys), MIN_DCLINE_COLS))
        for i, k in enumerate(dcline_keys):
            dc = network["dcline"][k]
            dcline[i, DC_F_BUS] = int(dc["f_bus"])
            dcline[i, DC_T_BUS] = int(dc["t_bus"])
            dcline[i, DC_BR_STATUS] = dc.get("br_status", 1)
            # PowerModels keeps "pf" at matpower's own value/sign (only "pt"/"qf"/"qt"
            # get sign-flipped by PowerModels.jl's parser); from_matpower's own dcline
            # converter is the one that negates it to match lightsim2grid's "received
            # at side 1" convention, so it must be passed through here unchanged.
            dcline[i, DC_PF] = dc["pf"]
            dcline[i, DC_VF] = dc.get("vf", 1.0)
            dcline[i, DC_VT] = dc.get("vt", 1.0)
            dcline[i, DC_QMINF] = dc["qminf"]
            dcline[i, DC_QMAXF] = dc["qmaxf"]
            dcline[i, DC_QMINT] = dc["qmint"]
            dcline[i, DC_QMAXT] = dc["qmaxt"]
            dcline[i, DC_LOSS0] = dc["loss0"]
            dcline[i, DC_LOSS1] = dc["loss1"]
            # DC_PT / DC_QF / DC_QT / DC_PMIN / DC_PMAX are not read by
            # from_matpower's dcline converter (see _aux_add_dc_line.py) -- left at 0.0
        mpc["dcline"] = dcline

    return mpc
