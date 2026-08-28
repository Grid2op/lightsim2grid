# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

import warnings

import numpy as np

from ._my_const import (
    BUS_I, BUS_TYPE, PD, QD, GS, BS, BASE_KV, VMAX, VMIN,
    GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN,
    F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, TAP, SHIFT, BR_STATUS,
    DC_F_BUS, DC_T_BUS, DC_BR_STATUS, DC_PF, DC_VF, DC_VT,
    DC_QMINF, DC_QMAXF, DC_QMINT, DC_QMAXT, DC_LOSS0, DC_LOSS1,
)


def mpc_to_powermodels(bus, gen, branch, dcline, baseMVA) -> dict:
    """
    Convert raw MATPOWER-shaped `bus`/`gen`/`branch`/`dcline` matrices into a
    PowerModels.jl-style network data dictionary (see
    https://lanl-ansi.github.io/PowerModels.jl/stable/network-data/), consumable by
    `lightsim2grid.network.from_powermodels.init`. This is the (verified) inverse of
    the conversion PowerModels.jl itself performs when parsing a matpower `.m` file
    (`PowerModels.jl/src/io/matpower.jl`): loads/shunts are split off the bus rows,
    `tap == 0` (matpower's "this is a plain line" sentinel) becomes `tap = 1.0` with
    `transformer = False`, and angles are converted from degrees to radians.
    """
    network = {"baseMVA": float(baseMVA), "bus": {}, "gen": {}, "branch": {}, "load": {}, "shunt": {}}

    for i in range(bus.shape[0]):
        bus_i = int(bus[i, BUS_I])
        key = str(i + 1)
        network["bus"][key] = {
            "bus_i": bus_i,
            "bus_type": int(bus[i, BUS_TYPE]),
            "base_kv": float(bus[i, BASE_KV]),
            "vmax": float(bus[i, VMAX]),
            "vmin": float(bus[i, VMIN]),
        }
        pd, qd = float(bus[i, PD]), float(bus[i, QD])
        if pd != 0. or qd != 0.:
            network["load"][key] = {"load_bus": bus_i, "pd": pd, "qd": qd, "status": 1}
        gs, bs = float(bus[i, GS]), float(bus[i, BS])
        if gs != 0. or bs != 0.:
            # matpower's GS/BS are MW/MVAr "at V=1pu"; PowerModels' gs/bs are per-unit
            network["shunt"][key] = {"shunt_bus": bus_i, "gs": gs / baseMVA, "bs": bs / baseMVA, "status": 1}

    for i in range(gen.shape[0]):
        key = str(i + 1)
        network["gen"][key] = {
            "gen_bus": int(gen[i, GEN_BUS]),
            "pg": float(gen[i, PG]),
            "qg": float(gen[i, QG]),
            "qmax": float(gen[i, QMAX]),
            "qmin": float(gen[i, QMIN]),
            "vg": float(gen[i, VG]),
            "mbase": float(gen[i, MBASE]),
            "gen_status": int(gen[i, GEN_STATUS]),
            "pmax": float(gen[i, PMAX]),
            "pmin": float(gen[i, PMIN]),
        }

    weird_shift_rows = []
    for i in range(branch.shape[0]):
        key = str(i + 1)
        tap = float(branch[i, TAP])
        transformer = tap != 0.
        shift_degree = float(branch[i, SHIFT])
        if not transformer and shift_degree != 0.:
            # matpower's own reference makeYbus.m applies `tap = tap .* exp(1j*pi/180*
            # branch(:, SHIFT))` unconditionally, to every branch, regardless of TAP --
            # a TAP == 0 ("plain line") branch with a non-zero SHIFT still gets a
            # magnitude-1 phase-shifting tap in matpower's own solver (verified against
            # a real MATPOWER/Octave run: Ybus is genuinely asymmetric here, not the
            # symmetric plain-line value dropping SHIFT would give). So the shift must
            # be kept; `classify_branches`/`is_transformer` in `from_powermodels`'s
            # `_aux_add_branch.py` already routes any branch with a non-zero shift
            # through `init_trafo` (ratio=1) regardless of this dict's own
            # `transformer` flag, so no further change is needed here besides not
            # zeroing the value.
            weird_shift_rows.append(i)
        network["branch"][key] = {
            "f_bus": int(branch[i, F_BUS]),
            "t_bus": int(branch[i, T_BUS]),
            "br_r": float(branch[i, BR_R]),
            "br_x": float(branch[i, BR_X]),
            "g_fr": 0., "g_to": 0.,  # matpower has no line-conductance concept
            "b_fr": float(branch[i, BR_B]) / 2., "b_to": float(branch[i, BR_B]) / 2.,
            "tap": tap if transformer else 1.0,
            "shift": np.deg2rad(shift_degree),
            "transformer": transformer,
            "rate_a": float(branch[i, RATE_A]),
            "br_status": int(branch[i, BR_STATUS]),
        }
    if weird_shift_rows:
        warnings.warn(f"{len(weird_shift_rows)} branch(es) (0-based mpc `branch` row ids: "
                      f"{weird_shift_rows}) have `TAP == 0` (so are treated as a plain "
                      "powerline for their ratio) but a non-zero `SHIFT`. Matching "
                      "matpower's own makeYbus.m, the SHIFT is kept and these branches "
                      "are modeled as a zero-ratio (tap=1) phase-shifting transformer.")

    if dcline.shape[0]:
        network["dcline"] = {}
        for i in range(dcline.shape[0]):
            key = str(i + 1)
            network["dcline"][key] = {
                "f_bus": int(dcline[i, DC_F_BUS]),
                "t_bus": int(dcline[i, DC_T_BUS]),
                "br_status": int(dcline[i, DC_BR_STATUS]),
                "pf": float(dcline[i, DC_PF]),  # matpower keeps "pf"'s own sign, no flip needed
                "vf": float(dcline[i, DC_VF]),
                "vt": float(dcline[i, DC_VT]),
                "qminf": float(dcline[i, DC_QMINF]),
                "qmaxf": float(dcline[i, DC_QMAXF]),
                "qmint": float(dcline[i, DC_QMINT]),
                "qmaxt": float(dcline[i, DC_QMAXT]),
                "loss0": float(dcline[i, DC_LOSS0]),
                "loss1": float(dcline[i, DC_LOSS1]),
            }

    return network
