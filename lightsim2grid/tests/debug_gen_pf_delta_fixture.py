# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
One-off (not collected by unittest / pytest) script used to (re)generate the committed
`pf_delta_case14.json` fixture used by `test_LSGrid_pf_delta.py`.

Real PFΔ data could not be downloaded for this repo (huggingface.co is blocked in the dev
sandbox this was written in), so this script builds a PFΔ-schema-shaped row from an
independent, non-lightsim2grid source of truth: pandapower's own internal per-unit "ppc"
representation (`net["_ppc"]`, which uses the exact same MATPOWER/PyPower column layout
`lightsim2grid.network.from_matpower` and PowerModels.jl both use) of its bundled IEEE
case14 network, solved with pandapower's own AC power flow. This sidesteps having to
hand-derive pandapower's own trafo T-model -> per-unit conversion.

Run manually (requires `pip install pandapower`):
    python lightsim2grid/tests/debug_gen_pf_delta_fixture.py
"""

import json
import os

import numpy as np
import pandapower as pp
import pandapower.networks as pn
from pandapower.pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BASE_KV, VMAX, VMIN, VM, VA
from pandapower.pypower.idx_gen import GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN
from pandapower.pypower.idx_brch import (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, TAP, SHIFT,
                                          BR_STATUS, PF, QF, PT, QT)


def build_pf_delta_row(net) -> dict:
    ppc = net["_ppc"]
    bus, gen, branch = ppc["bus"].real, ppc["gen"].real, ppc["branch"].real
    base_mva = float(ppc["baseMVA"])

    network = {"baseMVA": base_mva, "bus": {}, "gen": {}, "branch": {}, "load": {}, "shunt": {}}
    solution = {"bus": {}, "gen": {}, "branch": {}}

    load_idx, shunt_idx = 1, 1
    for i in range(bus.shape[0]):
        bus_i = int(bus[i, BUS_I]) + 1  # 1-based, mirrors PowerModels' own convention
        key = str(i + 1)
        network["bus"][key] = {
            "bus_i": bus_i,
            "bus_type": int(bus[i, BUS_TYPE]),
            "base_kv": float(bus[i, BASE_KV]),
            "vmax": float(bus[i, VMAX]),
            "vmin": float(bus[i, VMIN]),
        }
        solution["bus"][key] = {"vm": float(bus[i, VM]), "va": float(np.deg2rad(bus[i, VA]))}

        pd_, qd_ = float(bus[i, PD]), float(bus[i, QD])
        if pd_ != 0.0 or qd_ != 0.0:
            network["load"][str(load_idx)] = {"load_bus": bus_i, "pd": pd_, "qd": qd_, "status": 1}
            load_idx += 1

        gs_, bs_ = float(bus[i, GS]), float(bus[i, BS])
        if gs_ != 0.0 or bs_ != 0.0:
            # matpower's GS/BS are MW/MVAr "at V=1pu"; PowerModels' gs/bs are per-unit
            network["shunt"][str(shunt_idx)] = {"shunt_bus": bus_i, "gs": gs_ / base_mva,
                                                "bs": bs_ / base_mva, "status": 1}
            shunt_idx += 1

    for i in range(gen.shape[0]):
        key = str(i + 1)
        network["gen"][key] = {
            "gen_bus": int(gen[i, GEN_BUS]) + 1,
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
        solution["gen"][key] = {"pg": float(gen[i, PG]), "qg": float(gen[i, QG])}

    for i in range(branch.shape[0]):
        key = str(i + 1)
        tap = float(branch[i, TAP])
        shift_deg = float(branch[i, SHIFT])
        is_transformer = (abs(tap - 1.0) > 1e-9) or (shift_deg != 0.0)
        b_total = float(branch[i, BR_B])
        network["branch"][key] = {
            "f_bus": int(branch[i, F_BUS]) + 1,
            "t_bus": int(branch[i, T_BUS]) + 1,
            "br_r": float(branch[i, BR_R]),
            "br_x": float(branch[i, BR_X]),
            "g_fr": 0.0,
            "g_to": 0.0,
            "b_fr": b_total / 2.0,
            "b_to": b_total / 2.0,
            "tap": tap,
            "shift": np.deg2rad(shift_deg),  # PowerModels stores angles in radians
            "transformer": is_transformer,
            "rate_a": float(branch[i, RATE_A]),
            "br_status": int(branch[i, BR_STATUS]),
        }
        if int(branch[i, BR_STATUS]) != 0:
            solution["branch"][key] = {
                "pf": float(branch[i, PF]), "qf": float(branch[i, QF]),
                "pt": float(branch[i, PT]), "qt": float(branch[i, QT]),
            }

    return {
        "network": network,
        "solution": {"solution": solution, "termination_status": "SOLVED_PF", "objective": 0.0},
        "parent_seed": -1,
    }


if __name__ == "__main__":
    net = pn.case14()
    pp.runpp(net, numba=False)
    row = build_pf_delta_row(net)

    out_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pf_delta_case14.json")
    with open(out_path, "w") as f:
        json.dump(row, f, indent=2)

    n_trafo = sum(1 for br in row["network"]["branch"].values() if br["transformer"])
    n_shunt = len(row["network"]["shunt"])
    print(f"wrote {out_path}")
    print(f"  {len(row['network']['bus'])} buses, {len(row['network']['branch'])} branches "
          f"({n_trafo} transformers), {len(row['network']['gen'])} gens, "
          f"{len(row['network']['load'])} loads, {n_shunt} shunts")
    for br in row["network"]["branch"].values():
        if br["transformer"]:
            print(f"  sample transformer: tap={br['tap']}, shift={br['shift']} rad")
    for sh in row["network"]["shunt"].values():
        print(f"  sample shunt: gs={sh['gs']} pu, bs={sh['bs']} pu")
