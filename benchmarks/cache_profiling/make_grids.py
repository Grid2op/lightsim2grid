#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Dump a few pandapower reference grids as lightsim2grid binary (``.lsb``) files, so
that the C++ profiling driver (``profile_cached_pf.cpp``) can load them without
needing python, pandapower or the conversion code in the profile.

Two families are written:

* the plain pandapower cases, which have no remote voltage control and no SVC:
  every generator regulates its own bus, the ordinary PV path;
* a ``_fancy`` variant of the bigger ones, where a few *pairs* of generators are
  re-pointed at a common neighbouring load bus (a control GROUP: the bordered
  VoltageControl formulation, one voltage row and one reactive-sharing row per
  group) and a few voltage-mode SVCs are added. That is the configuration the
  voltage-control plan (``VoltageControlPlan``) is actually built for, and the
  only one where the size of its controller list is not zero.

    python make_grids.py [output_dir]
"""

import os
import sys

import numpy as np
import pandapower.networks as pn

from lightsim2grid.network import init_from_pandapower


# (file stem, pandapower factory).  Small + large, as asked: one IEEE case that
# fits in cache and one large pegase case that does not.
GRIDS = [
    ("case30", pn.case30),
    ("case118", pn.case118),
    ("case1354pegase", pn.case1354pegase),
    ("case9241pegase", pn.case9241pegase),
]

# (file stem, pandapower factory, nb of control groups, nb of SVCs) -- see
# `make_fancy` for what is actually built, and why the counts are upper bounds.
FANCY_GRIDS = [
    ("case118_fancy", pn.case118, 8, 4),
    ("case1354pegase_fancy", pn.case1354pegase, 20, 10),
    ("case9241pegase_fancy", pn.case9241pegase, 40, 30),
]

MAX_ITER = 10
TOL = 1e-8


def _adjacency(net, nb_bus):
    """bus -> set of buses one branch away, from the pandapower tables."""
    adj = {b: set() for b in range(nb_bus)}
    for f, t in zip(net.line.from_bus.values, net.line.to_bus.values):
        adj[int(f)].add(int(t))
        adj[int(t)].add(int(f))
    for f, t in zip(net.trafo.hv_bus.values, net.trafo.lv_bus.values):
        adj[int(f)].add(int(t))
        adj[int(t)].add(int(f))
    return adj


def _solves(grid):
    # from the SAME flat start the profiling driver uses (`grid.get_init_vm_pu()`,
    # not 1.0): a candidate accepted from a different initial point is a candidate
    # that can still diverge under the driver.
    v0 = np.full(grid.total_bus(), grid.get_init_vm_pu(), dtype=complex)
    return grid.ac_pf(v0, MAX_ITER, TOL).shape[0] > 0


def _build(net, groups, svc_buses, vm):
    """One candidate grid: `groups` = [(gen_a, gen_b, regulated_bus)], `svc_buses` a list."""
    grid = init_from_pandapower(net)
    for gen_a, gen_b, reg_bus in groups:
        for gen_id in (gen_a, gen_b):
            grid.set_gen_regulated_bus(gen_id, reg_bus)
            # both members of a group must agree on the setpoint to the last bit,
            # or fill_voltage_control_solver_data rejects the configuration
            grid.change_v_gen(gen_id, float(vm[reg_bus]))
    nb_svc = len(svc_buses)
    if nb_svc:
        buses = np.array(svc_buses, dtype=np.int32)
        grid.init_svcs([1] * nb_svc,                       # RegulationMode::VOLTAGE
                       np.array([vm[b] for b in svc_buses]),
                       np.zeros(nb_svc),                   # q setpoint (unused in voltage mode)
                       np.zeros(nb_svc),                   # slope: none
                       np.full(nb_svc, -50.), np.full(nb_svc, 50.),
                       buses,                              # regulated bus: its own
                       buses)
    return grid


def make_fancy(net, nb_groups, nb_svc):
    """
    A copy of `net` with up to `nb_groups` two-generator control groups and up to
    `nb_svc` voltage-mode SVCs, that still converges.

    Setpoints are the BASE case's own solved magnitudes, so the added controllers
    ask the grid for what it was doing anyway; without that the pegase cases
    diverge outright. Even so, not every pair can be moved onto a neighbouring
    bus and keep the grid solvable, so each candidate is accepted only if the
    grid still converges with it -- which is why the counts above are upper
    bounds and the function reports what it actually built. (Same chunk-and-verify
    workaround as ``benchmarks/make_exotic_grid.cpp``, for the same reason: see the
    remote-voltage-control entry in the changelog's TODO section.)
    """
    base = init_from_pandapower(net)
    nb_bus = base.total_bus()
    v0 = np.full(nb_bus, base.get_init_vm_pu(), dtype=complex)
    vm = np.abs(base.ac_pf(v0, MAX_ITER, TOL))
    if vm.size == 0:
        raise RuntimeError("the base case does not converge")

    gens = base.get_generators()
    gen_buses = {gens[i].bus_id for i in range(len(gens))}
    load_buses = {load.bus_id for load in base.get_loads()}
    adj = _adjacency(net, nb_bus)
    # a group needs two ACTIVE local regulators to enrol; the slack is left alone
    candidates = [i for i in range(len(gens))
                  if gens[i].connected and not gens[i].is_slack and gens[i].voltage_regulator_on]

    groups, taken, i = [], set(), 0
    while len(groups) < nb_groups and i + 1 < len(candidates):
        gen_a, gen_b = candidates[i], candidates[i + 1]
        i += 2
        # a load bus next to one of the two, not already regulated by a group
        common = [b for b in sorted(adj[gens[gen_a].bus_id] | adj[gens[gen_b].bus_id])
                  if b in load_buses and b not in gen_buses and b not in taken]
        if not common:
            continue
        reg_bus = common[0]
        if _solves(_build(net, groups + [(gen_a, gen_b, reg_bus)], [], vm)):
            groups.append((gen_a, gen_b, reg_bus))
            taken.add(reg_bus)

    svc_buses = []
    for bus in sorted(load_buses - gen_buses):
        if len(svc_buses) >= nb_svc:
            break
        if bus in taken:
            continue
        if _solves(_build(net, groups, svc_buses + [bus], vm)):
            svc_buses.append(bus)

    grid = _build(net, groups, svc_buses, vm)
    if not _solves(grid):
        raise RuntimeError("the fancy grid does not converge")
    return grid, len(groups), len(svc_buses)


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for name, factory in GRIDS:
        net = factory()
        grid = init_from_pandapower(net)
        path = os.path.join(out_dir, f"{name}.lsb")
        grid.save_binary(path)
        print(f"{name}: {grid.total_bus()} buses, "
              f"{len(net.line)} lines, {len(net.trafo)} trafos, "
              f"{len(net.load)} loads -> {path}")
    for name, factory, nb_groups, nb_svc in FANCY_GRIDS:
        net = factory()
        grid, got_groups, got_svc = make_fancy(net, nb_groups, nb_svc)
        path = os.path.join(out_dir, f"{name}.lsb")
        grid.save_binary(path)
        print(f"{name}: {grid.total_bus()} buses, {got_groups} control groups "
              f"({2 * got_groups} generators), {got_svc} voltage-mode SVCs -> {path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "grids")
