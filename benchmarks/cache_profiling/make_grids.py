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
* a ``_fancy`` variant of the bigger ones, where a few voltage-mode SVCs are
  added and a few *pairs* of generators are re-pointed at a common neighbouring
  load bus (a control GROUP: the bordered VoltageControl formulation, one voltage
  row and one reactive-sharing row per group). That is the configuration the
  voltage-control plan (``VoltageControlPlan``) is actually built for, and the
  only one where the size of its controller list is not zero. Every setpoint is
  a magnitude the grid already holds, so that each added controller is a no-op
  at the solution -- see `make_fancy` for the order this imposes.

With ``--two-busbars``, a third family is written instead: grids with two busbars per
substation (``_2bb``), the second one empty, for the phases that move elements
between busbars (``profile_batch``'s ``ss_topo_ac``). A synthetic grid (the 200-bus
Illinois case of the ACTIVSg series) and a larger one.

    python make_grids.py [output_dir] [--two-busbars]
"""

import os
import sys

import numpy as np
import pandapower as pp
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

# (file stem, pandapower factory): see the module docstring
TWO_BUSBAR_GRIDS = [
    ("illinois200_2bb", pn.case_illinois200),
    ("case1354pegase_2bb", pn.case1354pegase),
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


def _solved_vm(grid):
    """the solved magnitudes of `grid`, per bus (raises if it does not converge)"""
    v0 = np.full(grid.total_bus(), grid.get_init_vm_pu(), dtype=complex)
    vm = np.abs(grid.ac_pf(v0, MAX_ITER, TOL))
    if vm.size == 0:
        raise RuntimeError("the grid does not converge")
    return vm


def _build(net, groups, svc_buses, vm_groups, vm_svc):
    """
    One candidate grid: `groups` = [(gen_a, gen_b, regulated_bus)], `svc_buses` a
    list. `vm_svc` / `vm_groups` are the per-bus magnitudes the SVCs / the groups
    take their setpoints from (see make_fancy for which solve each comes from).
    """
    grid = init_from_pandapower(net)
    for gen_a, gen_b, reg_bus in groups:
        for gen_id in (gen_a, gen_b):
            grid.set_gen_regulated_bus(gen_id, reg_bus)
            # both members of a group must agree on the setpoint to the last bit,
            # or fill_voltage_control_solver_data rejects the configuration
            grid.change_v_gen(gen_id, float(vm_groups[reg_bus]))
    nb_svc = len(svc_buses)
    if nb_svc:
        buses = np.array(svc_buses, dtype=np.int32)
        grid.init_svcs([1] * nb_svc,                       # RegulationMode::VOLTAGE
                       np.array([vm_svc[b] for b in svc_buses]),
                       np.zeros(nb_svc),                   # q setpoint (unused in voltage mode)
                       np.zeros(nb_svc),                   # slope: none
                       np.full(nb_svc, -50.), np.full(nb_svc, 50.),
                       buses,                              # regulated bus: its own
                       buses)
    return grid


def make_fancy(net, nb_groups, nb_svc):
    """
    A copy of `net` with up to `nb_svc` voltage-mode SVCs and up to `nb_groups`
    two-generator control groups, that still converges.

    Every setpoint is a magnitude the grid already holds, so that each added
    controller is a no-op at the solution (the grid is asked for what it was doing
    anyway; without that the pegase cases diverge outright). That fixes the order
    things are built in:

    1. the SVCs first, with the plain case's solved magnitudes at their buses;
    2. then an AC solve of THAT grid -- every exotic element in place except the
       remote control -- gives the magnitudes the generator groups are pointed at:
       each group holds its target bus exactly where the grid with the SVCs
       already puts it.

    Even so, not every pair can be moved onto a neighbouring bus and keep the grid
    solvable, so each candidate is accepted only if the grid still converges with
    it -- which is why the counts above are upper bounds and the function reports
    what it actually built. (Same chunk-and-verify workaround as
    ``benchmarks/make_exotic_grid.cpp``, for the same reason: see the
    remote-voltage-control entry in the changelog's TODO section.)
    """
    base = init_from_pandapower(net)
    nb_bus = base.total_bus()
    vm_base = _solved_vm(base)

    gens = base.get_generators()
    gen_buses = {gens[i].bus_id for i in range(len(gens))}
    load_buses = {load.bus_id for load in base.get_loads()}
    adj = _adjacency(net, nb_bus)

    # 1. the SVCs, at load buses without a generator, held at the plain case's
    #    magnitude (so each injects nothing at the plain solution)
    svc_buses = []
    for bus in sorted(load_buses - gen_buses):
        if len(svc_buses) >= nb_svc:
            break
        if _solves(_build(net, [], svc_buses + [bus], vm_base, vm_base)):
            svc_buses.append(bus)

    # 2. the reference for the remote control: the grid with every exotic element
    #    in place except the remote control itself
    vm_ref = _solved_vm(_build(net, [], svc_buses, vm_base, vm_base))

    # 3. the groups: a group needs two ACTIVE local regulators to enrol (the slack
    #    is left alone), pointed at a load bus next to one of the two that no SVC
    #    holds and no other group regulates, and held exactly where the grid of
    #    step 2 already puts it
    candidates = [i for i in range(len(gens))
                  if gens[i].connected and not gens[i].is_slack and gens[i].voltage_regulator_on]
    groups, taken, i = [], set(svc_buses), 0
    while len(groups) < nb_groups and i + 1 < len(candidates):
        gen_a, gen_b = candidates[i], candidates[i + 1]
        i += 2
        common = [b for b in sorted(adj[gens[gen_a].bus_id] | adj[gens[gen_b].bus_id])
                  if b in load_buses and b not in gen_buses and b not in taken]
        if not common:
            continue
        reg_bus = common[0]
        if _solves(_build(net, groups + [(gen_a, gen_b, reg_bus)], svc_buses, vm_ref, vm_base)):
            groups.append((gen_a, gen_b, reg_bus))
            taken.add(reg_bus)

    grid = _build(net, groups, svc_buses, vm_ref, vm_base)
    if not _solves(grid):
        raise RuntimeError("the fancy grid does not converge")
    return grid, len(groups), len(svc_buses)


def two_busbars(net):
    """`net` with two busbars per substation: one substation per pandapower bus, its
    second busbar a copy of it, out of service, that no element stands on (the layout
    grid2op uses, which the converter expects)."""
    nb_bus = net.bus.shape[0]
    if list(net.bus.index) != list(range(nb_bus)):
        raise RuntimeError("the buses of this case are not numbered 0..n-1")
    pp.create_buses(net, nb_bus, vn_kv=net.bus.vn_kv.values, in_service=False)
    grid = init_from_pandapower(net, n_sub=nb_bus, n_busbar_per_sub=2)
    if not _solves(grid):
        raise RuntimeError("the two-busbar grid does not converge")
    return grid


def main(out_dir, with_two_busbars=False):
    os.makedirs(out_dir, exist_ok=True)
    if with_two_busbars:
        for name, factory in TWO_BUSBAR_GRIDS:
            net = factory()
            grid = two_busbars(net)
            path = os.path.join(out_dir, f"{name}.lsb")
            grid.save_binary(path)
            print(f"{name}: {grid.get_n_sub()} substations, two busbars each -> {path}")
        return
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
    args = [a for a in sys.argv[1:] if a != "--two-busbars"]
    main(args[0] if args else "grids", with_two_busbars="--two-busbars" in sys.argv[1:])
