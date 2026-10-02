#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Compare lightsim2grid with OpenLoadFlow on grid snapshots, outer loop by outer loop.

The development tool of ``docs/dev_notes/outer_loops_fixed_sparsity.md``: every outer loop
lightsim2grid implements is checked against OpenLoadFlow running **that loop only**, under
tight solves, before the next one is started. Both engines start from the same, unbaked
network; nothing is frozen by :func:`bake_outer_loops`.

Usage::

    python utils/olf_outer_compare.py SNAPSHOT_OR_DIR [...]                 # no loop at all
    python utils/olf_outer_compare.py DIR --loops DistributedSlack
    python utils/olf_outer_compare.py DIR --loops DistributedSlack --algo NROuter_KLU
    python utils/olf_outer_compare.py DIR --slack-mode most_meshed --csv out.csv

A directory is searched recursively for ``*.arc.xz`` / ``*.xiidm`` / ``*.iidm`` files.
Without a path, ``$LS2G_SNAPSHOTS`` is used. Snapshots are never part of the repository.

Selecting the loops on the OpenLoadFlow side
--------------------------------------------
OpenLoadFlow builds its loop list either from the parameters (each loop has a flag) or from
an explicit ``outerLoopNames`` list, run in the list's order. The explicit list cannot name
``AcHvdcAcEmulationLimits``: that loop only exists in the parameter-driven list, where
``hvdc_ac_emulation`` creates it. So:

* without ``AcHvdcAcEmulationLimits``, the explicit list is used, in OpenLoadFlow's default
  order, and the HVDC lines in AC emulation stay linear (no limit loop);
* with it, the parameter-driven list is used, every other loop being created by its flag.

``use_reactive_limits`` both creates the ``ReactiveLimits`` loop and turns on loading rules
(a unit with a too narrow reactive range does not regulate). It is set iff ``ReactiveLimits``
is requested, so that a run without it compares like with like.

The slack
---------
lightsim2grid's single slack is a generator: the reference unit of OpenLoadFlow's default
distributed slack (``_default_distributed_slack``). With ``--slack-mode name`` (default)
OpenLoadFlow is pinned on that generator's bus; with ``most_meshed`` it picks its own bus,
which only makes sense together with ``DistributedSlack`` (otherwise each engine's slack bus
absorbs the imbalance at a different place). Angles are always compared after rotating both
solutions so that they agree on lightsim2grid's slack bus.
"""

import argparse
import csv
import glob
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import pypowsybl as pp
import pypowsybl.loadflow as lf

from lightsim2grid.network.from_pypowsybl import init as init_from_pypowsybl
from lightsim2grid.network.from_pypowsybl import LightsimResultNetwork, OlfLoadingParameters
from lightsim2grid.network.from_pypowsybl._olf_compare import iidm_bus_voltages
from lightsim2grid.network.from_pypowsybl._aux_add_slack import _default_distributed_slack
from lightsim2grid.lightsim2grid_cpp import ScalingPolicyType
from lightsim2grid import algorithm as _algorithm

LOAD_PARAMETERS = {"iidm.die.with-extensions": "all"}

#: below this, a reactive power difference with OLF means nothing (see ``compare``)
Q_PRECISION_MVAR = 1e-3
SNAPSHOT_PATTERNS = ("*.arc.xz", "*.xiidm", "*.iidm")

#: OpenLoadFlow's default order of the loops in scope (``DefaultAcOuterLoopConfig``, v2.3.0)
OLF_ORDER = (
    "DistributedSlack",
    "AcHvdcAcEmulationLimits",
    "VoltageMonitoring",
    "ReactiveLimits",
    "PhaseControl",
    "TransformerVoltageControl",
    "ShuntVoltageControl",
)


def find_snapshots(paths):
    """Expand files and directories (searched recursively) into a sorted list of snapshots."""
    found = []
    for path in paths:
        if os.path.isdir(path):
            for pattern in SNAPSHOT_PATTERNS:
                found.extend(glob.glob(os.path.join(path, "**", pattern), recursive=True))
        else:
            found.append(path)
    return sorted(set(found))


def load(path, extra_load_mw=0.):
    """The snapshot, with ``extra_load_mw`` added to its largest load of the main synchronous
    component (a mismatch for the slack loops to share; the same on both engines)."""
    net = pp.network.load(path, LOAD_PARAMETERS)
    if extra_load_mw:
        loads = net.get_loads(attributes=["p0", "bus_id"])
        sync = net.get_buses(attributes=["synchronous_component"])["synchronous_component"]
        main = sync.value_counts().idxmax()
        in_main = loads["bus_id"].map(sync) == main
        big = loads.loc[in_main, "p0"].idxmax()
        net.update_loads(id=big, p0=float(loads.loc[big, "p0"]) + float(extra_load_mw))
    return net


def pick_slack(net):
    """The generator lightsim2grid uses as its single slack, and its bus id."""
    slack = _default_distributed_slack(net, net.get_generators())
    if not slack:
        raise RuntimeError("no generator takes part in the default distributed slack")
    gen_id = next(iter(slack))
    return gen_id, net.get_generators().loc[gen_id, "bus_id"]


def olf_parameters(loops, slack_bus_id=None, conv_eps=1e-9, max_nr_iter=50):
    """OpenLoadFlow parameters running exactly ``loops``, everything else at the installed
    pypowsybl's defaults (see the module docstring for how the loops are selected)."""
    unknown = set(loops) - set(OLF_ORDER)
    if unknown:
        raise ValueError(f"unknown loop(s) {sorted(unknown)}, use some of {OLF_ORDER}")
    loops = [name for name in OLF_ORDER if name in loops]

    params = lf.Parameters()
    params.distributed_slack = "DistributedSlack" in loops
    params.use_reactive_limits = "ReactiveLimits" in loops
    params.phase_shifter_regulation_on = "PhaseControl" in loops
    params.transformer_voltage_control_on = "TransformerVoltageControl" in loops
    params.shunt_compensator_voltage_control_on = "ShuntVoltageControl" in loops
    provider = dict(params.provider_parameters)
    provider["svcVoltageMonitoring"] = str("VoltageMonitoring" in loops).lower()
    provider["newtonRaphsonConvEpsPerEq"] = repr(conv_eps)
    provider["maxNewtonRaphsonIterations"] = str(max_nr_iter)
    if "AcHvdcAcEmulationLimits" in loops:
        provider.pop("outerLoopNames", None)
    else:
        params.hvdc_ac_emulation = True  # the lines keep their model, only the loop is dropped
        provider["outerLoopNames"] = ",".join(loops)
    if slack_bus_id is not None:
        params.read_slack_bus = False
        provider["slackBusSelectionMode"] = "NAME"
        provider["slackBusesIds"] = slack_bus_id
    params.provider_parameters = provider
    return params


def solve_olf(path, params, extra_load_mw=0.):
    net = load(path, extra_load_mw)
    res = lf.run_ac(net, params)
    return net, res[0]


def _enable_max_voltage_change(model, max_dva=1.0, max_dvm=0.4):
    """OpenLoadFlow's ``MAX_VOLTAGE_CHANGE`` step scaling, with this pypowsybl's bounds."""
    cfg = model.get_ac_algo_config()
    ip = list(cfg.int_params)
    ip[0] = int(ScalingPolicyType.MaxVoltageChange)
    cfg.int_params = ip
    rp = list(cfg.real_params)
    rp[0] = max_dva
    rp[1] = max_dvm
    cfg.real_params = rp
    model.set_ac_algo_config(cfg)


#: the lightsim2grid class of each OpenLoadFlow loop implemented so far
LIGHTSIM_LOOPS = {
    "DistributedSlack": lambda: _algorithm.DistributedSlack(),
}


def solve_lightsim(path, gen_slack_id, algo=None, max_iter=50, tol=1e-8, olf_rules=None, loops=(),
                   extra_load_mw=0.):
    """Build the grid from the unbaked snapshot and solve it from a DC start. ``olf_rules``
    (an ``OlfLoadingParameters``, or None for none) are OpenLoadFlow's loading rules: OLF
    always applies them, so a comparison needs them on this side too."""
    net = load(path, extra_load_mw)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        model = init_from_pypowsybl(net, gen_slack_id=gen_slack_id, sort_index=False,
                                    buses_for_sub=False, keep_half_open_lines=True,
                                    fuse_zero_impedance_branches=True, init_vm_pu=1.,
                                    olf_rules=olf_rules if olf_rules is not None else False)
    if loops:
        missing = [name for name in loops if name not in LIGHTSIM_LOOPS]
        if missing:
            raise ValueError(f"no lightsim2grid outer loop for {missing} yet")
        model.change_algorithm(algo if algo is not None else "NROuter_KLU")
        model.clear_outer_loops()
        for name in OLF_ORDER:
            if name in loops:
                model.add_outer_loop(LIGHTSIM_LOOPS[name]())
    elif algo is not None:
        model.change_algorithm(algo)
    model.set_keep_vinit_at_group_controlled_buses(True)  # before any dc_pf
    _enable_max_voltage_change(model)
    v_start = dc_start(model, distributed="DistributedSlack" in loops, max_iter=max_iter, tol=tol)
    V = model.ac_pf(v_start, max_iter, tol)
    return net, model, V


def dc_start(model, distributed, max_iter=50, tol=1e-8):
    """OpenLoadFlow's DC_VALUES start: the angles of a DC load flow, from a flat 1 pu (a
    bus regulated remotely keeps it, see set_keep_vinit_at_group_controlled_buses). With
    the distributed slack on, OpenLoadFlow's DC load flow already shares the imbalance on
    the participating units: the same here (set_dc_distribute_slack_on_can_participate),
    where a single-slack DC would leave it all on the slack generator -- a start far enough
    from OpenLoadFlow's to reach another root of a weak area."""
    flat = np.full(model.get_bus_vn_kv().shape[0], 1.0, dtype=complex)
    model.set_dc_distribute_slack_on_can_participate(distributed)
    v_dc = model.dc_pf(flat, max_iter, tol)
    return v_dc if v_dc.shape[0] > 0 else flat


def _rotated(angles_deg, ref):
    return angles_deg - angles_deg.loc[ref] if ref in angles_deg.index else angles_deg


def compare(net_olf, olf_result, net_ls, model, V, ref_bus, slack_gen_id):
    """Worst differences between the two solved states, on what both engines publish.

    The slack: OpenLoadFlow reports what is left on its slack bus in ``slack_bus_results``
    and leaves every generator at its target, lightsim2grid books it on its slack generator;
    both are compared as one number, and that generator is left out of the per-unit P.

    The reactive power is compared per bus (the sum over the bus' generators, which is what
    the solve decides) and per generator (which also depends on how each engine splits a
    bus' Q between its units). Never read a Q difference below ``Q_PRECISION_MVAR``: OLF
    dispatches nothing on a bus whose Q left to share is at most ``Q_DISPATCH_EPSILON``
    (a hard-coded 1e-5 pu, so 1e-3 MVar) and publishes 0 there, however exact the solve."""
    out = {}
    olf_v = iidm_bus_voltages(net_olf)
    # through the result view: it reads a bus merged by fuse_zero_impedance_branches at the
    # bus it was merged into (the merged one is left unsolved, at its starting voltage)
    ls_net = LightsimResultNetwork(model, net_ls)
    ls_bus = ls_net.get_buses()
    nominal = net_ls.get_voltage_levels()["nominal_v"]
    ls_v = pd.DataFrame({"vm_pu": ls_bus["v_mag"] / ls_bus["voltage_level_id"].map(nominal),
                         "va_deg": ls_bus["v_angle"]}).dropna()
    common = olf_v.index.intersection(ls_v.index)
    out["n_bus"] = len(common)
    dvm = (olf_v.loc[common, "vm_pu"] - ls_v.loc[common, "vm_pu"]).abs()
    dva = (_rotated(olf_v.loc[common, "va_deg"], ref_bus) - _rotated(ls_v.loc[common, "va_deg"], ref_bus)).abs()
    out["max_dvm_pu"] = float(dvm.max())
    out["worst_vm_bus"] = dvm.idxmax()
    out["max_dva_deg"] = float(dva.max())
    out["worst_va_bus"] = dva.idxmax()

    ls_gen = ls_net.get_generators()
    olf_slack = sum(r.active_power_mismatch for r in olf_result.slack_bus_results)
    # OLF publishes every unit at its (final) target and reports the leftover apart;
    # lightsim2grid books the leftover on its slack generator: compared on that generator
    ls_slack = -ls_gen.loc[slack_gen_id, "p"] + net_olf.get_generators().loc[slack_gen_id, "p"]
    out["olf_slack_mw"] = float(olf_slack)
    out["d_slack_mw"] = float(abs(olf_slack - ls_slack))

    for what, olf_df, ls_df in (("gen", net_olf.get_generators(), ls_gen),
                                ("vsc", net_olf.get_vsc_converter_stations(), ls_net.get_vsc_converter_stations())):
        olf_df = olf_df[olf_df["connected"] & olf_df["bus_id"].isin(common)] if len(olf_df) else olf_df
        ids = olf_df.index.intersection(ls_df.index)
        for col in ("p", "q"):
            keep = ids.drop(slack_gen_id, errors="ignore") if col == "p" else ids
            diff = (olf_df.loc[keep, col] - ls_df.loc[keep, col]).abs().dropna()
            out[f"max_d{col}_{what}"] = float(diff.max()) if len(diff) else 0.
            out[f"worst_{col}_{what}"] = diff.idxmax() if len(diff) else ""
        if what == "gen":
            q = pd.DataFrame({"olf": olf_df.loc[ids, "q"], "ls": ls_df.loc[ids, "q"],
                              "bus": olf_df.loc[ids, "bus_id"]}).groupby("bus").sum()
            diff = (q["olf"] - q["ls"]).abs()
            out["max_dq_bus"] = float(diff.max()) if len(diff) else 0.
            out["worst_q_bus"] = diff.idxmax() if len(diff) else ""
    return out


def _q(value):
    """A reactive power difference, as printed: nothing below OLF's dispatch threshold."""
    return f"<{Q_PRECISION_MVAR:g}" if value <= Q_PRECISION_MVAR else f"{value:.1e}"


def run_one(path, args):
    row = {"snapshot": os.path.basename(path)}
    net0 = load(path)
    gen_slack_id, slack_bus_id = pick_slack(net0)
    row["slack_gen"] = gen_slack_id

    t0 = time.perf_counter()
    params = olf_parameters(args.loops, slack_bus_id if args.slack_mode == "name" else None,
                            conv_eps=args.olf_eps, max_nr_iter=args.max_iter)
    net_olf, res = solve_olf(path, params, args.extra_load_mw)
    row["olf_status"] = res.status.name
    row["olf_status_text"] = res.status_text
    row["olf_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    # the loading rules depend on use_reactive_limits, set as on the OLF side
    rules = OlfLoadingParameters(reactive_limits=params.use_reactive_limits) if args.olf_rules else None
    net_ls, model, V = solve_lightsim(path, gen_slack_id, args.algo, args.max_iter, args.tol, rules, args.loops,
                                      args.extra_load_mw)
    row["ls_converged"] = V.shape[0] > 0
    row["ls_s"] = time.perf_counter() - t0
    stats = model.get_algo().get_linear_solver_stats() if hasattr(model.get_algo(), "get_linear_solver_stats") else None
    if stats is not None:
        row["ls_nb_analyze"] = stats.nb_analyze
        row["ls_nb_factorize"] = stats.nb_factorize
    if model.get_algo().supports_outer_loops():
        outer = model.get_algo().get_outer_loop_stats()
        row["ls_outer_status"] = outer.status.name
        row["ls_outer_iterations"] = outer.nb_outer_iterations
    if row["ls_converged"] and res.status == lf.ComponentStatus.CONVERGED:
        row.update(compare(net_olf, res, net_ls, model, V, slack_bus_id, gen_slack_id))
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="*", help="snapshots or directories (default: $LS2G_SNAPSHOTS)")
    parser.add_argument("--loops", default="", type=lambda s: [x for x in s.split(",") if x],
                        help=f"comma-separated OpenLoadFlow loops to run, among {','.join(OLF_ORDER)}")
    parser.add_argument("--algo", default=None, help="lightsim2grid algorithm name (default: the grid's)")
    parser.add_argument("--slack-mode", choices=("name", "most_meshed"), default="name")
    parser.add_argument("--olf-eps", type=float, default=1e-9, help="OpenLoadFlow newtonRaphsonConvEpsPerEq")
    parser.add_argument("--max-iter", type=int, default=50, help="Newton iterations, both engines")
    parser.add_argument("--tol", type=float, default=1e-8, help="lightsim2grid ac_pf tolerance")
    parser.add_argument("--no-olf-rules", dest="olf_rules", action="store_false",
                        help="build the lightsim2grid grid without OpenLoadFlow's loading rules")
    parser.add_argument("--extra-load-mw", type=float, default=0.,
                        help="added to the largest load of the main component, on both engines")
    parser.add_argument("--csv", default=None, help="write one row per snapshot there")
    parser.add_argument("--limit", type=int, default=None, help="only the first N snapshots")
    args = parser.parse_args(argv)

    paths = args.paths or ([os.environ["LS2G_SNAPSHOTS"]] if "LS2G_SNAPSHOTS" in os.environ else [])
    snapshots = find_snapshots(paths)[:args.limit]
    if not snapshots:
        parser.error("no snapshot found (give a path or set $LS2G_SNAPSHOTS)")

    rows = []
    for path in snapshots:
        try:
            row = run_one(path, args)
        except Exception as exc:  # one broken snapshot must not stop the sweep
            row = {"snapshot": os.path.basename(path), "error": f"{type(exc).__name__}: {exc}"}
        rows.append(row)
        if "error" in row:
            print(f"{row['snapshot']}: ERROR {row['error']}", flush=True)
        elif "max_dvm_pu" not in row:
            print(f"{row['snapshot']}: OLF {row['olf_status']} ({row['olf_status_text']}), "
                  f"lightsim2grid converged={row['ls_converged']}", flush=True)
        else:
            print(f"{row['snapshot']}: dVm={row['max_dvm_pu']:.1e} pu ({row['worst_vm_bus']}), "
                  f"dVa={row['max_dva_deg']:.1e} deg, slack {row['olf_slack_mw']:.2f} MW "
                  f"(d={row['d_slack_mw']:.1e}), dP gen={row['max_dp_gen']:.1e} MW, "
                  f"dQ bus={_q(row['max_dq_bus'])} / gen={_q(row['max_dq_gen'])} MVar, "
                  f"OLF {row['olf_s']:.1f}s / ls {row['ls_s']:.1f}s", flush=True)

    if args.csv:
        columns = sorted({k for row in rows for k in row}, key=lambda k: (k != "snapshot", k))
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
