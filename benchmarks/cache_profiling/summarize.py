#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Turn the callgrind output of run_profile.sh into the two tables the audit needs:

  * per grid / per phase, the instructions retired PER SOLVE, split between the
    powerflow algorithm (NR + the linear solver) and everything the LSGrid does
    around it;
  * inside LSGrid.cpp, the same instructions grouped by the function whose
    source lines they were spent on -- which is what says where the cached path
    still pays, and for what.

    python summarize.py <callgrind_out_dir> [src_root]
"""

import os
import re
import sys
from collections import OrderedDict, defaultdict

# the single-solve phases (run_profile.sh) ...
SINGLE_PHASES = ["cold", "idem", "inj", "inj_nores", "dcac", "nocache", "topo",
                 "inj_every2", "inj_every3", "idem__NRSing_KLU", "inj__NRSing_KLU"]
# ... and the batch ones (run_profile_batch.sh), where "per solve" means per ROW.
# Every compute phase comes with its `_1row` twin: the same phase over a single
# row, which is the fixed cost of a compute() -- see the marginal-cost table.
BATCH_COMPUTE_PHASES = ["ts_ac", "ts_dc", "ca_ac", "ca_dc", "ca_ac_mask", "ca_dc_mask"]
BATCH_PHASES = (BATCH_COMPUTE_PHASES
                + [p + "_1row" for p in BATCH_COMPUTE_PHASES]
                + ["ts_flows", "ca_flows", "ca_construct"])
# the phases the per-call-site ledger is printed for
LEDGER_PHASES = ["idem", "inj", "nocache", "ts_ac", "ts_dc", "ca_ac", "ca_ac_mask",
                 "ca_dc_mask", "ts_flows", "ca_flows"]
GRID_ORDER = ["case30", "case118", "case1354pegase", "case9241pegase",
              "case118_fancy", "case1354pegase_fancy", "case9241pegase_fancy"]

# inclusive costs worth naming, matched as substrings of the demangled symbol
INCLUSIVE_OF_INTEREST = OrderedDict([
    ("ac_pf",            "LSGrid::ac_pf("),
    ("dc_pf",            "LSGrid::dc_pf("),
    ("batch compute",    "BaseBatchSweep<"),
    ("build_solver_in",  "LSGrid::build_solver_input("),
    ("build_dc_sol_in",  "LSGrid::build_dc_solver_input("),
    ("NR compute_pf",    "NRAlgo<"),
    ("DC compute_pf",    "BaseDCAlgo<"),
    ("KLU analyze",      "KLULinearSolver::analyze("),
    ("KLU factorize",    "KLULinearSolver::factorize("),
    ("KLU refactorize",  "KLULinearSolver::refactorize("),
    ("KLU solve",        "KLULinearSolver::solve("),
    ("fill_J",           "::fill_J("),
    ("fill_internals",   "::fill_internal_variables("),
    ("check_invertible", "Contingency::check_invertible("),
    ("disconnected_bus", "::_disconnected_buses"),
    ("select_ref_slack", "::_select_ref_slack_and_masks"),
    ("Sbus assemble",    "SbusPolicy::Vary::assemble("),
    ("amps flows",       "::compute_amps_flows"),
    ("power flows",      "::compute_active_power_flows"),
    ("violations",       "batch_sweep_detail::check_"),
    ("LSGrid copy",      "LSGrid::LSGrid(ls2g::LSGrid const&)"),
])

# the source files the per-call-site ledger is printed for, relative to the
# repository root: what LSGrid.cpp does around a single solve, and what the
# batch classes do around each row
LEDGER_FILES = ["src/core/LSGrid.cpp",
                "src/core/batch_algorithm/BaseBatchSweep.cpp",
                "src/core/batch_algorithm/BaseBatchSweep.hpp",
                "src/core/batch_algorithm/BaseBatchSolverSynch.hpp",
                "src/core/batch_algorithm/BaseBatchSolverSynch.cpp",
                "src/core/batch_algorithm/SbusPolicy.cpp",
                "src/core/batch_algorithm/YbusPolicy.cpp"]

_NUM = re.compile(r"^\s*([\d,]+)\s")


def _read_total(path):
    with open(path) as fp:
        for line in fp:
            if "PROGRAM TOTALS" in line:
                m = _NUM.match(line)
                if m:
                    return int(m.group(1).replace(",", ""))
    return None


def _read_inclusive(path):
    """symbol -> inclusive Ir.

    callgrind_annotate lists the same function once per source file its code was
    attributed to; the true inclusive cost is the largest of them.
    """
    out = defaultdict(int)
    with open(path) as fp:
        for line in fp:
            m = _NUM.match(line)
            if not m or "file:function" in line or "PROGRAM TOTALS" in line:
                continue
            cost = int(m.group(1).replace(",", ""))
            rest = line[m.end():]
            # strip the "(xx.xx%)" column and the trailing "[/path/to.so]"
            rest = re.sub(r"^\s*\([^)]*\)\s*", "", rest).strip()
            rest = re.sub(r"\s*\[[^\]]*\]\s*$", "", rest).strip()
            if ":" in rest:
                rest = rest.split(":", 1)[1]
            out[rest] = max(out[rest], cost)
    return out


def _pick(inclusive, needle):
    best = 0
    for sym, cost in inclusive.items():
        if needle in sym:
            best = max(best, cost)
    return best


def _function_starts(src_path):
    """line number -> function name, for the definitions in a .cpp file.

    Definitions in LSGrid.cpp all start at column 0 with `<ret> LSGrid::name(`,
    which is enough to bucket annotated source lines by the function they are in
    -- including the code the compiler inlined, since an inlined callee's
    instructions are still attributed to ITS own source lines.
    """
    starts = []
    pat = re.compile(r"^[A-Za-z_].*?\b([A-Za-z_]\w*)::([~A-Za-z_]\w*)\s*\(")
    with open(src_path) as fp:
        for lineno, line in enumerate(fp, 1):
            m = pat.match(line)
            if m:
                starts.append((lineno, f"{m.group(1)}::{m.group(2)}"))
    return starts


def _enclosing(starts, lineno):
    name = "<file scope>"
    for start, fname in starts:
        if start <= lineno:
            name = fname
        else:
            break
    return name


def _present(out_dir, phases):
    """the phases of `phases` (in that order) that have at least one annotate file"""
    names = os.listdir(out_dir)
    return [p for p in phases
            if any(n.startswith("annotate.") and n.endswith(f".{p}.txt") for n in names)]


def main(out_dir, src_root):
    PHASES = _present(out_dir, SINGLE_PHASES + BATCH_PHASES)
    grids = [g for g in GRID_ORDER
             if any(os.path.exists(os.path.join(out_dir, f"annotate.{g}.{p}.txt")) for p in PHASES)]
    print("=" * 100)
    print("INSTRUCTIONS RETIRED (Ir) PER POWERFLOW -- collected region only "
          "(the ac_pf / dc_pf calls themselves; per ROW for the batch phases)")
    print("=" * 100)
    header = f"{'grid':<22}" + "".join(f"{p:>18}" for p in PHASES)
    print(header)
    totals = {}
    for grid in grids:
        row = f"{grid:<22}"
        for phase in PHASES:
            ann = os.path.join(out_dir, f"annotate.{grid}.{phase}.txt")
            nbf = os.path.join(out_dir, f"nb.{grid}.{phase}.txt")
            if not os.path.exists(ann):
                row += f"{'-':>18}"
                continue
            total = _read_total(ann)
            nb = int(open(nbf).read().strip()) if os.path.exists(nbf) else 1
            per = total / nb
            totals[(grid, phase)] = per
            totals[(grid, phase, "raw")] = (total, nb)
            row += f"{per:>18,.0f}"
        print(row)

    # ---- batch: the fixed cost of a compute() apart from the marginal row ------
    batch_compute = [p for p in BATCH_COMPUTE_PHASES if p in PHASES and p + "_1row" in PHASES]
    if batch_compute:
        print()
        print("=" * 100)
        print("BATCH: what a compute() pays BEFORE its first row (the 1-row run), and "
              "what each further row costs\n(the N-row run minus the 1-row run, over N-1)")
        print("=" * 100)
        header = f"{'grid':<22}" + "".join(f"{p + ' fixed':>18}{p + ' /row':>18}" for p in batch_compute)
        print(header)
        for grid in grids:
            row = f"{grid:<22}"
            for phase in batch_compute:
                n_run = totals.get((grid, phase, "raw"))
                one = totals.get((grid, phase + "_1row", "raw"))
                if n_run is None or one is None or n_run[1] <= 1:
                    row += f"{'-':>18}{'-':>18}"
                    continue
                marginal = (n_run[0] - one[0]) / (n_run[1] - one[1])
                row += f"{one[0]:>18,.0f}{marginal:>18,.0f}"
            print(row)

    print()
    print("=" * 100)
    print("WHERE THOSE INSTRUCTIONS GO (inclusive Ir per solve)")
    print("=" * 100)
    for grid in grids:
        for phase in PHASES:
            inc_path = os.path.join(out_dir, f"inclusive.{grid}.{phase}.txt")
            nbf = os.path.join(out_dir, f"nb.{grid}.{phase}.txt")
            if not os.path.exists(inc_path):
                continue
            nb = int(open(nbf).read().strip()) if os.path.exists(nbf) else 1
            inclusive = _read_inclusive(inc_path)
            per_total = totals.get((grid, phase), 0)
            print(f"\n--- {grid} / {phase} ({per_total:,.0f} Ir per solve) ---")
            ac = _pick(inclusive, INCLUSIVE_OF_INTEREST["ac_pf"]) / nb
            dc = _pick(inclusive, INCLUSIVE_OF_INTEREST["dc_pf"]) / nb
            for label, needle in INCLUSIVE_OF_INTEREST.items():
                val = _pick(inclusive, needle) / nb
                if val <= 0:
                    continue
                pct = 100. * val / per_total if per_total else 0.
                print(f"    {label:<18} {val:>14,.0f}  ({pct:5.1f}%)")
            algo = (_pick(inclusive, INCLUSIVE_OF_INTEREST["NR compute_pf"]) +
                    _pick(inclusive, INCLUSIVE_OF_INTEREST["DC compute_pf"])) / nb
            if phase in BATCH_PHASES:
                # a batch phase collects the batch call itself, not ac_pf / dc_pf:
                # "around" is then everything the batch does that is not a solve
                around = per_total - algo
                label = "--> batch around"
            else:
                around = (ac + dc) - algo
                label = "--> LSGrid around"
            if per_total:
                print(f"    {'--> algorithms':<18} {algo:>14,.0f}  "
                      f"({100. * algo / per_total:5.1f}%)")
                print(f"    {label:<18} {around:>14,.0f}  "
                      f"({100. * around / per_total:5.1f}%)")

    # ---- the ledger: what each call site inside the files of interest costs ----
    # A single-solve phase is read against LSGrid.cpp, a batch phase against the
    # batch sources: a batch phase never enters ac_pf, and the single solve never
    # enters the batch classes, so the other family's files would list nothing.
    lsgrid_cpp = os.path.join(src_root, "src", "core", "LSGrid.cpp")
    starts = _function_starts(lsgrid_cpp) if os.path.exists(lsgrid_cpp) else []
    for grid in grids:
        for phase in LEDGER_PHASES:
            path = os.path.join(out_dir, f"callgrind.{grid}.{phase}.out")
            nbf = os.path.join(out_dir, f"nb.{grid}.{phase}.txt")
            if not os.path.exists(path):
                continue
            nb = int(open(nbf).read().strip()) if os.path.exists(nbf) else 1
            self_costs, call_costs, self_by_fn = parse_callgrind(path)
            is_batch = phase in BATCH_PHASES
            files = [f for f in LEDGER_FILES if ("batch_algorithm" in f) == is_batch]
            basenames = tuple(os.path.basename(f) for f in files)

            print()
            print("=" * 100)
            print(f"{grid} / {phase}: what each CALL SITE in {', '.join(basenames)} costs "
                  f"(inclusive Ir per solve)")
            print("=" * 100)
            sites = [(cost, f, ln, callee)
                     for (f, ln, callee), cost in call_costs.items()
                     if f and f.endswith(basenames)]
            sites.sort(reverse=True)
            for cost, f, ln, callee in sites[:18]:
                print(f"    {cost / nb:>12,.0f}  {os.path.basename(f)[:22]:<22} L{ln:<5} "
                      f"{_src_line(src_root, f, ln)[:50]:<50} -> {_demangle_short(callee)[:40]}")

            if is_batch:
                # the batch code is mostly member templates defined in headers, so
                # bucket by the function callgrind attributes the instructions to
                # (an inlined callee's instructions are charged to the function it
                # was inlined into, on the callee's own source lines)
                print(f"    --- instructions executed IN those files themselves "
                      f"(no callee), by function ---")
                per_fn = defaultdict(int)
                for (f, fn), cost in self_by_fn.items():
                    if f and f.endswith(basenames):
                        per_fn[fn] += cost
                for fname, cost in sorted(per_fn.items(), key=lambda kv: -kv[1])[:12]:
                    print(f"    {cost / nb:>12,.0f}  {_demangle_short(fname)[:80]}")
            else:
                print(f"    --- instructions executed IN LSGrid.cpp itself "
                      f"(no callee), by function ---")
                per_fn = defaultdict(int)
                for (f, ln), cost in self_costs.items():
                    if f and f.endswith("LSGrid.cpp"):
                        per_fn[_enclosing(starts, ln)] += cost
                for fname, cost in sorted(per_fn.items(), key=lambda kv: -kv[1])[:10]:
                    print(f"    {cost / nb:>12,.0f}  {fname}")


def parse_callgrind(callgrind_out):
    """(self_costs, call_costs, self_by_fn): the first two keyed by (file, line),
    the third by (file, function name).

    The callgrind format writes, inside an `fn=` block, plain "<line> <Ir>" cost
    lines for the code of the function itself, and -- right after a `calls=`
    line -- one "<line> <Ir>" giving the INCLUSIVE cost of that call, charged to
    the line of the call site. Telling the two apart is the whole point: the
    first is where the instructions are executed, the second is the ledger of
    what each call site costs, which is what a breakdown of `ac_pf` is made of.

    `fl=` / `fi=` / `fe=` give the source file a cost line belongs to (`fi`/`fe`
    is how code inlined from another file is reported), and a position may be
    absolute ("1234"), relative ("+3", "-2") or "*" (same line as before).
    """
    self_costs = defaultdict(int)
    call_costs = defaultdict(int)
    self_by_fn = defaultdict(int)
    file_names = {}
    fn_names = {}
    cur_file = None
    cur_fn = None
    cur_callee = None
    last_line = 0
    pending_call = False

    def _resolve(table, val):
        m = re.match(r"^\((\d+)\)\s*(.*)$", val)
        if not m:
            return val
        fid, name = m.group(1), m.group(2)
        if name:
            table[fid] = name
        return table.get(fid)

    with open(callgrind_out) as fp:
        for raw in fp:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            key, sep, val = line.partition("=")
            if sep and key and not key[0].isdigit() and key[0] not in "+-*":
                if key in ("fl", "fi", "fe"):
                    cur_file = _resolve(file_names, val)
                    last_line = 0
                elif key in ("cfi", "cfl"):
                    _resolve(file_names, val)
                elif key == "fn":
                    cur_fn = _resolve(fn_names, val)
                    last_line = 0
                elif key == "cfn":
                    cur_callee = _resolve(fn_names, val)
                elif key == "calls":
                    pending_call = True
                continue
            parts = line.split()
            if len(parts) < 2 or not parts[1].lstrip("-").isdigit():
                continue
            pos, cost = parts[0], int(parts[1])
            if pos == "*":
                lineno = last_line
            elif pos[0] in "+-":
                lineno = last_line + int(pos)
            elif pos.isdigit():
                lineno = int(pos)
            else:
                continue
            last_line = lineno
            if pending_call:
                call_costs[(cur_file, lineno, cur_callee)] += cost
                pending_call = False
            else:
                self_costs[(cur_file, lineno)] += cost
                self_by_fn[(cur_file, cur_fn)] += cost
    return self_costs, call_costs, self_by_fn


def _src_line(src_root, path, lineno):
    full = path if os.path.isabs(path) else os.path.join(src_root, path)
    try:
        with open(full) as fp:
            for i, text in enumerate(fp, 1):
                if i == lineno:
                    return text.strip()
    except OSError:
        pass
    return ""


def _demangle_short(name):
    if not name:
        return "?"
    name = re.sub(r"<[^<>]*>", "<>", name)
    for _ in range(3):
        name = re.sub(r"<[^<>]*>", "<>", name)
    return name.split("(")[0]


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "callgrind_out",
         sys.argv[2] if len(sys.argv) > 2
         else os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
