# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Benchmark comparing LSGrid.save_binary / LSGrid.load_binary (the fast, additive
binary serialization path) against pickle.dump / pickle.load on the same grid.

Also reports, as a reference point, the "reading time" of building an LSGrid
from scratch from a pandapower network (init_from_pandapower) and from a
pypowsybl network (init_from_pypowsybl): these are NOT round trips of an
already-converted grid, they show how much is saved by re-loading a
lightsim2grid-native snapshot instead of re-running the full conversion.

In addition to the single l2rpn_idf_2023 fixture, this also scans a range of
grid sizes (the same pandapower cases used by benchmark_grid_size.py, from
14 to ~9241 buses) to see how the save_binary/load_binary speedup over pickle
evolves with grid size.

Usage: python benchmark_binary_serialization.py [--n_repeat N] [--n_repeat_convert N] [--cases case14.json,case118.json,...]
"""

import argparse
import os
import pickle
import tempfile
import time
import warnings

import grid2op
from lightsim2grid.lightSimBackend import LightSimBackend

TABULATE_AVAIL = False
try:
    from tabulate import tabulate
    TABULATE_AVAIL = True
except ImportError:
    print("The tabulate package is not installed. Output will use a plain format.")

ENV_NAME = "l2rpn_idf_2023"

# same graduated set of pandapower cases used by benchmark_grid_size.py, all
# already cached as json files in this directory
CASE_NAMES = [
    "case14.json",
    "case118.json",
    "case_illinois200.json",
    "case300.json",
    "case1354pegase.json",
    "case1888rte.json",
    "case2848rte.json",
    "case2869pegase.json",
    "case3120sp.json",
    "case6495rte.json",
    "case6515rte.json",
    "case9241pegase.json",
]


def _time_repeat(fun, n_repeat):
    """Run `fun()` `n_repeat` times, return (total_seconds, per_call_seconds)."""
    start = time.perf_counter()
    for _ in range(n_repeat):
        fun()
    total = time.perf_counter() - start
    return total, total / n_repeat


def benchmark_save_load(grid, n_repeat, tmpdir):
    binary_path = os.path.join(tmpdir, "benchmark.lsb")
    pickle_path = os.path.join(tmpdir, "benchmark.pickle")

    _, save_binary_s = _time_repeat(lambda: grid.save_binary(binary_path), n_repeat)
    _, load_binary_s = _time_repeat(lambda: type(grid).load_binary(binary_path), n_repeat)
    binary_size = os.path.getsize(binary_path)

    def _pickle_dump():
        with open(pickle_path, "wb") as f:
            pickle.dump(grid, f)

    def _pickle_load():
        with open(pickle_path, "rb") as f:
            pickle.load(f)

    _, pickle_dump_s = _time_repeat(_pickle_dump, n_repeat)
    _, pickle_load_s = _time_repeat(_pickle_load, n_repeat)
    pickle_size = os.path.getsize(pickle_path)

    return {
        "save_binary_s": save_binary_s,
        "load_binary_s": load_binary_s,
        "binary_size_bytes": binary_size,
        "pickle_dump_s": pickle_dump_s,
        "pickle_load_s": pickle_load_s,
        "pickle_size_bytes": pickle_size,
    }


def benchmark_init_from_pandapower(n_repeat):
    try:
        import pandapower.networks as pn
        from lightsim2grid.network import init_from_pandapower
    except ImportError:
        return None
    net = pn.case118()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        _, read_s = _time_repeat(lambda: init_from_pandapower(net), n_repeat)
    return read_s


def benchmark_init_from_pypowsybl(n_repeat):
    try:
        import pypowsybl.network as pp_network
        from lightsim2grid.network import init_from_pypowsybl
    except ImportError:
        return None
    net = pp_network.create_ieee118()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        _, read_s = _time_repeat(lambda: init_from_pypowsybl(net), n_repeat)
    return read_s


def run_fixed_grid_benchmark(n_repeat):
    """Original single-grid comparison: save_binary/load_binary vs pickle on
    the l2rpn_idf_2023 fixture (same grid used by test_binary_serialization.py
    and test_pickleable.py), plus pandapower/pypowsybl reference reading times.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        env = grid2op.make(ENV_NAME, test=True, backend=LightSimBackend())
    grid = env.backend._grid

    with tempfile.TemporaryDirectory() as tmpdir:
        res = benchmark_save_load(grid, n_repeat, tmpdir)

    pandapower_read_s = benchmark_init_from_pandapower(n_repeat)
    pypowsybl_read_s = benchmark_init_from_pypowsybl(n_repeat)

    rows = [
        ["save_binary (write)", f"{1000. * res['save_binary_s']:.3f}", f"{res['binary_size_bytes']} bytes"],
        ["load_binary (read)", f"{1000. * res['load_binary_s']:.3f}", "-"],
        ["pickle.dump (write)", f"{1000. * res['pickle_dump_s']:.3f}", f"{res['pickle_size_bytes']} bytes"],
        ["pickle.load (read)", f"{1000. * res['pickle_load_s']:.3f}", "-"],
    ]
    if pandapower_read_s is not None:
        rows.append(["init_from_pandapower (case118, reference)", f"{1000. * pandapower_read_s:.3f}", "-"])
    else:
        rows.append(["init_from_pandapower", "pandapower not installed, skipped", "-"])
    if pypowsybl_read_s is not None:
        rows.append(["init_from_pypowsybl (ieee118, reference)", f"{1000. * pypowsybl_read_s:.3f}", "-"])
    else:
        rows.append(["init_from_pypowsybl", "pypowsybl not installed, skipped", "-"])

    headers = [f"{ENV_NAME} ({n_repeat} repeats)", "time (ms/call)", "file size"]
    if TABULATE_AVAIL:
        print(tabulate(rows, headers=headers, tablefmt="rst"))
    else:
        print("\t".join(headers))
        for row in rows:
            print("\t".join(str(x) for x in row))

    speedup_write = res["pickle_dump_s"] / res["save_binary_s"]
    speedup_read = res["pickle_load_s"] / res["load_binary_s"]
    print(f"\nsave_binary is {speedup_write:.1f}x faster than pickle.dump")
    print(f"load_binary is {speedup_read:.1f}x faster than pickle.load")


def _load_pandapower_case(case_name):
    """Same load-or-fetch-and-cache pattern as benchmark_grid_size.py: the
    case*.json files are already committed in this directory, this is only a
    fallback."""
    import pandapower as pp
    if not os.path.exists(case_name):
        import pandapower.networks as pn
        case = getattr(pn, os.path.splitext(case_name)[0])()
        pp.to_json(case, case_name)
    return pp.from_json(case_name)


def _time_read_pandapower_file(case_name, n_repeat):
    """Time the full 'read a pandapower json FILE from disk and build an
    LSGrid from it' pipeline (pp.from_json + init_from_pandapower). This is
    what's directly comparable to load_binary/pickle.load: both of those also
    start from a file on disk, not an already-loaded in-memory object, so the
    pandapower reference needs to include the file read to be apples-to-apples.
    """
    import pandapower as pp
    from lightsim2grid.network import init_from_pandapower

    def _read_and_convert():
        net = pp.from_json(case_name)
        init_from_pandapower(net)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        _, read_s = _time_repeat(_read_and_convert, n_repeat)
    return read_s


def run_grid_size_scan(n_repeat, n_repeat_convert, case_names):
    """save_binary/load_binary vs pickle across a range of grid sizes, built
    directly via init_from_pandapower (no grid2op env needed, just an LSGrid)
    -- same pandapower cases benchmark_grid_size.py uses to scan solver speed.

    Also reports, as a reference, the time (and file size) to read the
    original pandapower json FILE and build an LSGrid from it from scratch --
    this is the "no snapshot at all" baseline that save_binary/load_binary and
    pickle are meant to be faster than for repeated re-loads of the same grid.
    """
    try:
        from lightsim2grid.network import init_from_pandapower
    except ImportError:
        print("\npandapower not installed, skipping the grid-size scan.")
        return

    rows = []
    for case_name in case_names:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            net = _load_pandapower_case(case_name)
            grid = init_from_pandapower(net)

        pp_read_s = _time_read_pandapower_file(case_name, n_repeat_convert)

        with tempfile.TemporaryDirectory() as tmpdir:
            res = benchmark_save_load(grid, n_repeat, tmpdir)

        speedup_write = res["pickle_dump_s"] / res["save_binary_s"]
        speedup_read = res["pickle_load_s"] / res["load_binary_s"]
        speedup_read_vs_pp_file = pp_read_s / res["load_binary_s"]

        rows.append([
            os.path.splitext(case_name)[0],
            grid.total_bus(),
            f"{1000. * res['save_binary_s']:.3f}",
            f"{1000. * res['pickle_dump_s']:.3f}",
            f"{speedup_write:.1f}x",
            f"{1000. * res['load_binary_s']:.3f}",
            f"{1000. * res['pickle_load_s']:.3f}",
            f"{speedup_read:.1f}x",
            res["binary_size_bytes"],
            res["pickle_size_bytes"],
            os.path.getsize(case_name),
            f"{1000. * pp_read_s:.2f}",
            f"{speedup_read_vs_pp_file:.1f}x",
        ])

    headers = [
        "grid", "nb bus",
        "save_binary (ms)", "pickle.dump (ms)", "write speedup",
        "load_binary (ms)", "pickle.load (ms)", "read speedup",
        "binary size (B)", "pickle size (B)", "pandapower json size (B)",
        "read pandapower file (ms, reference)", "load_binary speedup vs pandapower file",
    ]
    print(f"\nGrid size scan ({n_repeat} repeats for save/load, {n_repeat_convert} for the "
          f"'read pandapower file' reference):")
    if TABULATE_AVAIL:
        print(tabulate(rows, headers=headers, tablefmt="rst"))
    else:
        print("\t".join(headers))
        for row in rows:
            print("\t".join(str(x) for x in row))


def main(n_repeat=20, n_repeat_convert=5, case_names=CASE_NAMES):
    run_fixed_grid_benchmark(n_repeat)
    run_grid_size_scan(n_repeat, n_repeat_convert, case_names)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n_repeat", type=int, default=20,
                        help="number of times save_binary/load_binary/pickle are repeated")
    parser.add_argument("--n_repeat_convert", type=int, default=5,
                        help="number of times the pandapower/pypowsybl conversion (reference row) is repeated")
    parser.add_argument("--cases", type=str, default=None,
                        help="comma separated list of case*.json files to scan (default: the full CASE_NAMES list)")
    args = parser.parse_args()
    case_names = CASE_NAMES if args.cases is None else args.cases.split(",")
    main(n_repeat=args.n_repeat, n_repeat_convert=args.n_repeat_convert, case_names=case_names)
