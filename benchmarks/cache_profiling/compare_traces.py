#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Compare the A and B traces written by ab_test.sh.

`cmp` only answers "bit for bit or not", which is the wrong question for a change
that perturbs the Jacobian by a rounding: what matters is whether the ITERATION
COUNTS moved and how far apart the returned voltages are, in ulps of the values
themselves.

    python3 compare_traces.py <ab_out_dir>
"""

import glob
import os
import re
import sys


def read_trace(path):
    """[(nb_iter, [(re, im), ...]), ...], one entry per solve."""
    steps = []
    cur = None
    with open(path) as fp:
        for line in fp:
            m = re.match(r"^step (\d+) iter (\d+)$", line.strip())
            if m:
                cur = (int(m.group(2)), [])
                steps.append(cur)
                continue
            re_s, im_s = line.split()
            cur[1].append((float(re_s), float(im_s)))
    return steps


def compare(a_path, b_path):
    a, b = read_trace(a_path), read_trace(b_path)
    if len(a) != len(b):
        return None, None, f"different number of solves ({len(a)} vs {len(b)})"
    iter_a = [s[0] for s in a]
    iter_b = [s[0] for s in b]
    max_abs = 0.0
    max_rel = 0.0
    for (ia, va), (ib, vb) in zip(a, b):
        if len(va) != len(vb):
            return None, None, "different vector sizes"
        for (ar, ai), (br, bi) in zip(va, vb):
            d = abs(complex(ar, ai) - complex(br, bi))
            mag = abs(complex(ar, ai))
            max_abs = max(max_abs, d)
            if mag > 0.0:
                max_rel = max(max_rel, d / mag)
    note = "same iteration counts" if iter_a == iter_b else \
           f"ITERATIONS DIFFER: A={sum(iter_a)} B={sum(iter_b)} total"
    return max_abs, max_rel, note


def main(out_dir):
    print(f"{'grid':<16} {'phase':<10} {'max |dV| (pu)':>16} {'max rel':>12}   note")
    for a_path in sorted(glob.glob(os.path.join(out_dir, "trace.A.*.txt"))):
        b_path = a_path.replace("trace.A.", "trace.B.")
        if not os.path.exists(b_path):
            continue
        name = os.path.basename(a_path)[len("trace.A."):-len(".txt")]
        grid, phase = name.split(".", 1)
        max_abs, max_rel, note = compare(a_path, b_path)
        if max_abs is None:
            print(f"{grid:<16} {phase:<10} {'-':>16} {'-':>12}   {note}")
        else:
            # 1 ulp of a per-unit voltage near 1.0 is 2.2e-16
            ulps = max_rel / 2.220446049250313e-16
            print(f"{grid:<16} {phase:<10} {max_abs:>16.3e} {max_rel:>12.3e}   "
                  f"{note} ({ulps:.1f} ulp)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ab_out")
