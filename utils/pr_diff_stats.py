#!/usr/bin/env python3
# Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Split a branch's diff into C++ / Python / tests / docs and print it as a markdown table,
to put at the top of a pull request body.

A PR here routinely runs to one or two thousand lines, which is daunting to open and says
nothing about where the work actually is. The table says how much of it is code a reviewer
has to reason about and how much is tests, docs and changelog::

    python3 utils/pr_diff_stats.py                # against origin/dev_1.0.1 (the default)
    python3 utils/pr_diff_stats.py origin/main    # against another base branch

Counts the branch's own changes only (`git diff base...HEAD`), so commits that landed on
the base branch meanwhile are not attributed to the PR.
"""

import fnmatch
import subprocess
import sys

#: (bucket key, label), in the order the table prints them
BUCKETS = [
    ("cpp", "C++ (src/core, src/bindings)"),
    ("python", "Python (package)"),
    ("tests", "Tests"),
    ("docs", "Docs + changelog"),
    ("other", "Build / other"),
]

DEFAULT_BASE = "origin/dev_1.0.1"


def bucket_of(path):
    """Which row of the table a changed file belongs to.

    Tests come first on purpose: a test file is a test whatever its language, so
    ``src/tests/*.cpp`` must not fall into the C++ bucket.
    """
    if fnmatch.fnmatch(path, "src/tests/*") or fnmatch.fnmatch(path, "lightsim2grid/tests/*"):
        return "tests"
    if path.startswith("docs/") or path.endswith((".rst", ".md")):
        return "docs"
    if path.startswith(("src/core/", "src/bindings/")) and path.endswith((".cpp", ".hpp", ".tpp", ".h")):
        return "cpp"
    if path.startswith("lightsim2grid/") and path.endswith(".py"):
        return "python"
    return "other"


def diff_stats(base):
    """{bucket: [added, removed]} for `git diff base...HEAD`."""
    merge_base = subprocess.run(["git", "merge-base", base, "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    numstat = subprocess.run(["git", "diff", "--numstat", merge_base + "...HEAD"],
                             capture_output=True, text=True, check=True).stdout
    totals = {}
    for line in numstat.strip().splitlines():
        added, removed, path = line.split("\t")
        if added == "-":
            continue  # binary file: git reports no line counts
        entry = totals.setdefault(bucket_of(path), [0, 0])
        entry[0] += int(added)
        entry[1] += int(removed)
    return totals


def as_markdown(totals):
    lines = ["| | added | removed |", "|---|---|---|"]
    for key, label in BUCKETS:
        if key in totals:  # an empty bucket gets no row
            lines.append(f"| **{label}** | +{totals[key][0]} | -{totals[key][1]} |")
    added = sum(v[0] for v in totals.values())
    removed = sum(v[1] for v in totals.values())
    lines.append(f"| total | +{added} | -{removed} |")
    return "\n".join(lines)


def main():
    base = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BASE
    try:
        print(as_markdown(diff_stats(base)))
    except subprocess.CalledProcessError as exc:
        print(f"could not diff against '{base}': {exc.stderr.strip() or exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
