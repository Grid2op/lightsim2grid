#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Candidate change under A/B test (ab_test.sh): take 1/|V| from `Vm_` instead of
recomputing |V| with a `std::hypot` per bus per Newton iteration.

`Vm_` IS `|V_|`: the only two writers of `V_` keep it so. `update_state` sets
`Vm_ = |V_init|` next to `V_ = V_init`, and `apply_step` rebuilds
`V_ = Vm_ * exp(i.Va_)` and then makes `Vm_` non-negative with
`fix_negative_vm` (which flips the sign of `Vm_` and shifts `Va_` by half a
turn, leaving `V_` unchanged and equal to `Vm_new * exp(i.Va_new)`).

The two are therefore equal up to the rounding of `|cos + i.sin|`, not bit for
bit -- which is exactly what the A/B checks: same iteration counts, same
returned voltages.
"""

import os
import sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "..", "src", "core", "powerflow_algorithm",
                    "NRSystem.tpp")

OLD = "    inv_vm_cache_ = V_.array().abs().inverse();  // 1 / |V|"
NEW = ("    // 1 / |V|, taken from Vm_ (which IS |V_|) rather than a hypot pass over V_\n"
       "    inv_vm_cache_ = Vm_.array().inverse();")


def main():
    with open(PATH) as fp:
        src = fp.read()
    if NEW in src:
        print("inv_vm_from_vm: already applied")
        return 0
    if OLD not in src:
        print(f"inv_vm_from_vm: anchor not found in {PATH}", file=sys.stderr)
        return 1
    with open(PATH, "w") as fp:
        fp.write(src.replace(OLD, NEW))
    print("inv_vm_from_vm: applied")
    return 0


if __name__ == "__main__":
    sys.exit(main())
