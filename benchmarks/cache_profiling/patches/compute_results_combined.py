#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""Both compute_results candidates at once, for one A/B run (ab_test.sh)."""

import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

rc = 0
for name in ("bus_mismatch_no_temporary.py", "branch_results_hoist.py"):
    sys.argv = [os.path.join(HERE, name)]
    try:
        # each script ends in sys.exit(main()), which would stop this one
        runpy.run_path(sys.argv[0], run_name="__main__")
    except SystemExit as exc:
        rc = rc or (exc.code or 0)
sys.exit(rc)
