# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

# column indices, 0-based, matching MATPOWER's lib/idx_bus.m
(BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
 BASE_KV, ZONE, VMAX, VMIN) = range(13)

# bus type constants, matching MATPOWER's lib/idx_bus.m
PQ, PV, REF, NONE = 1, 2, 3, 4

# column indices, 0-based, matching MATPOWER's lib/idx_gen.m
(GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS,
 PMAX, PMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX,
 RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF) = range(21)

# column indices, 0-based, matching MATPOWER's lib/idx_brch.m
(F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C,
 TAP, SHIFT, BR_STATUS, ANGMIN, ANGMAX) = range(13)

# minimum number of columns lightsim2grid needs to find in each matpower
# matrix (matpower case files may have fewer optional trailing columns,
# e.g. no ANGMIN / ANGMAX, or more, e.g. solved-case result columns)
MIN_BUS_COLS = VMIN + 1
MIN_GEN_COLS = PMIN + 1
MIN_BRANCH_COLS = BR_STATUS + 1
