#!/usr/bin/env bash
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
#
# The wall-clock half of an A/B (see ab_test.sh for the instruction-count half).
# Instruction counts do not know about caches or the latency of a libm call, so a
# change that removes arithmetic has to be confirmed on a real run:
#
#   ./ab_wallclock.sh <grids_dir> <out_dir> <patch_script> [repeats] [phases...]
#
# Each (grid, phase, variant) is run `repeats` times and the BEST run is kept --
# the minimum is the least noisy statistic here, a slower run only ever means
# more interference.

set -euo pipefail

GRIDS_DIR=$1
OUT_DIR=$2
PATCH=$3
REPEATS=${4:-7}
shift 4 2>/dev/null || shift 3
PHASES=${*:-"idem inj"}

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "${HERE}/../.." && pwd)
BUILD="${OUT_DIR}/build_wall"
mkdir -p "${OUT_DIR}"

cleanup() { git -C "${REPO}" checkout -- src/core >/dev/null 2>&1 || true; }
trap cleanup EXIT

cmake -S "${HERE}" -B "${BUILD}" -DCMAKE_BUILD_TYPE=Release > /dev/null

nb_for() {
    case "$1" in
        case9241pegase) echo 50 ;;
        case1354pegase) echo 200 ;;
        *)              echo 2000 ;;
    esac
}

: > "${OUT_DIR}/wall.txt"
for variant in A B; do
    if [ "${variant}" = "B" ]; then python3 "${PATCH}"; fi
    cmake --build "${BUILD}" -j"$(nproc)" > /dev/null
    for grid_path in "${GRIDS_DIR}"/*.lsb; do
        grid=$(basename "${grid_path}" .lsb)
        nb=$(nb_for "${grid}")
        for phase in ${PHASES}; do
            best=""
            for _ in $(seq "${REPEATS}"); do
                ms=$("${BUILD}/profile_cached_pf" "${grid_path}" "${phase}" "${nb}" KLU \
                     | sed -n 's/.*, \([0-9.e+-]*\) ms\/solve.*/\1/p')
                best=$(awk -v a="${best}" -v b="${ms}" \
                       'BEGIN{ if(a=="" || b+0 < a+0) print b; else print a }')
            done
            echo "${grid} ${phase} ${variant} ${best}" | tee -a "${OUT_DIR}/wall.txt"
        done
    done
done

cleanup

echo
echo "============ A/B (wall clock, best of ${REPEATS}) ============"
printf "%-16s %-8s %12s %12s %9s\n" grid phase "A (ms)" "B (ms)" delta
while read -r grid phase _ a; do
    b=$(awk -v g="${grid}" -v p="${phase}" \
        '$1==g && $2==p && $3=="B" {print $4}' "${OUT_DIR}/wall.txt")
    printf "%-16s %-8s %12s %12s %9s\n" "${grid}" "${phase}" "${a}" "${b}" \
        "$(awk -v a="${a}" -v b="${b}" 'BEGIN{printf "%+.2f%%", 100.0*(b-a)/a}')"
done < <(awk '$3=="A"' "${OUT_DIR}/wall.txt")
