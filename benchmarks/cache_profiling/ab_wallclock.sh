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
#   DRIVER=batch ./ab_wallclock.sh <grids_dir> <out_dir> <patch_script> [repeats] [phases...]
#
# Each (grid, phase, variant) is run `repeats` times and the BEST run is kept --
# the minimum is the least noisy statistic here, a slower run only ever means
# more interference. DRIVER=batch measures profile_batch's phases (per row)
# instead of profile_cached_pf's (per solve) -- see ab_test.sh. The batch driver
# also prints its peak resident set size, kept alongside the timing: the memory
# items of the batch audit are invisible to an instruction count.

set -euo pipefail

GRIDS_DIR=$1
OUT_DIR=$2
PATCH=$3
REPEATS=${4:-7}
shift 4 2>/dev/null || shift 3
DRIVER=${DRIVER:-cached_pf}
if [ "${DRIVER}" = "batch" ]; then
    DEFAULT_PHASES="ts_ac ts_dc ca_ac ca_dc ts_flows ca_flows"
else
    DEFAULT_PHASES="idem inj"
fi
PHASES=${*:-${DEFAULT_PHASES}}

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "${HERE}/../.." && pwd)
BUILD="${OUT_DIR}/build_wall"
mkdir -p "${OUT_DIR}"

# The grids in a fixed, locale-independent order: every plain case first, then the
# `_fancy` ones, each family by increasing size (the number in the case name). A
# plain `*.lsb` glob sorts by locale, and a French one puts `case9241pegase_fancy`
# before `case9241pegase` -- which, with a phase that fails on the fancy grid, is
# how the plain one got skipped.
grids_in_order() {
    for f in "$1"/*.lsb; do
        local name fancy=0
        name=$(basename "${f}" .lsb)
        case "${name}" in *_fancy*) fancy=1 ;; esac
        printf '%s %s %s\n' "${fancy}" "$(echo "${name}" | sed -E 's/^[a-z]*([0-9]+).*/\1/')" "${f}"
    done | sort -k1,1n -k2,2n | awk '{print $3}'
}

cleanup() { git -C "${REPO}" checkout -- src/core >/dev/null 2>&1 || true; }
trap cleanup EXIT

cmake -S "${HERE}" -B "${BUILD}" -DCMAKE_BUILD_TYPE=Release > /dev/null

nb_for() {
    if [ "${DRIVER}" = "batch" ]; then
        # rows, each a whole solve: enough for the result matrices to leave the
        # cache on the big grids, which is the point of the flow phases
        case "$1" in
            case9241pegase*) echo 100 ;;
            case1354pegase*) echo 500 ;;
            *)               echo 2000 ;;
        esac
        return
    fi
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
    for grid_path in $(grids_in_order "${GRIDS_DIR}"); do
        grid=$(basename "${grid_path}" .lsb)
        nb=$(nb_for "${grid}")
        for phase in ${PHASES}; do
            best=""
            rss=""
            for _ in $(seq "${REPEATS}"); do
                line=$("${BUILD}/profile_${DRIVER}" "${grid_path}" "${phase}" "${nb}" KLU)
                ms=$(echo "${line}" | sed -n 's/.*, \([0-9.e+-]*\) ms\/solve.*/\1/p')
                rss=$(echo "${line}" | sed -n 's/.*peak rss \([0-9]*\) kB.*/\1/p')
                best=$(awk -v a="${best}" -v b="${ms}" \
                       'BEGIN{ if(a=="" || b+0 < a+0) print b; else print a }')
            done
            echo "${grid} ${phase} ${variant} ${best} ${rss:-0}" | tee -a "${OUT_DIR}/wall.txt"
        done
    done
done

cleanup

echo
echo "============ A/B (wall clock, best of ${REPEATS}) ============"
printf "%-16s %-11s %12s %12s %9s %12s %12s\n" grid phase "A (ms)" "B (ms)" delta "A rss (kB)" "B rss (kB)"
while read -r grid phase _ a rss_a; do
    b=$(awk -v g="${grid}" -v p="${phase}" \
        '$1==g && $2==p && $3=="B" {print $4}' "${OUT_DIR}/wall.txt")
    rss_b=$(awk -v g="${grid}" -v p="${phase}" \
        '$1==g && $2==p && $3=="B" {print $5}' "${OUT_DIR}/wall.txt")
    printf "%-16s %-11s %12s %12s %9s %12s %12s\n" "${grid}" "${phase}" "${a}" "${b}" \
        "$(awk -v a="${a}" -v b="${b}" 'BEGIN{printf "%+.2f%%", 100.0*(b-a)/a}')" \
        "${rss_a:-0}" "${rss_b:-0}"
done < <(awk '$3=="A"' "${OUT_DIR}/wall.txt")
