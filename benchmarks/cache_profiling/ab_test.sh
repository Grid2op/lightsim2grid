#!/usr/bin/env bash
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
#
# A/B one candidate change to src/core: build A (the tree as it is), build B (the
# tree with a patch applied), then run every (grid, phase) under BOTH and report
#   * the instructions retired per solve, A vs B;
#   * whether the two builds return the SAME answer -- every solve's iteration
#     count and full complex voltage vector, compared with 17 significant digits.
#
#   ./ab_test.sh <grids_dir> <out_dir> <patch_script> [phases...]
#
# <patch_script> is a python script that edits src/core in place; the tree is
# restored with `git checkout` afterwards, and again on exit.

set -euo pipefail

GRIDS_DIR=$1
OUT_DIR=$2
PATCH=$3
shift 3
PHASES=${*:-"idem inj dcac topo nocache cold"}

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "${HERE}/../.." && pwd)
BUILD="${OUT_DIR}/build"
mkdir -p "${OUT_DIR}"

cleanup() { git -C "${REPO}" checkout -- src/core >/dev/null 2>&1 || true; }
trap cleanup EXIT

cmake -S "${HERE}" -B "${BUILD}" -DCMAKE_BUILD_TYPE=Release > /dev/null

nb_for() {
    case "$1" in
        case9241pegase) echo 5 ;;
        case1354pegase) echo 10 ;;
        *)              echo 20 ;;
    esac
}

measure() {  # variant grid phase nb
    local variant=$1 grid=$2 phase=$3 nb=$4
    local out="${OUT_DIR}/cg.${variant}.${grid}.${phase}.out"
    local trace="${OUT_DIR}/trace.${variant}.${grid}.${phase}.txt"
    valgrind --tool=callgrind --instr-atstart=no --collect-atstart=no \
             --cache-sim=no --branch-sim=no --callgrind-out-file="${out}" \
             "${BUILD}/profile_cached_pf" "${GRIDS_DIR}/${grid}.lsb" "${phase}" \
             "${nb}" KLU always "${trace}" \
             > "${OUT_DIR}/run.${variant}.${grid}.${phase}.log" 2>&1
    callgrind_annotate --threshold=1 "${out}" 2>/dev/null \
        | sed -n 's/^ *\([0-9,]*\).*PROGRAM TOTALS.*/\1/p' | tr -d ,
}

for variant in A B; do
    if [ "${variant}" = "B" ]; then
        echo "--- applying ${PATCH} ---"
        python3 "${PATCH}"
    fi
    cmake --build "${BUILD}" -j"$(nproc)" > /dev/null
    for grid_path in "${GRIDS_DIR}"/*.lsb; do
        grid=$(basename "${grid_path}" .lsb)
        nb=$(nb_for "${grid}")
        for phase in ${PHASES}; do
            [ "${phase}" = "cold" ] && nb=1
            total=$(measure "${variant}" "${grid}" "${phase}" "${nb}")
            echo "${grid} ${phase} ${variant} $((total / nb))" \
                | tee -a "${OUT_DIR}/totals.txt"
            nb=$(nb_for "${grid}")
        done
    done
done

cleanup

echo
echo "================ A/B ================"
printf "%-16s %-10s %14s %14s %9s   %s\n" grid phase A B delta answer
while read -r grid phase _ a; do
    b=$(awk -v g="${grid}" -v p="${phase}" \
        '$1==g && $2==p && $3=="B" {print $4}' "${OUT_DIR}/totals.txt")
    delta=$(awk -v a="${a}" -v b="${b}" 'BEGIN{printf "%+.2f%%", 100.0*(b-a)/a}')
    if cmp -s "${OUT_DIR}/trace.A.${grid}.${phase}.txt" \
              "${OUT_DIR}/trace.B.${grid}.${phase}.txt"; then
        same="identical"
    else
        same="DIFFERS"
    fi
    printf "%-16s %-10s %14s %14s %9s   %s\n" "${grid}" "${phase}" "${a}" "${b}" "${delta}" "${same}"
done < <(awk '$3=="A"' "${OUT_DIR}/totals.txt")
