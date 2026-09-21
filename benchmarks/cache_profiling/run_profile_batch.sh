#!/usr/bin/env bash
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
#
# Instruction-count audit of the batch algorithms (TimeSeries, ContingencyAnalysis)
# -- the sibling of run_profile.sh, driving profile_batch instead.
#
#   ./run_profile_batch.sh <grids_dir> <build_dir> <out_dir> [algo]
#
# Produces, per (grid, phase): the raw callgrind output, a `callgrind_annotate`
# dump (exclusive and inclusive), and the number of ROWS the collected region
# covers, so summarize.py reports instructions per row. Every compute phase is
# run twice, over N rows and over 1 row (tag `<phase>_1row`): the difference,
# divided by N-1, is the marginal cost of a row, and the 1-row figure alone is
# what a compute() pays before its first row (the solver input rebuild and the
# "n" solve) -- two numbers a single per-row average would blur together.

set -euo pipefail

GRIDS_DIR=${1:-grids}
BUILD_DIR=${2:-../../build_profile}
OUT_DIR=${3:-callgrind_out_batch}
ALGO=${4:-KLU}

BIN="${BUILD_DIR}/profile_batch"
mkdir -p "${OUT_DIR}"

# the compute phases (run over N rows and over 1 row), then the ones that only
# make sense over N rows: reading flows back, and the per-step construction
COMPUTE_PHASES="ts_ac ts_dc ca_ac ca_dc ca_ac_mask ca_dc_mask"
OTHER_PHASES="ts_flows ca_flows ca_construct"

# same fixed, locale-independent order as run_profile.sh
grids_in_order() {
    for f in "$1"/*.lsb; do
        local name fancy=0
        name=$(basename "${f}" .lsb)
        case "${name}" in *_fancy*) fancy=1 ;; esac
        printf '%s %s %s\n' "${fancy}" "$(echo "${name}" | sed -E 's/^[a-z]*([0-9]+).*/\1/')" "${f}"
    done | sort -k1,1n -k2,2n | awk '{print $3}'
}

run_one() {  # grid_path grid phase tag nb
    local grid_path=$1 grid=$2 phase=$3 tag=$4 nb=$5
    local out="${OUT_DIR}/callgrind.${grid}.${tag}.out"

    echo "=== ${grid} / ${tag} (${nb} rows, ${ALGO}) ==="
    valgrind --tool=callgrind \
             --instr-atstart=no --collect-atstart=no \
             --cache-sim=no --branch-sim=no \
             --callgrind-out-file="${out}" \
             "${BIN}" "${grid_path}" "${phase}" "${nb}" "${ALGO}" always \
             > "${OUT_DIR}/run.${grid}.${tag}.log" \
             2> "${OUT_DIR}/valgrind.${grid}.${tag}.log" \
        || { echo "    driver failed (rc=$?): see valgrind.${grid}.${tag}.log -- phase skipped"; return 0; }
    echo "${nb}" > "${OUT_DIR}/nb.${grid}.${tag}.txt"
    callgrind_annotate --auto=no --threshold=99 "${out}" \
        > "${OUT_DIR}/annotate.${grid}.${tag}.txt"
    callgrind_annotate --inclusive=yes --threshold=99.9 "${out}" \
        > "${OUT_DIR}/inclusive.${grid}.${tag}.txt"
}

for grid_path in $(grids_in_order "${GRIDS_DIR}"); do
    grid=$(basename "${grid_path}" .lsb)
    # a row is a whole solve and the big grids are ~50x slower under callgrind
    case "${grid}" in
        case9241pegase*) nb=20 ;;
        case1354pegase*) nb=50 ;;
        *)               nb=200 ;;
    esac
    for phase in ${COMPUTE_PHASES}; do
        run_one "${grid_path}" "${grid}" "${phase}" "${phase}" "${nb}"
        run_one "${grid_path}" "${grid}" "${phase}" "${phase}_1row" 1
    done
    for phase in ${OTHER_PHASES}; do
        this_nb=${nb}
        # a construction is a copy of the grid plus 8 solves: a few are enough
        if [ "${phase}" = "ca_construct" ]; then this_nb=3; fi
        run_one "${grid_path}" "${grid}" "${phase}" "${phase}" "${this_nb}"
    done
done

echo
echo "Now run:  python3 summarize.py ${OUT_DIR}"
