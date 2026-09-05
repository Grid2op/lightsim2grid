#!/usr/bin/env bash
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
#
# Instruction-count audit of the cached powerflow path.
#
#   ./run_profile.sh <grids_dir> <build_dir> <out_dir> [algo]
#
# Produces, per (grid, phase): the raw callgrind output, a `callgrind_annotate`
# dump (exclusive and inclusive), and the number of solves the collected region
# covers. Numbers are Ir (instructions retired) for the collected region only --
# the powerflow calls themselves, never the grid loading, the warm-up solves, or
# the mutation applied between two steps. Feed the directory to summarize.py.

set -euo pipefail

GRIDS_DIR=${1:-grids}
BUILD_DIR=${2:-../../build_profile}
OUT_DIR=${3:-callgrind_out}
ALGO=${4:-KLU}

BIN="${BUILD_DIR}/profile_cached_pf"
mkdir -p "${OUT_DIR}"

# One entry per callgrind run, written as `phase[:refactor[:algo]]`:
#   * `cold` is a single solve by construction; the others average over several
#     so the per-solve figure is not one outlier;
#   * the `:every2` / `:every3` entries re-run the ordinary `inj` phase under a
#     refactorization policy that reuses the factorization across Newton
#     iterations -- what puts a number on refactorizing every iteration;
#   * the `::NRSing_KLU` entries re-run it with the single-slack Newton system,
#     which is the same solve without the distributed-slack extension.
PHASES="cold idem inj inj_nores dcac nocache topo
        inj:every2 inj:every3
        idem::NRSing_KLU inj::NRSing_KLU"

run_one() {
    local grid_path=$1 grid=$2 spec=$3 nb=$4
    local phase refactor algo tag
    phase=${spec%%:*}
    local rest=${spec#*:}
    if [ "${rest}" = "${spec}" ]; then rest=""; fi
    refactor=${rest%%:*}
    local algo_part=${rest#*:}
    if [ "${algo_part}" = "${rest}" ]; then algo_part=""; fi
    [ -z "${refactor}" ] && refactor=always
    algo=${algo_part:-${ALGO}}
    tag=${spec//:/_}
    tag=${tag%_}
    local out="${OUT_DIR}/callgrind.${grid}.${tag}.out"

    echo "=== ${grid} / ${tag} (${nb} solves, ${algo}, refactor=${refactor}) ==="
    valgrind --tool=callgrind \
             --instr-atstart=no --collect-atstart=no \
             --cache-sim=no --branch-sim=no \
             --callgrind-out-file="${out}" \
             "${BIN}" "${grid_path}" "${phase}" "${nb}" "${algo}" "${refactor}" \
             2> "${OUT_DIR}/valgrind.${grid}.${tag}.log"
    echo "${nb}" > "${OUT_DIR}/nb.${grid}.${tag}.txt"
    callgrind_annotate --auto=no --threshold=99 "${out}" \
        > "${OUT_DIR}/annotate.${grid}.${tag}.txt"
    callgrind_annotate --inclusive=yes --threshold=99.9 "${out}" \
        > "${OUT_DIR}/inclusive.${grid}.${tag}.txt"
}

for grid_path in "${GRIDS_DIR}"/*.lsb; do
    grid=$(basename "${grid_path}" .lsb)
    # the big grids are ~50x slower under callgrind: fewer repeats there
    case "${grid}" in
        case9241pegase) nb=5 ;;
        case1354pegase) nb=10 ;;
        *)              nb=20 ;;
    esac
    for spec in ${PHASES}; do
        # `cold` measures the very first solve, so it is one solve by definition
        this_nb=${nb}
        if [ "${spec}" = "cold" ]; then this_nb=1; fi
        run_one "${grid_path}" "${grid}" "${spec}" "${this_nb}"
    done
done

echo
echo "Now run:  python3 summarize.py ${OUT_DIR}"
