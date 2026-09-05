#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Candidate change under A/B test (ab_test.sh): stop paying for `std::vector<bool>`
eleven times per branch in the branch-flow loop.

`status1[el_id]` and `status2[el_id]` are read five times each per branch, and
each read is a bit-vector access -- a word offset, a shift and a mask, not a
load. On top of that, `get_bus_side_1_internal(el_id)` re-reads the SAME status
bit to decide whether to hand back `_deactivated_bus_id`, inside a branch that
has just established the side is connected, and it does not inline (it shows up
in the profile as `OneSideContainer::get_bus_internal`, 0.51M Ir per solve on
case9241pegase).

Reading each status bit once into a local, and taking the bus id straight from
the side's own `bus_id_` vector where the status is already known, is pure
bookkeeping: same values, same order, same arithmetic.
"""

import os
import sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "..", "src", "core", "element_container",
                    "TwoSidesContainer_rxh_A.hpp")

# (anchor, replacement) applied in order; each must match exactly once
EDITS = [
    (
        """            const std::vector<bool> & status1 = side_1_.get_status();
            const std::vector<bool> & status2 = side_2_.get_status();
            for(int el_id = 0; el_id < nb_element; ++el_id){

                // don't do anything if the element is disconnected
                if(!status_global_[el_id] || (!status1[el_id] && !status2[el_id])) {""",
        """            const std::vector<bool> & status1 = side_1_.get_status();
            const std::vector<bool> & status2 = side_2_.get_status();
            // Read once per branch, not once per use. These are std::vector<bool>,
            // so every `status1[el_id]` is a word offset, a shift and a mask rather
            // than a load, and the loop below used to ask each of them five times.
            const GlobalBusIdVect & buses1 = side_1_.get_buses();
            const GlobalBusIdVect & buses2 = side_2_.get_buses();
            for(int el_id = 0; el_id < nb_element; ++el_id){
                const bool st1 = status1[el_id];
                const bool st2 = status2[el_id];

                // don't do anything if the element is disconnected
                if(!status_global_[el_id] || (!st1 && !st2)) {""",
    ),
    (
        """                if(status1[el_id]){
                    bus_hv_id_me = get_bus_side_1_internal(el_id);""",
        """                if(st1){
                    // the side is connected: read the bus id itself rather than
                    // going through get_bus_side_1_internal, which re-reads the
                    // status bit we just tested (and does not inline)
                    bus_hv_id_me = GridModelBusId(buses1(el_id).cast_int());""",
    ),
    (
        """                if(status2[el_id]){
                    bus_lv_id_me = get_bus_side_2_internal(el_id);""",
        """                if(st2){
                    bus_lv_id_me = GridModelBusId(buses2(el_id).cast_int());""",
    ),
    (
        """                // retrieve voltages magnitude in kv instead of pu
                if(status1[el_id]){
                    real_type v_hv = Vm(bus_hv_solver_id.cast_int());""",
        """                // retrieve voltages magnitude in kv instead of pu
                if(st1){
                    real_type v_hv = Vm(bus_hv_solver_id.cast_int());""",
    ),
    (
        """                if(status2[el_id]){
                    real_type v_lv = Vm(bus_lv_solver_id.cast_int());""",
        """                if(st2){
                    real_type v_lv = Vm(bus_lv_solver_id.cast_int());""",
    ),
    (
        """                    const cplx_type Ehv = status1[el_id] ? V(bus_hv_solver_id.cast_int()) : cplx_type(0., 0.);
                    const cplx_type Elv = status2[el_id] ? V(bus_lv_solver_id.cast_int()) : cplx_type(0., 0.);""",
        """                    const cplx_type Ehv = st1 ? V(bus_hv_solver_id.cast_int()) : cplx_type(0., 0.);
                    const cplx_type Elv = st2 ? V(bus_lv_solver_id.cast_int()) : cplx_type(0., 0.);""",
    ),
    (
        """                    if(status1[el_id] && status2[el_id]){
                        real_type Va_hv = Va(bus_hv_solver_id.cast_int());""",
        """                    if(st1 && st2){
                        real_type Va_hv = Va(bus_hv_solver_id.cast_int());""",
    ),
]


def main():
    with open(PATH) as fp:
        src = fp.read()
    if "const bool st1 = status1[el_id];" in src:
        print("branch_results_hoist: already applied")
        return 0
    for old, new in EDITS:
        if src.count(old) != 1:
            print(f"branch_results_hoist: anchor found {src.count(old)} times, expected 1:\n{old[:120]}",
                  file=sys.stderr)
            return 1
        src = src.replace(old, new)
    with open(PATH, "w") as fp:
        fp.write(src)
    print("branch_results_hoist: applied")
    return 0


if __name__ == "__main__":
    sys.exit(main())
