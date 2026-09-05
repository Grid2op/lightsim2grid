#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

"""
Candidate change under A/B test (ab_test.sh): the per-bus power mismatch in
`LSGrid::compute_results`, without the heap temporary and without std::complex.

    mismatch = V.array() * (ac_cache_.mat * V).conjugate().array() - ac_cache_.inj.array();
    active_mismatch = mismatch.real() * sn_mva_;
    ...
    reactive_mismatch = mismatch.imag() * sn_mva_;

costs three things per solve that it does not have to. `ac_cache_.mat * V` is a
sparse-times-dense product inside a coefficient-wise expression, which Eigen can
only evaluate into a heap temporary -- NRSystem::_residual_into hit exactly this
and fixed it with a persistent buffer, and says so. The coefficient-wise complex
multiply then goes through `std::complex::operator*`, which follows every product
with a branch that re-derives the result if it came out NaN -- the same cost the
dS pass and the branch-flow loop both already spell out on real and imaginary
parts to avoid. And `mismatch` itself is a full complex vector that nothing reads:
its real and imaginary halves are immediately scaled into two real vectors, and
those are all anyone downstream sees.

One pass over the buses writing the two real vectors directly removes all three.
The arithmetic is grouped exactly as std::complex groups it (real = a.r*b.r -
a.i*b.i, imag = a.r*b.i + a.i*b.r, conj an exact sign flip), so the values are
bit-identical.
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
CPP = os.path.join(ROOT, "src", "core", "LSGrid.cpp")
HPP = os.path.join(ROOT, "src", "core", "LSGrid.hpp")

OLD_CPP = """    //handle_slack_bus active power
    CplxVect mismatch;  // power mismatch at each bus (SOLVER BUS !!!)
    RealVect reactive_mismatch;  // not used in dc mode (DO NOT ATTEMPT TO USE IT THERE)
    RealVect active_mismatch;
    if(ac){
        // In AC mode i am not forced to run through all the grid
        // auto tmp = (ac_cache_.mat * V).conjugate();
        mismatch = V.array() * (ac_cache_.mat * V).conjugate().array() - ac_cache_.inj.array();
        active_mismatch = mismatch.real() * sn_mva_;
    } else{"""

NEW_CPP = """    //handle_slack_bus active power
    RealVect reactive_mismatch;  // not used in dc mode (DO NOT ATTEMPT TO USE IT THERE)
    RealVect active_mismatch;
    if(ac){
        // The per-bus complex power mismatch, V .* conj(Ybus . V) - Sbus, written
        // straight onto the two real vectors that are the only thing anyone reads
        // (set_p_slack takes the active half, set_q the reactive one). It used to
        // be one Eigen expression, which cost three things it did not have to:
        //   - `ac_cache_.mat * V` inside a coefficient-wise expression is a
        //     sparse-times-dense product Eigen can only evaluate into a HEAP
        //     TEMPORARY, one allocation (and one full-vector zero-fill) per solve.
        //     NRSystem::_residual_into hit the same wall and answered it the same
        //     way, with a buffer that outlives the call;
        //   - the coefficient-wise complex multiply went through
        //     std::complex::operator*, which follows every product with a branch
        //     that re-derives the result if it came out NaN -- the cost the dS pass
        //     and the branch-flow loop already spell out on real and imaginary
        //     parts to avoid;
        //   - and `mismatch` was a full complex vector nothing ever looked at.
        // Bit-identical: the products are grouped exactly as std::complex groups
        // them and conj is an exact sign flip, so this is the same arithmetic in
        // the same order.
        ybus_v_res_.noalias() = ac_cache_.mat * V;   // noalias: no product temporary
        const Eigen::Index nb_bus_solver = V.size();
        active_mismatch.resize(nb_bus_solver);
        reactive_mismatch.resize(nb_bus_solver);
        for(Eigen::Index bus_id = 0; bus_id < nb_bus_solver; ++bus_id){
            const real_type v_r = std::real(V(bus_id)), v_i = std::imag(V(bus_id));
            // conj(Ybus . V)
            const real_type yv_r =  std::real(ybus_v_res_(bus_id));
            const real_type yv_i = -std::imag(ybus_v_res_(bus_id));
            active_mismatch(bus_id) =
                (v_r * yv_r - v_i * yv_i - std::real(ac_cache_.inj(bus_id))) * sn_mva_;
            reactive_mismatch(bus_id) =
                (v_r * yv_i + v_i * yv_r - std::imag(ac_cache_.inj(bus_id))) * sn_mva_;
        }
    } else{"""

OLD_TAIL = """    generators_.set_p_slack(active_mismatch, id_me_to_solver);

    if(ac) reactive_mismatch = mismatch.imag() * sn_mva_;
    // mainly to initialize the Q value of the generators in dc (just fill it with 0.)"""

NEW_TAIL = """    generators_.set_p_slack(active_mismatch, id_me_to_solver);

    // (in AC `reactive_mismatch` was filled in the same pass as the active half,
    // above; in DC it stays empty, which is what set_q expects there)
    // mainly to initialize the Q value of the generators in dc (just fill it with 0.)"""

OLD_HPP = """        // 5. generators
        RealVect total_q_min_per_bus_;  // TODO switches: move to BaseSubstation"""

NEW_HPP = """        /**
         * Ybus . V, scratch for the per-bus mismatch compute_results() derives from
         * it. A member rather than a local so the allocation happens once per
         * topology instead of once per solve: Eigen evaluates a sparse-times-dense
         * product into a heap temporary whenever it appears inside a larger
         * expression, which is what this used to be.
         *
         * AC only, and meaningless outside the compute_results() call that fills it.
         */
        CplxVect ybus_v_res_;

        // 5. generators
        RealVect total_q_min_per_bus_;  // TODO switches: move to BaseSubstation"""


def patch(path, old, new, name):
    with open(path) as fp:
        src = fp.read()
    if new in src:
        print(f"bus_mismatch_no_temporary: {name} already applied")
        return src, False
    if src.count(old) != 1:
        print(f"bus_mismatch_no_temporary: anchor '{name}' found {src.count(old)} times, "
              f"expected 1", file=sys.stderr)
        sys.exit(1)
    return src.replace(old, new), True


def main():
    src, changed_h = patch(HPP, OLD_HPP, NEW_HPP, "LSGrid.hpp member")
    if changed_h:
        with open(HPP, "w") as fp:
            fp.write(src)
    src, changed_1 = patch(CPP, OLD_CPP, NEW_CPP, "mismatch pass")
    if changed_1:
        with open(CPP, "w") as fp:
            fp.write(src)
    src, changed_2 = patch(CPP, OLD_TAIL, NEW_TAIL, "reactive tail")
    if changed_2:
        with open(CPP, "w") as fp:
            fp.write(src)
    print("bus_mismatch_no_temporary: applied")
    return 0


if __name__ == "__main__":
    sys.exit(main())
