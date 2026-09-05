// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASECONSTANTS_H
#define BASECONSTANTS_H

#include "ls2g_api.hpp"
#include "Utils.hpp"
#include "Eigen/Core"

namespace ls2g {

/**
Definition of some basic constants, because sometimes Eigen cannot deduce types.
Eg if I type "1.0" then Eigen cast it to "double" and i cannot use it with real_type = float for example
**/
class LS2G_API BaseConstants
{
    public:
        static const cplx_type my_i;
        static const real_type my_pi;
        static const real_type my_half_;
        static const real_type my_one_;
        static const real_type my_two_;
        static const real_type my_zero_;

        // ---- (magnitude, angle) canonical form, for every algorithm ----------
        //
        // Put a (Vm, Va) pair back in canonical form: magnitude >= 0, angle in
        // [-pi, pi]. Both can leave that form mid-solve, in either family: the
        // FDPF Q iteration (`Vm_(pq) -= q_`) and the Newton-Raphson step
        // (`Vm_(bus) += dx(col)`) can each drive a magnitude past zero, and both
        // accumulate into Va_ without ever wrapping it.
        //
        // One implementation, because the two families had grown two: the FDPF
        // had the cheap form below, while NRSystem::apply_step still repaired the
        // same overshoot with `Vm_ = V_.abs(); Va_ = V_.arg();` -- the hypot and
        // atan2 pair the cheap form exists to replace -- and, because atan2
        // returns in (-pi, pi], wrapped as a side effect where the FDPF wrapped
        // deliberately. Same intent, two implementations, two different answers.
        //
        // has_converged used to restore it by going through the complex
        // voltage -- `Vm_ = V_.abs(); Va_ = V_.arg();` -- which costs a hypot
        // and an atan2 per bus, twice per iteration, to rediscover two numbers
        // the solver already holds: V_ is built as Vm_ * exp(i.Va_), so |V_| is
        // just |Vm_| and arg(V_) is just Va_, plus half a turn where the
        // magnitude went negative. On a 121-bus grid that pair cost 2.5 us per
        // call out of a ~150 us solve.
        //
        // Equivalent to the abs()/arg() pair over the whole domain (checked in
        // test_fdpf_algorithm.cpp over Vm in [-3, 3] x Va in [-10, 10] rad),
        // with ONE exception: at Vm == 0 exactly, where the phase of a zero
        // phasor is undefined. arg() returned whichever of 0 / +pi / -pi
        // atan2's signed-zero rules produced for (+-0, +-0); this keeps the
        // incoming angle. Neither can reach a converged answer, because the
        // caller's next step divides the mismatch by that zero and fails its
        // allFinite() check.
        //
        // Kept as a named static (rather than inline in has_converged) so the
        // identity it rests on can be tested directly: the branches below are
        // unreachable from a converging solve -- a trajectory that overshoots a
        // magnitude ends up diverging, and a diverged solve clears its state.
        // The half that has to run on every iteration: a negative magnitude turned
        // into a positive one plus a half turn. `has_converged` divides the
        // mismatch by Vm_ immediately after (`mis_ /= Vm_`), so a negative
        // magnitude there does not merely look wrong -- it flips the sign of that
        // bus's P and Q and corrupts the convergence test.
        //
        // The `is there anything to do` test lives here, not at the call sites:
        // both of them (BaseFDPFAlgo::has_converged and NRSystem::apply_step)
        // want exactly this one, so duplicating it only gave the two families a
        // second chance to drift apart. It is what makes the function cheap --
        // a converging trajectory never drives a magnitude negative, so the
        // ordinary solve pays one vectorised read-only reduction over Vm instead
        // of the two full-length passes below. That is worth 5.6% of an FDPF
        // solve on case9241pegase (674.1M -> 636.3M Ir); an unguarded version of
        // this same function measured level with the code it replaced.
        static void fix_negative_vm(Eigen::Ref<RealVect> Vm, Eigen::Ref<RealVect> Va)
        {
            if (Vm.minCoeff() >= my_zero_) return;  // the ordinary case
            // a negative magnitude is the same phasor turned by pi
            Va.array() = (Vm.array() < my_zero_).select(Va.array() + my_pi, Va.array());
            Vm = Vm.cwiseAbs();
        }

        // The half that does not: bringing the angle back into [-pi, pi]. Nothing
        // inside the solve reads Va_ except the P iteration's own accumulation
        // (`Va_(pvpq) -= p_`) and the cos / sin that rebuild V, and both are
        // indifferent to a multiple of 2.pi -- V comes out bit for bit the same.
        // So this is about the pair the solver REPORTS, and once per solve is
        // enough. (The Newton-Raphson does not wrap at all: see
        // NRSystem::apply_step, whose fix-up covers only the negative magnitude.)
        //
        // Is any angle outside [-pi, pi]? One pass and a reduction, against the
        // four passes the wrap itself costs -- worth asking only if the answer is
        // usually no, which is what the benchmark says it is.
        static bool va_out_of_range(const Eigen::Ref<const RealVect> & Va)
        {
            return Va.cwiseAbs().maxCoeff() > my_pi;
        }

        // Guarded for the same reason fix_negative_vm is, and in the same place:
        // both call sites (the end of BaseFDPFAlgo::compute_pf and of
        // NRAlgo::compute_pf) asked the same question, so it is asked here once.
        // Without the guard this subtracts an exact zero from every angle, which
        // is correct and which the unguarded variant showed costs a little more
        // than the reduction that skips it.
        static void wrap_va(Eigen::Ref<RealVect> Va)
        {
            if (!va_out_of_range(Va)) return;  // the ordinary case
            Va.array() -= my_two_ * my_pi * (Va.array() * (my_one_ / (my_two_ * my_pi))).round();
        }

        // Both halves, in order. Kept as one entry point because the identity it
        // rests on -- that this is the abs()/arg() pair -- is what
        // test_fdpf_algorithm.cpp checks head-on, over the whole domain.
        static void canonicalise_vm_va(Eigen::Ref<RealVect> Vm, Eigen::Ref<RealVect> Va)
        {
            fix_negative_vm(Vm, Va);
            wrap_va(Va);
        }

        static const real_type my_180_pi_;
        static const real_type v_disco_el_;
        static const real_type theta_disco_el_;
        static const int _deactivated_bus_id;
        static const real_type _tol_equal_float;
        static const real_type _1_sqrt_3;
};

enum class LS2G_API FDPFMethod {XB, BX};  // Different type of FDPF powerflow
// FDPFMethod::XB => alg = 2 in pypower / pandapower
// FDPFMethod::BX => alg = 3 in pypower / pandapower


} // namespace ls2g

#endif // BASECONSTANTS_H