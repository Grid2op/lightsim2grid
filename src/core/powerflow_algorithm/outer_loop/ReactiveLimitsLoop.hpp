// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef REACTIVE_LIMITS_LOOP_H
#define REACTIVE_LIMITS_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "batch_algorithm/BusQCheck.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's ReactiveLimits outer loop (ReactiveLimitsOuterLoop, 2.3.0), bus by bus.
 *
 * A controller bus -- a bus whose voltage its units hold -- needs some reactive power from
 * them: the bus' generation. When that leaves the sum of their limits by more than
 * `max_reactive_power_mismatch`, the bus is switched PQ, its generation frozen at the limit
 * (PV -> PQ); a frozen bus whose voltage came back on the side of its set-point the limit
 * was stopping it from reaching is switched PV again (PQ -> PV), at most
 * `max_pq_pv_switch` times; and one whose limit moved since is frozen at the new one. If
 * every PV bus would switch, the strongest one stays PV (highest nominal voltage, then
 * highest target P).
 *
 * Fixed pattern. A controller bus holding its own voltage (the ordinary PV path) is
 * declared switchable (its Vm unknown and Q row reserved, pinned while it is PV); making
 * it PQ is a value edit of the algorithm's injection (OuterState::Sbus) and of the pinned
 * set (OuterState::pq_buses); making it PV again resets its magnitude to its set-point
 * (OuterState::vm_set). The bus' units share what is frozen as the results split any bus'
 * reactive power. A bus whose units are controllers of a voltage-control group (remote
 * regulation, an SVC) is frozen by holding each of them at its own limit
 * (OuterState::controller_hold_q, VoltageControl::set_held_controllers), the group's other
 * controllers regulating on; released, they regulate again.
 *
 * Trigger (detect): bus_q_check (BusQCheck.hpp) on the PV buses, compared with the
 * loop's tolerance, HIGH_Q / LOW_Q; LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q on the
 * frozen ones it would release.
 *
 * An SVC's limits are its susceptance range at the bus' voltage, so they move with the
 * solve. With `robust_mode`, a remote controller bus whose own voltage is unrealistic is
 * frozen at its units' target Q (MIN_REALISTIC_V / MAX_REALISTIC_V in OpenLoadFlow) and
 * restarted from 1 pu, as is one frozen at a limit with such a voltage. Not yet:
 * capability curves at a moving target P.
 */
class LS2G_API ReactiveLimitsLoop final : public BaseOuterLoop
{
    public:
        /// OpenLoadFlow's reactiveLimitsMaxPqPvSwitch
        int max_pq_pv_switch = 3;
        /// OpenLoadFlow's maxReactivePowerMismatch (its newtonRaphsonConvEpsPerEq with the
        /// default stopping criterion), pu of OpenLoadFlow's own 100 MVA base
        real_type max_reactive_power_mismatch = 1e-4;
        static constexpr real_type OLF_SB_MVA = 100.;
        /// OpenLoadFlow's voltageRemoteControlRobustMode: a remote controller whose own bus'
        /// voltage is unrealistic stops regulating, at its units' target Q, its bus back at 1 pu
        bool robust_mode = true;
        /// OpenLoadFlow's minRealisticVoltage / maxRealisticVoltage, pu (the robust mode
        /// compares with them, with a margin)
        real_type min_realistic_voltage = 0.8;
        real_type max_realistic_voltage = 1.2;
        /// OpenLoadFlow's REALISTIC_VOLTAGE_MARGIN
        static constexpr real_type REALISTIC_VOLTAGE_MARGIN = 1.02;

    protected:
        std::string _name() const override { return "ReactiveLimits"; }
        void _declare(const OuterContext & ctx, OuterDeclaration & decl) const override;
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        bool _can_fix_unrealistic_state() const override { return true; }
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new ReactiveLimitsLoop(*this));
        }
        AlgoConfig _get_params() const override;
        void _set_params(const AlgoConfig & params) override;

    private:
        /// A controller bus, and where the loop left it.
        struct ControllerBus {
            int entry = -1;          ///< its entry in plan_
            int bus_solver = -1;
            bool local = true;       ///< held through the PV path (else: by a group)
            real_type min_q = 0.;    ///< the sum of its units' fixed limits, MVar (an SVC's are not)
            real_type max_q = 0.;
            int reg_bus_solver = -1; ///< the bus it regulates
            real_type target_vm = 1.;    ///< its set-point, pu
            real_type nominal_v = 0.;    ///< kV, for the strongest PV bus
            real_type target_p = 0.;     ///< MW, its units' target P, for the same
            int state = 0;           ///< 0 PV, -1 frozen at its min, +1 at its max
            bool realistic = false;  ///< frozen by the robust mode (at its target Q), not a limit
            real_type target_q = 0.; ///< MVar, its units' target Q (the robust mode's)
            real_type frozen_q = 0.;     ///< MVar, while frozen
            int nb_pv_pq = 0;        ///< how many times it was switched PV -> PQ
        };

        /// the plan of every controller bus of `ctx` (bus_q_check)
        static bus_q_check::BusQPlan _plan(const OuterContext & ctx);
        /// the generation of a PV controller bus in this solve, MVar
        real_type _bus_q(const OuterContext & ctx, const ControllerBus & bus) const;
        /// its limits in this solve (an SVC's at the bus' voltage), MVar
        void _limits(const OuterContext & ctx, const ControllerBus & bus, real_type & q_min, real_type & q_max) const;
        /// one controller's limit in this solve, MVar
        real_type _controller_limit(const OuterContext & ctx, int ctrl_pos, bool max) const;
        /// one switch the current solve calls for: bus `k` of buses_, and why
        struct Switch {
            int k;
            LimitViolationType type;
            real_type value;     ///< its generation (MVar) / the voltage it regulates (kV)
            real_type limit;     ///< the limit it left (MVar) / its set-point (kV)
            bool realistic = false;  ///< the robust mode's (frozen at its target Q)
        };
        /// the switches the current solve calls for, as OpenLoadFlow's check computes them
        void _evaluate(const OuterContext & ctx, std::vector<Switch> & to_pq,
                       std::vector<Switch> & to_pv, std::vector<int> & moved,
                       int & remaining_pv) const;
        /// freeze / release a local bus in the algorithm's state
        void _freeze(OuterContext & ctx, ControllerBus & bus, real_type q_mvar, int state) const;
        void _release(OuterContext & ctx, ControllerBus & bus) const;

        // per solve (see _initialize)
        bus_q_check::BusQPlan plan_;
        std::vector<ControllerBus> buses_;
};

}  // namespace ls2g

#endif  // REACTIVE_LIMITS_LOOP_H
