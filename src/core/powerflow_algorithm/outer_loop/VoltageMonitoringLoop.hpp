// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef VOLTAGE_MONITORING_LOOP_H
#define VOLTAGE_MONITORING_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "batch_algorithm/SvcStandbyCheck.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's VoltageMonitoring outer loop (MonitoringVoltageOuterLoop, 2.3.0).
 *
 * An SVC whose standby automaton is in standby is a voltage monitor
 * (init_from_pypowsybl with svc_voltage_monitoring): idle, at Q = 0, while the voltage of the
 * bus it regulates stays inside its [low, high] thresholds. When it leaves them, the SVC is
 * switched on, once and for good: it regulates that bus at the automaton's low set-point
 * (below the low threshold) or its high one (above the high threshold).
 *
 * Such an SVC is idle in the grid (off, flagged standby with its set-points, see
 * LSGrid::set_svc_standby). With this loop in the grid's list, the voltage-control plan
 * holds it in a group of its own (VoltageControlPlan::build_controllers' `hold_monitors`),
 * so switching it on is a value edit (OuterState::svc_target_vm, then
 * NRSystem::release_held_svcs): one symbolic analysis for the whole solve.
 *
 * Trigger (detect): svc_standby_check (SvcStandbyCheck.hpp) on the monitors still idle,
 * compared strictly, LOW_VOLTAGE_SVC_STANDBY / HIGH_VOLTAGE_SVC_STANDBY. In detection mode,
 * on every idle SVC the grid flags standby.
 *
 * Needed (isNeeded) only with a monitor regulating its own bus, as OpenLoadFlow (which looks
 * at the controller buses that are also controlled): a monitor regulating a remote bus is
 * switched on only when a local one makes the loop run.
 */
class LS2G_API VoltageMonitoringLoop final : public BaseOuterLoop
{
    protected:
        std::string _name() const override { return "VoltageMonitoring"; }
        bool _holds_svc_monitors() const override { return true; }
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new VoltageMonitoringLoop(*this));
        }

    private:
        /// the SVCs the voltage-control plan of `grid` holds (its monitors), by SVC id
        static std::vector<int> _held_svcs(const LSGrid & grid);

        // per solve (see _initialize): the monitors, as svc_standby_check entries
        svc_standby_check::SvcStandbyPlan monitors_;
};

}  // namespace ls2g

#endif  // VOLTAGE_MONITORING_LOOP_H
