// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef HVDC_AC_EMULATION_LIMITS_LOOP_H
#define HVDC_AC_EMULATION_LIMITS_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "batch_algorithm/HvdcPCheck.hpp"

namespace ls2g {

/**
 * OpenLoadFlow's AcHvdcAcEmulationLimits outer loop (AcHvdcAcEmulationLimitsOuterLoop, 2.3.0).
 *
 * An angle-droop ("AC emulation") hvdc line in its linear regime whose droop flow leaves
 * what its converters can transmit in that direction is saturated there: the sending side
 * injects its maximum, the receiving side that minus the losses
 * (HvdcDroopSolverData::flows_pu). A saturated line whose droop flow comes back strictly
 * inside the limit of the direction it now flows in goes back to its linear regime; one whose
 * flow reversed beyond the other direction's limit is saturated on that side instead.
 *
 * The regime is the algorithm's (OuterState::hvdc_status), handed to the Newton's Hvdc
 * extension by value -- its entries are declared in every regime -- and used to publish the
 * line's flows; the line's own status_droop is not modified.
 *
 * Trigger (detect): hvdc_p_check (HvdcPCheck.hpp) on the loop's regimes, HIGH_P to saturate
 * and HVDC_AC_EMULATION_RELEASE to release. In detection mode, on the grid's own: its
 * linear lines and the ones the caller says an outer loop froze at a limit.
 */
class LS2G_API HvdcAcEmulationLimitsLoop final : public BaseOuterLoop
{
    protected:
        std::string _name() const override { return "AcHvdcAcEmulationLimits"; }
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new HvdcAcEmulationLimitsLoop(*this));
        }

    private:
        /// the lines in AC emulation (an active droop, in its linear regime) of `grid`
        static hvdc_p_check::HvdcPPlan _ac_emulation_lines(const OuterContext & ctx);

        // per solve (see _initialize): the lines, each with the loop's regime in frozen_dir
        hvdc_p_check::HvdcPPlan plan_;
};

}  // namespace ls2g

#endif  // HVDC_AC_EMULATION_LIMITS_LOOP_H
