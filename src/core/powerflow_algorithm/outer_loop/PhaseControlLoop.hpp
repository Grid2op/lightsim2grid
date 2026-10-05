// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef PHASE_CONTROL_LOOP_H
#define PHASE_CONTROL_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "element_container/TapChangers.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's PhaseControl outer loop in its CONTINUOUS_WITH_DISCRETISATION mode
 * (PhaseControlOuterLoop, 2.3.0), the one its phaseShifterRegulationOn creates. Not in the
 * default list (that parameter is off by default): add it to run it.
 *
 * The phase shifters: a transformer whose phase tap changer regulates (LSGrid::
 * set_trafo_phase_tap_regulation), connected at both ends, on two different buses, through
 * one of its own sides:
 *  - ACTIVE_POWER (OpenLoadFlow's CONTROLLER): the Newton first solves for the shift that
 *    gives the target active power on the regulated side (BranchControl); at the loop's first
 *    check the control is switched off and the shift rounded to the closest tap. A
 *    controller whose transformer is needed for the connectivity of the grid does not
 *    regulate (fixPhaseShifterNecessaryForConnectivity);
 *  - CURRENT_LIMITER (LIMITER): from the second check on, a current above the limit moves
 *    the tap one position, the way that lowers it (the sign of dI/da).
 * Every one of them is reserved up front (a PhaseShifterControl; BranchControl: a column and a
 * row per controller, the block of every one patched by value), so the whole solve is one
 * symbolic analysis.
 *
 * One known difference with OpenLoadFlow's current limiter. Its one-tap move
 * (PiModelArray.shiftOneTapPositionToChangeA1) steps to the next position, then compares
 * the previous one -- the starting tap -- with the shift the Newton wrote back, which is
 * that tap's alpha up to the rounding of its linear solve: on the wrong side of it, the
 * move is undone and the limiter stays where it was although above its limit. Which way the
 * rounding falls is not reproducible here; this loop compares with the tap's own alpha, so a
 * move towards the next position is always made.
 *
 * Trigger (detect): an active power off its target by more than the deadband (or the
 * detection tolerance, whichever is larger), PHASE_CONTROL_P; a current above its limit,
 * PHASE_LIMITER_CURRENT. Both ViolationCategory::CONTROL.
 */
class LS2G_API PhaseControlLoop final : public BaseOuterLoop
{
    protected:
        std::string _name() const override { return "PhaseControl"; }
        void _declare(const OuterContext & ctx) const override;
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new PhaseControlLoop(*this));
        }

    public:
        struct Shifter {
            int trafo;
            RegulationMode mode;  // ACTIVE_POWER or CURRENT_LIMITER
            int side;             // the regulated side
            bool regulates;       // false: an active power one needed for the connectivity
        };
        /// the phase shifters of `grid`, see the class comment
        static std::vector<Shifter> shifters(const LSGrid & grid);

    private:
        std::vector<Shifter> shifters_;  // per solve, see _initialize
};

}  // namespace ls2g

#endif  // PHASE_CONTROL_LOOP_H
