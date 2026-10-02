// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DISTRIBUTED_SLACK_LOOP_H
#define DISTRIBUTED_SLACK_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "element_container/SlackRedistribution.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's DistributedSlack outer loop (DistributedSlackOuterLoop, 2.3.0, balance
 * type PROPORTIONAL_TO_GENERATION_P_MAX), around a single-slack Newton.
 *
 * After a solve, the active power the slack bus absorbed (its mismatch) is shared on the
 * units taking part in the slack: every generator and storage unit with a "can participate
 * in the slack" weight (LSGrid::set_gen_can_participate_slack, the key OpenLoadFlow uses,
 * see _olf_rules.participation_weight), each bounded by its P limits (LSGrid::
 * set_gen_p_limits, OpenLoadFlow's target range) and never pushed across 0 MW. The sharing
 * is OpenLoadFlow's: not incremental, every pass starts again from the units' initial
 * targets with the cumulative mismatch (slack_redistribution::distribute). The new targets
 * go into the algorithm's injection (OuterState::Sbus) and are published as the units' P.
 *
 * Trigger (detect): the mismatch above `slack_bus_p_max_mismatch_mw` (SLACK_MISMATCH).
 *
 * Only on a grid with some unit flagged "can participate in the slack" (is_needed, and the
 * same rule in detection mode): a grid where no unit is flagged has no distributed slack set
 * up at all -- OpenLoadFlow with distributedSlack off, where the loop does not exist -- and
 * its single slack is kept on purpose.
 * UNSTABLE when the targets moved, FAILED when some of the mismatch could not be shared
 * (every unit at a bound) and `fail_on_residue` (OpenLoadFlow's FAIL behaviour; off, its
 * LEAVE_ON_SLACK_BUS: the slack bus keeps it).
 */
class LS2G_API DistributedSlackLoop final : public BaseOuterLoop
{
    public:
        /// OpenLoadFlow's slackBusPMaxMismatch, MW
        real_type slack_bus_p_max_mismatch_mw = 1.;
        /// OpenLoadFlow's slackDistributionFailureBehavior: FAIL (true) or LEAVE_ON_SLACK_BUS
        bool fail_on_residue = true;

        /// OpenLoadFlow's P_RESIDUE_EPS (1e-5 pu of its 100 MVA base), MW: what is left to
        /// share below it is nothing
        static constexpr real_type P_RESIDUE_EPS_MW = 1e-3;

    protected:
        std::string _name() const override { return "DistributedSlack"; }
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new DistributedSlackLoop(*this));
        }
        AlgoConfig _get_params() const override;
        void _set_params(const AlgoConfig & params) override;

    private:
        /// the mismatch of the slack bus (MW, > 0: the units must inject more), NaN without
        /// a solved slack bus
        real_type _mismatch_mw(const OuterContext & ctx) const;
        /// OpenLoadFlow's trigger: the mismatch (written in `mismatch_mw`) above the threshold
        bool _triggered(const OuterContext & ctx, real_type & mismatch_mw) const;
        /// whether some unit of `grid` takes part in the slack
        static bool _has_participant(const LSGrid & grid);

        // per solve (see _initialize): the participants, with their initial injection
        // (generator convention, MW), and where each one is now
        std::vector<slack_redistribution::Participant> units_;
        std::vector<int> unit_solver_bus_;
        std::vector<real_type> current_mw_;
};

}  // namespace ls2g

#endif  // DISTRIBUTED_SLACK_LOOP_H
