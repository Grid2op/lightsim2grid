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
 * go into the algorithm's injection (OuterInjections::Sbus) and are published as the units' P.
 *
 * Trigger (detect): the mismatch above `Params::slack_bus_p_max_mismatch_mw` (SLACK_MISMATCH).
 *
 * Only on a grid with some unit OFF the slack bus flagged "can participate in the slack"
 * (is_needed, and the same rule in detection mode): a grid where no unit is flagged has no
 * distributed slack set up at all -- OpenLoadFlow with distributedSlack off, where the loop
 * does not exist -- and its single slack is kept on purpose. Units on the slack bus alone (a
 * slack unit keeps its flag, see SlackParticipation::set_can_participate) are the same case:
 * they would only take what that bus already absorbs, the same solution, so the loop neither
 * runs nor fails on their limits -- such a unit past its active range is GenPCheck's LOW_P /
 * HIGH_P.
 * A batch row whose contingency loses active power (OuterRowView::lost_p_mw, a batch's
 * `redistribute_slack`) has it shared on the participants before its first solve, the same
 * way, and those targets are the initial ones the loop then starts from: OpenLoadFlow's
 * security analysis (its Default ContingencyActivePowerLossDistribution) before its
 * post-contingency load flow. Without it, the first solve puts the whole loss on the slack bus.
 *
 * UNSTABLE when the targets moved, FAILED when some of the mismatch could not be shared
 * (every unit at a bound) and `Params::fail_on_residue` (OpenLoadFlow's FAIL behaviour; off, its
 * LEAVE_ON_SLACK_BUS: the slack bus keeps it).
 *
 * Around a distributed-slack Newton (NROuterDist_*, OuterControls::supports_slack_share) the
 * Newton shares the slack itself, by weight and without bounds, on the same units: the loop
 * reserves their buses (declare) and hands it their weights (SlackShareControl). After a solve
 * the driver has moved each unit's target by its share; the loop then shares the same cumulative
 * mismatch from the initial targets as above. A unit that sharing saturates -- above `max_p`,
 * below `min_p`, across 0 MW -- is set at its bound and takes no share in the Newton any more
 * (weight 0); one that comes back inside its bounds gets its weight back. When no unit is pushed
 * past a bound the targets do not move: the loop is stable after the first solve. A solve resumed
 * from another (a batch row from its "n" solve, a time series' next step) starts with the units
 * that solve left out still out, and a batch row with the units the share of its lost power put
 * at a bound out as well. Both bet that the next solve pushes them the same way; when it does
 * not, the check gives them their weight back, at the cost of a re-solve. Every unit at a
 * bound: FAILED as above, or the reference bus takes the rest (LEAVE_ON_SLACK_BUS). The
 * dispatch differs from OpenLoadFlow's by at most `Params::slack_bus_p_max_mismatch_mw`, which
 * OpenLoadFlow leaves on the reference bus and the Newton shares out.
 */
class LS2G_API DistributedSlackLoop final : public BaseOuterLoop
{
    public:
        /// The loop's parameters, OpenLoadFlow's values by default; fixed at construction.
        struct Params {
            /// OpenLoadFlow's slackBusPMaxMismatch, MW
            real_type slack_bus_p_max_mismatch_mw = 1.;
            /// OpenLoadFlow's slackDistributionFailureBehavior: FAIL (true) or LEAVE_ON_SLACK_BUS
            bool fail_on_residue = true;
            /// OpenLoadFlow's P_RESIDUE_EPS (1e-5 pu of its 100 MVA base), MW: what is left to
            /// share below it is nothing
            real_type p_residue_eps_mw = 1e-3;
            /// OpenLoadFlow's PreviousStateInfo.moved: the units moved when they moved by more
            /// than this fraction of p_residue_eps_mw (its margin against rounding)
            real_type moved_fraction = 0.9;
        };

        DistributedSlackLoop();
        explicit DistributedSlackLoop(const Params & params);
        const Params & params() const { return params_; }

        /// What a row's DistributedSlack shares before its first solve (MW, generator
        /// convention), NaN when nothing is to be shared: the row's OuterRowView::lost_p_mw,
        /// counted at the units' row or grid targets, moved to the targets in force OpenLoadFlow
        /// counts the lost units at -- the ones `resumed_gen_p` / `resumed_sto_p` (a resumed
        /// solve's moved targets, MW, NaN untouched, storage units in the load convention; empty:
        /// none) give the generators and storage units the row disconnects or strands, unless
        /// the row sets their target itself.
        static real_type lost_p_to_share_mw(const OuterGridView & view,
                                            const std::vector<real_type> & resumed_gen_p,
                                            const std::vector<real_type> & resumed_sto_p);

    protected:
        std::string _name() const override { return "DistributedSlack"; }
        void _declare(const OuterContext & ctx) const override;
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new DistributedSlackLoop(*this));
        }
        AlgoConfig _get_params() const override;
        void _adjoint_terms(OuterAdjointTerms & out) const override;

    private:
        /// the mismatch of the slack bus (MW, > 0: the units must inject more), NaN without
        /// a solved slack bus
        real_type _mismatch_mw(const OuterContext & ctx) const;
        /// OpenLoadFlow's trigger: the mismatch (written in `mismatch_mw`) above the threshold
        bool _triggered(const OuterContext & ctx, real_type & mismatch_mw) const;
        /// whether some unit of the grid `view` shows takes part in the slack, on another
        /// solver bus than `skip_solver_bus` (-1: anywhere)
        static bool _has_participant(const OuterGridView & view, int skip_solver_bus = -1);
        /// the row's loss shared on the participants (see lost_p_to_share_mw), their targets
        /// moved and taken as the initial ones; nothing without a loss to share
        void _share_lost_p(OuterContext & ctx);
        /// whether the Newton shares the slack itself (see the class comment)
        static bool _in_newton(const OuterContext & ctx);
        /// the units' targets in force now, MW, generator convention
        void _targets_in_force(const OuterContext & ctx, std::vector<real_type> & out) const;
        /// hand the Newton the units' weights, 0 for the parked ones
        void _set_newton_weights(OuterContext & ctx) const;
        /// unit k's upper (`high`) or lower bound, MW, as slack_redistribution::distribute has it
        /// (0 MW on the side of 0 its initial target is not on); NaN without one
        real_type _bound(std::size_t k, bool high) const;
        /// park the units the solve this one resumes from parked (`was_parked`: kind, id, at the
        /// upper bound), their targets moved to that bound
        void _keep_parked(OuterContext & ctx,
                          const std::vector<std::pair<std::pair<slack_redistribution::UnitKind, int>, char> > & was_parked);
        /// _check around a distributed-slack Newton
        OuterLoopStatus _check_in_newton(OuterContext & ctx);

        // per solve (see _initialize): the participants, with their initial injection
        // (generator convention, MW), and where each one is now
        std::vector<slack_redistribution::Participant> units_;
        std::vector<int> unit_solver_bus_;
        std::vector<real_type> current_mw_;
        // around a distributed-slack Newton: the units at a bound, out of its share, and for each
        // whether at its upper one (else its lower one)
        std::vector<char> parked_;
        std::vector<char> parked_high_;
        // the units the row's lost-power share (_share_lost_p) put at a bound, and on which side
        std::vector<char> lost_p_saturated_;
        bool lost_p_high_ = false;
        real_type sn_mva_ = 0.;  // the grid's, at the last initialize

        const Params params_;
};

}  // namespace ls2g

#endif  // DISTRIBUTED_SLACK_LOOP_H
