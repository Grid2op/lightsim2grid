// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SHUNT_VOLTAGE_CONTROL_LOOP_H
#define SHUNT_VOLTAGE_CONTROL_LOOP_H

#include <set>
#include <vector>

#include "BaseOuterLoop.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's ShuntVoltageControl outer loop in its WITH_GENERATOR_VOLTAGE_CONTROL mode
 * (ShuntVoltageControlOuterLoop, 2.3.0), what its shuntCompensatorVoltageControlOn creates
 * with that mode. Not in the default list (that parameter is off by default): add it to run it.
 *
 * The controllers: as OpenLoadFlow, the shunts of one bus whose sections regulate a voltage
 * (LSGrid::set_shunt_section_regulation) are ONE controller, its susceptance the sum of theirs.
 * A controller regulates the bus its first shunt does; the controllers regulating one bus form
 * a group, the first one's target holding. A group whose bus a generator, an SVC, a VSC station
 * or a transformer (when TransformerVoltageControl comes before this loop) regulates is hidden
 * (OpenLoadFlow's control priorities) and never acts: neither solved for nor rounded
 * (OpenLoadFlow's getControllerElements keeps the visible controls).
 *
 * The first solve has every controller on: the Newton solves for the susceptances
 * (ShuntControl). At the loop's first check each controller is switched off and its
 * susceptance shared out between its shunts, the largest first, each rounded to its closest
 * section (LfShuntImpl.dispatchB), and the Newton runs again; then it is
 * stable. Everything is declared up front, so the whole solve is one symbolic analysis. The
 * section counts it leaves are in the results (ShuntInfo::res_section_count), the inputs are
 * not modified.
 *
 * Trigger (detect): a visible group off its target by more than half its deadband (by more
 * than the detection tolerance without one), SHUNT_VOLTAGE_CONTROL, ViolationCategory::CONTROL.
 */
class LS2G_API ShuntVoltageControlLoop final : public BaseOuterLoop
{
    public:
        struct Group {
            int bus_solver = -1;
            int bus_grid = -1;
            real_type target = 1.;          ///< pu
            real_type half_deadband = 0.;   ///< pu, 0 without one
            std::vector<int> controller_buses;           ///< solver ids
            std::vector<std::vector<int> > shunts;       ///< per controller bus
            bool hidden = false;
        };
        /// the shunt voltage controls of `grid`, `transformer_buses` the buses (solver ids)
        /// transformers regulate, see the class comment
        static std::vector<Group> groups(const LSGrid & grid, const std::set<int> & transformer_buses);

    protected:
        std::string _name() const override { return "ShuntVoltageControl"; }
        void _declare(const OuterContext & ctx, OuterDeclaration & decl) const override;
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new ShuntVoltageControlLoop(*this));
        }

    private:
        // the counts LfShuntImpl.dispatchB gives the shunts of one controller for `b` (pu)
        std::vector<int> _dispatch(const LSGrid & grid, const std::vector<int> & shunts, real_type b) const;

        // the buses of the ratio groups reserved before this loop (see _declare)
        mutable std::set<int> transformer_buses_;
        std::vector<Group> groups_;  // per solve, see _initialize
};

}  // namespace ls2g

#endif  // SHUNT_VOLTAGE_CONTROL_LOOP_H
