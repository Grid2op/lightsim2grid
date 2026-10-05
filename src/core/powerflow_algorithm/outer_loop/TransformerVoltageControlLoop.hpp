// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef TRANSFORMER_VOLTAGE_CONTROL_LOOP_H
#define TRANSFORMER_VOLTAGE_CONTROL_LOOP_H

#include <vector>

#include "BaseOuterLoop.hpp"
#include "batch_algorithm/BusQCheck.hpp"

namespace ls2g {

class LSGrid;

/**
 * OpenLoadFlow's TransformerVoltageControl outer loop in its AFTER_GENERATOR_VOLTAGE_CONTROL
 * mode (TransformerVoltageControlOuterLoop, 2.3.0), what its transformerVoltageControlOn
 * creates with that mode. Not in the default list (that parameter is off by default): add it to
 * run it.
 *
 * The controllers: a transformer whose ratio tap changer regulates a voltage (LSGrid::
 * set_trafo_ratio_tap_regulation), connected at both ends, on two different buses. They are
 * grouped by the bus they regulate, in transformer order: the first one's target holds, the
 * smallest deadband (none: 0.1 kV). A group whose bus a generator, an SVC or a VSC station
 * already regulates is hidden (OpenLoadFlow's control priorities) and never acts: neither
 * switched on nor rounded (OpenLoadFlow's getControllerElements keeps the visible controls).
 * Its transformers still keep a generator of their buses frozen (hasStepUpTransformers).
 *
 * Its three steps:
 *  - INITIAL (first check): every group off its target by more than half its deadband switches
 *    its transformers on. If none, the loop is done. Otherwise the generators regulating a bus
 *    of at most `Params::max_controlled_nominal_voltage` kV (step-up units aside) are frozen at the
 *    reactive power their bus injects -- OpenLoadFlow's GeneratorVoltageControlManager, which
 *    passes the bus' injection, its load included, as the generation; reproduced -- and a
 *    transformer with no PV bus left on its other side is switched off
 *    (fixTransformerVoltageControls). The Newton then solves for the ratios (BranchControl);
 *  - CONTROL: a transformer whose ratio left its group's range is rounded to its extreme tap
 *    and switched off, and the Newton runs again; once none does, each ratio is mapped back
 *    around its initial tap (TransformerRatioManager, `use_initial_tap_position`), every
 *    transformer is rounded to its closest tap and switched off, and the generators released;
 *  - COMPLETE: stable.
 * Everything is declared up front (a ratio column per transformer of a visible group, the
 * union of its rows; the switchable buses and held controllers of the freeze), so the whole
 * solve is one symbolic analysis. The taps it leaves are in the results
 * (TrafoInfo::res_ratio_tap_position), the inputs are not modified.
 *
 * Trigger (detect): a visible group off its target by more than half its deadband,
 * TRANSFORMER_VOLTAGE_DEADBAND, ViolationCategory::CONTROL.
 */
class LS2G_API TransformerVoltageControlLoop final : public BaseOuterLoop
{
    public:
        /// The loop's parameters, OpenLoadFlow's values by default; fixed at construction.
        struct Params {
            /// OpenLoadFlow's transformerVoltageControlUseInitialTapPosition
            bool use_initial_tap_position = true;
            /// OpenLoadFlow's generatorVoltageControlMinNominalVoltage (kV): the generators
            /// regulating a bus of at most this are frozen while the transformers act; < 0: the
            /// highest nominal voltage the transformers regulate (OpenLoadFlow's automatic value)
            real_type max_controlled_nominal_voltage = 120.;
            /// the deadband of a group none of whose transformers has one, kV
            /// (AbstractTransformerVoltageControlOuterLoop's MIN_TARGET_DEADBAND_KV)
            real_type min_target_deadband_kv = 0.1;
        };

        TransformerVoltageControlLoop();
        explicit TransformerVoltageControlLoop(const Params & params);
        const Params & params() const { return params_; }

        struct Group {
            int bus_solver = -1;
            int bus_grid = -1;
            real_type target = 1.;         ///< pu
            real_type half_deadband = 0.;  ///< pu
            std::vector<int> trafos;
            bool hidden = false;
        };
        /// the transformer voltage controls of `grid`, see the class comment
        std::vector<Group> groups(const LSGrid & grid) const;

    protected:
        std::string _name() const override { return "TransformerVoltageControl"; }
        void _declare(const OuterContext & ctx) const override;
        bool _is_needed(const OuterContext & ctx) const override;
        void _initialize(OuterContext & ctx) override;
        void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const override;
        OuterLoopStatus _check(OuterContext & ctx) override;
        bool _can_fix_unrealistic_state() const override { return true; }
        std::unique_ptr<BaseOuterLoop> _clone() const override {
            return std::unique_ptr<BaseOuterLoop>(new TransformerVoltageControlLoop(*this));
        }
        AlgoConfig _get_params() const override;

    private:
        enum class Step { INITIAL, CONTROL, COMPLETE };
        // TransformerRatioManager's, per transformer
        struct Ratio {
            real_type initial = 1.;  // its own, when switched on
            real_type min = 1., max = 1.;            // its own range
            real_type shared_min = 1., shared_max = 1., shared_initial = 1.;
        };
        // a generator bus frozen by INITIAL
        struct Frozen {
            int entry;
            int bus_solver;
            bool local;
            real_type target_vm;
        };

        real_type _limit(const LSGrid & grid) const;
        bool _step_up(const LSGrid & grid, const bus_q_check::BusQEntry & entry, real_type limit) const;
        void _freeze_generators(OuterContext & ctx, real_type limit);
        void _release_generators(OuterContext & ctx);
        void _fix_controls(OuterContext & ctx);
        int _closest_tap(const OuterContext & ctx, int trafo, real_type value) const;
        void _set_on(OuterContext & ctx, int trafo, bool on);

        // per solve, see _initialize
        std::vector<Group> groups_;
        bus_q_check::BusQPlan plan_;
        Step step_ = Step::INITIAL;
        std::vector<char> enabled_;        // per transformer
        std::vector<char> controller_;     // per transformer: in a group
        std::vector<Ratio> ratios_;        // per transformer
        std::vector<Frozen> frozen_;

        const Params params_;
};

}  // namespace ls2g

#endif  // TRANSFORMER_VOLTAGE_CONTROL_LOOP_H
