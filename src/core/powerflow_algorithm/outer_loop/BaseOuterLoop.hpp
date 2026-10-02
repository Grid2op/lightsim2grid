// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_OUTER_LOOP_H
#define BASE_OUTER_LOOP_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Utils.hpp"
#include "AlgoConfig.hpp"
#include "batch_algorithm/LimitViolation.hpp"

#include "Eigen/Core"

namespace ls2g {

class LSGrid;

/**
 * What an outer loop's check() concluded, as OpenLoadFlow's OuterLoopStatus:
 *  - STABLE: nothing to change, the loop is done for this pass;
 *  - UNSTABLE: the loop changed something, the Newton has to run again;
 *  - FAILED: the loop cannot go on (eg the slack could not be distributed); the whole
 *    solve stops there and is reported as failed.
 */
enum class OuterLoopStatus { STABLE, UNSTABLE, FAILED };

/**
 * Collects what the outer loops reserve in the Newton-Raphson's Jacobian, before its
 * sparsity is built. A loop asks for every slot ANY of its states may need -- that union
 * is what lets the whole solve run on one symbolic analysis; afterwards a loop only
 * rewrites values. Solver bus ids throughout.
 */
class LS2G_API OuterDeclaration final
{
    public:
        /// a bus that is PV in the labelling but may become PQ (or back) during the solve:
        /// it gets a Vm unknown and a Q equation, see Base::set_switchable_vm_buses
        void add_switchable_vm_bus(int solver_bus_id) { switchable_vm_buses_.push_back(solver_bus_id); }

        const std::vector<int> & switchable_vm_buses() const { return switchable_vm_buses_; }
        void clear() { switchable_vm_buses_.clear(); }

    private:
        std::vector<int> switchable_vm_buses_;
};

/**
 * What the outer-loop algorithm lets a loop edit between two Newton solves, all of it
 * private to the algorithm (the grid is never modified). Absent in detection mode.
 *
 * `Sbus` is the injection the next Newton solve reads (solver numbering, pu, generation
 * positive); `Sbus_init` the one the grid handed in, untouched by any loop.
 */
struct OuterState
{
    CplxVect * Sbus = nullptr;
    const CplxVect * Sbus_init = nullptr;
    /// the active power target (MW) a loop set for each generator / storage unit, NaN where
    /// it kept the grid's; empty until a loop sizes it. Storage units in the load convention,
    /// as the grid stores them. Published by LSGrid::compute_results.
    std::vector<real_type> gen_target_p;
    std::vector<real_type> storage_target_p;
};

/**
 * Everything a loop sees, after a Newton solve. Built by the outer-loop algorithm between
 * two solves (with `state`), or by the grid after any AC solve for detection (`state` is
 * then null: a loop reports what it WOULD do, from the grid and the solve alone).
 *
 * Solver bus numbering for the vectors; `grid` gives the element data and the
 * grid <-> solver bus maps.
 */
struct OuterContext
{
    const LSGrid * grid = nullptr;
    const CplxVect * V = nullptr;              ///< complex voltages, pu
    const CplxVect * bus_mismatch = nullptr;   ///< see BaseAlgo::get_bus_mismatch
    const RealVect * controller_q = nullptr;   ///< see BaseAlgo::get_controller_q
    int slack_bus = -1;                        ///< solver id of the (single) slack bus
    /// the distributed slack's unknown, pu (see BaseAlgo::get_slack_absorbed): 0 for a
    /// single-slack Newton, where the slack bus' mismatch carries the whole imbalance
    real_type slack_absorbed = 0.;
    /// how many times THIS loop was unstable so far in the solve (OpenLoadFlow's
    /// `context.getIteration()`: 0 means it has not changed anything yet)
    int iteration = 0;
    OuterState * state = nullptr;

    bool is_detection() const { return state == nullptr; }
};

/**
 * What the last outer-loop solve did, for tests and diagnostics.
 */
struct LS2G_API OuterLoopStats
{
    /// final status: STABLE, or UNSTABLE when the iteration cap was reached, or FAILED
    OuterLoopStatus status = OuterLoopStatus::STABLE;
    /// name of the loop that failed, empty otherwise
    std::string failed_loop;
    /// whether the unrealistic-voltage check failed (status is then FAILED)
    bool unrealistic_state = false;
    int nb_outer_iterations = 0;   ///< over every loop, capped by max_outer_iterations
    int nb_passes = 0;
    /// (loop name, number of UNSTABLE results), in the order of the loop list
    std::vector<std::pair<std::string, int> > loop_iterations;
    /// Newton iterations of every inner solve, the first solve included
    std::vector<int> nr_iterations;

    void clear() {
        status = OuterLoopStatus::STABLE;
        failed_loop.clear();
        unrealistic_state = false;
        nb_outer_iterations = 0;
        nb_passes = 0;
        loop_iterations.clear();
        nr_iterations.clear();
    }
};

/**
 * One OpenLoadFlow outer loop.
 *
 * Same contract as the element containers (GenericContainer): the public entry points are
 * non-virtual and never overridden, each forwards to one protected `_xxx` hook. A loop
 * carries its own per-solve state, (re)set by initialize(); the algorithm holds a clone per
 * instance, so a loop object is never shared between two solves running at once.
 *
 * The trigger of a loop is computed in ONE place, `_detect`: check() acts on what it
 * reports, and the grid's detection mode (no outer loop run) calls it as is.
 */
class LS2G_API BaseOuterLoop
{
    public:
        virtual ~BaseOuterLoop() = default;

        /// the OpenLoadFlow name of the loop
        std::string name() const { return _name(); }

        /// reserve, in the Jacobian, every slot any state of this loop may use. Called
        /// when the solver input is (re)built, with the grid data known.
        void declare(const OuterContext & ctx, OuterDeclaration & decl) const { _declare(ctx, decl); }

        /// OpenLoadFlow's isNeeded: whether the loop has anything to do on this grid
        bool is_needed(const OuterContext & ctx) const { return _is_needed(ctx); }

        /// before the first Newton solve
        void initialize(OuterContext & ctx) { _initialize(ctx); }

        /// the trigger, as violations; nothing is modified
        void detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const { _detect(ctx, out); }

        /// detect, then act on the algorithm's state
        OuterLoopStatus check(OuterContext & ctx) { return _check(ctx); }

        /// after the last Newton solve; called in the reverse order of the loop list
        void cleanup(OuterContext & ctx) { _cleanup(ctx); }

        /// whether this loop can bring an unrealistic voltage back (OpenLoadFlow's
        /// canFixUnrealisticState): the check is deferred until after the last such loop
        bool can_fix_unrealistic_state() const { return _can_fix_unrealistic_state(); }

        std::unique_ptr<BaseOuterLoop> clone() const { return _clone(); }

        /// the loop's parameters, flattened like an AlgoConfig (persistence, copies)
        AlgoConfig get_params() const { return _get_params(); }
        void set_params(const AlgoConfig & params) { _set_params(params); }

    protected:
        virtual std::string _name() const = 0;
        virtual void _declare(const OuterContext & /*ctx*/, OuterDeclaration & /*decl*/) const {}
        virtual bool _is_needed(const OuterContext & /*ctx*/) const { return true; }
        virtual void _initialize(OuterContext & /*ctx*/) {}
        virtual void _detect(const OuterContext & ctx, std::vector<LimitViolation> & out) const = 0;
        virtual OuterLoopStatus _check(OuterContext & ctx) = 0;
        virtual void _cleanup(OuterContext & /*ctx*/) {}
        virtual bool _can_fix_unrealistic_state() const { return false; }
        virtual std::unique_ptr<BaseOuterLoop> _clone() const = 0;
        virtual AlgoConfig _get_params() const { return AlgoConfig(); }
        virtual void _set_params(const AlgoConfig & /*params*/) {}
};

}  // namespace ls2g

#endif  // BASE_OUTER_LOOP_H
