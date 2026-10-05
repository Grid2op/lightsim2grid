// Copyright (c) 2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_OUTER_LOOP_H
#define BASE_OUTER_LOOP_H

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <set>
#include <utility>
#include <vector>

#include "Utils.hpp"
#include "TaggedIdVec.hpp"
#include "AlgoConfig.hpp"
#include "batch_algorithm/LimitViolation.hpp"
#include "OuterControls.hpp"

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
        /// a transformer (grid id) whose phase tap may move during the solve and, with
        /// `solves_shift`, whose shift the Newton solves for: see BranchControl
        void add_phase_shifter(int trafo_id, bool solves_shift) {
            phase_shifters_.push_back(trafo_id);
            phase_shifter_column_.push_back(solves_shift ? 1 : 0);
        }

        const std::vector<int> & phase_shifters() const { return phase_shifters_; }
        const std::vector<char> & phase_shifter_column() const { return phase_shifter_column_; }
        /// transformers (grid ids, in order) regulating the voltage of `bus_solver` at
        /// `target_vm` pu: with `solved`, the Newton solves for their ratios (BranchControl),
        /// otherwise only their taps may move
        void add_ratio_group(int bus_solver, real_type target_vm, const std::vector<int> & trafo_ids, bool solved) {
            ratio_groups_.push_back(RatioGroup{bus_solver, target_vm, trafo_ids, solved});
        }
        struct RatioGroup {
            int bus_solver;
            real_type target_vm;
            std::vector<int> trafos;
            bool solved;
        };
        const std::vector<RatioGroup> & ratio_groups() const { return ratio_groups_; }
        /// shunts regulating the voltage of `bus_solver` at `target_vm` pu, from the controller
        /// buses `controller_buses` (solver ids), each with its regulating shunts (grid ids):
        /// with `solved`, the Newton solves for their susceptances (ShuntControl), otherwise
        /// only their sections may move
        void add_shunt_group(int bus_solver, real_type target_vm, const std::vector<int> & controller_buses,
                             const std::vector<std::vector<int> > & shunts, bool solved) {
            shunt_groups_.push_back(ShuntGroup{bus_solver, target_vm, controller_buses, shunts, solved});
        }
        struct ShuntGroup {
            int bus_solver;
            real_type target_vm;
            std::vector<int> controller_buses;
            std::vector<std::vector<int> > shunts;
            bool solved;
        };
        const std::vector<ShuntGroup> & shunt_groups() const { return shunt_groups_; }
        void clear() {
            phase_shifters_.clear(); phase_shifter_column_.clear(); ratio_groups_.clear();
            shunt_groups_.clear();
        }

    private:
        std::vector<int> phase_shifters_;
        std::vector<char> phase_shifter_column_;
        std::vector<RatioGroup> ratio_groups_;
        std::vector<ShuntGroup> shunt_groups_;
};

/**
 * What the outer-loop algorithm lets a loop edit between two Newton solves, all of it
 * private to the algorithm (the grid is never modified). Absent in detection mode.
 *
 * `Sbus` is the injection the next Newton solve reads (solver numbering, pu, generation
 * positive); `Sbus_init` the one the grid handed in, untouched by any loop; `Sbus_target`
 * that one with the units' targets a loop moved (DistributedSlack's target P, and the target
 * Q that follows it, see LSGrid::set_gen_raw_target_q): a unit's own output, not a residual.
 */
struct OuterState
{
    CplxVect * Sbus = nullptr;
    const CplxVect * Sbus_init = nullptr;
    CplxVect * Sbus_target = nullptr;
    /// the active power target (MW) a loop set for each generator / storage unit, NaN where
    /// it kept the grid's; empty until a loop sizes it. Storage units in the load convention,
    /// as the grid stores them. Published by LSGrid::compute_results.
    std::vector<real_type> gen_target_p;
    std::vector<real_type> storage_target_p;
    /// the droop regime (0 linear, +1 saturated 1 -> 2, -1 saturated 2 -> 1) a loop set for
    /// each hvdc line (grid id), HVDC_KEEP where it kept the grid's; empty until a loop sizes it
    std::vector<int> hvdc_status;
    /// the set-point (pu) a loop switched each idle standby SVC on at (grid id), NaN where it
    /// is still held at Q = 0; empty until a loop sizes it
    std::vector<real_type> svc_target_vm;
    /// the position a loop moved each transformer's phase tap to (grid id), TAP_KEEP where it
    /// kept the solve's; empty until a loop sizes it
    std::vector<int> phase_tap;
    /// whether the Newton solves each declared transformer's shift for its active power
    /// (grid id): 1 on, 0 off, -1 kept as it is; empty until a loop sizes it
    std::vector<int> phase_control;
    /// the same for the ratio tap changers (moved, TAP_KEEP) and their voltage control
    /// (1 on, 0 off, -1 kept); the moves are applied by the next solve, then forgotten
    std::vector<int> ratio_tap;
    std::vector<int> ratio_control;
    /// the shunt controllers (by controller bus, solver id): their voltage control (1 on, 0 off,
    /// -1 kept; empty until a loop sizes it), and the section counts a loop switched their shunts
    /// to (applied by the next solve, then forgotten)
    std::vector<int> shunt_control;
    std::vector<std::pair<int, std::vector<int> > > shunt_sections;
    static constexpr int HVDC_KEEP = 2;
    static constexpr int TAP_KEEP = std::numeric_limits<int>::min();
};

class BranchControl;
class ShuntControl;

/**
 * One decision of an outer loop in a check: what it acted on, or a test it ran and did not act
 * on, with the two numbers it compared. Two runs that part ways show it here first, with the
 * margin that decided it (OuterLoopStats::decisions).
 *
 * The element and the reason use LimitViolation's conventions: a grid-model bus id for BUS,
 * the element's own id otherwise, -1 for GRID; `value` / `limit` in the unit of `reason`
 * (see LimitViolation), except where an action says otherwise. The actions:
 *
 * - DistributedSlack: DISTRIBUTE (every check: the slack mismatch, MW, against the threshold),
 *   UNITS_MOVED (how much the units moved, MW, against what counts as a move), FAIL_RESIDUE.
 * - ReactiveLimits, on a controller bus: PV_TO_PQ (PV_TO_PQ_UNREALISTIC, the robust mode's),
 *   KEPT_PV_STRONGEST; PQ_TO_PV, KEPT_PQ_MAX_SWITCH, and, every check for every frozen bus,
 *   KEPT_PQ / KEPT_PQ_GROUP_HOLDS (the release test that did not trigger: the regulated voltage
 *   the solve gave and the set-point, kV; GROUP_HOLDS when another controller of the group
 *   settled it); LIMIT_MOVED (the new limit and the one it was frozen at, MVar).
 * - AcHvdcAcEmulationLimits: SATURATE, RELEASE. VoltageMonitoring: SWITCH_ON.
 * - PhaseControl: ROUND_TAP (the shift and its tap's, rad), MOVE_TAP / KEPT_TAP (the current,
 *   A, against the limiter's). TransformerVoltageControl: SWITCH_ON (|target - v| against the
 *   half deadband, pu), ROUND_TO_RANGE (the ratio and the end it left), ROUND_TAP (the ratio and
 *   the tap position). ShuntVoltageControl: ROUND_SECTIONS (the susceptance solved for, pu).
 */
struct LS2G_API OuterDecision
{
    std::string loop;           ///< the loop's name (filled by the driver)
    int outer_iteration = 0;    ///< OuterLoopStats::nb_outer_iterations at the check
    std::string action;         ///< what the test decides, eg "PV_TO_PQ" (see each loop)
    bool taken = true;          ///< whether the loop acted on it
    ViolationElementType element_type = ViolationElementType::GRID;
    int element_id = -1;
    LimitViolationType reason = LimitViolationType::NOT_SIMULATED;
    real_type value = std::numeric_limits<real_type>::quiet_NaN();
    real_type limit = std::numeric_limits<real_type>::quiet_NaN();
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
    const RealVect * Va = nullptr;             ///< the angles, rad (a DC solve's are only here)
    /// the magnitudes, pu, as the Newton holds them (its own unknowns): read through vm(),
    /// not as |V|, which rebuilds them through cos / sin / hypot and can land an ulp away --
    /// enough to flip a strict comparison with the set-point a bus is held at, differently
    /// from one CPU to the next. Null outside a Newton solve (a detection, a batch).
    const RealVect * Vm = nullptr;
    const CplxVect * bus_mismatch = nullptr;   ///< see BaseAlgo::get_bus_mismatch
    const RealVect * controller_q = nullptr;   ///< see BaseAlgo::get_controller_q
    int slack_bus = -1;                        ///< solver id of the slack bus, -1 when several
    /// the distributed slack's unknown, pu (see BaseAlgo::get_slack_absorbed): 0 for a
    /// single-slack Newton, where the slack bus' mismatch carries the whole imbalance
    real_type slack_absorbed = 0.;
    /// the tolerance of a detection, MW (the caller's, eg compute_physical_violations'; 0
    /// in the outer-loop mode, where OpenLoadFlow compares strictly)
    real_type tol_mw = 0.;
    /// the same for a voltage comparison, pu (0 in the outer-loop mode)
    real_type tol_vm_pu = 0.;
    /// whether the voltage magnitudes are the solve's: false after a DC solve, or in a
    /// batch whose algorithm cannot feed the reactive checks
    bool vm_checks = true;
    /// the solver buses this solve left out (a batch row's stranded buses), sorted; may be null
    const std::vector<int> * masked = nullptr;
    /// the grid -> solver bus map of the solve (null: the grid's AC one); a DC solve, a
    /// batch's own labelling, have theirs
    const SolverBusIdVect * id_me_to_solver = nullptr;
    /// how many times THIS loop was unstable so far in the solve (OpenLoadFlow's
    /// `context.getIteration()`: 0 means it has not changed anything yet)
    int iteration = 0;
    OuterState * state = nullptr;
    /// what the loops reserve and act through (null in detection)
    OuterControls * controls = nullptr;
    /// the transformers whose phase the solve handles (their shift, tap, current), null
    /// when there are none or in detection
    const BranchControl * branch_control = nullptr;
    /// the shunts whose susceptance the solve handles, null when there are none or in detection
    const ShuntControl * shunt_control = nullptr;
    /// where a check records its decisions (null: not recorded, eg in detection)
    std::vector<OuterDecision> * trace = nullptr;

    bool is_detection() const { return state == nullptr; }
    /// the voltage magnitude of solver bus `bus`, pu: Vm's, |V| when there is none
    real_type vm(int bus) const { return Vm != nullptr ? (*Vm)(bus) : std::abs((*V)(bus)); }

    /// record a decision in `trace` (nothing without one)
    void record(const char * action, bool taken, ViolationElementType element_type, int element_id,
                LimitViolationType reason, real_type value, real_type limit) const
    {
        if (trace == nullptr) return;
        OuterDecision d;
        d.action = action;
        d.taken = taken;
        d.element_type = element_type;
        d.element_id = element_id;
        d.reason = reason;
        d.value = value;
        d.limit = limit;
        trace->push_back(d);
    }
    /// the same for a bus given by its SOLVER id (recorded as its grid-model id)
    void record_bus(const char * action, bool taken, int solver_bus, LimitViolationType reason,
                    real_type value, real_type limit) const;
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
    /// every decision the loops took, and the tests they ran without acting, in order
    std::vector<OuterDecision> decisions;

    void clear() {
        status = OuterLoopStatus::STABLE;
        failed_loop.clear();
        unrealistic_state = false;
        nb_outer_iterations = 0;
        nb_passes = 0;
        loop_iterations.clear();
        nr_iterations.clear();
        decisions.clear();
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

        /// whether this loop needs the idle standby SVCs held in their voltage-control
        /// groups (VoltageControlPlan::build_controllers' `hold_monitors`): read by the grid
        /// when it builds the solver input, before any declaration
        bool holds_svc_monitors() const { return _holds_svc_monitors(); }

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
        virtual bool _holds_svc_monitors() const { return false; }
        virtual std::unique_ptr<BaseOuterLoop> _clone() const = 0;
        virtual AlgoConfig _get_params() const { return AlgoConfig(); }
        virtual void _set_params(const AlgoConfig & /*params*/) {}
};

}  // namespace ls2g

#endif  // BASE_OUTER_LOOP_H
