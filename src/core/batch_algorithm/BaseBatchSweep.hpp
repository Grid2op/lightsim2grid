// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEBATCHSWEEP_H
#define BASEBATCHSWEEP_H

#include "BaseBatchSolverSynch.hpp"
#include "YbusPolicy.hpp"
#include "TopoPlan.hpp"
#include "SbusPolicy.hpp"
#include "LimitViolation.hpp"
#include "OperationalCheck.hpp"
#include "GenPCheck.hpp"
#include "GenPvReleaseCheck.hpp"
#include "HvdcPCheck.hpp"
#include "powerflow_algorithm/outer_loop/DefaultOuterLoops.hpp"
#include "BusGraph.hpp"
#include "BatchAdjoint.hpp"
#include "PowerFlowAdjoint.hpp"
#include "element_container/SlackRedistribution.hpp"

#include <set>
#include <map>
#include <vector>
#include <memory>
#include <queue>
#include <algorithm>
#include <iterator>
#include <exception>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <limits>
// isnan / isfinite / sqrt: the C header, not <cmath>. These are used unqualified below
// (isnan(...), not std::isnan(...)), and only math.h guarantees them in the global
// namespace across compilers -- cmath only guarantees std::isnan / std::isfinite;
// whether it ALSO exposes the unqualified global name is implementation-defined, which
// is exactly why this compiled on GCC and recent Clang (both happen to expose it
// transitively) but not on an older Clang/libc++ pairing. Mirrors the pre-refactor
// ContingencyAnalysis.cpp's own include of this same header for the same reason.
#include <math.h>
#include <cmath>

namespace ls2g {

/**
How each element of a batch is initialized, ie which complex voltage the solver starts
from when it solves step `i`. Only meaningful when SbusPolicy::supports_vary (it is what
distinguishes TimeSeries from InjectionSweep; ContingencyAnalysis and ScenarioSweep both
use FromSeed, matching their "every step restarts from the same base case" semantics).
 **/
enum class LS2G_API BatchInitKind
{
    /** warm start: step `i` starts from the solution of step `i-1` (see TimeSeries) **/
    FromPreviousStep,
    /** every step restarts from the same seed, so the steps are independent
        of one another and of their ordering **/
    FromSeed
};

// the operational limit checks live in OperationalCheck.hpp (shared with LSGrid::get_violations)

/**
 * What the outer loops did on one row of a batch whose algorithm runs them (NROuter_*):
 * how the row ended, and what it settled on that the grid does not say -- the active targets
 * DistributedSlack moved, the hvdc regimes, the taps and sections. Only what differs from the
 * grid is listed, by element id, sorted: a row the loops left alone holds nothing but its
 * status.
 */
struct LS2G_API OuterRowResult
{
    /// whether the row converged (as converged_mask()): false for a row whose Newton diverged
    /// or that was never solved, which `status` alone does not tell
    bool converged = false;
    /// how the loops ended -- theirs only: STABLE too on a row whose Newton diverged, or that
    /// was never solved (OpenLoadFlow keeps the solver's status apart from the loops' as well)
    OuterLoopStatus status = OuterLoopStatus::STABLE;
    /// the unrealistic-voltage check stopped the row (see OuterLoopStats::unrealistic_state)
    bool unrealistic_state = false;
    /// the loop that failed, empty otherwise
    std::string failed_loop;
    int nb_outer_iterations = 0;
    /// Newton iterations over every inner solve of the row
    int nb_nr_iterations = 0;
    /// (id, MW): generators and storage units (load convention) whose target a loop moved
    std::vector<std::pair<int, real_type> > gen_target_p;
    std::vector<std::pair<int, real_type> > storage_target_p;
    /// (hvdc line id, regime, see HvdcRegimeControl)
    std::vector<std::pair<int, int> > hvdc_regime;
    /// (transformer id, position), (shunt id, section count)
    std::vector<std::pair<int, int> > phase_tap;
    std::vector<std::pair<int, int> > ratio_tap;
    std::vector<std::pair<int, int> > shunt_sections;
    /// ("loop/ACTION", count) of every decision the loops acted on (OuterDecision::taken),
    /// sorted: how many switches the row took, by kind (eg "ReactiveLimits/PV_TO_PQ")
    std::vector<std::pair<std::string, int> > decisions_taken;
};

/**
Batch powerflow algorithm, policy-parameterized on what varies from step to step:

- `YbusPolicy` (`YbusPolicy::NOOP` / `YbusPolicy::Contingency`): whether the admittance
  matrix varies (a per-step contingency -- line/trafo disconnection).
- `SbusPolicy` (`SbusPolicy::NOOP` / `SbusPolicy::Vary`): whether the injection varies
  (per-step generator / static generator / load values).
- `BatchInitKind`: only relevant when SbusPolicy::supports_vary -- whether each step is
  seeded from the previous step's result (chained, single-threaded) or independently
  from the same seed (parallelizable).

Only 3 of the 4 (YbusPolicy, SbusPolicy) combinations are instantiated -- "nothing
varies" is degenerate and excluded. See the aliases at the bottom of this file:
`TimeSeries`, `InjectionSweep`, `ContingencyAnalysis` and `ScenarioSweep` are all the
same class template, differing only in these compile-time policies.

Methods that only make sense for a subset of instantiations (eg `add_n1` -- a
ContingencyAnalysis-only concept) are member function *templates*, SFINAE-gated via
`std::enable_if` on the relevant policy's `supports_*` flag, so they are simply absent
(a compile error to call, "no matching member function") on instantiations they do not
apply to -- NOT a plain `static_assert` in a non-template method body, which would be
forced to compile by this file's explicit class-template instantiation (see the bottom
of `BaseBatchSweep.cpp`) even for instantiations that never call the method.

A member function TEMPLATE's definition must be visible in every translation unit that
calls it (ordinary template-instantiation rules; `extern template class` at the bottom
of this file only suppresses re-instantiation of the class's *non-template* members).
Since pybind11's bindings (`binding_batch.cpp`) call most of the SFINAE-gated methods
below directly, essentially all of them are defined inline, right here, rather than
split out to `BaseBatchSweep.cpp` -- only `compute()`, `_run_range()` and
`_run_one_step()` are ordinary (non-template) members and live in the .cpp.
 **/
template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
class LS2G_API BaseBatchSweep: public BaseBatchSolverSynch
{
    public:
        // per-instance boolean mask matrix type, nested here (not at namespace
        // scope) for the same reason SbusPolicy::Vary::RealMat/CplxMat are nested.
        using BoolMat = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

        explicit BaseBatchSweep(const LSGrid & init_grid_model):
            BaseBatchSolverSynch(init_grid_model),
            ybus_policy_(),
            sbus_policy_()
            {}

        // ContingencyAnalysis-only: constructing with `compute_limit_violations` set
        // up front (mirrors the pre-refactor ContingencyAnalysis(grid, bool) ctor).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        explicit BaseBatchSweep(const LSGrid & init_grid_model, bool compute_limit_violations):
            BaseBatchSweep(init_grid_model)
            { _compute_limit_violations_ = compute_limit_violations; }

        ~BaseBatchSweep() noexcept override = default;
        BaseBatchSweep(const BaseBatchSweep&) = delete;
        BaseBatchSweep(BaseBatchSweep&&) = delete;
        BaseBatchSweep & operator=(BaseBatchSweep&&) = delete;
        BaseBatchSweep & operator=(const BaseBatchSweep&) = delete;

        // name of this instantiation, used in error messages. A 3-way dispatch since
        // the 4th (NOOP,NOOP) combination is never instantiated.
        static const char * algo_name() {
            if(YbusPolicy::supports_contingency && SbusPolicy::supports_vary) return "ScenarioSweep";
            if(YbusPolicy::supports_contingency) return "ContingencyAnalysis";
            return INIT == BatchInitKind::FromPreviousStep ? "TimeSeries" : "InjectionSweep";
        }

        // see BaseBatchSolverSynch::supports_multithread. TimeSeries (INIT ==
        // FromPreviousStep) chains its steps and cannot be split; every other
        // instantiation restarts each step from the same seed (or has no notion of
        // "steps" chaining at all, eg ContingencyAnalysis).
        bool supports_multithread() const override {
            return INIT != BatchInitKind::FromPreviousStep;
        }

        // ================= the three cache levels ================================
        // See BaseBatchSolverSynch's block comment for what each level holds and why
        // dropping one drops every level below it. Each override extends its base with
        // the state THIS class adds, and nothing else: no level reaches sideways into
        // another's.

        // L1 -- what was read off the grid. On top of the solver caches the base
        // drops, the per-contingency Ybus coefficients: they are values read off that
        // matrix, so they do not survive it (a change_algorithm() that switches AC <->
        // DC is exactly this case -- same contingencies, coefficients read off a
        // different matrix).
        void clear_grid_results() override {
            if(_grid_cache_valid_){
                _reset_ybus_coeffs();
                // the union layout the topological actions asked for, and each row's
                // action resolved against the grid (see _maybe_topo_union)
                _topo_used_buses_me_.clear();
                _extra_buses_me_.clear();
                _union_edges_me_.clear();
                _row_topo_.clear();
            }
            BaseBatchSolverSynch::clear_grid_results();
        }

        // L2 -- what was built for this batch from that grid: what the graph walk
        // settled about each contingency, the reference-slack choice and the skip
        // mask that came with it, and the PV <-> PQ structure a generator contingency
        // reserved. The base drops the algorithm, whose sparsity that structure is.
        void clear_batch_inputs() override {
            // not under the test below: the actions are resolved before the "n" solve,
            // so a batch whose "n" case diverged (never valid) holds a resolution too --
            // what the actions were resolved into (the actions themselves are a
            // registration and stay)
            _topo_resolved_ = false;
            _topology_active_ = false;
            _base_masked_.clear();
            _row_emptied_buses_.clear();
            _row_placed_buses_.clear();
            _reset_topo_policy_state();
            if(_batch_inputs_valid_){
                // the workers go with the member algorithm the base drops: each holds a
                // ledger, a Jacobian sparsity and a factorization built from exactly the
                // batch inputs being dropped here (see _compute_threaded)
                _thread_algos_.clear();
                _li_masked.clear();
                _cont_connected_.clear();
                _skip_mask.clear();
                _gen_contingency_active_ = false;
                _switchable_buses_.clear();
                _row_pv_to_pq_.clear();
                _row_pv_pinned_.clear();
                _row_vm_seed_.clear();
                _pv_pinning_active_ = false;
                _row_slack_gens_off_.clear();
                _row_slack_storages_off_.clear();
                _li_defaults_vect_cache_.clear();
                _physical_violations_n_.clear();
                _clear_slack_redistribution();
            }
            BaseBatchSolverSynch::clear_batch_inputs();
        }

        // L3 -- what the last run produced. On top of the voltages and the flows the
        // base drops: convergence, limit violations, and the Jacobians the adjoint
        // kept.
        void clear_batch_outputs() override {
            _adjoint_.clear();
            _adjoint_row_ok_.clear();
            _adjoint_info_.clear();
            _converged.clear();
            _converged_mask_.clear();
            _row_solve_time_.clear();
            _row_nb_iter_.clear();
            _row_lin_stats_.clear();
            _row_J_dim_.clear();
            _row_J_nnz_.clear();
            _outer_results_.clear();
            _gen_results_ = RealMat();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _physical_violations_.clear();
            _gen_p_plan_.clear();
            // NOT _physical_violations_n_: it describes the BASE CASE, which is L2 (see
            // clear_batch_inputs). A compute() that reuses the base case does not re-solve
            // it -- so the report of the solve that built it is the only correct one
            // there, and re-deriving it from a stale algorithm state is exactly the bug
            // this split avoids.
            _timer_modif_Ybus = 0.;
            BaseBatchSolverSynch::clear_batch_outputs();
        }

        // Everything: the three levels, plus what sits on a different axis from them
        // -- the REGISTRATIONS (which contingencies, which injections) and the
        // settings derived from them. A cache level answers "is what I built still
        // right?"; this answers "start over".
        void clear() override {
            BaseBatchSolverSynch::clear();   // -> L1 -> L2 -> L3
            _reset_ybus_policy();
            sbus_policy_.clear();
            topo_actions_.clear();
            _status = 1;
            _nb_steps_locked = -1;
            _handle_disconnected_grid = false;
        }

        // Contingency-varying instantiations only (ContingencyAnalysis AND
        // ScenarioSweep): drop the results, and only the results (mirrors the
        // pre-refactor ContingencyAnalysis::clear_results_only). Kept as the name the
        // python binding has always used; clear_batch_outputs() is the same call.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void clear_results_only() { clear_batch_outputs(); }

        // ================= new setter API (SbusPolicy::supports_vary only) =======
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void modify_gen_p(const Eigen::Ref<const typename S::RealMat> & gen_p) {
            _check_cols(gen_p, static_cast<Eigen::Index>(_grid_model.get_generators_as_data().nb()), "modify_gen_p");
            _lock_or_check_nb_steps(gen_p.rows(), "modify_gen_p");
            sbus_policy_.gen_p = gen_p;
        }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void modify_sgen_p(const Eigen::Ref<const typename S::RealMat> & sgen_p) {
            _check_cols(sgen_p, static_cast<Eigen::Index>(_grid_model.get_static_generators_as_data().nb()), "modify_sgen_p");
            _lock_or_check_nb_steps(sgen_p.rows(), "modify_sgen_p");
            sbus_policy_.sgen_p = sgen_p;
        }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void modify_load_p(const Eigen::Ref<const typename S::RealMat> & load_p) {
            _check_cols(load_p, static_cast<Eigen::Index>(_grid_model.get_loads_as_data().nb()), "modify_load_p");
            _lock_or_check_nb_steps(load_p.rows(), "modify_load_p");
            sbus_policy_.load_p = load_p;
        }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void modify_load_q(const Eigen::Ref<const typename S::RealMat> & load_q) {
            _check_cols(load_q, static_cast<Eigen::Index>(_grid_model.get_loads_as_data().nb()), "modify_load_q");
            _lock_or_check_nb_steps(load_q.rows(), "modify_load_q");
            sbus_policy_.load_q = load_q;
        }
        // Per-step generator target voltage magnitude, in pu (vm_pu), NOT kV. Unlike
        // the four setters above, this does NOT feed the injection (Sbus): it only
        // re-seeds |V| at each voltage-regulating generator's regulated bus before
        // that step's solve (see _apply_step_gen_v / GeneratorContainer::set_vm).
        // Left unset (0 rows), every row keeps using the grid's own target_vm_pu,
        // exactly as before this setter existed.
        //
        // A bus has ONE magnitude, so a row cannot ask one bus for two. That constrains
        // this setter in two ways, and a row that breaks either is not solved -- it
        // reports as a row skipped before the solver (see _row_gen_v_conflicts) rather
        // than silently taking whichever element set_vm happens to visit last:
        //
        //  - several generators regulating the same bus must be given the same target;
        //  - a generator sharing its regulated bus with an SVC or an hvdc converter
        //    station must be given THAT element's target, because this setter cannot
        //    move it. ONLY GENERATOR set-points vary per row: the batch has no
        //    modify_svc_v / modify_hvdc_v, so such a generator's set-point is in
        //    practice fixed for the whole sweep. See the TODO at the top of the
        //    changelog.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void modify_gen_v(const Eigen::Ref<const typename S::RealMat> & gen_v) {
            _check_cols(gen_v, static_cast<Eigen::Index>(_grid_model.get_generators_as_data().nb()), "modify_gen_v");
            _lock_or_check_nb_steps(gen_v.rows(), "modify_gen_v");
            sbus_policy_.gen_v = gen_v;
        }

        // ============ contingency mask setters (Contingency && Vary only) ========
        // ScenarioSweep only: row i, column j True means "deactivate this branch for
        // step i". See YbusPolicy::Contingency::line_mask/trafo_mask.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_lines(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(n_line_), "set_contingency_lines");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_lines");
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();

            ybus_policy_.line_mask = mask;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_trafos(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(n_trafos_), "set_contingency_trafos");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_trafos");
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();

            ybus_policy_.trafo_mask = mask;
        }
        /**
         * ScenarioSweep only: row i, column g True means "disconnect generator g for
         * step i". Unlike the two branch masks above this does NOT touch Ybus -- a
         * generator carries no admittance -- it changes two other things:
         *
         *  - the injection: that step's Sbus loses the generator's active power (and,
         *    if it does not regulate voltage, its reactive setpoint too), and the
         *    distributed slack is re-weighted without it. The lost MW is picked up by
         *    the slack; a caller wanting redispatch says so with modify_gen_p.
         *  - the labelling: when the LAST generator regulating its own bus' voltage is
         *    taken out, that bus stops being PV and becomes PQ for the step. Its
         *    magnitude is then solved for instead of held at the setpoint.
         *
         * The second is what would normally cost a fresh symbolic factorization per
         * row. It does not here: every bus that can flip is given a Vm unknown and a Q
         * equation once, up front, so the Jacobian sparsity is the union over all the
         * steps, and each step merely pins the Q row of the buses that are still PV
         * (see _maybe_prepare_gen_contingency / NRSystem::set_pv_pinned_buses). One
         * analysis, one factorization, refactorizations after that -- as before.
         *
         * Only generators that regulate their OWN bus are supported. A generator that
         * regulates a remote bus, or one whose bus is held by an SVC or an HVDC
         * converter station, is rejected by compute(): those go through the
         * VoltageControl extension, whose own Jacobian rows and columns this does not
         * yet know how to reserve and mask.
         */
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_gens(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(_grid_model.get_generators_as_data().nb()),
                             "set_contingency_gens");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_gens");
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();

            sbus_policy_.gen_off = mask;
        }
        /**
         * ScenarioSweep only: row i, column s True means "disconnect shunt s for step
         * i". A shunt is an admittance: the row loses its diagonal coefficient (AC; in DC
         * its active part leaves Pbus), as a branch mask edits Ybus. Refused by compute()
         * for a shunt that regulates a voltage when the grid's outer loops hold a
         * ShuntVoltageControl one (the shunt's sections are that loop's, which does not
         * know of a row taking it out).
         */
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_shunts(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(_grid_model.get_shunts().nb()), "set_contingency_shunts");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_shunts");
            // a Ybus edit, registered like the branch masks (L2)
            clear_batch_inputs();

            ybus_policy_.shunt_mask = mask;
        }
        /**
         * ScenarioSweep only: row i, column l True means "disconnect load l for step i".
         * An injection only: the row's load power (its own modify_load_p / modify_load_q
         * row, or the grid's) leaves its Sbus; what that changes in the balance is made up
         * by the distributed slack (shared OLF-like before the solve with
         * redistribute_slack). A row naming a load in both this mask and its topological
         * action is refused.
         */
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_loads(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(_grid_model.get_loads_as_data().nb()), "set_contingency_loads");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_loads");
            // an injection, as modify_load_p: rebuilt at every compute(), nothing settled
            // about the batch depends on it
            sbus_policy_.load_off = mask;
        }

        /**
         * ScenarioSweep only: one topological action per row (a `TopoAction`, grid2op
         * semantics: set_bus on loads / generators / storage units / line and trafo
         * ends, set_line_status), played by that row on top of its injections and
         * masks. An action that changes nothing is a plain row.
         *
         * Every action is checked against the grid here (the element and the busbar
         * exist, no contradiction -- see TopoAction::check_validity) and a
         * `std::invalid_argument` naming the row is raised, nothing registered, if one
         * is invalid. What a row plays is settled at compute() (see
         * _maybe_resolve_topology): disconnections -- a branch (set_line_status -1, or
         * set_bus -1 on one of its ends), a generator, a load or a storage unit (set_bus
         * -1), a branch or a generator exactly as the set_contingency_* masks do it --
         * and reconnections and moves between busbars, which create a bus or merge two.
         * The whole sweep keeps one symbolic analysis: the solver labelling holds the
         * union of the buses the rows use, the admittance entries a row writes are
         * reserved up front, and a row is value edits only. A row naming an element in
         * both a mask and its action is refused, and so, for now, are the cases listed
         * at _maybe_resolve_topology (TODO in the changelog).
         */
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_topo_actions(const std::vector<TopoAction> & actions) {
            std::vector<TopoAction> checked;
            checked.reserve(actions.size());
            for(size_t i = 0; i < actions.size(); ++i){
                TopoAction act = actions[i];
                try{
                    act.check_validity(_grid_model);
                }catch(const std::exception & exc_){
                    std::ostringstream msg;
                    msg << algo_name() << "::set_topo_actions: action " << i << " is invalid: " << exc_.what();
                    throw std::invalid_argument(msg.str());
                }
                checked.push_back(act);
            }
            _lock_or_check_nb_steps(static_cast<Eigen::Index>(actions.size()), "set_topo_actions");
            // a registration that reaches L1, unlike the masks: the buses the actions
            // use are part of the solver labelling (see _maybe_topo_union), so what was
            // built for the OLD actions goes, the grid cache included
            clear_grid_results();
            topo_actions_ = checked;
        }
        // the (checked) actions registered with set_topo_actions
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        const std::vector<TopoAction> & get_topo_actions() const { return topo_actions_; }

        // the branches (gridmodel numbering, powerlines then transformers, sorted) row
        // `row` disconnects, by its masks and by its topological action together --
        // what the row's flows are zeroed for and its current checks skip. Needs a
        // compute() (the actions are resolved there).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        std::vector<int> get_row_disconnected_branches(Eigen::Index row) const {
            if(row < 0 || row >= _nb_steps_locked){
                std::ostringstream exc_;
                exc_ << algo_name() << "::get_row_disconnected_branches: row " << row
                     << " does not exist, the batch has " << _nb_steps_locked << " rows.";
                throw std::out_of_range(exc_.str());
            }
            // resolved, not converged: a batch whose "n" case diverged still knows what
            // each row disconnects (ScenarioSweep.run() reports it for that case too)
            if(!topo_actions_.empty() && !_topo_resolved_){
                std::ostringstream exc_;
                exc_ << algo_name() << "::get_row_disconnected_branches: the topological actions "
                        "are resolved by compute(); call it first.";
                throw std::runtime_error(exc_.str());
            }
            return ybus_policy_.branch_ids_for_row(row, n_line_);
        }

        // ========== legacy contingency-list API (Contingency && !Vary only) ======
        // ContingencyAnalysis only -- deliberately NOT unified with
        // set_contingency_lines/trafos above (two genuinely different usages -- a
        // shared base case with an arbitrary set of distinct scenarios here, vs. one
        // contingency per row paired with that row's own injection on ScenarioSweep).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_all_n1(){
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            for(int l_id = 0; l_id < static_cast<int>(n_total_); ++l_id){
                ybus_policy_.li_defaults.insert(_branch_contingency({l_id}));
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_n1(int line_id){
            _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(_branch_contingency({line_id}));
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s) _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            for(const auto line_id : vect_n1s){
                ybus_policy_.li_defaults.insert(_branch_contingency({line_id}));
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_nk(const std::vector<int> & vect_nk){
            for(const auto line_id : vect_nk) _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(_branch_contingency(vect_nk));
        }
        /**
         * ContingencyAnalysis only: one contingency tripping, together, the branches
         * (lines then transformers), the shunts, the loads and the generators given (any
         * of the four may be empty, not all). A contingency registers once: adding one
         * already registered changes nothing.
         *
         *  - a branch is taken out of Ybus, as add_nk does;
         *  - a shunt too: it is an admittance, one diagonal coefficient (AC; in DC its
         *    active part leaves Pbus). Refused by compute() when the shunt regulates a
         *    voltage and the grid's outer loops hold a ShuntVoltageControl one (the
         *    shunt's sections are that loop's, which does not know of a row taking it out);
         *  - a load leaves the row's injection, its active and reactive power;
         *  - a generator leaves it too (its active power, and the reactive setpoint of
         *    one that does not regulate voltage), and the distributed slack is
         *    re-weighted without it. When it is the last one regulating its own bus,
         *    that bus turns PQ for the row, as set_contingency_gens does it on a
         *    ScenarioSweep (the same machinery, see _maybe_prepare_gen_contingency, and
         *    the same scope: a generator regulating a remote bus, or on a bus a control
         *    group holds, is refused by compute()).
         *  - a branch of `branch_ids_side1` (`branch_ids_side2`) is opened at side 1 (2)
         *    only: its other end stays connected and keeps the Kron-reduced part of the
         *    branch (its charging), as a branch open at one end in the grid does. That
         *    is OpenLoadFlow's model of a branch whose breaker opens at one end only.
         *    For the flows, the current checks and the outer loops it counts as open.
         *    Named at both ends, or also in `branch_ids`, it is opened at both.
         * What it loses is made up by the distributed slack (redistribute_slack shares
         * it OLF-like before the solve), and a bus it strands is handled as a branch
         * contingency strands one (handle_disconnected_grid).
         */
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_contingency(const std::vector<int> & branch_ids,
                             const std::vector<int> & shunt_ids,
                             const std::vector<int> & load_ids,
                             const std::vector<int> & gen_ids,
                             const std::vector<int> & branch_ids_side1 = std::vector<int>(),
                             const std::vector<int> & branch_ids_side2 = std::vector<int>()){
            const ContingencyElements cont = _checked_contingency(branch_ids, shunt_ids, load_ids, gen_ids, "add_contingency",
                                                                  branch_ids_side1, branch_ids_side2);
            // a contingency registration (L2), as add_nk
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(cont);
        }
        // shortcuts: one contingency tripping these shunts / loads / generators only
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_contingency_shunt(const std::vector<int> & shunt_ids){
            const ContingencyElements cont = _checked_contingency({}, shunt_ids, {}, {}, "add_contingency_shunt");
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(cont);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_contingency_load(const std::vector<int> & load_ids){
            const ContingencyElements cont = _checked_contingency({}, {}, load_ids, {}, "add_contingency_load");
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(cont);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_contingency_gen(const std::vector<int> & gen_ids){
            const ContingencyElements cont = _checked_contingency({}, {}, {}, gen_ids, "add_contingency_gen");
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(cont);
        }
        // removes the contingency tripping exactly these elements; whether there was one
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_contingency(const std::vector<int> & branch_ids,
                                const std::vector<int> & shunt_ids,
                                const std::vector<int> & load_ids,
                                const std::vector<int> & gen_ids,
                                const std::vector<int> & branch_ids_side1 = std::vector<int>(),
                                const std::vector<int> & branch_ids_side2 = std::vector<int>()){
            const ContingencyElements cont = _checked_contingency(branch_ids, shunt_ids, load_ids, gen_ids, "remove_contingency",
                                                                  branch_ids_side1, branch_ids_side2);
            clear_batch_inputs();
            return ybus_policy_.li_defaults.erase(cont) >= 1;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_n1(int line_id){
            _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            auto nb_removed = ybus_policy_.li_defaults.erase(_branch_contingency({line_id}));
            return nb_removed >= 1;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        size_t remove_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s) _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            size_t nb_removed = 0;
            for(const auto line_id : vect_n1s){
                nb_removed += ybus_policy_.li_defaults.erase(_branch_contingency({line_id}));
            }
            return nb_removed;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_nk(const std::vector<int> & vect_nk){
            for(const auto line_id : vect_nk) _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            auto nb_removed = ybus_policy_.li_defaults.erase(_branch_contingency(vect_nk));
            return nb_removed >= 1;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::set<ContingencyElements> & my_defaults() const {return ybus_policy_.li_defaults;}
        // the branches each contingency trips, at both ends, in row order (a contingency
        // tripping only shunts, loads or generators, or opening branches at one end, lists
        // none of those: see my_contingencies)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<std::vector<int> > my_defaults_vect() const {
            std::vector<std::vector<int> > res;
            res.reserve(ybus_policy_.li_defaults.size());
            for(const auto & cont: ybus_policy_.li_defaults){
                res.push_back(std::vector<int>(cont.branches.begin(), cont.branches.end()));
            }
            return res;
        }
        // every branch each contingency opens, at both ends or at one, in row order: what
        // the flows, the current checks and the outer loops see as open (see add_contingency)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<std::vector<int> > _branches_out_vect() const {
            std::vector<std::vector<int> > res;
            res.reserve(ybus_policy_.li_defaults.size());
            for(const auto & cont: ybus_policy_.li_defaults) res.push_back(cont.branches_out());
            return res;
        }
        // everything each contingency trips, in row order: (branches, shunts, loads,
        // generators, branches opened at side 1 only, at side 2 only), each sorted
        using ContingencyTuple = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>,
                                            std::vector<int>, std::vector<int> >;
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<ContingencyTuple> my_contingencies() const {
            std::vector<ContingencyTuple> res;
            res.reserve(ybus_policy_.li_defaults.size());
            for(const auto & cont: ybus_policy_.li_defaults){
                res.push_back(ContingencyTuple(std::vector<int>(cont.branches.begin(), cont.branches.end()),
                                               std::vector<int>(cont.shunts.begin(), cont.shunts.end()),
                                               std::vector<int>(cont.loads.begin(), cont.loads.end()),
                                               std::vector<int>(cont.gens.begin(), cont.gens.end()),
                                               std::vector<int>(cont.branches_side1.begin(), cont.branches_side1.end()),
                                               std::vector<int>(cont.branches_side2.begin(), cont.branches_side2.end())));
            }
            return res;
        }

        // "handle disconnected grid" mode (ContingencyAnalysis AND ScenarioSweep --
        // both edit Ybus per step, so a contingency can strand a component either
        // way; TimeSeries/InjectionSweep have no notion of a contingency at all)
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool get_handle_disconnected_grid() const {return _handle_disconnected_grid;}
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void set_handle_disconnected_grid(bool val) {
            if(val == _handle_disconnected_grid) return;
            // L2: this decides whether a contingency that strands a component is
            // masked out or skipped outright, and it is _maybe_prepare_masks() that
            // settles the reference slack and the skip mask from it.
            clear_batch_inputs();
            _handle_disconnected_grid = val;
        }

        // OLF-style bounded redistribution of the active power a row loses (a generator
        // contingency, an island cut off in "handle disconnected grid" mode) BEFORE its
        // solve, the saturated units leaving that row's distributed slack -- see
        // _prepare_slack_redistribution and SlackRedistribution.hpp. Off by default.
        // Not L2: the pre-pass is rebuilt at every compute() (it depends on the
        // per-row injections, which a second compute() comes to change).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool get_redistribute_slack() const {return _redistribute_slack_;}
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void set_redistribute_slack(bool val) {_redistribute_slack_ = val;}

        // limit violations (ContingencyAnalysis AND ScenarioSweep -- see
        // get_violations()/get_violations_n() below; deliberately NO
        // converged()/converged_n() equivalent here: a non-converged row's
        // get_violations() entry already carries a GRID/NOT_SIMULATED-or-DIVERGENCE
        // sentinel LimitViolation, so a separate convergence flag would be
        // redundant. converged()/converged_n() stay ContingencyAnalysis-only,
        // unchanged, below.)
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool get_compute_limit_violations() const noexcept {return _compute_limit_violations_;}
        // What this flag invalidates is only the RESULTS (L3): whether violations are
        // recorded changes nothing about the grid, the Jacobian or the "n" solve. It
        // nonetheless clear()s the whole object, registrations included -- a contract
        // that predates the cache levels and is explicitly tested (see
        // test_ContingencyAnalysis_limit_violations.py's
        // test_setter_toggle_clears_and_raises_when_off), so it stays. The established convention is therefore to set
        // this flag FIRST, before handle_disconnected_grid / registering any
        // contingency or injection -- ContingencyAnalysis's 2-arg constructor enforces
        // this implicitly; ScenarioSweep has no such constructor and must follow the
        // same order explicitly (see its own tests / scenarioSweep.py).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void set_compute_limit_violations(bool val){
            if(val == _compute_limit_violations_) return;
            _compute_limit_violations_ = val;
            clear();
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        real_type get_violation_threshold() const noexcept {return _violation_threshold_;}
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void set_violation_threshold(real_type val){
            if(!(val > 0. && val <= 1.)){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_violation_threshold: the threshold should be "
                        "a real number in the range ]0., 1.] (got " << val << ").";
                throw std::runtime_error(exc_.str());
            }
            // L3: a tighter threshold reports violations the recorded ones do not
            // contain, so those have to be recomputed. Nothing else changes -- the
            // same grid, the same batch, the same solve.
            if(val < _violation_threshold_) clear_batch_outputs();
            _violation_threshold_ = val;
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        real_type get_violation_rel_tol() const noexcept {return _violation_rel_tol_;}
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void set_violation_rel_tol(real_type val){
            const std::string fun_name = std::string(algo_name()) + "::set_violation_rel_tol";
            batch_sweep_detail::check_violation_rel_tol(val, fun_name.c_str());
            if(val == _violation_rel_tol_) return;
            // L3, like set_violation_threshold -- but either direction: a larger tolerance
            // drops recorded violations as surely as a smaller one adds some
            clear_batch_outputs();
            _violation_rel_tol_ = val;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::vector<char> & converged() const {
            _check_limit_violations_enabled("converged");
            return _converged;
        }
        // per-row converged mask, available on EVERY instantiation (TimeSeries /
        // InjectionSweep / ContingencyAnalysis / ScenarioSweep) unconditionally --
        // unlike converged() above (ContingencyAnalysis-only, and gated behind
        // compute_limit_violations=True), this is always sized to nb_steps() and
        // always populated by compute(), no setup required. Row i is 1 iff that row
        // was both invertible and reported convergence by the solver; a row an
        // aborted TimeSeries chain never reaches (see _run_range's
        // BatchInitKind::FromPreviousStep early return) reads as 0, same as an
        // outright divergence -- indistinguishable from the caller's point of view,
        // which matches nb_solved()/nb_converged()'s own "never attempted" == "did
        // not converge" convention.
        const std::vector<char> & converged_mask() const { return _converged_mask_; }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        const std::vector<std::vector<LimitViolation> > & get_violations() const {
            _check_limit_violations_enabled("get_violations");
            return _violations;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool converged_n() const {
            _check_limit_violations_enabled("converged_n");
            return _converged_n_;
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        const std::vector<LimitViolation> & get_violations_n() const {
            _check_limit_violations_enabled("get_violations_n");
            return _violations_n_;
        }

        // ---- physical-limit checks (EVERY instantiation) ---------------------------
        // Opt in to the checks whose violation says the converged row is not a state the
        // grid can reach at all -- ViolationCategory::PHYSICAL -- or that it leaves a control
        // an outer loop would apply -- ViolationCategory::CONTROL --, as opposed to the
        // operational limits `compute_limit_violations` reports (a voltage band, a thermal
        // rating: states the grid does reach and should not sit in). Each one a
        // condition an OpenLoadFlow outer loop acts on -- reported, and only enforced when the
        // algorithm runs the loops (NROuter_*, after which they have nothing left to report):
        //
        //   * the REACTIVE CAPABILITY of each bus whose voltage is held by machines
        //     (LOW_Q / HIGH_Q on the BUS): did it need more reactive power than the sum of
        //     what its voltage-regulating generators, hvdc converter stations and
        //     voltage-mode SVCs can produce? `ReactiveLimits`' own trigger;
        //   * the RELEASE of each PQ generator the caller flagged as pinned at a reactive
        //     limit (LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q on the GENERATOR, see
        //     GenPvReleaseCheck.hpp and LSGrid::set_gen_can_be_pv): does the bus it would
        //     regulate sit below its target while the machine absorbs all it can (or above
        //     while it produces all it can)? The other direction of the same
        //     `ReactiveLimits` loop, PQ -> PV;
        //   * the SWITCH ON of each idle SVC the caller flagged as carrying a standby
        //     automaton (LOW_VOLTAGE_SVC_STANDBY / HIGH_VOLTAGE_SVC_STANDBY on the SVC, see
        //     SvcStandbyCheck.hpp and LSGrid::set_svc_standby): does the bus it regulates
        //     sit outside the automaton's voltage thresholds? OpenLoadFlow's
        //     `MonitoringVoltageOuterLoop`;
        //   * the OWN BUS of each generator regulating a remote bus
        //     (LOW_VOLTAGE_REMOTE_CONTROL / HIGH_VOLTAGE_REMOTE_CONTROL on the GENERATOR):
        //     does holding the remote target take it outside the realistic range?
        //     `ReactiveLimits`' trigger in its robust remote voltage control mode;
        //   * the ACTIVE POWER of each angle-droop ("AC emulation") hvdc line still in the
        //     linear regime (HIGH_P on the HVDC, see HvdcPCheck.hpp): did it transmit more
        //     than `pmax_1to2_mw` / `pmax_2to1_mw` allow in that direction? OpenLoadFlow's
        //     `HvdcAcEmulationLimits`;
        //   * the RELEASE of each angle-droop hvdc line the caller flagged as frozen at its
        //     active power limit (HVDC_AC_EMULATION_RELEASE on the HVDC, see HvdcPCheck.hpp and
        //     LSGrid::set_hvdc_ac_emulation_frozen): would its droop ask for less than that
        //     limit? The other half of the same `HvdcAcEmulationLimits` loop;
        //   * the ACTIVE POWER of each generator and each storage unit carrying the
        //     DISTRIBUTED SLACK (LOW_P / HIGH_P on the GENERATOR / STORAGE, see
        //     GenPCheck.hpp): the slack is solved inside the Jacobian by participation
        //     factors that know nothing about limits, so `target_p + share` can land beyond
        //     `min_p_mw` / `max_p_mw`. OpenLoadFlow's `DistributedSlack`, whose whole job is
        //     to take a saturated machine out of the distribution and re-share what is left.
        //
        // Available on all four instantiations (unlike compute_limit_violations, which is
        // contingency-only): a time series that walks a load curve is exactly as likely to
        // ask a bus for reactive power it does not have, or a machine for active power it
        // cannot deliver, as a contingency is.
        //
        // WHAT EACH CHECK NEEDS. The hvdc ones need only the bus angles and the active-power
        // one only the slack the row distributed -- both of which every algorithm leaves
        // behind, AC and DC alike. The reactive one needs an AC algorithm that publishes
        // its per-bus mismatch (BaseAlgo::fills_bus_mismatch -- every built-in AC family
        // does; a plugin solver has to opt in), and compute() raises for one that does not
        // rather than reporting nothing. In DC it is simply not applicable: a DC powerflow
        // has no reactive power at all, so a DC batch reports the active-power checks
        // and nothing is hidden by it. The release check needs a voltage magnitude, so it
        // is AC only too, and the flags themselves (LSGrid::set_gen_can_be_pv): a grid with
        // none reports nothing there. Same for the standby SVC check (voltage magnitudes,
        // LSGrid::set_svc_standby) and the remote voltage control one (voltage magnitudes).
        // The active-power check also needs the limits
        // themselves, which are optional (LSGrid::set_gen_p_limits /
        // set_storage_p_limits): a grid that has none simply reports nothing there.
        //
        // Setting this drops this batch's base case and results, but NOT its registrations
        // (the contingencies / injections), so unlike compute_limit_violations -- which
        // clear()s the whole object -- it can be set at any point before compute().
        // ---- generator results per row (opt-in) ----------------------------------------
        // Every converged row records the active and reactive output of every generator,
        // by the rules LSGrid::compute_results publishes them with after a single solve
        // (LSGrid::generator_results: the setpoint plus the share of the distributed slack,
        // the reactive output the algorithm solved or the bus' residual split by reactive
        // range). A generator the row disconnects, or strands on a bus it masks, reads 0.
        // Off by default: it reads the row's mismatch and allocates per row. Takes effect
        // at the next compute(); changes nothing a kept base case is made of.
        bool get_compute_gen_results() const noexcept {return _compute_gen_results_;}
        void set_compute_gen_results(bool val){ _compute_gen_results_ = val; }
        // (nb_rows, 2 * nb_gen), row-major: generator g of row i at columns 2g (P, MW) and
        // 2g + 1 (Q, MVAr). Empty if the last compute() did not record them.
        const RealMat & get_gen_results() const { return _gen_results_; }

        bool get_compute_physical_violations() const noexcept {return _compute_physical_violations_;}
        void set_compute_physical_violations(bool val){
            if(val == _compute_physical_violations_) return;
            _compute_physical_violations_ = val;
            // L2, not L3: the base case's own report (get_physical_violations_n) can only be
            // read off the solve that built it, and a kept base case is not re-solved. So
            // a change of what is reported has to drop that base case with the results.
            clear_batch_inputs();
        }
        // Absolute slack on every comparison these checks make, so that an element resting
        // exactly on its limit is not reported over solver noise: a violation needs
        // `value > limit + tol` (or, for LOW_Q, `value < limit - tol`). Defaults to 1e-4,
        // well above what a converged solve leaves behind and well below anything
        // meaningful. In MVA: one noise floor for both halves, since MW and MVAr are the
        // same scale and this is a numerical threshold rather than a physical parameter.
        real_type get_physical_violation_tol_mva() const noexcept {return _physical_tol_mva_;}
        void set_physical_violation_tol_mva(real_type val){
            // isfinite unqualified, like everywhere else in this header (see the note on
            // the <math.h> include at the top)
            if(!(val >= 0.) || !isfinite(val)){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_physical_violation_tol_mva: the tolerance should "
                        "be a finite, non-negative number of MVA (got " << val << ").";
                throw std::runtime_error(exc_.str());
            }
            // any change, in either direction: a smaller tolerance reports violations the
            // recorded ones do not contain, a larger one leaves recorded ones that should
            // no longer be reported. L2 for the same reason as the flag above.
            if(val != _physical_tol_mva_) clear_batch_inputs();
            _physical_tol_mva_ = val;
        }
        // The same, for the comparisons made on a voltage: the PQ -> PV release check
        // (GenPvReleaseCheck.hpp) reports a flagged machine whose regulated voltage is
        // below (at min_q) or above (at max_q) its target by more than this, in pu, and the
        // standby SVC check (SvcStandbyCheck.hpp) a flagged SVC whose regulated voltage is
        // outside its automaton's thresholds by more than this. A separate knob because a
        // voltage and a power are not the same scale.
        real_type get_physical_violation_tol_vm_pu() const noexcept {return _physical_tol_vm_pu_;}
        void set_physical_violation_tol_vm_pu(real_type val){
            if(!(val >= 0.) || !isfinite(val)){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_physical_violation_tol_vm_pu: the tolerance should "
                        "be a finite, non-negative number of pu (got " << val << ").";
                throw std::runtime_error(exc_.str());
            }
            if(val != _physical_tol_vm_pu_) clear_batch_inputs();
            _physical_tol_vm_pu_ = val;
        }
        /**
         * Per row: the physical limits this row's solution leaves. A row that did not
         * converge (or that was never simulated) has an EMPTY entry rather than a sentinel
         * -- ask converged_mask() to tell that apart from "converged, no violation". Every
         * entry has category PHYSICAL, and one of these shapes:
         *
         *   - element_type BUS, element_id the grid bus id, violation_type LOW_Q / HIGH_Q,
         *     `value` the reactive power the machines holding that bus had to produce
         *     (MVAr) and `limit` their summed capability;
         *   - element_type GENERATOR, element_id the generator id, violation_type
         *     LOW_VOLTAGE_AT_MIN_Q / HIGH_VOLTAGE_AT_MAX_Q, `value` the voltage of the bus
         *     that flagged PQ machine would regulate and `limit` its target, both in kV;
         *   - element_type SVC, element_id the svc id, violation_type
         *     LOW_VOLTAGE_SVC_STANDBY / HIGH_VOLTAGE_SVC_STANDBY, `value` the voltage of the
         *     bus that flagged idle standby SVC regulates and `limit` the automaton's
         *     threshold, both in kV;
         *   - element_type GENERATOR, element_id the generator id, violation_type
         *     LOW_VOLTAGE_REMOTE_CONTROL / HIGH_VOLTAGE_REMOTE_CONTROL, `value` the voltage of
         *     that remote controller's own bus and `limit` the realistic bound, both in kV;
         *   - element_type HVDC, element_id the hvdc line id, violation_type HIGH_P, `side`
         *     the direction (1 for 1 -> 2), `value` the active power leaving that side (MW,
         *     positive) and `limit` that direction's pmax;
         *   - element_type HVDC, element_id the hvdc line id, violation_type
         *     HVDC_AC_EMULATION_RELEASE, `side` the direction it is frozen in, `value` the
         *     flow its droop asks for in that direction and `limit` the limit it is frozen
         *     at, both in MW;
         *   - element_type GENERATOR, element_id the generator id, violation_type LOW_P /
         *     HIGH_P, `value` its converged active power (MW, target plus its share of the
         *     distributed slack) and `limit` its min_p_mw / max_p_mw.
         *
         * Requires compute_physical_violations = true.
         */
        const std::vector<std::vector<LimitViolation> > & get_physical_violations() const {
            _check_physical_violations_enabled("get_physical_violations");
            return _physical_violations_;
        }
        /// The same, for the base ("n") case every row is solved from (no injection
        /// change, no contingency). Empty if that solve did not converge.
        const std::vector<LimitViolation> & get_physical_violations_n() const {
            _check_physical_violations_enabled("get_physical_violations_n");
            return _physical_violations_n_;
        }

        // contingency-list inspection (ContingencyAnalysis only). Public (directly
        // bound to Python) -- see the class-level note on why this needs a fully
        // inline body, including everything it transitively calls.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        IntVect is_grid_connected_after_contingency(){
            const bool ac_solver_used = _algo.ac_solver_used();
            if(!_grid_cache_valid_ || ybus_policy_.li_coeffs.size() != ybus_policy_.li_defaults.size()){
                const size_t nb_total_bus = _grid_model.total_bus();
                CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                                    {_grid_model.get_init_vm_pu(), 0.});
                prepare_solver_input_base(Vinit, ac_solver_used);
                ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, active_layout().id_me_to_solver, n_line_);
            }
            _prepare_connectivity();
            IntVect res = IntVect::Constant(static_cast<Eigen::Index>(_cont_connected_.size()), 0);
            for(size_t cont_id = 0; cont_id < _cont_connected_.size(); ++cont_id){
                res(static_cast<Eigen::Index>(cont_id)) = _cont_connected_[cont_id] ? 1 : 0;
            }
            return res;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        int pick_reference_slack(){
            // this reorders the slack, which every row of the next batch is solved
            // against: L2 state, changed here outside compute()
            clear_batch_inputs();
            const bool ac_solver_used = _algo.ac_solver_used();
            if(!_grid_cache_valid_ || ybus_policy_.li_coeffs.size() != ybus_policy_.li_defaults.size()){
                const size_t nb_total_bus = _grid_model.total_bus();
                CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                                    {_grid_model.get_init_vm_pu(), 0.});
                prepare_solver_input_base(Vinit, ac_solver_used);
                ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, active_layout().id_me_to_solver, n_line_);
            }
            _prepare_connectivity();
            // a suggestion only: the automatic choice, whatever reference the grid
            // forces, and the slack order the next compute() starts from is left as
            // it was (compute() makes its own choice, see _maybe_prepare_masks)
            SolverBusLayout & layout = active_layout();
            const GlobalBusIdVect slack_me = layout.slack_bus_id_me;
            const SolverBusIdVect slack_solver = layout.slack_bus_id_solver;
            _select_ref_slack_and_masks(false);
            const int res = (layout.slack_bus_id_me.size() == 0) ? -1 : layout.slack_bus_id_me[0].cast_int();
            layout.slack_bus_id_me = slack_me;
            layout.slack_bus_id_solver = slack_solver;
            return res;
        }

        // ================= legacy bundled compute_Vs / get_sbuses ================
        // TimeSeries/InjectionSweep only: NOT available on ScenarioSweep (it never
        // had this call), hence the extra `&& !Y::supports_contingency` on top of
        // `S::supports_vary`.
        template<class S = SbusPolicy, class Y = YbusPolicy,
                 typename std::enable_if<S::supports_vary && !Y::supports_contingency, int>::type = 0>
        int compute_Vs(const Eigen::Ref<const typename S::RealMat> & gen_p,
                       const Eigen::Ref<const typename S::RealMat> & sgen_p,
                       const Eigen::Ref<const typename S::RealMat> & load_p,
                       const Eigen::Ref<const typename S::RealMat> & load_q,
                       const Eigen::Ref<const CplxVect> & Vinit,
                       const int max_iter,
                       const real_type tol)
        {
            modify_gen_p(gen_p);
            modify_sgen_p(sgen_p);
            modify_load_p(load_p);
            modify_load_q(load_q);
            try {
                compute(Vinit, max_iter, tol);
            } catch(const std::exception &) {
                // base-case ("n") non-convergence: compute() throws; this legacy
                // wrapper instead preserves its own historical -1 return convention.
                _status = 0;
                return -1;
            }
            return _status;
        }
        // the whole nb_steps x nb_bus injection matrix, built on request: the row
        // loop itself never holds it (see SbusPolicy::Vary::fill_row)
        template<class S = SbusPolicy, class Y = YbusPolicy,
                 typename std::enable_if<S::supports_vary && !Y::supports_contingency, int>::type = 0>
        Eigen::Ref<const typename S::CplxMat> get_sbuses() const {return sbus_policy_.materialize();}

        // aggregate status (TimeSeries/InjectionSweep/ScenarioSweep -- everything
        // with SbusPolicy::supports_vary; ContingencyAnalysis reports per-row via
        // converged() instead, it has no aggregate status).
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        int get_status() const { return _status; }

        // ================= base-case reuse ========================================
        // Every compute() has a BASE CASE to establish before any row is solved: read
        // the grid (admittance matrix, bus labelling, pv/pq split), walk the graph to
        // settle what each contingency strands, solve one "n" powerflow, and analyze
        // and factorize the Jacobian. None of it depends on the injections, so a
        // second compute() on the same object repeats all of it for nothing -- which
        // is exactly what a caller in a loop does, and what training a model on a
        // batch does thousands of times.
        //
        // That base case is levels L1 and L2 (see BaseBatchSolverSynch's block comment
        // on clear_grid_results / clear_batch_inputs / clear_batch_outputs): it is kept
        // between calls, and rebuilt exactly when a modifier drops the level it belongs
        // to. Every modifier of this class names its level, so there is no separate
        // list to keep in sync. The grid cannot change underneath: this class holds its
        // OWN copy of it, taken at construction, and offers no way to modify it.
        //
        // That includes the WORKER algorithms of a multi-threaded batch. The base case
        // lives on the member algorithm, which is what the single-threaded path runs the
        // rows with -- so at nb_thread == 1 keeping the base case already meant keeping
        // the factorization. At nb_thread > 1 the rows are run by one algorithm per
        // thread instead, and those were rebuilt every call: a fresh AlgorithmSelector
        // each, and a first row under tell_all_changed(), so build_J_sparsity + analyze
        // + factorize rather than a refactorize -- the one thing a batch exists to pay
        // once. They are kept on exactly the same terms now.
        //
        // ON by default. Turn it off to have every compute() rebuild everything, as it
        // did before this existed -- worth doing to tell a suspected caching bug from a
        // real one, and the reason this is a public switch rather than an internal
        // detail.
        void set_reuse_base_case(bool val) {
            if(val == _reuse_base_case_) return;
            _reuse_base_case_ = val;
            // L1: with reuse off, nothing kept may be believed -- and compute() starts
            // from a clean slate every time while it stays off.
            if(!val) clear_grid_results();
        }
        bool get_reuse_base_case() const { return _reuse_base_case_; }

        // Throw the kept base case away, so the next compute() builds a fresh one.
        // A synonym for clear_grid_results(), the top of the hierarchy. This class
        // calls the right level itself whenever it changes anything a kept base case
        // is made of, so there is normally nothing to do -- this is public for the one
        // case this class cannot see: a grid handed to it that later grows a way of
        // being modified in place (see the change_gridmodel TODO in
        // BaseBatchSolverSynch).
        void invalidate_base_case() { clear_grid_results(); }

        // Whether the last compute() kept a base case instead of building one. For
        // tests, and for anyone measuring where a batch's time goes.
        bool base_case_was_reused() const { return _base_case_was_reused_; }

        // Whether the last compute() also kept the WORKER algorithms of the
        // multi-threaded path, rather than building and analyzing one per thread. False
        // for a single-threaded batch, which has no workers (the member algorithm does
        // the rows itself, and keeping that is base_case_was_reused above).
        bool thread_algos_were_reused() const { return _thread_algos_were_reused_; }

        // ================= reverse-mode differentiation ============================
        // See BatchAdjoint for the maths and the cost. In short: with this on, every
        // row's converged Jacobian is kept, and solve_JT() then answers
        // `J_i^T lambda_i = xbar_i` for the whole batch at one transposed solve per
        // row -- which is the gradient of a scalar loss with respect to every
        // injection of every row.

        // Keep each row's converged Jacobian during compute(), so solve_JT() can run
        // afterwards. Off by default: it costs nb_rows * nnz_J reals (see
        // adjoint_memory_bytes(), answerable before compute() through
        // BatchAdjoint::memory_bytes) plus one Jacobian fill per row, and a caller
        // who only wants flows should pay neither. Only the Newton-Raphson family of
        // AC algorithms has a Jacobian to keep; compute() raises if this is on with
        // any other.
        void set_keep_jacobian(bool val) {
            if(val == _keep_jacobian_) return;
            _keep_jacobian_ = val;
            // L3: the kept Jacobians are an output, and nothing else this flag touches
            // outlives a compute().
            clear_batch_outputs();
        }
        bool get_keep_jacobian() const { return _keep_jacobian_; }

        // Bytes the kept Jacobians occupy after a compute() (0 when none were kept).
        std::size_t adjoint_memory_bytes() const { return _adjoint_.memory_bytes(); }

        // Dimension of the augmented Newton-Raphson system: the length of one
        // cotangent, and of one lambda. 0 before a compute() that kept the Jacobians.
        int dim_J() const { return static_cast<int>(_adjoint_.dim_J()); }

        // Where each GRID bus sits in the Jacobian, -1 where it owns no such
        // row / column. Keyed by grid bus id -- the same numbering as the columns of
        // get_voltages() -- rather than by the solver's own, which is an internal
        // labelling a caller has no business reconstructing:
        //   *_col_of_bus: the unknown (a voltage angle / magnitude) -> where a
        //                 cotangent of that bus goes in xbar;
        //   *_row_of_bus: the equation (an active / reactive mismatch) -> where that
        //                 bus's injection gradient is read out of lambda.
        IntVect get_theta_col_of_bus() const { return _bus_map_to_grid(_algo.get_theta_to_J_col_python()); }
        IntVect get_vm_col_of_bus() const { return _bus_map_to_grid(_algo.get_vm_to_J_col_python()); }
        IntVect get_p_row_of_bus() const { return _bus_map_to_grid(_algo.get_p_to_J_row_python()); }
        IntVect get_q_row_of_bus() const { return _bus_map_to_grid(_algo.get_q_to_J_row_python()); }

        // Solve the adjoint system of every row. `xbar` is (nb_rows, k * dim_J): one
        // row per batch row, holding k cotangents of dim_J coefficients laid end to
        // end (k = 1 for the gradient of a scalar loss). Returns lambda, same shape.
        // Rows that did not converge come back zero -- see adjoint_row_ok().
        BatchAdjoint::RealMatRM solve_JT(const Eigen::Ref<const BatchAdjoint::RealMatRM> & xbar) {
            if(!_adjoint_.is_allocated()){
                std::ostringstream exc_;
                exc_ << algo_name() << "::solve_JT: no Jacobian was kept for this batch. Set "
                        "`keep_jacobian` to True BEFORE calling compute().";
                throw std::runtime_error(exc_.str());
            }
            return _adjoint_.solve_JT(xbar, _adjoint_identity_rows(), _nb_thread, _adjoint_row_ok_);
        }

        // ---- the generator voltage setpoint (modify_gen_v) ------------------------
        // A bus whose magnitude the Newton-Raphson does NOT solve for -- a PV bus, the
        // slack -- holds it at whatever the starting voltage had, and modify_gen_v is
        // what puts it there (GeneratorContainer::set_vm). So for those buses gen_v is
        // a PARAMETER of the system, and a loss differentiates through it in two parts:
        //
        //   direct    L depends on V_k = v_k . e^{j.theta_k} explicitly, at fixed
        //             unknowns: Re(conj(V_k/|V_k|) . dL/dV_k). Note this is the very
        //             same quantity a PQ bus contributes to xbar -- a bus whose |V| is
        //             an unknown hands it to the adjoint, a bus whose |V| is a
        //             parameter hands it to that parameter's gradient.
        //   indirect  v_k moves the solution: -lambda^T dF/dv_k, which is
        //             gen_v_indirect_grad() below. dF/dv_k is the dS/d|V| column the
        //             Jacobian does not store for such a bus, precisely because its
        //             magnitude is not an unknown.
        //
        // The caller does the direct half (it is two lines of the same arithmetic that
        // builds xbar) and asks for the indirect one here, where Ybus is.

        // The GRID bus whose magnitude each generator's gen_v set-point helps fix, -1
        // where it fixes none -- which is the case for a generator that is disconnected,
        // not regulating or treated as off; for one regulating a bus a voltage-control
        // group holds (the set-point is then that group's v_set, not a fixed |V| -- see
        // get_gen_v_vc_row); and for one whose bus is pinned by an SVC or an hvdc
        // converter station, whose set-point no batch can move, so the generator's is not
        // free either. A bus the solve gives a magnitude unknown (a switchable bus) is a
        // target too: its set-point reaches a row while the bus is pinned PV there, and
        // nothing while it is PQ -- which get_gen_v_share_rows says, row by row.
        IntVect get_gen_v_target_bus() const {
            const power_flow_adjoint::GenVLayout layout = _gen_v_layout();
            const auto solver_to_me = active_layout().id_solver_to_me.as_eigen();
            IntVect res = IntVect::Constant(layout.target_bus.size(), -1);
            for(Eigen::Index g = 0; g < res.size(); ++g){
                const int b = layout.target_bus[g];
                if(b >= 0 && b < solver_to_me.size()) res[g] = solver_to_me[b];
            }
            return res;
        }

        // How much of its bus' derivative each generator carries: 1 where it is the only
        // regulator of that bus, 1/n where n of them share it, 0 where it carries none
        // (the -1 entries of get_gen_v_target_bus).
        //
        // WHY A SHARE, AND NOT A DERIVATIVE EACH. Generators regulating one bus must be
        // given the same set-point -- a bus has one magnitude -- so the loss is a
        // function only ON the diagonal v_1 = ... = v_n. Off it there is no value to
        // compare against: such a row is refused, not solved differently. So the partial
        // derivative of one set-point with the others held fixed does not exist, and
        // what these n numbers are is NOT a gradient in the usual sense. What does exist
        // is the derivative along the tie, and that is what they sum to.
        //
        // Splitting it equally is the choice that behaves, for the two things a caller
        // actually does. Tie them (one parameter driving the group, which is how the
        // degree of freedom really looks) and the chain rule adds the shares back to the
        // true derivative. Treat them as separate parameters and step on all of them,
        // and they move together, so the iterate stays where the function is defined --
        // which giving the whole derivative to one of them does not do: the next row
        // would be refused. And nothing here depends on the order the generators happen
        // to sit in their container, which is what "whichever set_vm writes last" would
        // have made the answer depend on.
        RealVect get_gen_v_share() const { return _gen_v_layout().share; }

        // get_gen_v_share, row by row: 0 on a row where the generator's set-point reaches
        // nothing -- its bus switched PQ there, its group's voltage row lost its set-point
        // (every controller held, or stranded), the row not solved. What the DIRECT half of
        // the gen_v gradient is weighted by; gen_v_indirect_grad already is. Needs the
        // Jacobians kept (keep_jacobian), which is when it is asked for.
        BatchAdjoint::RealMatRM get_gen_v_share_rows() const;

        // The Jacobian row whose adjoint is each generator's gen_v gradient, -1 where
        // it has none there. A generator regulating a bus a voltage-control group holds
        // (a remote regulator, or a local one on a group-controlled bus) fixes no |V|:
        // the group's bordered row  F_v = |V_reg| + sum s.Q_c - v_set  keeps |V_reg| an
        // unknown, and the set-point is v_set. dF_v/dv_set = -1, so its gradient is
        // lambda at that row -- no direct half, and no dS/d|V| column. -1 as well for a
        // group holding an SVC or an hvdc station: their set-point no batch moves, so a
        // row moving the generator's is refused (_row_gen_v_conflicts) and there is no
        // derivative to give.
        IntVect get_gen_v_vc_row() const {
            const power_flow_adjoint::GenVLayout layout = _gen_v_layout();
            const IntVect v_row = _algo.get_group_v_row();
            IntVect res = IntVect::Constant(layout.vc_group.size(), -1);
            for(Eigen::Index g = 0; g < layout.vc_group.size(); ++g){
                const int grp = layout.vc_group[g];
                if(grp >= 0 && grp < v_row.size()) res[g] = v_row[grp];
            }
            return res;
        }

        // The indirect half of the gen_v gradient: `-lambda^T dF/dv` per row and per
        // generator -- for a generator of a voltage-control group, all of it: lambda at
        // the group's voltage row (get_gen_v_vc_row), zero on a row where that row lost
        // its set-point (stranded by handle_disconnected_grid, every controller held) --
        // otherwise keyed like get_gen_v_target_bus() (zero wherever that says -1, and
        // on a row that did not converge) and already carrying each generator's share of
        // its bus in that row (see get_gen_v_share_rows -- the caller must weight the
        // DIRECT half by the same thing). `lambda` is what solve_JT returned, so the
        // adjoint system is solved once and both halves of the gradient read it.
        //
        // dF/dv_k is the dS/d|V_k| column, taken on THIS row's own system -- its
        // contingency edits, and the taps and sections its outer loops settled on -- as
        // the algorithm handed it over with the row's Jacobian (AdjointSolveInfo).
        BatchAdjoint::RealMatRM gen_v_indirect_grad(const Eigen::Ref<const BatchAdjoint::RealMatRM> & lambda) const;

        // What the gen_p gradient owes beyond lambda at the generator's P row, per row and
        // per generator (AdjointSolveInfo::gen_p_terms): a unit an outer loop froze at a
        // limit of its capability curve injects a reactive power that follows its target
        // P; a unit the distributed slack (NROuterDist_*) holds at a bound is not moved by
        // its target at all, and this takes its P row back off. Zero everywhere without
        // such a unit. `lambda` is solve_JT's.
        BatchAdjoint::RealMatRM gen_p_indirect_grad(const Eigen::Ref<const BatchAdjoint::RealMatRM> & lambda) const;

        // The transformers each row's outer loops left at another tap than the grid's, with
        // the pi block (y11, y12, y21, y22, pu) they have there -- what that row's flows are
        // computed with: three aligned columns (row, transformer id, block). Empty without
        // outer loops, or when no loop moved a tap. For a caller recomputing flows from the
        // voltages (a differentiable flow), which would otherwise read the grid's tap.
        std::tuple<IntVect, IntVect, Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >
        get_outer_trafo_admittances() const {
            std::size_t n = 0;
            for(const auto & row : _row_trafo_blocks_) n += row.size();
            IntVect rows(static_cast<Eigen::Index>(n)), trafos(static_cast<Eigen::Index>(n));
            Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> blocks(static_cast<Eigen::Index>(n), 4);
            Eigen::Index k = 0;
            for(std::size_t i = 0; i < _row_trafo_blocks_.size(); ++i){
                for(const auto & tb : _row_trafo_blocks_[i]){
                    rows[k] = static_cast<int>(i);
                    trafos[k] = tb.first;
                    for(int c = 0; c < 4; ++c) blocks(k, c) = tb.second[static_cast<std::size_t>(c)];
                    ++k;
                }
            }
            return std::make_tuple(rows, trafos, blocks);
        }

        // Per row: 1 where solve_JT() actually solved that row's adjoint system.
        const std::vector<char> & adjoint_row_ok() const { return _adjoint_row_ok_; }

        // Linear-solver counters and timings of the last solve_JT: how much of the
        // backward went into refactorizing versus into the transposed solves.
        const LinearSolverStats & adjoint_solver_stats() const { return _adjoint_.get_linear_solver_stats(); }

        // ================= unified compute() ======================================
        void compute(const Eigen::Ref<const CplxVect> & Vinit, int max_iter, real_type tol);

        Eigen::Ref<RealMat > compute_flows() {
            _maybe_check_results_match_defaults("compute_flows");
            compute_flows_from_Vs();
            _maybe_clean_flows(true);
            return _amps_flows;
        }
        Eigen::Ref<RealMat > compute_power_flows() {
            _maybe_check_results_match_defaults("compute_power_flows");
            compute_flows_from_Vs(false);
            _maybe_clean_flows(false);
            return _active_power_flows;
        }

        // P and Q at both ends of every branch, every row: see get_branch_results for the
        // layout. A branch this row disconnects reads 0, like one out of service in the
        // grid, and so does every branch of a row that did not converge.
        const RealMat & compute_branch_results() {
            _maybe_check_results_match_defaults("compute_branch_results");
            compute_branch_results_from_Vs();
            const std::vector<std::vector<int> > tripped = _rows_tripped_branch_ids();
            const Eigen::Index nb_rows = _branch_results.rows();
            for(Eigen::Index i = 0; i < nb_rows; ++i){
                const size_t row = static_cast<size_t>(i);
                if(row >= _converged_mask_.size() || !_converged_mask_[row]){
                    _branch_results.row(i).setZero();
                    continue;
                }
                if(row >= tripped.size()) continue;
                for(int br_id : tripped[row]){
                    _branch_results.row(i).segment(4 * static_cast<Eigen::Index>(br_id), 4).setZero();
                }
            }
            return _branch_results;
        }

        // ---- the admittance matrix of one row -----------------------------------------
        // The matrix row `row` was solved with: the base one, minus the branches that row
        // disconnects. An entry those branches were the only contributors to is removed
        // from the pattern rather than left as an explicit zero, so the result has the
        // pattern a one-at-a-time solve of that topology builds (its values may differ in
        // the last bits: removing a branch is a subtraction, not a rebuild). Grid bus
        // numbering (size total_bus, like LSGrid.get_Ybus) unless `solver_numbering`.
        // AC only for get_Ybus, DC only for get_dcYbus: the batch builds the one matrix
        // its algorithm needs.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        Eigen::SparseMatrix<cplx_type> get_Ybus(Eigen::Index row, bool solver_numbering = false) const {
            _maybe_check_results_match_defaults("get_Ybus");
            _check_row_matrix_available(row, true, "get_Ybus");
            return _row_matrix(ac_cache_.mat, static_cast<size_t>(row), solver_numbering);
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        Eigen::SparseMatrix<real_type> get_dcYbus(Eigen::Index row, bool solver_numbering = false) const {
            _maybe_check_results_match_defaults("get_dcYbus");
            _check_row_matrix_available(row, false, "get_dcYbus");
            return _row_matrix(dc_cache_.mat, static_cast<size_t>(row), solver_numbering);
        }

        // ---- per-row solve statistics -------------------------------------------------
        // What the algorithm reported for each row of the last compute(): its time (s, the
        // algorithm's own get_computation_time) and its number of iterations. 0 for a row
        // that was not handed to the solver (skipped before it, or never reached).
        RealVect get_row_solve_times() const {
            return RealVect::Map(_row_solve_time_.data(), static_cast<Eigen::Index>(_row_solve_time_.size()));
        }
        IntVect get_row_nb_iter() const {
            return IntVect::Map(_row_nb_iter_.data(), static_cast<Eigen::Index>(_row_nb_iter_.size()));
        }
        // What the linear solver ran for each row: its counters are the row's own (what they
        // grew by during its solve), its timers and min_rgrowth_ the row's too (they cover
        // one solve). All 0 for a row never handed to the solver.
        const std::vector<LinearSolverStats> & get_row_linear_solver_stats() const { return _row_lin_stats_; }
        // The size of the Jacobian each row's last Newton iteration solved (0 for a row never
        // handed to the solver, or an algorithm without one): the same on every row for a
        // fixed sparsity pattern, the row's own for an algorithm that rebuilds it.
        IntVect get_row_J_dim() const {
            return IntVect::Map(_row_J_dim_.data(), static_cast<Eigen::Index>(_row_J_dim_.size()));
        }
        IntVect get_row_J_nnz() const {
            return IntVect::Map(_row_J_nnz_.data(), static_cast<Eigen::Index>(_row_J_nnz_.size()));
        }
        // What the outer loops did on each row of the last compute() (see OuterRowResult),
        // when the algorithm runs them (NROuter_*); empty otherwise.
        const std::vector<OuterRowResult> & get_outer_results() const { return _outer_results_; }

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}
        // both ContingencyAnalysis and ScenarioSweep edit Ybus per step.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        double modif_Ybus_time() const {return _timer_modif_Ybus;}
        // solve_time(): only ever bound for ContingencyAnalysis today (kept as a
        // synonym, same as the pre-refactor code -- see binding_batch.cpp).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        double solve_time() const {return _timer_solver;}

    protected:
        // ----- base-case reuse ----------------------------------------------------

        // The starting voltage, mapped onto the kept grid cache's bus labelling. This is
        // LSGrid::_pre_process_own_cache's own two steps -- gather the caller's V onto
        // the solver's buses, then snap every regulated bus to its target magnitude --
        // and it is redone on every call rather than cached, because a call is free to
        // start anywhere and doing it costs microseconds (a gather over the buses and
        // three passes over the regulating elements) against the tens of milliseconds
        // that keeping the grid cache saves.
        CplxVect _vinit_on_grid_cache(const Eigen::Ref<const CplxVect> & Vinit) const {
            CplxVect res = Vinit(active_layout().id_solver_to_me.as_eigen());
            const SolverBusIdVect & me_to_solver = active_layout().id_me_to_solver;
            seed_vm_keeping(res, _vm_held_buses(), [this, &me_to_solver](CplxVect & V){
                _grid_model.get_generators().set_vm(V, me_to_solver);
                _grid_model.get_dclines().set_vm(V, me_to_solver);
                _grid_model.get_svcs().set_vm(V, me_to_solver);
            });
            return res;
        }

        // ----- reverse-mode differentiation helpers -------------------------------

        // Re-key a solver-bus-keyed map onto grid bus ids (see get_p_row_of_bus).
        IntVect _bus_map_to_grid(const IntVect & solver_keyed) const {
            return power_flow_adjoint::bus_map_to_grid(active_layout(), solver_keyed);
        }

        // which part of the gen_v gradient each generator has (see PowerFlowAdjoint)
        power_flow_adjoint::GenVLayout _gen_v_layout() const {
            return power_flow_adjoint::gen_v_layout(_grid_model, active_layout(), _grid_model.get_ac_voltage_control_plan());
        }

        // Keep this row's converged Jacobian for the adjoint, and what the algorithm says
        // of the system it ends on (its identity rows, the dF/d|V| columns of the gen_v
        // buses, ...: AdjointSolveInfo). MUST be called while the row's own state is
        // still installed on the algorithm -- its Ybus edits, its masked buses, its PV
        // pinning, its outer loops' last decisions -- because refresh_J_at_solution()
        // re-fills J from exactly that, and a J refreshed after the loop restored its
        // resting state would describe a system this row never solved. Nothing happens
        // without keep_jacobian: the row loop pays one test.
        void _maybe_store_jacobian(size_t i, AlgorithmSelector & algo) {
            if(!_adjoint_.is_allocated()) return;
            algo.refresh_J_at_solution();
            _adjoint_.store_row(static_cast<Eigen::Index>(i), algo.get_J());
            if(i < _adjoint_info_.size()) algo.get_adjoint_solve_info(_adjoint_gen_v_.vm_buses, _adjoint_info_[i]);
        }

        // Which equations each row froze to the identity, as Jacobian row indices --
        // what solve_JT drops from lambda (see BatchAdjoint::solve_JT). Two sources,
        // and a row can have both:
        //   - a bus a contingency stranded: masked, so BOTH its P and its Q equation;
        //   - a bus still PV in this row (a reserved switchable bus whose generator
        //     is still on, the batch's or an outer loop's): pinned, so its Q equation
        //     only -- as the algorithm reported it with the row's Jacobian.
        // Empty when no row does either, which is the common case.
        std::vector<std::vector<int> > _adjoint_identity_rows() const {
            const Eigen::Index nb_rows = _adjoint_.nb_rows();
            // _li_masked is filled whenever connectivity was analysed, but only the
            // masked mode ever hands it to the algorithm: without it a contingency
            // that strands a bus makes the row skipped entirely, not masked
            const bool has_masking = _mask_mode() && !_li_masked.empty();
            bool has_pinning = false;
            for(const AdjointSolveInfo & info : _adjoint_info_) has_pinning = has_pinning || !info.pinned_buses.empty();
            if(nb_rows <= 0 || (!has_masking && !has_pinning)) return std::vector<std::vector<int> >();

            const IntVect p_row = _algo.get_p_to_J_row_python();
            const IntVect q_row = _algo.get_q_to_J_row_python();
            const AdjointSolveInfo none;
            std::vector<std::vector<int> > res(static_cast<size_t>(nb_rows));
            for(Eigen::Index i = 0; i < nb_rows; ++i){
                const size_t row = static_cast<size_t>(i);
                const std::vector<int> * masked = has_masking && row < _li_masked.size() ? &_li_masked[row] : nullptr;
                power_flow_adjoint::identity_rows(row < _adjoint_info_.size() ? _adjoint_info_[row] : none,
                                                  masked, p_row, q_row, res[row]);
            }
            return res;
        }

        // ----- bookkeeping shared by every instantiation (plain, always compiled) -
        void _check_ok_el(Eigen::Index el){
            if(el < 0 || el >= static_cast<Eigen::Index>(n_total_)){
                std::ostringstream exc_;
                exc_ << algo_name() << ": cannot add the contingency with id " << el
                     << "; the grid counts only " << n_total_ << " powerlines / trafos.";
                throw std::runtime_error(exc_.str());
            }
        }
        void _check_ok_id(int el, int nb, const char * kind, const char * fun_name) const {
            if(el < 0 || el >= nb){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": no " << kind << " with id " << el
                     << "; the grid counts " << nb << " of them.";
                throw std::runtime_error(exc_.str());
            }
        }
        // a branch-only contingency (add_n1 / add_nk / remove_*), ids already checked
        static ContingencyElements _branch_contingency(const std::vector<int> & branch_ids){
            ContingencyElements res;
            res.branches.insert(branch_ids.begin(), branch_ids.end());
            return res;
        }
        // what add_contingency* / remove_contingency are given, checked: every id exists,
        // and something is tripped at all. A branch opened at both ends one way or another
        // (named at both sides, or also among `branch_ids`) is in `branches` only, as
        // OpenLoadFlow merges the two ends of a branch it opens.
        ContingencyElements _checked_contingency(const std::vector<int> & branch_ids,
                                                 const std::vector<int> & shunt_ids,
                                                 const std::vector<int> & load_ids,
                                                 const std::vector<int> & gen_ids,
                                                 const char * fun_name,
                                                 const std::vector<int> & branch_ids_side1 = std::vector<int>(),
                                                 const std::vector<int> & branch_ids_side2 = std::vector<int>()) const {
            ContingencyElements res;
            for(int el : branch_ids){ _check_ok_id(el, static_cast<int>(n_total_), "powerline / trafo", fun_name); res.branches.insert(el); }
            for(int el : branch_ids_side1){ _check_ok_id(el, static_cast<int>(n_total_), "powerline / trafo", fun_name); res.branches_side1.insert(el); }
            for(int el : branch_ids_side2){
                _check_ok_id(el, static_cast<int>(n_total_), "powerline / trafo", fun_name);
                if(res.branches_side1.erase(el) > 0) res.branches.insert(el);
                else res.branches_side2.insert(el);
            }
            for(int el : res.branches){ res.branches_side1.erase(el); res.branches_side2.erase(el); }
            for(int el : shunt_ids){ _check_ok_id(el, _grid_model.get_shunts().nb(), "shunt", fun_name); res.shunts.insert(el); }
            for(int el : load_ids){ _check_ok_id(el, _grid_model.get_loads().nb(), "load", fun_name); res.loads.insert(el); }
            for(int el : gen_ids){ _check_ok_id(el, _grid_model.get_generators().nb(), "generator", fun_name); res.gens.insert(el); }
            if(res.empty()){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": the contingency trips nothing (every list is empty).";
                throw std::runtime_error(exc_.str());
            }
            return res;
        }
        void _check_cols(const Eigen::Ref<const RealMat> & data, Eigen::Index nb_expected_cols, const char * setter_name) const {
            if(data.cols() != nb_expected_cols){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << setter_name << ": got " << data.cols()
                     << " columns, the grid has " << nb_expected_cols << " such elements.";
                throw std::runtime_error(exc_.str());
            }
        }
        void _check_cols_bool(const Eigen::Ref<const BoolMat> & data, Eigen::Index nb_expected_cols, const char * setter_name) const {
            if(data.cols() != nb_expected_cols){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << setter_name << ": got " << data.cols()
                     << " columns, the grid has " << nb_expected_cols << " such elements.";
                throw std::runtime_error(exc_.str());
            }
        }
        void _lock_or_check_nb_steps(Eigen::Index rows, const char * setter_name){
            _results_stale_ = true;
            if(_nb_steps_locked < 0){ _nb_steps_locked = rows; return; }
            if(rows != _nb_steps_locked){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << setter_name << ": got " << rows
                     << " rows, but a previous call already fixed the number of "
                        "simulations to " << _nb_steps_locked << ". Every modify_* / "
                        "set_contingency_* call must agree on the number of rows.";
                throw std::runtime_error(exc_.str());
            }
        }
        void _check_limit_violations_enabled(const std::string & fun_name) const {
            if(!_compute_limit_violations_){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": limit violations were not "
                        "requested. Construct this object with `compute_limit_violations=True` "
                        "to use this feature.";
                throw std::runtime_error(exc_.str());
            }
        }
        void _check_physical_violations_enabled(const std::string & fun_name) const {
            if(!_compute_physical_violations_){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": the physical-limit checks were not "
                        "requested. Set `compute_physical_violations = True` before compute() "
                        "to use this feature.";
                throw std::runtime_error(exc_.str());
            }
        }
        // `base_w` is the weight vector to mask, so a row that already re-weighted the
        // slack for its own generator contingency (see _row_slack_weights) is masked
        // on top of that rather than back on the layout's untouched weights. Returns
        // `base_w` itself when nothing is masked -- the common row -- and otherwise
        // the masked copy written into `scratch` (which `base_w` may already be).
        const RealVect & _masked_slack_weights(const std::vector<int> & masked,
                                               const RealVect & base_w,
                                               RealVect & scratch) const {
            if(masked.empty()) return base_w;
            if(&scratch != &base_w) scratch = base_w;
            const real_type orig_sum = scratch.sum();
            for(int b : masked) if(b >= 0 && b < scratch.size()) scratch(b) = 0.;
            const real_type new_sum = scratch.sum();
            if(new_sum > 1e-12 && orig_sum > 1e-12) scratch *= (orig_sum / new_sum);
            return scratch;
        }

        // BFS connected-component labelling: returns the solver bus ids NOT part of
        // the largest connected component (empty if the whole matrix is connected).
        // The search _prepare_connectivity falls back to for the contingencies the
        // base graph's DFS tree cannot settle on its own (an N-k that removes several
        // tree edges, a base graph that is not connected to begin with): its answer
        // is the reference the tree's answer must agree with, bus for bus. Reached
        // from pick_reference_slack() (public) hence this plain (non-SFINAE, but
        // member-template-on-T) helper must also be fully inline.
        template<typename T>
        std::vector<int> _disconnected_buses(const Eigen::SparseMatrix<T> & mat) const {
            const int n = static_cast<int>(mat.cols());
            std::vector<int> comp_of_bus(n, -1);
            int nb_comp = 0;
            std::queue<int> neighborhood;
            for(int start = 0; start < n; ++start){
                if(comp_of_bus[start] != -1) continue;
                comp_of_bus[start] = nb_comp;
                neighborhood.push(start);
                while(!neighborhood.empty()){
                    const int col_id = neighborhood.front();
                    neighborhood.pop();
                    for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, col_id); it; ++it){
                        const int row = static_cast<int>(it.row());
                        if(comp_of_bus[row] == -1 && std::abs(it.value()) > 1e-8){
                            comp_of_bus[row] = nb_comp;
                            neighborhood.push(row);
                        }
                    }
                }
                ++nb_comp;
            }

            if(nb_comp <= 1) return std::vector<int>();

            std::vector<int> nb_bus_per_comp(nb_comp, 0);
            for(int bus_id = 0; bus_id < n; ++bus_id) nb_bus_per_comp[comp_of_bus[bus_id]] += 1;
            const int main_comp = static_cast<int>(std::distance(
                nb_bus_per_comp.begin(),
                std::max_element(nb_bus_per_comp.begin(), nb_bus_per_comp.end())));

            std::vector<int> masked;
            masked.reserve(n - nb_bus_per_comp[main_comp]);
            for(int bus_id = 0; bus_id < n; ++bus_id){
                if(comp_of_bus[bus_id] != main_comp) masked.push_back(bus_id);
            }
            return masked;
        }

        // Connectivity of every contingency, settled once per compute() before the
        // row loop: _li_masked[i] is the (sorted) list of solver buses contingency i
        // strands -- empty when the grid stays in one piece -- and _cont_connected_[i]
        // says just that. One depth-first search of the base graph (BusGraph) answers
        // an N-1 in constant time and lists the stranded side in time proportional to
        // its size; the contingencies it cannot settle (an N-k removing several tree
        // edges, a base graph that is not connected) get the breadth-first labelling
        // this used to run for every one of them, on a patched copy of the matrix
        // made only if such a contingency exists. Keyed off ybus_policy_.li_coeffs,
        // representation-agnostic (li_defaults on ContingencyAnalysis, the masks on
        // ScenarioSweep). Reached from the public is_grid_connected_after_contingency()
        // / pick_reference_slack(), hence fully inline.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _prepare_connectivity(){
            const size_t nb_cont = ybus_policy_.li_coeffs.size();
            _li_masked.assign(nb_cont, std::vector<int>());
            _cont_connected_.assign(nb_cont, 1);
            const bool ac_solver_used = _algo.ac_solver_used();
            const real_type threshold = BusGraph::default_threshold();
            // the union layout's extra buses (_base_masked_) are isolated in the base
            // matrix: left out of the tree, which would otherwise never be connected and
            // send every row to the search below as soon as one action creates a bus
            if(ac_solver_used) bus_graph_.build(ac_cache_.mat, threshold, _base_masked_);
            else bus_graph_.build(dc_cache_.mat, threshold, _base_masked_);

            std::vector<std::pair<int, int> > edges;
            // the patched copy the fallback search walks, made on first need
            Eigen::SparseMatrix<cplx_type> ybus_work;
            Eigen::SparseMatrix<real_type> bbus_work;
            bool work_ready = false;
            for(size_t cont_id = 0; cont_id < nb_cont; ++cont_id){
                const std::vector<Coeff> & coeffs = ybus_policy_.li_coeffs[cont_id];
                BusGraph::Cut cut;
                cut.verdict = BusGraph::Verdict::Unknown;
                cut.child = -1;
                // a row that places a branch itself (set_topo_actions: an end moved, a
                // branch reconnected) ADDS edges the base tree cannot represent: the
                // search settles it
                const bool places = cont_id < ybus_policy_.topo_branches_moved.size() &&
                                    !ybus_policy_.topo_branches_moved[cont_id].empty();
                const bool settled = !places && (ac_solver_used
                    ? BusGraph::removed_edges(ac_cache_.mat, coeffs, threshold, edges)
                    : BusGraph::removed_edges(dc_cache_.mat, coeffs, threshold, edges));
                if(settled) cut = bus_graph_.cut(edges);
                // the tree answers for the buses it holds; the ones it leaves out are cut
                // off in every row it settles (a row that places a branch is searched),
                // stranding whatever the row puts on them
                if(cut.verdict == BusGraph::Verdict::Connected){
                    _li_masked[cont_id] = _base_masked_;
                    _cont_connected_[cont_id] = _only_masked_by_design(cont_id, _li_masked[cont_id]) ? 1 : 0;
                    continue;
                }
                if(cut.verdict == BusGraph::Verdict::Split){
                    bus_graph_.stranded_buses(cut, _li_masked[cont_id]);
                    if(!_base_masked_.empty()){
                        std::vector<int> merged;
                        std::set_union(_li_masked[cont_id].begin(), _li_masked[cont_id].end(),
                                       _base_masked_.begin(), _base_masked_.end(), std::back_inserter(merged));
                        _li_masked[cont_id].swap(merged);
                    }
                } else if(ac_solver_used){
                    if(!work_ready){ ybus_work = ac_cache_.mat; work_ready = true; }
                    for(const auto & c: coeffs) ybus_work.coeffRef(c.row_id, c.col_id) -= c.value;
                    _li_masked[cont_id] = _disconnected_buses(ybus_work);
                    for(const auto & c: coeffs) ybus_work.coeffRef(c.row_id, c.col_id) += c.value;
                } else {
                    if(!work_ready){ bbus_work = dc_cache_.mat; work_ready = true; }
                    for(const auto & c: coeffs) bbus_work.coeffRef(c.row_id, c.col_id) -= std::real(c.value);
                    _li_masked[cont_id] = _disconnected_buses(bbus_work);
                    for(const auto & c: coeffs) bbus_work.coeffRef(c.row_id, c.col_id) += std::real(c.value);
                }
                // connected: nothing stranded beyond the buses masked by design (see
                // _only_masked_by_design)
                _cont_connected_[cont_id] = _only_masked_by_design(cont_id, _li_masked[cont_id]) ? 1 : 0;
            }
        }

        // whether `masked` (sorted), the buses row `row` cuts off, holds nothing but buses
        // masked by design: the extra buses of the union layout (_base_masked_, in every
        // row that does not use them) and the base buses the row's action leaves with no
        // element at all (_row_emptied_buses_: a merge). Nothing stands on either, so
        // nothing is stranded. An extra bus the row puts an element on
        // (_row_placed_buses_) is not one of them: that element is stranded.
        bool _only_masked_by_design(size_t row, const std::vector<int> & masked) const {
            const std::vector<int> * emptied = row < _row_emptied_buses_.size() ? &_row_emptied_buses_[row] : nullptr;
            const std::vector<int> * placed = row < _row_placed_buses_.size() ? &_row_placed_buses_[row] : nullptr;
            for(int b : masked){
                if(std::binary_search(_base_masked_.begin(), _base_masked_.end(), b) &&
                   (placed == nullptr || !std::binary_search(placed->begin(), placed->end(), b))) continue;
                if(emptied != nullptr && std::binary_search(emptied->begin(), emptied->end(), b)) continue;
                return false;
            }
            return true;
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _prepare_connectivity(){}

        // pre-pass run before the (n-)powerflow in "handle disconnected grid" mode:
        // from _li_masked (see _prepare_connectivity, which must have run), chooses
        // the reference slack, moves it to index 0 of BOTH slack_bus_id_solver (the
        // one the solver reads) and slack_bus_id_me (kept aligned with it), and fills
        // _skip_mask with the contingencies that strand it.
        // The reference is, with `respect_forced`, the one forced on the grid
        // (LSGrid::set_reference_slack_bus): it is chosen for the whole batch, so it
        // is kept even if some contingencies strand it (those are skipped). Otherwise,
        // or if the forced bus is not a slack bus with a positive weight, it is the
        // slack bus stranded by the fewest contingencies.
        // Everything is compared in SOLVER numbering (_li_masked and slack_weights
        // are solver-indexed), the slack list being walked by position.
        // Called from _maybe_prepare_masks() (compute()-only, both instantiations,
        // respect_forced) and pick_reference_slack() (public, ContingencyAnalysis-only,
        // a suggestion: not respect_forced; hence this too must be fully inline).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _select_ref_slack_and_masks(bool respect_forced){
            const size_t nb_cont = ybus_policy_.li_coeffs.size();
            _skip_mask.assign(nb_cont, 0);

            SolverBusLayout & layout = active_layout();
            const int nb_slack = static_cast<int>(layout.slack_bus_id_solver.size());
            if(nb_slack == 0 || static_cast<int>(layout.slack_bus_id_me.size()) != nb_slack) return;
            const RealVect & weights = layout.slack_weights;
            const auto weight_of = [&weights](int bus_solver){
                return (bus_solver >= 0 && bus_solver < weights.size()) ? weights(bus_solver) : real_type(0.);
            };

            // positions (in the slack list) of the candidates
            std::vector<int> candidates;
            for(int k = 0; k < nb_slack; ++k){
                if(weight_of(layout.slack_bus_id_solver[k].cast_int()) > 0.) candidates.push_back(k);
            }
            if(candidates.empty()){
                for(int k = 0; k < nb_slack; ++k) candidates.push_back(k);
            }

            // every _li_masked entry is sorted (see _prepare_connectivity)
            auto is_masked = [](const std::vector<int> & masked, int bus){
                return std::binary_search(masked.begin(), masked.end(), bus);
            };

            int best_k = -1;
            const int forced = respect_forced ? _grid_model.get_reference_slack_bus() : -1;
            if(forced >= 0){
                for(int k : candidates){
                    if(layout.slack_bus_id_me[k].cast_int() == forced){ best_k = k; break; }
                }
            }
            if(best_k < 0){
                int best_strand = -1;
                real_type best_weight = -1.;
                int best_bus = -1;
                for(int k : candidates){
                    const int bus = layout.slack_bus_id_solver[k].cast_int();
                    int strand = 0;
                    for(const auto & masked : _li_masked) if(is_masked(masked, bus)) ++strand;
                    const real_type weight = weight_of(bus);
                    const bool better = (best_strand < 0) ||
                                        (strand < best_strand) ||
                                        (strand == best_strand && weight > best_weight) ||
                                        (strand == best_strand && weight == best_weight && bus < best_bus);
                    if(better){ best_strand = strand; best_weight = weight; best_bus = bus; best_k = k; }
                }
            }

            if(best_k > 0){
                // same permutation on both lists: best_k to the front, the others in order
                std::vector<int> me = layout.slack_bus_id_me.to_int_vector();
                std::vector<int> solver = layout.slack_bus_id_solver.to_int_vector();
                std::rotate(me.begin(), me.begin() + best_k, me.begin() + best_k + 1);
                std::rotate(solver.begin(), solver.begin() + best_k, solver.begin() + best_k + 1);
                layout.slack_bus_id_me = GlobalBusIdVect(me);
                layout.slack_bus_id_solver = SolverBusIdVect(solver);
            }

            const int ref_solver = layout.slack_bus_id_solver[0].cast_int();
            for(size_t cont_id = 0; cont_id < nb_cont; ++cont_id){
                if(is_masked(_li_masked[cont_id], ref_solver)) _skip_mask[cont_id] = 1;
                // the masked loop runs for the topological actions' sake as well: without
                // handle_disconnected_grid a row that strands a bus keeps its legacy
                // answer, NOT_SIMULATED
                if(!_handle_disconnected_grid && !_only_masked_by_design(cont_id, _li_masked[cont_id])) _skip_mask[cont_id] = 1;
            }
        }

        // ----- number of steps: SFINAE-dispatched (ContingencyAnalysis counts its
        // distinct scenarios; every other instantiation reads the setter-locked
        // count) -----
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        size_t _nb_steps() const { return ybus_policy_.li_defaults.size(); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        size_t _nb_steps() const {
            if(_nb_steps_locked < 0){
                std::ostringstream exc_;
                exc_ << algo_name() << "::compute: nothing was ever set (no modify_* / "
                        "set_contingency_* call) -- nothing to compute. Call the relevant "
                        "setter(s) first, or use ac_pf()/dc_pf() directly for a single powerflow.";
                throw std::runtime_error(exc_.str());
            }
            return static_cast<size_t>(_nb_steps_locked);
        }

        // ----- Ybus-varying preparation: 3-way (NOOP / Contingency-list / masks) --
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _prepare_ybus_varying(bool, Eigen::Index) {}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _prepare_ybus_varying(bool ac_solver_used, Eigen::Index) {
            ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, active_layout().id_me_to_solver, n_line_);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _prepare_ybus_varying(bool ac_solver_used, Eigen::Index nb_steps) {
            ybus_policy_.init_li_coeffs_from_masks(_grid_model, ac_solver_used, active_layout().id_me_to_solver, n_line_, nb_steps);
        }

        // ----- Sbus-varying preparation: 2-way --------------------------------
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _prepare_sbus_varying(bool, Eigen::Index) {}
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _prepare_sbus_varying(bool ac_solver_used, Eigen::Index nb_steps) {
            const CplxVect complete = ac_solver_used ? CplxVect(ac_cache_.inj) : CplxVect(dc_cache_.inj.template cast<cplx_type>());
            sbus_policy_.prepare(_grid_model, ac_solver_used, nb_buses_solver_, active_layout().id_me_to_solver,
                                 complete, _grid_model.get_sn_mva(), nb_steps, algo_name());
        }

        // ----- reset ybus_policy_'s own state: 2-way --------------------------
        // The COEFFICIENTS only (L1): values read off Ybus / Bbus, so they go with the
        // matrix they were read from. The contingencies themselves are registrations
        // and stay -- _prepare_ybus_varying rebuilds the coefficients for them.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _reset_ybus_coeffs() {
            ybus_policy_.li_coeffs.clear();
            ybus_policy_.li_sbus_dp.clear();
            ybus_policy_.li_rows.clear();
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _reset_ybus_coeffs() {}

        // everything, the registrations included: see clear().
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _reset_ybus_policy() {
            ybus_policy_.li_defaults.clear();
            ybus_policy_.li_coeffs.clear();
            ybus_policy_.li_sbus_dp.clear();
            ybus_policy_.li_rows.clear();
            ybus_policy_.line_mask = BoolMat();
            ybus_policy_.trafo_mask = BoolMat();
            ybus_policy_.shunt_mask = BoolMat();
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _reset_ybus_policy() {}

        // ----- per-step Ybus edit (invertibility as the return value): 2-way -----
        // AC: whether the grid stays connected was settled for every contingency
        // before the loop (_prepare_connectivity). DC: the solver takes the
        // responsibility, so the matrix is always reported "connected", as before.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool _remove_step_coeffs(Eigen::SparseMatrix<cplx_type> & Ybus, size_t i, bool ac_solver_used, AlgorithmSelector & algo) {
            YbusPolicy::Contingency::remove_from_Ybus(Ybus, ybus_policy_.li_coeffs[i], ac_solver_used, algo);
            return !ac_solver_used || (i < _cont_connected_.size() && _cont_connected_[i] != 0);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        bool _remove_step_coeffs(Eigen::SparseMatrix<cplx_type> &, size_t, bool, AlgorithmSelector &) { return true; }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _readd_step_coeffs(Eigen::SparseMatrix<cplx_type> & Ybus, size_t i, bool ac_solver_used, AlgorithmSelector & algo) {
            YbusPolicy::Contingency::readd_to_Ybus(Ybus, ybus_policy_.li_coeffs[i], ac_solver_used, algo);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _readd_step_coeffs(Eigen::SparseMatrix<cplx_type> &, size_t, bool, AlgorithmSelector &) {}

        // ----- per-row Ybus, for get_Ybus / get_dcYbus -------------------------------
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _check_row_matrix_available(Eigen::Index row, bool ac, const char * fun_name) const {
            std::ostringstream exc_;
            if(_algo.ac_solver_used() != ac){
                exc_ << algo_name() << "::" << fun_name << ": the last compute() used a "
                     << (ac ? "DC" : "AC") << " algorithm, which builds no "
                     << (ac ? "AC Ybus: use get_dcYbus." : "DC Bbus: use get_Ybus.");
            } else if(row < 0 || static_cast<size_t>(row) >= ybus_policy_.li_coeffs.size()){
                exc_ << algo_name() << "::" << fun_name << ": no row " << row << " (rows of the last "
                     "compute(): " << ybus_policy_.li_coeffs.size() << "). Call compute() first.";
            } else {
                return;
            }
            throw std::runtime_error(exc_.str());
        }
        static real_type _as_scalar(const cplx_type & val, real_type) { return std::real(val); }
        static cplx_type _as_scalar(const cplx_type & val, cplx_type) { return val; }

        template<typename Scalar, class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        Eigen::SparseMatrix<Scalar> _row_matrix(const Eigen::SparseMatrix<Scalar> & base, size_t row,
                                                bool solver_numbering) const {
            using index_type = typename Eigen::SparseMatrix<Scalar>::StorageIndex;
            // an entry the row's branches were the only contributors to comes out at 0, or
            // within rounding of it: (a + b) - b need not be exactly a. Relative to what the
            // entry held, that residual is at the level of the machine epsilon; a genuine
            // remaining contribution is not.
            const real_type rel_tol = 1e-12;
            Eigen::SparseMatrix<Scalar> mat = base;
            std::vector<std::pair<std::pair<int, int>, real_type> > touched;  // (row, col), |base value|
            for(const Coeff & c : ybus_policy_.li_coeffs[row]){
                const int r = static_cast<int>(c.row_id);
                const int col = static_cast<int>(c.col_id);
                touched.push_back({{r, col}, std::abs(base.coeff(r, col))});
                mat.coeffRef(r, col) -= _as_scalar(c.value, Scalar());
            }
            std::sort(touched.begin(), touched.end());

            const GlobalBusIdVect & id_solver_to_me = active_layout().id_solver_to_me;
            const Eigen::Index size = solver_numbering ? mat.rows() : static_cast<Eigen::Index>(_grid_model.total_bus());
            std::vector<Eigen::Triplet<Scalar> > triplets;
            triplets.reserve(static_cast<size_t>(mat.nonZeros()));
            for(Eigen::Index col_ = 0; col_ < mat.outerSize(); ++col_){
                for(typename Eigen::SparseMatrix<Scalar>::InnerIterator it(mat, col_); it; ++it){
                    const std::pair<int, int> key(static_cast<int>(it.row()), static_cast<int>(it.col()));
                    const auto found = std::lower_bound(touched.begin(), touched.end(),
                                                        std::make_pair(key, real_type(-1.)));
                    if(found != touched.end() && found->first == key &&
                       std::abs(it.value()) <= rel_tol * found->second) continue;
                    const index_type r = solver_numbering ? static_cast<index_type>(it.row())
                                                          : static_cast<index_type>(id_solver_to_me[static_cast<size_t>(it.row())].cast_int());
                    const index_type c = solver_numbering ? static_cast<index_type>(it.col())
                                                          : static_cast<index_type>(id_solver_to_me[static_cast<size_t>(it.col())].cast_int());
                    triplets.push_back({r, c, it.value()});
                }
            }
            Eigen::SparseMatrix<Scalar> res(size, size);
            res.setFromTriplets(triplets.begin(), triplets.end());
            res.makeCompressed();
            return res;
        }

        // ----- per-row branches a row disconnects, all rows at once: 3-way ------------
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        std::vector<std::vector<int> > _rows_tripped_branch_ids() const { return {}; }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<std::vector<int> > _rows_tripped_branch_ids() const { return _branches_out_vect(); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        std::vector<std::vector<int> > _rows_tripped_branch_ids() const {
            std::vector<std::vector<int> > res(static_cast<size_t>(_nb_result_rows()));
            for(size_t row = 0; row < res.size(); ++row){
                res[row] = ybus_policy_.branch_ids_for_row(static_cast<Eigen::Index>(row), n_line_);
            }
            return res;
        }

        // ----- DC: the buses a row releases get the caller's magnitude back ----------
        // (see _dc_row_vm_reset_). A bus is released by a row that disconnects every
        // generator writing its magnitude when nothing else (SVC, HVDC station, storage
        // unit) writes it too: the same walk set_vm does, so the magnitude is the one a
        // dc_pf of that topology starts from. Rebuilt every compute(): the masks are a
        // batch input, the starting voltage is this call's. A generator contingency
        // (ScenarioSweep's masks, ContingencyAnalysis' contingencies) is what releases one.
        void _prepare_dc_row_vm_reset(const Eigen::Ref<const CplxVect> & Vinit, bool ac_solver_used, size_t nb_steps){
            _dc_row_vm_reset_.clear();
            if(ac_solver_used || !_any_gen_off()) return;

            std::map<int, std::vector<int> > gen_writers;  // grid bus -> generators writing it
            std::set<int> other_writers;                   // grid buses something else writes
            _grid_model.get_generators().for_each_vm_target(
                [&gen_writers](int gen_id, int bus, real_type, const auto &){ gen_writers[bus].push_back(gen_id); });
            const auto mark_other = [&other_writers](int, int bus, real_type, const auto &){ other_writers.insert(bus); };
            _grid_model.get_svcs().for_each_vm_target(mark_other);
            _grid_model.get_dclines().for_each_vm_target(mark_other);
            _grid_model.get_storages().for_each_vm_target(mark_other);

            const SolverBusIdVect & me_to_solver = active_layout().id_me_to_solver;
            _dc_row_vm_reset_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            for(size_t row = 0; row < nb_steps; ++row){
                for(const auto & bus_gens : gen_writers){
                    if(other_writers.count(bus_gens.first)) continue;
                    bool all_off = true;
                    for(int gen_id : bus_gens.second){
                        if(!_gen_off_in_row(row, gen_id)){ all_off = false; break; }
                    }
                    if(!all_off) continue;
                    const int bus_solver = me_to_solver[static_cast<size_t>(bus_gens.first)].cast_int();
                    if(bus_solver == BaseConstants::_deactivated_bus_id) continue;
                    _dc_row_vm_reset_[row].push_back({bus_solver, std::abs(Vinit(bus_gens.first))});
                }
            }
        }
        void _apply_step_dc_vm_reset(size_t i, CplxVect & V) const {
            if(i >= _dc_row_vm_reset_.size()) return;
            for(const std::pair<int, real_type> & bus_vm : _dc_row_vm_reset_[i]){
                V(bus_vm.first) = std::polar(bus_vm.second, std::arg(V(bus_vm.first)));
            }
        }

        // ----- per-row solve statistics: written by the row loops, see get_row_nb_iter.
        // One slot per row, pre-sized by compute(): each worker thread writes its own
        // rows only, so no two threads touch one element.
        // `lin_before`: the algorithm's linear-solver stats right before the row's solve.
        void _record_row_solve_stats(size_t i, const AlgorithmSelector & algo, const LinearSolverStats & lin_before){
            if(i >= _row_solve_time_.size()) return;
            _row_solve_time_[i] = algo.get_computation_time();
            _row_nb_iter_[i] = algo.get_nb_iter();
            _record_row_lin_stats(i, algo, lin_before);
            _record_row_outer(i, algo);
        }

        void _record_row_lin_stats(size_t i, const AlgorithmSelector & algo, const LinearSolverStats & lin_before){
            if(i >= _row_lin_stats_.size()) return;
            // the counters cover the algorithm's lifetime: the row's is what they grew by;
            // the timers (and min_rgrowth_) are reset at the start of each solve already
            LinearSolverStats r = algo.get_linear_solver_stats();
            r.nb_reset -= lin_before.nb_reset;
            r.nb_analyze -= lin_before.nb_analyze;
            r.nb_factorize -= lin_before.nb_factorize;
            r.nb_refactorize -= lin_before.nb_refactorize;
            r.nb_refactorize_failed -= lin_before.nb_refactorize_failed;
            r.nb_fallback_factorize -= lin_before.nb_fallback_factorize;
            r.nb_fallback_factorize_failed -= lin_before.nb_fallback_factorize_failed;
            r.nb_solve -= lin_before.nb_solve;
            r.nb_solve_transpose -= lin_before.nb_solve_transpose;
            _row_lin_stats_[i] = r;
            if(algo.supports_jacobian()){
                const auto J = algo.get_J();
                _row_J_dim_[i] = static_cast<int>(J.rows());
                _row_J_nnz_[i] = static_cast<int>(J.nonZeros());
            }
        }

        // the pi block of each transformer whose taps row i's loops moved, at the row's
        // positions (the grid's for the changer they left alone): what the row's flows read
        void _record_row_trafo_blocks(size_t i, const OuterRowResult & r){
            if(i >= _row_trafo_blocks_.size()) return;
            auto & out = _row_trafo_blocks_[i];
            out.clear();
            if(r.phase_tap.empty() && r.ratio_tap.empty()) return;
            out = _grid_model.get_trafos().moved_tap_blocks(r.ratio_tap, r.phase_tap);
        }

        // what the outer loops did on row i (an algorithm that runs them only, see
        // _outer_results_), read while the row's solve is still the algorithm's last
        void _record_row_outer(size_t i, const AlgorithmSelector & algo){
            if(i >= _outer_results_.size()) return;
            OuterRowResult & r = _outer_results_[i];
            r = OuterRowResult();
            const OuterLoopStats stats = algo.get_outer_loop_stats();
            r.status = stats.status;
            r.unrealistic_state = stats.unrealistic_state;
            r.failed_loop = stats.failed_loop;
            r.nb_outer_iterations = stats.nb_outer_iterations;
            for(int n : stats.nr_iterations) r.nb_nr_iterations += n;
            std::map<std::string, int> taken;
            for(const OuterDecision & d : stats.decisions){
                if(d.taken) ++taken[d.loop + "/" + d.action];
            }
            r.decisions_taken.assign(taken.begin(), taken.end());
            std::vector<real_type> gen_p, storage_p;
            algo.get_outer_target_p(gen_p, storage_p);
            auto moved = [](const std::vector<real_type> & v, std::vector<std::pair<int, real_type> > & out){
                for(std::size_t k = 0; k < v.size(); ++k) if(!std::isnan(v[k])) out.emplace_back(static_cast<int>(k), v[k]);
            };
            moved(gen_p, r.gen_target_p);
            moved(storage_p, r.storage_target_p);
            std::vector<int> ints;
            algo.get_outer_hvdc_status(ints);
            for(std::size_t k = 0; k < ints.size(); ++k){
                if(ints[k] != HvdcRegimeControl::KEEP) r.hvdc_regime.emplace_back(static_cast<int>(k), ints[k]);
            }
            // a position the loops handled, where it differs from the grid's
            const TrafoContainer & trafos = _grid_model.get_trafos();
            algo.get_outer_phase_tap(ints);
            r.phase_tap = trafos.moved_taps(true, ints);
            algo.get_outer_ratio_tap(ints);
            r.ratio_tap = trafos.moved_taps(false, ints);
            _record_row_trafo_blocks(i, r);
            algo.get_outer_shunt_sections(ints);
            const ShuntContainer & shunts = _grid_model.get_shunts();
            for(std::size_t k = 0; k < ints.size(); ++k){
                if(ints[k] == std::numeric_limits<int>::min()) continue;
                const int el = static_cast<int>(k);
                if(el < shunts.nb() && shunts.has_sections(el) && shunts.get_section_count(el) == ints[k]) continue;
                r.shunt_sections.emplace_back(el, ints[k]);
            }
        }

        // ----- per-step Sbus: 2-way (row of sbus_policy_ vs. the fixed member) ---
        // Row i's injection: built into `scratch` (one buffer per range: no row
        // allocates, no two threads share it) where it varies, the fixed vector
        // itself where it does not. Nothing the size of nb_steps x nb_bus exists.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        const CplxVect & _step_sbus(size_t i, CplxVect & scratch) const {
            sbus_policy_.fill_row(static_cast<Eigen::Index>(i), scratch);
            return scratch;
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        const CplxVect & _step_sbus(size_t, CplxVect &) const { return ac_cache_.inj; }

        // ----- per-step generator vm seeding: 2-way (SbusPolicy::Vary::gen_v row
        // applied via GeneratorContainer::set_vm, vs. a no-op leaving V's magnitude
        // alone). Unlike _step_sbus above this mutates V in place rather than
        // returning a value: set_vm's contract is "rescale |V| at each PV
        // generator's regulated bus", there is nothing to hand back. -----
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _apply_step_gen_v(size_t, CplxVect &) const {}
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _apply_step_gen_v(size_t i, CplxVect & V) const {
            if(sbus_policy_.gen_v.rows() == 0) return;  // never set: the grid's own
                // target_vm_pu_ was already seeded once by _finish_preprocessing,
                // before the very first step -- nothing to redo here, no per-row
                // regression.
            // a row of the row-major gen_v is contiguous: a view binds to set_vm's
            // Eigen::Ref<const RealVect> parameter without a copy (as _step_sbus)
            const Eigen::Map<const RealVect> target_vm_pu_row(sbus_policy_.gen_v.row(static_cast<Eigen::Index>(i)).data(),
                                                              sbus_policy_.gen_v.cols());
            seed_vm_keeping(V, _vm_held_buses(), [this, i, &target_vm_pu_row](CplxVect & V_seed){
                _grid_model.get_generators().set_vm(V_seed, active_layout().id_me_to_solver, target_vm_pu_row);
                _reseat_shared_buses(i, target_vm_pu_row, V_seed);
            });
        }

        // ----- modify_gen_v on a voltage-control group ------------------------------
        // A generator regulating a bus a group holds fixes no |V| there -- the bus keeps
        // its magnitude unknown, pinned by the group's bordered voltage row -- so
        // set_vm's re-seed above only moves the starting point, and the row puts the
        // magnitude back at the grid's own v_set. The set-point that row reads is what a
        // per-row gen_v has to change: handed to the algorithm, group by group, before
        // each row's solve. The (generator, group) pairs depend on the grid only, found
        // once per compute() by _prepare_gen_v_vc.
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _apply_step_vc_v_set(size_t, AlgorithmSelector &) const {}
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _apply_step_vc_v_set(size_t i, AlgorithmSelector & algo) const {
            if(_gen_v_vc_pairs_.empty()){
                // nothing to vary: the grid's own set-points (and a worker a previous
                // compute() left with an override gets it dropped)
                algo.set_voltage_control_v_set(RealVect());
                return;
            }
            RealVect v_set = RealVect::Constant(_gen_v_vc_n_groups_,
                                                std::numeric_limits<real_type>::quiet_NaN());
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(row < sbus_policy_.gen_v.rows()){
                for(const auto & gen_and_group : _gen_v_vc_pairs_){
                    const real_type val = sbus_policy_.gen_v(row, gen_and_group.first);
                    // generators of one group agree on a solved row (_row_gen_v_conflicts)
                    if(std::isfinite(val)) v_set(gen_and_group.second) = val;
                }
            }
            algo.set_voltage_control_v_set(v_set);
        }

        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _prepare_gen_v_vc() {
            _gen_v_vc_pairs_.clear();
            _gen_v_vc_n_groups_ = 0;
            if(sbus_policy_.gen_v.rows() == 0) return;   // never set: the grid's own
            const IntVect vc_group = _gen_v_vc_group(/*include_pinned=*/true);
            for(Eigen::Index g = 0; g < vc_group.size(); ++g){
                if(vc_group[g] >= 0) _gen_v_vc_pairs_.push_back(std::make_pair(static_cast<int>(g), vc_group[g]));
            }
            if(!_gen_v_vc_pairs_.empty())
                _gen_v_vc_n_groups_ = static_cast<Eigen::Index>(_grid_model.get_ac_voltage_control_plan().controllers().n_groups());
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _prepare_gen_v_vc() {}

        // Per generator, the voltage-control group whose v_set its gen_v is (the group
        // holding the bus it regulates), -1 for none. `include_pinned` false also
        // answers -1 for a group holding an SVC or an hvdc station (see
        // get_gen_v_vc_row). Same writers as set_vm (vm_writers_by_bus), keyed on the
        // REGULATED bus.
        IntVect _gen_v_vc_group(bool include_pinned = false) const {
            IntVect res = IntVect::Constant(
                static_cast<Eigen::Index>(_grid_model.get_generators_as_data().nb()), -1);
            if(!_algo.ac_solver_used()) return res;
            const VoltageControlSolverData & ctrl = _grid_model.get_ac_voltage_control_plan().controllers();
            if(ctrl.n_groups() == 0) return res;
            std::map<int, int> group_of_bus;                 // regulated solver bus -> group
            for(int grp = 0; grp < ctrl.n_groups(); ++grp){
                bool pinned = false;
                for(int j = ctrl.grp_start(grp); j < ctrl.grp_start(grp) + ctrl.grp_count(grp); ++j){
                    if(ctrl.kind(j) != VoltageControlSolverData::GEN) pinned = true;
                }
                if(pinned && !include_pinned) continue;
                group_of_bus[ctrl.reg_bus(grp)] = grp;
            }
            std::map<int, std::vector<int> > gens_of_bus;
            _grid_model.get_generators().vm_writers_by_bus(active_layout().id_me_to_solver, gens_of_bus);
            for(const auto & bus_and_gens : gens_of_bus){
                const std::map<int, int>::const_iterator it = group_of_bus.find(bus_and_gens.first);
                if(it == group_of_bus.end()) continue;
                for(int gen_id : bus_and_gens.second) res[static_cast<Eigen::Index>(gen_id)] = it->second;
            }
            return res;
        }

        // ----- modify_gen_v: set-points that contradict each other ----------------
        // A bus has ONE voltage magnitude. Two generators regulating it with two
        // different targets state two constraints it cannot both satisfy, and set_vm
        // resolves that by applying whichever generator it visits last -- an answer to
        // a question nobody asked, and silent. A row whose set-points do that has an
        // input that cannot be met, so it is not solved: it is reported exactly like a
        // row skipped before the solver ever ran (converged_mask 0, a GRID /
        // NOT_SIMULATED violation, zero voltages, and so no gradient).
        //
        // Only the groups matter -- the buses with more than one regulating generator,
        // which most grids have few of and many have none -- and which generator writes
        // which bus does not depend on the row. So the groups are found once per
        // compute() and each row then only compares numbers.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _prepare_gen_v_constraints() {
            _gen_v_constraints_.clear();
            if(sbus_policy_.gen_v.rows() == 0) return;   // never set: the grid's own targets
            const SolverBusIdVect & me_to_solver = active_layout().id_me_to_solver;

            // which generators write which bus ...
            std::map<int, std::vector<int> > gens_of_bus;
            _grid_model.get_generators().vm_writers_by_bus(me_to_solver, gens_of_bus);
            // ... and what the elements whose set-point modify_gen_v CANNOT move ask of
            // those same buses. A converter station and a voltage-mode SVC carry their
            // own target and the batch has no per-row setter for either, so a row that
            // moves a generator sharing their bus contradicts a value it cannot change.
            std::map<int, real_type> fixed_of_bus;
            _grid_model.get_svcs().vm_targets_by_bus(me_to_solver, fixed_of_bus);
            _grid_model.get_dclines().vm_targets_by_bus(me_to_solver, fixed_of_bus);

            for(const auto & bus_and_gens : gens_of_bus){
                const std::map<int, real_type>::const_iterator fixed = fixed_of_bus.find(bus_and_gens.first);
                const bool has_fixed = (fixed != fixed_of_bus.end());
                // a bus one generator writes, with nothing else on it, can be given
                // whatever this row likes: there is nothing to contradict
                if(bus_and_gens.second.size() < 2 && !has_fixed) continue;
                GenVConstraint c;
                c.bus_solver = bus_and_gens.first;
                c.gens = bus_and_gens.second;
                c.has_fixed = has_fixed;
                c.fixed_vm = has_fixed ? fixed->second : 0.;
                _gen_v_constraints_.push_back(c);
            }
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _prepare_gen_v_constraints() {}

        // The groups are those of the base placement: a generator the row takes off the
        // bus (disconnected, or moved elsewhere -- see _maybe_resolve_topology) no longer
        // writes it, so it has nothing to agree with there. Where it lands is checked
        // against what stands there once, at compute() (_maybe_prepare_gen_contingency).
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _row_gen_v_conflicts(size_t i) const {
            if(_gen_v_constraints_.empty()) return false;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(row >= sbus_policy_.gen_v.rows()) return false;
            for(size_t c = 0; c < _gen_v_constraints_.size(); ++c){
                const GenVConstraint & con = _gen_v_constraints_[c];
                // whatever else is on the bus pins the value; otherwise the generators
                // left there only have to agree among themselves
                bool has_ref = con.has_fixed;
                real_type ref = con.fixed_vm;
                for(int gen_id : con.gens){
                    if(_gen_off_in_row(i, gen_id)) continue;
                    const real_type vm = sbus_policy_.gen_v(row, gen_id);
                    if(!has_ref){
                        ref = vm;
                        has_ref = true;
                        continue;
                    }
                    if(std::abs(vm - ref) > BaseConstants::_tol_equal_float) return true;
                }
            }
            return false;
        }

        // set_vm writes every generator at the bus of the base placement, the last one
        // visited winning: on a shared bus row i takes one of them off, that may be the
        // one that left. The bus is put back at what stays there asks for -- the
        // generators left (they agree, _row_gen_v_conflicts), or the element whose
        // set-point no row moves. A bus they all leave is the PV -> PQ path's.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _reseat_shared_buses(size_t i, const Eigen::Ref<const RealVect> & target_vm_pu_row, CplxVect & V) const {
            for(const GenVConstraint & con : _gen_v_constraints_){
                bool some_off = false;
                int live = -1;
                for(int gen_id : con.gens){
                    if(_gen_off_in_row(i, gen_id)) some_off = true;
                    else if(live < 0) live = gen_id;
                }
                if(!some_off) continue;
                const real_type target = con.has_fixed ? con.fixed_vm
                                                       : (live >= 0 ? target_vm_pu_row(live) : std::numeric_limits<real_type>::quiet_NaN());
                const int b = con.bus_solver;
                if(!std::isfinite(target) || b < 0 || b >= V.size()) continue;
                const real_type vm = std::abs(V(b));
                V(b) = vm > 0. ? V(b) * (target / vm) : cplx_type(target, 0.);
            }
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _row_gen_v_conflicts(size_t) const { return false; }

        // ----- sbus_policy_.gen_v, generically (2-way, same split as above): feeds
        // BaseBatchSolverSynch's generic (SbusPolicy-agnostic) DC fast-path magnitude
        // reconstruction (_dc_gen_v_ / _dc_vm_row_grid) -- empty wherever
        // modify_gen_v does not exist (SbusPolicy::NOOP, ie ContingencyAnalysis) or
        // was never called. -----
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        const RealMat & _sbus_gen_v() const { return sbus_policy_.gen_v; }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        RealMat _sbus_gen_v() const { return RealMat(); }

        // ----- per-row "which branches did THIS row itself disconnect" ------------
        // used by _record_row_violations below to exclude a row's own disconnected
        // branches from that row's current-limit checks (they still look "connected"
        // via get_status_global(), only the Ybus coefficients were edited). Two
        // sources of truth, per instantiation: ContingencyAnalysis has a cache built
        // once per compute() from _branches_out_vect() (the sparse, add_n1-driven
        // representation); ScenarioSweep has no such cache (no li_defaults at all)
        // and instead reads the row directly off line_mask/trafo_mask.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<int> _row_skip_branch_ids(size_t i) const {
            return (i < _li_defaults_vect_cache_.size()) ? _li_defaults_vect_cache_[i] : std::vector<int>();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        std::vector<int> _row_skip_branch_ids(size_t i) const {
            return ybus_policy_.branch_ids_for_row(static_cast<Eigen::Index>(i), n_line_);
        }

        // ----- the row as the outer loops see it (see OuterRowView) --------------
        // Everything this row changes compared with the grid, read off the same per-row
        // data the row's Sbus / Ybus are built from: its injections and generator outages
        // (SbusPolicy::Vary), its topological action, its contingency's branches. Dense per
        // element kind, and only for the kinds the row changes; `v` is the worker's own,
        // refilled for every row (its vectors keep their capacity).
        void _fill_outer_row_view(size_t i, const std::vector<int> * masked, OuterRowView & v) const {
            v.layout = &active_layout();
            v.masked = (masked != nullptr && !masked->empty()) ? masked : nullptr;
            v.resume = false;
            v.resume_from = nullptr;
            v.moved_trafo_blocks = nullptr;
            v.lost_p_mw = i < _row_lost_p_mw_.size() ? _row_lost_p_mw_[i] : std::numeric_limits<real_type>::quiet_NaN();
            v.clear_overrides();
            _fill_outer_row_injections(i, v);
            _fill_outer_row_trips(i, v);
            _fill_outer_row_topo(i, v);
            _fill_outer_row_branches(i, v);
        }
        // the generators a ContingencyAnalysis contingency trips (a ScenarioSweep's come
        // with its injections above) and the loads a row's contingency trips, both
        // classes: out of the row, as a topological action takes them out
        void _fill_outer_row_trips(size_t i, OuterRowView & v) const {
            const std::set<int> * gens = _row_tripped_gens_ca(i);
            if(gens != nullptr){
                const int nb_gen = _grid_model.get_generators().nb();
                if(v.gens.off.empty()) v.gens.off.assign(static_cast<std::size_t>(nb_gen), 0);
                for(int g : *gens) if(g >= 0 && g < nb_gen) v.gens.off[static_cast<std::size_t>(g)] = 1;
            }
            const std::vector<int> loads = _row_tripped_loads(i);
            if(!loads.empty()){
                const int nb_load = _grid_model.get_loads().nb();
                if(v.loads.off.empty()) v.loads.off.assign(static_cast<std::size_t>(nb_load), 0);
                for(int l : loads) if(l >= 0 && l < nb_load) v.loads.off[static_cast<std::size_t>(l)] = 1;
            }
        }
        // where the row's loops start: a row chained on the previous one (TimeSeries) from its
        // equations -- the buses its loops froze, the taps they moved, ... (the loops undo what
        // this row does not call for); a row seeded from the "n" solve (init_from_n_powerflow,
        // a TimeSeries' first one included) from that solve's; any other from the grid
        void _set_outer_row_resume(size_t i, OuterRowView & v) const {
            const bool chained = INIT == BatchInitKind::FromPreviousStep && i > 0;
            v.resume_from = chained ? nullptr : _seeded_row_resume_point();
            v.resume = chained || v.resume_from != nullptr;
        }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _fill_outer_row_injections(size_t i, OuterRowView & v) const {
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            auto fill = [row](const RealMat & mat, std::vector<real_type> & out){
                if(mat.rows() == 0 || row >= mat.rows()) return;
                out.resize(static_cast<std::size_t>(mat.cols()));
                for(Eigen::Index c = 0; c < mat.cols(); ++c) out[static_cast<std::size_t>(c)] = mat(row, c);
            };
            fill(sbus_policy_.gen_p, v.gens.target_p);
            fill(sbus_policy_.gen_v, v.gens.target_vm);
            fill(sbus_policy_.load_q, v.loads.target_q);
            const int nb_gen = _grid_model.get_generators().nb();
            if(sbus_policy_.has_gen_off()){
                v.gens.off.assign(static_cast<std::size_t>(nb_gen), 0);
                for(int g = 0; g < nb_gen; ++g) v.gens.off[static_cast<std::size_t>(g)] = _gen_off_in_row(i, g) ? 1 : 0;
            }
            if(i < sbus_policy_.topo_gens_on.size() && !sbus_policy_.topo_gens_on[i].empty()){
                v.gens.bus_solver.assign(static_cast<std::size_t>(nb_gen), OuterRowView::Units::GRID_BUS);
                for(const auto & gen_on : sbus_policy_.topo_gens_on[i]){
                    if(gen_on.gen_id >= 0 && gen_on.gen_id < nb_gen) v.gens.bus_solver[static_cast<std::size_t>(gen_on.gen_id)] = gen_on.bus_solver;
                    // a placed generator is on in its row (topo_gen_off may list its base bus)
                    if(!v.gens.off.empty()) v.gens.off[static_cast<std::size_t>(gen_on.gen_id)] = 0;
                }
            }
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _fill_outer_row_injections(size_t, OuterRowView &) const {}
        // loads and storage units a topological action moves or disconnects (generators come
        // with the injections above)
        void _fill_outer_row_topo(size_t i, OuterRowView & v) const {
            if(i >= _row_topo_.size()) return;
            const RowTopoPlan & plan = _row_topo_[i];
            const SolverBusIdVect & to_solver = active_layout().id_me_to_solver;
            auto place = [&to_solver](const std::vector<TopoElPlacement> & placements, int nb, OuterRowView::Units & u){
                if(placements.empty()) return;
                // the loads the row's contingency trips may have set `off` already
                if(u.off.empty()) u.off.assign(static_cast<std::size_t>(nb), 0);
                if(u.bus_solver.empty()) u.bus_solver.assign(static_cast<std::size_t>(nb), OuterRowView::Units::GRID_BUS);
                for(const TopoElPlacement & p : placements){
                    if(p.el_id < 0 || p.el_id >= nb) continue;
                    if(p.row_bus_me < 0){ u.off[static_cast<std::size_t>(p.el_id)] = 1; continue; }
                    const int b = to_solver[p.row_bus_me].cast_int();
                    u.bus_solver[static_cast<std::size_t>(p.el_id)] = b == BaseConstants::_deactivated_bus_id ? -1 : b;
                }
            };
            place(plan.loads, _grid_model.get_loads().nb(), v.loads);
            place(plan.storages, _grid_model.get_storages().nb(), v.storages);
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _fill_outer_row_branches(size_t i, OuterRowView & v) const {
            const int nb_lines = static_cast<int>(n_line_);
            const int nb_trafos = static_cast<int>(n_trafos_);
            for(int br_id : _row_skip_branch_ids(i)){
                std::vector<char> & off = br_id < nb_lines ? v.lines.off : v.trafos.off;
                const int nb = br_id < nb_lines ? nb_lines : nb_trafos;
                const int el = br_id < nb_lines ? br_id : br_id - nb_lines;
                if(el < 0 || el >= nb) continue;
                if(off.empty()) off.assign(static_cast<std::size_t>(nb), 0);
                off[static_cast<std::size_t>(el)] = 1;
            }
            if(i < ybus_policy_.topo_branches_moved.size()){
                for(const auto & placement : ybus_policy_.topo_branches_moved[i]){
                    const bool is_line = placement.branch_id < nb_lines;
                    OuterRowView::Branches & b = is_line ? v.lines : v.trafos;
                    const int nb = is_line ? nb_lines : nb_trafos;
                    const int el = is_line ? placement.branch_id : placement.branch_id - nb_lines;
                    if(el < 0 || el >= nb) continue;
                    if(b.bus1_solver.empty()){
                        b.bus1_solver.assign(static_cast<std::size_t>(nb), OuterRowView::Units::GRID_BUS);
                        b.bus2_solver.assign(static_cast<std::size_t>(nb), OuterRowView::Units::GRID_BUS);
                    }
                    b.bus1_solver[static_cast<std::size_t>(el)] = placement.bus1_solver < 0 ? -1 : placement.bus1_solver;
                    b.bus2_solver[static_cast<std::size_t>(el)] = placement.bus2_solver < 0 ? -1 : placement.bus2_solver;
                }
            }
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _fill_outer_row_branches(size_t, OuterRowView &) const {}

        // ----- per-row violation recording ---------------------------------------
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _record_row_violations(size_t i, const CplxVect & V, const std::vector<int> * masked_ids){
            const bool ac_solver_used = _algo.ac_solver_used();
            const real_type sn_mva = _grid_model.get_sn_mva();

            std::vector<int> skip_lines, skip_trafos;
            for(int br_id : _row_skip_branch_ids(i)){
                if(static_cast<size_t>(br_id) < n_line_) skip_lines.push_back(br_id);
                else skip_trafos.push_back(static_cast<int>(br_id - static_cast<int>(n_line_)));
            }
            const std::vector<int> * masked_ids_use = (masked_ids != nullptr && !masked_ids->empty()) ? masked_ids : nullptr;

            batch_sweep_detail::check_bus_voltage_violations(V, active_layout().id_me_to_solver, _grid_model.get_bus_vmin_kv(), _grid_model.get_bus_vmax_kv(),
                                         _grid_model.get_bus_vn_kv(), _grid_model.get_substations(),
                                         _violation_threshold_, _violation_rel_tol_, masked_ids_use, _violations[i]);
            const std::vector<BranchBusOverride> * overrides =
                (i < _row_branch_overrides_.size() && !_row_branch_overrides_[i].empty()) ? &_row_branch_overrides_[i] : nullptr;
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, _violation_rel_tol_, skip_lines, _violations[i],
                                     overrides, 0);
            // the taps this row's outer loops left its transformers at (_record_row_trafo_blocks)
            const auto * blocks = (i < _row_trafo_blocks_.size() && !_row_trafo_blocks_[i].empty())
                                  ? &_row_trafo_blocks_[i] : nullptr;
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, _violation_rel_tol_, skip_trafos, _violations[i],
                                     overrides, n_line_, blocks);
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _record_row_violations_dispatch(size_t i, const CplxVect & V, const std::vector<int> * masked_ids) {
            _record_row_violations(i, V, masked_ids);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _record_row_violations_dispatch(size_t, const CplxVect &, const std::vector<int> *) {}

        // ----- per-row reactive-capability recording ------------------------------
        // Whether THIS row disconnected that generator: a ScenarioSweep generator
        // contingency (its mask or its topological action), a ContingencyAnalysis
        // contingency tripping it; always false where Ybus does not vary.
        template<class Y = YbusPolicy, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _gen_off_in_row(size_t i, int gen_id) const {
            return sbus_policy_.gen_off_in(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(gen_id));
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && Y::supports_contingency, int>::type = 0>
        bool _gen_off_in_row(size_t i, int gen_id) const { return ybus_policy_.gen_off_in_row(i, gen_id); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && !Y::supports_contingency, int>::type = 0>
        bool _gen_off_in_row(size_t, int) const { return false; }
        // ... whether some row disconnects a generator at all
        template<class Y = YbusPolicy, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _any_gen_off() const { return sbus_policy_.has_gen_off(); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && Y::supports_contingency, int>::type = 0>
        bool _any_gen_off() const { return ybus_policy_.has_gen_off(); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && !Y::supports_contingency, int>::type = 0>
        bool _any_gen_off() const { return false; }
        // ... whether some row's topological action puts a base-off generator back
        // (ScenarioSweep only), and the generators each row puts back
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _has_gens_on() const { return sbus_policy_.has_gens_on(); }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _has_gens_on() const { return false; }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        const std::vector<std::vector<::ls2g::SbusPolicy::Vary::GenOn> > & _topo_gens_on_rows() const { return sbus_policy_.topo_gens_on; }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        const std::vector<std::vector<::ls2g::SbusPolicy::Vary::GenOn> > & _topo_gens_on_rows() const {
            static const std::vector<std::vector<::ls2g::SbusPolicy::Vary::GenOn> > none;
            return none;
        }

        // Whether THIS row disconnected that load / shunt: its contingency
        // (ContingencyAnalysis) or its set_contingency_loads / set_contingency_shunts mask
        // (ScenarioSweep; a load its topological action takes off is _load_leaves_base_bus')
        template<class Y = YbusPolicy, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _load_off_in_row(size_t i, int load_id) const {
            return sbus_policy_.load_off_in(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(load_id));
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && Y::supports_contingency, int>::type = 0>
        bool _load_off_in_row(size_t i, int load_id) const { return ybus_policy_.load_off_in_row(i, load_id); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && !Y::supports_contingency, int>::type = 0>
        bool _load_off_in_row(size_t, int) const { return false; }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool _shunt_off_in_row(size_t i, int shunt_id) const { return ybus_policy_.shunt_off_in_row(i, shunt_id); }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        bool _shunt_off_in_row(size_t, int) const { return false; }
        // the loads / shunts row i's contingency trips (sorted ids), for the few places
        // that walk them rather than test them one by one
        template<class Y = YbusPolicy, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        std::vector<int> _row_tripped_loads(size_t i) const {
            std::vector<int> res;
            const auto & mask = sbus_policy_.load_off;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(mask.rows() == 0 || row >= mask.rows()) return res;
            for(Eigen::Index l = 0; l < mask.cols(); ++l) if(mask(row, l)) res.push_back(static_cast<int>(l));
            return res;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && Y::supports_contingency, int>::type = 0>
        std::vector<int> _row_tripped_loads(size_t i) const {
            if(i >= ybus_policy_.li_rows.size()) return std::vector<int>();
            return std::vector<int>(ybus_policy_.li_rows[i].loads.begin(), ybus_policy_.li_rows[i].loads.end());
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!S::supports_vary && !Y::supports_contingency, int>::type = 0>
        std::vector<int> _row_tripped_loads(size_t) const { return std::vector<int>(); }
        template<class Y = YbusPolicy, class S = SbusPolicy, typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        std::vector<int> _row_tripped_shunts(size_t i) const { return ybus_policy_.shunt_ids_for_row(static_cast<Eigen::Index>(i)); }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<int> _row_tripped_shunts(size_t i) const {
            if(i >= ybus_policy_.li_rows.size()) return std::vector<int>();
            return std::vector<int>(ybus_policy_.li_rows[i].shunts.begin(), ybus_policy_.li_rows[i].shunts.end());
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        std::vector<int> _row_tripped_shunts(size_t) const { return std::vector<int>(); }
        // the generators row i's ContingencyAnalysis contingency trips (a ScenarioSweep
        // reads its mask through _gen_off_in_row)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::set<int> * _row_tripped_gens_ca(size_t i) const {
            if(i >= ybus_policy_.li_rows.size() || ybus_policy_.li_rows[i].gens.empty()) return nullptr;
            return &ybus_policy_.li_rows[i].gens;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        const std::set<int> * _row_tripped_gens_ca(size_t) const { return nullptr; }

        // the generators row i has out, ascending, in `out`: exactly those _gen_off_in_row says
        // so of, listed. A ContingencyAnalysis row's are its contingency's trip list (asking
        // _gen_off_in_row for every generator would be a set lookup each); elsewhere the masks
        // are scanned.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _row_gens_off(size_t i, int nb_gen, std::vector<int> & out) const {
            out.clear();
            const std::set<int> * gens = _row_tripped_gens_ca(i);
            if(gens == nullptr) return;
            for(int gen_id : *gens) if(gen_id >= 0 && gen_id < nb_gen) out.push_back(gen_id);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _row_gens_off(size_t i, int nb_gen, std::vector<int> & out) const {
            out.clear();
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id) if(_gen_off_in_row(i, gen_id)) out.push_back(gen_id);
        }

        // This row's own active set-point for that generator, in MW -- the number the
        // distributed slack's share is added ON TOP of (see GenPCheck.hpp). Exactly what
        // SbusPolicy::Vary::fill_row stamped into this row's injection: its own gen_p row
        // where modify_gen_p was given one, the grid's target otherwise. Where the
        // injection does not vary at all (ContingencyAnalysis), the grid's target IS the
        // row's.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _gen_target_p_in_row(size_t i, int gen_id) const {
            const auto & mat = sbus_policy_.gen_p;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(mat.rows() > 0 && row < mat.rows() && gen_id < mat.cols()) return mat(row, gen_id);
            return _grid_target_p(gen_id);
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _gen_target_p_in_row(size_t, int gen_id) const { return _grid_target_p(gen_id); }

        // This row's own voltage target for that generator, in pu: its own gen_v row where
        // modify_gen_v was given one, the grid's target_vm_pu otherwise. For a PQ machine
        // the solve never reads it (modify_gen_v re-seeds the PV buses only); the PQ -> PV
        // release check does, as "the target it would hold if released" (see
        // GenPvReleaseCheck.hpp) -- so a sweep that moves a pinned machine's target moves
        // what it is checked against.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _gen_target_vm_in_row(size_t i, int gen_id) const {
            const auto & mat = sbus_policy_.gen_v;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(mat.rows() > 0 && row < mat.rows() && gen_id < mat.cols()) return mat(row, gen_id);
            return _grid_model.get_generators().get_target_vm_pu(gen_id);
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _gen_target_vm_in_row(size_t, int gen_id) const {
            return _grid_model.get_generators().get_target_vm_pu(gen_id);
        }

        real_type _grid_target_p(int gen_id) const {
            const Eigen::Ref<const RealVect> tgt = _grid_model.get_gen_target_p();
            return (gen_id >= 0 && gen_id < tgt.size()) ? tgt(gen_id) : 0.;
        }

        // The active power imbalance a DC row leaves for the slack machines to make up
        // (MW), the same `-sum(Pbus)` LSGrid::_fill_bus_mismatch_dc shares out -- read off
        // the injection the row was actually solved with (empty means the DC entry point
        // fell back to the member dc_cache_.inj, see compute_one_powerflow).
        // Summed over the SOLVED buses only: a bus `masked` (stranded by this row's
        // contingency, sorted solver ids, may be nullptr) is out of the DC balance
        // (BaseDCAlgo builds its right-hand side over the live buses only).
        real_type _dc_imbalance_mw(const Eigen::Ref<const CplxVect> & sbus_solver,
                                   const std::vector<int> * masked) const {
            const real_type sn_mva = _grid_model.get_sn_mva();
            const bool has_masked = (masked != nullptr) && !masked->empty();
            real_type sum = 0.;
            if(sbus_solver.size() > 0){
                for(Eigen::Index b = 0; b < sbus_solver.size(); ++b){
                    if(has_masked && std::binary_search(masked->begin(), masked->end(), static_cast<int>(b))) continue;
                    sum += std::real(sbus_solver(b));
                }
            } else {
                for(Eigen::Index b = 0; b < dc_cache_.inj.size(); ++b){
                    if(has_masked && std::binary_search(masked->begin(), masked->end(), static_cast<int>(b))) continue;
                    sum += dc_cache_.inj(b);
                }
            }
            return -sum * sn_mva;
        }

        // This row's masked (stranded) solver buses, or nullptr when it masks none.
        // Only the "handle disconnected grid" mode ever hands _li_masked to the
        // algorithm; without it a contingency that strands a bus makes the row skipped
        // outright, so nothing is masked (see _adjoint_identity_rows, same reasoning).
        const std::vector<int> * _row_masked_ids(size_t i) const {
            if(!_mask_mode()) return nullptr;
            if(i >= _li_masked.size()) return nullptr;
            return _li_masked[i].empty() ? nullptr : &_li_masked[i];
        }

        // the (generator, solver bus) pairs row i's topological action places, restricted
        // to the generators that pin their bus there (a non-regulating one is an
        // injection, out of the reactive check); empty where the row places none
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        std::vector<std::pair<int, int> > _row_gens_placed(size_t i) const {
            std::vector<std::pair<int, int> > res;
            if(i >= sbus_policy_.topo_gens_on.size()) return res;
            const auto & generators = _grid_model.get_generators();
            for(const auto & gen_on : sbus_policy_.topo_gens_on[i]){
                if(!generators.would_be_local_voltage_controller(gen_on.gen_id)) continue;
                res.push_back(std::make_pair(gen_on.gen_id, gen_on.bus_solver));
            }
            return res;
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        std::vector<std::pair<int, int> > _row_gens_placed(size_t) const { return std::vector<std::pair<int, int> >(); }
        // every (generator, solver bus) pair row i's topological action places, whether or
        // not the generator pins its bus there (see _record_row_gen_results)
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        std::vector<std::pair<int, int> > _row_gens_placed_all(size_t i) const {
            std::vector<std::pair<int, int> > res;
            if(i >= sbus_policy_.topo_gens_on.size()) return res;
            for(const auto & gen_on : sbus_policy_.topo_gens_on[i]) res.push_back(std::make_pair(gen_on.gen_id, gen_on.bus_solver));
            return res;
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        std::vector<std::pair<int, int> > _row_gens_placed_all(size_t) const { return std::vector<std::pair<int, int> >(); }

        // ---- row i's topological action, as the slack pre-pass reads it (see
        // _prepare_slack_redistribution; nothing outside ScenarioSweep) ----------------
        // whether the action places that generator on a bus (moved, or reactivated)
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _gen_placed_in_row(size_t i, int gen_id) const {
            if(i >= sbus_policy_.topo_gens_on.size()) return false;
            for(const auto & gen_on : sbus_policy_.topo_gens_on[i]) if(gen_on.gen_id == gen_id) return true;
            return false;
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _gen_placed_in_row(size_t, int) const { return false; }
        // whether the action takes that load / storage unit off the bus it holds in the
        // base grid (disconnected, or moved)
        // (or its contingency trips it: set_contingency_loads, a ContingencyAnalysis load)
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _load_leaves_base_bus(size_t i, int load_id) const {
            if(_load_off_in_row(i, load_id)) return true;
            if(i >= sbus_policy_.topo_loads_off.size()) return false;
            const std::vector<int> & off = sbus_policy_.topo_loads_off[i];
            return std::find(off.begin(), off.end(), load_id) != off.end();
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _load_leaves_base_bus(size_t i, int load_id) const { return _load_off_in_row(i, load_id); }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _storage_leaves_base_bus(size_t i, int storage_id) const {
            if(i >= sbus_policy_.topo_storages_off.size()) return false;
            const std::vector<int> & off = sbus_policy_.topo_storages_off[i];
            return std::find(off.begin(), off.end(), storage_id) != off.end();
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _storage_leaves_base_bus(size_t, int) const { return false; }

        // What the action changes in the active balance, generator convention, as power
        // LOST (what the remaining units make up): for each generator, load or storage
        // unit it places, disconnects or reactivates, what the base balance counted for
        // it minus what the row's main component gets from it. A unit moved inside the
        // main component nets out; one disconnected, or put on a bus the row masks, is
        // lost; one reactivated is a gain. `in_island(bus_solver)`: the row masks that bus.
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _row_topo_lost_mw(size_t i, const InIsland & in_island) const {
            real_type lost = 0.;
            const auto & generators = _grid_model.get_generators();
            const std::vector<bool> & gen_status = generators.get_status();
            const GlobalBusIdVect & gen_bus = generators.get_bus_id();
            if(i < sbus_policy_.topo_gens_on.size()){
                for(const auto & gen_on : sbus_policy_.topo_gens_on[i]){
                    const int gen_id = gen_on.gen_id;
                    const real_type p = _gen_target_p_in_row(i, gen_id);
                    const bool base_on = gen_status[gen_id] &&
                                         gen_bus(gen_id).cast_int() != BaseConstants::_deactivated_bus_id;
                    if(base_on) lost += p;
                    if(!in_island(gen_on.bus_solver)) lost -= p;
                }
            }
            // a load / storage unit leaving its base bus was connected in the base grid
            // (see _maybe_resolve_topology): its consumption leaves the balance, and comes
            // back where the row places it in the main component
            if(i < sbus_policy_.topo_loads_off.size()){
                for(int load_id : sbus_policy_.topo_loads_off[i]) lost -= _load_target_p_in_row(i, load_id);
            }
            if(i < sbus_policy_.topo_loads_on.size()){
                for(const auto & load_bus : sbus_policy_.topo_loads_on[i]){
                    if(!in_island(load_bus.second)) lost += _load_target_p_in_row(i, load_bus.first);
                }
            }
            const Eigen::Ref<const RealVect> storage_p = _grid_model.get_storages().get_target_p();
            if(i < sbus_policy_.topo_storages_off.size()){
                for(int storage_id : sbus_policy_.topo_storages_off[i]) lost -= storage_p(storage_id);
            }
            if(i < sbus_policy_.topo_storages_on.size()){
                for(const auto & storage_bus : sbus_policy_.topo_storages_on[i]){
                    if(!in_island(storage_bus.second)) lost += storage_p(storage_bus.first);
                }
            }
            return lost;
        }
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _row_topo_lost_mw(size_t, const InIsland &) const { return 0.; }

        // The participating units the action places in the main component, on the bus it
        // gives them (the ones it moves are left out of their base bus by the caller).
        // The list is re-sorted (family, then id) when one is added: each family's
        // per-row list of moved set-points has to come out sorted, see _row_target_p.
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _append_topo_gen_participants(size_t i, const InIsland & in_island,
                                           std::vector<slack_redistribution::Participant> & units) const {
            if(i >= sbus_policy_.topo_gens_on.size() || sbus_policy_.topo_gens_on[i].empty()) return;
            const size_t before = units.size();
            for(const auto & gen_on : sbus_policy_.topo_gens_on[i]){
                if(in_island(gen_on.bus_solver)) continue;
                slack_redistribution::append_participant(
                    _grid_model.get_generators(), slack_redistribution::UnitKind::GENERATOR, gen_on.gen_id,
                    gen_on.bus_me, _gen_target_p_in_row(i, gen_on.gen_id), units);
            }
            if(units.size() != before) _sort_participants(units);
        }
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _append_topo_gen_participants(size_t, const InIsland &, std::vector<slack_redistribution::Participant> &) const {}
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _append_topo_storage_participants(size_t i, const InIsland & in_island,
                                               std::vector<slack_redistribution::Participant> & units) const {
            if(i >= sbus_policy_.topo_storages_on.size() || sbus_policy_.topo_storages_on[i].empty()) return;
            const auto & storages = _grid_model.get_storages();
            const auto solver_to_me = active_layout().id_solver_to_me.as_eigen();
            const size_t before = units.size();
            for(const auto & storage_bus : sbus_policy_.topo_storages_on[i]){
                if(in_island(storage_bus.second)) continue;
                slack_redistribution::append_participant(
                    storages, slack_redistribution::UnitKind::STORAGE, storage_bus.first,
                    solver_to_me[storage_bus.second], -storages.get_target_p()(storage_bus.first), units);
            }
            if(units.size() != before) _sort_participants(units);
        }
        template<class InIsland, class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _append_topo_storage_participants(size_t, const InIsland &, std::vector<slack_redistribution::Participant> &) const {}
        static void _sort_participants(std::vector<slack_redistribution::Participant> & units){
            std::stable_sort(units.begin(), units.end(),
                             [](const slack_redistribution::Participant & a, const slack_redistribution::Participant & b){
                                 if(a.kind != b.kind) return a.kind < b.kind;
                                 return a.el_id < b.el_id; });
        }


        // MUST be called while the row's own state is still installed on the algorithm
        // -- its Ybus edits, its masked buses, its PV pinning -- for the same reason
        // _maybe_store_jacobian must: what is read here (the per-bus mismatch, the
        // controllers' reactive output, the bus angles) belongs to the system THIS row
        // solved, and the `algo` handed in is the one that solved it (the member one, or
        // this thread's). `V_solver` is this row's converged complex voltage: a
        // voltage-mode SVC's capability is a susceptance range, so what it is worth in MVAr
        // depends on it.
        // MUST be called while the row's own state is still installed on the algorithm,
        // like _record_row_physical: the mismatch, the slack the solve absorbed and the
        // controllers' reactive output are the ones of the system THIS row solved. `Ybus`
        // is the row's matrix (its branch edits in place), read only when the algorithm
        // leaves no mismatch behind.
        void _record_row_gen_results(size_t i, AlgorithmSelector & algo,
                                     const Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const Eigen::Ref<const CplxVect> & V_solver,
                                     const Eigen::Ref<const RealVect> & slack_weights,
                                     const Eigen::Ref<const CplxVect> & sbus_solver){
            if(!_compute_gen_results_) return;
            if(static_cast<Eigen::Index>(i) >= _gen_results_.rows()) return;
            const real_type sn_mva = _grid_model.get_sn_mva();
            const bool ac = algo.ac_solver_used();
            const std::vector<int> * masked = _row_masked_ids(i);
            const auto & generators = _grid_model.get_generators();
            const int nb_gen = static_cast<int>(generators.nb());
            const int nb_storage = static_cast<int>(_grid_model.get_storages().nb());
            const SolverBusIdVect & me_to_solver = active_layout().id_me_to_solver;

            // the generators this row has out -- its contingency, and the ones standing on
            // a bus it strands -- and the units its slack pre-pass took out of the slack
            std::vector<bool> gen_off(static_cast<std::size_t>(nb_gen), false);
            std::vector<bool> no_share(static_cast<std::size_t>(nb_gen), false);
            RealVect target_p(nb_gen);
            const GlobalBusIdVect & gen_buses = generators.get_bus_id();
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                bool off = _gen_off_in_row(i, gen_id);
                if(!off && masked != nullptr){
                    const int bus_me = gen_buses(gen_id).cast_int();
                    const int bus_solver = bus_me == BaseConstants::_deactivated_bus_id ?
                                           BaseConstants::_deactivated_bus_id : me_to_solver[bus_me].cast_int();
                    off = bus_solver != BaseConstants::_deactivated_bus_id &&
                          std::binary_search(masked->begin(), masked->end(), bus_solver);
                }
                gen_off[gen_id] = off;
                no_share[gen_id] = _row_takes_no_share(i, ViolationElementType::GENERATOR, gen_id);
                target_p(gen_id) = _row_target_p(i, ViolationElementType::GENERATOR, gen_id);
            }
            // the generators this row's action places (moved, or reactivated): off their
            // base bus, which is what gen_off said of a moved one, but on the row's --
            // unless that bus is one the row masks
            const std::vector<std::pair<int, int> > gen_placed = _row_gens_placed_all(i);
            const GlobalBusIdVect & solver_to_me = active_layout().id_solver_to_me;
            std::vector<std::pair<int, int> > gen_placed_me;
            gen_placed_me.reserve(gen_placed.size());
            for(const auto & gen_bus : gen_placed){
                const bool in_island = gen_bus.second < 0 ||
                                       (masked != nullptr && std::binary_search(masked->begin(), masked->end(), gen_bus.second));
                gen_off[gen_bus.first] = in_island;
                if(!in_island) gen_placed_me.push_back({gen_bus.first, solver_to_me[static_cast<size_t>(gen_bus.second)].cast_int()});
            }
            std::vector<bool> storage_no_share;
            for(int storage_id = 0; storage_id < nb_storage; ++storage_id){
                if(!_row_takes_no_share(i, ViolationElementType::STORAGE, storage_id)) continue;
                if(storage_no_share.empty()) storage_no_share.assign(static_cast<std::size_t>(nb_storage), false);
                storage_no_share[storage_id] = true;
            }

            // what each bus produced beyond its setpoints: the same two sources
            // LSGrid::_fill_bus_mismatch_ac / _dc read, for this row
            RealVect active_mismatch, reactive_mismatch, ctrl_q;
            IntVect ctrl_kind, ctrl_elem;
            if(ac){
                if(algo.fills_bus_mismatch()){
                    const Eigen::Ref<const CplxVect> mis = algo.get_bus_mismatch();
                    const real_type slack_absorbed = algo.get_slack_absorbed();
                    reactive_mismatch = mis.imag() * sn_mva;
                    active_mismatch = (mis.real().array() - slack_absorbed * slack_weights.array()) * sn_mva;
                } else {
                    const CplxVect ybus_v = Ybus * V_solver;
                    const CplxVect & sb = sbus_solver.size() > 0 ? CplxVect(sbus_solver) : CplxVect(ac_cache_.inj);
                    const CplxVect mismatch = V_solver.array() * ybus_v.array().conjugate() - sb.array();
                    active_mismatch = mismatch.real() * sn_mva;
                    reactive_mismatch = mismatch.imag() * sn_mva;
                }
                ctrl_q = algo.get_controller_q();
                ctrl_kind = algo.get_controller_kind();
                ctrl_elem = algo.get_controller_elem_id();
            } else {
                const Eigen::Index nb_bus = slack_weights.size();
                active_mismatch = RealVect::Zero(nb_bus);
                const real_type imbalance = _dc_imbalance_mw(sbus_solver, masked);
                for(Eigen::Index k = 0; k < nb_bus; ++k){
                    if(slack_weights(k) <= BaseConstants::my_zero_) continue;
                    active_mismatch(k) = slack_weights(k) * imbalance;
                }
            }

            RealVect p_mw, q_mvar;
            _grid_model.generator_results(ac, target_p, gen_off, no_share, storage_no_share, me_to_solver,
                                          active_mismatch, reactive_mismatch, ctrl_q, ctrl_kind, ctrl_elem,
                                          p_mw, q_mvar, gen_placed_me);
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                _gen_results_(static_cast<Eigen::Index>(i), 2 * gen_id) = p_mw(gen_id);
                _gen_results_(static_cast<Eigen::Index>(i), 2 * gen_id + 1) = q_mvar(gen_id);
            }
        }

        void _record_row_physical(size_t i, AlgorithmSelector & algo,
                                  const Eigen::Ref<const CplxVect> & V_solver,
                                  const Eigen::Ref<const RealVect> & slack_weights,
                                  const Eigen::Ref<const CplxVect> & sbus_solver){
            if(!_compute_physical_violations_) return;
            if(i >= _physical_violations_.size()) return;
            const std::vector<int> * masked = _row_masked_ids(i);
            // the outer loops' own triggers, on the row as it is: its own machines, targets,
            // branches and stranded buses (see OuterRowView)
            OuterRowView row;
            _fill_outer_row_view(i, masked, row);
            _take_moved_targets(algo, row);
            // the taps this row's outer loops left its transformers at (_record_row_trafo_blocks)
            if(i < _row_trafo_blocks_.size() && !_row_trafo_blocks_[i].empty()) row.moved_trafo_blocks = &_row_trafo_blocks_[i];
            _record_outer_detect(algo, V_solver, row, _physical_violations_[i]);
            if(_gen_pv_release_check_on_ && !_gen_pv_release_plan_.empty()){
                // the row's converged voltage, against the target this row gives the machine
                gen_pv_release_check::check_gen_pv_release_violations(
                    _gen_pv_release_plan_, V_solver, _physical_tol_vm_pu_, masked,
                    [this, i](int gen_id){ return this->_gen_target_vm_in_row(i, gen_id); },
                    [this, i](int gen_id){ return this->_gen_off_in_row(i, gen_id); },
                    _physical_violations_[i]);
            }
            if(!_gen_p_plan_.empty()){
                // the weights are THIS row's (a generator contingency re-derives them),
                // and so is the slack the solve absorbed -- both belong to the system the
                // algorithm handed in has just solved, which is why this runs here and not
                // after the row's state is put back
                gen_p_check::SlackShareInputs slack(algo.get_bus_mismatch(), slack_weights,
                                                    _grid_model.get_sn_mva(),
                                                    algo.ac_solver_used());
                if(slack.ac) slack.slack_absorbed = algo.get_slack_absorbed();
                else slack.dc_imbalance_mw = _dc_imbalance_mw(sbus_solver, masked);
                // the targets are this row's own -- the redistributed ones where the
                // slack pre-pass moved them -- and a unit that pre-pass saturated took
                // no share of what the solve had left to distribute
                gen_p_check::check_gen_p_violations(
                    _gen_p_plan_, slack, _physical_tol_mva_, masked,
                    [this, i](ViolationElementType el_type, int el_id){ return this->_row_target_p(i, el_type, el_id); },
                    [this, i](ViolationElementType el_type, int el_id){
                        return (el_type == ViolationElementType::GENERATOR) && this->_gen_off_in_row(i, el_id); },
                    [this, i](ViolationElementType el_type, int el_id){ return this->_row_takes_no_share(i, el_type, el_id); },
                    _physical_violations_[i]);
            }
        }

        void _store_row_status(size_t i, bool conv, bool invertible, const CplxVect & V_solver){
            _converged_mask_[i] = conv ? 1 : 0;
            if(!_compute_limit_violations_) return;
            _converged[i] = conv ? 1 : 0;
            if(!conv){
                _violations[i].push_back(LimitViolation{
                    ViolationElementType::GRID, -1, 0,
                    invertible ? LimitViolationType::DIVERGENCE : LimitViolationType::NOT_SIMULATED,
                    std::numeric_limits<real_type>::quiet_NaN(),
                    std::numeric_limits<real_type>::quiet_NaN()});
            } else {
                _record_row_violations_dispatch(i, V_solver, nullptr);
            }
        }

        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _refresh_defaults_vect_cache(){
            // read by the limit checks, and by the outer loops' view of each row
            _li_defaults_vect_cache_ = (_compute_limit_violations_ || _compute_physical_violations_ || _algo.supports_outer_loops())
                                       ? _branches_out_vect() : std::vector<std::vector<int> >();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _refresh_defaults_vect_cache(){}

        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _record_n_case_violations(const CplxVect & V_n){
            if(!_compute_limit_violations_) return;
            _converged_n_ = true;
            const std::vector<int> no_skip;
            // the buses the "n" solve masks (the union layout's extra busbars, which no
            // element of the base grid stands on) read their seed, not a solved voltage:
            // skipped, as in a row that does not use them
            const std::vector<int> * masked = _base_masked_.empty() ? nullptr : &_base_masked_;
            batch_sweep_detail::check_bus_voltage_violations(V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vmin_kv(), _grid_model.get_bus_vmax_kv(),
                                         _grid_model.get_bus_vn_kv(), _grid_model.get_substations(),
                                         _violation_threshold_, _violation_rel_tol_, masked, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, _violation_rel_tol_, no_skip, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, _violation_rel_tol_, no_skip, _violations_n_);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _record_n_case_violations(const CplxVect &){}

        // Once per compute(), before anything is solved: the one capability check (made
        // here, once, rather than discovered row by row) and this call's result buffers.
        // Sized even though the plans are not built yet, so that a batch whose base case
        // diverges answers get_physical_violations() with one empty entry per row rather
        // than with an empty list.
        void _prepare_physical_check(size_t nb_steps, bool ac_solver_used){
            _gen_pv_release_plan_.clear();
            _gen_p_plan_.clear();
            _physical_violations_.clear();
            _gen_pv_release_check_on_ = false;
            if(!_compute_physical_violations_) return;
            // the reactive half needs a per-bus mismatch; the hvdc half needs only the bus
            // angles, which every algorithm solves (see the flag's own doc)
            if(ac_solver_used){
                if(!_algo.fills_bus_mismatch()){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::compute: `compute_physical_violations` needs an "
                            "algorithm that publishes its per-bus mismatch -- that mismatch IS "
                            "the reactive power the machines pinning each bus had to produce. "
                            "Every built-in AC algorithm does; the active one ("
                         << _algo.get_name() << ") does not, so it is a plugin solver that has "
                            "not opted in (see BaseAlgo::fills_bus_mismatch). Pick a built-in AC "
                            "algorithm (change_algorithm), or turn "
                            "`compute_physical_violations` off.";
                    throw std::runtime_error(exc_.str());
                }
                // the reactive and voltage checks (the loops' and the release one) read the
                // mismatch and compare voltage magnitudes: AC only
                _gen_pv_release_check_on_ = true;
            }
            _physical_violations_.assign(nb_steps, std::vector<LimitViolation>());
        }

        // ... and the routing itself, once the labelling and the control plan this batch
        // solves against are settled (`active_layout()`, which the L1 / L2 preparation
        // above has built or kept). A plan means nothing against another labelling, so it
        // is never carried across a compute().
        void _build_physical_plans(){
            if(!_compute_physical_violations_) return;
            if(_gen_pv_release_check_on_){
                // which limit a flagged machine sits at is decided here, once: no batch axis
                // varies a reactive setpoint
                gen_pv_release_check::build_gen_pv_release_plan(
                    _grid_model, active_layout().id_me_to_solver, _gen_pv_release_plan_);
            }
            gen_p_check::build_gen_p_plan(_grid_model, active_layout().id_me_to_solver,
                                          _gen_p_plan_);
            _outer_detect_loops_ = detection_loops(_grid_model.get_outer_loops());
        }

        // the targets the solve's own loops moved (DistributedSlack), laid over the row's: a
        // unit's reactive limits follow its target P, and the solve held it within the ones
        // at that P (see LSGrid::get_physical_violations)
        static void _take_moved_targets(const AlgorithmSelector & algo, OuterRowView & row){
            if(!algo.supports_outer_loops()) return;
            std::vector<real_type> gen_target_p, storage_target_p;
            algo.get_outer_target_p(gen_target_p, storage_target_p);
            row.take_moved_targets(gen_target_p, storage_target_p);
        }

        // The grid's outer loops' own triggers (detection mode, see outer_detect), on the row
        // `row` describes, as LSGrid::get_physical_violations reads them on a grid (a bus out
        // of its reactive range, a remote controller at an unrealistic voltage, the hvdc droop
        // limits, the slack left on a single slack bus, the standby SVCs' thresholds -- the
        // voltage checks only with an algorithm that can feed the reactive ones).
        void _record_outer_detect(const AlgorithmSelector & algo,
                                  const Eigen::Ref<const CplxVect> & V_solver,
                                  const OuterRowView & row,
                                  std::vector<LimitViolation> & out){
            if(_outer_detect_loops_.empty()) return;
            const SolverBusLayout & layout = active_layout();
            const bool ac = algo.ac_solver_used();
            const CplxVect V(V_solver);
            const RealVect Va = algo.get_Va();
            const CplxVect mismatch = ac ? CplxVect(algo.get_bus_mismatch()) : CplxVect();
            const RealVect ctrl_q = ac ? algo.get_controller_q() : RealVect();
            OuterContext ctx;
            ctx.grid = &_grid_model;
            ctx.V = &V;
            ctx.Va = &Va;
            ctx.bus_mismatch = ac ? &mismatch : nullptr;
            ctx.controller_q = &ctrl_q;
            ctx.slack_bus = layout.slack_bus_id_solver.size() == 1 ? layout.slack_bus_id_solver[0].cast_int() : -1;
            ctx.slack_absorbed = ac ? algo.get_slack_absorbed() : 0.;
            ctx.tol_mw = _physical_tol_mva_;
            ctx.tol_vm_pu = _physical_tol_vm_pu_;
            ctx.vm_checks = ac && _gen_pv_release_check_on_;
            ctx.ac = ac;
            ctx.row = &row;
            outer_detect(_outer_detect_loops_, ctx, out);
        }

        // The base ("n") case's own report, read off the solve _finish_preprocessing just
        // ran on the member algorithm. Only called when that solve actually ran: a
        // compute() that REUSED the base case leaves the member algorithm holding the state
        // of the last row of the PREVIOUS call, and the report already stored (of the same
        // base case, under the same settings -- both setters drop it) is the right one to
        // keep.
        void _record_n_case_physical(){
            if(!_compute_physical_violations_) return;
            _physical_violations_n_.clear();
            // the base case is the grid's own: its targets, no contingency
            OuterRowView n_view = _base_row_view_;
            _take_moved_targets(_algo, n_view);
            std::vector<std::pair<int, std::array<cplx_type, 4> > > n_blocks;
            if(_algo.supports_outer_loops()){
                const TrafoContainer & trafos = _grid_model.get_trafos();
                std::vector<int> phase_tap, ratio_tap;
                _algo.get_outer_phase_tap(phase_tap);
                _algo.get_outer_ratio_tap(ratio_tap);
                n_blocks = trafos.moved_tap_blocks(trafos.moved_taps(false, ratio_tap), trafos.moved_taps(true, phase_tap));
                if(!n_blocks.empty()) n_view.moved_trafo_blocks = &n_blocks;
            }
            _record_outer_detect(_algo, _algo.get_V(), n_view, _physical_violations_n_);
            if(_gen_pv_release_check_on_ && !_gen_pv_release_plan_.empty()){
                gen_pv_release_check::check_gen_pv_release_violations(
                    _gen_pv_release_plan_, _algo.get_V(), _physical_tol_vm_pu_, nullptr,
                    [this](int gen_id){ return this->_grid_model.get_generators().get_target_vm_pu(gen_id); },
                    [](int){ return false; },
                    _physical_violations_n_);
            }
            if(!_gen_p_plan_.empty()){
                // the base case is the grid's own: its targets, its weights, no
                // contingency -- and, in DC, the injection _solve_n_case handed the solver
                const bool ac = _algo.ac_solver_used();
                gen_p_check::SlackShareInputs slack(_algo.get_bus_mismatch(),
                                                    active_layout().slack_weights,
                                                    _grid_model.get_sn_mva(), ac);
                if(ac) slack.slack_absorbed = _algo.get_slack_absorbed();
                else slack.dc_imbalance_mw = -dc_cache_.inj.sum() * _grid_model.get_sn_mva();
                gen_p_check::check_gen_p_violations(
                    _gen_p_plan_, slack, _physical_tol_mva_, nullptr,
                    [this](ViolationElementType el_type, int el_id){ return this->_grid_target_p_of(el_type, el_id); },
                    [](ViolationElementType, int){ return false; },  // the base case disconnects no generator
                    [](ViolationElementType, int){ return false; },  // ... and redistributes nothing
                    _physical_violations_n_);
            }
        }

        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_masks(){
            if(!_mask_mode()) return;
            // Masking is a value-level edit at constant sparsity, and one that can
            // move a pivot: the voltage row of a stranded group goes from "Vm(reg)
            // = Vset" to "Q_first = 0", so the entry KLU pivoted at in the base
            // factorization is a zero in that row's matrix, and klu_refactor halts
            // on it (SparseLU re-pivots and never noticed). Let the linear solver
            // fall back to a numeric factorize on such a row -- and on the next
            // ordinary row, whose pivots the stranded matrix moved back. A bool
            // store; the fallback itself only ever runs on a failure.
            _algo.set_refactor_fallback(true);
            if(!_voltage_control_may_mask_wired_){
                // first compute() with masking on for this object (or the first ever):
                // _algo's VoltageControl extension needs its stranded-controller
                // Jacobian slot, which must exist in the sparsity BEFORE the next
                // build_J_sparsity() -- force one. A no-op the next time (this stays
                // true until the object is destroyed; turning handle_disconnected_grid
                // back off does not need a matching rebuild, it just leaves the slot
                // reserved and unused).
                _algo.set_may_mask_voltage_control(true);
                _algo_controler.tell_pv_changed();
                _voltage_control_may_mask_wired_ = true;
            }
            if(!_algo.supports_bus_masking()){
                // name what asked for the masking: topological actions turn it on without
                // handle_disconnected_grid (see _mask_mode)
                std::ostringstream exc_;
                if(_handle_disconnected_grid){
                    exc_ << algo_name() << ": the `handle_disconnected_grid` mode requires a "
                            "Newton-Raphson algorithm (AC) or the DC solver (the active algorithm "
                            "does not support bus masking). Use `change_algorithm` to select an NR "
                            "solver (e.g. NR_KLU / NR_SLU) or the DC solver.";
                } else {
                    exc_ << algo_name() << "::set_topo_actions: topological actions require a "
                            "Newton-Raphson algorithm (the active one, " << _algo.get_name()
                         << ", cannot mask the buses a row leaves unused). Use `change_algorithm` "
                            "to select an NR solver (e.g. NR_KLU / NR_SLU).";
                }
                throw std::runtime_error(exc_.str());
            }
            _select_ref_slack_and_masks(true);
            _skip_rows_stranding_regulated_bus();
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_masks(){}

        // ================= topological actions (ScenarioSweep only) ==============
        // Turns topo_actions_ into what the row loop consumes, once per compute()
        // (L2, before _prepare_ybus_varying): the branches each row disconnects
        // (ybus_policy_.topo_branches_off -- read wherever the masks are, so the Ybus
        // edit, the connectivity, the current checks and the flow cleaning all see
        // them) and the ones it places (topo_branches_moved: an end moved, a branch
        // reconnected), the generators it takes off their bus (sbus_policy_.topo_gen_off
        // -- read wherever gen_off is, so the injection, the slack re-weighting and the
        // PV -> PQ pinning all see them) and the ones it places (topo_gens_on), the loads
        // / storage units it takes off or places (topo_loads_off / _on, topo_storages_off
        // / _on, injections only), and the base buses it leaves empty.
        //
        // What is refused here, and why, is the frontier of this version: a row naming
        // an element in both a mask and its action is ambiguous; a slack participant
        // moved or reactivated would change the slack set with the row, a regulating
        // storage unit placed the PV set, with no slot reserved for either;
        // keep_jacobian's gen_v gradient maps a generator to its base bus; the DC path
        // is not wired. Refused elsewhere: a half-open branch put on (_maybe_topo_union),
        // a generator reactivated on a slack bus or held remotely
        // (_maybe_prepare_gen_contingency).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _maybe_resolve_topology(size_t nb_steps){
            _topology_active_ = false;
            _topo_resolved_ = false;
            _base_masked_.clear();
            _row_emptied_buses_.clear();
            _row_placed_buses_.clear();
            _row_slack_storages_off_.clear();
            _row_branch_overrides_.clear();
            _reset_topo_policy_state();
            if(topo_actions_.empty()){
                _topo_resolved_ = true;
                return;
            }
            if(topo_actions_.size() != nb_steps){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_topo_actions: " << topo_actions_.size()
                     << " actions were registered for a batch of " << nb_steps << " rows.";
                throw std::runtime_error(exc_.str());
            }
            const auto & generators = _grid_model.get_generators();
            const auto & storages = _grid_model.get_storages();
            const int nb_gen = static_cast<int>(generators.nb());
            const int nb_line = static_cast<int>(n_line_);
            const auto & id_me_to_solver = active_layout().id_me_to_solver;
            // every bus a row uses is in the layout (see _maybe_topo_union): a -1 here
            // is a bug, not an input error
            auto solver_of = [&](int bus_me, size_t row) -> int {
                const int res = bus_me == BaseConstants::_deactivated_bus_id ? BaseConstants::_deactivated_bus_id
                                                                             : id_me_to_solver[bus_me].cast_int();
                if(res == BaseConstants::_deactivated_bus_id){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::set_topo_actions: internal error, the action of row " << row
                         << " uses bus " << bus_me << ", which the solver labelling does not hold.";
                    throw std::runtime_error(exc_.str());
                }
                return res;
            };
            if(_row_topo_.size() != nb_steps){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_topo_actions: internal error, " << _row_topo_.size()
                     << " actions were resolved with the grid for a batch of " << nb_steps << " rows.";
                throw std::runtime_error(exc_.str());
            }
            _row_emptied_buses_.assign(nb_steps, std::vector<int>());
            _row_placed_buses_.assign(nb_steps, std::vector<int>());
            _row_slack_storages_off_.assign(nb_steps, std::vector<int>());
            const std::vector<std::size_t> & base_count = _grid_model.get_substations().get_nb_elements_per_bus();
            ybus_policy_.topo_branches_off.assign(nb_steps, std::vector<int>());
            ybus_policy_.topo_branches_moved.assign(nb_steps, std::vector<typename YbusPolicy::Contingency::BranchPlacement>());
            sbus_policy_.topo_loads_off.assign(nb_steps, std::vector<int>());
            sbus_policy_.topo_loads_on.assign(nb_steps, std::vector<std::pair<int, int> >());
            sbus_policy_.topo_storages_off.assign(nb_steps, std::vector<int>());
            sbus_policy_.topo_storages_on.assign(nb_steps, std::vector<std::pair<int, int> >());
            sbus_policy_.topo_gens_on.assign(nb_steps, std::vector<typename SbusPolicy::Vary::GenOn>());
            _row_branch_overrides_.assign(nb_steps, std::vector<BranchBusOverride>());
            bool any_gen_off = false;
            auto gen_off_in_row = [&](size_t row, int gen_id){
                if(!any_gen_off){
                    sbus_policy_.topo_gen_off = BoolMat::Constant(static_cast<Eigen::Index>(nb_steps), nb_gen, false);
                    any_gen_off = true;
                }
                sbus_policy_.topo_gen_off(static_cast<Eigen::Index>(row), gen_id) = true;
            };
            for(size_t row = 0; row < nb_steps; ++row){
                // resolved by _maybe_topo_union, with the grid
                const RowTopoPlan & plan = _row_topo_[row];
                if(plan.empty()) continue;
                _topology_active_ = true;
                _row_emptied_buses_[row] = _emptied_buses(plan, base_count, solver_of, row);
                _row_placed_buses_[row] = _placed_buses(plan, solver_of, row);
                if(!_algo.ac_solver_used()){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::set_topo_actions: row " << row << " carries a topological "
                            "action, which the DC algorithm does not support yet. Use an AC "
                            "Newton-Raphson algorithm (change_algorithm), or an action that changes nothing.";
                    throw std::runtime_error(exc_.str());
                }
                // ---- branches: off (as a mask entry), or placed by the row ----------
                for(const TopoBranchPlacement & br : plan.branches){
                    const bool is_trafo = br.branch_id >= nb_line;
                    const int local_id = is_trafo ? br.branch_id - nb_line : br.branch_id;
                    const char * kind = is_trafo ? "transformer" : "powerline";
                    const bool masked = is_trafo
                        ? (ybus_policy_.trafo_mask.rows() > 0 && ybus_policy_.trafo_mask(static_cast<Eigen::Index>(row), local_id))
                        : (ybus_policy_.line_mask.rows() > 0 && ybus_policy_.line_mask(static_cast<Eigen::Index>(row), local_id));
                    if(masked) _throw_topo_overlap(row, kind, local_id, is_trafo ? "set_contingency_trafos" : "set_contingency_lines");
                    BranchBusOverride ov;
                    ov.branch_id = br.branch_id;
                    ov.from_me = br.row_bus1_me;
                    ov.to_me = br.row_bus2_me;
                    ov.connected = br.row_on;
                    _row_branch_overrides_[row].push_back(ov);
                    if(!br.row_on){
                        ybus_policy_.topo_branches_off[row].push_back(br.branch_id);
                    }else{
                        typename YbusPolicy::Contingency::BranchPlacement placement;
                        placement.branch_id = br.branch_id;
                        placement.bus1_solver = solver_of(br.row_bus1_me, row);
                        placement.bus2_solver = solver_of(br.row_bus2_me, row);
                        placement.base_on = br.base_on;
                        ybus_policy_.topo_branches_moved[row].push_back(placement);
                    }
                }
                std::sort(ybus_policy_.topo_branches_off[row].begin(), ybus_policy_.topo_branches_off[row].end());
                std::sort(_row_branch_overrides_[row].begin(), _row_branch_overrides_[row].end(),
                          [](const BranchBusOverride & a, const BranchBusOverride & b){ return a.branch_id < b.branch_id; });
                // ---- generators: off (the generator-contingency path), or placed ------
                for(const TopoElPlacement & g : plan.gens){
                    const bool was_on = g.base_bus_me != BaseConstants::_deactivated_bus_id;
                    if(sbus_policy_.gen_off.rows() > 0 && sbus_policy_.gen_off(static_cast<Eigen::Index>(row), g.el_id)){
                        if(was_on) _throw_topo_overlap(row, "generator", g.el_id, "set_contingency_gens");
                        // off in the base grid: the mask alone changes nothing, but the
                        // action puts it back on
                        std::ostringstream exc_;
                        exc_ << algo_name() << "::set_topo_actions: row " << row << " reactivates generator "
                             << g.el_id << " through its topological action while set_contingency_gens "
                                "takes it out. Name it in one of the two only.";
                        throw std::runtime_error(exc_.str());
                    }
                    if(g.row_bus_me == BaseConstants::_deactivated_bus_id){
                        gen_off_in_row(row, g.el_id);
                        continue;
                    }
                    // placed on a bus (back on its last one, or moved): what a bus's
                    // label costs per row is handled by _maybe_prepare_gen_contingency;
                    // what cannot move per row is refused here
                    if(abs(generators.get_gen_slack_weight(g.el_id)) > BaseConstants::_tol_equal_float){
                        std::ostringstream exc_;
                        exc_ << algo_name() << "::set_topo_actions: the action of row " << row
                             << (was_on ? " moves" : " reactivates") << " generator " << g.el_id
                             << ", which takes part in the distributed slack. The slack set would change "
                                "with the row: not supported yet (see the TODO in the changelog).";
                        throw std::runtime_error(exc_.str());
                    }
                    if(_keep_jacobian_){
                        std::ostringstream exc_;
                        exc_ << algo_name() << "::set_topo_actions: the action of row " << row
                             << (was_on ? " moves" : " reactivates") << " generator " << g.el_id
                             << ", which `keep_jacobian` does not support yet (the gen_v gradient maps "
                                "every generator to the bus it holds in the base grid, see "
                                "get_gen_v_target_bus). Turn it off.";
                        throw std::runtime_error(exc_.str());
                    }
                    if(was_on) gen_off_in_row(row, g.el_id);   // gone from its base bus
                    typename SbusPolicy::Vary::GenOn gen_on;
                    gen_on.gen_id = g.el_id;
                    gen_on.bus_me = g.row_bus_me;
                    gen_on.bus_solver = solver_of(g.row_bus_me, row);
                    sbus_policy_.topo_gens_on[row].push_back(gen_on);
                }
                // ---- loads and storage units: injections only ------------------------
                for(const TopoElPlacement & l : plan.loads){
                    // (a load also in set_contingency_loads: _check_topo_load_overlap)
                    if(l.base_bus_me != BaseConstants::_deactivated_bus_id) sbus_policy_.topo_loads_off[row].push_back(l.el_id);
                    if(l.row_bus_me != BaseConstants::_deactivated_bus_id){
                        sbus_policy_.topo_loads_on[row].push_back(std::make_pair(l.el_id, solver_of(l.row_bus_me, row)));
                    }
                }
                for(const TopoElPlacement & st : plan.storages){
                    if(st.row_bus_me != BaseConstants::_deactivated_bus_id && storages.get_voltage_regulator_on(st.el_id)){
                        std::ostringstream exc_;
                        exc_ << algo_name() << "::set_topo_actions: the action of row " << row
                             << " places storage unit " << st.el_id << ", which regulates voltage: not supported yet.";
                        throw std::runtime_error(exc_.str());
                    }
                    // a share of the distributed slack: taken out, the row's weights are
                    // re-derived without it (_row_slack_weights); placed, it would take
                    // its share on another bus -- the slack set would change with the row,
                    // refused as for a generator
                    if(abs(storages.get_slack_weight(st.el_id)) > BaseConstants::_tol_equal_float){
                        if(st.row_bus_me != BaseConstants::_deactivated_bus_id){
                            const bool was_on = st.base_bus_me != BaseConstants::_deactivated_bus_id;
                            std::ostringstream exc_;
                            exc_ << algo_name() << "::set_topo_actions: the action of row " << row
                                 << (was_on ? " moves" : " reactivates") << " storage unit " << st.el_id
                                 << ", which takes part in the distributed slack. The slack set would change "
                                    "with the row: not supported yet (see the TODO in the changelog).";
                            throw std::runtime_error(exc_.str());
                        }
                        _row_slack_storages_off_[row].push_back(st.el_id);
                    }
                    if(st.base_bus_me != BaseConstants::_deactivated_bus_id) sbus_policy_.topo_storages_off[row].push_back(st.el_id);
                    if(st.row_bus_me != BaseConstants::_deactivated_bus_id){
                        sbus_policy_.topo_storages_on[row].push_back(std::make_pair(st.el_id, solver_of(st.row_bus_me, row)));
                    }
                }
            }
            if(!_topology_active_){
                // every action is a no-op: nothing for the row loop to read
                _reset_topo_policy_state();
                _row_branch_overrides_.clear();
            }
            // the base mask: the buses of the union layout that are empty in the base
            // grid, masked in the "n" case and in every row that does not use them
            for(int bus_me : _extra_buses_me_){
                const int b = id_me_to_solver[bus_me].cast_int();
                if(b != BaseConstants::_deactivated_bus_id) _base_masked_.push_back(b);
            }
            std::sort(_base_masked_.begin(), _base_masked_.end());
            // on the member algorithm from here on: compute() reset it just before this,
            // and the "n" solve runs with it (the workers set it themselves)
            _algo.set_masked_buses(_base_masked_);
            _topo_resolved_ = true;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && S::supports_vary), int>::type = 0>
        void _maybe_resolve_topology(size_t){}

        // ================= the union layout (L1) ==================================
        // The solver labelling is built once for the batch, so it has to hold every bus
        // any row uses -- a busbar no element stands on in the base grid included --
        // and the Ybus pattern every entry any row writes: a branch end moved to
        // another busbar, a branch reconnected. Both are read off the actions here,
        // before the grid cache is built; the buses go into the labelling
        // (LSGrid::init_converter_bus_id's `keep_also`), the entries are reserved as
        // stored zeros afterwards (_maybe_reserve_union_pattern), before the "n" solve
        // builds the Jacobian's sparsity from that pattern.
        void _maybe_topo_union(){
            _topo_used_buses_me_.clear();
            _union_edges_me_.clear();
            _extra_buses_me_.clear();
            _row_topo_.clear();
            if(topo_actions_.empty()) return;
            std::set<int> used;
            std::set<std::pair<int, int> > edges;
            // kept: _maybe_resolve_topology reads the same plans (both levels are rebuilt
            // with the grid or the actions, see set_topo_actions)
            _row_topo_.resize(topo_actions_.size());
            for(size_t row = 0; row < topo_actions_.size(); ++row){
                RowTopoPlan & plan = _row_topo_[row];
                resolve_row_topo(topo_actions_[row], _grid_model, plan);
                for(const TopoElPlacement & el : plan.loads) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoElPlacement & el : plan.gens) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoElPlacement & el : plan.storages) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoBranchPlacement & br : plan.branches){
                    if(!br.row_on) continue;
                    _check_branch_placement(row, br);
                    used.insert(br.row_bus1_me);
                    used.insert(br.row_bus2_me);
                    edges.insert(std::make_pair(std::min(br.row_bus1_me, br.row_bus2_me),
                                                std::max(br.row_bus1_me, br.row_bus2_me)));
                }
            }
            _topo_used_buses_me_.assign(used.begin(), used.end());
            _union_edges_me_.assign(edges.begin(), edges.end());
        }

        // What a row cannot put a branch on yet. A branch with one end open in the base
        // grid ("half-open"): the placement would stamp the closed-branch block, so the
        // row would close that end, and whether it should is not decided -- the one-off
        // TopoAction::apply_to_gridmodel closes it on the bus it was last on, even for a
        // set_line_status the plan would otherwise read as a no-op. Taking such a branch
        // out is well defined (its effective block goes) and stays allowed. And a branch
        // with no bus recorded on an end it goes back to: there is no bus to stamp.
        void _check_branch_placement(size_t row, const TopoBranchPlacement & br) const {
            const bool is_trafo = br.branch_id >= static_cast<int>(n_line_);
            const int local_id = is_trafo ? br.branch_id - static_cast<int>(n_line_) : br.branch_id;
            const char * kind = is_trafo ? "transformer" : "powerline";
            if(br.base_half_open()){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_topo_actions: the action of row " << row << " puts " << kind << " "
                     << local_id << " on, one end of which is open in the base grid: closing that end is not "
                        "supported yet. Disconnecting it is.";
                throw std::runtime_error(exc_.str());
            }
            if(br.row_bus1_me < 0 || br.row_bus2_me < 0){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_topo_actions: the action of row " << row << " reconnects " << kind << " "
                     << local_id << ", which has no bus recorded on side " << (br.row_bus1_me < 0 ? 1 : 2)
                     << ": give that end a busbar in the action.";
                throw std::runtime_error(exc_.str());
            }
        }

        // The base buses (solver ids, sorted) a row's action leaves with no element at all:
        // the busbar a merge empties. Counted off the base grid's per-bus element counts,
        // each element the plan takes off a bus and puts on another moved accordingly (a
        // branch end counts where it is closed). Only the action's own moves count, so a
        // row without one keeps the answer the masks always gave.
        template<class SolverOf>
        std::vector<int> _emptied_buses(const RowTopoPlan & plan, const std::vector<std::size_t> & base_count,
                                        const SolverOf & solver_of, size_t row) const {
            std::map<int, long> change;
            auto move = [&change](int bus_me, long delta){ if(bus_me >= 0) change[bus_me] += delta; };
            for(const TopoBranchPlacement & br : plan.branches){
                if(br.base_on && br.base_on1) move(br.base_bus1_me, -1);
                if(br.base_on && br.base_on2) move(br.base_bus2_me, -1);
                if(br.row_on){
                    move(br.row_bus1_me, 1);
                    move(br.row_bus2_me, 1);
                }
            }
            for(const std::vector<TopoElPlacement> * els : {&plan.loads, &plan.gens, &plan.storages}){
                for(const TopoElPlacement & el : *els){
                    move(el.base_bus_me, -1);
                    move(el.row_bus_me, 1);
                }
            }
            std::vector<int> res;
            for(const auto & bus_change : change){
                const std::size_t bus = static_cast<std::size_t>(bus_change.first);
                if(bus >= base_count.size() || base_count[bus] == 0) continue;   // not a bus of the base grid
                if(static_cast<long>(base_count[bus]) + bus_change.second == 0) res.push_back(solver_of(bus_change.first, row));
            }
            std::sort(res.begin(), res.end());
            return res;
        }

        // The buses (solver ids, sorted) a row's action puts an element or a branch end on:
        // an extra bus of the union layout is masked by design only in the rows that do
        // not use it (see _only_masked_by_design).
        template<class SolverOf>
        std::vector<int> _placed_buses(const RowTopoPlan & plan, const SolverOf & solver_of, size_t row) const {
            std::vector<int> res;
            for(const TopoBranchPlacement & br : plan.branches){
                if(!br.row_on) continue;
                res.push_back(solver_of(br.row_bus1_me, row));
                res.push_back(solver_of(br.row_bus2_me, row));
            }
            for(const std::vector<TopoElPlacement> * els : {&plan.loads, &plan.gens, &plan.storages}){
                for(const TopoElPlacement & el : *els){
                    if(el.row_bus_me != BaseConstants::_deactivated_bus_id) res.push_back(solver_of(el.row_bus_me, row));
                }
            }
            std::sort(res.begin(), res.end());
            res.erase(std::unique(res.begin(), res.end()), res.end());
            return res;
        }

        // once the cache is built (the element counts are then settled): which of the
        // buses the rows use are empty in the base grid -- the ones to mask when unused
        void _maybe_finalize_extra_buses(){
            _extra_buses_me_.clear();
            const SubstationContainer & subs = _grid_model.get_substations();
            for(int bus_me : _topo_used_buses_me_){
                if(!subs.is_bus_connected(GlobalBusId(bus_me))) _extra_buses_me_.push_back(bus_me);
            }
        }

        // reserve, as stored zeros, the entries of the admittance matrix a row may
        // write that the base pattern does not hold: the diagonal of every extra bus,
        // and the four entries of every branch placement. Done at L1, right after the
        // cache is built and before anything reads its pattern (the algorithm is reset
        // at L2, the graph and the coefficient lists are built there).
        void _maybe_reserve_union_pattern(bool ac_solver_used){
            if(_extra_buses_me_.empty() && _union_edges_me_.empty()) return;
            if(!ac_solver_used) return;   // a DC row with an action is refused at L2
            const auto & id_me_to_solver = active_layout().id_me_to_solver;
            Eigen::SparseMatrix<cplx_type> & mat = ac_cache_.mat;
            auto stored = [&mat](int r, int c){
                for(Eigen::SparseMatrix<cplx_type>::InnerIterator it(mat, c); it; ++it){
                    if(it.row() == r) return true;
                }
                return false;
            };
            auto reserve = [&](int r, int c){
                if(r < 0 || c < 0) return;
                if(!stored(r, c)) mat.coeffRef(r, c) = cplx_type(0., 0.);
            };
            for(int bus_me : _extra_buses_me_){
                const int b = id_me_to_solver[bus_me].cast_int();
                reserve(b, b);
            }
            for(const auto & edge : _union_edges_me_){
                const int b1 = id_me_to_solver[edge.first].cast_int();
                const int b2 = id_me_to_solver[edge.second].cast_int();
                reserve(b1, b1);
                reserve(b2, b2);
                reserve(b1, b2);
                reserve(b2, b1);
            }
            mat.makeCompressed();
        }

        // an extra bus starts from a finite, nonzero voltage even where it is masked:
        // the Newton-Raphson divides by |V| at every bus, and the caller's starting
        // voltage is often exactly 0 on a busbar the grid does not use
        void _maybe_seed_extra_buses(Eigen::Ref<CplxVect> V) const {
            if(_extra_buses_me_.empty()) return;
            const auto & id_me_to_solver = active_layout().id_me_to_solver;
            const SubstationContainer & subs = _grid_model.get_substations();
            for(int bus_me : _extra_buses_me_){
                const int b = id_me_to_solver[bus_me].cast_int();
                if(b < 0 || b >= V.size()) continue;
                if(std::abs(V(b)) > 0.) continue;
                cplx_type seed(_grid_model.get_init_vm_pu(), 0.);
                const int first_busbar = subs.local_to_gridmodel(subs.sub_id_of_bus(bus_me), LocalBusId(1)).cast_int();
                const int b1 = id_me_to_solver[first_busbar].cast_int();
                if(b1 >= 0 && b1 < V.size() && std::abs(V(b1)) > 0.) seed = V(b1);
                V(b) = seed;
            }
        }

        // the masked row loop runs where handle_disconnected_grid says so, and where a
        // row's topological action does something (its extra buses are masked when
        // unused). Actions that all change nothing leave plain rows, on any algorithm.
        // Read from _maybe_resolve_topology on, which is what settles _topology_active_
        // (every caller runs after it: the connectivity, the masks, the row loops).
        bool _mask_mode() const { return _handle_disconnected_grid || _topology_active_; }

        // A row naming a load both in set_contingency_loads and in its topological action
        // is refused. Not with the other overlaps in _maybe_resolve_topology: the load
        // mask is an injection, set without dropping the batch inputs, so a mask set after
        // a compute() never goes through that resolution again. Checked at every compute().
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _check_topo_load_overlap() const {
            if(topo_actions_.empty()) return;
            for(size_t row = 0; row < _row_topo_.size(); ++row){
                for(const TopoElPlacement & l : _row_topo_[row].loads){
                    if(sbus_policy_.load_off_in(static_cast<Eigen::Index>(row), l.el_id))
                        _throw_topo_overlap(row, "load", l.el_id, "set_contingency_loads");
                }
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && S::supports_vary), int>::type = 0>
        void _check_topo_load_overlap() const {}

        [[noreturn]] void _throw_topo_overlap(size_t row, const char * el_kind, int el_id, const char * setter) const {
            std::ostringstream exc_;
            exc_ << algo_name() << "::set_topo_actions: row " << row << " disconnects " << el_kind << " "
                 << el_id << " both through " << setter << " and through its topological action. "
                    "Name it in one of the two only.";
            throw std::runtime_error(exc_.str());
        }

        // what _maybe_resolve_topology wrote into the two policies (2-way: the
        // fields only exist on the (Contingency, Vary) instantiation)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _reset_topo_policy_state(){
            ybus_policy_.topo_branches_off.clear();
            // the placements too: init_li_coeffs_from_masks replays whatever is left here,
            // on a batch that registered no action as well
            ybus_policy_.topo_branches_moved.clear();
            sbus_policy_.clear_topo();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && S::supports_vary), int>::type = 0>
        void _reset_topo_policy_state(){}

        // "handle disconnected grid" mode, AC: a row that strands the regulated bus of a
        // voltage-control group while some of its controllers stay in the main component
        // cannot be solved -- the group's voltage row asks a masked (frozen) magnitude to
        // be v_set, which none of the live controllers' reactive unknowns can reach -- so
        // it is skipped (NOT_SIMULATED) rather than left to diverge. The one contingency
        // at a time path refuses it as well ("regulates a disconnected bus"). The other
        // cases need nothing: a group whose controllers are all stranded is pinned to
        // Q = 0 by the VoltageControl extension (see VoltageControl::set_masked_buses),
        // and in a group with only some of them stranded the live ones hold the bus and
        // share among themselves, the stranded ones' reactive power landing on masked
        // rows. Called after _select_ref_slack_and_masks, which (re)fills _skip_mask.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _skip_rows_stranding_regulated_bus(){
            if(!_algo.ac_solver_used()) return;  // no voltage control in DC
            const VoltageControlSolverData & ctrl = active_layout().voltage_control.controllers();
            const int ng = ctrl.n_groups();
            if(ng == 0) return;
            // every _li_masked entry is sorted (see _prepare_connectivity)
            auto is_masked = [](const std::vector<int> & masked, int bus){
                return std::binary_search(masked.begin(), masked.end(), bus);
            };
            const size_t nb_rows = std::min(_skip_mask.size(), _li_masked.size());
            for(size_t row = 0; row < nb_rows; ++row){
                const std::vector<int> & masked = _li_masked[row];
                if(_skip_mask[row] || masked.empty()) continue;
                for(int g = 0; g < ng; ++g){
                    if(!is_masked(masked, ctrl.reg_bus(g))) continue;
                    const int first = ctrl.grp_start(g);
                    // a held controller (LSGrid::set_hold_frozen_regulators) holds no bus
                    bool some_live = false;
                    for(int off = 0; off < ctrl.grp_count(g) && !some_live; ++off){
                        if(ctrl.is_held(first + off)) continue;
                        some_live = !is_masked(masked, ctrl.bus(first + off));
                    }
                    if(some_live){ _skip_mask[row] = 1; break; }
                }
            }
        }

        // ================= generator contingencies ===============================
        // ScenarioSweep (set_contingency_gens, set_topo_actions) and ContingencyAnalysis
        // (a contingency tripping a generator), read alike through _gen_off_in_row.
        // Turns the generators each row takes out into the two things the per-row loop needs, once
        // per compute(): which buses can lose their voltage pinning at all
        // (_switchable_buses_, handed to the algorithm so the Jacobian reserves a Vm
        // unknown + Q equation for each), and which of them actually do, row by row
        // (_row_pv_to_pq_). Also collects the slack-participating generators each row
        // takes out, for _row_slack_weights.
        //
        // Runs BEFORE _finish_preprocessing, so the "n" warm-up solve -- the one that
        // builds the sparsity, under its own tell_all_changed() -- already sees the
        // enlarged layout. That is what makes the whole sweep share ONE symbolic
        // analysis, the base case included.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_gen_contingency(size_t nb_steps){
            _switchable_buses_.clear();
            _row_pv_to_pq_.clear();
            _row_pv_pinned_.clear();
            _row_vm_seed_.clear();
            _pv_pinning_active_ = false;
            _row_slack_gens_off_.clear();
            _gen_contingency_active_ = false;
            // a generator is taken out by the set_contingency_gens mask or by the row's
            // topological action (set_bus -1, see _maybe_resolve_topology): the two
            // axes are read together through _gen_off_in_row (on a ContingencyAnalysis, the
            // generators its contingencies trip). The row's action
            // may also put a base-off generator back (topo_gens_on): its bus turns PV.
            // And it may take out a storage unit that regulates its bus: one controller
            // fewer there too.
            const bool has_gens_on = _has_gens_on();
            if(!_any_gen_off() && !has_gens_on && !_has_regulating_storage_off()){
                // No mask this time. _algo is a member and outlives the compute() that
                // reserved slots for it, so it has to be told the set is empty again --
                // otherwise a bus a PREVIOUS compute() made switchable would keep its Vm
                // unknown + Q equation with nothing pinning them, and silently solve as
                // PQ. (_finish_preprocessing's tell_all_changed() rebuilds the sparsity.)
                _push_switchable_to_algo();
                return;
            }
            _gen_contingency_active_ = true;

            const auto & generators = _grid_model.get_generators();
            const int nb_gen = static_cast<int>(generators.nb());
            const auto & id_me_to_solver = active_layout().id_me_to_solver;
            const Eigen::Index nb_rows = static_cast<Eigen::Index>(nb_steps);

            // ---- the distributed slack, AC and DC alike ---------------------------
            // Which participating machines each row takes out. Needed on both families
            // (DC has a distributed slack too), and the only part of this that applies
            // in DC at all -- see the early return below.
            // the generators each row has out, listed once (see _row_gens_off): what the
            // loops below read, rather than asking _gen_off_in_row for every (row, generator)
            std::vector<std::vector<int> > gens_off(nb_steps);
            for(size_t step = 0; step < nb_steps; ++step) _row_gens_off(step, nb_gen, gens_off[step]);
            _row_slack_gens_off_.assign(nb_steps, std::vector<int>());
            for(size_t step = 0; step < nb_steps; ++step){
                for(int gen_id : gens_off[step]){
                    if(abs(generators.get_gen_slack_weight(gen_id)) < BaseConstants::_tol_equal_float) continue;
                    _row_slack_gens_off_[step].push_back(gen_id);
                }
            }

            // DC knows no PV / PQ distinction: |V| is not a variable there, so nothing
            // is relabelled and no Jacobian slot has to be reserved. The injection
            // correction (SbusPolicy::Vary::fill_row) and the slack re-weighting above
            // are the whole of the DC story.
            if(!_algo.ac_solver_used()){
                _push_switchable_to_algo();  // empty here: nothing is switchable in DC
                return;
            }

            if(!_algo.supports_pv_pinning()){
                std::ostringstream exc_;
                exc_ << algo_name() << ": generator contingencies (" << _gen_trip_api() << ") require a "
                        "Newton-Raphson algorithm (the active algorithm cannot relabel a bus PV -> PQ "
                        "without rebuilding the Jacobian). Use `change_algorithm` to select an NR "
                        "solver (e.g. NR_KLU / NR_SLU), or the DC solver.";
                throw std::runtime_error(exc_.str());
            }

            // ---- 1. reject what this version cannot model --------------------------
            // Only a generator pinning its OWN bus is in scope. A remote controller --
            // or a bus a group (remote generator, SVC, HVDC station) controls -- lives
            // in the VoltageControl extension, which owns Jacobian columns and rows of
            // its own that nothing here reserves or masks.
            const std::set<int> & group_buses = _grid_model.get_ac_voltage_control_plan().group_controlled_buses();
            std::vector<char> ever_off(static_cast<std::size_t>(nb_gen), 0);
            for(const std::vector<int> & off : gens_off) for(int gen_id : off) ever_off[static_cast<std::size_t>(gen_id)] = 1;
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                if(!ever_off[static_cast<std::size_t>(gen_id)]) continue;
                const int bus_id_me = generators.get_bus_id()(gen_id).cast_int();
                const bool remote = generators.is_remote_voltage_controller(gen_id);
                const bool in_group = bus_id_me != BaseConstants::_deactivated_bus_id &&
                                      group_buses.find(bus_id_me) != group_buses.end();
                if(remote || in_group){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::" << _gen_trip_api() << ": generator " << gen_id
                         << (remote ? " regulates the voltage of a remote bus"
                                    : " stands on a bus whose voltage a control group holds (a remote "
                                      "generator, an SVC or an HVDC converter station)")
                         << ". Only generators regulating their own bus can be disconnected for now: "
                            "remote voltage control is not supported by this feature yet.";
                    throw std::runtime_error(exc_.str());
                }
            }

            // ---- 2. which bus does each local voltage controller pin? --------------
            // -1 for a generator that pins nothing (disconnected, not regulating,
            // regulating remotely, or "pseudo off" when that rule is active).
            std::vector<int> pinned_bus_of_gen(nb_gen, -1);
            std::map<int, std::vector<int> > gens_of_bus;  // solver bus -> its local controllers
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                if(!generators.is_local_voltage_controller(gen_id)) continue;
                const int bus_id_me = generators.get_bus_id()(gen_id).cast_int();
                if(bus_id_me == BaseConstants::_deactivated_bus_id) continue;
                const int bus_solver = id_me_to_solver[bus_id_me].cast_int();
                if(bus_solver == BaseConstants::_deactivated_bus_id) continue;
                pinned_bus_of_gen[gen_id] = bus_solver;
                gens_of_bus[bus_solver].push_back(gen_id);
            }

            // ... and the storage units regulating their bus, which pin it the same way
            // (a storage unit only ever regulates its own, see LSGrid::check_grid). No row
            // places one (refused by _maybe_resolve_topology): a row only takes one out.
            const auto & storages = _grid_model.get_storages();
            std::map<int, std::vector<int> > storages_of_bus;  // solver bus -> its regulating storage units
            for(int storage_id = 0; storage_id < static_cast<int>(storages.nb()); ++storage_id){
                if(!storages.is_local_voltage_controller(storage_id)) continue;
                const int bus_id_me = storages.get_bus_id()(storage_id).cast_int();
                if(bus_id_me == BaseConstants::_deactivated_bus_id) continue;
                const int bus_solver = id_me_to_solver[bus_id_me].cast_int();
                if(bus_solver == BaseConstants::_deactivated_bus_id) continue;
                storages_of_bus[bus_solver].push_back(storage_id);
            }

            // ---- 3. per row: the buses that lose EVERY one of their controllers ----
            // Only a bus one of whose controllers the row takes out can: the buses of the
            // generators it has out and of the storage units it takes off, in ascending order
            _row_pv_to_pq_.assign(nb_steps, std::vector<int>());
            std::vector<int> candidate_buses;
            for(Eigen::Index step = 0; step < nb_rows; ++step){
                candidate_buses.clear();
                for(int gen_id : gens_off[static_cast<size_t>(step)]){
                    if(pinned_bus_of_gen[gen_id] >= 0) candidate_buses.push_back(pinned_bus_of_gen[gen_id]);
                }
                for(const auto & bus_storages : storages_of_bus){
                    for(int storage_id : bus_storages.second){
                        if(_storage_leaves_base_bus(static_cast<size_t>(step), storage_id)){
                            candidate_buses.push_back(bus_storages.first);
                            break;
                        }
                    }
                }
                std::sort(candidate_buses.begin(), candidate_buses.end());
                candidate_buses.erase(std::unique(candidate_buses.begin(), candidate_buses.end()), candidate_buses.end());
                for(int bus_solver : candidate_buses){
                    bool all_off = true;
                    const auto bus_gens = gens_of_bus.find(bus_solver);
                    if(bus_gens != gens_of_bus.end()){
                        for(int gen_id : bus_gens->second){
                            if(!_gen_off_in_row(static_cast<size_t>(step), gen_id)){
                                all_off = false;
                                break;
                            }
                        }
                    }
                    const auto bus_storages = storages_of_bus.find(bus_solver);
                    if(all_off && bus_storages != storages_of_bus.end()){
                        for(int storage_id : bus_storages->second){
                            if(!_storage_leaves_base_bus(static_cast<size_t>(step), storage_id)){
                                all_off = false;
                                break;
                            }
                        }
                    }
                    if(all_off) _row_pv_to_pq_[static_cast<size_t>(step)].push_back(bus_solver);
                }
            }

            // ---- 3b. per row: the buses a reactivated generator pins ----------------
            // A base-off generator the row puts back on its bus (a bus of the layout,
            // see _maybe_resolve_topology) regulating that bus locally makes it PV for
            // the row. A base-PQ bus already owns a Vm unknown and a Q equation, so
            // nothing is reserved: the row pins its Q row (set_pv_pinned_buses) and
            // seeds |V| there at the set-point (_row_vm_seed_). A base-PV bus stays PV
            // through it, held by the controllers the row leaves there -- provided they
            // ask the same magnitude -- or, when the row takes them all off, at the
            // incoming set-point.
            std::vector<std::set<int> > extra_pv(nb_steps);
            _row_vm_seed_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            // the magnitude a generator asks in a row: the row's own (modify_gen_v) or
            // the grid's
            auto gen_vm_in_row = [&](size_t step, int gen_id){
                const real_type own = _gen_target_vm_in_row(step, gen_id);
                return std::isfinite(own) ? own : generators.get_target_vm_pu(gen_id);
            };
            if(has_gens_on){
                const std::vector<int> & slack_solver = active_layout().slack_bus_id_solver.to_int_vector();
                const auto & topo_gens_on = _topo_gens_on_rows();
                for(size_t step = 0; step < nb_steps && step < topo_gens_on.size(); ++step){
                    std::vector<int> & to_pq = _row_pv_to_pq_[step];
                    for(const auto & gen_on : topo_gens_on[step]){
                        const int gen_id = gen_on.gen_id;
                        const int bus_solver = gen_on.bus_solver;
                        if(!generators.would_be_local_voltage_controller(gen_id)){
                            // does not regulate (or regulates remotely): an injection, the
                            // bus keeps its label. A remote one lives in the VoltageControl
                            // extension, out of scope here.
                            if(generators.regulates_remote(gen_id) && generators.get_voltage_regulator_on(gen_id)){
                                std::ostringstream exc_;
                                exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                     << " reactivates generator " << gen_id << ", which regulates the voltage "
                                        "of a remote bus. Only generators regulating their own bus can be "
                                        "reactivated for now.";
                                throw std::runtime_error(exc_.str());
                            }
                            continue;
                        }
                        const int bus_me = gen_on.bus_me;
                        if(group_buses.find(bus_me) != group_buses.end()){
                            std::ostringstream exc_;
                            exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                 << " reactivates generator " << gen_id << " on a bus whose voltage a control "
                                    "group holds (a remote generator, an SVC or an HVDC converter station): "
                                    "not supported yet.";
                            throw std::runtime_error(exc_.str());
                        }
                        if(std::find(slack_solver.begin(), slack_solver.end(), bus_solver) != slack_solver.end()){
                            std::ostringstream exc_;
                            exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                 << " reactivates generator " << gen_id << " on a slack bus: not supported yet.";
                            throw std::runtime_error(exc_.str());
                        }
                        // the magnitude this generator asks of the bus
                        const real_type vm = gen_vm_in_row(step, gen_id);
                        const auto base_gens = gens_of_bus.find(bus_solver);
                        const auto base_storages = storages_of_bus.find(bus_solver);
                        const bool base_pv = base_gens != gens_of_bus.end() || base_storages != storages_of_bus.end();
                        if(base_pv){
                            // a base-PV bus: kept PV whatever the row does to its own
                            // controllers. The first one the row leaves there holds it
                            // (they agree among themselves, _row_gen_v_conflicts), a
                            // generator or else a regulating storage unit; one the row
                            // takes off or moves away holds nothing.
                            std::vector<int>::iterator it = std::find(to_pq.begin(), to_pq.end(), bus_solver);
                            if(it != to_pq.end()) to_pq.erase(it);
                            const char * holder_kind = "generator";
                            int holder = -1;
                            real_type bus_vm = 0.;
                            if(base_gens != gens_of_bus.end()){
                                for(int other : base_gens->second){
                                    if(_gen_off_in_row(step, other)) continue;
                                    holder = other;
                                    bus_vm = gen_vm_in_row(step, other);
                                    break;
                                }
                            }
                            if(holder < 0 && base_storages != storages_of_bus.end()){
                                for(int storage_id : base_storages->second){
                                    if(_storage_leaves_base_bus(step, storage_id)) continue;
                                    holder_kind = "storage unit";
                                    holder = storage_id;
                                    bus_vm = storages.get_target_vm_pu(storage_id);
                                    break;
                                }
                            }
                            if(holder >= 0){
                                if(std::abs(bus_vm - vm) > BaseConstants::_tol_equal_float){
                                    std::ostringstream exc_;
                                    exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                         << " reactivates generator " << gen_id << " (set-point " << vm
                                         << " pu) on a bus " << holder_kind << " " << holder << " holds at " << bus_vm
                                         << " pu. A bus has one magnitude: give them the same set-point.";
                                    throw std::runtime_error(exc_.str());
                                }
                                continue;
                            }
                        }
                        // a bus nothing else holds in this row: at this magnitude. A
                        // base-PQ bus turns PV for the row; a base-PV bus whose own
                        // controllers the row all takes off is pinned already (kept off
                        // to_pq above), only its |V| is seeded
                        bool seeded = false;
                        for(const auto & seed : _row_vm_seed_[step]){
                            if(seed.first != bus_solver) continue;
                            seeded = true;
                            if(std::abs(seed.second - vm) > BaseConstants::_tol_equal_float){
                                std::ostringstream exc_;
                                exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                     << " reactivates two generators on one bus with different set-points ("
                                     << seed.second << " and " << vm << " pu). A bus has one magnitude.";
                                throw std::runtime_error(exc_.str());
                            }
                        }
                        if(!seeded) _row_vm_seed_[step].push_back(std::make_pair(bus_solver, vm));
                        if(!base_pv) extra_pv[step].insert(bus_solver);
                    }
                }
            }

            // ---- 3c. the union, and each row's pinned set ---------------------------
            std::set<int> switchable;
            for(size_t step = 0; step < nb_steps; ++step) switchable.insert(_row_pv_to_pq_[step].begin(), _row_pv_to_pq_[step].end());
            _switchable_buses_.assign(switchable.begin(), switchable.end());  // sorted (std::set)
            _row_pv_pinned_.assign(nb_steps, std::vector<int>());
            for(size_t step = 0; step < nb_steps; ++step){
                std::vector<int> & pinned = _row_pv_pinned_[step];
                std::sort(_row_pv_to_pq_[step].begin(), _row_pv_to_pq_[step].end());
                std::set_difference(_switchable_buses_.begin(), _switchable_buses_.end(),
                                    _row_pv_to_pq_[step].begin(), _row_pv_to_pq_[step].end(), std::back_inserter(pinned));
                pinned.insert(pinned.end(), extra_pv[step].begin(), extra_pv[step].end());
                std::sort(pinned.begin(), pinned.end());
                if(pinned != _switchable_buses_) _pv_pinning_active_ = true;
            }
            if(!_switchable_buses_.empty()) _pv_pinning_active_ = true;

            // ---- 4. reserve the Jacobian slots, and start pinned ------------------
            _push_switchable_to_algo();
        }

        // what the loads and shunts row i's contingency trips take out of the balance, as
        // power LOST (generator convention): minus their consumption, each counted once,
        // where it stood in the base grid (the slack pre-pass' term (e))
        real_type _row_trip_lost_mw(size_t i) const {
            real_type lost = 0.;
            const auto & loads = _grid_model.get_loads();
            for(int load_id : _row_tripped_loads(i)){
                if(!loads.get_status()[static_cast<size_t>(load_id)]) continue;
                if(loads.get_bus_id()(load_id).cast_int() == BaseConstants::_deactivated_bus_id) continue;
                lost -= _load_target_p_in_row(i, load_id);
            }
            const auto & shunts = _grid_model.get_shunts();
            for(int shunt_id : _row_tripped_shunts(i)){
                if(!shunts.get_status()[static_cast<size_t>(shunt_id)]) continue;
                if(shunts.get_bus_id()(shunt_id).cast_int() == BaseConstants::_deactivated_bus_id) continue;
                lost -= shunts.get_target_p()(shunt_id);
            }
            return lost;
        }

        // A shunt a row trips must not be one the ShuntVoltageControl outer loop moves:
        // its sections are that loop's (NRSystem's ShuntControl patches the diagonal by
        // the difference with what it believes is stamped there), and it does not know
        // of a row taking the shunt out. Refused, with an algorithm running outer loops
        // whose list holds that loop (OpenLoadFlow's default list does not: the sections
        // are then nobody's), once per batch (_prepare_ybus_varying has just listed the
        // rows' shunts).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _check_shunt_trips(size_t nb_steps) const {
            if(!_algo.supports_outer_loops()) return;
            bool shunt_loop = false;
            for(const auto & loop : _grid_model.get_outer_loops()){
                if(loop != nullptr && loop->name() == "ShuntVoltageControl") shunt_loop = true;
            }
            if(!shunt_loop) return;
            const ShuntContainer & shunts = _grid_model.get_shunts();
            for(size_t row = 0; row < nb_steps; ++row){
                for(int shunt_id : _row_tripped_shunts(row)){
                    if(!shunts.has_sections(shunt_id) || !shunts.get_section_regulating(shunt_id) ||
                       !shunts.get_status(shunt_id)) continue;
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::compute: row " << row << " trips shunt " << shunt_id
                         << ", which regulates a voltage, with an algorithm running outer loops ("
                         << _algo.get_name() << "). The ShuntVoltageControl loop owns its sections and does "
                            "not support a row taking it out yet: leave it in, or use an algorithm without "
                            "outer loops.";
                    throw std::runtime_error(exc_.str());
                }
            }
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _check_shunt_trips(size_t) const {}

        // whether some row's action takes out a storage unit that regulates its bus
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _has_regulating_storage_off() const {
            const auto & storages = _grid_model.get_storages();
            for(const std::vector<int> & off : sbus_policy_.topo_storages_off){
                for(int storage_id : off) if(storages.is_local_voltage_controller(storage_id)) return true;
            }
            return false;
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _has_regulating_storage_off() const { return false; }

        // Hand _switchable_buses_ to the member algorithm, pinned. Every switchable bus
        // is PV in the "n" case -- that is where its controller still stands -- so the
        // base solve pins all of them and the per-row loop releases the ones that flip.
        // Called on EVERY path out of _maybe_prepare_gen_contingency, the empty ones
        // included: _algo is a member, and what a previous compute() reserved on it has
        // to be taken back or the next one solves a bus as PQ that nothing pins.
        void _push_switchable_to_algo(){
            _algo.set_switchable_vm_buses(_switchable_buses_);
            _algo.set_pv_pinned_buses(_switchable_buses_);
            // same reason as in _maybe_prepare_masks: a Q row pinned to identity in
            // one row and live in the next is a pivot that may move
            if(_pv_pinning_active_) _algo.set_refactor_fallback(true);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_gen_contingency(size_t){}
        // the setter / method a generator contingency comes from, for the messages
        static const char * _gen_trip_api() {
            return SbusPolicy::supports_vary ? "set_contingency_gens / set_topo_actions" : "add_contingency_gen / add_contingency";
        }

        // true when this compute() has generator contingencies to apply at all; every
        // per-row helper below is skipped (and the loop stays bit-identical to before
        // this feature) when it is false.
        bool _has_gen_contingency() const { return _gen_contingency_active_; }
        // ... and true when some bus can actually flip PV -> PQ, so the per-row
        // set_pv_pinned_buses calls are worth making. False in DC, and false when the
        // masked generators happen never to leave a bus without a controller.
        bool _has_pv_switching() const { return _pv_pinning_active_; }

        // whether the grid has voltage-control groups (a remote regulator, an SVC): their
        // bordered rows put zeros on the diagonal, so the pivots KLU fixed at the "n"
        // factorization can be zero in a row's edited Jacobian (a contingency changes
        // Ybus), whatever the batch masks. Only meaningful once the plan is built.
        bool _has_voltage_control_groups() const {
            return _grid_model.get_ac_voltage_control_plan().controllers().n_groups() > 0;
        }

        // whether some row solves with slack weights of its own (_row_slack_weights): a
        // participant it disconnects, or one its slack pre-pass saturates
        bool _row_slack_weights_vary() const {
            const auto any_row = [](const std::vector<std::vector<int> > & rows){
                for(const std::vector<int> & row : rows) if(!row.empty()) return true;
                return false;
            };
            return any_row(_row_slack_gens_off_) || any_row(_row_slack_storages_off_) ||
                   any_row(_row_sat_gens_) || any_row(_row_sat_storages_);
        }

        // the buses to pin PV in row i: the switchable buses that are STILL PV there
        // (the complement, _row_pv_to_pq_[i], is what becomes PQ) plus the base-PQ
        // buses a reactivated generator pins. Sorted, precomputed per row.
        const std::vector<int> & _row_pv_pinned(size_t i) const {
            if(i >= _row_pv_pinned_.size()) return _switchable_buses_;
            return _row_pv_pinned_[i];
        }

        // |V| of the buses row i turns PV with a reactivated generator, seeded at that
        // generator's set-point before the solve (a pinned bus keeps its starting
        // magnitude). A no-op on every other row.
        void _apply_step_topo_seed(size_t i, Eigen::Ref<CplxVect> V) const {
            if(i >= _row_vm_seed_.size()) return;
            for(const auto & seed : _row_vm_seed_[i]){
                const int b = seed.first;
                if(b < 0 || b >= V.size()) continue;
                const real_type vm = std::abs(V(b));
                V(b) = vm > 0. ? V(b) * (seed.second / vm) : cplx_type(seed.second, 0.);
            }
        }

        // row i's distributed-slack weights: the layout's own unless the row takes a
        // participating generator or storage unit out, in which case they are
        // re-derived without it and renormalised (LSGrid::get_slack_weights_solver_without).
        // Should a row somehow leave no participant at all, the reference slack bus
        // keeps the whole share -- the angle reference is a property of the batch,
        // picked once, and must not move from row to row.
        // Returns the layout's own vector by reference on the common row; the
        // re-derived one is written into the caller's `scratch` (one per range, so
        // no row allocates and no two threads share it).
        const RealVect & _row_slack_weights(size_t i, RealVect & scratch) const {
            const RealVect & base_w = active_layout().slack_weights;
            const bool has_gens_off = (i < _row_slack_gens_off_.size()) && !_row_slack_gens_off_[i].empty();
            // ... and the units the slack pre-pass saturated take no share either
            // (see _prepare_slack_redistribution)
            const bool has_sat_gens = (i < _row_sat_gens_.size()) && !_row_sat_gens_[i].empty();
            const bool has_sat_storages = (i < _row_sat_storages_.size()) && !_row_sat_storages_[i].empty();
            const bool has_storages_off = (i < _row_slack_storages_off_.size()) && !_row_slack_storages_off_[i].empty();
            if(!has_gens_off && !has_sat_gens && !has_sat_storages && !has_storages_off) return base_w;
            const auto & generators = _grid_model.get_generators();
            std::vector<bool> gen_off(generators.nb(), false);
            if(has_gens_off) for(int gen_id : _row_slack_gens_off_[i]) gen_off[gen_id] = true;
            if(has_sat_gens) for(int gen_id : _row_sat_gens_[i]) gen_off[gen_id] = true;
            std::vector<bool> storage_off;
            if(has_sat_storages || has_storages_off){
                storage_off.assign(_grid_model.get_storages().nb(), false);
                if(has_sat_storages) for(int storage_id : _row_sat_storages_[i]) storage_off[storage_id] = true;
                if(has_storages_off) for(int storage_id : _row_slack_storages_off_[i]) storage_off[storage_id] = true;
            }
            scratch = _grid_model.get_slack_weights_solver_without(
                static_cast<size_t>(base_w.size()), active_layout().id_me_to_solver, gen_off, storage_off);
            if(abs(scratch.sum()) < BaseConstants::_tol_equal_float){
                scratch.setZero();
                if(active_layout().slack_bus_id_solver.size() > 0){
                    const int ref = active_layout().slack_bus_id_solver[static_cast<int>(0)].cast_int();
                    if(ref >= 0 && ref < scratch.size()) scratch(ref) = 1.;
                }
            }
            return scratch;
        }

        // ---- OLF-style bounded slack redistribution (option `redistribute_slack`) ----
        //
        // The distributed slack of the solve shares whatever imbalance the row leaves
        // by fixed per-bus weights, with no limit: a row that loses a big generator (a
        // generator contingency, or an island cut off in "handle disconnected grid"
        // mode) pushes the remaining machines past their max_p, where OpenLoadFlow's
        // DistributedSlack outer loop stops each one at its bound and re-shares the
        // excess. The part of that imbalance known BEFORE the solve -- the set-points of
        // what the row takes out -- is shared here the way OLF does it
        // (slack_redistribution::distribute), once per compute() for every row:
        //   _row_slack_dp_pu_  : per row, the correction to add to the row's injection
        //                        (solver bus, dP in pu), one entry per participant bus
        //   _row_gen_new_p_ /
        //   _row_sto_new_p_    : per row, the units whose set-point moved (id, new
        //                        injection in MW, generator convention), sorted by id --
        //                        what the p-limit check reads as the row's target
        //   _row_sat_gens_ /
        //   _row_sat_storages_ : per row, the units that reached a bound (sorted): out
        //                        of that row's distributed slack (_row_slack_weights),
        //                        so the solve only shares what is left (the change in
        //                        the losses) on the units that can still move
        // With outer loops (NROuter_*) none of these: the Newton has a single slack, the
        // DistributedSlack loop owns the participants and shares the loss itself, before
        // the row's first solve, from the targets in force (OuterRowView::lost_p_mw):
        //   _row_lost_p_mw_    : per row, what it loses (MW), NaN where nothing is shared
        // Nothing the size of nb_rows x nb_bus exists; a row that loses nothing has
        // every entry empty and solves exactly as before.
        void _clear_slack_redistribution(){
            _row_slack_dp_pu_.clear();
            _row_gen_new_p_.clear();
            _row_sto_new_p_.clear();
            _row_sat_gens_.clear();
            _row_sat_storages_.clear();
            _row_lost_p_mw_.clear();
        }

        // this row's own active set-point of a static generator / a load (MW): its
        // own row where modify_sgen_p / modify_load_p was given one, the grid's target
        // otherwise (see _gen_target_p_in_row)
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _sgen_target_p_in_row(size_t i, int sgen_id) const {
            const auto & mat = sbus_policy_.sgen_p;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(mat.rows() > 0 && row < mat.rows() && sgen_id < mat.cols()) return mat(row, sgen_id);
            return _grid_model.get_sgen_target_p()(sgen_id);
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _sgen_target_p_in_row(size_t, int sgen_id) const { return _grid_model.get_sgen_target_p()(sgen_id); }
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _load_target_p_in_row(size_t i, int load_id) const {
            const auto & mat = sbus_policy_.load_p;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(mat.rows() > 0 && row < mat.rows() && load_id < mat.cols()) return mat(row, load_id);
            return _grid_model.get_load_target_p()(load_id);
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _load_target_p_in_row(size_t, int load_id) const { return _grid_model.get_load_target_p()(load_id); }

        // what row i's own injections (modify_gen_p / modify_sgen_p / modify_load_p)
        // take out of the active balance of the grid's targets, in MW and generator
        // convention: minus SbusPolicy::Vary::row_p_change_mw. 0 where the injection
        // does not vary (ContingencyAnalysis).
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        real_type _row_injection_change_lost_mw(size_t i) const {
            return -sbus_policy_.row_p_change_mw(static_cast<Eigen::Index>(i));
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        real_type _row_injection_change_lost_mw(size_t) const { return 0.; }

        // the grid's own target of a slack unit, generator convention (a storage unit's
        // target is in load convention)
        real_type _grid_target_p_of(ViolationElementType el_type, int el_id) const {
            if(el_type == ViolationElementType::GENERATOR) return _grid_target_p(el_id);
            const Eigen::Ref<const RealVect> tgt = _grid_model.get_storage_target_p();
            return (el_id >= 0 && el_id < tgt.size()) ? -tgt(el_id) : 0.;
        }

        static const real_type * _find_row_new_p(const std::vector<std::vector<std::pair<int, real_type> > > & rows,
                                                 size_t i, int el_id){
            if(i >= rows.size()) return nullptr;
            const std::vector<std::pair<int, real_type> > & row = rows[i];
            auto it = std::lower_bound(row.begin(), row.end(), std::make_pair(el_id, real_type(0.)),
                                       [](const std::pair<int, real_type> & a, const std::pair<int, real_type> & b){
                                           return a.first < b.first; });
            if(it == row.end() || it->first != el_id) return nullptr;
            return &it->second;
        }

        // this row's active set-point of a slack unit (MW, generator convention): the
        // one the slack pre-pass wrote where it moved it, the row's own otherwise
        real_type _row_target_p(size_t i, ViolationElementType el_type, int el_id) const {
            // the target this row's outer loops moved it to (DistributedSlack), first
            if(i < _outer_results_.size()){
                const OuterRowResult & r = _outer_results_[i];
                const auto & moved_by_loops = (el_type == ViolationElementType::GENERATOR) ? r.gen_target_p : r.storage_target_p;
                auto it = std::lower_bound(moved_by_loops.begin(), moved_by_loops.end(), std::make_pair(el_id, real_type(0.)),
                                           [](const std::pair<int, real_type> & a, const std::pair<int, real_type> & b){ return a.first < b.first; });
                if(it != moved_by_loops.end() && it->first == el_id){
                    // a storage unit's, in the load convention: its injection is the opposite
                    return (el_type == ViolationElementType::GENERATOR) ? it->second : -it->second;
                }
            }
            if(el_type == ViolationElementType::GENERATOR){
                const real_type * moved = _find_row_new_p(_row_gen_new_p_, i, el_id);
                return moved != nullptr ? *moved : _gen_target_p_in_row(i, el_id);
            }
            const real_type * moved = _find_row_new_p(_row_sto_new_p_, i, el_id);
            return moved != nullptr ? *moved : _grid_target_p_of(el_type, el_id);
        }

        // whether the slack pre-pass took that unit out of this row's distributed slack
        bool _row_takes_no_share(size_t i, ViolationElementType el_type, int el_id) const {
            const std::vector<std::vector<int> > & rows =
                (el_type == ViolationElementType::GENERATOR) ? _row_sat_gens_ : _row_sat_storages_;
            if(i >= rows.size()) return false;
            return std::binary_search(rows[i].begin(), rows[i].end(), el_id);
        }

        // Row i's injection WITH the slack pre-pass correction: `_step_sbus` (the
        // row's own, or the fixed base vector by reference) plus this row's dP on the
        // participants' buses, written into `scratch`. In DC on a fixed-injection sweep
        // the base vector is the (real) dc_cache_.inj, which compute_one_powerflow
        // would otherwise fall back to on an empty Sbus.
        const CplxVect & _step_sbus_row(size_t i, CplxVect & scratch) const {
            const CplxVect & base = _step_sbus(i, scratch);
            const bool has_slack_dp = i < _row_slack_dp_pu_.size() && !_row_slack_dp_pu_[i].empty();
            const std::vector<std::pair<int, cplx_type> > * trip_dp = _row_sbus_dp(i);
            if(!has_slack_dp && trip_dp == nullptr) return base;
            if(&base != &scratch){
                if(base.size() > 0) scratch = base;
                else scratch = dc_cache_.inj.template cast<cplx_type>();
            }
            if(has_slack_dp){
                for(const std::pair<int, real_type> & bus_dp : _row_slack_dp_pu_[i]){
                    if(bus_dp.first < 0 || bus_dp.first >= scratch.size()) continue;
                    scratch(bus_dp.first) += cplx_type(bus_dp.second, 0.);
                }
            }
            if(trip_dp != nullptr){
                for(const std::pair<int, cplx_type> & bus_dp : *trip_dp){
                    if(bus_dp.first < 0 || bus_dp.first >= scratch.size()) continue;
                    scratch(bus_dp.first) += bus_dp.second;
                }
            }
            return scratch;
        }

        // what row i's contingency adds to its injection (see
        // YbusPolicy::Contingency::li_sbus_dp: the DC phase shifters and shunts it
        // disconnects, the loads and generators a ContingencyAnalysis contingency
        // trips), or nullptr if nothing -- always the case where Ybus does not vary
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        const std::vector<std::pair<int, cplx_type> > * _row_sbus_dp(size_t i) const {
            if(i >= ybus_policy_.li_sbus_dp.size() || ybus_policy_.li_sbus_dp[i].empty()) return nullptr;
            return &ybus_policy_.li_sbus_dp[i];
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        const std::vector<std::pair<int, cplx_type> > * _row_sbus_dp(size_t) const { return nullptr; }

        // The pre-pass itself, once per compute() and after the per-row injections are
        // known (_prepare_sbus_varying) and the masks / generator contingencies settled.
        // For each row: what it loses, in MW and generator convention, is
        //   (a) the generators it disconnects (a ScenarioSweep generator contingency,
        //       slack or not): their ROW set-point (SbusPolicy::Vary::fill_row already
        //       takes them out of the row's injection, so their power IS what the
        //       remaining machines have to make up), plus
        //   (b) the elements stranded on the buses it masks, element by element
        //       (generators not already counted in (a), static generators, minus loads,
        //       storage units and shunts, plus the HVDC converters stranded there: the
        //       in-main converter of their line keeps injecting, see HvdcLineContainer,
        //       so what the balance loses is the stranded station's own setpoint) --
        //       their rows are masked in the solve, so the balance loses their net
        //       injection. Element-wise rather than the real part of the row's
        //       injection over the masked buses, which in DC would count the
        //       phase-shifter term, plus
        //   (c) what the row's own injections take out of the balance of the grid's
        //       targets (a ScenarioSweep's modify_gen_p / modify_sgen_p /
        //       modify_load_p): sum(target - row) over the generators and static
        //       generators minus the same over the loads, every element of the solved
        //       system counted -- those (a) and (b) then take out included, at their
        //       row value, so nothing is counted twice, plus
        //   (d) what the row's topological action changes (set_topo_actions, see
        //       _row_topo_lost_mw): every unit it touches is counted where the ROW puts
        //       it, so (a) and (b) leave out a generator it moves and a load / storage
        //       unit it takes off its base bus. A unit moved inside the main component
        //       loses nothing, nor does a bus the row leaves empty: what stood there
        //       stands elsewhere.
        // A line or transformer that leaves the grid connected loses no injection:
        // the change in the losses it causes stays with the solve's distributed slack.
        // The participants are the slack units left in the main component and not
        // disconnected by the row, and the units flagged "can participate in the slack"
        // (see SlackParticipation::set_can_participate), each on the bus the row gives
        // it; slack_redistribution::distribute does the rest.
        void _prepare_slack_redistribution(size_t nb_steps){
            _clear_slack_redistribution();
            if(!_redistribute_slack_) return;
            using slack_redistribution::Participant;
            using slack_redistribution::UnitKind;

            _row_slack_dp_pu_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            _row_gen_new_p_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            _row_sto_new_p_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            _row_sat_gens_.assign(nb_steps, std::vector<int>());
            _row_sat_storages_.assign(nb_steps, std::vector<int>());
            const bool outer = _algo.supports_outer_loops();
            if(outer) _row_lost_p_mw_.assign(nb_steps, std::numeric_limits<real_type>::quiet_NaN());

            const auto & generators = _grid_model.get_generators();
            const auto & sgens = _grid_model.get_static_generators();
            const auto & loads = _grid_model.get_loads();
            const auto & storages = _grid_model.get_storages();
            const auto & shunts = _grid_model.get_shunts();
            const SolverBusIdVect & id2s = active_layout().id_me_to_solver;
            const real_type sn_mva = _grid_model.get_sn_mva();
            const real_type eps_mw = slack_redistribution::default_eps_mw;
            const int nb_gen = generators.nb();
            const std::vector<bool> & gen_status = generators.get_status();
            const GlobalBusIdVect & gen_bus = generators.get_bus_id();

            std::vector<Participant> units;
            std::vector<real_type> new_inj;
            std::vector<char> saturated;
            for(size_t i = 0; i < nb_steps; ++i){
                if(i < _skip_mask.size() && _skip_mask[i]) continue;
                const std::vector<int> * masked = _row_masked_ids(i);
                const auto in_island_solver = [masked](int bus_solver){
                    if(masked == nullptr || bus_solver < 0) return false;
                    return std::binary_search(masked->begin(), masked->end(), bus_solver);
                };
                const auto in_island = [&](int bus_me){ return in_island_solver(id2s[bus_me].cast_int()); };
                const auto gen_off = [this, i](int gen_id){ return this->_gen_off_in_row(i, gen_id); };
                const auto gen_p_row = [this, i](int gen_id){ return this->_gen_target_p_in_row(i, gen_id); };

                // (a) the generators this row disconnects (one its action moves is off its
                // base bus too, and counted in (d) where it lands)
                real_type lost_mw = 0.;
                for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                    if(!gen_status[gen_id]) continue;
                    if(gen_bus(gen_id).cast_int() == BaseConstants::_deactivated_bus_id) continue;
                    if(!gen_off(gen_id)) continue;
                    if(_gen_placed_in_row(i, gen_id)) continue;
                    lost_mw += gen_p_row(gen_id);
                }
                // (b) the elements stranded on the masked buses
                if(masked != nullptr){
                    lost_mw += slack_redistribution::sum_setpoints_if(
                        generators, 1., in_island,
                        [&](int gen_id){ return gen_off(gen_id) ? 0. : gen_p_row(gen_id); });
                    lost_mw += slack_redistribution::sum_setpoints_if(
                        sgens, 1., in_island, [this, i](int sgen_id){ return this->_sgen_target_p_in_row(i, sgen_id); });
                    lost_mw += slack_redistribution::sum_setpoints_if(
                        loads, -1., in_island, [this, i](int load_id){
                            return this->_load_leaves_base_bus(i, load_id) ? 0. : this->_load_target_p_in_row(i, load_id); });
                    lost_mw += slack_redistribution::sum_setpoints_if(
                        storages, -1., in_island, [this, i, &storages](int storage_id){
                            return this->_storage_leaves_base_bus(i, storage_id) ? 0. : storages.get_target_p()(storage_id); });
                    lost_mw += slack_redistribution::sum_setpoints_if(
                        shunts, -1., in_island, [this, i, &shunts](int shunt_id){
                            return this->_shunt_off_in_row(i, shunt_id) ? 0. : shunts.get_target_p()(shunt_id); });
                    lost_mw += slack_redistribution::sum_hvdc_station_setpoints_if(
                        _grid_model.get_dclines(), in_island);
                }
                // (c) the row's own injection change
                lost_mw += _row_injection_change_lost_mw(i);
                // (d) what the row's topological action moves, disconnects or reactivates
                lost_mw += _row_topo_lost_mw(i, in_island_solver);
                // (e) the loads and shunts the row's contingency trips (its generators are
                // (a)): their consumption leaves the balance, stranded or not ((b) leaves
                // them out)
                lost_mw += _row_trip_lost_mw(i);
                if(outer){
                    // the DistributedSlack loop shares it (see OuterRowView::lost_p_mw). Even
                    // when it is nothing at the grid's targets: the loop counts the units the
                    // row loses at the targets in force, which a row resumed from the "n"
                    // solve takes from that solve (a unit at 0 MW in the grid that the "n"
                    // solve's DistributedSlack moved up loses that much); it then decides
                    // itself whether anything is left to share
                    _row_lost_p_mw_[i] = lost_mw;
                    continue;
                }
                if(std::abs(lost_mw) <= eps_mw) continue;

                units.clear();
                const auto keep_bus = [&in_island](int bus_me){ return !in_island(bus_me); };
                slack_redistribution::append_participants(
                    generators, UnitKind::GENERATOR, 1., keep_bus, gen_off, gen_p_row, units);
                _append_topo_gen_participants(i, in_island_solver, units);
                slack_redistribution::append_participants(
                    storages, UnitKind::STORAGE, -1., keep_bus,
                    [this, i](int storage_id){ return this->_storage_leaves_base_bus(i, storage_id); },
                    [&storages](int storage_id){ return storages.get_target_p()(storage_id); }, units);
                _append_topo_storage_participants(i, in_island_solver, units);
                if(units.empty()) continue;  // the row's own fallback stays (see _row_slack_weights)

                const slack_redistribution::Report report = slack_redistribution::distribute(
                    units, lost_mw, eps_mw, new_inj, saturated);

                std::vector<std::pair<int, real_type> > & dp_row = _row_slack_dp_pu_[i];
                for(size_t k = 0; k < units.size(); ++k){
                    const real_type dp_mw = new_inj[k] - units[k].injection_mw;
                    if(std::abs(dp_mw) > BaseConstants::_tol_equal_float){
                        // the participants are listed by id (generators, then storage
                        // units), so each per-family list comes out sorted
                        if(units[k].kind == UnitKind::GENERATOR) _row_gen_new_p_[i].push_back(std::make_pair(units[k].el_id, new_inj[k]));
                        else _row_sto_new_p_[i].push_back(std::make_pair(units[k].el_id, new_inj[k]));
                        const int bus_solver = id2s[units[k].bus].cast_int();
                        if(bus_solver >= 0) dp_row.push_back(std::make_pair(bus_solver, dp_mw / sn_mva));
                    }
                    // a unit only flagged "can participate in the slack" was never in it
                    if(saturated[k] && units[k].in_slack){
                        if(units[k].kind == UnitKind::GENERATOR) _row_sat_gens_[i].push_back(units[k].el_id);
                        else _row_sat_storages_[i].push_back(units[k].el_id);
                    }
                }
                // every unit of the slack saturated: `saturated` is all zero, they all stay in
                // this row's slack (see slack_redistribution::distribute)
                (void) report;
                // one entry per bus: merge the units sharing one
                if(dp_row.size() > 1){
                    std::sort(dp_row.begin(), dp_row.end(),
                              [](const std::pair<int, real_type> & a, const std::pair<int, real_type> & b){
                                  return a.first < b.first; });
                    std::vector<std::pair<int, real_type> > merged;
                    merged.reserve(dp_row.size());
                    for(const std::pair<int, real_type> & bus_dp : dp_row){
                        if(!merged.empty() && merged.back().first == bus_dp.first) merged.back().second += bus_dp.second;
                        else merged.push_back(bus_dp);
                    }
                    dp_row.swap(merged);
                }
            }
        }

        // whether row i turns some switchable bus PQ -- only then is the algorithm's
        // pinning (all switchable buses, its resting state) worth touching: every
        // set_pv_pinned_buses call marks the mask positions dirty, which costs a pass
        // over the Jacobian's nonzeros at the next fill.
        bool _row_flips_pv(size_t i) const {
            return _has_pv_switching() && i < _row_pv_pinned_.size() && _row_pv_pinned_[i] != _switchable_buses_;
        }

        // per-range worker: NON mask-mode path, shared by every instantiation.
        void _run_range(size_t step_begin, size_t step_end,
                        AlgorithmSelector & algo, AlgoControl & control,
                        Eigen::SparseMatrix<cplx_type> & Ybus,
                        const Eigen::Ref<const CplxVect> & Vinit_solver,
                        bool ac_solver_used, int max_iter, real_type tol_solver,
                        int & nb_solved, int & nb_converged, double & timer_solver, double & timer_modif_ybus,
                        int & first_diverging_step, std::exception_ptr & err,
                        bool needs_solver_init);
        void _run_one_step(size_t i, AlgorithmSelector & algo, AlgoControl & control,
                           Eigen::SparseMatrix<cplx_type> & Ybus, CplxVect & V,
                           CplxVect & sbus_scratch, RealVect & sw_scratch,
                           bool ac_solver_used, int max_iter, real_type tol_solver,
                           int & nb_solved, int & nb_converged, double & timer_solver, double & timer_modif_ybus,
                           bool & conv, bool & invertible, OuterRowView * row_view);

        // per-range worker: mask-mode path (ContingencyAnalysis AND ScenarioSweep --
        // shared, using _step_sbus(cont_id) instead of a hardcoded fixed ac_cache_.inj so it
        // picks up per-row injections on ScenarioSweep; a no-op change for
        // ContingencyAnalysis, whose overload of _step_sbus still returns ac_cache_.inj).
        // Only ever called from _maybe_run_range_masked below, itself only ever
        // called from compute() (non-template, .cpp-confined) -- safe to stay
        // declared here / defined out-of-line in BaseBatchSweep.cpp.
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _run_range_masked(size_t cont_begin, size_t cont_end,
                               AlgorithmSelector & algo, AlgoControl & control,
                               Eigen::SparseMatrix<cplx_type> & Ybus,
                               const Eigen::Ref<const CplxVect> & Vinit_solver,
                               bool ac_solver_used, int max_iter, real_type tol, real_type sn_mva,
                               double & timer_modif_ybus, int & nb_solved, int & nb_converged, double & timer_solver,
                               int & first_diverging_step, std::exception_ptr & err, bool needs_solver_init)
        {
            // the row as the outer loops see it, for an algorithm that runs them
            OuterRowView row_view;
            const bool outer = algo.supports_outer_loops();
            try {
                CplxVect V;
                CplxVect sbus_scratch;   // this row's injection, where it varies
                RealVect sw_scratch;
                if(needs_solver_init) control.tell_all_changed();
                // the loop's invariant: between two rows the algorithm masks nothing
                // and pins every switchable bus, so a row that strands nothing and
                // flips nothing leaves both alone (each call would cost a pass over the
                // Jacobian's nonzeros at the next fill). Established once here rather
                // than assumed of a member algorithm a previous compute() may have left
                // mid-row.
                algo.set_masked_buses(_base_masked_);
                if(_has_pv_switching()) algo.set_pv_pinned_buses(_switchable_buses_);

                for(size_t cont_id = cont_begin; cont_id < cont_end; ++cont_id){
                    const std::vector<Coeff> & coeffs_modif = ybus_policy_.li_coeffs[cont_id];
                    bool conv = false;
                    bool do_store = false;
                    LimitViolationType div_reason = LimitViolationType::NOT_SIMULATED;

                    // same as _run_one_step: a row whose gen_v set-points contradict
                    // each other asks for a magnitude no solution can hold, so it is
                    // skipped before the solver rather than silently resolved by
                    // whichever generator set_vm visits last (_row_gen_v_conflicts)
                    if(!_skip_mask[cont_id] && !_row_gen_v_conflicts(cont_id)){
                        const std::vector<int> & masked = _li_masked[cont_id];
                        auto t1 = CustTimer();
                        YbusPolicy::Contingency::remove_from_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                        timer_modif_ybus += t1.duration();

                        const bool own_mask = masked != _base_masked_;
                        if(own_mask) algo.set_masked_buses(masked);
                        // generator contingencies: release the pinning of the buses
                        // this row turns PQ, keep it on the others
                        const bool flips = _row_flips_pv(cont_id);
                        if(flips) algo.set_pv_pinned_buses(_row_pv_pinned(cont_id));
                        V = Vinit_solver;
                        _apply_step_gen_v(cont_id, V);
                        _apply_step_topo_seed(cont_id, V);
                        _apply_step_dc_vm_reset(cont_id, V);
                        _apply_step_vc_v_set(cont_id, algo);
                        const RealVect & sw = _masked_slack_weights(masked, _row_slack_weights(cont_id, sw_scratch), sw_scratch);
                        const CplxVect & sb = _step_sbus_row(cont_id, sbus_scratch);
                        // a DC row's Pbus is its own (injections, slack pre-pass, phase
                        // shifters it disconnects), and the DC algorithm re-reads it only
                        // when told to: before EVERY row, the first one included -- the
                        // algorithm last saw the base case's
                        if(!ac_solver_used) control.tell_recompute_sbus();
                        if(outer){
                            _fill_outer_row_view(cont_id, &masked, row_view);
                            _set_outer_row_resume(cont_id, row_view);
                            algo.set_outer_row_view(&row_view);
                        }
                        const LinearSolverStats lin_before = algo.get_linear_solver_stats_before_solve();
                        conv = compute_one_powerflow(algo, control, nb_solved, nb_converged, timer_solver, Ybus, V, sb,
                                                     active_layout().slack_bus_id_solver.as_eigen(), sw,
                                                     active_layout().bus_pv.as_eigen(), active_layout().bus_pq.as_eigen(), max_iter, tol / sn_mva);
                        _record_row_solve_stats(cont_id, algo, lin_before);
                        if(needs_solver_init){ control.tell_none_changed(); needs_solver_init = false; }
                        // before the two restores below, and before the Ybus is put
                        // back: see _maybe_store_jacobian (and _record_row_physical, which
                        // reads the mismatch of the system this row solved)
                        if(conv){
                            _maybe_store_jacobian(cont_id, algo);
                            _record_row_physical(cont_id, algo, V, sw, sb);
                            _record_row_gen_results(cont_id, algo, Ybus, V, sw, sb);
                        }
                        if(own_mask) algo.set_masked_buses(_base_masked_);
                        if(flips) algo.set_pv_pinned_buses(_switchable_buses_);

                        auto t2 = CustTimer();
                        YbusPolicy::Contingency::readd_to_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                        timer_modif_ybus += t2.duration();

                        if(conv){
                            for(int b : masked) if(b >= 0 && b < V.size()) V(b) = cplx_type(0., 0.);
                            do_store = true;
                        } else {
                            div_reason = LimitViolationType::DIVERGENCE;
                        }
                    } else {
                        div_reason = LimitViolationType::NOT_SIMULATED;
                    }

                    if(do_store) _voltages.row(cont_id)(active_layout().id_solver_to_me.as_eigen()) = V.array();
                    else if(first_diverging_step < 0) first_diverging_step = static_cast<int>(cont_id);

                    _converged_mask_[cont_id] = do_store ? 1 : 0;
                    if(_compute_limit_violations_){
                        _converged[cont_id] = do_store ? 1 : 0;
                        if(!do_store){
                            _violations[cont_id].push_back(LimitViolation{
                                ViolationElementType::GRID, -1, 0, div_reason,
                                std::numeric_limits<real_type>::quiet_NaN(),
                                std::numeric_limits<real_type>::quiet_NaN()});
                        } else {
                            _record_row_violations(cont_id, V, &_li_masked[cont_id]);
                        }
                    }
                }
            } catch(...) {
                err = std::current_exception();
            }
            // the view dies here: the algorithm must not keep pointing at it
            if(outer) algo.set_outer_row_view(&_base_row_view_);
        }
        // dispatch into _run_range_masked only where it exists (ContingencyAnalysis
        // AND ScenarioSweep); unreachable on TimeSeries/InjectionSweep
        // (_handle_disconnected_grid can never be true there -- no setter reaches it).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _maybe_run_range_masked(size_t cont_begin, size_t cont_end,
                                     AlgorithmSelector & algo, AlgoControl & control,
                                     Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const Eigen::Ref<const CplxVect> & Vinit_solver,
                                     bool ac_solver_used, int max_iter, real_type tol, real_type sn_mva,
                                     double & timer_modif_ybus, int & nb_solved, int & nb_converged, double & timer_solver,
                                     int & first_diverging_step, std::exception_ptr & err, bool needs_solver_init){
            _run_range_masked(cont_begin, cont_end, algo, control, Ybus, Vinit_solver, ac_solver_used,
                              max_iter, tol, sn_mva, timer_modif_ybus, nb_solved, nb_converged, timer_solver,
                              first_diverging_step, err, needs_solver_init);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_run_range_masked(size_t, size_t, AlgorithmSelector &, AlgoControl &,
                                     Eigen::SparseMatrix<cplx_type> &, const Eigen::Ref<const CplxVect> &,
                                     bool, int, real_type, real_type, double &, int &, int &, double &,
                                     int &, std::exception_ptr &, bool){
            throw std::logic_error("unreachable: handle_disconnected_grid cannot be set on this instantiation");
        }

        // ----- threading spawn: 2-way (per-thread Ybus copy vs shared ac_cache_.mat) ------
        // Only ever called from compute() (non-template, .cpp-confined): safe to
        // stay declared here / defined out-of-line in BaseBatchSweep.cpp.
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _compute_threaded(size_t nb_steps, const CplxVect & Vinit_solver, bool ac_solver_used,
                               int max_iter, real_type tol, real_type sn_mva, double & timer_thread_init);
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _compute_threaded(size_t nb_steps, const CplxVect & Vinit_solver, bool ac_solver_used,
                               int max_iter, real_type tol, real_type sn_mva, double & timer_thread_init);

        // ----- flow-cleanup dispatch (zero the flow through a disconnected branch)-
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_check_results_match_defaults(const std::string &) {}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _maybe_check_results_match_defaults(const std::string & fun_name) const {
            const size_t nb_res = static_cast<size_t>(_nb_result_rows());
            if(nb_res != ybus_policy_.li_defaults.size()){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": the results were computed for "
                     << nb_res << " contingency(ies) but the object now holds " << ybus_policy_.li_defaults.size()
                     << ". The contingency set was modified (add_n1 / remove_n1 / ...) after the "
                        "last compute(); call compute() again before reading the flows.";
                throw std::runtime_error(exc_.str());
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _maybe_check_results_match_defaults(const std::string & fun_name) const {
            if(_results_stale_){
                std::ostringstream exc_;
                exc_ << algo_name() << "::" << fun_name << ": the injections and/or the contingency "
                        "masks were modified (modify_gen_p / modify_sgen_p / modify_load_p / "
                        "modify_load_q / set_contingency_lines / set_contingency_trafos / "
                        "set_contingency_gens / set_contingency_shunts / set_contingency_loads / "
                        "set_topo_actions) after the "
                        "last compute(); call compute() again before reading the flows.";
                throw std::runtime_error(exc_.str());
            }
        }

        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_clean_flows(bool) {}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _maybe_clean_flows(bool is_amps){
            auto timer = CustTimer();
            size_t cont_id = 0;
            for(const auto & this_cont: ybus_policy_.li_defaults){
                // a branch opened at one end only reads 0 too (see add_contingency)
                for(auto l_id : this_cont.branches_out()){
                    real_type & el = is_amps ? _amps_flows(cont_id, l_id): _active_power_flows(cont_id, l_id);
                    if(isfinite(el)) el = 0.;
                }
                ++cont_id;
            }
            if (is_amps) _timer_compute_A += timer.duration();
            else _timer_compute_P += timer.duration();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _maybe_clean_flows(bool is_amps){
            auto timer = CustTimer();
            const Eigen::Index nb_steps = _nb_result_rows();
            // the masks and the row's topological action alike (branch_ids_for_row
            // merges them): a branch the row disconnected carries no flow
            for(Eigen::Index row = 0; row < nb_steps; ++row){
                for(int br_id : ybus_policy_.branch_ids_for_row(row, n_line_)){
                    real_type & el = is_amps ? _amps_flows(row, br_id) : _active_power_flows(row, br_id);
                    if(isfinite(el)) el = 0.;
                }
            }
            if (is_amps) _timer_compute_A += timer.duration();
            else _timer_compute_P += timer.duration();
        }

    private:
        YbusPolicy ybus_policy_;
        SbusPolicy sbus_policy_;

        // aggregate status (1: success, 0: failure); meaningful only where
        // SbusPolicy::supports_vary (get_status() is gated accordingly). Plain,
        // always-present -- a single int costs nothing on the other instantiations.
        int _status = 1;

        // row-count lock shared by every setter (modify_* on the Sbus side,
        // set_contingency_lines/trafos on the Ybus side): -1 = not yet established.
        Eigen::Index _nb_steps_locked = -1;

        // set by every row-locked setter (_lock_or_check_nb_steps -- modify_gen_p /
        // modify_sgen_p / modify_load_p / modify_load_q / set_contingency_lines /
        // set_contingency_trafos), cleared on a successful compute(): guards
        // compute_flows()/compute_power_flows() on ScenarioSweep (Contingency && Vary)
        // against reading _voltages left over from a PREVIOUS compute() while
        // _maybe_clean_flows zeroes flows according to a mask/injection that was
        // re-set (same row count, different content) since -- unlike
        // ContingencyAnalysis's own _maybe_check_results_match_defaults, a row-count
        // comparison cannot catch this here: set_contingency_lines/trafos can only
        // ever replace a row's content, the row-count lock already forbids changing
        // the row count itself. Plain, always-present; only read on
        // (Contingency, Vary) (see _maybe_check_results_match_defaults below).
        bool _results_stale_ = false;

        // "handle disconnected grid" mode + limit violations: plain, always-present
        // state (SFINAE-gated methods above restrict who can reach it); empty /
        // false on every instantiation but ContingencyAnalysis.
        bool _handle_disconnected_grid = false;
        // whether _algo (the single-threaded / member algo -- each per-thread one is
        // freshly spawned and told directly, see _compute_threaded) has already been
        // told may_mask_voltage_control(true): sparsity only needs a forced rebuild
        // (tell_pv_changed()) the first time a compute() actually turns masking on,
        // see _maybe_prepare_masks().
        bool _voltage_control_may_mask_wired_ = false;
        // generator contingencies (ScenarioSweep only; plain, always-present state,
        // like everything else here -- the SFINAE-gated methods above restrict who
        // can fill it, and it stays empty on every other instantiation). Rebuilt by
        // _maybe_prepare_gen_contingency at each compute().
        //   _switchable_buses_    : solver bus ids that can flip PV -> PQ in SOME row,
        //                           sorted; each owns a reserved Vm unknown + Q
        //                           equation for the whole sweep.
        //   _row_pv_to_pq_        : per row, the switchable buses that actually flip
        //                           (sorted). Empty vector == that row keeps them all PV.
        //   _row_slack_gens_off_  : per row, the disconnected generators that carried a
        //                           non-zero distributed-slack weight. Usually empty.
        //   _row_pv_pinned_       : per row, the buses to pin PV (sorted): the switchable
        //                           ones still PV plus the base-PQ buses a reactivated
        //                           generator pins (set_topo_actions).
        //   _row_vm_seed_         : per row, (solver bus, |V|) of the latter, seeded before
        //                           the solve.
        //   _pv_pinning_active_   : whether any row pins differently from the resting
        //                           state (every switchable bus pinned, nothing else).
        // the grid's outer loops, for their triggers (see _record_outer_detect)
        std::vector<std::shared_ptr<BaseOuterLoop> > _outer_detect_loops_;
        bool _gen_contingency_active_ = false;
        std::vector<int> _switchable_buses_;
        std::vector<std::vector<int> > _row_pv_to_pq_;
        std::vector<std::vector<int> > _row_pv_pinned_;
        std::vector<std::vector<std::pair<int, real_type> > > _row_vm_seed_;
        bool _pv_pinning_active_ = false;
        std::vector<std::vector<int> > _row_slack_gens_off_;
        // per row, the storage units taking part in the slack its topological action
        // disconnects (_maybe_resolve_topology), for _row_slack_weights as well
        std::vector<std::vector<int> > _row_slack_storages_off_;
        // topological actions (ScenarioSweep only; plain, always-present state like
        // the above). topo_actions_ is a registration (one checked TopoAction per
        // row, see set_topo_actions); _row_topo_ is each action resolved against the
        // grid (L1, _maybe_topo_union) and _topology_active_ what
        // _maybe_resolve_topology made of them for this batch (L2).
        std::vector<TopoAction> topo_actions_;
        std::vector<RowTopoPlan> _row_topo_;
        bool _topology_active_ = false;
        // whether the policies hold the resolution of the registered actions (see
        // get_row_disconnected_branches): set at the end of _maybe_resolve_topology,
        // dropped with the batch inputs
        bool _topo_resolved_ = false;
        // L2, per row: the base buses its action leaves with no element (see
        // _emptied_buses), masked without stranding anything
        std::vector<std::vector<int> > _row_emptied_buses_;
        // L2, per row: the buses its action puts something on (see _placed_buses)
        std::vector<std::vector<int> > _row_placed_buses_;
        // L1: the buses the actions use (gridmodel ids, sorted), those of them empty in
        // the base grid, and the (bus, bus) pairs a row's branch placement writes
        std::vector<int> _topo_used_buses_me_;
        std::vector<int> _extra_buses_me_;
        std::vector<std::pair<int, int> > _union_edges_me_;
        // L2: the extra buses in solver numbering (sorted): masked in the "n" case and in
        // every row that does not use them -- the masked loop's resting state
        std::vector<int> _base_masked_;
        // OLF-style bounded slack redistribution (see _prepare_slack_redistribution):
        // the option, and the per-row data it builds at each compute()
        bool _redistribute_slack_ = false;
        std::vector<std::vector<std::pair<int, real_type> > > _row_slack_dp_pu_;

        // per-row solve statistics of the last compute(): see get_row_nb_iter
        std::vector<double> _row_solve_time_;
        std::vector<int> _row_nb_iter_;
        std::vector<LinearSolverStats> _row_lin_stats_;
        std::vector<int> _row_J_dim_;
        std::vector<int> _row_J_nnz_;
        // what the outer loops did on each row, when the algorithm runs them (empty otherwise)
        std::vector<OuterRowResult> _outer_results_;
        std::vector<std::vector<std::pair<int, real_type> > > _row_gen_new_p_;
        std::vector<std::vector<std::pair<int, real_type> > > _row_sto_new_p_;
        std::vector<std::vector<int> > _row_sat_gens_;
        std::vector<std::vector<int> > _row_sat_storages_;
        std::vector<real_type> _row_lost_p_mw_;  // outer loops: see _clear_slack_redistribution
        // per contingency, the solver buses it strands (sorted, empty if none) and
        // whether the grid stays connected -- both settled by _prepare_connectivity
        // from bus_graph_, the DFS tree of the base graph built there.
        std::vector<std::vector<int> > _li_masked;
        std::vector<char> _cont_connected_;
        BusGraph bus_graph_;
        std::vector<char> _skip_mask;
        bool _compute_limit_violations_ = false;
        real_type _violation_threshold_ = 1.0;
        real_type _violation_rel_tol_ = DEFAULT_VIOLATION_REL_TOL;
        std::vector<char> _converged;
        // twin of _converged above, but unconditional (every instantiation, not just
        // ContingencyAnalysis; not gated behind compute_limit_violations) -- see
        // converged_mask().
        std::vector<char> _converged_mask_;
        std::vector<std::vector<LimitViolation> > _violations;
        bool _converged_n_ = false;
        std::vector<LimitViolation> _violations_n_;
        // physical-limit checks (every instantiation, see
        // set_compute_physical_violations). The plans are the routing the per-row checks
        // need and are rebuilt once per compute(), from the labelling the batch solves in
        // -- like every other solver-keyed thing here they mean nothing against another
        // labelling, so they are dropped with the results. The outer loops' own triggers need
        // no plan: each row's are read through its OuterRowView (see _record_outer_detect).
        bool _compute_physical_violations_ = false;
        // see set_compute_gen_results
        bool _compute_gen_results_ = false;
        RealMat _gen_results_;
        real_type _physical_tol_mva_ = 1e-4;
        real_type _physical_tol_vm_pu_ = 1e-4;
        // `_gen_pv_release_check_on_`: whether the algorithm of THIS compute() can feed the
        // reactive and voltage checks -- the loops' own and the PQ -> PV release one -- at
        // all (see _prepare_physical_check): false in DC
        bool _gen_pv_release_check_on_ = false;
        gen_pv_release_check::GenPvReleasePlan _gen_pv_release_plan_;
        gen_p_check::GenPPlan _gen_p_plan_;
        std::vector<std::vector<LimitViolation> > _physical_violations_;
        std::vector<LimitViolation> _physical_violations_n_;
        // per-row branch-id cache (ContingencyAnalysis only), refreshed once per
        // compute() call -- avoids rebuilding my_defaults_vect() per row / per thread.
        std::vector<std::vector<int> > _li_defaults_vect_cache_;

        double _timer_modif_Ybus = 0.;

        // base-case reuse (see set_reuse_base_case): whether L1 / L2 may be kept at
        // all, and whether the last compute() actually kept them.
        bool _reuse_base_case_ = true;
        bool _base_case_was_reused_ = false;

        // The worker algorithms of a multi-threaded compute(), kept between calls --
        // part of L2, same as the member _algo (see clear_batch_inputs), and kept on
        // the same terms (set_reuse_base_case governs both: there is no separate switch,
        // because there is no state in which keeping one and rebuilding the other would
        // be the right answer). Empty on the single-threaded path, and whenever the
        // batch inputs were dropped. Indexed by thread, and every call gives thread t
        // the same one back, so nothing is shared between threads.
        std::vector<std::unique_ptr<AlgorithmSelector> > _thread_algos_;
        bool _thread_algos_were_reused_ = false;

        // What a row's gen_v has to agree with, at the buses where it is not free:
        // `gens` are the generators writing one bus, and `fixed_vm` the value something
        // ELSE on that bus asks for -- an SVC or an hvdc converter station, whose
        // set-point modify_gen_v cannot move. Rebuilt by _prepare_gen_v_constraints at
        // each compute(), and empty unless modify_gen_v was used (the grid's own
        // target_vm_pu_ is LSGrid's to validate, not this class's). One entry only for a
        // bus that can actually be contradicted. See _row_gen_v_conflicts.
        struct GenVConstraint
        {
            int bus_solver;   // the bus they all write
            std::vector<int> gens;
            bool has_fixed;
            real_type fixed_vm;
        };
        std::vector<GenVConstraint> _gen_v_constraints_;
        // (generator, voltage-control group) whose v_set its gen_v is -- see
        // _apply_step_vc_v_set. Rebuilt by _prepare_gen_v_vc at each compute().
        std::vector<std::pair<int, int> > _gen_v_vc_pairs_;
        Eigen::Index _gen_v_vc_n_groups_ = 0;

        // reverse-mode differentiation (see the public block above and BatchAdjoint)
        bool _keep_jacobian_ = false;
        BatchAdjoint _adjoint_;
        std::vector<char> _adjoint_row_ok_;
        // per row, what its algorithm said of the system it ended on (filled with its
        // Jacobian), and the gen_v layout the columns were asked for
        std::vector<AdjointSolveInfo> _adjoint_info_;
        power_flow_adjoint::GenVLayout _adjoint_gen_v_;
};

/**
Time series: the same grid topology is used along with time series of injections
(productions and loads) to compute powerflows. Each step is warm started with the
solution of the previous one, so the steps are solved in order, on a single thread.
 **/
using TimeSeries = BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromPreviousStep>;

/**
Injection sweep: same inputs, same outputs, same interface as TimeSeries, but every step
restarts from the same seed -- independent scenarios rather than consecutive instants,
splittable over several threads.
 **/
using InjectionSweep = BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromSeed>;

/**
Contingency (security) analysis: the same base injection is used along with a set of
distinct topology changes (line/trafo disconnections) to compute powerflows.
 **/
using ContingencyAnalysis = BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::NOOP, BatchInitKind::FromSeed>;

/**
Scenario sweep: varies both the injection AND a contingency per step, independently,
row-aligned (row i pairs its own injection with its own contingency mask). The 4th
instantiation of this template -- see set_contingency_lines/set_contingency_trafos.
 **/
using ScenarioSweep = BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::Vary, BatchInitKind::FromSeed>;

#ifndef LS2G_BUILDING_CORE
    // all 4 instantiations are compiled once, into the core library (see the bottom
    // of BaseBatchSweep.cpp); every other translation unit links against them
    // instead of instantiating the template again (this only covers the class's
    // NON-template members -- compute()/_run_range()/_run_one_step() -- everything
    // else here is a member template, instantiated independently wherever it is
    // actually called, per ordinary template rules).
    extern template class LS2G_API BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromPreviousStep>;
    extern template class LS2G_API BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromSeed>;
    extern template class LS2G_API BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::NOOP, BatchInitKind::FromSeed>;
    extern template class LS2G_API BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::Vary, BatchInitKind::FromSeed>;
#endif  // LS2G_BUILDING_CORE

} // namespace ls2g

#endif  // BASEBATCHSWEEP_H
