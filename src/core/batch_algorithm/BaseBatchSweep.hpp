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
#include "BusQCheck.hpp"
#include "GenPCheck.hpp"
#include "HvdcPCheck.hpp"
#include "BusGraph.hpp"
#include "BatchAdjoint.hpp"

#include <set>
#include <map>
#include <vector>
#include <memory>
#include <queue>
#include <algorithm>
#include <iterator>
#include <exception>
#include <sstream>
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

namespace batch_sweep_detail {

// appends BUS limit violations for a single, already-solved voltage vector `V` (solver
// numbering). `masked_solver_ids` is nullptr for the base ("n") case; in "handle
// disconnected grid" mode it is this contingency's masked solver-bus-id list -- those
// buses are forced to V=0 post-solve and must never be checked. `subs` resolves the
// violating bus's substation name. `threshold` (in ]0., 1.]) tightens both checks; 1.
// is the "report exactly at the configured limit" default. Moved verbatim (as a free
// function, `inline` so it can live in this header without an ODR violation) from the
// pre-refactor ContingencyAnalysis.cpp's anonymous namespace.
inline void check_bus_voltage_violations(
    const Eigen::Ref<const CplxVect> & V,
    const SolverBusIdVect & id_me_to_solver,
    const Eigen::Ref<const RealVect> & bus_vmin_kv,
    const Eigen::Ref<const RealVect> & bus_vmax_kv,
    const Eigen::Ref<const RealVect> & bus_vn_kv,
    const SubstationContainer & subs,
    real_type threshold,
    const std::vector<int> * masked_solver_ids,
    std::vector<LimitViolation> & out)
{
    const Eigen::Index nb_bus = bus_vmin_kv.size();
    if(nb_bus == 0) return;  // voltage limits never configured on this grid

    const std::vector<std::string> & sub_names = subs.get_sub_names();  // empty if never set

    std::vector<bool> is_masked;
    if(masked_solver_ids != nullptr && !masked_solver_ids->empty()){
        is_masked.assign(static_cast<size_t>(V.size()), false);
        for(int b : *masked_solver_ids){
            if(b >= 0 && static_cast<size_t>(b) < is_masked.size()) is_masked[static_cast<size_t>(b)] = true;
        }
    }

    for(Eigen::Index grid_id = 0; grid_id < nb_bus; ++grid_id){
        const real_type vmin = bus_vmin_kv(grid_id);
        const real_type vmax = bus_vmax_kv(grid_id);
        if(isnan(vmin) && isnan(vmax)) continue;  // no limit configured for this specific bus

        const int solver_id = id_me_to_solver(static_cast<int>(grid_id)).cast_int();
        if(solver_id == BaseConstants::_deactivated_bus_id) continue;
        if(!is_masked.empty() && is_masked[static_cast<size_t>(solver_id)]) continue;

        const real_type vm_kv = std::abs(V(solver_id)) * bus_vn_kv(grid_id);
        const bool has_min = !isnan(vmin);
        const bool has_max = !isnan(vmax);

        const real_type vnom = bus_vn_kv(grid_id);
        real_type anchor = vnom;
        if(has_min && (anchor < vmin)) anchor = vmin;
        if(has_max && (anchor > vmax)) anchor = vmax;
        const real_type shrink = 1. - threshold;
        const real_type low_eff = vmin + shrink * (anchor - vmin);
        const real_type high_eff = vmax - shrink * (vmax - anchor);

        if(has_min && has_max && low_eff > high_eff){
            std::ostringstream exc_;
            exc_ << "ContingencyAnalysis: bus " << grid_id << " has inconsistent voltage "
                    "limits: its effective minimum (" << low_eff << " kV, from vmin = "
                 << vmin << " kV) is above its effective maximum (" << high_eff
                 << " kV, from vmax = " << vmax << " kV) at violation_threshold = "
                 << threshold << ". Check the vmin_kv / vmax_kv passed to "
                    "LSGrid::set_bus_voltage_limits (they must satisfy vmin <= vmax, and "
                    "should bracket the bus nominal voltage of " << vnom << " kV).";
            throw std::runtime_error(exc_.str());
        }

        if(has_min && vm_kv <= low_eff){
            const int sub_id = subs.sub_id_of_bus(static_cast<int>(grid_id));
            const std::string sub_name = static_cast<size_t>(sub_id) < sub_names.size() ? sub_names[sub_id] : std::string();
            out.push_back(LimitViolation{ViolationElementType::BUS, static_cast<int>(grid_id), 0,
                                          LimitViolationType::LOW_VOLTAGE, vm_kv, vmin, sub_name});
        } else if(has_max && vm_kv >= high_eff){
            const int sub_id = subs.sub_id_of_bus(static_cast<int>(grid_id));
            const std::string sub_name = static_cast<size_t>(sub_id) < sub_names.size() ? sub_names[sub_id] : std::string();
            out.push_back(LimitViolation{ViolationElementType::BUS, static_cast<int>(grid_id), 0,
                                          LimitViolationType::HIGH_VOLTAGE, vm_kv, vmax, sub_name});
        }
    }
}

// appends CURRENT limit violations (both sides) for a single, already-solved voltage
// vector `V` (solver numbering), for every element of `structure_data` (LineContainer
// or TrafoContainer). `skip_ids` are the LOCAL (own-type) element ids disconnected BY
// THIS CONTINGENCY. Moved verbatim from the pre-refactor
// ContingencyAnalysis.cpp's anonymous namespace.
template<class T>
inline void check_current_violations(
    const T & structure_data,
    ViolationElementType el_type,
    const Eigen::Ref<const CplxVect> & V,
    const SolverBusIdVect & id_me_to_solver,
    const Eigen::Ref<const RealVect> & bus_vn_kv,
    bool ac_solver_used,
    real_type sn_mva,
    const Eigen::Ref<const RealVect> & limit1,
    const Eigen::Ref<const RealVect> & limit2,
    real_type threshold,
    const std::vector<int> & skip_ids,
    std::vector<LimitViolation> & out,
    // the row's own placement of some branches (see BaseBatchSolverSynch::
    // _row_branch_overrides_), sorted by branch id, and the offset of this container's
    // ids in that numbering (0 for lines, n_line for trafos); nullptr for none
    const std::vector<BranchBusOverride> * overrides = nullptr,
    size_t lag_id = 0)
{
    if(limit1.size() == 0 && limit2.size() == 0) return;  // thermal limits never configured

    const auto & el_status = structure_data.get_status_global();
    const auto & status1 = structure_data.get_status_side_1();
    const auto & status2 = structure_data.get_status_side_2();
    const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
    const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
    const std::vector<std::string> & el_names = structure_data.get_names();  // empty if never set

    Eigen::Ref<const CplxVect> yac_eff_11 = structure_data.yac_eff_11();
    Eigen::Ref<const CplxVect> yac_eff_12 = structure_data.yac_eff_12();
    Eigen::Ref<const CplxVect> yac_eff_21 = structure_data.yac_eff_21();
    Eigen::Ref<const CplxVect> yac_eff_22 = structure_data.yac_eff_22();
    Eigen::Ref<const CplxVect> yac_raw_11 = structure_data.yac_11();
    Eigen::Ref<const CplxVect> yac_raw_12 = structure_data.yac_12();
    Eigen::Ref<const CplxVect> yac_raw_21 = structure_data.yac_21();
    Eigen::Ref<const CplxVect> yac_raw_22 = structure_data.yac_22();
    Eigen::Ref<const RealVect> ydc_11 = structure_data.ydc_11();
    Eigen::Ref<const RealVect> ydc_12 = structure_data.ydc_12();
    Eigen::Ref<const RealVect> ydc_21 = structure_data.ydc_21();
    Eigen::Ref<const RealVect> ydc_22 = structure_data.ydc_22();
    Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift();
    const bool has_tau_shift = dc_x_tau_shift.size() > 0;  // trafo only, empty for lines

    const size_t nb_el = structure_data.nb();
    const real_type sqrt_3 = sqrt(3.);

    std::vector<bool> skip;
    if(!skip_ids.empty()){
        skip.assign(nb_el, false);
        for(int id : skip_ids){
            if(id >= 0 && static_cast<size_t>(id) < nb_el) skip[static_cast<size_t>(id)] = true;
        }
    }

    size_t next_override = 0;
    if(overrides != nullptr){
        while(next_override < overrides->size() && (*overrides)[next_override].branch_id < static_cast<int>(lag_id)) ++next_override;
    }
    for(size_t el_id = 0; el_id < nb_el; ++el_id){
        bool on = el_status[el_id];
        bool s1 = status1[el_id];
        bool s2 = status2[el_id];
        int bus_from_me = BaseConstants::_deactivated_bus_id;
        int bus_to_me = BaseConstants::_deactivated_bus_id;
        bool own_placement = false;
        if(overrides != nullptr && next_override < overrides->size() &&
           (*overrides)[next_override].branch_id == static_cast<int>(el_id + lag_id)){
            const BranchBusOverride & ov = (*overrides)[next_override];
            ++next_override;
            on = ov.connected;
            if(on){ s1 = true; s2 = true; bus_from_me = ov.from_me; bus_to_me = ov.to_me; own_placement = true; }
        }
        if(!on) continue;
        if(!skip.empty() && skip[el_id]) continue;

        const Eigen::Index el_idx = static_cast<Eigen::Index>(el_id);
        const bool has_lim1 = limit1.size() > 0 && !isnan(limit1(el_idx));
        const bool has_lim2 = limit2.size() > 0 && !isnan(limit2(el_idx));
        if(!has_lim1 && !has_lim2) continue;

        int solver_from = BaseConstants::_deactivated_bus_id;
        int solver_to = BaseConstants::_deactivated_bus_id;
        if(s1){
            if(!own_placement) bus_from_me = bus_from(static_cast<int>(el_id)).cast_int();
            solver_from = id_me_to_solver(bus_from_me).cast_int();
            if(solver_from == BaseConstants::_deactivated_bus_id) continue;
        }
        if(s2){
            if(!own_placement) bus_to_me = bus_to(static_cast<int>(el_id)).cast_int();
            solver_to = id_me_to_solver(bus_to_me).cast_int();
            if(solver_to == BaseConstants::_deactivated_bus_id) continue;
        }

        const cplx_type Efrom = s1 ? V(solver_from) : cplx_type(0., 0.);
        const cplx_type Eto = s2 ? V(solver_to) : cplx_type(0., 0.);
        const real_type v_from_kv = s1 ? std::abs(Efrom) * bus_vn_kv(bus_from_me) : real_type(1.);
        const real_type v_to_kv = s2 ? std::abs(Eto) * bus_vn_kv(bus_to_me) : real_type(1.);

        real_type amps1 = 0.;
        real_type amps2 = 0.;
        if(ac_solver_used){
            // a branch the row placed itself has both ends closed: its raw block
            const cplx_type y11 = own_placement ? yac_raw_11(el_idx) : yac_eff_11(el_idx);
            const cplx_type y12 = own_placement ? yac_raw_12(el_idx) : yac_eff_12(el_idx);
            const cplx_type y21 = own_placement ? yac_raw_21(el_idx) : yac_eff_21(el_idx);
            const cplx_type y22 = own_placement ? yac_raw_22(el_idx) : yac_eff_22(el_idx);
            const cplx_type I_from = std::conj(y11 * Efrom + y12 * Eto);
            const cplx_type S_from = Efrom * I_from;
            amps1 = std::abs(S_from) * sn_mva / (sqrt_3 * v_from_kv);

            const cplx_type I_to = std::conj(y22 * Eto + y21 * Efrom);
            const cplx_type S_to = Eto * I_to;
            amps2 = std::abs(S_to) * sn_mva / (sqrt_3 * v_to_kv);
        } else if(s1 && s2){
            const real_type theta_from = std::arg(Efrom);
            const real_type theta_to = std::arg(Eto);
            real_type p_from = (ydc_11(el_idx) * theta_from + ydc_12(el_idx) * theta_to) * sn_mva;
            if(has_tau_shift) p_from -= dc_x_tau_shift(el_idx);
            amps1 = std::abs(p_from) / (sqrt_3 * v_from_kv);

            real_type p_to = (ydc_22(el_idx) * theta_to + ydc_21(el_idx) * theta_from) * sn_mva;
            if(has_tau_shift) p_to += dc_x_tau_shift(el_idx);
            amps2 = std::abs(p_to) / (sqrt_3 * v_to_kv);
        }

        // the name is copied only into a violation: one std::string per element
        // examined was a malloc per branch per row on a grid with long names
        if(has_lim1 && amps1 >= threshold * limit1(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 1,
                                          LimitViolationType::CURRENT, amps1, limit1(el_idx),
                                          el_id < el_names.size() ? el_names[el_id] : std::string()});
        }
        if(has_lim2 && amps2 >= threshold * limit2(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 2,
                                          LimitViolationType::CURRENT, amps2, limit2(el_idx),
                                          el_id < el_names.size() ? el_names[el_id] : std::string()});
        }
    }
}

}  // namespace batch_sweep_detail

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
                // the union layout the topological actions asked for (see _maybe_topo_union)
                _topo_used_buses_me_.clear();
                _extra_buses_me_.clear();
                _union_edges_me_.clear();
            }
            BaseBatchSolverSynch::clear_grid_results();
        }

        // L2 -- what was built for this batch from that grid: what the graph walk
        // settled about each contingency, the reference-slack choice and the skip
        // mask that came with it, and the PV <-> PQ structure a generator contingency
        // reserved. The base drops the algorithm, whose sparsity that structure is.
        void clear_batch_inputs() override {
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
                _li_defaults_vect_cache_.clear();
                _physical_violations_n_.clear();
                // what the topological actions were resolved into for this batch
                // (the actions themselves are a registration and stay)
                _row_topo_.clear();
                _topology_active_ = false;
                _base_masked_.clear();
                _reset_topo_policy_state();
            }
            BaseBatchSolverSynch::clear_batch_inputs();
        }

        // L3 -- what the last run produced. On top of the voltages and the flows the
        // base drops: convergence, limit violations, and the Jacobians the adjoint
        // kept.
        void clear_batch_outputs() override {
            _adjoint_.clear();
            _adjoint_row_ok_.clear();
            _converged.clear();
            _converged_mask_.clear();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _physical_violations_.clear();
            _bus_q_plan_.clear();
            _hvdc_p_plan_.clear();
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
         * ScenarioSweep only: one topological action per row (a `TopoAction`, grid2op
         * semantics: set_bus on loads / generators / storage units / line and trafo
         * ends, set_line_status), played by that row on top of its injections and
         * masks. An action that changes nothing is a plain row.
         *
         * Every action is checked against the grid here (the element and the busbar
         * exist, no contradiction -- see TopoAction::check_validity) and a
         * `std::invalid_argument` naming the row is raised, nothing registered, if one
         * is invalid. What THIS version can play is settled at compute() (see
         * _maybe_resolve_topology): disconnections only -- a branch (set_line_status
         * -1, or set_bus -1 on one of its ends), a generator, a load or a storage unit
         * (set_bus -1). A branch or a generator is disconnected exactly as the
         * set_contingency_* masks do it -- same Ybus edit, same PV -> PQ handling, one
         * symbolic analysis for the whole sweep -- and a row naming an element in both
         * a mask and its action is refused. Moving an element to a busbar, reconnecting
         * one, and the DC algorithm are refused for now (TODO in the changelog).
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
            if(!topo_actions_.empty() && !_batch_inputs_valid_){
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
                std::set<int> this_default = {l_id};
                ybus_policy_.li_defaults.insert(this_default);
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_n1(int line_id){
            _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            std::set<int> this_default = {line_id};
            ybus_policy_.li_defaults.insert(this_default);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s) _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            for(const auto line_id : vect_n1s){
                std::set<int> this_default = {line_id};
                ybus_policy_.li_defaults.insert(this_default);
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk){
                _check_ok_el(line_id);
                this_default.insert(line_id);
            }
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            ybus_policy_.li_defaults.insert(this_default);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_n1(int line_id){
            _check_ok_el(line_id);
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            std::set<int> this_default = {line_id};
            auto nb_removed = ybus_policy_.li_defaults.erase(this_default);
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
                std::set<int> this_default = {line_id};
                nb_removed += ybus_policy_.li_defaults.erase(this_default);
            }
            return nb_removed;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk){
                _check_ok_el(line_id);
                this_default.insert(line_id);
            }
            // a contingency registration (L2): what the graph walk settled about the
            // OLD set, the masks and the reserved pv/pq structure all go with it
            clear_batch_inputs();
            auto nb_removed = ybus_policy_.li_defaults.erase(this_default);
            return nb_removed >= 1;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::set<std::set<int> > & my_defaults() const {return ybus_policy_.li_defaults;}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        std::vector<std::vector<int> > my_defaults_vect() const {
            std::vector<std::vector<int> > res;
            res.reserve(ybus_policy_.li_defaults.size());
            for(const auto & set_int: ybus_policy_.li_defaults){
                std::vector<int> this_def(set_int.begin(), set_int.end());
                res.push_back(this_def);
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
        // grid can reach at all -- ViolationCategory::PHYSICAL, as opposed to the
        // operational limits `compute_limit_violations` reports (a voltage band, a thermal
        // rating: states the grid does reach and should not sit in). Three today, each a
        // condition an OpenLoadFlow outer loop acts on, and none enforced here -- no bus is
        // switched PV -> PQ, no droop is clamped, no machine leaves the slack distribution,
        // no row is re-solved:
        //
        //   * the REACTIVE CAPABILITY of each bus whose voltage is held by machines
        //     (LOW_Q / HIGH_Q on the BUS, see BusQCheck.hpp): did it need more reactive
        //     power than the sum of what its voltage-regulating generators, hvdc converter
        //     stations and voltage-mode SVCs can produce? OpenLoadFlow's `ReactiveLimits`;
        //   * the ACTIVE POWER of each angle-droop ("AC emulation") hvdc line still in the
        //     linear regime (HIGH_P on the HVDC, see HvdcPCheck.hpp): did it transmit more
        //     than `pmax_1to2_mw` / `pmax_2to1_mw` allow in that direction? OpenLoadFlow's
        //     `HvdcAcEmulationLimits`;
        //   * the ACTIVE POWER of each generator carrying the DISTRIBUTED SLACK (LOW_P /
        //     HIGH_P on the GENERATOR, see GenPCheck.hpp): the slack is solved inside the
        //     Jacobian by participation factors that know nothing about limits, so
        //     `target_p + share` can land beyond `min_p_mw` / `max_p_mw`. OpenLoadFlow's
        //     `DistributedSlack`, whose whole job is to take a saturated machine out of the
        //     distribution and re-share what is left.
        //
        // Available on all four instantiations (unlike compute_limit_violations, which is
        // contingency-only): a time series that walks a load curve is exactly as likely to
        // ask a bus for reactive power it does not have, or a machine for active power it
        // cannot deliver, as a contingency is.
        //
        // WHAT EACH CHECK NEEDS. The hvdc one needs only the bus angles and the generator
        // one only the slack the row distributed -- both of which every algorithm leaves
        // behind, AC and DC alike. The reactive one needs an AC algorithm that publishes
        // its per-bus mismatch (BaseAlgo::fills_bus_mismatch -- every built-in AC family
        // does; a plugin solver has to opt in), and compute() raises for one that does not
        // rather than reporting nothing. In DC it is simply not applicable: a DC powerflow
        // has no reactive power at all, so a DC batch reports the two active-power checks
        // and nothing is hidden by it. The generator check also needs the limits
        // themselves, which are optional (LSGrid::set_gen_p_limits): a grid that has none
        // simply reports nothing there.
        //
        // Setting this drops this batch's base case and results, but NOT its registrations
        // (the contingencies / injections), so unlike compute_limit_violations -- which
        // clear()s the whole object -- it can be set at any point before compute().
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
        /**
         * Per row: the physical limits this row's solution leaves. A row that did not
         * converge (or that was never simulated) has an EMPTY entry rather than a sentinel
         * -- ask converged_mask() to tell that apart from "converged, no violation". Every
         * entry has category PHYSICAL, and one of two shapes:
         *
         *   - element_type BUS, element_id the grid bus id, violation_type LOW_Q / HIGH_Q,
         *     `value` the reactive power the machines holding that bus had to produce
         *     (MVAr) and `limit` their summed capability;
         *   - element_type HVDC, element_id the hvdc line id, violation_type HIGH_P, `side`
         *     the direction (1 for 1 -> 2), `value` the active power leaving that side (MW,
         *     positive) and `limit` that direction's pmax;
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
            _select_ref_slack_and_masks();
            if(active_layout().slack_bus_id_me.size() == 0) return -1;
            return active_layout().slack_bus_id_me[static_cast<int>(0)].cast_int();
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
        // not regulating or treated as off; for one regulating a bus the solve gives a
        // magnitude unknown to (a bus a voltage-control group holds: the set-point is
        // then that group's v_set, not a fixed |V| -- see get_gen_v_vc_row); and for one
        // whose bus is pinned by an SVC or an hvdc converter station, whose set-point no
        // batch can move, so the generator's is not free either.
        IntVect get_gen_v_target_bus() const {
            std::map<int, std::vector<int> > gens_of_bus;
            _grid_model.get_generators().vm_writers_by_bus(active_layout().id_me_to_solver, gens_of_bus);
            std::map<int, real_type> pinned;
            _grid_model.get_svcs().vm_targets_by_bus(active_layout().id_me_to_solver, pinned);
            _grid_model.get_dclines().vm_targets_by_bus(active_layout().id_me_to_solver, pinned);

            const IntVect vm_col = _algo.get_vm_to_J_col_python();
            const auto solver_to_me = active_layout().id_solver_to_me.as_eigen();
            IntVect res = IntVect::Constant(
                static_cast<Eigen::Index>(_grid_model.get_generators_as_data().nb()), -1);
            for(const auto & bus_and_gens : gens_of_bus){
                const int b = bus_and_gens.first;
                if(b < 0 || b >= solver_to_me.size()) continue;
                if(b < vm_col.size() && vm_col[b] >= 0) continue;      // |V| is an unknown there
                if(pinned.find(b) != pinned.end()) continue;           // held by something else
                for(size_t k = 0; k < bus_and_gens.second.size(); ++k){
                    res[static_cast<Eigen::Index>(bus_and_gens.second[k])] = solver_to_me[b];
                }
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
        RealVect get_gen_v_share() const {
            const IntVect target = get_gen_v_target_bus();
            const IntVect vc_group = _gen_v_vc_group();
            std::map<int, int> count, count_vc;   // per fixed bus / per group
            for(Eigen::Index g = 0; g < target.size(); ++g){
                if(target[g] >= 0) count[target[g]] += 1;
                else if(g < vc_group.size() && vc_group[g] >= 0) count_vc[vc_group[g]] += 1;
            }
            RealVect res = RealVect::Zero(target.size());
            for(Eigen::Index g = 0; g < target.size(); ++g){
                if(target[g] >= 0) res[g] = 1. / static_cast<real_type>(count[target[g]]);
                else if(g < vc_group.size() && vc_group[g] >= 0)
                    res[g] = 1. / static_cast<real_type>(count_vc[vc_group[g]]);
            }
            return res;
        }

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
            const IntVect vc_group = _gen_v_vc_group();
            const IntVect v_row = _algo.get_group_v_row();
            IntVect res = IntVect::Constant(vc_group.size(), -1);
            for(Eigen::Index g = 0; g < vc_group.size(); ++g){
                const int grp = vc_group[g];
                if(grp >= 0 && grp < v_row.size()) res[g] = v_row[grp];
            }
            return res;
        }

        // The indirect half of the gen_v gradient: `-lambda^T dF/dv` per row and per
        // generator -- for a generator of a voltage-control group, all of it: lambda at
        // the group's voltage row (get_gen_v_vc_row), zero on a row where
        // handle_disconnected_grid stranded the group -- otherwise
        // keyed like get_gen_v_target_bus() (zero wherever that says -1, and
        // on a row that did not converge) and already carrying each generator's share of
        // its bus (see get_gen_v_share -- the caller must weight the DIRECT half by the
        // same thing). `lambda` is what solve_JT returned, so the
        // adjoint system is solved once and both halves of the gradient read it.
        //
        // dF/dv_k is the dS/d|V_k| column, taken on THIS row's admittance matrix -- the
        // base one with the row's own contingency edits applied, which is why this
        // lives here and not in the caller.
        BatchAdjoint::RealMatRM gen_v_indirect_grad(const Eigen::Ref<const BatchAdjoint::RealMatRM> & lambda);

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
            _grid_model.get_generators().set_vm(res, me_to_solver);
            _grid_model.get_dclines().set_vm(res, me_to_solver);
            _grid_model.get_svcs().set_vm(res, me_to_solver);
            return res;
        }

        // ----- reverse-mode differentiation helpers -------------------------------

        // Re-key a solver-bus-keyed map onto grid bus ids (see get_p_row_of_bus).
        IntVect _bus_map_to_grid(const IntVect & solver_keyed) const {
            const auto me_to_solver = active_layout().id_me_to_solver.as_eigen();
            IntVect res = IntVect::Constant(me_to_solver.size(), -1);
            for(Eigen::Index bus_me = 0; bus_me < me_to_solver.size(); ++bus_me){
                const int bus_solver = me_to_solver[bus_me];
                if(bus_solver < 0) continue;   // bus not in the solver (deactivated, isolated)
                if(bus_solver >= solver_keyed.size()) continue;
                res[bus_me] = solver_keyed[bus_solver];
            }
            return res;
        }

        // Keep this row's converged Jacobian for the adjoint. MUST be called while the
        // row's own state is still installed on the algorithm -- its Ybus edits, its
        // masked buses, its PV pinning -- because refresh_J_at_solution() re-fills J
        // from exactly that, and a J refreshed after the loop restored its resting
        // state would describe a system this row never solved.
        void _maybe_store_jacobian(size_t i, AlgorithmSelector & algo) {
            if(!_adjoint_.is_allocated()) return;
            algo.refresh_J_at_solution();
            _adjoint_.store_row(static_cast<Eigen::Index>(i), algo.get_J());
        }

        // Which equations each row froze to the identity, as Jacobian row indices --
        // what solve_JT drops from lambda (see BatchAdjoint::solve_JT). Two sources,
        // and a row can have both:
        //   - a bus a contingency stranded: masked, so BOTH its P and its Q equation;
        //   - a bus still PV in this row (a reserved switchable bus whose generator
        //     is still on): pinned, so its Q equation only.
        // Empty when the batch does neither, which is the common case.
        std::vector<std::vector<int> > _adjoint_identity_rows() const {
            const Eigen::Index nb_rows = _adjoint_.nb_rows();
            // _li_masked is filled whenever connectivity was analysed, but only the
            // masked mode ever hands it to the algorithm: without it a contingency
            // that strands a bus makes the row skipped entirely, not masked
            const bool has_masking = _mask_mode() && !_li_masked.empty();
            const bool has_pinning = _has_pv_switching();
            if(nb_rows <= 0 || (!has_masking && !has_pinning)) return std::vector<std::vector<int> >();

            const IntVect p_row = _algo.get_p_to_J_row_python();
            const IntVect q_row = _algo.get_q_to_J_row_python();
            auto push = [](const IntVect & map, int bus, std::vector<int> & out){
                if(bus >= 0 && bus < map.size() && map[bus] >= 0) out.push_back(map[bus]);
            };

            std::vector<std::vector<int> > res(static_cast<size_t>(nb_rows));
            for(Eigen::Index i = 0; i < nb_rows; ++i){
                std::vector<int> & row = res[static_cast<size_t>(i)];
                if(has_masking && static_cast<size_t>(i) < _li_masked.size()){
                    for(int bus : _li_masked[static_cast<size_t>(i)]){
                        push(p_row, bus, row);
                        push(q_row, bus, row);
                    }
                }
                if(has_pinning){
                    // _row_pv_pinned falls back to the whole switchable set for a row
                    // that flips nothing -- which is right: such a row leaves every
                    // reserved bus pinned, exactly as the loop's resting state has it.
                    for(int bus : _row_pv_pinned(static_cast<size_t>(i))) push(q_row, bus, row);
                }
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
            if(ac_solver_used) bus_graph_.build(ac_cache_.mat, threshold);
            else bus_graph_.build(dc_cache_.mat, threshold);

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
                if(cut.verdict == BusGraph::Verdict::Connected) continue;
                if(cut.verdict == BusGraph::Verdict::Split){
                    bus_graph_.stranded_buses(cut, _li_masked[cont_id]);
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
                // connected: nothing stranded beyond the buses masked in every row (the
                // extra buses of the union layout a row does not use, see _base_masked_)
                _cont_connected_[cont_id] = _only_base_masked(_li_masked[cont_id]) ? 1 : 0;
            }
        }

        // whether `masked` (sorted) holds nothing but buses of _base_masked_ (sorted)
        bool _only_base_masked(const std::vector<int> & masked) const {
            if(masked.empty()) return true;
            return std::includes(_base_masked_.begin(), _base_masked_.end(), masked.begin(), masked.end());
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _prepare_connectivity(){}

        // pre-pass run before the (n-)powerflow in "handle disconnected grid" mode:
        // from _li_masked (see _prepare_connectivity, which must have run), chooses
        // the reference slack that minimises the number of skipped contingencies
        // (reordering active_layout().slack_bus_id_me so it is index 0) and fills
        // _skip_mask. Called from _maybe_prepare_masks() (compute()-only, both
        // instantiations) and pick_reference_slack() (public, ContingencyAnalysis-only,
        // hence this too must be fully inline).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _select_ref_slack_and_masks(){
            const size_t nb_cont = ybus_policy_.li_coeffs.size();
            _skip_mask.assign(nb_cont, 0);

            const Eigen::Index nb_slack = active_layout().slack_bus_id_me.size();
            std::vector<int> candidates;
            for(Eigen::Index i = 0; i < nb_slack; ++i){
                const int b = active_layout().slack_bus_id_me[static_cast<int>(i)].cast_int();
                if(b >= 0 && b < active_layout().slack_weights.size() && active_layout().slack_weights(b) > 0.) candidates.push_back(b);
            }
            if(candidates.empty()){
                for(Eigen::Index i = 0; i < nb_slack; ++i) candidates.push_back(active_layout().slack_bus_id_me[static_cast<int>(i)].cast_int());
            }
            if(candidates.empty()) return;

            // every _li_masked entry is sorted (see _prepare_connectivity)
            auto is_masked = [](const std::vector<int> & masked, int bus){
                return std::binary_search(masked.begin(), masked.end(), bus);
            };
            int best_bus = candidates[0];
            int best_strand = -1;
            real_type best_weight = -1.;
            for(int bus : candidates){
                int strand = 0;
                for(const auto & masked : _li_masked) if(is_masked(masked, bus)) ++strand;
                const real_type weight = (bus < active_layout().slack_weights.size()) ? active_layout().slack_weights(bus) : 0.;
                const bool better = (best_strand < 0) ||
                                    (strand < best_strand) ||
                                    (strand == best_strand && weight > best_weight) ||
                                    (strand == best_strand && weight == best_weight && bus < best_bus);
                if(better){ best_strand = strand; best_weight = weight; best_bus = bus; }
            }

            {
                std::vector<int> reordered;
                reordered.reserve(static_cast<size_t>(nb_slack));
                reordered.push_back(best_bus);
                for(Eigen::Index i = 0; i < nb_slack; ++i){
                    const int b = active_layout().slack_bus_id_me[static_cast<int>(i)].cast_int();
                    if(b != best_bus) reordered.push_back(b);
                }
                active_layout().slack_bus_id_me = GlobalBusIdVect(reordered);
            }

            for(size_t cont_id = 0; cont_id < nb_cont; ++cont_id){
                if(is_masked(_li_masked[cont_id], best_bus)) _skip_mask[cont_id] = 1;
                // the masked loop runs for the topological actions' sake as well: without
                // handle_disconnected_grid a row that strands a bus keeps its legacy
                // answer, NOT_SIMULATED
                if(!_handle_disconnected_grid && !_only_base_masked(_li_masked[cont_id])) _skip_mask[cont_id] = 1;
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
        void _reset_ybus_coeffs() { ybus_policy_.li_coeffs.clear(); }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _reset_ybus_coeffs() {}

        // everything, the registrations included: see clear().
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _reset_ybus_policy() {
            ybus_policy_.li_defaults.clear();
            ybus_policy_.li_coeffs.clear();
            ybus_policy_.line_mask = BoolMat();
            ybus_policy_.trafo_mask = BoolMat();
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

        // ----- per-row Ybus edit, matrix only: 2-way ------------------------------
        // The pair above goes through YbusPolicy, which also tells the ALGORITHM its
        // matrix changed. These touch nothing but the coefficients, for a reader that
        // wants row i's matrix without disturbing a solver (gen_v_indirect_grad).
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _patch_ybus_values(Eigen::SparseMatrix<cplx_type> & Ybus, size_t i, bool undo) const {
            if(i >= ybus_policy_.li_coeffs.size()) return;
            for(const auto & c : ybus_policy_.li_coeffs[i]){
                if(undo) Ybus.coeffRef(c.row_id, c.col_id) += c.value;
                else     Ybus.coeffRef(c.row_id, c.col_id) -= c.value;
            }
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _patch_ybus_values(Eigen::SparseMatrix<cplx_type> &, size_t, bool) const {}

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
            _grid_model.get_generators().set_vm(V, active_layout().id_me_to_solver, target_vm_pu_row);
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
                c.gens = bus_and_gens.second;
                c.has_fixed = has_fixed;
                c.fixed_vm = has_fixed ? fixed->second : 0.;
                _gen_v_constraints_.push_back(c);
            }
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _prepare_gen_v_constraints() {}

        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _row_gen_v_conflicts(size_t i) const {
            if(_gen_v_constraints_.empty()) return false;
            const Eigen::Index row = static_cast<Eigen::Index>(i);
            if(row >= sbus_policy_.gen_v.rows()) return false;
            for(size_t c = 0; c < _gen_v_constraints_.size(); ++c){
                const GenVConstraint & con = _gen_v_constraints_[c];
                // whatever else is on the bus pins the value; otherwise the generators
                // only have to agree among themselves
                const real_type ref = con.has_fixed ? con.fixed_vm
                                                    : sbus_policy_.gen_v(row, con.gens[0]);
                for(size_t k = 0; k < con.gens.size(); ++k){
                    if(std::abs(sbus_policy_.gen_v(row, con.gens[k]) - ref) >
                       BaseConstants::_tol_equal_float) return true;
                }
            }
            return false;
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
        // once per compute() from my_defaults_vect() (the sparse, add_n1-driven
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
                                         _violation_threshold_, masked_ids_use, _violations[i]);
            const std::vector<BranchBusOverride> * overrides =
                (i < _row_branch_overrides_.size() && !_row_branch_overrides_[i].empty()) ? &_row_branch_overrides_[i] : nullptr;
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, skip_lines, _violations[i], overrides, 0);
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, skip_trafos, _violations[i], overrides, n_line_);
        }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _record_row_violations_dispatch(size_t i, const CplxVect & V, const std::vector<int> * masked_ids) {
            _record_row_violations(i, V, masked_ids);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _record_row_violations_dispatch(size_t, const CplxVect &, const std::vector<int> *) {}

        // ----- per-row reactive-capability recording ------------------------------
        // Whether THIS row disconnected that generator (a ScenarioSweep generator
        // contingency); always false where no such mask exists.
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        bool _gen_off_in_row(size_t i, int gen_id) const {
            return sbus_policy_.gen_off_in(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(gen_id));
        }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        bool _gen_off_in_row(size_t, int) const { return false; }

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

        real_type _grid_target_p(int gen_id) const {
            const Eigen::Ref<const RealVect> tgt = _grid_model.get_gen_target_p();
            return (gen_id >= 0 && gen_id < tgt.size()) ? tgt(gen_id) : 0.;
        }

        // The active power imbalance a DC row leaves for the slack machines to make up
        // (MW), the same `-sum(Pbus)` LSGrid::_fill_bus_mismatch_dc shares out -- read off
        // the injection the row was actually solved with (empty means the DC entry point
        // fell back to the member dc_cache_.inj, see compute_one_powerflow).
        real_type _dc_imbalance_mw(const Eigen::Ref<const CplxVect> & sbus_solver) const {
            const real_type sn_mva = _grid_model.get_sn_mva();
            if(sbus_solver.size() > 0) return -sbus_solver.real().sum() * sn_mva;
            return -dc_cache_.inj.sum() * sn_mva;
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

        // MUST be called while the row's own state is still installed on the algorithm
        // -- its Ybus edits, its masked buses, its PV pinning -- for the same reason
        // _maybe_store_jacobian must: what is read here (the per-bus mismatch, the
        // controllers' reactive output, the bus angles) belongs to the system THIS row
        // solved, and the `algo` handed in is the one that solved it (the member one, or
        // this thread's). `V_solver` is this row's converged complex voltage: a
        // voltage-mode SVC's capability is a susceptance range, so what it is worth in MVAr
        // depends on it.
        void _record_row_physical(size_t i, AlgorithmSelector & algo,
                                  const Eigen::Ref<const CplxVect> & V_solver,
                                  const Eigen::Ref<const RealVect> & slack_weights,
                                  const Eigen::Ref<const CplxVect> & sbus_solver){
            if(!_compute_physical_violations_) return;
            if(i >= _physical_violations_.size()) return;
            const std::vector<int> * masked = _row_masked_ids(i);
            if(_bus_q_check_on_ && !_bus_q_plan_.empty()){
                // get_controller_q() returns by value: asked for only where a bus' reactive
                // power actually needs it (see BusQPlan::needs_controller_q), so an ordinary
                // grid of local PV machines pays no per-row allocation.
                const RealVect ctrl_q = _bus_q_plan_.needs_controller_q ? algo.get_controller_q()
                                                                       : RealVect();
                bus_q_check::check_bus_q_violations(
                    _bus_q_plan_, _grid_model, algo.get_bus_mismatch(), V_solver, ctrl_q,
                    _grid_model.get_sn_mva(), _physical_tol_mva_, masked,
                    [this, i](int gen_id){ return this->_gen_off_in_row(i, gen_id); },
                    _physical_violations_[i]);
            }
            if(!_hvdc_p_plan_.empty()){
                hvdc_p_check::check_hvdc_p_violations(_hvdc_p_plan_, algo.get_Va(),
                                                      _physical_tol_mva_, masked,
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
                else slack.dc_imbalance_mw = _dc_imbalance_mw(sbus_solver);
                gen_p_check::check_gen_p_violations(
                    _gen_p_plan_, slack, _physical_tol_mva_, masked,
                    [this, i](int gen_id){ return this->_gen_target_p_in_row(i, gen_id); },
                    [this, i](int gen_id){ return this->_gen_off_in_row(i, gen_id); },
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
            _li_defaults_vect_cache_ = _compute_limit_violations_ ? my_defaults_vect() : std::vector<std::vector<int> >();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _refresh_defaults_vect_cache(){}

        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _record_n_case_violations(const CplxVect & V_n){
            if(!_compute_limit_violations_) return;
            _converged_n_ = true;
            const std::vector<int> no_skip;
            batch_sweep_detail::check_bus_voltage_violations(V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vmin_kv(), _grid_model.get_bus_vmax_kv(),
                                         _grid_model.get_bus_vn_kv(), _grid_model.get_substations(),
                                         _violation_threshold_, nullptr, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, no_skip, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V_n, active_layout().id_me_to_solver, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, no_skip, _violations_n_);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _record_n_case_violations(const CplxVect &){}

        // Once per compute(), before anything is solved: the one capability check (made
        // here, once, rather than discovered row by row) and this call's result buffers.
        // Sized even though the plans are not built yet, so that a batch whose base case
        // diverges answers get_physical_violations() with one empty entry per row rather
        // than with an empty list.
        void _prepare_physical_check(size_t nb_steps, bool ac_solver_used){
            _bus_q_plan_.clear();
            _hvdc_p_plan_.clear();
            _gen_p_plan_.clear();
            _physical_violations_.clear();
            _bus_q_check_on_ = false;
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
                _bus_q_check_on_ = true;
            }
            _physical_violations_.assign(nb_steps, std::vector<LimitViolation>());
        }

        // ... and the routing itself, once the labelling and the control plan this batch
        // solves against are settled (`active_layout()`, which the L1 / L2 preparation
        // above has built or kept). A plan means nothing against another labelling, so it
        // is never carried across a compute().
        void _build_physical_plans(){
            if(!_compute_physical_violations_) return;
            if(_bus_q_check_on_){
                bus_q_check::build_bus_q_plan(_grid_model, active_layout().id_me_to_solver,
                                              active_layout().voltage_control.controllers(),
                                              _bus_q_plan_);
            }
            hvdc_p_check::build_hvdc_p_plan(_grid_model, active_layout().id_me_to_solver,
                                            _hvdc_p_plan_);
            gen_p_check::build_gen_p_plan(_grid_model, active_layout().id_me_to_solver,
                                          _gen_p_plan_);
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
            if(_bus_q_check_on_ && !_bus_q_plan_.empty()){
                const RealVect ctrl_q = _bus_q_plan_.needs_controller_q ? _algo.get_controller_q()
                                                                       : RealVect();
                bus_q_check::check_bus_q_violations(
                    _bus_q_plan_, _grid_model, _algo.get_bus_mismatch(), _algo.get_V(), ctrl_q,
                    _grid_model.get_sn_mva(), _physical_tol_mva_, nullptr,
                    [](int){ return false; },  // the base case disconnects no generator
                    _physical_violations_n_);
            }
            if(!_hvdc_p_plan_.empty()){
                hvdc_p_check::check_hvdc_p_violations(_hvdc_p_plan_, _algo.get_Va(),
                                                      _physical_tol_mva_, nullptr,
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
                    [this](int gen_id){ return this->_grid_target_p(gen_id); },
                    [](int){ return false; },  // the base case disconnects no generator
                    _physical_violations_n_);
            }
        }

        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_masks(){
            if(!_mask_mode()) return;
            // Masking is a value-level edit at constant sparsity, and one that can
            // move a pivot: the row of a stranded lone controller goes from "Vm(reg)
            // = Vset" to "Q_c = 0", so the entry KLU pivoted at in the base
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
                std::ostringstream exc_;
                exc_ << algo_name() << ": the `handle_disconnected_grid` mode requires a "
                        "Newton-Raphson algorithm (AC) or the DC solver (the active algorithm "
                        "does not support bus masking). Use `change_algorithm` to select an NR "
                        "solver (e.g. NR_KLU / NR_SLU) or the DC solver.";
                throw std::runtime_error(exc_.str());
            }
            _select_ref_slack_and_masks();
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _maybe_prepare_masks(){}

        // ================= topological actions (ScenarioSweep only) ==============
        // Turns topo_actions_ into what the row loop consumes, once per compute()
        // (L2, before _prepare_ybus_varying): the branches each row disconnects
        // (ybus_policy_.topo_branches_off -- read wherever the masks are, so the Ybus
        // edit, the connectivity, the current checks and the flow cleaning all see
        // them), the generators (sbus_policy_.topo_gen_off -- read wherever gen_off
        // is, so the injection, the slack re-weighting and the PV -> PQ pinning all
        // see them) and the loads / storage units (sbus_policy_.topo_loads_off /
        // topo_storages_off, injection only).
        //
        // What is refused here, and why, is the frontier of this version: a row that
        // MOVES an element (set_bus > 0) or reconnects one changes the bus set or the
        // Ybus pattern, which the fixed solver layout cannot take yet (TODO in the
        // changelog); a row naming an element in both a mask and its action is
        // ambiguous; the DC path is not wired.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _maybe_resolve_topology(size_t nb_steps){
            _row_topo_.clear();
            _topology_active_ = false;
            _base_masked_.clear();
            _row_branch_overrides_.clear();
            _reset_topo_policy_state();
            if(topo_actions_.empty()) return;
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
            _row_topo_.resize(nb_steps);
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
                RowTopoPlan & plan = _row_topo_[row];
                resolve_row_topo(topo_actions_[row], _grid_model, plan);
                if(plan.empty()) continue;
                _topology_active_ = true;
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
                    if(was_on && sbus_policy_.gen_off.rows() > 0 && sbus_policy_.gen_off(static_cast<Eigen::Index>(row), g.el_id)){
                        _throw_topo_overlap(row, "generator", g.el_id, "set_contingency_gens");
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
                    if(_keep_jacobian_ || _compute_physical_violations_){
                        std::ostringstream exc_;
                        exc_ << algo_name() << "::set_topo_actions: the action of row " << row
                             << (was_on ? " moves" : " reactivates") << " generator " << g.el_id
                             << ", which `keep_jacobian` and `compute_physical_violations` do not support yet "
                                "(their plans are keyed on the base placement of the generators). Turn them off.";
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
            if(topo_actions_.empty()) return;
            std::set<int> used;
            std::set<std::pair<int, int> > edges;
            RowTopoPlan plan;
            for(const TopoAction & act : topo_actions_){
                resolve_row_topo(act, _grid_model, plan);
                for(const TopoElPlacement & el : plan.loads) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoElPlacement & el : plan.gens) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoElPlacement & el : plan.storages) if(el.row_bus_me >= 0) used.insert(el.row_bus_me);
                for(const TopoBranchPlacement & br : plan.branches){
                    if(!br.row_on) continue;
                    used.insert(br.row_bus1_me);
                    used.insert(br.row_bus2_me);
                    edges.insert(std::make_pair(std::min(br.row_bus1_me, br.row_bus2_me),
                                                std::max(br.row_bus1_me, br.row_bus2_me)));
                }
            }
            _topo_used_buses_me_.assign(used.begin(), used.end());
            _union_edges_me_.assign(edges.begin(), edges.end());
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
        void _maybe_seed_extra_buses(CplxVect & V) const {
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

        // the masked row loop runs where handle_disconnected_grid says so, and where
        // topological actions are registered (their extra buses are masked when unused)
        bool _mask_mode() const { return _handle_disconnected_grid || !topo_actions_.empty(); }

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
            sbus_policy_.clear_topo();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && S::supports_vary), int>::type = 0>
        void _reset_topo_policy_state(){}

        // ================= generator contingencies (ScenarioSweep only) ==========
        // Turns sbus_policy_.gen_off into the two things the per-row loop needs, once
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
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
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
            // axes are read together through sbus_policy_.gen_off_in. The row's action
            // may also put a base-off generator back (topo_gens_on): its bus turns PV.
            const bool has_gens_on = sbus_policy_.has_gens_on();
            if(!sbus_policy_.has_gen_off() && !has_gens_on){
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
            _row_slack_gens_off_.assign(nb_steps, std::vector<int>());
            for(Eigen::Index step = 0; step < nb_rows; ++step){
                for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                    if(!sbus_policy_.gen_off_in(step, gen_id)) continue;
                    if(abs(generators.get_gen_slack_weight(gen_id)) < BaseConstants::_tol_equal_float) continue;
                    _row_slack_gens_off_[static_cast<size_t>(step)].push_back(gen_id);
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
                exc_ << algo_name() << ": generator contingencies (`set_contingency_gens` / "
                        "`set_topo_actions`) require a "
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
            for(int gen_id = 0; gen_id < nb_gen; ++gen_id){
                bool ever_off = false;
                for(Eigen::Index step = 0; step < nb_rows; ++step){
                    if(sbus_policy_.gen_off_in(step, gen_id)){ ever_off = true; break; }
                }
                if(!ever_off) continue;
                const int bus_id_me = generators.get_bus_id()(gen_id).cast_int();
                const bool remote = generators.is_remote_voltage_controller(gen_id);
                const bool in_group = bus_id_me != BaseConstants::_deactivated_bus_id &&
                                      group_buses.find(bus_id_me) != group_buses.end();
                if(remote || in_group){
                    std::ostringstream exc_;
                    exc_ << algo_name() << "::set_contingency_gens: generator " << gen_id
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

            // ---- 3. per row: the buses that lose EVERY one of their controllers ----
            _row_pv_to_pq_.assign(nb_steps, std::vector<int>());
            for(Eigen::Index step = 0; step < nb_rows; ++step){
                for(const auto & bus_gens : gens_of_bus){
                    bool all_off = true;
                    for(int gen_id : bus_gens.second){
                        if(!sbus_policy_.gen_off_in(step, gen_id)){
                            all_off = false;
                            break;
                        }
                    }
                    if(all_off) _row_pv_to_pq_[static_cast<size_t>(step)].push_back(bus_gens.first);
                }
            }

            // ---- 3b. per row: the buses a reactivated generator pins ----------------
            // A base-off generator the row puts back on its bus (a bus of the layout,
            // see _maybe_resolve_topology) regulating that bus locally makes it PV for
            // the row. A base-PQ bus already owns a Vm unknown and a Q equation, so
            // nothing is reserved: the row pins its Q row (set_pv_pinned_buses) and
            // seeds |V| there at the set-point (_row_vm_seed_). A base-PV bus whose own
            // controllers the row takes out stays PV through it -- provided both ask
            // the same magnitude.
            std::vector<std::set<int> > extra_pv(nb_steps);
            _row_vm_seed_.assign(nb_steps, std::vector<std::pair<int, real_type> >());
            if(has_gens_on){
                const std::vector<int> & slack_solver = active_layout().slack_bus_id_solver.to_int_vector();
                for(size_t step = 0; step < nb_steps && step < sbus_policy_.topo_gens_on.size(); ++step){
                    std::vector<int> & to_pq = _row_pv_to_pq_[step];
                    for(const auto & gen_on : sbus_policy_.topo_gens_on[step]){
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
                        // the magnitude this generator asks of the bus: the row's own
                        // (modify_gen_v) or the grid's
                        real_type vm = generators.get_target_vm_pu(gen_id);
                        if(sbus_policy_.gen_v.rows() > 0 && static_cast<Eigen::Index>(step) < sbus_policy_.gen_v.rows()){
                            const real_type own = sbus_policy_.gen_v(static_cast<Eigen::Index>(step), gen_id);
                            if(std::isfinite(own)) vm = own;
                        }
                        const auto base_gens = gens_of_bus.find(bus_solver);
                        if(base_gens != gens_of_bus.end()){
                            // a base-PV bus: kept PV whatever the row does to its own
                            // controllers, whose set-point it must share
                            real_type bus_vm = generators.get_target_vm_pu(base_gens->second[0]);
                            if(sbus_policy_.gen_v.rows() > 0 && static_cast<Eigen::Index>(step) < sbus_policy_.gen_v.rows()){
                                const real_type own = sbus_policy_.gen_v(static_cast<Eigen::Index>(step), base_gens->second[0]);
                                if(std::isfinite(own)) bus_vm = own;
                            }
                            if(std::abs(bus_vm - vm) > BaseConstants::_tol_equal_float){
                                std::ostringstream exc_;
                                exc_ << algo_name() << "::set_topo_actions: the action of row " << step
                                     << " reactivates generator " << gen_id << " (set-point " << vm
                                     << " pu) on a bus another generator holds at " << bus_vm
                                     << " pu. A bus has one magnitude: give them the same set-point.";
                                throw std::runtime_error(exc_.str());
                            }
                            std::vector<int>::iterator it = std::find(to_pq.begin(), to_pq.end(), bus_solver);
                            if(it != to_pq.end()) to_pq.erase(it);
                            continue;
                        }
                        // a base-PQ bus: PV for this row, at this magnitude
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
                        extra_pv[step].insert(bus_solver);
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
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && S::supports_vary), int>::type = 0>
        void _maybe_prepare_gen_contingency(size_t){}

        // true when this compute() has generator contingencies to apply at all; every
        // per-row helper below is skipped (and the loop stays bit-identical to before
        // this feature) when it is false.
        bool _has_gen_contingency() const { return _gen_contingency_active_; }
        // ... and true when some bus can actually flip PV -> PQ, so the per-row
        // set_pv_pinned_buses calls are worth making. False in DC, and false when the
        // masked generators happen never to leave a bus without a controller.
        bool _has_pv_switching() const { return _pv_pinning_active_; }

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
        void _apply_step_topo_seed(size_t i, CplxVect & V) const {
            if(i >= _row_vm_seed_.size()) return;
            for(const auto & seed : _row_vm_seed_[i]){
                const int b = seed.first;
                if(b < 0 || b >= V.size()) continue;
                const real_type vm = std::abs(V(b));
                V(b) = vm > 0. ? V(b) * (seed.second / vm) : cplx_type(seed.second, 0.);
            }
        }

        // row i's distributed-slack weights: the layout's own unless the row takes a
        // participating generator out, in which case they are re-derived without it
        // and renormalised (LSGrid::get_slack_weights_solver_without -- the storage
        // units taking part in the slack stay in: no row disconnects one).
        // Should a row somehow leave no participant at all, the reference slack bus
        // keeps the whole share -- the angle reference is a property of the batch,
        // picked once, and must not move from row to row.
        // Returns the layout's own vector by reference on the common row; the
        // re-derived one is written into the caller's `scratch` (one per range, so
        // no row allocates and no two threads share it).
        const RealVect & _row_slack_weights(size_t i, RealVect & scratch) const {
            const RealVect & base_w = active_layout().slack_weights;
            if(i >= _row_slack_gens_off_.size() || _row_slack_gens_off_[i].empty()) return base_w;
            const auto & generators = _grid_model.get_generators();
            std::vector<bool> gen_off(generators.nb(), false);
            for(int gen_id : _row_slack_gens_off_[i]) gen_off[gen_id] = true;
            scratch = _grid_model.get_slack_weights_solver_without(
                static_cast<size_t>(base_w.size()), active_layout().id_me_to_solver, gen_off);
            if(abs(scratch.sum()) < BaseConstants::_tol_equal_float){
                scratch.setZero();
                if(active_layout().slack_bus_id_solver.size() > 0){
                    const int ref = active_layout().slack_bus_id_solver[static_cast<int>(0)].cast_int();
                    if(ref >= 0 && ref < scratch.size()) scratch(ref) = 1.;
                }
            }
            return scratch;
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
                           bool & conv, bool & invertible);

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
                        _apply_step_vc_v_set(cont_id, algo);
                        const RealVect & sw = _masked_slack_weights(masked, _row_slack_weights(cont_id, sw_scratch), sw_scratch);
                        const CplxVect & sb = _step_sbus(cont_id, sbus_scratch);
                        conv = compute_one_powerflow(algo, control, nb_solved, nb_converged, timer_solver, Ybus, V, sb,
                                                     active_layout().slack_bus_id_solver.as_eigen(), sw,
                                                     active_layout().bus_pv.as_eigen(), active_layout().bus_pq.as_eigen(), max_iter, tol / sn_mva);
                        if(needs_solver_init){ control.tell_none_changed(); needs_solver_init = false; }
                        // before the two restores below, and before the Ybus is put
                        // back: see _maybe_store_jacobian (and _record_row_physical, which
                        // reads the mismatch of the system this row solved)
                        if(conv){
                            _maybe_store_jacobian(cont_id, algo);
                            _record_row_physical(cont_id, algo, V, sw, sb);
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
                        "set_contingency_gens / set_topo_actions) after the "
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
            for(const auto & l_id_this_cont: ybus_policy_.li_defaults){
                for(auto l_id : l_id_this_cont){
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
        bool _gen_contingency_active_ = false;
        std::vector<int> _switchable_buses_;
        std::vector<std::vector<int> > _row_pv_to_pq_;
        std::vector<std::vector<int> > _row_pv_pinned_;
        std::vector<std::vector<std::pair<int, real_type> > > _row_vm_seed_;
        bool _pv_pinning_active_ = false;
        std::vector<std::vector<int> > _row_slack_gens_off_;
        // topological actions (ScenarioSweep only; plain, always-present state like
        // the above). topo_actions_ is a registration (one checked TopoAction per
        // row, see set_topo_actions); _row_topo_ / _topology_active_ are what
        // _maybe_resolve_topology made of it for this batch (L2).
        std::vector<TopoAction> topo_actions_;
        std::vector<RowTopoPlan> _row_topo_;
        bool _topology_active_ = false;
        // L1: the buses the actions use (gridmodel ids, sorted), those of them empty in
        // the base grid, and the (bus, bus) pairs a row's branch placement writes
        std::vector<int> _topo_used_buses_me_;
        std::vector<int> _extra_buses_me_;
        std::vector<std::pair<int, int> > _union_edges_me_;
        // L2: the extra buses in solver numbering (sorted): masked in the "n" case and in
        // every row that does not use them -- the masked loop's resting state
        std::vector<int> _base_masked_;
        // per contingency, the solver buses it strands (sorted, empty if none) and
        // whether the grid stays connected -- both settled by _prepare_connectivity
        // from bus_graph_, the DFS tree of the base graph built there.
        std::vector<std::vector<int> > _li_masked;
        std::vector<char> _cont_connected_;
        BusGraph bus_graph_;
        std::vector<char> _skip_mask;
        bool _compute_limit_violations_ = false;
        real_type _violation_threshold_ = 1.0;
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
        // labelling, so they are dropped with the results. `_bus_q_check_on_` is whether
        // the algorithm of THIS compute() can feed the reactive half at all (see
        // _prepare_physical_check): false in DC, where there is no reactive power to check.
        bool _compute_physical_violations_ = false;
        real_type _physical_tol_mva_ = 1e-4;
        bool _bus_q_check_on_ = false;
        bus_q_check::BusQPlan _bus_q_plan_;
        hvdc_p_check::HvdcPPlan _hvdc_p_plan_;
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
