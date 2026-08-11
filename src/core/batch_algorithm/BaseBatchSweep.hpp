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
#include "SbusPolicy.hpp"
#include "LimitViolation.hpp"

#include <set>
#include <vector>
#include <queue>
#include <algorithm>
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
    std::vector<LimitViolation> & out)
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

    for(size_t el_id = 0; el_id < nb_el; ++el_id){
        if(!el_status[el_id]) continue;
        if(!skip.empty() && skip[el_id]) continue;

        const Eigen::Index el_idx = static_cast<Eigen::Index>(el_id);
        const bool has_lim1 = limit1.size() > 0 && !isnan(limit1(el_idx));
        const bool has_lim2 = limit2.size() > 0 && !isnan(limit2(el_idx));
        if(!has_lim1 && !has_lim2) continue;

        const bool s1 = status1[el_id];
        const bool s2 = status2[el_id];
        int bus_from_me = BaseConstants::_deactivated_bus_id;
        int bus_to_me = BaseConstants::_deactivated_bus_id;
        int solver_from = BaseConstants::_deactivated_bus_id;
        int solver_to = BaseConstants::_deactivated_bus_id;
        if(s1){
            bus_from_me = bus_from(static_cast<int>(el_id)).cast_int();
            solver_from = id_me_to_solver(bus_from_me).cast_int();
            if(solver_from == BaseConstants::_deactivated_bus_id) continue;
        }
        if(s2){
            bus_to_me = bus_to(static_cast<int>(el_id)).cast_int();
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
            const cplx_type I_from = std::conj(yac_eff_11(el_idx) * Efrom + yac_eff_12(el_idx) * Eto);
            const cplx_type S_from = Efrom * I_from;
            amps1 = std::abs(S_from) * sn_mva / (sqrt_3 * v_from_kv);

            const cplx_type I_to = std::conj(yac_eff_22(el_idx) * Eto + yac_eff_21(el_idx) * Efrom);
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

        const std::string el_name = el_id < el_names.size() ? el_names[el_id] : std::string();
        if(has_lim1 && amps1 >= threshold * limit1(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 1,
                                          LimitViolationType::CURRENT, amps1, limit1(el_idx), el_name});
        }
        if(has_lim2 && amps2 >= threshold * limit2(el_idx)){
            out.push_back(LimitViolation{el_type, static_cast<int>(el_id), 2,
                                          LimitViolationType::CURRENT, amps2, limit2(el_idx), el_name});
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

        void clear() override {
            BaseBatchSolverSynch::clear();
            sbus_policy_.clear();
            _reset_ybus_policy();
            _status = 1;
            _nb_steps_locked = -1;
            _handle_disconnected_grid = false;
            _li_masked.clear();
            _skip_mask.clear();
            _converged.clear();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _li_defaults_vect_cache_.clear();
            _timer_modif_Ybus = 0.;
            _timer_total = 0.;
            _timer_pre_proc = 0.;
        }

        // ContingencyAnalysis-only: like clear() but keeps the registered
        // contingencies, only drops computed results (mirrors the pre-refactor
        // ContingencyAnalysis::clear_results_only).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void clear_results_only() {
            BaseBatchSolverSynch::clear();
            _li_masked.clear();
            _skip_mask.clear();
            _converged.clear();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _li_defaults_vect_cache_.clear();
            _timer_total = 0.;
            _timer_modif_Ybus = 0.;
            _timer_thread_init = 0.;
            _timer_pre_proc = 0.;
        }

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

        // ============ contingency mask setters (Contingency && Vary only) ========
        // ScenarioSweep only: row i, column j True means "deactivate this branch for
        // step i". See YbusPolicy::Contingency::line_mask/trafo_mask.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_lines(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(n_line_), "set_contingency_lines");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_lines");
            ybus_policy_.line_mask = mask;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void set_contingency_trafos(const Eigen::Ref<const BoolMat> & mask) {
            _check_cols_bool(mask, static_cast<Eigen::Index>(n_trafos_), "set_contingency_trafos");
            _lock_or_check_nb_steps(mask.rows(), "set_contingency_trafos");
            ybus_policy_.trafo_mask = mask;
        }

        // ========== legacy contingency-list API (Contingency && !Vary only) ======
        // ContingencyAnalysis only -- deliberately NOT unified with
        // set_contingency_lines/trafos above (two genuinely different usages -- a
        // shared base case with an arbitrary set of distinct scenarios here, vs. one
        // contingency per row paired with that row's own injection on ScenarioSweep).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_all_n1(){
            clear_results_only();
            for(int l_id = 0; l_id < static_cast<int>(n_total_); ++l_id){
                std::set<int> this_default = {l_id};
                ybus_policy_.li_defaults.insert(this_default);
            }
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_n1(int line_id){
            _check_ok_el(line_id);
            clear_results_only();
            std::set<int> this_default = {line_id};
            ybus_policy_.li_defaults.insert(this_default);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void add_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s) _check_ok_el(line_id);
            clear_results_only();
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
            clear_results_only();
            ybus_policy_.li_defaults.insert(this_default);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool remove_n1(int line_id){
            _check_ok_el(line_id);
            clear_results_only();
            std::set<int> this_default = {line_id};
            auto nb_removed = ybus_policy_.li_defaults.erase(this_default);
            return nb_removed >= 1;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        size_t remove_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s) _check_ok_el(line_id);
            clear_results_only();
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
            clear_results_only();
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

        // "handle disconnected grid" mode (ContingencyAnalysis only)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool get_handle_disconnected_grid() const {return _handle_disconnected_grid;}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void set_handle_disconnected_grid(bool val) {_handle_disconnected_grid = val;}

        // limit violations (ContingencyAnalysis only)
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        bool get_compute_limit_violations() const noexcept {return _compute_limit_violations_;}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void set_compute_limit_violations(bool val){
            if(val == _compute_limit_violations_) return;
            _compute_limit_violations_ = val;
            clear();
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        real_type get_violation_threshold() const noexcept {return _violation_threshold_;}
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void set_violation_threshold(real_type val){
            if(!(val > 0. && val <= 1.)){
                std::ostringstream exc_;
                exc_ << algo_name() << "::set_violation_threshold: the threshold should be "
                        "a real number in the range ]0., 1.] (got " << val << ").";
                throw std::runtime_error(exc_.str());
            }
            if(val < _violation_threshold_) clear_results_only();
            _violation_threshold_ = val;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::vector<char> & converged() const {
            _check_limit_violations_enabled("converged");
            return _converged;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
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
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        const std::vector<LimitViolation> & get_violations_n() const {
            _check_limit_violations_enabled("get_violations_n");
            return _violations_n_;
        }

        // contingency-list inspection (ContingencyAnalysis only). Public (directly
        // bound to Python) -- see the class-level note on why this needs a fully
        // inline body, including everything it transitively calls.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        IntVect is_grid_connected_after_contingency(){
            const bool ac_solver_used = _algo.ac_solver_used();
            const bool inputs_ready = (ac_solver_used ? Ybus_.cols() != 0 : Bbus_.cols() != 0);
            if(!inputs_ready || ybus_policy_.li_coeffs.size() != ybus_policy_.li_defaults.size()){
                const size_t nb_total_bus = _grid_model.total_bus();
                CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                                    {_grid_model.get_init_vm_pu(), 0.});
                prepare_solver_input_base(Vinit, ac_solver_used);
                ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, id_me_to_solver_, n_line_);
            }
            IntVect res = IntVect::Constant(ybus_policy_.li_coeffs.size(), 0);
            if(ac_solver_used){
                Eigen::SparseMatrix<cplx_type> Ybus = Ybus_;
                int cont_id = 0;
                for(const auto & coeffs_modif: ybus_policy_.li_coeffs){
                    if(YbusPolicy::Contingency::remove_from_Ybus(Ybus, coeffs_modif, true, _algo)) res(cont_id) = 1;
                    else res(cont_id) = 0;
                    YbusPolicy::Contingency::readd_to_Ybus(Ybus, coeffs_modif, true, _algo);
                    ++cont_id;
                }
            } else {
                Eigen::SparseMatrix<real_type> Bbus = Bbus_;
                int cont_id = 0;
                for(const auto & coeffs_modif: ybus_policy_.li_coeffs){
                    for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) -= std::real(c.value);
                    res(cont_id) = _disconnected_buses(Bbus).empty() ? 1 : 0;
                    for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) += std::real(c.value);
                    ++cont_id;
                }
            }
            return res;
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        int pick_reference_slack(){
            const bool ac_solver_used = _algo.ac_solver_used();
            const bool inputs_ready = (ac_solver_used ? Ybus_.cols() != 0 : Bbus_.cols() != 0);
            if(!inputs_ready || ybus_policy_.li_coeffs.size() != ybus_policy_.li_defaults.size()){
                const size_t nb_total_bus = _grid_model.total_bus();
                CplxVect Vinit = CplxVect::Constant(static_cast<Eigen::Index>(nb_total_bus),
                                                    {_grid_model.get_init_vm_pu(), 0.});
                prepare_solver_input_base(Vinit, ac_solver_used);
                ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, id_me_to_solver_, n_line_);
            }
            _select_ref_slack_and_masks();
            if(slack_ids_me_.size() == 0) return -1;
            return slack_ids_me_[static_cast<int>(0)].cast_int();
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
        template<class S = SbusPolicy, class Y = YbusPolicy,
                 typename std::enable_if<S::supports_vary && !Y::supports_contingency, int>::type = 0>
        Eigen::Ref<const typename S::CplxMat> get_sbuses() const {return sbus_policy_.sbuses;}

        // aggregate status (TimeSeries/InjectionSweep/ScenarioSweep -- everything
        // with SbusPolicy::supports_vary; ContingencyAnalysis reports per-row via
        // converged() instead, it has no aggregate status).
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        int get_status() const { return _status; }

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
        RealVect _masked_slack_weights(const std::vector<int> & masked) const {
            RealVect w = slack_weights_;
            if(masked.empty()) return w;
            const real_type orig_sum = w.sum();
            for(int b : masked) if(b >= 0 && b < w.size()) w(b) = 0.;
            const real_type new_sum = w.sum();
            if(new_sum > 1e-12 && orig_sum > 1e-12) w *= (orig_sum / new_sum);
            return w;
        }

        // BFS connected-component labelling: returns the solver bus ids NOT part of
        // the largest connected component (empty if the whole matrix is connected).
        // Called from is_grid_connected_after_contingency()/pick_reference_slack(),
        // both public, hence this plain (non-SFINAE, but member-template-on-T)
        // helper must also be fully inline.
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

        // pre-pass run before the (n-)powerflow: fills _li_masked (the masked bus
        // set of each contingency), chooses the reference slack that minimises the
        // number of skipped contingencies (reordering slack_ids_me_ so it is index 0)
        // and fills _skip_mask. ContingencyAnalysis only; called from
        // _maybe_prepare_masks() (compute()-only) and pick_reference_slack() (public,
        // hence this too must be fully inline).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _select_ref_slack_and_masks(){
            const size_t nb_cont = ybus_policy_.li_coeffs.size();
            _li_masked.assign(nb_cont, std::vector<int>());
            _skip_mask.assign(nb_cont, 0);

            if(_algo.ac_solver_used()){
                Eigen::SparseMatrix<cplx_type> Ybus = Ybus_;
                size_t cont_id = 0;
                for(const auto & coeffs_modif: ybus_policy_.li_coeffs){
                    for(const auto & c: coeffs_modif) Ybus.coeffRef(c.row_id, c.col_id) -= c.value;
                    _li_masked[cont_id] = _disconnected_buses(Ybus);
                    for(const auto & c: coeffs_modif) Ybus.coeffRef(c.row_id, c.col_id) += c.value;
                    ++cont_id;
                }
            } else {
                Eigen::SparseMatrix<real_type> Bbus = Bbus_;
                size_t cont_id = 0;
                for(const auto & coeffs_modif: ybus_policy_.li_coeffs){
                    for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) -= std::real(c.value);
                    _li_masked[cont_id] = _disconnected_buses(Bbus);
                    for(const auto & c: coeffs_modif) Bbus.coeffRef(c.row_id, c.col_id) += std::real(c.value);
                    ++cont_id;
                }
            }

            const Eigen::Index nb_slack = slack_ids_me_.size();
            std::vector<int> candidates;
            for(Eigen::Index i = 0; i < nb_slack; ++i){
                const int b = slack_ids_me_[static_cast<int>(i)].cast_int();
                if(b >= 0 && b < slack_weights_.size() && slack_weights_(b) > 0.) candidates.push_back(b);
            }
            if(candidates.empty()){
                for(Eigen::Index i = 0; i < nb_slack; ++i) candidates.push_back(slack_ids_me_[static_cast<int>(i)].cast_int());
            }
            if(candidates.empty()) return;

            auto is_masked = [](const std::vector<int> & masked, int bus){
                return std::find(masked.begin(), masked.end(), bus) != masked.end();
            };
            int best_bus = candidates[0];
            int best_strand = -1;
            real_type best_weight = -1.;
            for(int bus : candidates){
                int strand = 0;
                for(const auto & masked : _li_masked) if(is_masked(masked, bus)) ++strand;
                const real_type weight = (bus < slack_weights_.size()) ? slack_weights_(bus) : 0.;
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
                    const int b = slack_ids_me_[static_cast<int>(i)].cast_int();
                    if(b != best_bus) reordered.push_back(b);
                }
                slack_ids_me_ = GlobalBusIdVect(reordered);
            }

            for(size_t cont_id = 0; cont_id < nb_cont; ++cont_id){
                if(is_masked(_li_masked[cont_id], best_bus)) _skip_mask[cont_id] = 1;
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
            ybus_policy_.init_li_coeffs(_grid_model, ac_solver_used, id_me_to_solver_, n_line_);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && S::supports_vary, int>::type = 0>
        void _prepare_ybus_varying(bool ac_solver_used, Eigen::Index nb_steps) {
            ybus_policy_.init_li_coeffs_from_masks(_grid_model, ac_solver_used, id_me_to_solver_, n_line_, n_trafos_, nb_steps);
        }

        // ----- Sbus-varying preparation: 2-way --------------------------------
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        void _prepare_sbus_varying(bool, Eigen::Index) {}
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        void _prepare_sbus_varying(bool ac_solver_used, Eigen::Index nb_steps) {
            const CplxVect complete = ac_solver_used ? CplxVect(Sbus_) : CplxVect(Pbus_.template cast<cplx_type>());
            sbus_policy_.assemble(_grid_model, ac_solver_used, nb_buses_solver_, id_me_to_solver_,
                                  complete, _grid_model.get_sn_mva(), nb_steps, algo_name());
        }

        // ----- reset ybus_policy_'s own state: 2-way --------------------------
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
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        bool _remove_step_coeffs(Eigen::SparseMatrix<cplx_type> & Ybus, size_t i, bool ac_solver_used, AlgorithmSelector & algo) {
            return YbusPolicy::Contingency::remove_from_Ybus(Ybus, ybus_policy_.li_coeffs[i], ac_solver_used, algo);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        bool _remove_step_coeffs(Eigen::SparseMatrix<cplx_type> &, size_t, bool, AlgorithmSelector &) { return true; }
        template<class Y = YbusPolicy, typename std::enable_if<Y::supports_contingency, int>::type = 0>
        void _readd_step_coeffs(Eigen::SparseMatrix<cplx_type> & Ybus, size_t i, bool ac_solver_used, AlgorithmSelector & algo) {
            YbusPolicy::Contingency::readd_to_Ybus(Ybus, ybus_policy_.li_coeffs[i], ac_solver_used, algo);
        }
        template<class Y = YbusPolicy, typename std::enable_if<!Y::supports_contingency, int>::type = 0>
        void _readd_step_coeffs(Eigen::SparseMatrix<cplx_type> &, size_t, bool, AlgorithmSelector &) {}

        // ----- per-step Sbus: 2-way (row of sbus_policy_ vs. the fixed member) ---
        template<class S = SbusPolicy, typename std::enable_if<S::supports_vary, int>::type = 0>
        CplxVect _step_sbus(size_t i) const { return sbus_policy_.sbuses.row(i); }
        template<class S = SbusPolicy, typename std::enable_if<!S::supports_vary, int>::type = 0>
        CplxVect _step_sbus(size_t) const { return Sbus_; }

        // ----- per-row violation recording ---------------------------------------
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _record_row_violations(size_t i, const CplxVect & V, const std::vector<int> * masked_ids){
            const bool ac_solver_used = _algo.ac_solver_used();
            const real_type sn_mva = _grid_model.get_sn_mva();

            std::vector<int> skip_lines, skip_trafos;
            if(i < _li_defaults_vect_cache_.size()){
                for(int br_id : _li_defaults_vect_cache_[i]){
                    if(static_cast<size_t>(br_id) < n_line_) skip_lines.push_back(br_id);
                    else skip_trafos.push_back(static_cast<int>(br_id - static_cast<int>(n_line_)));
                }
            }
            const std::vector<int> * masked_ids_use = (masked_ids != nullptr && !masked_ids->empty()) ? masked_ids : nullptr;

            batch_sweep_detail::check_bus_voltage_violations(V, id_me_to_solver_, _grid_model.get_bus_vmin_kv(), _grid_model.get_bus_vmax_kv(),
                                         _grid_model.get_bus_vn_kv(), _grid_model.get_substations(),
                                         _violation_threshold_, masked_ids_use, _violations[i]);
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V, id_me_to_solver_, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, skip_lines, _violations[i]);
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V, id_me_to_solver_, _grid_model.get_bus_vn_kv(), ac_solver_used, sn_mva,
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, skip_trafos, _violations[i]);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _record_row_violations_dispatch(size_t i, const CplxVect & V, const std::vector<int> * masked_ids) {
            _record_row_violations(i, V, masked_ids);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _record_row_violations_dispatch(size_t, const CplxVect &, const std::vector<int> *) {}

        void _store_row_status(size_t i, bool conv, bool invertible, const CplxVect & V_solver){
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

        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _record_n_case_violations(const CplxVect & V_n){
            if(!_compute_limit_violations_) return;
            _converged_n_ = true;
            const std::vector<int> no_skip;
            batch_sweep_detail::check_bus_voltage_violations(V_n, id_me_to_solver_, _grid_model.get_bus_vmin_kv(), _grid_model.get_bus_vmax_kv(),
                                         _grid_model.get_bus_vn_kv(), _grid_model.get_substations(),
                                         _violation_threshold_, nullptr, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_powerlines_as_data(), ViolationElementType::LINE,
                                     V_n, id_me_to_solver_, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_powerlines_as_data().get_limit_a1_ka(),
                                     _grid_model.get_powerlines_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, no_skip, _violations_n_);
            batch_sweep_detail::check_current_violations(_grid_model.get_trafos_as_data(), ViolationElementType::TRAFO,
                                     V_n, id_me_to_solver_, _grid_model.get_bus_vn_kv(), _algo.ac_solver_used(), _grid_model.get_sn_mva(),
                                     _grid_model.get_trafos_as_data().get_limit_a1_ka(),
                                     _grid_model.get_trafos_as_data().get_limit_a2_ka(),
                                     _violation_threshold_, no_skip, _violations_n_);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _record_n_case_violations(const CplxVect &){}

        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _maybe_prepare_masks(){
            if(!_handle_disconnected_grid) return;
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
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _maybe_prepare_masks(){}

        // per-range worker: NON mask-mode path, shared by every instantiation.
        void _run_range(size_t step_begin, size_t step_end,
                        AlgorithmSelector & algo, AlgoControl & control,
                        Eigen::SparseMatrix<cplx_type> & Ybus,
                        const Eigen::Ref<const CplxVect> & Vinit_solver,
                        bool ac_solver_used, int max_iter, real_type tol_solver,
                        int & nb_solved, double & timer_solver, double & timer_modif_ybus,
                        int & first_diverging_step, std::exception_ptr & err,
                        bool needs_solver_init);
        void _run_one_step(size_t i, AlgorithmSelector & algo, AlgoControl & control,
                           Eigen::SparseMatrix<cplx_type> & Ybus, CplxVect & V,
                           bool ac_solver_used, int max_iter, real_type tol_solver,
                           int & nb_solved, double & timer_solver, double & timer_modif_ybus,
                           bool & conv, bool & invertible);

        // per-range worker: mask-mode path (ContingencyAnalysis only). Only ever
        // called from _maybe_run_range_masked below, itself only ever called from
        // compute() (non-template, .cpp-confined) -- safe to stay declared here /
        // defined out-of-line in BaseBatchSweep.cpp.
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _run_range_masked(size_t cont_begin, size_t cont_end,
                               AlgorithmSelector & algo, AlgoControl & control,
                               Eigen::SparseMatrix<cplx_type> & Ybus,
                               const Eigen::Ref<const CplxVect> & Vinit_solver,
                               bool ac_solver_used, int max_iter, real_type tol, real_type sn_mva,
                               double & timer_modif_ybus, int & nb_solved, double & timer_solver,
                               std::exception_ptr & err, bool needs_solver_init)
        {
            try {
                CplxVect V;
                if(needs_solver_init) control.tell_all_changed();

                for(size_t cont_id = cont_begin; cont_id < cont_end; ++cont_id){
                    const std::vector<Coeff> & coeffs_modif = ybus_policy_.li_coeffs[cont_id];
                    bool conv = false;
                    bool do_store = false;
                    LimitViolationType div_reason = LimitViolationType::NOT_SIMULATED;

                    if(!_skip_mask[cont_id]){
                        const std::vector<int> & masked = _li_masked[cont_id];
                        auto t1 = CustTimer();
                        YbusPolicy::Contingency::remove_from_Ybus(Ybus, coeffs_modif, ac_solver_used, algo);
                        timer_modif_ybus += t1.duration();

                        algo.set_masked_buses(masked);
                        V = Vinit_solver;
                        if(masked.empty()){
                            conv = compute_one_powerflow(algo, control, nb_solved, timer_solver, Ybus, V, Sbus_,
                                                         slack_ids_solver_.as_eigen(), slack_weights_,
                                                         bus_pv_.as_eigen(), bus_pq_.as_eigen(), max_iter, tol / sn_mva);
                        } else {
                            const RealVect sw = _masked_slack_weights(masked);
                            conv = compute_one_powerflow(algo, control, nb_solved, timer_solver, Ybus, V, Sbus_,
                                                         slack_ids_solver_.as_eigen(), sw,
                                                         bus_pv_.as_eigen(), bus_pq_.as_eigen(), max_iter, tol / sn_mva);
                        }
                        if(needs_solver_init){ control.tell_none_changed(); needs_solver_init = false; }
                        algo.set_masked_buses(std::vector<int>());

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

                    if(do_store) _voltages.row(cont_id)(id_solver_to_me_.as_eigen()) = V.array();

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
        // dispatch into _run_range_masked only where it exists; unreachable
        // elsewhere (_handle_disconnected_grid can never be true there).
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<Y::supports_contingency && !S::supports_vary, int>::type = 0>
        void _maybe_run_range_masked(size_t cont_begin, size_t cont_end,
                                     AlgorithmSelector & algo, AlgoControl & control,
                                     Eigen::SparseMatrix<cplx_type> & Ybus,
                                     const Eigen::Ref<const CplxVect> & Vinit_solver,
                                     bool ac_solver_used, int max_iter, real_type tol, real_type sn_mva,
                                     double & timer_modif_ybus, int & nb_solved, double & timer_solver,
                                     std::exception_ptr & err, bool needs_solver_init){
            _run_range_masked(cont_begin, cont_end, algo, control, Ybus, Vinit_solver, ac_solver_used,
                              max_iter, tol, sn_mva, timer_modif_ybus, nb_solved, timer_solver, err, needs_solver_init);
        }
        template<class Y = YbusPolicy, class S = SbusPolicy,
                 typename std::enable_if<!(Y::supports_contingency && !S::supports_vary), int>::type = 0>
        void _maybe_run_range_masked(size_t, size_t, AlgorithmSelector &, AlgoControl &,
                                     Eigen::SparseMatrix<cplx_type> &, const Eigen::Ref<const CplxVect> &,
                                     bool, int, real_type, real_type, double &, int &, double &,
                                     std::exception_ptr &, bool){
            throw std::logic_error("unreachable: handle_disconnected_grid cannot be set on this instantiation");
        }

        // ----- threading spawn: 2-way (per-thread Ybus copy vs shared Ybus_) ------
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
            const size_t nb_res = static_cast<size_t>(_voltages.rows());
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
        void _maybe_check_results_match_defaults(const std::string &) const {}

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
            const Eigen::Index nb_steps = _voltages.rows();
            for(Eigen::Index row = 0; row < nb_steps; ++row){
                if(ybus_policy_.line_mask.rows() > 0){
                    for(Eigen::Index col = 0; col < ybus_policy_.line_mask.cols(); ++col){
                        if(!ybus_policy_.line_mask(row, col)) continue;
                        real_type & el = is_amps ? _amps_flows(row, col) : _active_power_flows(row, col);
                        if(isfinite(el)) el = 0.;
                    }
                }
                if(ybus_policy_.trafo_mask.rows() > 0){
                    for(Eigen::Index col = 0; col < ybus_policy_.trafo_mask.cols(); ++col){
                        if(!ybus_policy_.trafo_mask(row, col)) continue;
                        const Eigen::Index l_id = static_cast<Eigen::Index>(n_line_) + col;
                        real_type & el = is_amps ? _amps_flows(row, l_id) : _active_power_flows(row, l_id);
                        if(isfinite(el)) el = 0.;
                    }
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

        // "handle disconnected grid" mode + limit violations: plain, always-present
        // state (SFINAE-gated methods above restrict who can reach it); empty /
        // false on every instantiation but ContingencyAnalysis.
        bool _handle_disconnected_grid = false;
        std::vector<std::vector<int> > _li_masked;
        std::vector<char> _skip_mask;
        bool _compute_limit_violations_ = false;
        real_type _violation_threshold_ = 1.0;
        std::vector<char> _converged;
        std::vector<std::vector<LimitViolation> > _violations;
        bool _converged_n_ = false;
        std::vector<LimitViolation> _violations_n_;
        // per-row branch-id cache (ContingencyAnalysis only), refreshed once per
        // compute() call -- avoids rebuilding my_defaults_vect() per row / per thread.
        std::vector<std::vector<int> > _li_defaults_vect_cache_;

        double _timer_modif_Ybus = 0.;
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
