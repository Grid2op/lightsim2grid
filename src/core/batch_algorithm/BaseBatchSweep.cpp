// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseBatchSweep.hpp"

#include <algorithm>
#include <memory>
#include <sstream>
#include <thread>
#include <vector>

namespace ls2g {

template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
void BaseBatchSweep<YbusPolicy, SbusPolicy, INIT>::_run_one_step(
    size_t i, AlgorithmSelector & algo, AlgoControl & control,
    Eigen::SparseMatrix<cplx_type> & Ybus, CplxVect & V,
    bool ac_solver_used, int max_iter, real_type tol_solver,
    int & nb_solved, int & nb_converged, double & timer_solver, double & timer_modif_ybus,
    bool & conv, bool & invertible)
{
    // re-seed |V| at every voltage-regulating generator's regulated bus with this
    // row's own target (modify_gen_v; a no-op if it was never set, or on a
    // SbusPolicy::NOOP instantiation). Applied unconditionally, AC or DC, matching
    // the existing precedent in BaseBatchSolverSynch::_finish_preprocessing (which
    // also seeds vm before branching AC/DC) -- DC's own equations ignore |V|
    // magnitude, so this is harmless there. Covers both BatchInitKind variants: V
    // already holds "this row's starting point" by the time _run_one_step runs
    // (either freshly reset to Vinit_solver for FromSeed, or carried over from the
    // previous row's convergence for FromPreviousStep / TimeSeries), so a single
    // unconditional call here re-applies the *current* row's target every time.
    _apply_step_gen_v(i, V);

    auto t1 = CustTimer();
    invertible = _remove_step_coeffs(Ybus, i, ac_solver_used, algo);
    timer_modif_ybus += t1.duration();

    if(invertible){
        conv = compute_one_powerflow(algo, control, nb_solved, nb_converged, timer_solver,
                                     Ybus, V, _step_sbus(i),
                                     active_layout().slack_bus_id_solver.as_eigen(), active_layout().slack_weights,
                                     active_layout().bus_pv.as_eigen(), active_layout().bus_pq.as_eigen(),
                                     max_iter, tol_solver);
    } else {
        conv = false;
    }

    auto t2 = CustTimer();
    _readd_step_coeffs(Ybus, i, ac_solver_used, algo);
    timer_modif_ybus += t2.duration();
}

template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
void BaseBatchSweep<YbusPolicy, SbusPolicy, INIT>::_run_range(
    size_t step_begin, size_t step_end,
    AlgorithmSelector & algo, AlgoControl & control,
    Eigen::SparseMatrix<cplx_type> & Ybus,
    const Eigen::Ref<const CplxVect> & Vinit_solver,
    bool ac_solver_used, int max_iter, real_type tol_solver,
    int & nb_solved, int & nb_converged, double & timer_solver, double & timer_modif_ybus,
    int & first_diverging_step, std::exception_ptr & err,
    bool needs_solver_init)
{
    try {
        CplxVect V = Vinit_solver;

        if(needs_solver_init) control.tell_all_changed();
        if(!ac_solver_used) control.tell_recompute_sbus();

        for(size_t i = step_begin; i < step_end; ++i){
            // this single line is the whole difference between the FromSeed and
            // FromPreviousStep instantiations (chained vs. independent steps) --
            // INIT is a compile time constant, only one branch is ever taken.
            if(INIT == BatchInitKind::FromSeed) V = Vinit_solver;

            bool conv = false;
            bool invertible = true;
            _run_one_step(i, algo, control, Ybus, V, ac_solver_used, max_iter, tol_solver,
                          nb_solved, nb_converged, timer_solver, timer_modif_ybus, conv, invertible);

            control.tell_none_changed();
            if(!ac_solver_used) control.tell_recompute_sbus();

            if(conv && _dc_lazy_storage_used_ && _compute_limit_violations_){
                // V's angle is stale in the DC fast path (never refreshed from the
                // solve -- see compute_one_powerflow / BaseAlgo::set_lazy_v):
                // current-limit checks need the real per-bus angle, so rebuild it
                // here for just this row, from the fresh theta (algo.get_Va()) and
                // V's own magnitude (still correct -- DC never changes |V|, see
                // BaseDCAlgo::compute_pf_dc). Paid only when limit violations are
                // actually requested.
                const RealVect Vm_local = V.array().abs();
                const auto & Va_local = algo.get_Va();
                CplxVect V_violations(V.size());
                for(Eigen::Index k = 0; k < V_violations.size(); ++k){
                    V_violations[k] = std::polar(Vm_local[k], Va_local[k]);
                }
                _store_row_status(i, conv, invertible, V_violations);
            } else {
                _store_row_status(i, conv, invertible, V);
            }

            if(!conv){
                first_diverging_step = static_cast<int>(i);
                // A diverging/non-invertible step stops this whole range only when
                // steps are chained (TimeSeries, INIT == FromPreviousStep): a later
                // row's Vinit depends on this one, so there is nothing meaningful left
                // to compute past the failure. Every other instantiation (InjectionSweep,
                // ContingencyAnalysis, ScenarioSweep) uses BatchInitKind::FromSeed --
                // its rows are independent of one another and of their ordering (see
                // eg InjectionSweep's own docstring), so one failing row must not
                // prevent the rest of the range from being computed. Aborting there
                // would silently leave every later row as an unattempted zero row,
                // indistinguishable from an actual divergence. A skip is just one more
                // zero row / NOT_SIMULATED-or-DIVERGENCE violation, already recorded
                // above by _store_row_status.
                if(INIT == BatchInitKind::FromPreviousStep) return;
                continue;
            }
            if(_dc_lazy_storage_used_){
                _thetas.row(i)(active_layout().id_solver_to_me.as_eigen()) = algo.get_Va().array();
                // marks this row as actually solved -- see _dc_row_solved_ / _dc_vm_row_grid:
                // a row that never reaches here (eg an islanding DC contingency, which
                // diverges and is skipped above) must reconstruct as exact complex 0, not
                // as its (never written, so 0) theta paired with a nonzero base magnitude.
                if(static_cast<size_t>(i) < _dc_row_solved_.size()) _dc_row_solved_[i] = 1;
            } else {
                _voltages.row(i)(active_layout().id_solver_to_me.as_eigen()) = V.array();
            }
        }
    } catch(...) {
        err = std::current_exception();
    }
}

template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
template<class Y, typename std::enable_if<!Y::supports_contingency, int>::type>
void BaseBatchSweep<YbusPolicy, SbusPolicy, INIT>::_compute_threaded(
    size_t nb_steps, const CplxVect & Vinit_solver, bool ac_solver_used,
    int max_iter, real_type tol, real_type sn_mva, double & timer_thread_init)
{
    const real_type tol_ = tol / sn_mva;
    const int nb_thread = std::min(static_cast<int>(nb_steps), std::max(1, _nb_thread));
    // DC theta-only fast path (see BaseAlgo::set_lazy_v): _handle_disconnected_grid
    // never applies on this (!Y::supports_contingency) instantiation.
    const bool use_dc_lazy_v = !ac_solver_used;

    if(nb_thread <= 1){
        // single-threaded path: reuse the (already warmed up) member solver and
        // the member accumulators -> identical to the legacy code.
        _algo.set_lazy_v(use_dc_lazy_v);
        int step_diverge = -1;
        std::exception_ptr err;
        _run_range(0, nb_steps, _algo, _algo_controler, ac_cache_.mat, Vinit_solver,
                  ac_solver_used, max_iter, tol_, _nb_solved, _nb_converged, _timer_solver, _timer_modif_Ybus,
                  step_diverge, err, false);
        if(err) std::rethrow_exception(err);
        _status = step_diverge < 0 ? 1 : 0;
        return;
    }

    // multi-threaded path: every thread owns its solver and its solver control.
    // The admittance matrix is NOT copied per thread: the topology is fixed here,
    // so ac_cache_.mat/dc_cache_.mat stay read-only for the whole loop and are simply shared.
    auto timer_thread = CustTimer();
    std::vector<std::unique_ptr<AlgorithmSelector> > algos(nb_thread);
    std::vector<AlgoControl> controls(nb_thread);
    std::vector<int> th_nb_solved(nb_thread, 0);
    std::vector<int> th_nb_converged(nb_thread, 0);
    std::vector<double> th_timer_solver(nb_thread, 0.);
    std::vector<double> th_timer_modif_ybus(nb_thread, 0.);
    std::vector<int> th_diverge(nb_thread, -1);
    std::vector<std::exception_ptr> th_err(nb_thread);

    auto init_thread = [&](int t){
        algos[t] = make_thread_algo();
        algos[t]->set_lazy_v(use_dc_lazy_v);
        controls[t] = _algo_controler;
    };

    std::vector<std::thread> threads;
    threads.reserve(nb_thread - 1);
    for(int t = 1; t < nb_thread; ++t){
        size_t b, e;
        split_range(nb_steps, nb_thread, t, b, e);
        threads.emplace_back([this, t, b, e, ac_solver_used, max_iter, tol_,
                              &init_thread, &algos, &controls, &th_nb_solved, &th_nb_converged,
                              &th_timer_solver, &th_timer_modif_ybus, &th_diverge, &th_err, &Vinit_solver](){
            init_thread(t);
            _run_range(b, e, *algos[t], controls[t], ac_cache_.mat, Vinit_solver, ac_solver_used, max_iter, tol_,
                      th_nb_solved[t], th_nb_converged[t], th_timer_solver[t], th_timer_modif_ybus[t], th_diverge[t], th_err[t], true);
        });
    }
    timer_thread_init = timer_thread.duration();

    {
        size_t b, e;
        split_range(nb_steps, nb_thread, 0, b, e);
        init_thread(0);
        _run_range(b, e, *algos[0], controls[0], ac_cache_.mat, Vinit_solver, ac_solver_used, max_iter, tol_,
                  th_nb_solved[0], th_nb_converged[0], th_timer_solver[0], th_timer_modif_ybus[0], th_diverge[0], th_err[0], true);
    }

    for(auto & th : threads) th.join();

    for(int t = 0; t < nb_thread; ++t){
        if(th_err[t]) std::rethrow_exception(th_err[t]);
    }
    bool all_conv = true;
    for(int t = 0; t < nb_thread; ++t){
        _nb_solved += th_nb_solved[t];
        _nb_converged += th_nb_converged[t];
        _timer_solver += th_timer_solver[t];
        _timer_modif_Ybus += th_timer_modif_ybus[t];
        if(th_diverge[t] >= 0) all_conv = false;
    }
    _status = all_conv ? 1 : 0;
}

template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
template<class Y, typename std::enable_if<Y::supports_contingency, int>::type>
void BaseBatchSweep<YbusPolicy, SbusPolicy, INIT>::_compute_threaded(
    size_t nb_steps, const CplxVect & Vinit_solver, bool ac_solver_used,
    int max_iter, real_type tol, real_type sn_mva, double & timer_thread_init)
{
    const real_type tol_ = tol / sn_mva;
    const bool mask_mode = _handle_disconnected_grid;
    // DC theta-only fast path (see BaseAlgo::set_lazy_v): never on together with
    // mask_mode -- that path stays on the legacy, always-eager _voltages
    // accumulation (see _run_range_masked).
    const bool use_dc_lazy_v = !ac_solver_used && !mask_mode;

    if(std::min(static_cast<int>(nb_steps), std::max(1, _nb_thread)) <= 1){
        // single-threaded path: reuse the (already warmed-up) member solver, member
        // ac_cache_.mat and the member accumulators -> identical to the legacy code.
        _algo.set_lazy_v(use_dc_lazy_v);
        std::exception_ptr err;
        int step_diverge = -1;
        if(mask_mode){
            _maybe_run_range_masked(0, nb_steps, _algo, _algo_controler, ac_cache_.mat, Vinit_solver,
                                    ac_solver_used, max_iter, tol, sn_mva,
                                    _timer_modif_Ybus, _nb_solved, _nb_converged, _timer_solver, step_diverge, err, false);
        } else {
            _run_range(0, nb_steps, _algo, _algo_controler, ac_cache_.mat, Vinit_solver, ac_solver_used,
                      max_iter, tol_, _nb_solved, _nb_converged, _timer_solver, _timer_modif_Ybus, step_diverge, err, false);
        }
        if(err) std::rethrow_exception(err);
        // meaningful only where SbusPolicy::supports_vary (ContingencyAnalysis has no
        // aggregate _status / get_status() -- each row is reported independently via
        // converged()/get_violations() instead); harmless write otherwise.
        if(SbusPolicy::supports_vary) _status = step_diverge < 0 ? 1 : 0;
        return;
    }

    // multi-threaded path: every thread owns its solver / control / Ybus copy
    // (unlike the !Y::supports_contingency overload above: emulating a
    // disconnection means editing Ybus, so it cannot be shared read-only).
    const int nb_thread = std::min(static_cast<int>(nb_steps), std::max(1, _nb_thread));
    std::vector<std::unique_ptr<AlgorithmSelector> > algos(nb_thread);
    std::vector<AlgoControl> controls(nb_thread);
    std::vector<Eigen::SparseMatrix<cplx_type> > ybus_copies(nb_thread);
    std::vector<double> th_timer_modif(nb_thread, 0.);
    std::vector<int> th_nb_solved(nb_thread, 0);
    std::vector<int> th_nb_converged(nb_thread, 0);
    std::vector<double> th_timer_solver(nb_thread, 0.);
    std::vector<int> th_diverge(nb_thread, -1);
    std::vector<std::exception_ptr> th_err(nb_thread);

    auto init_thread = [&](int t){
        algos[t] = make_thread_algo();
        algos[t]->set_lazy_v(use_dc_lazy_v);
        controls[t] = _algo_controler;
        ybus_copies[t] = ac_cache_.mat;
    };

    auto timer_thread = CustTimer();
    std::vector<std::thread> threads;
    threads.reserve(nb_thread - 1);
    for(int t = 1; t < nb_thread; ++t){
        size_t b, e;
        split_range(nb_steps, nb_thread, t, b, e);
        threads.emplace_back([this, t, b, e, ac_solver_used, max_iter, tol, tol_, sn_mva, mask_mode,
                              &init_thread, &algos, &controls, &ybus_copies, &th_timer_modif,
                              &th_nb_solved, &th_nb_converged, &th_timer_solver, &th_diverge, &th_err, &Vinit_solver](){
            init_thread(t);
            if(mask_mode){
                _maybe_run_range_masked(b, e, *algos[t], controls[t], ybus_copies[t], Vinit_solver, ac_solver_used,
                                        max_iter, tol, sn_mva, th_timer_modif[t], th_nb_solved[t], th_nb_converged[t], th_timer_solver[t],
                                        th_diverge[t], th_err[t], true);
            } else {
                _run_range(b, e, *algos[t], controls[t], ybus_copies[t], Vinit_solver, ac_solver_used, max_iter, tol_,
                          th_nb_solved[t], th_nb_converged[t], th_timer_solver[t], th_timer_modif[t], th_diverge[t], th_err[t], true);
            }
        });
    }
    timer_thread_init = timer_thread.duration();

    {
        size_t b, e;
        split_range(nb_steps, nb_thread, 0, b, e);
        init_thread(0);
        if(mask_mode){
            _maybe_run_range_masked(b, e, *algos[0], controls[0], ybus_copies[0], Vinit_solver, ac_solver_used,
                                    max_iter, tol, sn_mva, th_timer_modif[0], th_nb_solved[0], th_nb_converged[0], th_timer_solver[0],
                                    th_diverge[0], th_err[0], true);
        } else {
            _run_range(b, e, *algos[0], controls[0], ybus_copies[0], Vinit_solver, ac_solver_used, max_iter, tol_,
                      th_nb_solved[0], th_nb_converged[0], th_timer_solver[0], th_timer_modif[0], th_diverge[0], th_err[0], true);
        }
    }

    for(auto & th : threads) th.join();

    for(int t = 0; t < nb_thread; ++t){
        if(th_err[t]) std::rethrow_exception(th_err[t]);
    }
    bool all_conv = true;
    for(int t = 0; t < nb_thread; ++t){
        _nb_solved += th_nb_solved[t];
        _nb_converged += th_nb_converged[t];
        _timer_solver += th_timer_solver[t];
        _timer_modif_Ybus += th_timer_modif[t];
        if(th_diverge[t] >= 0) all_conv = false;
    }
    if(SbusPolicy::supports_vary) _status = all_conv ? 1 : 0;
}

template<class YbusPolicy, class SbusPolicy, BatchInitKind INIT>
void BaseBatchSweep<YbusPolicy, SbusPolicy, INIT>::compute(
    const Eigen::Ref<const CplxVect> & Vinit, int max_iter, real_type tol)
{
    auto timer = CustTimer();
    auto timer_preproc = CustTimer();

    // perform some initial checks and reset timers
    size_t nb_total_bus = _reset_data_and_check_vinit(Vinit);
    _status = 0;
    _timer_modif_Ybus = 0.;
    _timer_thread_init = 0.;

    const auto & sn_mva = _grid_model.get_sn_mva();
    const bool ac_solver_used = _algo.ac_solver_used();

    const size_t nb_steps = _nb_steps();

    // per-row converged mask (converged_mask()): every instantiation, always on,
    // regardless of compute_limit_violations -- see BaseBatchSweep.hpp's comment
    // next to it. Rows never reached this compute() (eg a TimeSeries chain that
    // aborts, or the diverging-"n"-case early return below) stay 0.
    _converged_mask_.assign(nb_steps, 0);

    // limit-violation bookkeeping (ContingencyAnalysis only; no-op elsewhere)
    _refresh_defaults_vect_cache();
    if(_compute_limit_violations_){
        _converged.assign(nb_steps, 0);
        _violations.assign(nb_steps, std::vector<LimitViolation>());
        _converged_n_ = false;
        _violations_n_.clear();
    }

    // prepare the gridmodel (compute Ybus, Sbus etc.)
    CplxVect Vinit_solver = prepare_solver_input_base(Vinit, ac_solver_used);

    // initialize whatever varies (Ybus and/or Sbus -- each a no-op where the
    // corresponding policy is NOOP)
    _prepare_ybus_varying(ac_solver_used, static_cast<Eigen::Index>(nb_steps));
    _prepare_sbus_varying(ac_solver_used, static_cast<Eigen::Index>(nb_steps));

    // "handle disconnected grid" mode pre-pass (ContingencyAnalysis only; no-op
    // elsewhere -- and _handle_disconnected_grid can never be true elsewhere, since
    // no setter exists to set it there)
    _maybe_prepare_masks();

    // DC theta-only fast path (see BaseAlgo::set_lazy_v): every DC compute() except
    // the "handle disconnected grid" masked one (which stays on the always-eager
    // legacy path -- see _run_range_masked). _sbus_gen_v() is the generic,
    // SbusPolicy-agnostic copy of sbus_policy_.gen_v (BaseBatchSolverSynch's
    // magnitude-reconstruction helpers do not need to know about SbusPolicy at all).
    const bool use_dc_lazy_v = !ac_solver_used && !_handle_disconnected_grid;
    if(use_dc_lazy_v) _dc_gen_v_ = _sbus_gen_v();

    bool n_powerflow_has_conv = _finish_preprocessing(
        nb_steps, nb_total_bus, Vinit_solver, max_iter, tol, timer_preproc, use_dc_lazy_v
    );

    if(!n_powerflow_has_conv){
        // a diverging base case: a Vary-carrying instantiation (TimeSeries /
        // InjectionSweep / ScenarioSweep) reports it via the aggregate status --
        // _status is set to 0 and compute() returns normally, so get_status() (and
        // the legacy compute_Vs wrapper, which just forwards it) reads 0, not the
        // pre-refactor BaseInjectionSweep::compute_Vs's historical -1 (every existing
        // caller only ever tests `!= 1`, so this is a difference in the exact value,
        // not in whether the failure is reported). ContingencyAnalysis instead
        // throws: every contingency is solved relative to this base case, so a
        // diverging one makes the whole analysis meaningless. Both branches below are
        // plain, generic code (no Y/S-specific member access), so this runtime `if`
        // on a compile-time constant is safe to compile for every instantiation.
        if(SbusPolicy::supports_vary){
            _status = 0;
            // With compute_limit_violations set, a diverging "n" powerflow must NOT
            // leave get_violations()/get_violations_n() looking like empty lists --
            // that reads as "converged, nothing found" (ScenarioSweep has no
            // converged()/converged_n() escape hatch; a violation list is meant to
            // fully encode row status by itself, see BaseBatchSweep.hpp's note next
            // to get_compute_limit_violations()). Stamp the same GRID/DIVERGENCE
            // sentinel _store_row_status() already uses for a per-row divergence.
            if(_compute_limit_violations_){
                const LimitViolation sentinel{
                    ViolationElementType::GRID, -1, 0, LimitViolationType::DIVERGENCE,
                    std::numeric_limits<real_type>::quiet_NaN(),
                    std::numeric_limits<real_type>::quiet_NaN()};
                _violations_n_.push_back(sentinel);
                for(auto & row : _violations) row.push_back(sentinel);
            }
            _results_stale_ = false;
            _timer_total = timer.duration();
            return;
        }
        std::ostringstream exc_;
        exc_ << algo_name() << "::compute: the pre-contingency (\"n\", no disconnection) "
                "powerflow did not converge (error: " << _algo.get_error() << "). No contingency "
                "can be simulated from a non-converged base case; check the input `Vinit` / "
                "`tol` / `max_iter`, or the grid state itself.";
        throw std::runtime_error(exc_.str());
    }

    // pre-contingency ("n") case violations (ContingencyAnalysis only; no-op
    // elsewhere)
    _record_n_case_violations(_algo.get_V());

    // compute the powerflows, possibly split across several threads
    _compute_threaded(nb_steps, Vinit_solver, ac_solver_used, max_iter, tol, sn_mva, _timer_thread_init);

    _results_stale_ = false;
    _timer_total = timer.duration();
}

// Compile all 4 instantiations once, here, into the core library: every other
// translation unit picks them up through the `extern template` declarations at the
// bottom of BaseBatchSweep.hpp instead of instantiating the template again. This
// only instantiates the class's NON-template members (compute/_run_range/
// _run_one_step) plus, transitively, whichever member templates those bodies
// themselves call (eg _compute_threaded, _prepare_ybus_varying, ...) -- exactly the
// mechanism that makes SFINAE-gating the *other*, Python-bound member templates
// necessary to define inline in BaseBatchSweep.hpp instead (see the class-level
// comment there).
template class LS2G_API BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromPreviousStep>;  // TimeSeries
template class LS2G_API BaseBatchSweep<YbusPolicy::NOOP, SbusPolicy::Vary, BatchInitKind::FromSeed>;          // InjectionSweep
template class LS2G_API BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::NOOP, BatchInitKind::FromSeed>;   // ContingencyAnalysis
template class LS2G_API BaseBatchSweep<YbusPolicy::Contingency, SbusPolicy::Vary, BatchInitKind::FromSeed>;   // ScenarioSweep

} // namespace ls2g
