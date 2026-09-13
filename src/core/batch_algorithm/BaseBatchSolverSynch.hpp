// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEMULTIPLEPOWERFLOW_H
#define BASEMULTIPLEPOWERFLOW_H

#include "LSGrid.hpp"

#include <memory>

namespace ls2g {

/**
This is a utility class, used for TimeSeries and SecurityAnalysis that abstract some computations when
the same solver is re used multiple times.

It allows to perform "batch" powerflow one a time in a synchronous manner.

The "solver" of the gridmodel is never really used to perform powerflows.

**/
class LS2G_API BaseBatchSolverSynch : protected BaseConstants
{
    public:
        using RealMat = Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using CplxMat =  Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        
        explicit BaseBatchSolverSynch(const LSGrid & init_grid_model):
            _grid_model(init_grid_model),
            n_line_(init_grid_model.nb_powerline()),
            n_trafos_(init_grid_model.nb_trafo()),
            n_total_(n_line_ + n_trafos_),
            _algo(),
            _voltages(),
            _amps_flows(),
            _active_power_flows()
            {
                // give the internal solvers access to the (copied) grid: the
                // hvdc angle-droop data flows through this pointer (the member
                // _grid_model has a stable address, the class is not copyable)
                _algo.set_lsgrid(&_grid_model);
                // inherit the source grid's AC algorithm (type + config, eg a
                // ScalingPolicyType::MaxVoltageChange damping set via
                // set_ac_algo_config): _algo is otherwise a fresh, independent,
                // default (undamped NR_SparseLU) solver, so without this the n /
                // n-1 powerflows run here can converge to a different root than
                // (or fail to converge unlike) init_grid_model.ac_pf() -- same class
                // of bug as the one just fixed in LSGrid's copy constructor.
                // Inherit by *name*, not AlgorithmType: the enum collapses every
                // plugin (and every built-in with no dedicated enum member, eg
                // NRRefactorRetry_*) onto AlgorithmType::Custom, for which
                // AlgorithmSelector::change_algorithm(AlgorithmType) unconditionally
                // throws -- so a grid using such a solver could never be handed to
                // TimeSeries / ContingencyAnalysis / SecurityAnalysis at all.
                _algo.change_algorithm(init_grid_model.get_algo().get_name());
                _algo.set_config(init_grid_model.get_ac_algo_config());
            }
        virtual ~BaseBatchSolverSynch() noexcept = default;  // to avoid warning about overload virtual
        BaseBatchSolverSynch(const BaseBatchSolverSynch&) = delete;
        BaseBatchSolverSynch(BaseBatchSolverSynch&&) = delete;
        BaseBatchSolverSynch & operator=(BaseBatchSolverSynch&&) = delete;
        BaseBatchSolverSynch & operator=(const BaseBatchSolverSynch&) = delete;
    
        bool get_init_from_n_powerflow() const noexcept {return _init_from_n_powerflow;}
        void set_init_from_n_powerflow(bool do_it) {
            if(do_it == _init_from_n_powerflow) return;
            // L2: this decides what every row starts from -- the "n" solve's answer or
            // the caller's own Vinit.
            clear_batch_inputs();
            _init_from_n_powerflow = do_it;
        }

        // Whether the elements of this batch may be solved concurrently. This is a
        // property of the algorithm, not a user setting: it is false exactly when one
        // element is initialized with the result of the element before it, because
        // splitting such a batch into per-thread ranges would break that chain and make
        // the results depend on the number of threads. Algorithms whose elements all
        // start from the same seed (ContingencyAnalysis, InjectionSweep) are unaffected
        // and return the default.
        static constexpr bool _default_supports_multithread = true;
        virtual bool supports_multithread() const {return _default_supports_multithread;}

        // Number of OS threads used to solve the batch (default: 1). With nb_thread == 1
        // the behaviour is identical to the legacy sequential path. With nb_thread > 1
        // the batch is split into contiguous ranges, each solved by its own thread (own
        // solver), writing to disjoint rows of the result matrix -- so the results do not
        // depend on the number of threads. Values < 1 are clamped to 1.
        int get_nb_thread() const {return _nb_thread;}
        void set_nb_thread(int n) {
            // rejected here rather than at compute time: an algorithm that chains its
            // elements can never become parallelisable later, so there is no state in
            // which a stored value > 1 could turn out to be legal.
            if(n > 1 && !supports_multithread()){
                std::ostringstream exc_;
                exc_ << "BaseBatchSolverSynch::set_nb_thread: this algorithm initializes each "
                        "computation with the result of the previous one, so its computations "
                        "cannot be split over several threads (got nb_thread = " << n << "). Use "
                        "InjectionSweep instead: it computes the very same injections, but starts "
                        "every computation from the same voltage, which makes them independent "
                        "of one another (and of their ordering).";
                throw std::runtime_error(exc_.str());
            }
            const int val = (n < 1 ? 1 : n);
            if(val == _nb_thread) return;
            // L2: the batch is split differently, and each worker gets its own solver
            // warmed up from the base case -- so the base case is rebuilt with it.
            clear_batch_inputs();
            _nb_thread = val;
        }

        // solver "control"
        // L1, on purpose: a different algorithm may be a different FAMILY (AC <-> DC),
        // and then not one of the three levels survives -- the matrix, the injections
        // and the bus labelling are all the other family's. Telling the two cases apart
        // would buy a cache hit on a call nobody makes in a loop, so this stays the
        // blunt, safe drop.
        virtual void change_algorithm(const AlgorithmType & type){
            _algo.change_algorithm(type);
            this->clear_grid_results();
        }
        // String-based overload: looks up the solver by registry name, so it
        // works for plugin solvers (and built-ins with no dedicated AlgorithmType
        // member, eg NRRefactorRetry_*) which change_algorithm(AlgorithmType)
        // cannot reach (AlgorithmType::Custom is rejected there).
        virtual void change_algorithm(const std::string & name){
            _algo.change_algorithm(name);
            this->clear_grid_results();
        }

        // Returns the enum-typed solvers available in this build.
        // Does not include plugin (Custom) solvers; use AlgorithmRegistry::available_default_algorithms()
        // for the full list.
        std::vector<AlgorithmType> available_default_algorithms() const {return _algo.available_default_algorithms(); }

        // Returns all solver names currently registered (built-in + plugins).
        std::vector<std::string> available_algorithm_names() const {
            return AlgorithmRegistry::instance().available_algorithm_names();
        }

        AlgorithmType get_algo_type() const {return _algo.get_type(); }
        // Registry name of the currently selected algorithm. Unlike get_algo_type(),
        // stays meaningful for plugin / no-enum-member solvers (see AlgorithmSelector::get_name()).
        const std::string & get_algo_name() const {return _algo.get_name(); }

        // the config (eg ScalingPolicyType / damping parameters) of the internal
        // solver used for every n / n-1 powerflow here. This is a fresh, independent
        // solver (see the constructor): the source LSGrid's own config is copied in
        // once at construction time, but is NOT kept in sync afterwards, so re-apply
        // it here if you change it on the grid model after building this object (or
        // after a change_algorithm(), which resets to that new algorithm's defaults).
        AlgoConfig get_algo_config() const {return _algo.get_config(); }
        // L1, for the same reason as change_algorithm above: a config change (damping,
        // scaling policy, ...) changes what the "n" solve converges to, and possibly
        // whether it converges at all. Not a call anyone makes per batch, so it is not
        // worth reasoning about which levels a given field really touches.
        void set_algo_config(const AlgoConfig & cfg) {
            _algo.set_config(cfg);
            this->clear_grid_results();
        }

        // // TODO
        // void change_gridmodel(const GridModel & new_grid_model){
        //     CplxVect V = CplxVect::Constant(new_grid_model.total_bus(), 1.04);
        //     _grid_model = new_grid_model.copy(); // not implemented !!!
        // }

        // utlities informations
        double amps_computation_time() const {return _timer_compute_A;}
        double solver_time() const {return _timer_solver;}
        // time spent building the per-thread solvers (0. when nb_thread == 1)
        double thread_init_time() const {return _timer_thread_init;}
        int nb_solved() const {return _nb_solved;}
        // number of those nb_solved() attempts that actually converged (ie a row
        // that was invertible AND whose solver reported convergence). Always
        // <= nb_solved(): a row skipped outright (eg a non-invertible / islanding
        // Ybus) never reaches compute_one_powerflow, so it counts towards neither.
        int nb_converged() const {return _nb_converged;}
        /**
         * ================= the three cache levels =================================
         *
         * A batch holds three tiers of state, each built from the one above it, and
         * each kept across compute() calls so a caller in a loop does not pay for it
         * again (see BaseBatchSweep::set_reuse_base_case). Dropping a level drops
         * every level below it -- that nesting IS the contract, and it is why these
         * three are the only entry points: a modifier names the highest level it
         * invalidates and never has to enumerate the rest.
         *
         *   L1  clear_grid_results()   what was read off the GRID: Ybus / Bbus, the
         *                              injections, the grid<->solver bus labelling,
         *                              the pv/pq split, the slack and its weights.
         *                              Only a different grid invalidates this (see
         *                              the change_gridmodel TODO below), or the
         *                              caller opting out of reuse entirely.
         *
         *   L2  clear_batch_inputs()   what was built for THIS BATCH from that grid:
         *                              the per-row Ybus edit lists, what the graph
         *                              walk settled about each contingency, the
         *                              masks, the reserved PV<->PQ structure, the
         *                              "n" powerflow that seeds the rows -- and the
         *                              algorithm itself, whose ledger, Jacobian
         *                              sparsity and factorization are exactly what
         *                              that "n" solve built.
         *
         *   L3  clear_batch_outputs()  what the last run PRODUCED: voltages, flows,
         *                              convergence, limit violations, the kept
         *                              Jacobians -- plus the counters and timers
         *                              that describe that run.
         *
         * Each is idempotent by design: dropping a level that is already invalid
         * costs nothing but the (cheap) recursion downwards. That matters -- the
         * python wrapper registers an N-1 sweep one contingency at a time, so L2
         * runs once per line, and only the first of them has anything to do.
         */

        // L1 -- see above.
        virtual void clear_grid_results() {
            if(_grid_cache_valid_){
                _grid_cache_valid_ = false;
                ac_cache_.clear();
                dc_cache_.clear();
                nb_buses_solver_ = -1;
            }
            clear_batch_inputs();
        }

        // L2 -- see above. The algorithm goes with it: an analysis and a
        // factorization are not results, they are what the "n" solve built and what
        // the next batch would otherwise redo. Dropping one without the other is the
        // one combination that cannot be made safe -- a kept base case solving
        // against a freshly reset (so default-constructed) system is a segfault, not
        // a wrong answer.
        virtual void clear_batch_inputs() {
            if(_batch_inputs_valid_){
                _batch_inputs_valid_ = false;
                _algo.reset();
                _base_case_v_solver_ = CplxVect();
                _prepared_nb_steps_ = _nb_steps_none;
            }
            clear_batch_outputs();
        }

        // L3 -- see above. THE RESULTS, and nothing else: the voltages and the two
        // flow matrices derived from them, plus the counters and the timers of the
        // run that produced them (compute() zeroes those again as the next one
        // starts). Not the algorithm, not the caches.
        virtual void clear_batch_outputs() {
            _voltages = CplxMat();
            _amps_flows = RealMat();
            _active_power_flows = RealMat();
            // the DC fast path holds the same voltages as an angle plus a shared
            // magnitude (see _dc_vm_row_grid): one result in several pieces, so the
            // pieces go together -- a stale one would be reconstructed into a voltage
            // that was never solved
            _thetas = RealMat();
            _dc_row_solved_.clear();
            _dc_lazy_storage_used_ = false;
            _dc_base_vm_solver_ = RealVect();
            _dc_base_vm_grid_ = RealVect();
            _dc_gen_v_ = RealMat();
            _nb_solved = 0;
            _nb_converged = 0;
            _timer_compute_A = 0.;
            _timer_compute_P = 0.;
            _timer_solver = 0.;
            _timer_thread_init = 0.;
            _timer_total = 0.;
            _timer_pre_proc = 0.;
            _thread_solver_stats_.clear();
        }

        // Drop everything: the three levels above, and (in BaseBatchSweep) the
        // registrations and settings on top of them -- a different axis, which is why
        // this is not a fourth level.
        virtual void clear() {
            clear_grid_results();
            // NB: _nb_thread is deliberately NOT reset -- it is a setting, not a result.
        }
    public:

        // field-wise +=, so get_linear_solver_stats() can report the whole compute()
        // rather than whichever algorithm happened to be the member one.
        static void _accumulate_solver_stats(LinearSolverStats & into, const LinearSolverStats & from){
            into.nb_reset += from.nb_reset;
            into.nb_analyze += from.nb_analyze;
            into.nb_factorize += from.nb_factorize;
            into.nb_refactorize += from.nb_refactorize;
            into.nb_refactorize_failed += from.nb_refactorize_failed;
            into.nb_fallback_factorize += from.nb_fallback_factorize;
            into.nb_fallback_factorize_failed += from.nb_fallback_factorize_failed;
            into.nb_solve += from.nb_solve;
            into.timer_initialize_ += from.timer_initialize_;
            into.timer_factor_ += from.timer_factor_;
            into.timer_refactor_ += from.timer_refactor_;
            into.timer_solve_ += from.timer_solve_;
        }

        /**
         * The linear-solver counters of the whole compute(): how many symbolic analyses,
         * numeric factorizations and refactorizations it took, summed over the member
         * algorithm AND every worker thread's own (see _thread_solver_stats_).
         *
         * What a batch is FOR is that the expensive half stays flat however many rows it
         * runs: the Jacobian's sparsity pattern is fixed for the whole batch, so each
         * algorithm analyzes ONCE and every row after its first only refactorizes.
         * Counted per algorithm rather than per batch, because that is what the work
         * actually is: single-threaded that means nb_analyze == 1; with nb_thread > 1 it
         * means **one analyze per algorithm used**, ie one per worker plus the member
         * one that solved the "n" warm-up case. Reading more than that is a row changing
         * the sparsity pattern, which is a bug, not a tuning question.
         */
        LinearSolverStats get_linear_solver_stats() const {
            LinearSolverStats res = _algo.get_linear_solver_stats();
            for(const auto & th : _thread_solver_stats_) _accumulate_solver_stats(res, th);
            return res;
        }

        /**
         * The same counters, kept apart: index 0 is the member algorithm (which solves
         * the "n" warm-up case, and the whole batch when single-threaded), then one entry
         * per worker thread of the last multi-threaded compute(). Empty of worker entries
         * after a single-threaded run. This is what to look at to tell "every algorithm
         * analyzed once" from "one algorithm analyzed several times" -- the sum alone
         * cannot distinguish them.
         */
        std::vector<LinearSolverStats> get_linear_solver_stats_per_algo() const {
            std::vector<LinearSolverStats> res;
            res.reserve(1 + _thread_solver_stats_.size());
            res.push_back(_algo.get_linear_solver_stats());
            for(const auto & th : _thread_solver_stats_) res.push_back(th);
            return res;
        }

        // results
        // this should not be const, see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html#pass-by-reference
        // tl;dr: const can make copies ! OR NOT I AM LOST
        const RealMat & get_flows() const {return _amps_flows;}
        const RealMat & get_power_flows() const {return _active_power_flows;}
        // DC theta-only fast path (see BaseAlgo::set_lazy_v): when a compute() used
        // it, _voltages is left empty and _thetas (+ the small _dc_* inputs) holds
        // everything needed to rebuild it -- done here, lazily, on first request, and
        // cached (_voltages stays mutable for exactly this reason). A compute() that
        // did not use the fast path (AC, or DC with handle_disconnected_grid) already
        // leaves _voltages fully populated, so this is then a no-op.
        const CplxMat & get_voltages() const {
            if(_dc_lazy_storage_used_ && _voltages.size() == 0 && _thetas.size() != 0){
                _rebuild_voltages_from_thetas();
            }
            return _voltages;
        }

    protected:
        /**
         * The flows of ONE row through the branches of `structure_data` (lines or
         * trafos): amps at the "from" side when `amps`, active power there otherwise,
         * into `out.row(i)`. The whole batch used to be walked branch by branch, each
         * branch reading a column of the row-major voltage matrix -- a strided pass
         * over every row per branch, two heap temporaries the size of the batch per
         * branch, and a matrix that does not fit in cache re-read once per branch.
         * Row by row, the row's voltages are read once, contiguously, and the row of
         * flows is written once. The arithmetic per branch is the one it was, in the
         * same order, so a flow here is the flow it used to be, bit for bit.
         *
         * Exactly one of `V_row` (the row's complex voltages, grid numbering: AC, or
         * DC without the fast path) and `theta_row` + `vm_row` (the DC fast path: the
         * row's angles and its reconstructed magnitude, grid numbering) is non-null.
         */
        template<class T>
        void _flows_of_row(const T & structure_data,
                           Eigen::Index i,
                           size_t lag_id,
                           bool is_trafo,
                           bool amps,
                           bool is_ac,
                           real_type sn_mva,
                           const cplx_type * V_row,
                           const real_type * theta_row,
                           const RealVect * vm_row,
                           RealMat & out) const
        {
            const auto & bus_vn_kv = _grid_model.get_bus_vn_kv();
            const auto & el_status = structure_data.get_status_global();
            const auto & status1 = structure_data.get_status_side_1();
            const auto & status2 = structure_data.get_status_side_2();
            const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
            const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();

            // AC uses complex (Kron-reduced) coefficients, DC uses real susceptance coefficients
            Eigen::Ref<const CplxVect> vect_yac_ff = structure_data.yac_eff_11();
            Eigen::Ref<const CplxVect> vect_yac_ft = structure_data.yac_eff_12();
            Eigen::Ref<const RealVect> vect_ydc_ff = structure_data.ydc_11();
            Eigen::Ref<const RealVect> vect_ydc_ft = structure_data.ydc_12();
            Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift(); // not used in AC nor if it's powerline anyway

            const size_t nb_el = structure_data.nb();
            const real_type sqrt_3 = sqrt(3.);
            const bool dc_lazy = (V_row == nullptr);

            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                const bool s1 = status1[el_id];
                const bool s2 = status2[el_id];
                const Eigen::Index col = static_cast<Eigen::Index>(el_id + lag_id);

                // a half-open branch (see keep_half_open_lines) has bus_id ==
                // _deactivated_bus_id on its open side and must not be used to index
                // the voltages / bus_vn_kv -- an open-end voltage of exactly 0 is
                // substituted (both sides treated the same way). For AC, yac_eff_* is
                // already Kron-reduced for whichever side is open, so this alone gives
                // the correct "or"-side (side 1) flow either way; DC has no such
                // reduction (handled explicitly below).
                const int from_me = bus_from(el_id).cast_int();
                const int to_me = bus_to(el_id).cast_int();
                // vn_kv base for the amps conversion: whichever side is actually
                // energized. If side 1 (the one being measured) is open the numerator
                // is exactly 0 regardless, so the base only has to avoid a 0/0.
                const real_type bus_vn_kv_f = s1 ? bus_vn_kv(from_me) : (s2 ? bus_vn_kv(to_me) : real_type(1.));

                real_type res;
                if(is_ac){
                    const cplx_type Efrom = s1 ? V_row[from_me] : cplx_type(0., 0.);
                    const cplx_type Eto = s2 ? V_row[to_me] : cplx_type(0., 0.);
                    const cplx_type y_ff = vect_yac_ff(el_id);
                    const cplx_type y_ft = vect_yac_ft(el_id);
                    // trafo equations (to get the power at the "from" side)
                    cplx_type I_ft = y_ff * Efrom + y_ft * Eto;
                    I_ft = std::conj(I_ft);
                    const cplx_type S_ft = Efrom * I_ft;
                    if(amps){
                        const real_type v_f_kv = (s1 ? std::abs(Efrom) : std::abs(Eto)) * bus_vn_kv_f;
                        res = std::abs(S_ft) * sn_mva;
                        res /= sqrt_3 * v_f_kv;
                    } else {
                        res = S_ft.real() * sn_mva;
                    }
                } else {
                    // unlike yac_eff_*, ydc_11/ydc_12 are NOT Kron-reduced for a
                    // half-open branch: DC treats one side open as fully disconnected
                    // (see fillBdc: "disco on one side == disco on both sides"), so
                    // report 0 rather than mixing ydc_ff/ydc_ft with a meaningless
                    // open-end angle.
                    if(!(s1 && s2)){
                        out(i, col) = 0.;
                        continue;
                    }
                    // DC active flow from the bus angles (theta) directly, like the
                    // gridmodel results: P = ydc_ff . theta_from + ydc_ft . theta_to
                    const real_type theta_from = dc_lazy ? theta_row[from_me] : std::arg(V_row[from_me]);
                    const real_type theta_to = dc_lazy ? theta_row[to_me] : std::arg(V_row[to_me]);
                    const real_type y_ff = vect_ydc_ff(el_id);
                    const real_type y_ft = vect_ydc_ft(el_id);
                    res = (y_ff * theta_from + y_ft * theta_to) * sn_mva;
                    if(is_trafo) res -= dc_x_tau_shift(el_id);
                    if(amps){
                        res = std::abs(res);
                        // the magnitude: the row's reconstructed one on the fast path,
                        // |V| otherwise (both sides are closed here)
                        const real_type vm = dc_lazy ? (*vm_row)(from_me) : std::abs(V_row[from_me]);
                        const real_type v_f_kv = vm * bus_vn_kv_f;
                        res /= sqrt_3 * v_f_kv;
                    }
                }
                out(i, col) = res;
            }
        }

        // member version: forwards to the explicit overload below using the
        // member solver / control / accumulators (single-threaded path).
        // V stays a plain reference (not Eigen::Ref): it forwards into both
        // algo.compute_pf(Ybus, V, ...) (AC) and algo.compute_pf_dc(Bbus, V, ...) (DC)
        // depending on ac_solver_used(), and the DC path resizes/reassigns V -- so it
        // can't become Eigen::Ref even though the AC path alone would allow it.
        bool compute_one_powerflow(const EigenRefConstCplxSpMat     & Ybus,
                                   CplxVect & V,
                                   const Eigen::Ref<const CplxVect> & Sbus,
                                   const Eigen::Ref<const IntVect> & slack_ids,
                                   const Eigen::Ref<const RealVect> & slack_weights,
                                   const Eigen::Ref<const IntVect> & bus_pv,
                                   const Eigen::Ref<const IntVect> & bus_pq,
                                   int max_iter,
                                   double tol
                                   );

        // explicit version: operates on the passed solver / control and writes
        // its book-keeping into the passed accumulators. This is what the
        // multi-threaded ContingencyAnalysis uses (one solver per thread). The
        // read-only member dc_cache_.mat is only read here (safe to share across threads).
        // V stays plain -- same AC/DC forwarding reason as the member overload above.
        bool compute_one_powerflow(AlgorithmSelector & algo,
                                   AlgoControl & control,
                                   int & nb_solved,
                                   int & nb_converged,
                                   double & timer_solver,
                                   const EigenRefConstCplxSpMat     &  Ybus,
                                   CplxVect & V,
                                   const Eigen::Ref<const CplxVect> & Sbus,
                                   const Eigen::Ref<const IntVect> & slack_ids,
                                   const Eigen::Ref<const RealVect> & slack_weights,
                                   const Eigen::Ref<const IntVect> & bus_pv,
                                   const Eigen::Ref<const IntVect> & bus_pq,
                                   int max_iter,
                                   double tol
                                   );

        // Warm up a (freshly built) solver so its symbolic factorization / DC
        // internal Ybus / sparsity pattern match the member solver after the
        // n-powerflow. Mirrors the n-powerflow block of _finish_preprocessing
        // (works on the member ac_cache_.mat / dc_cache_.mat / ac_cache_.inj / dc_cache_.inj — all read-only).
        // `control` is left in the "nothing changed" state on return so the
        // subsequent per-contingency solves reuse the factorization.
        bool warmup_solver(AlgorithmSelector & algo,
                           AlgoControl & control,
                           const Eigen::Ref<const CplxVect> & Vinit_solver,
                           int max_iter,
                           real_type tol);

        void compute_flows_from_Vs(bool amps=true);

        // ----- multi-threading helpers (shared by every batch algorithm) ------------
        // Build a solver for one worker thread. Deliberately a FRESH AlgorithmSelector
        // rather than a copy of `_algo`: only the algorithm identity and its config are
        // inherited, never the member solver's (already warmed up) internal state.
        // Selects by *name*, not by AlgorithmType: a plugin -- or a built-in with no
        // dedicated enum member, eg NRRefactorRetry_* -- reports AlgorithmType::Custom,
        // which change_algorithm(AlgorithmType) unconditionally rejects.
        // Only reads shared state (_grid_model, _algo's name / config), so several
        // threads may call it concurrently.
        std::unique_ptr<AlgorithmSelector> make_thread_algo() const;

        // Contiguous split of [0, nb_items) over `nb_thread` threads: the first
        // `nb_items % nb_thread` threads get one extra item. Writes thread `t`'s share
        // into [begin, end).
        static void split_range(size_t nb_items, int nb_thread, int t,
                                size_t & begin, size_t & end);

        CplxVect extract_Vsolver_from_Vinit(const Eigen::Ref<const CplxVect> & Vinit,
                                            size_t nb_buses_solver,
                                            size_t nb_total_bus,
                                            const SolverBusIdVect & id_me_to_ac_solver){
            // extract V solver from the given V
            CplxVect Vinit_solver = CplxVect::Constant(nb_buses_solver, {_grid_model.get_init_vm_pu(), 0.});
            int tmp;
            for(size_t bus_id_grid = 0; bus_id_grid < nb_total_bus; ++bus_id_grid){
                tmp = static_cast<int>(id_me_to_ac_solver[bus_id_grid]);
                if(tmp == BaseConstants::_deactivated_bus_id) continue;
                Vinit_solver[tmp] = Vinit[bus_id_grid];
            }
            return Vinit_solver;
        }
    protected:

        /**
         * The family-agnostic half of whichever cache this sweep is building into.
         *
         * A batch runs AC or DC for its whole life, so exactly one of ac_cache_ /
         * dc_cache_ is ever populated. Most of this class does not care which --
         * the bus labelling, the slack and the pv-pq split have the same type
         * either way -- so it reads them through here rather than branching at
         * every use. `_ac_solver_used` is set by prepare_solver_input_base, and
         * defaults to the AC side so that a read before any prep sees an empty
         * layout rather than a dangling one.
         */
        [[nodiscard]] SolverBusLayout & active_layout() noexcept {
            return _ac_solver_used ? static_cast<SolverBusLayout &>(ac_cache_)
                                   : static_cast<SolverBusLayout &>(dc_cache_);
        }
        [[nodiscard]] const SolverBusLayout & active_layout() const noexcept {
            return _ac_solver_used ? static_cast<const SolverBusLayout &>(ac_cache_)
                                   : static_cast<const SolverBusLayout &>(dc_cache_);
        }

        CplxVect prepare_solver_input_base(const Eigen::Ref<const CplxVect> & Vinit, bool ac_solver_used){
            // Which family this batch runs, for active_layout(). Fixed for the whole
            // sweep (it comes from the algorithm), but recorded here rather than
            // asked of _algo on every access.
            _ac_solver_used = ac_solver_used;

            // clear previous data: the whole cache, in one call, because it is one
            // object -- there is no way to forget a field
            ac_cache_.clear();
            dc_cache_.clear();
            nb_buses_solver_ = -1;

            // fill the data correctly. One call fills the ENTIRE cache -- labelling,
            // matrix, injections, slack, pv-pq split, slack weights -- into vectors
            // this batch owns. Nothing of it is left in _grid_model to be collected
            // afterwards.
            _algo_controler.tell_all_changed();
            CplxVect res;
            if(ac_solver_used){
                res = _grid_model.build_solver_input(Vinit, ac_cache_, _algo_controler);
                nb_buses_solver_ = static_cast<int>(ac_cache_.mat.cols());
            } else {
                // native real DC path: build the real Bbus / Pbus
                res = _grid_model.build_dc_solver_input(Vinit, dc_cache_, _algo_controler);
                nb_buses_solver_ = static_cast<int>(dc_cache_.mat.cols());
            }
            // L1 is now built: every later compute() may map its own starting voltage
            // onto this labelling instead of reading the grid again (see
            // BaseBatchSweep::_vinit_on_grid_cache).
            _grid_cache_valid_ = true;
            return res;
        }

        size_t _reset_data_and_check_vinit(const Eigen::Ref<const CplxVect> & Vinit){
            const size_t nb_total_bus = _grid_model.total_bus();
            if(static_cast<size_t>(Vinit.size()) != nb_total_bus){
                std::ostringstream exc_;
                exc_ << "TimeSeries::compute_Sbuses: Size of the Vinit should be the same as the total number of buses. Currently:  ";
                exc_ << "Vinit: " << Vinit.size() << " and there are " << nb_total_bus << " buses.";
                exc_ << "(fyi: Components of Vinit corresponding to deactivated bus will be ignored anyway, so you can put whatever you want there).";
                throw std::runtime_error(exc_.str());
            }

            // reset timers
            _nb_solved = 0;
            _nb_converged = 0;
            _timer_pre_proc = 0.;
            _timer_total = 0.;
            _timer_solver = 0.;
            return nb_total_bus;
        }

        // Size this run's result buffers. Purely per-call (L3): a kept base case says
        // nothing about how many rows the NEXT batch has -- it only guarantees the
        // grid and the batch inputs behind them are still the right ones.
        void _size_result_buffers(size_t nb_steps, size_t nb_total_bus, bool use_dc_lazy_v){
            // the DC theta-only fast path accumulates into _thetas (real) instead of
            // _voltages (complex) -- see get_voltages() / _flows_of_row.
            _dc_lazy_storage_used_ = use_dc_lazy_v;
            if(use_dc_lazy_v){
                _thetas = RealMat::Zero(nb_steps, nb_total_bus);
                _dc_row_solved_.assign(nb_steps, 0);
                _voltages = CplxMat();
            } else {
                _voltages = BaseBatchSolverSynch::CplxMat::Zero(nb_steps, nb_total_bus);
                _thetas = RealMat();
                _dc_row_solved_.clear();
            }
            _amps_flows = RealMat::Zero(0, n_total_);
            _active_power_flows = RealMat::Zero(0, n_total_);
        }

        // The "n" powerflow: the batch's base case (L2). Besides its own answer --
        // kept in _base_case_v_solver_, the seed a set_init_from_n_powerflow() batch
        // starts every row from -- this is what builds the ledger, the Jacobian
        // sparsity and the factorization every row afterwards refactorizes into. That
        // is why it belongs to L2 and not to the results: keeping it is the whole
        // point of keeping a base case.
        //
        // Vinit_solver as Eigen::Ref relies on the reassignment below
        // (Vinit_solver = _algo.get_V()) always being same-size as the caller's
        // Vinit_solver: both trace back to the same ac_cache_.mat/dc_cache_.mat solver-space
        // dimension (nb_buses_solver_) for the duration of one call, so this holds
        // structurally, not by luck. No virtual dispatch here to enforce it --
        // if that invariant is ever broken, Eigen::Ref's operator= will assert
        // (debug) or corrupt memory (release), same risk as any Eigen::Ref sink.
        bool _solve_n_case(
            Eigen::Ref<CplxVect> Vinit_solver,  // is modified if _init_from_n_powerflow is true !
            size_t max_iter,
            real_type tol
        ){
            // The solver is NOT reset here: the caller did, before the preparation
            // hooks that configure it for the batch (BaseBatchSweep::compute) -- a
            // reset drops the PV pinning a generator-contingency sweep hands the
            // algorithm, which this solve must run with.
            _algo_controler.tell_all_changed();
            _algo.tell_solver_control(_algo_controler);
            _grid_model.get_generators().set_vm(Vinit_solver, active_layout().id_me_to_solver);
            CplxVect Vinit_solver2 = Vinit_solver;
            bool conv;
            // the "n" powerflow warm-up solve always needs the full complex V: its
            // result may seed every row (_init_from_n_powerflow below) and, for
            // ContingencyAnalysis / ScenarioSweep, feeds _record_n_case_violations --
            // so it must never be lazy, even when the per-row sweep that follows will be.
            _algo.set_lazy_v(false);
            if(_algo.ac_solver_used()){
                // ac_cache_.inj is already per-unit (pre_process_solver / fillSbus_me divides by
                // sn_mva when != 1), same convention as LSGrid::ac_pf's acSbus_ -- so tol
                // (a physical MW/MVAr tolerance) must be converted the same way LSGrid::ac_pf
                // does (`tol / sn_mva_`), or this initial solve accepts a per-unit mismatch
                // up to sn_mva times looser than what the caller asked for.
                conv = _algo.compute_pf(
                    ac_cache_.mat,
                    Vinit_solver2,
                    ac_cache_.inj,
                    active_layout().slack_bus_id_solver.as_eigen(),
                    active_layout().slack_weights,
                    active_layout().bus_pv.as_eigen(),
                    active_layout().bus_pq.as_eigen(),
                    max_iter,
                    tol / _grid_model.get_sn_mva());
            } else {
                conv = _algo.compute_pf_dc(
                    dc_cache_.mat,
                    Vinit_solver2,
                    dc_cache_.inj,
                    active_layout().slack_bus_id_solver.as_eigen(),
                    active_layout().slack_weights,
                    active_layout().bus_pv.as_eigen(),
                    active_layout().bus_pq.as_eigen());
            }
            if(conv) _base_case_v_solver_ = _algo.get_V();
            return conv;
        }

        // DC theta-only fast path: the shared magnitude every row's |V| is
        // reconstructed from (see _dc_vm_row_grid). Read off THIS call's starting
        // voltage, so it is redone on every compute() -- a kept base case says
        // nothing about where the caller chose to start.
        //
        // Magnitude is a pure echo of the input in DC (see BaseDCAlgo::compute_pf_dc)
        // and never changes across the sweep except where a row explicitly re-seeds it
        // (BaseBatchSweep::_apply_step_gen_v). Grid buses outside
        // active_layout().id_solver_to_me (eg an unused substation's second bus)
        // default to 0, not 1: they are never part of the solved system, and the
        // legacy (eager) _voltages was CplxMat::Zero(...)-initialized and never wrote
        // them -- 0 magnitude reproduces that "untouched column reads back as exact
        // complex 0" contract regardless of theta (also 0 there, for the same reason).
        void _init_dc_base_vm(const Eigen::Ref<const CplxVect> & Vinit_solver, size_t nb_total_bus){
            _dc_base_vm_solver_ = Vinit_solver.array().abs();
            _dc_base_vm_grid_ = RealVect::Zero(static_cast<Eigen::Index>(nb_total_bus));
            _dc_base_vm_grid_(active_layout().id_solver_to_me.as_eigen()) = _dc_base_vm_solver_;
        }

        // The two above, in the order compute() needs them, plus the L2 short circuit:
        // with the batch inputs still valid there is no "n" solve to redo -- the
        // sparsity, the ledger and the factorization on the algorithm are the ones
        // this batch needs, and the caller has already put this call's starting
        // voltage on the kept labelling (BaseBatchSweep::_vinit_on_grid_cache).
        bool _finish_preprocessing(
            size_t nb_steps,
            size_t nb_total_bus,
            Eigen::Ref<CplxVect> Vinit_solver,  // is modified if _init_from_n_powerflow is true !
            size_t max_iter,
            real_type tol,
            CustTimer  & timer_preproc,  // non const because double duration() is not const
            bool use_dc_lazy_v = false   // see BaseBatchSweep::compute()
        ){
            _size_result_buffers(nb_steps, nb_total_bus, use_dc_lazy_v);

            bool conv;
            if(_batch_inputs_valid_){
                _algo_controler.tell_none_changed();
                conv = true;   // it converged when it was built, or it was not kept
            } else {
                conv = _solve_n_case(Vinit_solver, max_iter, tol);
            }
            // check if we init the n-1 cases with results from the n cases or not
            if(conv && _init_from_n_powerflow) Vinit_solver = _base_case_v_solver_;
            if(use_dc_lazy_v && conv) _init_dc_base_vm(Vinit_solver, nb_total_bus);

            // everything init from n-case above
            _algo_controler.tell_none_changed();

            // end of pre processing
            _timer_pre_proc = timer_preproc.duration();
            return conv;
        }

        // number of computed rows, whichever accumulator (_voltages or, DC fast
        // path, _thetas) actually holds them.
        Eigen::Index _nb_result_rows() const {
            return _dc_lazy_storage_used_ ? _thetas.rows() : _voltages.rows();
        }

        // DC fast path only: row i's magnitude (grid-space, one entry per bus), used
        // to rebuild get_voltages()'s complex matrix and by the amps of a row. |V|
        // never changes across a DC solve (see BaseDCAlgo::compute_pf_dc) -- it is
        // either the shared base (_dc_gen_v_ never set: the common case, returned by
        // reference, nothing computed) or, when a row re-seeds generator voltage
        // targets (BaseBatchSweep::modify_gen_v), the base rescaled at those
        // generators' regulated buses -- delegated to GeneratorContainer::set_vm, the
        // exact function BaseBatchSweep::_apply_step_gen_v already uses live, so this
        // reconstruction is guaranteed consistent with what that row's solve actually
        // saw; written into the caller's `scratch`.
        const RealVect & _dc_vm_row_grid(Eigen::Index i, RealVect & scratch) const {
            if(static_cast<size_t>(i) < _dc_row_solved_.size() && !_dc_row_solved_[static_cast<size_t>(i)]){
                // never actually solved (eg an islanding contingency that diverges,
                // see _dc_row_solved_) -- 0 magnitude reproduces the legacy "untouched
                // row reads back as exact complex 0" contract regardless of theta.
                scratch = RealVect::Zero(_dc_base_vm_grid_.size());
                return scratch;
            }
            if(_dc_gen_v_.rows() == 0) return _dc_base_vm_grid_;
            CplxVect tmp = _dc_base_vm_solver_.cast<cplx_type>();
            const RealVect row = _dc_gen_v_.row(i);
            _grid_model.get_generators().set_vm(tmp, active_layout().id_me_to_solver, row);
            scratch = _dc_base_vm_grid_;
            scratch(active_layout().id_solver_to_me.as_eigen()) = tmp.array().abs();
            return scratch;
        }

        // DC fast path only: rebuilds the full complex _voltages from _thetas + the
        // reconstructed per-row magnitude -- see get_voltages(). Same total number of
        // std::polar calls as the pre-optimization eager path, but now genuinely
        // optional: paid only if/when a caller actually asks for complex voltages.
        void _rebuild_voltages_from_thetas() const {
            const Eigen::Index nb_steps = _thetas.rows();
            const Eigen::Index nb_bus = _thetas.cols();
            _voltages = CplxMat::Zero(nb_steps, nb_bus);
            RealVect scratch;
            for(Eigen::Index i = 0; i < nb_steps; ++i){
                const RealVect & vm_row = _dc_vm_row_grid(i, scratch);
                for(Eigen::Index j = 0; j < nb_bus; ++j){
                    _voltages(i, j) = std::polar(vm_row(j), _thetas(i, j));
                }
            }
        }

    protected:
        // ----- the three cache levels (see clear_grid_results / clear_batch_inputs /
        // clear_batch_outputs above) -------------------------------------------------
        // Each says "the corresponding tier of state is built and still describes this
        // batch". Raised by compute() as it builds each tier, lowered by the matching
        // clear_*(). There is no flag for L3: results are sized by every compute() and
        // never reused, so there is nothing to remember about them.
        bool _grid_cache_valid_ = false;
        bool _batch_inputs_valid_ = false;
        // the row count the kept batch inputs were built for: a batch of a different
        // size has different per-row state to prepare, so compute() drops L2 on a
        // change. _nb_steps_none (never a legal row count) means "nothing prepared".
        static const size_t _nb_steps_none = static_cast<size_t>(-1);
        size_t _prepared_nb_steps_ = _nb_steps_none;
        // the "n" solve's voltages, in solver space -- the seed a
        // set_init_from_n_powerflow() batch starts every row from. Part of L2: it is
        // that solve's answer, and it is kept with the factorization that produced it.
        CplxVect _base_case_v_solver_;

        bool _init_from_n_powerflow = false;
        // number of OS threads used to solve the batch (see set_nb_thread)
        int _nb_thread = 1;
        //timers
        double _timer_total = 0.;
        double _timer_pre_proc = 0.;
        double _timer_thread_init = 0.;
        // one entry per worker thread of the last multi-threaded compute(), harvested
        // after the joins (see BaseBatchSweep::_compute_threaded). Empty after a
        // single-threaded run: there the member _algo did everything.
        std::vector<LinearSolverStats> _thread_solver_stats_;

        // inputs
        LSGrid _grid_model;

        // properties of the grid
        const size_t n_line_;
        const size_t n_trafos_;
        const size_t n_total_;

        // algo
        AlgorithmSelector _algo;

        // outputs
        // mutable: get_voltages() rebuilds this lazily, from _thetas, on first request
        // when the DC fast path was used (see _rebuild_voltages_from_thetas) -- a pure
        // caching side effect of an otherwise-const accessor.
        mutable CplxMat _voltages;
        RealMat _amps_flows;
        RealMat _active_power_flows;

        // ----- DC theta-only fast path (see BaseAlgo::set_lazy_v) -----------------
        // true for the duration of one compute() call using it: DC, and not the
        // "handle disconnected grid" masked path (which stays on the legacy, always-
        // eager _voltages accumulation -- see BaseBatchSweep::_run_range_masked).
        bool _dc_lazy_storage_used_ = false;
        // grid-space theta accumulator (nb_steps x nb_total_bus), filled instead of
        // _voltages when _dc_lazy_storage_used_.
        RealMat _thetas;
        // one entry per row: whether _thetas.row(i) actually holds a converged solve
        // (0 initially, set by _run_range only where it writes that row). A row that
        // never converges (eg a contingency that islands the grid with no slack
        // reference in the stranded piece, in DC) is skipped there -- its theta stays
        // 0 (indistinguishable from a genuinely converged, exactly-flat-angle row) --
        // so the magnitude reconstruction below must consult this explicitly to
        // reproduce the legacy contract of an unwritten row (CplxMat::Zero(...),
        // never touched) reading back as exact complex 0, not (0-angle magnitude).
        std::vector<char> _dc_row_solved_;
        // solver-space / grid-space magnitude shared by every row that does not
        // re-seed a generator voltage target (see _dc_vm_row_grid).
        RealVect _dc_base_vm_solver_;
        RealVect _dc_base_vm_grid_;
        // copy of SbusPolicy::Vary::gen_v (empty if never set / not applicable, eg
        // ContingencyAnalysis) -- kept here, generic, so the magnitude reconstruction
        // helpers above do not need to know about SbusPolicy at all.
        RealMat _dc_gen_v_;


        // timers
        int _nb_solved = 0;
        int _nb_converged = 0;
        double _timer_compute_A = 0.;
        double _timer_compute_P = 0.;
        double _timer_solver = 0.;

        // solver control
        AlgoControl _algo_controler;

        // ---- what this batch solves, one object per family --------------------
        // A batch runs AC or DC, never both, so only one of these is ever built --
        // `active_layout()` picks it. They belong to the BATCH, not to
        // `_grid_model`: LSGrid::pre_process_solver fills whichever cache it is
        // handed, all of it, so nothing this batch builds can end up half in the
        // grid's cache and half in ours. That used to be exactly what happened --
        // the pv-pq split and the slack weights were not parameters, so they were
        // written into the grid's members while everything else came here, and had
        // to be fetched back out afterwards (three copies, and a grid left holding
        // a mixture). See SolverSideCache.hpp.
        AcSolverCache ac_cache_;
        DcSolverCache dc_cache_;
        /// which of the two above this sweep is building into; see active_layout()
        bool _ac_solver_used = true;
        int nb_buses_solver_;

};


} // namespace ls2g

#endif // BASEMULTIPLEPOWERFLOW_H