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
        void set_init_from_n_powerflow(bool do_it) noexcept {_init_from_n_powerflow = do_it;}

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
            _nb_thread = (n < 1 ? 1 : n);
        }

        // solver "control"
        virtual void change_algorithm(const AlgorithmType & type){
            _algo.change_algorithm(type);
            this->clear();
        }
        // String-based overload: looks up the solver by registry name, so it
        // works for plugin solvers (and built-ins with no dedicated AlgorithmType
        // member, eg NRRefactorRetry_*) which change_algorithm(AlgorithmType)
        // cannot reach (AlgorithmType::Custom is rejected there).
        virtual void change_algorithm(const std::string & name){
            _algo.change_algorithm(name);
            this->clear();
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
        void set_algo_config(const AlgoConfig & cfg) {_algo.set_config(cfg); }

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
        virtual void clear() {
            _algo.reset();
            _amps_flows = RealMat();
            _active_power_flows = RealMat();
            _voltages = CplxMat();
            _thetas = RealMat();
            _dc_row_solved_.clear();
            _dc_lazy_storage_used_ = false;
            _dc_base_vm_solver_ = RealVect();
            _dc_base_vm_grid_ = RealVect();
            _dc_gen_v_ = RealMat();
            _dc_vm_cache_ = RealMat();
            _nb_solved = 0;
            _nb_converged = 0;
            _timer_compute_A = 0.;
            _timer_compute_P = 0.;
            _timer_solver = 0.;
            _timer_thread_init = 0.;
            // NB: _nb_thread is deliberately NOT reset -- it is a setting, not a result.
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
        template<class T>
        void compute_amps_flows(const T & structure_data,
                                real_type sn_mva,
                                size_t lag_id,
                                bool is_trafo)
        {
            const auto & bus_vn_kv = _grid_model.get_bus_vn_kv();
            const auto & el_status = structure_data.get_status_global();
            const auto & status1 = structure_data.get_status_side_1();
            const auto & status2 = structure_data.get_status_side_2();
            const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
            const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
            bool is_ac = _algo.ac_solver_used();

            // AC uses complex (Kron-reduced) coefficients, DC uses real susceptance coefficients
            Eigen::Ref<const CplxVect> vect_yac_ff = structure_data.yac_eff_11();
            Eigen::Ref<const CplxVect> vect_yac_ft = structure_data.yac_eff_12();
            Eigen::Ref<const RealVect> vect_ydc_ff = structure_data.ydc_11();
            Eigen::Ref<const RealVect> vect_ydc_ft = structure_data.ydc_12();
            Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift(); // not used in AC nor if it's powerline anyway

            size_t nb_el = structure_data.nb();
            real_type sqrt_3 = sqrt(3.);
            // DC theta-only fast path (see BaseAlgo::set_lazy_v / BaseBatchSweep::compute):
            // when active, _voltages is empty and the per-bus values are read straight from
            // _thetas (real, no .arg() needed) + the small _dc_* reconstruction inputs
            // (no .abs() on a full complex column needed either) instead.
            const bool dc_lazy = !is_ac && _dc_lazy_storage_used_;
            const Eigen::Index nb_steps = _nb_result_rows();

            RealVect res;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                const bool s1 = status1[el_id];
                const bool s2 = status2[el_id];

                // retrieve which buses are used; a half-open branch (see
                // keep_half_open_lines) has bus_id == _deactivated_bus_id on
                // its open side and must not be used to index _voltages /
                // bus_vn_kv -- substitute an open-end voltage of exactly 0
                // instead (both sides treated the same way). For AC, yac_eff_*
                // is already Kron-reduced for whichever side is open, so this
                // alone gives the correct "or"-side (side 1) flow either way;
                // DC has no such reduction (handled explicitly below).
                GlobalBusId bus_from_me = bus_from(el_id);
                GlobalBusId bus_to_me = bus_to(el_id);

                const real_type bus_vn_kv_f = s1 ? bus_vn_kv(bus_from_me.cast_int())
                                                  : (s2 ? bus_vn_kv(bus_to_me.cast_int()) : real_type(1.));

                // TODO speed: this copies a full nb_steps-sized column per line/trafo,
                // every call. _voltages/_thetas are RowMajor, so a column isn't contiguous
                // and can't bind to Eigen::Ref; the ternary also needs a common type with
                // Zero(nb_steps). Avoiding the copy would need a control-flow restructure
                // (separate open-side / closed-side code paths) rather than a
                // reference-type change -- see CHANGELOG [TODO].
                CplxVect Efrom, Eto;  // AC, or DC without the fast path (_voltages-backed)
                RealVect theta_from, theta_to;  // DC fast path (_thetas-backed)
                RealVect v_f_kv;
                if(dc_lazy){
                    theta_from = s1 ? RealVect(_thetas.col(bus_from_me.cast_int())) : RealVect::Zero(nb_steps);
                    theta_to   = s2 ? RealVect(_thetas.col(bus_to_me.cast_int()))   : RealVect::Zero(nb_steps);
                    v_f_kv = (s1 ? _dc_vm_col(bus_from_me.cast_int(), nb_steps) : _dc_vm_col(bus_to_me.cast_int(), nb_steps)) * bus_vn_kv_f;
                } else {
                    Efrom = s1 ? CplxVect(_voltages.col(bus_from_me.cast_int())) : CplxVect::Zero(nb_steps);
                    Eto   = s2 ? CplxVect(_voltages.col(bus_to_me.cast_int()))   : CplxVect::Zero(nb_steps);
                    // vn_kv base for the amps conversion: use whichever side is
                    // actually energized. If side 1 (the one being measured) is
                    // open the numerator below is exactly 0 regardless (AC: Efrom
                    // == 0 forces S_ft == 0; DC: forced explicitly), so the choice
                    // of base here only has to avoid a spurious 0/0 division.
                    v_f_kv = (s1 ? Efrom.array().abs() : Eto.array().abs()) * bus_vn_kv_f;
                }

                if(is_ac){
                    // retrieve physical parameters (complex)
                    const cplx_type y_ff = vect_yac_ff(el_id);  // scalar
                    const cplx_type y_ft = vect_yac_ft(el_id);
                    // trafo equations (to get the power at the "from" side)
                    CplxVect I_ft =  y_ff * Efrom + y_ft * Eto;
                    I_ft = I_ft.array().conjugate();
                    const CplxVect S_ft = Efrom.array() * I_ft.array();

                    // now compute the current flow
                    res = S_ft.array().abs() * sn_mva;
                }else{
                    // unlike yac_eff_*, ydc_11/ydc_12 are NOT Kron-reduced for a
                    // half-open branch: DC treats one side open as fully
                    // disconnected (see fillBdc: "disco on one side == disco on
                    // both sides"), so report 0 rather than mixing ydc_ff/ydc_ft
                    // with a meaningless open-end angle.
                    if(!(s1 && s2)){
                        _amps_flows.col(el_id + lag_id).setZero();
                        continue;
                    }
                    // DC active flow from the bus angles (theta) directly, like the gridmodel
                    // results: P = ydc_ff . theta_from + ydc_ft . theta_to  (theta = bus voltage angle)
                    const real_type y_ff = vect_ydc_ff(el_id);
                    const real_type y_ft = vect_ydc_ft(el_id);
                    if(!dc_lazy){
                        theta_from = Efrom.array().arg();
                        theta_to = Eto.array().arg();
                    }
                    res = (y_ff * theta_from.array() + y_ft * theta_to.array()) * sn_mva;
                    if(is_trafo) res.array() -= dc_x_tau_shift(el_id);
                    res.array() = res.array().abs();
                }
                res.array() /= sqrt_3 * v_f_kv.array();
                _amps_flows.col(el_id + lag_id) = res;
            }
        }
        template<class T>
        void compute_active_power_flows(const T & structure_data,
                                        real_type sn_mva,
                                        size_t lag_id,
                                        bool is_trafo)
        {
            const auto & el_status = structure_data.get_status_global();
            const auto & status1 = structure_data.get_status_side_1();
            const auto & status2 = structure_data.get_status_side_2();
            const GlobalBusIdVect & bus_from = structure_data.get_bus_id_side_1();
            const GlobalBusIdVect & bus_to = structure_data.get_bus_id_side_2();
            const bool is_ac = _algo.ac_solver_used();

            // AC uses complex (Kron-reduced) coefficients, DC uses real susceptance coefficients
            Eigen::Ref<const CplxVect> vect_yac_ff = structure_data.yac_eff_11();
            Eigen::Ref<const CplxVect> vect_yac_ft = structure_data.yac_eff_12();
            Eigen::Ref<const RealVect> vect_ydc_ff = structure_data.ydc_11();
            Eigen::Ref<const RealVect> vect_ydc_ft = structure_data.ydc_12();
            Eigen::Ref<const RealVect> dc_x_tau_shift = structure_data.dc_x_tau_shift(); // not used in AC nor if it's powerline anyway

            size_t nb_el = structure_data.nb();
            // see compute_amps_flows above -- same DC theta-only fast path
            const bool dc_lazy = !is_ac && _dc_lazy_storage_used_;
            const Eigen::Index nb_steps = _nb_result_rows();

            RealVect res;
            for(size_t el_id = 0; el_id < nb_el; ++el_id){
                if(!el_status[el_id]) continue;

                const bool s1 = status1[el_id];
                const bool s2 = status2[el_id];

                // a half-open branch (see keep_half_open_lines) has bus_id ==
                // _deactivated_bus_id on its open side; substitute an open-end
                // voltage of exactly 0 rather than indexing _voltages with it
                // (both sides treated the same way). For AC, yac_eff_* is
                // already Kron-reduced for whichever side is open, so this
                // alone gives the correct "or"-side (side 1) flow either way
                // (0 when side 1 itself is open, matching
                // TwoSidesContainer_rxh_A::compute_results_tsc_rxha_no_amps);
                // DC has no such reduction (handled explicitly below).
                GlobalBusId bus_from_me = bus_from(el_id);
                GlobalBusId bus_to_me = bus_to(el_id);
                // TODO speed: same structural copy as in compute_amps_flows above (see
                // the TODO there for why this isn't a simple Eigen::Ref fix) -- see
                // CHANGELOG [TODO].
                CplxVect Efrom, Eto;
                RealVect theta_from, theta_to;
                if(dc_lazy){
                    // DC active power needs only theta -- no magnitude, no .arg() round trip
                    theta_from = s1 ? RealVect(_thetas.col(bus_from_me.cast_int())) : RealVect::Zero(nb_steps);
                    theta_to   = s2 ? RealVect(_thetas.col(bus_to_me.cast_int()))   : RealVect::Zero(nb_steps);
                } else {
                    Efrom = s1 ? CplxVect(_voltages.col(bus_from_me.cast_int())) : CplxVect::Zero(nb_steps);
                    Eto   = s2 ? CplxVect(_voltages.col(bus_to_me.cast_int()))   : CplxVect::Zero(nb_steps);
                }

                // trafo equations (to get the power at the "from" side)
                if(is_ac){
                    // retrieve physical parameters (complex)
                    const cplx_type y_ff = vect_yac_ff(el_id);  // scalar
                    const cplx_type y_ft = vect_yac_ft(el_id);
                    CplxVect I_ft = y_ff * Efrom + y_ft * Eto;
                    I_ft = I_ft.array().conjugate();
                    const CplxVect S_ft = Efrom.array() * I_ft.array();

                    // now compute the active flow
                    res = S_ft.array().real() * sn_mva;
                }else{
                    // unlike yac_eff_*, ydc_11/ydc_12 are NOT Kron-reduced for a
                    // half-open branch: DC treats one side open as fully
                    // disconnected (see fillBdc: "disco on one side == disco on
                    // both sides"), so report 0 rather than mixing ydc_ff/ydc_ft
                    // with a meaningless open-end angle.
                    if(!(s1 && s2)){
                        _active_power_flows.col(el_id + lag_id).setZero();
                        continue;
                    }
                    // DC active flow from the bus angles (theta) directly, like the gridmodel
                    // results: P = ydc_ff . theta_from + ydc_ft . theta_to  (theta = bus voltage angle)
                    const real_type y_ff = vect_ydc_ff(el_id);
                    const real_type y_ft = vect_ydc_ft(el_id);
                    if(!dc_lazy){
                        theta_from = Efrom.array().arg();
                        theta_to = Eto.array().arg();
                    }
                    res = (y_ff * theta_from.array() + y_ft * theta_to.array()) * sn_mva;
                    if(is_trafo) res.array() -= dc_x_tau_shift(el_id);
                }
                _active_power_flows.col(el_id + lag_id) = res;
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
        // read-only member Bbus_ is only read here (safe to share across threads).
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
        // (works on the member Ybus_ / Bbus_ / Sbus_ / Pbus_ — all read-only).
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

        CplxVect prepare_solver_input_base(const Eigen::Ref<const CplxVect> & Vinit, bool ac_solver_used){
            // clear previous data
            Sbus_ = CplxVect();
            Ybus_ = Eigen::SparseMatrix<cplx_type>();
            Pbus_ = RealVect();
            Bbus_ = Eigen::SparseMatrix<real_type>();
            id_solver_to_me_.clear();
            id_me_to_solver_.clear();
            slack_ids_solver_ = SolverBusIdVect();
            slack_ids_me_ = GlobalBusIdVect();
            nb_buses_solver_ = -1;

            // fill the data correctly
            _algo_controler.tell_all_changed();
            CplxVect res;
            if(ac_solver_used){
                res = _grid_model.pre_process_solver(
                    Vinit,
                    Sbus_,
                    Ybus_,
                    id_me_to_solver_,
                    id_solver_to_me_,
                    slack_ids_me_,
                    slack_ids_solver_,
                    ac_solver_used,
                    _algo_controler);
                nb_buses_solver_ = static_cast<int>(Ybus_.cols());
            } else {
                // native real DC path: build the real Bbus / Pbus
                res = _grid_model.pre_process_dc_solver(
                    Vinit,
                    Pbus_,
                    Bbus_,
                    id_me_to_solver_,
                    id_solver_to_me_,
                    slack_ids_me_,
                    slack_ids_solver_,
                    _algo_controler);
                nb_buses_solver_ = static_cast<int>(Bbus_.cols());
            }

            // extract relevant information
            // ask for the family this batch actually solved with: the grid keeps
            // one pv-pq split and one set of slack weights PER family, and the
            // family-less accessors answer for AC whenever an AC powerflow has run
            // on this grid -- which is not necessarily the one we just built.
            const SolverBusIdVect & gm_bus_pv = ac_solver_used ? _grid_model.get_ac_pv_solver() : _grid_model.get_dc_pv_solver();
            const SolverBusIdVect & gm_bus_pq = ac_solver_used ? _grid_model.get_ac_pq_solver() : _grid_model.get_dc_pq_solver();
            const RealVect gm_bus_sw = ac_solver_used ? _grid_model.get_ac_slack_weights_solver() : _grid_model.get_dc_slack_weights_solver();

            // TODO copies are made here, which is not ideal
            bus_pv_ = gm_bus_pv;  // was _to_intvect()
            bus_pq_ = gm_bus_pq;  // was _to_intvect()
            slack_weights_ = gm_bus_sw;
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

        // Vinit_solver as Eigen::Ref relies on the reassignment below
        // (Vinit_solver = _algo.get_V()) always being same-size as the caller's
        // Vinit_solver: both trace back to the same Ybus_/Bbus_ solver-space
        // dimension (nb_buses_solver_) for the duration of one call, so this holds
        // structurally, not by luck. No virtual dispatch here to enforce it --
        // if that invariant is ever broken, Eigen::Ref's operator= will assert
        // (debug) or corrupt memory (release), same risk as any Eigen::Ref sink.
        bool _finish_preprocessing(
            size_t nb_steps,
            size_t nb_total_bus,
            Eigen::Ref<CplxVect> Vinit_solver,  // is modified if _init_from_n_powerflow is true !
            size_t max_iter,
            real_type tol,
            CustTimer  & timer_preproc,  // non const because double duration() is not const
            bool use_dc_lazy_v = false  // see BaseBatchSweep::compute()
        ){

                // init the results matrices: the DC theta-only fast path accumulates into
                // _thetas (real) instead of _voltages (complex) -- see get_voltages() /
                // compute_amps_flows / compute_active_power_flows.
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
                _dc_vm_cache_ = RealMat();
                _amps_flows = RealMat::Zero(0, n_total_);
                _active_power_flows = RealMat::Zero(0, n_total_);

                // reset the solver
                _algo.reset();

                // perform the initial powerflow / "powerflow in n"
                // (needed to init the underlying solver with the correct sparsity pattern in particular)
                _algo_controler.tell_all_changed();
                _algo.tell_solver_control(_algo_controler);
                _grid_model.get_generators().set_vm(Vinit_solver, id_me_to_solver_);
                CplxVect Vinit_solver2 = Vinit_solver;
                bool conv;
                // the "n" powerflow warm-up solve always needs the full complex V: its
                // result may seed every row (_init_from_n_powerflow below) and, for
                // ContingencyAnalysis / ScenarioSweep, feeds _record_n_case_violations --
                // so it must never be lazy, even when the per-row sweep that follows will be.
                _algo.set_lazy_v(false);
                if(_algo.ac_solver_used()){
                    // Sbus_ is already per-unit (pre_process_solver / fillSbus_me divides by
                    // sn_mva when != 1), same convention as LSGrid::ac_pf's acSbus_ -- so tol
                    // (a physical MW/MVAr tolerance) must be converted the same way LSGrid::ac_pf
                    // does (`tol / sn_mva_`), or this initial solve accepts a per-unit mismatch
                    // up to sn_mva times looser than what the caller asked for.
                    conv = _algo.compute_pf(
                        Ybus_,
                        Vinit_solver2,
                        Sbus_,
                        slack_ids_solver_.as_eigen(),
                        slack_weights_,
                        bus_pv_.as_eigen(),
                        bus_pq_.as_eigen(),
                        max_iter,
                        tol / _grid_model.get_sn_mva());
                } else {
                    conv = _algo.compute_pf_dc(
                        Bbus_,
                        Vinit_solver2,
                        Pbus_,
                        slack_ids_solver_.as_eigen(),
                        slack_weights_,
                        bus_pv_.as_eigen(),
                        bus_pq_.as_eigen());
                }

                // check if we init the n-1 cases with results from the n cases
                // or not
                if(_init_from_n_powerflow) Vinit_solver = _algo.get_V();

                if(use_dc_lazy_v && conv){
                    // magnitude is a pure echo of the input in DC (see BaseDCAlgo::compute_pf_dc)
                    // and never changes across the sweep except where a row explicitly re-seeds it
                    // (BaseBatchSweep::_apply_step_gen_v) -- this is the shared base every row's
                    // magnitude is reconstructed from, see _dc_vm_row_grid.
                    // Grid buses outside id_solver_to_me_ (eg an unused substation's second bus)
                    // default to 0, not 1: they are never part of the solved system, and the legacy
                    // (eager) _voltages was CplxMat::Zero(...)-initialized and never wrote them --
                    // 0 magnitude reproduces that "untouched column reads back as exact complex 0"
                    // contract regardless of theta (also 0 there, for the same reason).
                    _dc_base_vm_solver_ = Vinit_solver.array().abs();
                    _dc_base_vm_grid_ = RealVect::Zero(static_cast<Eigen::Index>(nb_total_bus));
                    _dc_base_vm_grid_(id_solver_to_me_.as_eigen()) = _dc_base_vm_solver_;
                }

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
        // to rebuild get_voltages()'s complex matrix. |V| never changes across a DC
        // solve (see BaseDCAlgo::compute_pf_dc) -- it is either the shared base
        // (_dc_gen_v_ never set: the common case, nothing to redo per row) or, when a
        // row re-seeds generator voltage targets (BaseBatchSweep::modify_gen_v), the
        // base rescaled at those generators' regulated buses -- delegated to
        // GeneratorContainer::set_vm, the exact function BaseBatchSweep::
        // _apply_step_gen_v already uses live, so this reconstruction is guaranteed
        // consistent with what that row's solve actually saw.
        RealVect _dc_vm_row_grid(Eigen::Index i) const {
            if(static_cast<size_t>(i) < _dc_row_solved_.size() && !_dc_row_solved_[static_cast<size_t>(i)]){
                // never actually solved (eg an islanding contingency that diverges,
                // see _dc_row_solved_) -- 0 magnitude reproduces the legacy "untouched
                // row reads back as exact complex 0" contract regardless of theta.
                return RealVect::Zero(_dc_base_vm_grid_.size());
            }
            if(_dc_gen_v_.rows() == 0) return _dc_base_vm_grid_;
            CplxVect tmp = _dc_base_vm_solver_.cast<cplx_type>();
            const RealVect row = _dc_gen_v_.row(i);
            _grid_model.get_generators().set_vm(tmp, id_me_to_solver_, row);
            RealVect res = _dc_base_vm_grid_;
            res(id_solver_to_me_.as_eigen()) = tmp.array().abs();
            return res;
        }

        // DC fast path only: all `nb_steps` rows' magnitude at a single grid bus
        // (one column of what get_voltages() would reconstruct). The common case
        // (_dc_gen_v_ never set) is a plain broadcast, no reconstruction needed; the
        // row-varying case is served from a cache built once (not once per branch).
        // Either way, a never-solved row (_dc_row_solved_) is forced back to 0 -- see
        // _dc_vm_row_grid.
        RealVect _dc_vm_col(Eigen::Index bus_id, Eigen::Index nb_steps) const {
            if(_dc_gen_v_.rows() == 0){
                RealVect res = RealVect::Constant(nb_steps, _dc_base_vm_grid_(bus_id));
                for(Eigen::Index i = 0; i < nb_steps; ++i){
                    if(static_cast<size_t>(i) < _dc_row_solved_.size() && !_dc_row_solved_[static_cast<size_t>(i)]) res(i) = 0.;
                }
                return res;
            }
            _ensure_dc_vm_cache(nb_steps);
            return _dc_vm_cache_.col(bus_id);
        }

        void _ensure_dc_vm_cache(Eigen::Index nb_steps) const {
            if(_dc_vm_cache_.rows() == nb_steps) return;
            const Eigen::Index nb_bus = _dc_base_vm_grid_.size();
            _dc_vm_cache_ = RealMat::Zero(nb_steps, nb_bus);
            for(Eigen::Index i = 0; i < nb_steps; ++i){
                _dc_vm_cache_.row(i) = _dc_vm_row_grid(i).transpose();
            }
        }

        // DC fast path only: rebuilds the full complex _voltages from _thetas + the
        // reconstructed per-row magnitude -- see get_voltages(). Same total number of
        // std::polar calls as the pre-optimization eager path, but now genuinely
        // optional: paid only if/when a caller actually asks for complex voltages.
        void _rebuild_voltages_from_thetas() const {
            const Eigen::Index nb_steps = _thetas.rows();
            const Eigen::Index nb_bus = _thetas.cols();
            _voltages = CplxMat::Zero(nb_steps, nb_bus);
            for(Eigen::Index i = 0; i < nb_steps; ++i){
                const RealVect vm_row = _dc_vm_row_grid(i);
                for(Eigen::Index j = 0; j < nb_bus; ++j){
                    _voltages(i, j) = std::polar(vm_row(j), _thetas(i, j));
                }
            }
        }

    protected:
        bool _init_from_n_powerflow = false;
        // number of OS threads used to solve the batch (see set_nb_thread)
        int _nb_thread = 1;
        //timers
        double _timer_total = 0.;
        double _timer_pre_proc = 0.;
        double _timer_thread_init = 0.;

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
        // lazily-built cache of _dc_vm_row_grid for every row, used by compute_amps_flows
        // so a magnitude that does vary per row (_dc_gen_v_ set) is reconstructed once,
        // not once per branch endpoint that reads it.
        mutable RealMat _dc_vm_cache_;


        // timers
        int _nb_solved = 0;
        int _nb_converged = 0;
        double _timer_compute_A = 0.;
        double _timer_compute_P = 0.;
        double _timer_solver = 0.;

        // solver control
        AlgoControl _algo_controler;

        // internal data
        CplxVect Sbus_;
        Eigen::SparseMatrix<cplx_type> Ybus_;
        // DC counterparts (real): the DC solver works on the real system Bbus . theta = Pbus
        RealVect Pbus_;
        Eigen::SparseMatrix<real_type> Bbus_;
        GlobalBusIdVect id_solver_to_me_;
        SolverBusIdVect id_me_to_solver_;
        SolverBusIdVect slack_ids_solver_;
        GlobalBusIdVect slack_ids_me_;
        int nb_buses_solver_;

        // TODO everything is copied here
        // not optimal...
        SolverBusIdVect bus_pv_;
        SolverBusIdVect bus_pq_;
        RealVect slack_weights_;

};


} // namespace ls2g

#endif // BASEMULTIPLEPOWERFLOW_H