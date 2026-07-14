// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SECURITYANALYSIS_H
#define SECURITYANALYSIS_H

#include "BaseBatchSolverSynch.hpp"
#include "LimitViolation.hpp"
#include <set>
#include <exception>

namespace ls2g {

/**
Class to perform a contingency analysis (security analysis), which consist of performing some powerflow after some powerlines
have been disconnected 
 **/
class LS2G_API ContingencyAnalysis final: public BaseBatchSolverSynch
{
    public:
        explicit  ContingencyAnalysis(const LSGrid & init_grid_model,
                                       bool compute_limit_violations = false):
                            BaseBatchSolverSynch(init_grid_model),
                            _li_defaults(),
                            _li_coeffs(),
                            _compute_limit_violations_(compute_limit_violations),
                            _timer_modif_Ybus(0.),
                            _timer_thread_init(0.)
                            { }

        // whether limit violations are computed (see `converged`, `get_violations`,
        // `converged_n`, `get_violations_n` below). Defaults to `False` (settable at
        // construction too): computing violations means an extra per-element current/voltage
        // check inline in every contingency's solve, so users who only need
        // `compute_flows()`/`get_flows()` should not pay for it. Changing this flag clears any
        // previously-computed results (limit-violation results as well as voltages/flows,
        // same as `clear()`) so stale results computed under the other setting can never be
        // returned.
        bool get_compute_limit_violations() const noexcept {return _compute_limit_violations_;}
        void set_compute_limit_violations(bool val){
            if(val == _compute_limit_violations_) return;
            _compute_limit_violations_ = val;
            clear();
        }

        ~ContingencyAnalysis() noexcept = default;
        ContingencyAnalysis(const ContingencyAnalysis&) = delete;
        ContingencyAnalysis(ContingencyAnalysis&&) = delete;
        ContingencyAnalysis & operator=(ContingencyAnalysis&&) = delete;
        ContingencyAnalysis & operator=(const ContingencyAnalysis&) = delete;

        // utilities to add defaults to simulate
        void add_all_n1(){
            for(int l_id = 0; l_id < static_cast<int>(n_total_); ++l_id){
                std::set<int> this_default = {l_id};
                _li_defaults.insert(this_default);
            }
        }
        void add_n1(int line_id){
            check_ok_el(line_id);
            std::set<int> this_default = {line_id};
            _li_defaults.insert(this_default);
        }
        void add_multiple_n1(const std::vector<int> & vect_n1s){
            for(const auto line_id : vect_n1s){
                check_ok_el(line_id);
                std::set<int> this_default = {line_id};
                _li_defaults.insert(this_default);
            }
        }
        void add_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk)
            {
                check_ok_el(line_id);
                this_default.insert(line_id);
            }
            _li_defaults.insert(this_default);
        }

        // utilities to remove defaults to simulate (TODO)
        virtual void clear(){
            BaseBatchSolverSynch::clear();
            _li_defaults.clear();
            _li_coeffs.clear();
            _li_masked.clear();
            _skip_mask.clear();
            _converged.clear();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _timer_total = 0.;
            _timer_modif_Ybus = 0.;
            _timer_thread_init = 0.;
            _timer_pre_proc = 0.;
        }
        void clear_results_only(){
            BaseBatchSolverSynch::clear();
            _li_masked.clear();
            _skip_mask.clear();
            _converged.clear();
            _violations.clear();
            _converged_n_ = false;
            _violations_n_.clear();
            _timer_total = 0.;
            _timer_modif_Ybus = 0.;
            _timer_thread_init = 0.;
            _timer_pre_proc = 0.;
        }
        
        bool remove_n1(int line_id){
            check_ok_el(line_id);
            std::set<int> this_default = {line_id};
            auto nb_removed = _li_defaults.erase(this_default);
            return nb_removed >= 1;
        }
        size_t remove_multiple_n1(const std::vector<int> & vect_n1s){
            size_t nb_removed = 0;
            for(const auto line_id : vect_n1s){
                check_ok_el(line_id);
                std::set<int> this_default = {line_id};
                nb_removed += _li_defaults.erase(this_default);
            }
            return nb_removed;
        }
        bool remove_nk(const std::vector<int> & vect_nk){
            std::set<int> this_default;
            for(const auto line_id : vect_nk)
            {
                check_ok_el(line_id);
                this_default.insert(line_id);
            }
            auto nb_removed = _li_defaults.erase(this_default);
            return nb_removed >= 1;
        }

        // "handle disconnected grid" mode: when ON, a contingency that splits the grid
        // is no longer skipped. Instead the largest connected component is solved while
        // the buses of the other component(s) are "masked" (their NR equations are
        // forced to identity, their voltage reported as 0). Default: OFF (legacy
        // behaviour: such contingencies are skipped). Only the NR family supports it.
        bool get_handle_disconnected_grid() const {return _handle_disconnected_grid;}
        void set_handle_disconnected_grid(bool val) {_handle_disconnected_grid = val;}

        // number of OS threads used to solve the contingencies (default: 1).
        // With nb_thread == 1 the behaviour is identical to the legacy sequential
        // path. With nb_thread > 1 the contingency list is split into contiguous
        // ranges, each solved by its own thread (own solver + own Ybus copy),
        // writing to disjoint rows of the result matrix. Values < 1 are clamped to 1.
        int get_nb_thread() const {return _nb_thread;}
        void set_nb_thread(int n) {_nb_thread = (n < 1 ? 1 : n);}

        // make the computation. Throws if the pre-contingency ("n", no disconnection) powerflow
        // itself does not converge -- every contingency is solved starting from / relative to
        // that base case, so a diverging base case makes the whole analysis meaningless.
        void compute(const Eigen::Ref<const CplxVect> & Vinit, int max_iter, real_type tol);
        IntVect is_grid_connected_after_contingency();

        // Choose, over the currently-registered contingencies, the reference slack
        // bus (gridmodel numbering) that is stranded by the FEWEST of them — so the
        // "handle disconnected grid" mode can keep its (untouched) Jacobian reused
        // for as many contingencies as possible. Returns -1 if there is no slack.
        // Side-effect free w.r.t. the registered defaults; prepares the solver
        // inputs on demand (like is_grid_connected_after_contingency). Feed the
        // result to LSGrid::set_reference_slack_bus() before the base ac_pf so the
        // GPU companion (gpusim2grid) inherits this reference.
        int pick_reference_slack();

        Eigen::Ref<RealMat > compute_flows() {
            compute_flows_from_Vs();
            clean_flows();
            return _amps_flows;
        }

        Eigen::Ref<RealMat > compute_power_flows() {
            compute_flows_from_Vs(false);
            clean_flows(false);
            return _active_power_flows;
        }

        // python cannot handle set of sets... (c++ can because set are ordered set, but in python
        // they are unordered one, with hash functions)
        const std::set<std::set<int> > & my_defaults() const {return _li_defaults;}
        std::vector<std::vector<int> > my_defaults_vect() const {
            std::vector<std::vector<int> > res;
            res.reserve(_li_defaults.size());
            for(const auto & set_int: _li_defaults){
                std::vector<int> this_def(set_int.begin(), set_int.end());
                res.push_back(this_def);
            }
            return res;
        }

        // limit violations (only available if constructed with `compute_limit_violations=True`,
        // see get_compute_limit_violations() above; computed inline, per contingency, during
        // `compute()`). Row order matches `my_defaults_vect()` / `_voltages` rows.
        const std::vector<char> & converged() const {
            _check_limit_violations_enabled("converged");
            return _converged;
        }
        const std::vector<std::vector<LimitViolation> > & get_violations() const {
            _check_limit_violations_enabled("get_violations");
            return _violations;
        }
        // pre-contingency ("n") case
        bool converged_n() const {
            _check_limit_violations_enabled("converged_n");
            return _converged_n_;
        }
        const std::vector<LimitViolation> & get_violations_n() const {
            _check_limit_violations_enabled("get_violations_n");
            return _violations_n_;
        }

        // timers
        double total_time() const {return _timer_total;}
        double preprocessing_time() const {return _timer_pre_proc;}
        double modif_Ybus_time() const {return _timer_modif_Ybus;}
        double thread_init_time() const {return _timer_thread_init;}
        double solve_time() const {return _timer_solver;}

    protected:
        // prevent the insertion of "out of range" elements
        void check_ok_el(Eigen::Index el){
            if(el < 0){
                std::ostringstream exc_;
                exc_ << "SecurityAnalysis: cannot add the contingency with id ";
                exc_ << el << " contingency id should be > 0";
                throw std::runtime_error(exc_.str());
            }
            if(el >= static_cast<Eigen::Index>(n_total_)){
                std::ostringstream exc_;
                exc_ << "SecurityAnalysis: cannot add the contingency with id ";
                exc_ << el << " because the grid counts only " << n_total_ << " powerlines / trafos.";
                throw std::runtime_error(exc_.str());
            }
        }
        // raises if the user never asked for limit violations to be computed (see
        // get_compute_limit_violations())
        void _check_limit_violations_enabled(const std::string & fun_name) const {
            if(!_compute_limit_violations_){
                std::ostringstream exc_;
                exc_ << "ContingencyAnalysis::" << fun_name << ": limit violations were not "
                        "requested. Construct this object with `compute_limit_violations=True` "
                        "to use this feature.";
                throw std::runtime_error(exc_.str());
            }
        }

        void init_li_coeffs(bool ac_solver_used, const SolverBusIdVect &id_me_to_solver);
        // remove the line parameters from Ybus, this is to emulate its disconnection.
        // `algo` is the solver whose internal Ybus is updated in DC mode (one per thread).
        bool remove_from_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus, const std::vector<Coeff> & coeffs, bool ac_solver_used, AlgorithmSelector & algo);
        // after the coefficient has been removed with "remove_from_Ybus", add it back to Ybus
        void readd_to_Ybus(Eigen::SparseMatrix<cplx_type> & Ybus, const std::vector<Coeff> & coeffs, bool ac_solver_used, AlgorithmSelector & algo);

        // by default the flows are not 0 when the powerline is connected in the original topology
        // this function sorts this out
        void clean_flows(bool is_amps=true);

        // sometimes, when i perform some disconnection, I make the graph non connexe
        // in this case, well, i don't use the results of the simulation
        bool check_invertible(const Eigen::SparseMatrix<cplx_type> & Ybus) const;

        // ----- "handle disconnected grid" mode helpers (mask mode only) -----------
        // connected-component labelling of the admittance matrix: returns the solver
        // bus ids that are NOT part of the largest connected component (empty if it is
        // connected). Templated on the scalar type so it works on both the AC (complex
        // Ybus_) and the DC (real Bbus_) matrices. Explicitly instantiated in the .cpp.
        template<typename T>
        std::vector<int> disconnected_buses(const Eigen::SparseMatrix<T> & mat) const;
        // pre-pass run before the (n-)powerflow: fills _li_masked (the masked bus set
        // of each contingency), chooses the reference slack that minimises the number
        // of skipped contingencies (reordering slack_ids_me_ so it is index 0) and
        // fills _skip_mask (contingencies that still strand the chosen reference).
        void select_ref_slack_and_masks();
        // a copy of slack_weights_ with the masked slack buses zeroed and the vector
        // rescaled so its total is preserved (slack power stays among live slacks).
        RealVect masked_slack_weights(const std::vector<int> & masked) const;

        // solve the contingencies in [cont_begin, cont_end) using the passed solver
        // (one per thread), its own AlgoControl and a thread-local Ybus copy. Writes
        // directly to the (disjoint) rows of _voltages; book-keeping goes to the
        // passed accumulators (merged by the caller after all threads join). Any
        // thrown exception is captured into `err` instead of crossing the thread
        // boundary. This is the body shared by the single- and multi-threaded paths.
        void run_contingency_range(size_t cont_begin,
                                   size_t cont_end,
                                   AlgorithmSelector & algo,
                                   AlgoControl & control,
                                   Eigen::SparseMatrix<cplx_type> & Ybus,
                                   const CplxVect & Vinit_solver,
                                   bool ac_solver_used,
                                   bool mask_mode,
                                   int max_iter,
                                   real_type tol,
                                   real_type sn_mva,
                                   double & timer_modif_Ybus,
                                   int & nb_solved,
                                   double & timer_solver,
                                   std::exception_ptr & err,
                                   bool needs_solver_init,
                                   // limit-violation checking: only non-null when
                                   // _compute_limit_violations_ is true. Built once (before
                                   // the threads are spawned, see compute()) since it walks
                                   // the whole `_li_defaults` set -- avoids every contingency
                                   // (or every thread) rebuilding it via my_defaults_vect().
                                   // The per-bus / per-branch limit vectors themselves are NOT
                                   // threaded through here: they are trivial reference-return
                                   // accessors on the (already shared, read-only during this
                                   // phase) _grid_model, so run_contingency_range fetches them
                                   // directly where needed instead of adding 7 more parameters.
                                   const std::vector<std::vector<int> > * li_defaults_vect);

    private:
        // li_default
        std::set<std::set<int> > _li_defaults;  // do not use unordered_set here, we rely on the order for different functions !
        std::vector<std::vector<Coeff> > _li_coeffs;  // for each n-k, stores the coefficients I need to modify in the Ybus

        // number of OS threads used to solve the contingencies (see set_nb_thread)
        int _nb_thread = 1;

        // "handle disconnected grid" mode (see set_handle_disconnected_grid)
        bool _handle_disconnected_grid = false;
        std::vector<std::vector<int> > _li_masked;  // per contingency: masked solver bus ids (mask mode)
        std::vector<char> _skip_mask;               // per contingency: 1 => skip (strands the ref slack)

        // limit violations (see get_compute_limit_violations() / converged() / get_violations())
        bool _compute_limit_violations_ = false;
        std::vector<char> _converged;                          // per contingency, aligned with my_defaults_vect()
        std::vector<std::vector<LimitViolation> > _violations;  // per contingency, aligned with my_defaults_vect()
        bool _converged_n_ = false;                              // pre-contingency ("n") case
        std::vector<LimitViolation> _violations_n_;

        //timers
        double _timer_modif_Ybus;  // time to update the Ybus between the defaults simulation
        double _timer_thread_init;
};

} // namespace ls2g

#endif  //COMPUTERS_H
