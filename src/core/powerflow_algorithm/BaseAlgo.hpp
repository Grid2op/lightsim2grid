// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASEALGO_H
#define BASEALGO_H

#include "ls2g_api.hpp"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI
#include <memory>

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.hpp"
#include "CustTimer.hpp"
#include "BaseConstants.hpp"
#include "AlgoConfig.hpp"
#include "linear_solvers/LinearSolverStats.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

namespace ls2g {

class LSGrid;

struct TimerJac {
    double timer_Fx_         = -1.;
    double timer_solve_      = -1.;
    double timer_factor_     = -1.;
    double timer_refactor_   = -1.;
    double timer_initialize_ = -1.;
    double timer_check_      = -1.;
    double timer_dSbus_      = -1.;
    double timer_fillJ_      = -1.;
    double timer_Va_Vm_      = -1.;
    double timer_pre_proc_   = -1.;
    double timer_scale_      = -1.;
    double timer_mismatch_   = -1.;
    double timer_total_nr_   = -1.;
};
using TimerPTDFLODFType = std::tuple<double, double, double>;

/**
This class represents a algorithm to compute powerflow.

It can be derived for different usecase, for example for DC powerflow, AC powerflow using Newton Raphson method etc.
**/
class LS2G_API BaseAlgo : public BaseConstants
{
    public:
        const bool IS_AC;  // should be static ideally...

    public:
        explicit BaseAlgo(bool is_ac=true) noexcept:
            BaseConstants(),
            IS_AC(is_ac),
            n_(-1),
            err_(ErrorType::NotInitError),
            timer_Fx_(0.),
            timer_solve_(0.),
            timer_check_(0.),
            timer_total_nr_(0.),
            lsgrid_ptr_(nullptr){};

        virtual ~BaseAlgo() noexcept = default;

        // no copy allowed
        BaseAlgo(const BaseAlgo&) = delete;
        BaseAlgo(BaseAlgo&&) = delete;
        BaseAlgo & operator=(BaseAlgo&&) = delete;
        BaseAlgo & operator=(const BaseAlgo&) = delete;

        virtual void set_lsgrid(const LSGrid * gridmodel){
            lsgrid_ptr_ = gridmodel;
        }

        // ---- input validation, used ONLY by the *_with_input_validation
        // wrappers below (never by compute_pf / compute_pf_dc themselves) ----
        // compute_pf / compute_pf_dc are called, unchecked, on every hot
        // C++-internal path (LSGrid::ac_pf/dc_pf, and BaseBatchSolverSynch --
        // ie every contingency in ContingencyAnalysis, every timestep of
        // TimeSeries / SecurityAnalysis): those callers build Ybus / Sbus /
        // pv / pq / slack_ids themselves and can be assumed correct, so
        // paying an O(n) validation pass on each of those calls would be
        // pure overhead. The check is instead done once, at the actual
        // trust boundary: python calling a solver directly (see
        // compute_pf_with_input_validation below, bound as both
        // Solver.compute_pf and Solver.solve in binding_solvers.cpp).
        //
        // Throws std::runtime_error for structural problems: non-square
        // Ybus, V / Sbus (or Pbus) / slack_weights not of size n (the
        // slack weights are indexed *by bus id*, see eg SlackPolicies.tpp),
        // no slack bus at all, or a bus listed twice / in more than one of
        // slack_ids / pv / pq.
        // Throws std::out_of_range (python IndexError) for any id of
        // slack_ids / pv / pq outside [0, n).
        // `caller` names the solver in the error message (eg "NRAlgo::compute_pf").
        static void check_pf_inputs(
            const std::string & caller,
            Eigen::Index ybus_rows,
            Eigen::Index ybus_cols,
            Eigen::Index v_size,
            Eigen::Index sbus_size,
            const Eigen::Ref<const IntVect> & slack_ids,
            const Eigen::Ref<const RealVect> & slack_weights,
            const Eigen::Ref<const IntVect> & pv,
            const Eigen::Ref<const IntVect> & pq);

        // Throws std::runtime_error unless max_iter >= 0 (0 is a legitimate,
        // well-defined call: every solver builds its initial state / first
        // Jacobian before the iteration loop, so max_iter=0 deliberately
        // returns that pre-iteration state without taking any NR/GS step --
        // used eg by test_jacobian.py to check J iteration-by-iteration; only
        // negative values are nonsensical) and tol is a finite
        // strictly positive number (NaN and infinity are rejected). AC only:
        // DC is a single linear solve without iterations / tolerance.
        static void check_iter_tol(
            const std::string & caller,
            int max_iter,
            real_type tol);

        // Checked wrapper around the (virtual) compute_pf: validates with
        // check_pf_inputs + check_iter_tol, then dispatches to compute_pf on
        // whichever concrete solver `this` is. Not virtual: one shared
        // implementation (BaseAlgo.cpp) is enough since it only needs to
        // validate and then call the already-virtual compute_pf. This is
        // the method actually bound to python's Solver.compute_pf / .solve
        // (see binding_solvers.cpp) -- the raw compute_pf keeps its name and
        // stays reachable only from C++.
        // V is read-only here (see compute_pf's declaration above): this changes the
        // python-facing Solver.compute_pf/.solve parameter from a mutable numpy array
        // to a read-only one -- nothing in-tree relied on in-place mutation.
        // Ybus is a bare (non-Ref) Eigen::SparseMatrix, unlike every other Ybus/Bbus
        // parameter in this file: this is the one signature pybind11 actually binds
        // (see binding_solvers.cpp), and pybind11's Eigen support has no type_caster
        // for Eigen::Ref<const SparseMatrix<...>> (only for a concrete, owning
        // SparseMatrix -- see TODO in CHANGELOG.rst about writing one to avoid this
        // copy). Internally this is still cheap: it just forwards into the
        // Ref-taking compute_pf below, a single C++-side (not python-side) bind.
        bool compute_pf_with_input_validation(
            const Eigen::SparseMatrix<cplx_type> & Ybus,
            const Eigen::Ref<const CplxVect> & V,
            const Eigen::Ref<const CplxVect> & Sbus,
            const Eigen::Ref<const IntVect>  & slack_ids,
            const Eigen::Ref<const RealVect> & slack_weights,
            const Eigen::Ref<const IntVect>  & pv,
            const Eigen::Ref<const IntVect>  & pq,
            int                              max_iter,
            real_type                        tol);

        // Same idea for the DC entry point (validates, then calls the
        // virtual compute_pf_dc). Not bound to python today -- compute_pf_dc
        // itself isn't exposed under any name -- kept symmetric with the AC
        // wrapper above and ready if that binding gap is ever closed.
        // V stays plain: see the comment on the virtual compute_pf_dc above (resized).
        // Bbus is a bare (non-Ref) Eigen::SparseMatrix for the same pybind11 reason as
        // Ybus in compute_pf_with_input_validation above.
        bool compute_pf_dc_with_input_validation(
            const Eigen::SparseMatrix<real_type> & Bbus,
            const Eigen::Ref<const CplxVect> & V,
            const Eigen::Ref<const RealVect> & Pbus,
            const Eigen::Ref<const IntVect> & slack_ids,
            const Eigen::Ref<const RealVect> & slack_weights,
            const Eigen::Ref<const IntVect> & pv,
            const Eigen::Ref<const IntVect> & pq);

        virtual EigenRefConstRealSpMat get_J() const {
            throw std::runtime_error("AlgorithmSelector::get_J: There is not Jacobian matrix for this solver type.");
        }

        // bus_id -> Jacobian column of that bus' theta / vm / q unknown (-1 if none).
        // Only Newton-Raphson solvers define a Jacobian, hence these throw by default.
        virtual IntVect get_theta_to_J_col_python() const {
            throw std::runtime_error("get_theta_to_J_col: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_to_J_col_python() const {
            throw std::runtime_error("get_vm_to_J_col: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_to_J_col_python() const {
            throw std::runtime_error("get_q_to_J_col: only available for Newton-Raphson solvers.");
        }

        // bus_id -> Jacobian row of that bus' P / Q mismatch equation (-1 if none).
        // The row counterpart of the *_to_J_col maps. Only Newton-Raphson solvers
        // define a Jacobian, hence these throw by default.
        virtual IntVect get_p_to_J_row_python() const {
            throw std::runtime_error("get_p_to_J_row: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_to_J_row_python() const {
            throw std::runtime_error("get_q_to_J_row: only available for Newton-Raphson solvers.");
        }

        // Compact (bus, row/col) registration pair lists -- the row/col
        // counterpart of the *_to_J_col / *_to_J_row bus-keyed maps, but
        // preserving every registration (a bus may appear more than once, or
        // be absent from the bus-keyed map's CURRENT value if a later
        // registration shadowed it there). Only Newton-Raphson solvers define
        // a Jacobian, hence these throw by default.
        virtual IntVect get_p_buses_python() const {
            throw std::runtime_error("get_p_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_p_rows_python() const {
            throw std::runtime_error("get_p_rows: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_buses_python() const {
            throw std::runtime_error("get_q_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_q_rows_python() const {
            throw std::runtime_error("get_q_rows: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_theta_buses_python() const {
            throw std::runtime_error("get_theta_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_theta_cols_python() const {
            throw std::runtime_error("get_theta_cols: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_buses_python() const {
            throw std::runtime_error("get_vm_buses: only available for Newton-Raphson solvers.");
        }
        virtual IntVect get_vm_cols_python() const {
            throw std::runtime_error("get_vm_cols: only available for Newton-Raphson solvers.");
        }

        // VoltageControl (remote gen + SVC) converged reactive injection per
        // controller (pu) and its (kind, element id) identity, in controller
        // registration order. Empty for algorithms without the extension (DC,
        // Gauss-Seidel, NR without any voltage-mode controller).
        virtual RealVect get_controller_q()       const { return RealVect(); }
        virtual IntVect  get_controller_kind()    const { return IntVect(); }
        virtual IntVect  get_controller_elem_id() const { return IntVect(); }
        // J column of each controller's own Q unknown, in controller registration
        // order -- NOT the bus-keyed q_to_J_col (see NRSystem::controller_q_col's
        // own doc): needed whenever two controllers share a bus.
        virtual IntVect  get_controller_q_col()   const { return IntVect(); }

        // MultiSlack: J column of the slack_absorbed unknown (-1 when the
        // distributed-slack-in-Jacobian extension is not active).
        virtual int      get_slack_col()          const { return -1; }
        // MultiSlack: converged VALUE of the slack_absorbed unknown (pu; 0
        // when the extension is not active). NOT the same as the per-solve
        // initial guess of 0 -- this is the ground truth after convergence.
        virtual real_type get_slack_absorbed()     const { return static_cast<real_type>(0.); }

        // ---- capability flags -------------------------------------------------
        // Each concrete solver family declares its own `static constexpr bool`
        // (queryable at compile time even without an instance, e.g.
        // `NRAlgo<SparseLULinearSolver, SingleSlackNRSystem>::SUPPORTS_HVDC_DROOP`)
        // and overrides the matching virtual accessor below with `return
        // ITS_OWN_CONSTANT;` so that code holding only a `unique_ptr<BaseAlgo>`
        // (AlgorithmSelector, including plugin/Custom solvers) can query the
        // capability polymorphically. Defaults here are the conservative
        // "capability absent" value -- a plugin that does not override stays
        // maximally restricted, exactly like today's built-in Gauss-Seidel.
        static constexpr bool IS_DC = false;
        static constexpr bool SUPPORTS_HVDC_DROOP = false;
        static constexpr bool IS_FDPF = false;
        static constexpr bool SUPPORTS_REMOTE_VOLTAGE_CONTROL = false;
        static constexpr bool FILLS_BUS_MISMATCH = false;

        virtual bool is_dc() const noexcept { return IS_DC; }
        // Only the Newton-Raphson algorithms implement the hvdc angle-droop
        // ("AC emulation") equations (their NRSystem includes the Hvdc
        // extension); Gauss-Seidel / Fast-Decoupled do not, and a plugin
        // (Custom) solver defaults to `false` unless it overrides this --
        // see LSGrid::ac_pf's angle-droop guard.
        virtual bool supports_hvdc_droop() const noexcept { return SUPPORTS_HVDC_DROOP; }
        virtual bool is_fdpf() const noexcept { return IS_FDPF; }
        // Only the Newton-Raphson algorithms implement remote voltage control
        // (generators/SVCs regulating a bus other than their own -- their
        // NRSystem includes the VoltageControl extension, see
        // LSGrid::fill_voltage_control_solver_data). Gauss-Seidel /
        // Fast-Decoupled / DC do not consume that data at all.
        virtual bool supports_remote_voltage_control() const noexcept { return SUPPORTS_REMOTE_VOLTAGE_CONTROL; }
        /**
         * Does this algorithm leave a usable per-bus mismatch behind
         * (`get_bus_mismatch()`, size nb_bus in SOLVER numbering, at the voltage it
         * converged to)?
         *
         * LSGrid::compute_results reads it to work out what the generators and the
         * converter stations actually produced, so an algorithm that answers true
         * and then leaves the buffer empty (or the wrong size) would take the whole
         * result publication down with it. Built-in algorithms are covered by this
         * suite; an external (plugin) solver claiming the capability is CHECKED, once
         * per solve, in LSGrid::process_results -- see _check_solver_output.
         *
         * Defaults to false, so a plugin written against an older header keeps the
         * behaviour it has today: LSGrid falls back to deriving the mismatch itself.
         */
        virtual bool fills_bus_mismatch() const noexcept { return FILLS_BUS_MISMATCH; }

        Eigen::Ref<const RealVect> get_Va() const{
            return Va_;
        }
        Eigen::Ref<const RealVect> get_Vm() const{
            return Vm_;
        }
        Eigen::Ref<const CplxVect> get_V() const{
            return V_;
        }
        ErrorType get_error() const {
            return err_;
        }
        // Force the reported error state. Used by LSGrid to mark a solve as
        // non-converged when a (possibly plugin) solver returned non-finite
        // voltages while still claiming convergence -- see
        // LSGrid::process_results / LSGrid::_check_solver_output.
        /**
         * Restore the iteration count a solve reached. Same purpose as set_error:
         * LSGrid::process_results resets a diverged algorithm and puts the diagnosis
         * back, so that "it gave up after N iterations" survives the reset.
         */
        void set_nb_iter(int nb_iter) { nr_iter_ = nb_iter; }

        void set_error(ErrorType error) {
            err_ = error;
        }
        int get_nb_iter() const {
            return nr_iter_;
        }

        bool converged() const{
            return err_ == ErrorType::NoError;
        }

        std::tuple<double, double, double, double> get_timers() const
        {
            // TODO change the order of the timers here!
            auto res = std::tuple<double, double, double, double>(
              timer_Fx_, timer_solve_, timer_check_, timer_total_nr_);
            return res;
        }
        
        virtual TimerJac get_timers_jacobian() const
        {
            TimerJac res;
            res.timer_Fx_       = timer_Fx_;
            res.timer_solve_    = timer_solve_;
            res.timer_check_    = timer_check_;
            res.timer_total_nr_ = timer_total_nr_;
            return res;
        }
        
        virtual TimerPTDFLODFType get_timers_ptdf_lodf() const
        {
            TimerPTDFLODFType res = {
                -1.,  // not available for non NR solver, so I put -1
                -1.,  // not available for non NR solver, so I put -1
                -1.,  // not available for non NR solver, so I put -1
            };
            return res;
        }

        // Per-call counters and timings for the underlying linear solver (see
        // LinearSolverPolicy / RefactorRetryLinearSolver). Default: all-zero, meaning
        // "not tracked" -- overridden by NRAlgo/BaseDCAlgo (a single LinearSolver).
        // BaseFDPFAlgo has two linear solvers (B'/B'') and exposes them separately
        // (get_linear_solver_stats_bp/_bpp on the concrete FDPF_* Python type instead),
        // so it does not override this generic single-solver accessor.
        virtual LinearSolverStats get_linear_solver_stats() const { return LinearSolverStats{}; }

        // Complex AC entry point: solves V . (Ybus . V)* = Sbus.
        // Every AC solver overrides this; the DC solver does not (it uses `compute_pf_dc`) and
        // therefore inherits this throwing default (symmetric with `compute_pf_dc` below).
        // Ybus stays a plain reference (not Eigen::Ref): NRSystem::update_state caches
        // its address (Ybus_ptr_) across phases (update_state -> build_J_sparsity),
        // which needs the caller's actual, long-lived matrix, not a call-site Ref
        // temporary. Since compute_pf is one shared virtual signature (BaseAlgo,
        // NRAlgo, GaussSeidelAlgo, BaseFDPFAlgo all share it), every Ybus/Ybus-forwarding
        // parameter in this family must stay plain for the same reason.
        // V is the Vinit on input; none of the AC overriders write back into it (the
        // converged result is retrieved via get_V()), so it's read-only in practice --
        // unlike compute_pf_dc's V (BaseDCAlgo resizes/reassigns it), which stays plain.
        virtual
        bool compute_pf(const EigenRefConstCplxSpMat     & /*Ybus*/,
                        const Eigen::Ref<const CplxVect> & /*V*/,
                        const Eigen::Ref<const CplxVect> & /*Sbus*/,
                        const Eigen::Ref<const IntVect>  & /*slack_ids*/,
                        const Eigen::Ref<const RealVect> & /*slack_weights*/,
                        const Eigen::Ref<const IntVect>  & /*pv*/,
                        const Eigen::Ref<const IntVect>  & /*pq*/,
                        int /*max_iter*/,
                        real_type /*tol*/
                        ){
            throw std::runtime_error("compute_pf (complex AC entry point) is not available for this solver (DC solvers use compute_pf_dc).");
        }

        // Native real-valued DC entry point: DC only needs `Bbus . theta = Pbus` (all real).
        // Only the DC solver overrides this; every other solver type throws.
        // `V` carries the complex initial voltage (slack angle + voltage setpoints) on input
        // and the complex result on output. There is no max_iter / tol: DC is a single linear solve.
        // V stays a plain reference (not Eigen::Ref): BaseDCAlgo::compute_pf_dc
        // resizes/reassigns it (V.resize(...), V = CplxVect() on the early-error
        // path) -- unlike compute_pf's V (AC family), which is read-only and was
        // converted to Eigen::Ref. compute_pf/compute_pf_dc are independent virtual
        // chains (different overrider sets), so there's no need for them to match.
        virtual
        bool compute_pf_dc(const EigenRefConstRealSpMat     & /*Bbus*/,
                           const Eigen::Ref<const CplxVect> & /*V*/,
                           const Eigen::Ref<const RealVect> & /*Pbus*/,
                           const Eigen::Ref<const IntVect>  & /*slack_ids*/,
                           const Eigen::Ref<const RealVect> & /*slack_weights*/,
                           const Eigen::Ref<const IntVect>  & /*pv*/,
                           const Eigen::Ref<const IntVect>  & /*pq*/){
            throw std::runtime_error("compute_pf_dc is only available for DC solvers.");
        }

        void tell_solver_control(const AlgoControl & solver_control){
            _solver_control = solver_control;
        }

        // Batch DC sweeps (BaseBatchSweep's theta-only fast path) call this before a
        // run of per-step compute_pf_dc() calls: when true, the DC solver only builds
        // Va_ (theta) and skips Vm_/V_ entirely. DC magnitude never changes across a
        // solve (it is a pure echo of the caller's V, see BaseDCAlgo::compute_pf_dc),
        // so the batch caller reconstructs it itself, from its own (much smaller)
        // input data, only if/when something actually asks for it (get_voltages(),
        // current-limit checks) -- see BaseBatchSolverSynch. No-op default: AC
        // solvers always need V_ built by the solve itself, and a solver that does
        // not override this (eg a plugin) keeps building it eagerly, exactly as
        // before this existed.
        virtual void set_lazy_v(bool) {}
        virtual bool lazy_v() const { return false; }
        virtual void reset();
        // TODO speed: prevent copy and use Eigen::Ref here
        virtual RealMat get_ptdf(){
            throw std::runtime_error("Impossible to get the PTDF matrix with this solver type.");
        }
        // TODO speed: prevent copy and use Eigen::Ref here
        virtual RealMat get_lodf(const Eigen::Ref<const IntVect> & /*from_bus*/,
                                 const Eigen::Ref<const IntVect> & /*to_bus*/){  // TODO interface is likely to change
            throw std::runtime_error("Impossible to get the LODF matrix with this solver type.");
        }
        virtual Eigen::SparseMatrix<real_type> get_bsdf(){  // TODO interface is likely to change
            throw std::runtime_error("Impossible to get the BSDF matrix with this solver type.");
        }

        virtual void update_internal_Ybus(const Coeff & /*new_coeffs*/, bool /*add*/){
            throw std::runtime_error("Function update_internal_Ybus not implemented in general.");
        }

        // bus masking (ContingencyAnalysis "handle disconnected grid" mode): forces
        // the given solver buses' equations to identity so an isolated island does
        // not make the system singular, without changing the matrix sparsity. Only
        // the Newton-Raphson family supports it (see supports_bus_masking); the
        // default is a no-op so other algorithms are unaffected.
        virtual bool supports_bus_masking() const { return false; }
        virtual void set_masked_buses(const std::vector<int> & /*solver_bus_ids*/) {}

        virtual AlgoConfig get_config() const { return AlgoConfig{}; }
        virtual void set_config(const AlgoConfig&) {}
        
    protected:
        virtual void reset_timer(){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_check_ = 0.;
            timer_total_nr_ = 0.;
        }

        bool is_linear_solver_valid() const {
            // bool res = true;
            // if((err_ == ErrorType::NotInitError) || (err_ == ErrorType::LicenseError)) res = false;  // cannot use a non intialize solver
            // return res;
            return (err_ != ErrorType::LicenseError);
        }
        RealVect _evaluate_Fx(const EigenRefConstCplxSpMat     & Ybus,
                              const Eigen::Ref<const CplxVect> & V,
                              const Eigen::Ref<const CplxVect> & Sbus,
                              size_t                           slack_id,  // id of the slack bus
                              real_type                        slack_absorbed,
                              const Eigen::Ref<const RealVect> & slack_weights,
                              const Eigen::Ref<const IntVect>  & pv,
                              const Eigen::Ref<const IntVect>  & pq);

        RealVect _evaluate_Fx(const EigenRefConstCplxSpMat     & Ybus,
                              const Eigen::Ref<const CplxVect> & V,
                              const Eigen::Ref<const CplxVect> & Sbus,
                              const Eigen::Ref<const IntVect>  & pv,
                              const Eigen::Ref<const IntVect>  & pq);

        bool _check_for_convergence(const Eigen::Ref<const RealVect> & F,
                                    real_type tol);

        bool _check_for_convergence(const Eigen::Ref<const RealVect> & p,
                                    const Eigen::Ref<const RealVect> & q,
                                    real_type tol);

        Eigen::VectorXi extract_slack_bus_id(const Eigen::Ref<const IntVect> & pv,
                                             const Eigen::Ref<const IntVect> & pq,
                                             unsigned int nb_bus);

        /**
        When there are multiple slacks, add the other "slack buses" in the PV buses indexes
        (behaves as if only the first element is used for the slack !!!, called "ref slack")

        `slack_ids(0)` MUST be the intended reference / angle slack whenever this is called --
        every algorithm family (DC via BaseDCAlgo, Fast-Decoupled via BaseFDPFAlgo, Gauss-Seidel)
        funnels through this one method for that choice, and the (single-slack) Newton-Raphson
        path makes the same assumption independently in NRSystem's MultiSlack extension
        (`ref_slack_id_ = slack_ids[0]`, see NRSystem.hpp). Reordering `slack_ids` so a different
        bus ends up at index 0 (eg ContingencyAnalysis::select_ref_slack_and_masks) silently
        changes which bus is removed from the reduced linear system / used as the angle datum --
        not a crash, just wrong (or non-reproducible) results that are miserable to debug back to
        this assumption.
        **/
        Eigen::VectorXi retrieve_pv_with_slack(const Eigen::Ref<const IntVect> & slack_ids,
                                               const Eigen::Ref<const IntVect> & pv) const {
            if(slack_ids.size() > 1){
                const auto nb_slack_added = slack_ids.size() - 1;
                Eigen::VectorXi my_pv = Eigen::VectorXi(pv.size() + nb_slack_added);
                for(auto i = 0; i < nb_slack_added; ++i){
                    my_pv(i) = slack_ids(i+1);
                }
                for(auto i = 0; i < pv.size(); ++i){
                    my_pv(i + nb_slack_added) = pv(i);
                }
                return my_pv;
            }else{
                return pv;
            }
        }

        /**
        When there are multiple slacks, add the other "slack buses" in the PV buses indexes
        **/
        Eigen::VectorXi add_slack_to_pv(const Eigen::Ref<const IntVect> & slack_ids,
                                        const Eigen::Ref<const IntVect> & pv) const {
            Eigen::VectorXi my_pv = Eigen::VectorXi(slack_ids.size() + pv.size());
            my_pv << slack_ids, pv;
            return my_pv;
        }
        
        // terribly inefficient way to know if an element is in a vector
        bool isin(int k, const Eigen::Ref<const Eigen::VectorXi> & vect) const{
            for(auto el : vect){
                if(el == k) return true;
            }
            return false;
        }

        /**
         * The per-bus complex power mismatch left by the last solve, in solver bus
         * numbering. Empty on an algorithm that has never run (or has been reset),
         * and on one whose `fills_bus_mismatch()` is false. See mis_bus_ for
         * exactly what it contains.
         *
         * PUBLIC on purpose: LSGrid::compute_results reads it back through
         * AlgorithmSelector to publish what the generators and converter stations
         * produced, which is the whole point of the capability flag above. It sits
         * in a `protected:` block otherwise, hence the explicit `public:` here and
         * the one that closes it.
         */
    public:
        Eigen::Ref<const CplxVect> get_bus_mismatch() const {return mis_bus_;}
    protected:

        // Bf / Bf_T are resized and filled from scratch by fillBf_for_PTDF: a real
        // reference is needed, Eigen::Ref<SparseMatrix> can't resize/reserve.
        void get_Bf(Eigen::SparseMatrix<real_type> & Bf) const;
        void get_Bf_transpose(Eigen::SparseMatrix<real_type> & Bf_T) const;
        
    protected:
        // solver initialization
        int n_;

        /**
         * The per-bus complex power mismatch of the last evaluation:
         * V .* conj(Ybus * V) - Sbus, plus whatever the algorithm's components inject
         * (the distributed slack's share, the hvdc droop flows, a voltage controller's
         * reactive output). Size nb_bus, in SOLVER bus numbering.
         *
         * Held here, and not in each algorithm, because both families compute exactly
         * this and both need it to persist across iterations -- and because it is the
         * one intermediate a caller may legitimately want after the solve (see
         * get_bus_mismatch). Filled by the concrete algorithm: BaseFDPFAlgo writes it
         * directly, NRAlgo hands its address to the NRSystem it owns.
         *
         * It is the RAW mismatch in both families. The FDPF used to divide it by Vm in
         * place -- the `mis / Vm` its P and Q rows need -- which left the same member
         * meaning two different things depending on who last wrote it; that division
         * now happens where the rows are extracted.
         */
        CplxVect mis_bus_;
        /// Ybus * V of the last evaluation, kept out of the expression above so Eigen
        /// does not have to allocate a temporary for the sparse * dense product.
        CplxVect ybus_v_;

        // solution of the problem
        RealVect Vm_;  // voltage magnitude
        RealVect Va_;  // voltage angle
        CplxVect V_;  // complex voltage

        int nr_iter_;  // number of iteration performs by the solver (may vary depending on the solver)
        ErrorType err_; //error message:
        // -1 : the solver has not been initialized (call initialize in this case)
        // 0 everything ok
        // 1: i can't factorize the matrix (klu_factor)
        // 2: i can't refactorize the matrix (klu_refactor)
        // 3: i can't solve the system (klu_solve)
        // 4: end of possible iterations (divergence because nr_iter_ >= max_iter)

        // timers
        double timer_Fx_;
        double timer_solve_;
        double timer_check_;
        double timer_total_nr_;

        const LSGrid * lsgrid_ptr_;  // does not have ownership so that's fine (pointer to the base gridmodel, can be used for some powerflow)
        AlgoControl _solver_control;

};


} // namespace ls2g

#endif // BASEALGO_H