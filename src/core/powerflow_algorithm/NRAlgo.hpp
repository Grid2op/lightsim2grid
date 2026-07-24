// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_ALGO_H
#define NR_ALGO_H

#include <cmath>     // std::isfinite (policy parameter validation)
#include <sstream>
#include <stdexcept>

#include "BaseAlgo.hpp"
#include "NRSystem.hpp"
#include "ScalingPolicies.hpp"
#include "RefactorPolicies.hpp"
#include "linear_solvers/LinearSolverStats.hpp"

namespace ls2g {
  
/**
 * Unified Newton-Raphson solver parameterised by LinearSolver and NRSystem.
 *
 * NRSystem is either MultiSlackNRSystem (distributed slack, lag=1) or
 * SingleSlackNRSystem (traditional single slack, lag=0).
 *
 * Step-scaling and refactorization strategies are runtime-configurable enums,
 * not template parameters, so the same binary can switch strategy at run time.
 */
template<class LinearSolver, class NRSystem>
class NRAlgo final : public BaseAlgo
{
public:
    NRAlgo() noexcept :
        BaseAlgo(true),
        need_factorize_(true),
        scaling_policy_(create_scaling_policy<NRSystem>(ScalingPolicyType::NoScaling)),
        refactor_policy_(RefactorPolicyType::AlwaysRefactor),
        max_dVa_(static_cast<real_type>(0.5)),
        max_dVm_(static_cast<real_type>(0.1)),
        ls_c_(static_cast<real_type>(1e-4)),
        ls_rho_(static_cast<real_type>(0.5)),
        ls_max_iter_(10),
        iw_mu_min_(static_cast<real_type>(1e-4)),
        iw_mu_max_(static_cast<real_type>(1.0)),
        refactor_every_n_(4),
        timer_dSbus_(0.),
        timer_fillJ_(0.),
        timer_Va_Vm_(0.),
        timer_pre_proc_(0.),
        timer_scale_(0.),
        timer_mismatch_(0.) {}

    ~NRAlgo() noexcept override = default;

    // Both SingleSlackNRSystem and MultiSlackNRSystem include the Hvdc
    // extension (see NRSystem.hpp), so every NRAlgo instantiation supports
    // the angle-droop ("AC emulation") equations regardless of NRSystem.
    static constexpr bool SUPPORTS_HVDC_DROOP = true;
    bool supports_hvdc_droop() const noexcept override { return SUPPORTS_HVDC_DROOP; }

    // Both SingleSlackNRSystem and MultiSlackNRSystem also include the
    // VoltageControl extension (remote gen / SVC regulation), so every
    // NRAlgo instantiation supports it regardless of NRSystem.
    static constexpr bool SUPPORTS_REMOTE_VOLTAGE_CONTROL = true;
    bool supports_remote_voltage_control() const noexcept override { return SUPPORTS_REMOTE_VOLTAGE_CONTROL; }

    // ----- Jacobian accessor ---------------------------------------------------

    Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const override {
        return _system.J();
    }

    Eigen::SparseMatrix<real_type> get_J_python() const {
        Eigen::SparseMatrix<real_type> res = get_J();
        return res;
    }

    // ----- column (unknown) -> bus-id converters (NR-only) ---------------------
    // bus_id -> Jacobian column of that bus' theta / vm / q unknown (-1 if none).
    // Only meaningful after a topology has been initialised (i.e. after at least
    // one compute_pf), same lifetime as get_J.
    IntVect get_theta_to_J_col_python() const override { return _to_intvect(_system.theta_to_J_col()); }
    IntVect get_vm_to_J_col_python()    const override { return _to_intvect(_system.vm_to_J_col()); }
    IntVect get_q_to_J_col_python()     const override { return _to_intvect(_system.q_to_J_col()); }

    // ----- row (equation) -> bus-id converters (NR-only) -----------------------
    // bus_id -> Jacobian row of that bus' P / Q mismatch equation (-1 if none).
    // Same lifetime / semantics as the *_to_J_col converters above.
    IntVect get_p_to_J_row_python() const override { return _to_intvect(_system.p_to_J_row()); }
    IntVect get_q_to_J_row_python() const override { return _to_intvect(_system.q_to_J_row()); }

    // ----- compact (bus, row/col) registration pair lists (NR-only) ------------
    // Row/col counterpart of the *_to_J_col / *_to_J_row bus-keyed maps above,
    // but preserving every registration (see NRSystem::p_buses() etc.).
    IntVect get_p_buses_python()     const override { return _to_intvect(_system.p_buses()); }
    IntVect get_p_rows_python()      const override { return _to_intvect(_system.p_rows()); }
    IntVect get_q_buses_python()     const override { return _to_intvect(_system.q_buses()); }
    IntVect get_q_rows_python()      const override { return _to_intvect(_system.q_rows()); }
    IntVect get_theta_buses_python() const override { return _to_intvect(_system.theta_buses()); }
    IntVect get_theta_cols_python()  const override { return _to_intvect(_system.theta_cols()); }
    IntVect get_vm_buses_python()    const override { return _to_intvect(_system.vm_buses()); }
    IntVect get_vm_cols_python()     const override { return _to_intvect(_system.vm_cols()); }

    // ----- VoltageControl (remote gen + SVC) converged results -----------------
    RealVect get_controller_q()       const override { return _system.controller_q(); }
    IntVect  get_controller_kind()    const override { return _system.controller_kind(); }
    IntVect  get_controller_elem_id() const override { return _system.controller_elem_id(); }
    IntVect  get_controller_q_col()   const override { return _system.controller_q_col(); }
    int      get_slack_col()          const override { return _system.slack_col(); }
    real_type get_slack_absorbed()    const override { return _system.slack_absorbed(); }

    // ----- timers --------------------------------------------------------------

    TimerJac get_timers_jacobian() const override
    {
        const LinearSolverStats lsstats = detail::get_stats_impl(_linear_solver, 0);
        TimerJac res;
        res.timer_Fx_         = timer_Fx_;
        res.timer_solve_      = lsstats.timer_solve_;
        res.timer_factor_     = lsstats.timer_factor_;
        res.timer_refactor_   = lsstats.timer_refactor_;
        res.timer_initialize_ = lsstats.timer_initialize_;
        res.timer_check_      = timer_check_;
        res.timer_dSbus_      = timer_dSbus_;
        res.timer_fillJ_      = timer_fillJ_;
        res.timer_Va_Vm_      = timer_Va_Vm_;
        res.timer_pre_proc_   = timer_pre_proc_;
        res.timer_scale_      = timer_scale_;
        res.timer_mismatch_   = timer_mismatch_;
        res.timer_total_nr_   = timer_total_nr_;
        return res;
    }

    // Per-call counters and timings for the underlying linear solver (factor/refactor/
    // analyze/solve counts, refactor-failure and fallback-factor counts if the wrapped
    // LinearSolver tracks them -- see LinearSolverPolicy/RefactorRetryLinearSolver).
    // Returns an all-zero LinearSolverStats if LinearSolver doesn't track any.
    LinearSolverStats get_linear_solver_stats() const override {
        return detail::get_stats_impl(_linear_solver, 0);
    }

    // ----- powerflow -----------------------------------------------------------

    bool compute_pf(
        const EigenRefConstCplxSpMat     & Ybus,
        const Eigen::Ref<const CplxVect> & V,
        const Eigen::Ref<const CplxVect> & Sbus,
        const Eigen::Ref<const IntVect>  & slack_ids,
        const Eigen::Ref<const RealVect> & slack_weights,
        const Eigen::Ref<const IntVect>  & pv,
        const Eigen::Ref<const IntVect>  & pq,
        int                              max_iter,
        real_type                        tol
    ) override;

    void reset() override;

    // ----- bus masking ---------------------------------------------------------
    bool supports_bus_masking() const override { return true; }
    void set_masked_buses(const std::vector<int> & solver_bus_ids) override {
        _system.set_masked_buses(solver_bus_ids);
    }
    
    // ----- scaling policy ------------------------------------------------------
    ScalingPolicyType get_scaling_policy_type()  const { return scaling_policy_->type(); }
    void set_scaling_policy(ScalingPolicyType t)  { 
        scaling_policy_ = create_scaling_policy<NRSystem>(t);
        update_scaling_policy_params<NRSystem>(
            scaling_policy_.get(),
            max_dVa_, max_dVm_,
            ls_c_, ls_rho_, ls_max_iter_,
            iw_mu_min_, iw_mu_max_
        );
    }

    // ----- policy parameter validation ----------------------------------------
    // These parameters reach the solver from two untrusted-ish directions: the
    // python setters below (exposed by `lightsim2grid.algorithm`, ie the
    // `newtonpf` entry point) and set_config(), which is fed by an AlgoConfig --
    // part of the serialized grid state, so it also arrives from a pickle or a
    // binary file. They must be validated on the way in: `refactor_every_n_` in
    // particular is used as `iter % refactor_every_n_` (see
    // should_refactor_policy), where 0 is an integer division by zero -- a
    // SIGFPE that kills the process instead of raising a catchable error.
    static int _checked_refactor_every_n(int v) {
        if (v < 1) {
            std::ostringstream exc_;
            exc_ << "NRAlgo: refactor_every_n must be >= 1, got " << v
                 << " (it is used as `iter % refactor_every_n`, so 0 would be a division by zero).";
            throw std::runtime_error(exc_.str());
        }
        return v;
    }
    static real_type _checked_in_0_1(real_type v, const char * name) {
        if (!std::isfinite(v) || (v <= 0.) || (v >= 1.)) {
            std::ostringstream exc_;
            exc_ << "NRAlgo: " << name << " must be a finite number in (0, 1), got " << v << ".";
            throw std::runtime_error(exc_.str());
        }
        return v;
    }
    static real_type _checked_positive(real_type v, const char * name) {
        if (!std::isfinite(v) || (v <= 0.)) {
            std::ostringstream exc_;
            exc_ << "NRAlgo: " << name << " must be a finite, strictly positive number, got " << v << ".";
            throw std::runtime_error(exc_.str());
        }
        return v;
    }
    static real_type _checked_finite(real_type v, const char * name) {
        if (!std::isfinite(v)) {
            std::ostringstream exc_;
            exc_ << "NRAlgo: " << name << " must be a finite number, got " << v << ".";
            throw std::runtime_error(exc_.str());
        }
        return v;
    }
    static int _checked_non_negative(int v, const char * name) {
        if (v < 0) {
            std::ostringstream exc_;
            exc_ << "NRAlgo: " << name << " must be >= 0, got " << v << ".";
            throw std::runtime_error(exc_.str());
        }
        return v;
    }

    // MaxVoltageChange params
    real_type get_max_dVa() const { return max_dVa_; }
    void      set_max_dVa(real_type v) { max_dVa_ = _checked_positive(v, "max_dVa"); }
    real_type get_max_dVm() const { return max_dVm_; }
    void      set_max_dVm(real_type v) { max_dVm_ = _checked_positive(v, "max_dVm"); }

    // LineSearch (Armijo) params
    real_type get_ls_c()       const { return ls_c_; }
    void      set_ls_c(real_type v)  { ls_c_ = _checked_in_0_1(v, "ls_c"); }
    real_type get_ls_rho()     const { return ls_rho_; }
    void      set_ls_rho(real_type v){ ls_rho_ = _checked_in_0_1(v, "ls_rho"); }
    int       get_ls_max_iter()const { return ls_max_iter_; }
    void      set_ls_max_iter(int v) { ls_max_iter_ = _checked_non_negative(v, "ls_max_iter"); }

    // Iwamoto params
    real_type get_iw_mu_min() const { return iw_mu_min_; }
    void      set_iw_mu_min(real_type v) { iw_mu_min_ = _checked_finite(v, "iw_mu_min"); }
    real_type get_iw_mu_max() const { return iw_mu_max_; }
    void      set_iw_mu_max(real_type v) { iw_mu_max_ = _checked_finite(v, "iw_mu_max"); }

    // ----- refactor policy -----------------------------------------------------

    RefactorPolicyType get_refactor_policy()  const { return refactor_policy_; }
    void set_refactor_policy(RefactorPolicyType t)  { refactor_policy_ = t; }

    int  get_refactor_every_n() const { return refactor_every_n_; }
    void set_refactor_every_n(int v)  { refactor_every_n_ = _checked_refactor_every_n(v); }

    // ----- AlgoConfig serialization -------------------------------------------

    AlgoConfig get_config() const override {
        AlgoConfig cfg;
        cfg.int_params  = { static_cast<int>(scaling_policy_->type()),
                            static_cast<int>(refactor_policy_),
                            ls_max_iter_,
                            refactor_every_n_ };
        cfg.real_params = { static_cast<double>(max_dVa_),
                            static_cast<double>(max_dVm_),
                            static_cast<double>(ls_c_),
                            static_cast<double>(ls_rho_),
                            static_cast<double>(iw_mu_min_),
                            static_cast<double>(iw_mu_max_) };
        return cfg;
    }

    void set_config(const AlgoConfig& cfg) override {
        if (cfg.int_params.size()  < 4) throw std::runtime_error("NRAlgo::set_config: int_params must have at least 4 elements");
        if (cfg.real_params.size() < 6) throw std::runtime_error("NRAlgo::set_config: real_params must have at least 6 elements");
        // NB an AlgoConfig is part of the serialized grid state, so these values
        // can come straight from a pickle or a binary file: validate them exactly
        // like the setters do (see the note on _checked_refactor_every_n).
        max_dVa_         = _checked_positive(static_cast<real_type>(cfg.real_params[0]), "max_dVa");
        max_dVm_         = _checked_positive(static_cast<real_type>(cfg.real_params[1]), "max_dVm");
        ls_c_            = _checked_in_0_1(static_cast<real_type>(cfg.real_params[2]), "ls_c");
        ls_rho_          = _checked_in_0_1(static_cast<real_type>(cfg.real_params[3]), "ls_rho");
        iw_mu_min_       = _checked_finite(static_cast<real_type>(cfg.real_params[4]), "iw_mu_min");
        iw_mu_max_       = _checked_finite(static_cast<real_type>(cfg.real_params[5]), "iw_mu_max");
        ls_max_iter_     = _checked_non_negative(cfg.int_params[2], "ls_max_iter");
        refactor_every_n_= _checked_refactor_every_n(cfg.int_params[3]);
        refactor_policy_ = static_cast<RefactorPolicyType>(cfg.int_params[1]);
        set_scaling_policy(static_cast<ScalingPolicyType>(cfg.int_params[0]));
    }

    // ----- debug ---------------------------------------------------------------

    // TODO make it work !
    // Eigen::SparseMatrix<real_type>
    // create_jacobian_matrix_test(const Eigen::SparseMatrix<cplx_type>& Ybus,
    //                             const CplxVect& V,
    //                             const RealVect& slack_weights,
    //                             const Eigen::VectorXi& pq,
    //                             const Eigen::VectorXi& pvpq)
    // {
    //     // DO NOT USE, FOR DEBUG ONLY
    //     (void)pvpq; (void)slack_weights;
    //     CplxVect Sbus_dummy(V.size()); Sbus_dummy.setZero();
    //     Eigen::VectorXi pv_dummy(0);
    //     IntVect slack_ids_dummy(1); slack_ids_dummy(0) = 0;
    //     RealVect sw_dummy = RealVect::Ones(V.size());
    //     sw_dummy /= sw_dummy.sum();
    //     _system.update_state(BaseAlgo::lsgrid_ptr_, Ybus, V, Sbus, slack_weights);
    //     _system.init_topology(slack_ids, slack_weights, pv_dummy, pq);
    //     // _system.init_topology(Ybus, Sbus_dummy, slack_ids_dummy, sw_dummy, pv_dummy, pq);
    //     // _system.update_state(Ybus, V, Sbus_dummy);
    //     _system.build_J_sparsity();
    //     _system.fill_J();
    //     return _system.J();
    // }

protected:
    void reset_timer() override {
        BaseAlgo::reset_timer();
        detail::reset_stats_timers_impl(_linear_solver, 0);
        timer_dSbus_      = 0.;
        timer_fillJ_      = 0.;
        timer_Va_Vm_      = 0.;
        timer_pre_proc_   = 0.;
        timer_scale_      = 0.;
        timer_mismatch_      = 0.;
        _system.reset_timers();
    }

    bool should_refactor_policy(int iter) const {
        switch (refactor_policy_) {
            case RefactorPolicyType::AlwaysRefactor: return true;
            case RefactorPolicyType::EveryN:         return (iter % refactor_every_n_) == 1;
            case RefactorPolicyType::Chord:          return iter == 1;
            default:                                 return true;
        }
    }

private:
    static IntVect _to_intvect(const std::vector<int>& v) {
        return Eigen::Map<const IntVect>(v.data(), static_cast<Eigen::Index>(v.size()));
    }

    LinearSolver _linear_solver;
    NRSystem     _system;

    bool need_factorize_;

    // Runtime policy enums
    std::unique_ptr<ScalingPolicy<NRSystem> >  scaling_policy_;
    RefactorPolicyType refactor_policy_;

    // MaxVoltageChange params
    real_type max_dVa_;
    real_type max_dVm_;

    // LineSearch (Armijo) params
    real_type ls_c_;
    real_type ls_rho_;
    int       ls_max_iter_;

    // Iwamoto params
    real_type iw_mu_min_;
    real_type iw_mu_max_;

    // EveryN param
    int refactor_every_n_;

    // Timers
    double timer_dSbus_;
    double timer_fillJ_;
    double timer_Va_Vm_;
    double timer_pre_proc_;
    double timer_scale_;
    double timer_mismatch_;

    // No copy
    NRAlgo(const NRAlgo&) = delete;
    NRAlgo(NRAlgo&&) = delete;
    NRAlgo& operator=(NRAlgo&&) = delete;
    NRAlgo& operator=(const NRAlgo&) = delete;
};

#include "NRAlgo.tpp"

} // namespace ls2g

#endif // NR_ALGO_H
