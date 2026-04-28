// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef NR_ALGO_H
#define NR_ALGO_H

#include "BaseAlgo.hpp"
#include "NRSystem.hpp"
#include "ScalingPolicies.hpp"
#include "RefactorPolicies.hpp"

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
class NRAlgo : public BaseAlgo
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
        timer_refactor_(0.),
        timer_initialize_(0.),
        timer_dSbus_(0.),
        timer_fillJ_(0.),
        timer_Va_Vm_(0.),
        timer_pre_proc_(0.) {}

    virtual ~NRAlgo() noexcept = default;

    // ----- Jacobian accessor ---------------------------------------------------

    virtual
    Eigen::Ref<const Eigen::SparseMatrix<real_type>> get_J() const override {
        return _system.J();
    }

    virtual
    Eigen::SparseMatrix<real_type> get_J_python() const {
        Eigen::SparseMatrix<real_type> res = get_J();
        return res;
    }

    // ----- timers --------------------------------------------------------------

    virtual
    TimerJacType get_timers_jacobian() const override
    {
        return TimerJacType(timer_Fx_,
                            timer_solve_,
                            timer_refactor_,
                            timer_initialize_,
                            timer_check_,
                            timer_dSbus_,
                            timer_fillJ_,
                            timer_Va_Vm_,
                            timer_pre_proc_,
                            timer_total_nr_);
    }

    // ----- powerflow -----------------------------------------------------------

    virtual
    bool compute_pf(const Eigen::SparseMatrix<cplx_type>& Ybus,
                    CplxVect& V,
                    const CplxVect& Sbus,
                    Eigen::Ref<const IntVect> slack_ids,
                    const RealVect& slack_weights,
                    Eigen::Ref<const IntVect> pv,
                    Eigen::Ref<const IntVect> pq,
                    int max_iter,
                    real_type tol
                    ) override;

    virtual void reset() override;

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

    // MaxVoltageChange params
    real_type get_max_dVa() const { return max_dVa_; }
    void      set_max_dVa(real_type v) { max_dVa_ = v; }
    real_type get_max_dVm() const { return max_dVm_; }
    void      set_max_dVm(real_type v) { max_dVm_ = v; }

    // LineSearch (Armijo) params
    real_type get_ls_c()       const { return ls_c_; }
    void      set_ls_c(real_type v)  { ls_c_ = v; }
    real_type get_ls_rho()     const { return ls_rho_; }
    void      set_ls_rho(real_type v){ ls_rho_ = v; }
    int       get_ls_max_iter()const { return ls_max_iter_; }
    void      set_ls_max_iter(int v) { ls_max_iter_ = v; }

    // Iwamoto params
    real_type get_iw_mu_min() const { return iw_mu_min_; }
    void      set_iw_mu_min(real_type v) { iw_mu_min_ = v; }
    real_type get_iw_mu_max() const { return iw_mu_max_; }
    void      set_iw_mu_max(real_type v) { iw_mu_max_ = v; }

    // ----- refactor policy -----------------------------------------------------

    RefactorPolicyType get_refactor_policy()  const { return refactor_policy_; }
    void set_refactor_policy(RefactorPolicyType t)  { refactor_policy_ = t; }

    int  get_refactor_every_n() const { return refactor_every_n_; }
    void set_refactor_every_n(int v)  { refactor_every_n_ = v; }

    // ----- AlgoConfig serialization -------------------------------------------

    virtual AlgoConfig get_config() const override {
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

    virtual void set_config(const AlgoConfig& cfg) override {
        if (cfg.int_params.size()  < 4) throw std::runtime_error("NRAlgo::set_config: int_params must have at least 4 elements");
        if (cfg.real_params.size() < 6) throw std::runtime_error("NRAlgo::set_config: real_params must have at least 6 elements");
        max_dVa_         = static_cast<real_type>(cfg.real_params[0]);
        max_dVm_         = static_cast<real_type>(cfg.real_params[1]);
        ls_c_            = static_cast<real_type>(cfg.real_params[2]);
        ls_rho_          = static_cast<real_type>(cfg.real_params[3]);
        iw_mu_min_       = static_cast<real_type>(cfg.real_params[4]);
        iw_mu_max_       = static_cast<real_type>(cfg.real_params[5]);
        ls_max_iter_     = cfg.int_params[2];
        refactor_every_n_= cfg.int_params[3];
        refactor_policy_ = static_cast<RefactorPolicyType>(cfg.int_params[1]);
        set_scaling_policy(static_cast<ScalingPolicyType>(cfg.int_params[0]));
    }

    // ----- debug ---------------------------------------------------------------

    Eigen::SparseMatrix<real_type>
    create_jacobian_matrix_test(const Eigen::SparseMatrix<cplx_type>& Ybus,
                                const CplxVect& V,
                                const RealVect& slack_weights,
                                const Eigen::VectorXi& pq,
                                const Eigen::VectorXi& pvpq)
    {
        // DO NOT USE, FOR DEBUG ONLY
        (void)pvpq; (void)slack_weights;
        CplxVect Sbus_dummy(V.size()); Sbus_dummy.setZero();
        Eigen::VectorXi pv_dummy(0);
        IntVect slack_ids_dummy(1); slack_ids_dummy(0) = 0;
        RealVect sw_dummy = RealVect::Ones(V.size());
        sw_dummy /= sw_dummy.sum();
        _system.init_topology(Ybus, Sbus_dummy, slack_ids_dummy, sw_dummy, pv_dummy, pq);
        _system.update_state(Ybus, V, Sbus_dummy);
        _system.build_J_sparsity();
        _system.fill_J();
        return _system.J();
    }

protected:
    virtual void reset_timer() override {
        BaseAlgo::reset_timer();
        timer_refactor_   = 0.;
        timer_dSbus_      = 0.;
        timer_fillJ_      = 0.;
        timer_Va_Vm_      = 0.;
        timer_pre_proc_   = 0.;
        timer_initialize_ = 0.;
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
    double timer_refactor_;
    double timer_initialize_;
    double timer_dSbus_;
    double timer_fillJ_;
    double timer_Va_Vm_;
    double timer_pre_proc_;

    // No copy
    NRAlgo(const NRAlgo&) = delete;
    NRAlgo(NRAlgo&&) = delete;
    NRAlgo& operator=(NRAlgo&&) = delete;
    NRAlgo& operator=(const NRAlgo&) = delete;
};

#include "NRAlgo.tpp"

} // namespace ls2g

#endif // NR_ALGO_H
