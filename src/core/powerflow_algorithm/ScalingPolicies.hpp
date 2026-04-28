// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SCALING_POLICIES_H
#define SCALING_POLICIES_H


#include "ls2g_api.hpp"


namespace ls2g {

/**
 * Runtime step-scaling strategies for the Newton-Raphson loop.
 *
 *   NoScaling       — full Newton step (alpha = 1), zero overhead.
 *   MaxVoltageChange — clamp step so max|dVa| <= max_dVa and max|dVm| <= max_dVm.
 *   LineSearch       — Armijo backtracking line search.
 *   Iwamoto          — Iwamoto optimal multiplier step.
 */
enum class ScalingPolicyType : int {
    Unknown          = -1,
    NoScaling        = 0,
    MaxVoltageChange = 1,
    LineSearch       = 2,
    Iwamoto          = 3
};


template<class NRSystem>
class LS2G_API ScalingPolicy
{
    public:
        virtual ~ScalingPolicy() noexcept = default;
        virtual ScalingPolicyType type() const = 0;
        virtual real_type scale(const NRSystem& system, const RealVect & F) = 0;

    // call to update the policy at each iteration
    // nothing to do in general, is used for 
    // LineSearch and Iwamoto
    // virtual void update(const NRSystem& system, const RealVect & F) {} 
};


template<class NRSystem>
class LS2G_API NoScalingPolicy final : public ScalingPolicy<NRSystem>
{
    public:
        virtual ScalingPolicyType type() const final {return ScalingPolicyType::NoScaling;}
        virtual real_type scale(const NRSystem& system, const RealVect & F) final
        {
            return 1.;
        }
};


template<class NRSystem>
class LS2G_API MaxVoltageChangeScalingPolicy final : public ScalingPolicy<NRSystem>
{
    public:
        virtual ScalingPolicyType type() const final {return ScalingPolicyType::MaxVoltageChange;}
        virtual real_type scale(const NRSystem& system, const RealVect & F) final
        {
            real_type alpha = static_cast<real_type>(1.0);
            if (system.theta_size() > 0) {
                const real_type max_abs_dtheta = system.theta(F).array().abs().maxCoeff();
                if (max_abs_dtheta > max_dVa_)
                    alpha = std::min(alpha, max_dVa_ / max_abs_dtheta);
            }
            if (system.vm_size() > 0) {
                const real_type max_abs_dvm = system.vm(F).array().abs().maxCoeff();
                if (max_abs_dvm > max_dVm_)
                    alpha = std::min(alpha, max_dVm_ / max_abs_dvm);
            }
            return alpha;
        }

        // getter / setters
        real_type max_dVa() const { return max_dVa_;}
        real_type max_dVm() const { return max_dVm_;}

        void max_dVa(real_type val) { max_dVa_ = val;}
        void max_dVm(real_type val) { max_dVm_ = val;}

    private:
        real_type max_dVa_ = static_cast<real_type>(0.5);
        real_type max_dVm_ = static_cast<real_type>(0.1);
};


template<class NRSystem>
class LS2G_API LineSearchScalingPolicy final : public ScalingPolicy<NRSystem>
{
    public:
        virtual ScalingPolicyType type() const final {return ScalingPolicyType::LineSearch;}

        virtual real_type scale(const NRSystem& system, const RealVect & F) final
        {
            // TODO speed Save current merit (||mismatch||^2 before step)
            const real_type F_norm_sq_0 = F.squaredNorm();
            real_type alpha = static_cast<real_type>(1.0);
            for (int k = 0; k < ls_max_iter_; ++k) {
                const real_type threshold =
                    (static_cast<real_type>(1.0)
                        - static_cast<real_type>(2.0) * ls_c_ * alpha) * F_norm_sq_0;
                if (system.mismatch_sq_norm_at(alpha * F) <= threshold) break;
                alpha *= ls_rho_;
            }
            return alpha;
        }

    // virtual void update(const NRSystem& system, const RealVect & F) final {
    //     F_norm_sq_0_ = F.squaredNorm();
    // } 

    // should be called each iteration of NR
    // real_type F_norm_sq_0() const {return F_norm_sq_0_;}
    // void F_norm_sq_0(real_type val) {F_norm_sq_0_ = val;}  // F_norm_sq_0 = F.squaredNorm()

    // getter / setter
    int ls_max_iter() const {return ls_max_iter_;}
    void ls_max_iter(int val) {ls_max_iter_ = val;}
    real_type ls_c() const {return ls_c_;}
    void ls_c(real_type val) {ls_c_ = val;}
    real_type ls_rho() const {return ls_rho_;}
    void ls_rho(real_type val) {ls_rho_ = val;}
    
    private:
        int ls_max_iter_ = 20;
        real_type ls_c_ = static_cast<real_type>(1e-4);
        real_type ls_rho_ = static_cast<real_type>(0.5);

        // needs update at each iteration of NR
        // real_type F_norm_sq_0_ = 0.;
};


template<class NRSystem>
class IwamotoScalingPolicy final : public ScalingPolicy<NRSystem>
{
    public:
        virtual ScalingPolicyType type() const final {return ScalingPolicyType::Iwamoto;}

        virtual real_type scale(const NRSystem& system, const RealVect & F) final
        {
            // TODO speed this is already computed in NR (in previous step)
            real_type g0 = F.squaredNorm();

            real_type g1 = system.mismatch_sq_norm_at(F);
            real_type mu = (g0 + g1 > static_cast<real_type>(0.))
                            ? g0 / (g0 + g1)
                            : static_cast<real_type>(1.0);
            mu = std::max(iw_mu_min_, std::min(iw_mu_max_, mu));
            return mu;
        }

    // virtual void update(const NRSystem& system, const RealVect & F) final {
    //     g0_ = F.squaredNorm();
    //     g1_ = system.mismatch_sq_norm_at(F);
    // } 

    // should be called each iteration of NR
    // real_type g0() const {return g0_;}  //  F.squaredNorm();
    // void g0(real_type val) {g0_ = val;}

    // real_type g1() const {return g1_;}  // _system.mismatch_sq_norm_at(F);
    // void g1(real_type val) {g1_ = val;}

        // getter / setter
        real_type iw_mu_min() const {return iw_mu_min_;}
        void iw_mu_min(real_type val) {iw_mu_min_ = val;}
        real_type iw_mu_max() const {return iw_mu_max_;}
        void iw_mu_max(real_type val) {iw_mu_max_ = val;}

    private:
       real_type iw_mu_min_ = static_cast<real_type>(1e-4);
       real_type iw_mu_max_ = static_cast<real_type>(1.0);

        // needs update at each iteration of NR
        // real_type g0_ = 0.;
        // real_type g1_ = 0.;

};


template <class NRSystem>
std::unique_ptr<ScalingPolicy<NRSystem> > create_scaling_policy(ScalingPolicyType type){
    switch (type)
    {
    case ScalingPolicyType::NoScaling:
        return std::make_unique<NoScalingPolicy<NRSystem> >();
    case ScalingPolicyType::MaxVoltageChange:
        return std::make_unique<MaxVoltageChangeScalingPolicy<NRSystem> >();
    case ScalingPolicyType::LineSearch:
        return std::make_unique<LineSearchScalingPolicy<NRSystem> >();
    case ScalingPolicyType::Iwamoto:
        return std::make_unique<IwamotoScalingPolicy<NRSystem> >();
    default:
        throw std::runtime_error("Unknown scaling policy used");
    }
}

template<class NRSystem>
void update_scaling_policy_params(
    ScalingPolicy<NRSystem> * ptr_basse_pol,
    // MaxVoltageChange
    real_type max_dVa,
    real_type max_dVm,
    // LineSearch (Armijo) params
    real_type ls_c,
    real_type ls_rho,
    int       ls_max_iter,
    // Iwamoto params
    real_type iw_mu_min,
    real_type iw_mu_max
)
{
    switch (ptr_basse_pol->type())
    {
    case ScalingPolicyType::NoScaling:
        break;
    case ScalingPolicyType::MaxVoltageChange:
        {        
            MaxVoltageChangeScalingPolicy<NRSystem> *ptr_pol = dynamic_cast<MaxVoltageChangeScalingPolicy<NRSystem>* >(ptr_basse_pol);
            ptr_pol->max_dVa(max_dVa);
            ptr_pol->max_dVm(max_dVm);
        }
        return;
    case ScalingPolicyType::LineSearch:
        {
            LineSearchScalingPolicy<NRSystem> *ptr_pol = dynamic_cast<LineSearchScalingPolicy<NRSystem>* >(ptr_basse_pol);
            ptr_pol->ls_c(ls_c);
            ptr_pol->ls_max_iter(ls_max_iter);
            ptr_pol->ls_rho(ls_rho);
        }
        return;
    case ScalingPolicyType::Iwamoto:
        {
            IwamotoScalingPolicy<NRSystem> *ptr_pol = dynamic_cast<IwamotoScalingPolicy<NRSystem>* >(ptr_basse_pol);
            ptr_pol->iw_mu_min(iw_mu_min);
            ptr_pol->iw_mu_max(iw_mu_max);
        }
        return;
    default:
        throw std::runtime_error("Unknown scaling policy used");
    }
}

} // namespace ls2g

#endif // SCALING_POLICIES_H
