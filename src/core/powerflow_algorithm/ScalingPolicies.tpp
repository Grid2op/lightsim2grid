// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

// ---- NoScalingPolicy -------------------------------------------------------

template<class SlackPolicyT, class NRAlgoT>
void NoScalingPolicy::apply_step(
    NRAlgoT& algo,
    RealVect& F,
    real_type /*F_norm_sq_0*/,
    real_type& slack_absorbed,
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect& Sbus,
    size_t slack_bus_id,
    const RealVect& slack_weights,
    const Eigen::VectorXi& my_pv,
    Eigen::Ref<const IntVect> pq) const
{
    SlackPolicyT::update_slack_absorbed(algo, F, slack_absorbed);
    algo._do_voltage_update(F, my_pv, pq);
    F = SlackPolicyT::evaluate_Fx(algo, Ybus, algo.V_, Sbus,
                                   slack_bus_id, slack_absorbed, slack_weights,
                                   my_pv, pq);
}

// ---- MaxVoltageChangePolicy ------------------------------------------------

template<class NRAlgoT>
real_type MaxVoltageChangePolicy::_compute_alpha(NRAlgoT& algo, const RealVect& dx) const
{
    real_type alpha = static_cast<real_type>(1.0);
    if(algo._layout.theta_size() > 0){
        const real_type max_abs_dtheta = algo._layout.theta(dx).array().abs().maxCoeff();
        if(max_abs_dtheta > max_dVa_) alpha = std::min(alpha, max_dVa_ / max_abs_dtheta);
    }
    if(algo._layout.vm_size() > 0){
        const real_type max_abs_dvm = algo._layout.vm(dx).array().abs().maxCoeff();
        if(max_abs_dvm > max_dVm_) alpha = std::min(alpha, max_dVm_ / max_abs_dvm);
    }
    return alpha;
}

template<class SlackPolicyT, class NRAlgoT>
void MaxVoltageChangePolicy::apply_step(
    NRAlgoT& algo,
    RealVect& F,
    real_type /*F_norm_sq_0*/,
    real_type& slack_absorbed,
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect& Sbus,
    size_t slack_bus_id,
    const RealVect& slack_weights,
    const Eigen::VectorXi& my_pv,
    Eigen::Ref<const IntVect> pq) const
{
    const real_type alpha = _compute_alpha(algo, F);
    F *= alpha;
    SlackPolicyT::update_slack_absorbed(algo, F, slack_absorbed);
    algo._do_voltage_update(F, my_pv, pq);
    F = SlackPolicyT::evaluate_Fx(algo, Ybus, algo.V_, Sbus,
                                   slack_bus_id, slack_absorbed, slack_weights,
                                   my_pv, pq);
}

// ---- LineSearchPolicy ------------------------------------------------------

template<class NRAlgoT>
CplxVect LineSearchPolicy::_compute_trial_V(
    NRAlgoT& algo,
    const RealVect& dx,
    real_type alpha,
    const Eigen::VectorXi& my_pv,
    Eigen::Ref<const IntVect> pq) const
{
    RealVect Va_trial = algo.Va_;
    RealVect Vm_trial = algo.Vm_;
    if(algo._layout.nb_pv() > 0)
        Va_trial(my_pv) -= alpha * algo._layout.theta(dx).segment(0, algo._layout.nb_pv());
    if(algo._layout.nb_pq() > 0){
        Va_trial(pq) -= alpha * algo._layout.theta(dx).segment(algo._layout.nb_pv(),
                                                                 algo._layout.nb_pq());
        Vm_trial(pq) -= alpha * algo._layout.vm(dx);
    }
    const cplx_type m_i = BaseConstants::my_i;
    return Vm_trial.array() * (Va_trial.array().cos().template cast<cplx_type>()
                               + m_i * Va_trial.array().sin().template cast<cplx_type>());
}

template<class SlackPolicyT, class NRAlgoT>
void LineSearchPolicy::apply_step(
    NRAlgoT& algo,
    RealVect& F,
    real_type F_norm_sq_0,
    real_type& slack_absorbed,
    const Eigen::SparseMatrix<cplx_type>& Ybus,
    const CplxVect& Sbus,
    size_t slack_bus_id,
    const RealVect& slack_weights,
    const Eigen::VectorXi& my_pv,
    Eigen::Ref<const IntVect> pq) const
{
    real_type alpha = static_cast<real_type>(1.0);
    for(int k = 0; k < max_iter_; ++k){
        real_type sa_trial = slack_absorbed;
        RealVect dx_trial = alpha * F;
        SlackPolicyT::update_slack_absorbed(algo, dx_trial, sa_trial);
        CplxVect V_trial = _compute_trial_V(algo, F, alpha, my_pv, pq);
        RealVect F_trial = SlackPolicyT::evaluate_Fx(algo, Ybus, V_trial, Sbus,
                                                      slack_bus_id, sa_trial, slack_weights,
                                                      my_pv, pq);
        const real_type threshold = (static_cast<real_type>(1.0)
                                     - static_cast<real_type>(2.0) * c_ * alpha) * F_norm_sq_0;
        if(F_trial.squaredNorm() <= threshold)
            break;
        alpha *= rho_;
    }
    // commit the accepted step
    F *= alpha;
    SlackPolicyT::update_slack_absorbed(algo, F, slack_absorbed);
    algo._do_voltage_update(F, my_pv, pq);
    F = SlackPolicyT::evaluate_Fx(algo, Ybus, algo.V_, Sbus,
                                   slack_bus_id, slack_absorbed, slack_weights,
                                   my_pv, pq);
}
