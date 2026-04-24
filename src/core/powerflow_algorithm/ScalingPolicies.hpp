// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SCALING_POLICIES_H
#define SCALING_POLICIES_H

#include "Utils.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <algorithm>

namespace ls2g {

/**
 * Identity scaling: applies the full Newton step without modification.
 * Used as the default ScalingPolicy; inlines to zero overhead.
 */
struct NoScalingPolicy {
    template<class SlackPolicyT, class NRAlgoT>
    void apply_step(NRAlgoT& algo,
                    RealVect& F,
                    real_type F_norm_sq_0,
                    real_type& slack_absorbed,
                    const Eigen::SparseMatrix<cplx_type>& Ybus,
                    const CplxVect& Sbus,
                    size_t slack_bus_id,
                    const RealVect& slack_weights,
                    const Eigen::VectorXi& my_pv,
                    Eigen::Ref<const IntVect> pq) const;
};

/**
 * Limits the Newton step so that no voltage angle changes by more than max_dVa
 * radians and no voltage magnitude changes by more than max_dVm pu in one
 * iteration.
 *
 * Both limits are instance members, settable at runtime from C++ or Python.
 */
class MaxVoltageChangePolicy {
public:
    MaxVoltageChangePolicy() noexcept
        : max_dVa_(static_cast<real_type>(0.5)),
          max_dVm_(static_cast<real_type>(0.1)) {}

    real_type get_max_dVa() const { return max_dVa_; }
    void      set_max_dVa(real_type v) { max_dVa_ = v; }
    real_type get_max_dVm() const { return max_dVm_; }
    void      set_max_dVm(real_type v) { max_dVm_ = v; }

    template<class SlackPolicyT, class NRAlgoT>
    void apply_step(NRAlgoT& algo,
                    RealVect& F,
                    real_type F_norm_sq_0,
                    real_type& slack_absorbed,
                    const Eigen::SparseMatrix<cplx_type>& Ybus,
                    const CplxVect& Sbus,
                    size_t slack_bus_id,
                    const RealVect& slack_weights,
                    const Eigen::VectorXi& my_pv,
                    Eigen::Ref<const IntVect> pq) const;

private:
    real_type max_dVa_;
    real_type max_dVm_;

    template<class NRAlgoT>
    real_type _compute_alpha(NRAlgoT& algo, const RealVect& dx) const;
};

/**
 * Armijo backtracking line search.
 *
 * At each Newton iteration, the step alpha is halved (by factor rho) until the
 * sufficient-decrease condition ||F(x + alpha*dx)||^2 <= (1 - 2*c*alpha)*||F(x)||^2
 * is satisfied or max_iter line-search steps have been taken.
 *
 * All three parameters are settable at runtime from C++ or Python.
 *   c       : Armijo constant      (default 1e-4, range (0,0.5))
 *   rho     : backtracking factor  (default 0.5,  range (0,1))
 *   max_iter: max LS iterations    (default 10)
 */
class LineSearchPolicy {
public:
    LineSearchPolicy() noexcept
        : c_(static_cast<real_type>(1e-4)),
          rho_(static_cast<real_type>(0.5)),
          max_iter_(10) {}

    real_type get_c()    const { return c_; }
    void      set_c(real_type v) { c_ = v; }
    real_type get_rho()  const { return rho_; }
    void      set_rho(real_type v) { rho_ = v; }
    int       get_max_iter() const { return max_iter_; }
    void      set_max_iter(int v)  { max_iter_ = v; }

    template<class SlackPolicyT, class NRAlgoT>
    void apply_step(NRAlgoT& algo,
                    RealVect& F,
                    real_type F_norm_sq_0,
                    real_type& slack_absorbed,
                    const Eigen::SparseMatrix<cplx_type>& Ybus,
                    const CplxVect& Sbus,
                    size_t slack_bus_id,
                    const RealVect& slack_weights,
                    const Eigen::VectorXi& my_pv,
                    Eigen::Ref<const IntVect> pq) const;

private:
    real_type c_;
    real_type rho_;
    int       max_iter_;

    template<class NRAlgoT>
    CplxVect _compute_trial_V(NRAlgoT& algo,
                               const RealVect& dx,
                               real_type alpha,
                               const Eigen::VectorXi& my_pv,
                               Eigen::Ref<const IntVect> pq) const;
};

#include "ScalingPolicies.tpp"

} // namespace ls2g

#endif // SCALING_POLICIES_H
