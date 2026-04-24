// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef SLACK_POLICIES_H
#define SLACK_POLICIES_H

#include "Utils.hpp"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#include <vector>

namespace ls2g {

/**
 * Policy for multi-slack Newton-Raphson: adds a distributed-slack row/col to the Jacobian.
 * lag == 1 means J has an extra first row and column for the slack equation.
 */
struct MultiSlackPolicy {
    static constexpr int lag = 1;

    template<class NRAlgoT>
    static Eigen::VectorXi get_my_pv(NRAlgoT& algo,
                                      Eigen::Ref<const IntVect> slack_ids,
                                      Eigen::Ref<const IntVect> pv);

    static real_type initial_slack_absorbed(const CplxVect& Sbus);

    template<class NRAlgoT>
    static RealVect evaluate_Fx(NRAlgoT& algo,
                                 const Eigen::SparseMatrix<cplx_type>& Ybus,
                                 const CplxVect& V,
                                 const CplxVect& Sbus,
                                 size_t slack_bus_id,
                                 real_type slack_absorbed,
                                 const RealVect& slack_weights,
                                 const Eigen::VectorXi& my_pv,
                                 Eigen::Ref<const IntVect> pq);

    template<class NRAlgoT>
    static void update_slack_absorbed(NRAlgoT& algo,
                                       const RealVect& dx,
                                       real_type& slack_absorbed);

    template<class NRAlgoT>
    static void fill_jacobian_matrix(NRAlgoT& algo,
                                      const Eigen::SparseMatrix<cplx_type>& Ybus,
                                      const CplxVect& V,
                                      size_t slack_bus_id,
                                      const RealVect& slack_weights,
                                      const Eigen::VectorXi& pq,
                                      const Eigen::VectorXi& pvpq,
                                      const std::vector<int>& pq_inv,
                                      const std::vector<int>& pvpq_inv);

private:
    template<class NRAlgoT>
    static void fill_jac_unknown_sparsity(NRAlgoT& algo,
                                           const Eigen::SparseMatrix<cplx_type>& Ybus,
                                           const CplxVect& V,
                                           size_t slack_bus_id,
                                           const RealVect& slack_weights,
                                           const Eigen::VectorXi& pq,
                                           const Eigen::VectorXi& pvpq,
                                           const std::vector<int>& pq_inv,
                                           const std::vector<int>& pvpq_inv);

    template<class NRAlgoT>
    static void fill_jac_known_sparsity(NRAlgoT& algo,
                                         size_t slack_bus_id,
                                         const Eigen::VectorXi& pq,
                                         const Eigen::VectorXi& pvpq);

    template<class NRAlgoT>
    static void fill_value_map_impl(NRAlgoT& algo,
                                     size_t slack_bus_id,
                                     const Eigen::VectorXi& pq,
                                     const Eigen::VectorXi& pvpq,
                                     bool reset_J);
};

/**
 * Policy for single-slack Newton-Raphson: no distributed-slack row/col in the Jacobian.
 * lag == 0 means J is exactly (n_pvpq + n_pq) x (n_pvpq + n_pq).
 */
struct SingleSlackPolicy {
    static constexpr int lag = 0;

    template<class NRAlgoT>
    static Eigen::VectorXi get_my_pv(NRAlgoT& algo,
                                      Eigen::Ref<const IntVect> slack_ids,
                                      Eigen::Ref<const IntVect> pv);

    static real_type initial_slack_absorbed(const CplxVect& Sbus);

    template<class NRAlgoT>
    static RealVect evaluate_Fx(NRAlgoT& algo,
                                 const Eigen::SparseMatrix<cplx_type>& Ybus,
                                 const CplxVect& V,
                                 const CplxVect& Sbus,
                                 size_t slack_bus_id,
                                 real_type slack_absorbed,
                                 const RealVect& slack_weights,
                                 const Eigen::VectorXi& my_pv,
                                 Eigen::Ref<const IntVect> pq);

    template<class NRAlgoT>
    static void update_slack_absorbed(NRAlgoT& algo,
                                       const RealVect& dx,
                                       real_type& slack_absorbed);

    template<class NRAlgoT>
    static void fill_jacobian_matrix(NRAlgoT& algo,
                                      const Eigen::SparseMatrix<cplx_type>& Ybus,
                                      const CplxVect& V,
                                      size_t slack_bus_id,
                                      const RealVect& slack_weights,
                                      const Eigen::VectorXi& pq,
                                      const Eigen::VectorXi& pvpq,
                                      const std::vector<int>& pq_inv,
                                      const std::vector<int>& pvpq_inv);

private:
    template<class NRAlgoT>
    static void fill_jac_unknown_sparsity(NRAlgoT& algo,
                                           const Eigen::SparseMatrix<cplx_type>& Ybus,
                                           const CplxVect& V,
                                           const Eigen::VectorXi& pq,
                                           const Eigen::VectorXi& pvpq,
                                           const std::vector<int>& pq_inv,
                                           const std::vector<int>& pvpq_inv);

    template<class NRAlgoT>
    static void fill_jac_known_sparsity(NRAlgoT& algo,
                                         const Eigen::VectorXi& pq,
                                         const Eigen::VectorXi& pvpq);

    template<class NRAlgoT>
    static void fill_value_map_impl(NRAlgoT& algo,
                                     const Eigen::VectorXi& pq,
                                     const Eigen::VectorXi& pvpq,
                                     bool reset_J);
};

#include "SlackPolicies.tpp"

} // namespace ls2g

#endif // SLACK_POLICIES_H
