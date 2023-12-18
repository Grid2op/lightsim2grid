// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASE_NR_SINGLESLACK_ALGO_H
#define BASE_NR_SINGLESLACK_ALGO_H

#include "BaseNRAlgo.h"

/**
Base class for Newton Raphson based solver (only interesting for single slack)
**/
template<class LinearSolver>
class BaseNRSingleSlackAlgo : public BaseNRAlgo<LinearSolver>
{
    public:
        BaseNRSingleSlackAlgo():BaseNRAlgo<LinearSolver>(){}

        ~BaseNRSingleSlackAlgo(){}

        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                        CplxVect & V,
                        const CplxVect & Sbus,
                        const Eigen::VectorXi & slack_ids,
                        const RealVect & slack_weights,
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        );


    protected:
        void fill_jacobian_matrix(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                  const CplxVect & V,
                                  const Eigen::VectorXi & pq,
                                  const Eigen::VectorXi & pvpq,
                                  const std::vector<int> & pq_inv,
                                  const std::vector<int> & pvpq_inv);

        void fill_jacobian_matrix_kown_sparsity_pattern(
                 const Eigen::VectorXi & pq,
                 const Eigen::VectorXi & pvpq
                 );

        void fill_jacobian_matrix_unkown_sparsity_pattern(
                 const Eigen::SparseMatrix<cplx_type> & Ybus,
                 const CplxVect & V,
                 const Eigen::VectorXi & pq,
                 const Eigen::VectorXi & pvpq,
                 const std::vector<int> & pq_inv,
                 const std::vector<int> & pvpq_inv
                 );

        void fill_value_map(const Eigen::VectorXi & pq,
                            const Eigen::VectorXi & pvpq,
                            bool reset_J);

};

#include "BaseNRSingleSlackAlgo.tpp"

#endif // BASE_NR_SINGLESLACK_ALGO_H
