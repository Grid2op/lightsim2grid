// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GAUSSSEIDEL_ALGO_H
#define GAUSSSEIDEL_ALGO_H

#include "BaseAlgo.h"

class GaussSeidelAlgo : public BaseAlgo
{
    public:
        GaussSeidelAlgo():BaseAlgo(true) {};

        ~GaussSeidelAlgo(){}

        // todo  can be factorized
        Eigen::SparseMatrix<real_type> get_J(){
            throw std::runtime_error("get_J: There is no jacobian in the Gauss Seidel method");
        }

        // todo change the name!
        virtual
        bool compute_pf(const Eigen::SparseMatrix<cplx_type> & Ybus,
                        CplxVect & V,
                        const CplxVect & Sbus,
                        const Eigen::VectorXi & slack_ids,
                        const RealVect & slack_weights,  // currently unused
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        real_type tol
                        ) ;

    protected:

        virtual
        void one_iter(CplxVect & tmp_Sbus,
                      const Eigen::SparseMatrix<cplx_type> & Ybus,
                      const Eigen::VectorXi & pv,
                      const Eigen::VectorXi & pq
                      );

    private:
        // no copy allowed
        GaussSeidelAlgo( const GaussSeidelAlgo & ) =delete;
        GaussSeidelAlgo & operator=( const GaussSeidelAlgo & ) =delete;

};

#endif // GAUSSSEIDEL_ALGO_H
