// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef DCSOLVER_H
#define DCSOLVER_H

#include "BaseSolver.h"
// TODO make err_ more explicit: use an enum
template<class LinearSolver>
class BaseDCSolver: public BaseSolver
{
    public:
        BaseDCSolver():BaseSolver(false), _linear_solver(), need_factorize_(true){};

        ~BaseDCSolver(){}

        virtual void reset();

        // TODO SLACK : this should be handled in Sbus by the gridmodel maybe ?
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
                        );

    private:
        // no copy allowed
        BaseDCSolver( const BaseSolver & ) =delete ;
        BaseDCSolver & operator=( const BaseSolver & ) =delete;

    protected:
        LinearSolver  _linear_solver;
        bool need_factorize_;

};

#include "DCSolver.tpp"

#endif // DCSOLVER_H
