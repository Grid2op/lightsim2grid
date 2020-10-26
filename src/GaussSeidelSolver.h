// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef GAUSSSEIDELSOLVER_H
#define GAUSSSEIDELSOLVER_H

#include "BaseSolver.h"

class GaussSeidelSolver : public BaseSolver
{
    public:
        GaussSeidelSolver():BaseSolver() {};

        ~GaussSeidelSolver(){}

        // todo  can be factorized
        Eigen::SparseMatrix<double> get_J(){
            throw std::runtime_error("get_J: There is no jacobian in the Gauss Seidel method");
        }

        // todo change the name!
        bool compute_pf(const Eigen::SparseMatrix<cdouble> & Ybus,
                        Eigen::VectorXcd & V,
                        const Eigen::VectorXcd & Sbus,
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        double tol
                        ) ;


    protected:


        void one_iter_all_at_once(Eigen::VectorXcd & tmp_Sbus,
                                  const Eigen::SparseMatrix<cdouble> & Ybus,
                                  const Eigen::VectorXi & pv,
                                  const Eigen::VectorXi & pq
                                  );
        void one_iter(Eigen::VectorXcd & tmp_Sbus,
                      const Eigen::SparseMatrix<cdouble> & Ybus,
                      const Eigen::VectorXi & pv,
                      const Eigen::VectorXi & pq
                      );

    private:
        // no copy allowed
        GaussSeidelSolver( const GaussSeidelSolver & ) ;
        GaussSeidelSolver & operator=( const GaussSeidelSolver & ) ;

};

#endif // GAUSSSEIDELSOLVER_H
