// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef CKTSO_SOLVER_AVAILABLE
#ifndef CKTSOSOLVER_H
#define CKTSOSOLVER_H

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

// import cktso package
#include "cktso.h"

/**
class to handle the solver using newton-raphson method, using CKTSO algorithm and sparse matrices.
CKTSO, according to https://github.com/BDonnot/lightsim2grid/issues/52
is the successor of NICSLU.

As long as the admittance matrix of the system does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

NB: the code of CKTSO is not included in this repository. This class is only compiled if the "setup.py"
can find a version of `https://github.com/chenxm1986/cktso`. Be careful though, this code is under some
specific license.

**/
class CKTSOLinearSolver
{
    public:
        CKTSOLinearSolver():
            solver_(nullptr),
            nb_thread_(1),
            ai_(nullptr), 
            ap_(nullptr),
            iparm_(nullptr),
            oparm_(nullptr)
            {}

        ~CKTSOLinearSolver()
         {
            if(solver_!= nullptr) solver_->DestroySolver();
            if(ai_!= nullptr) delete [] ai_;
            if(ap_!= nullptr) delete [] ap_;
            
            // should not be deleted, see https://github.com/BDonnot/lightsim2grid/issues/52#issuecomment-1333565959
            // if(iparm_!= nullptr) delete iparm_;
            // if(oparm_!= nullptr) delete oparm_;
            
         }

        // public api
        ErrorType reset();
        ErrorType initialize(Eigen::SparseMatrix<real_type> & J);
        ErrorType solve(Eigen::SparseMatrix<real_type> & J, RealVect & b, bool doesnt_need_refactor);

        // can this linear solver solve problem where RHS is a matrix
        static const bool CAN_SOLVE_MAT;
        
        // prevent copy and assignment
        CKTSOLinearSolver(const CKTSOLinearSolver & other) = delete;
        CKTSOLinearSolver & operator=( const CKTSOLinearSolver & ) = delete;
        
    private:
        // solver initialization
        ICktSo solver_;
        const unsigned int nb_thread_;
        int * ai_;
        int * ap_;
        int * iparm_;
        long long * oparm_;

};

#endif // CKTSOSOLVER_H
#endif  // CKTSO_SOLVER_AVAILABLE
