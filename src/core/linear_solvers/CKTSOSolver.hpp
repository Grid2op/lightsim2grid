// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef CKTSO_SOLVER_AVAILABLE
#ifndef CKTSOSOLVER_H
#define CKTSOSOLVER_H

#include <memory>

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.hpp"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

// import cktso package
#include "cktso.h"

namespace ls2g {

/**
class to handle the solver using newton-raphson method, using CKTSO algorithm and sparse matrices.
CKTSO, according to https://github.com/Grid2Op/lightsim2grid/issues/52
is the successor of NICSLU.

As long as the admittance matrix of the system does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

NB: the code of CKTSO is not included in this repository. This class is only compiled if the "setup.py"
can find a version of `https://github.com/chenxm1986/cktso`. Be careful though, this code is under some
specific license.

**/
class LS2G_API CKTSOLinearSolver final
{
    public:
        CKTSOLinearSolver() noexcept:
            solver_(nullptr),
            nb_thread_(1),
            ai_(nullptr), 
            ap_(nullptr),
            iparm_(nullptr),
            oparm_(nullptr)
            {}

        ~CKTSOLinearSolver() noexcept
        {
           if(solver_!= nullptr) solver_->DestroySolver();
           // ai_ / ap_ are std::unique_ptr and free themselves.

           // should not be deleted, see https://github.com/Grid2Op/lightsim2grid/issues/52#issuecomment-1333565959
           // if(iparm_!= nullptr) delete iparm_;
           // if(oparm_!= nullptr) delete oparm_;

        }
        // CKTSOLinearSolver(CKTSOLinearSolver && other) noexcept: nb_thread_(ohter.nb_thread_){
        // TODO !
        //     std::swap(solver_, other.solver_);

        //     if(ai_!= nullptr) delete [] ai_;
        //     ai_ = other.ai_;
        //     other.ai_ = nullptr;

        //     if(ap_!= nullptr) delete [] ap_;
        //     ap_ = other.ap_;
        //     other.ap_ = nullptr;

        //     iparm_ = other.iparm_;
        //     other.iparm_ = nullptr;

        //     oparm_ = other.oparm_;
        //     other.oparm_ = nullptr;
        // }
        

        // public api
        ErrorType reset();
        ErrorType analyze(const Eigen::Ref<const Eigen::SparseMatrix<real_type>> & J);   // creates solver + reordering + symbolic factorization
        ErrorType factorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type>> & J); // numeric factorization (requires values)
        ErrorType refactorize(const Eigen::Ref<const Eigen::SparseMatrix<real_type>> & J);  // re-numeric factorization, reuses symbolic
        ErrorType solve(RealVect & b);

        // can this linear solver solve problem where RHS is a matrix
        static const bool CAN_SOLVE_MAT;
        
        // prevent copy and assignment
        CKTSOLinearSolver(const CKTSOLinearSolver&) = delete;
        CKTSOLinearSolver(CKTSOLinearSolver&&) = delete;
        CKTSOLinearSolver & operator=(CKTSOLinearSolver&&) = delete;
        CKTSOLinearSolver & operator=(const CKTSOLinearSolver&) = delete;
        
    private:
        // solver initialization
        ICktSo solver_;
        const unsigned int nb_thread_;
        // ai_ / ap_ own the CSC index buffers passed to CKTSO's Analyze; they must
        // stay alive until reset() / destruction, which std::unique_ptr handles.
        std::unique_ptr<int[]> ai_;
        std::unique_ptr<int[]> ap_;
        int * iparm_;   // owned by CKTSO, not freed here
        long long * oparm_;  // owned by CKTSO, not freed here

};

#endif // CKTSOSOLVER_H

} // namespace ls2g

#endif  // CKTSO_SOLVER_AVAILABLE