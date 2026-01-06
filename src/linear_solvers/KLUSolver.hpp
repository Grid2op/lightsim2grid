// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef KLU_SOLVER_AVAILABLE
#ifndef KLSOLVER_H
#define KLSOLVER_H
#include <iostream>

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.hpp"

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

// import klu package
extern "C" {
    #include "cs.h"
    #include "klu.h"
}

/**
class to handle the solver using newton-raphson method, using KLU algorithm and sparse matrices.

As long as the admittance matrix of the sytem does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

**/
class KLULinearSolver final
{
    public:
        KLULinearSolver() noexcept :symbolic_(nullptr),numeric_(nullptr),common_(){}

        ~KLULinearSolver() noexcept
        {
            // std::cout << "KLULinearSolver destructor 1" << std::endl;
            if(symbolic_ != nullptr) klu_free_symbolic(&symbolic_, &common_);
            // std::cout << "KLULinearSolver destructor 2" << std::endl;
            if(numeric_ != nullptr) klu_free_numeric(&numeric_, &common_);
            // std::cout << "KLULinearSolver destructor 3" << std::endl;
        }

        // KLULinearSolver(KLULinearSolver && other) noexcept {
        //     TODO
        //     if(symbolic_ != nullptr) klu_free_symbolic(&symbolic_, &common_);
        //     symbolic_ = other.symbolic_;
        //     other.symbolic_ = nullptr;
            
        //     if(numeric_ != nullptr) klu_free_numeric(&numeric_, &common_);
        //     numeric_ = other.numeric_;
        //     other.numeric_ = nullptr;

        //     std::swap(common_, other.common_);
        // }

        // public api
        ErrorType reset();
        ErrorType initialize(const Eigen::SparseMatrix<real_type>& J);
        ErrorType solve(const Eigen::SparseMatrix<real_type>& J, RealVect & b, bool doesnt_need_refactor);

        // can this linear solver solve problem where RHS is a matrix
        static const bool CAN_SOLVE_MAT;
        
    private:
        // solver initialization
        klu_symbolic* symbolic_;
        klu_numeric* numeric_;
        klu_common common_;

        // no copy allowed
        KLULinearSolver(const KLULinearSolver&) = delete;
        KLULinearSolver(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(KLULinearSolver&&) = delete;
        KLULinearSolver & operator=(const KLULinearSolver&) = delete;
};

#endif // KLSOLVER_H
#endif // KLU_SOLVER_AVAILABLE
