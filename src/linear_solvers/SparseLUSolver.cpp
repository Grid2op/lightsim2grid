// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SparseLUSolver.h"

#include <iostream>

const bool SparseLULinearSolver::CAN_SOLVE_MAT = true;

ErrorType SparseLULinearSolver::initialize(const Eigen::SparseMatrix<real_type> & J){
    // default Eigen representation: column major, which is good for klu !
    // J is const here
    ErrorType res = ErrorType::NoError; 
    solver_.analyzePattern(J);
    // do not check here for "solver_.info" it is not set to "Success"
    solver_.factorize(J);
    if(solver_.info() != Eigen::Success) res = ErrorType::SolverFactor; 
    return res;
}

ErrorType SparseLULinearSolver::solve(const Eigen::SparseMatrix<real_type> & J, RealVect & b, bool doesnt_need_refactor){
    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call sparselu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    ErrorType err = ErrorType::NoError;
    bool stop = false;
    if(!doesnt_need_refactor){
        // if the call to "klu_factor" has been made this iteration, there is no need
        // to re factor again the matrix
        // i'm in the case where it has not
        solver_.factorize(J);
        if (solver_.info() != Eigen::Success) {
            err = ErrorType::SolverFactor;
            stop = true;
        }
    }
    if(!stop){
        RealVect Va = solver_.solve(b);
        if (solver_.info() != Eigen::Success) {
            err = ErrorType::SolverSolve;
        }
        b = Va;
    }
    return err;
}
