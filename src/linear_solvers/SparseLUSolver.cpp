// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SparseLUSolver.hpp"

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

ErrorType SparseLULinearSolver::refactor(const Eigen::SparseMatrix<real_type> & J){
    solver_.factorize(J);
    if(solver_.info() != Eigen::Success) return ErrorType::SolverFactor;
    return ErrorType::NoError;
}

ErrorType SparseLULinearSolver::solve(RealVect & b){
    ErrorType err = ErrorType::NoError;
    RealVect Va = solver_.solve(b);
    // std::cout << "\t\tSparseLUSolver.cpp: solver_.info: " << solver_.info() << std::endl;  // TODO DEBUG WINDOWS
    if(solver_.info() != Eigen::Success) err = ErrorType::SolverSolve;
    b = Va;
    return err;
}
