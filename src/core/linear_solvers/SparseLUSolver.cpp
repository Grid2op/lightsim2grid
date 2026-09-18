// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SparseLUSolver.hpp"

#include <iostream>

namespace ls2g {

ErrorType SparseLULinearSolver::analyze(const EigenRefConstRealSpMat & J){
    solver_.analyzePattern(J);
    // analyzePattern does not set solver_.info() to Success, so no check here
    return ErrorType::NoError;
}

ErrorType SparseLULinearSolver::factorize(const EigenRefConstRealSpMat & J){
    solver_.factorize(J);
    if(solver_.info() != Eigen::Success) return ErrorType::SolverFactor;
    return ErrorType::NoError;
}

ErrorType SparseLULinearSolver::refactorize(const EigenRefConstRealSpMat & J){
    solver_.factorize(J);
    if(solver_.info() != Eigen::Success) return ErrorType::SolverFactor;
    return ErrorType::NoError;
}

ErrorType SparseLULinearSolver::solve_transpose(Eigen::Ref<RealVect> b){
    ErrorType err = ErrorType::NoError;
    // Eigen reuses the L / U factors of J and walks them the other way round
    // (transpose(), not adjoint(): both agree here since real_type is real, but
    // adjoint() would conjugate if this were ever instantiated on a complex type).
    RealVect Va = solver_.transpose().solve(b);
    if(solver_.info() != Eigen::Success) err = ErrorType::SolverSolve;
    b = Va;
    return err;
}

ErrorType SparseLULinearSolver::solve(Eigen::Ref<RealVect> b) const{
    ErrorType err = ErrorType::NoError;
    RealVect Va = solver_.solve(b);
    // std::cout << "\t\tSparseLUSolver.cpp: solver_.info: " << solver_.info() << std::endl;  // TODO DEBUG WINDOWS
    if(solver_.info() != Eigen::Success) err = ErrorType::SolverSolve;
    b = Va;
    return err;
}

} // namespace ls2g
