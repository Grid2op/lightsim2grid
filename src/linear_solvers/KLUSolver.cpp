// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "KLUSolver.h"

#include <iostream>

const bool KLULinearSolver::CAN_SOLVE_MAT = false;

ErrorType KLULinearSolver::reset(){
    klu_free_symbolic(&symbolic_, &common_);
    klu_free_numeric(&numeric_, &common_);
    common_ = klu_common();
    symbolic_ = nullptr;
    numeric_ = nullptr;
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::initialize(Eigen::SparseMatrix<real_type>&  J){
    // default Eigen representation: column major, which is good for klu !
    // J is const here, even if it's not said in klu_analyze
    const auto n = J.cols();
    common_ = klu_common();
    ErrorType res = ErrorType::NoError; 
    symbolic_ = klu_analyze(n, J.outerIndexPtr(), J.innerIndexPtr(), &common_);
    if(common_.status != KLU_OK){
        res = ErrorType::SolverAnalyze; 
    }else{
        numeric_ = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, &common_);
        if(common_.status != KLU_OK) res = ErrorType::SolverFactor; 
    }
    return res;
}

ErrorType KLULinearSolver::solve(Eigen::SparseMatrix<real_type>& J, RealVect & b, bool doesnt_need_refactor){
    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call klu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    int ok;
    ErrorType err = ErrorType::NoError;
    bool stop = false;
    if(!doesnt_need_refactor){
        // if the call to "klu_factor" has been made this iteration, there is no need
        // to re factor again the matrix
        // i'm in the case where it has not
        ok = klu_refactor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, numeric_, &common_);
        if (ok != 1) {
            // std::cout << "\t KLU: refactor error" << std::endl;
            err = ErrorType::SolverReFactor;
            stop = true;
        }
    }
    if(!stop){
        const auto n = J.cols();
        ok = klu_solve(symbolic_, numeric_, n, 1, &b(0), &common_);
        if (ok != 1) {
            // std::cout << "\t KLU: klu_solve error" << std::endl;
            err = ErrorType::SolverSolve;
        }
    }
    return err;
}
