// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NICSLUSolver.h"
#include <iostream>

const bool NICSLULinearSolver::CAN_SOLVE_MAT = false;

ErrorType NICSLULinearSolver::reset(){
    // free everything
    solver_.Free();
    if(ai_ != nullptr) delete [] ai_;  // created in NICSLUSolver::initialize
    if(ap_ != nullptr) delete [] ap_;  // created in NICSLUSolver::initialize
    ai_ = nullptr;
    ap_ = nullptr;

    solver_ = CNicsLU();
    int ret = solver_.Initialize();
    if (ret < 0){   // fail
        if (ret == -100){
            // fail to initialize because of a license issue
            std::string msg = "Fail to initilize the NICSLU solver because we cannot find the nicslu.lic file.";
            msg += "Please copy this file at the location you want to use the NICSLU solver.";
            std::cout << msg << std::endl;
        }
        // std::cout << "NICSLULinearSolver::reset() solver_.Initialize error: "  << ret << std::endl;
        return ErrorType::LicenseError;
    }
    // solver_.SetConfiguration(0, 1.); //enable timer
    return ErrorType::NoError;
}

ErrorType NICSLULinearSolver::initialize(Eigen::SparseMatrix<real_type> & J){
    // default Eigen representation: column major, which is good for klu !
    const auto n = J.cols(); // should be equal to J_.nrows()
    int ret;
    ErrorType err = ErrorType::NoError; // reset error message
    const unsigned int nnz = J.nonZeros();

    ai_ = new unsigned int [nnz];  // deleted in destructor and NICSLUSolver::reset
    const int * ref_ai = J.innerIndexPtr();
    for(unsigned int i = 0; i < nnz; ++i){
        ai_[i] = static_cast<unsigned int>(ref_ai[i]);
    }

    ap_ = new unsigned int [n+1];  // deleted in destructor and NICSLUSolver::reset
    const int * ref_ap = J.outerIndexPtr();
    for(int i = 0; i < n+1; ++i){
        ap_[i] = static_cast<unsigned int>(ref_ap[i]);
    }
    ret = solver_.Analyze(n,
                          J.valuePtr(),
                          ai_,
                          ap_,
                          MATRIX_COLUMN_REAL);  // MATRIX_COLUMN_REAL, MATRIX_ROW_REAL
    if (ret < 0){
        err = ErrorType::SolverAnalyze;
        // std::cout << "NICSLULinearSolver::initialize() solver_.Analyze error: "  << ret << std::endl;
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return err;
    }
    // solver.FactorizeMatrix(ax, 0); //use all created threads
    ret = solver_.FactorizeMatrix(J.valuePtr(), nb_thread_);
    if (ret < 0){
        err = ErrorType::SolverFactor;
        // std::cout << "NICSLULinearSolver::initialize() solver_.FactorizeMatrix error: "  << ret << std::endl;
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return err;
    }
    return err;
}

ErrorType NICSLULinearSolver::solve(Eigen::SparseMatrix<real_type> & J, RealVect & b, bool doesnt_need_refactor){
    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call klu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    int ret;
    bool stop = false;
    RealVect x;
    ErrorType err = ErrorType::NoError;
    if(!doesnt_need_refactor){
        ret  = solver_.FactorizeMatrix(J.valuePtr(), nb_thread_);  // TODO maybe 0 instead of nb_thread_ here, see https://github.com/chenxm1986/nicslu/blob/master/nicslu202110/demo/demo2.cpp
        if (ret < 0) {
            // std::cout << "NICSLULinearSolver::solve solver_.FactorizeMatrix error: " << ret << std::endl;
            // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
            err = ErrorType::SolverReFactor;
            stop = true;
        }
    }
    if(!stop){
        const auto n = J.cols(); // should be equal to J_.nrows()
        x = RealVect(n);
        ret = solver_.Solve(&b(0), &x(0));
        if (ret < 0) {
            // std::cout << "NICSLULinearSolver::solve solver_.Solve error: " << ret << std::endl;
            err = ErrorType::SolverSolve;
        }
        b = x;
    }
    return err;
}
