// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// This file has been inspired from https://github.com/chenxm1986/cktso/blob/master/demo/demo.cpp
#include "CKTSOSolver.h"
#include <iostream>

const bool CKTSOLinearSolver::CAN_SOLVE_MAT = false;


ErrorType CKTSOLinearSolver::reset(){
    // free everything
    if(solver_ != nullptr) solver_->DestroySolver();
    if(ai_ != nullptr) delete [] ai_;
    if(ap_ != nullptr) delete [] ap_;

    // should not be deleted, see https://github.com/BDonnot/lightsim2grid/issues/52#issuecomment-1333565959
    // if(iparm_!= nullptr) delete iparm_;
    // if(oparm_!= nullptr) delete oparm_;

    ai_ = nullptr;
    ap_ = nullptr;
    iparm_ = nullptr;
    oparm_ = nullptr;

    solver_ = nullptr;
    return ErrorType::NoError;
}

ErrorType CKTSOLinearSolver::initialize(Eigen::SparseMatrix<real_type> & J){
    const long long *oparm = oparm_;
    int ret_ = CKTSO_CreateSolver(&solver_, &iparm_, &oparm);
    if (ret_ < 0){   // fail
        if (ret_ == -99){
            // fail to initialize because of a license issue
            std::string msg = "Fail to initilize the CKTSO solver because we cannot find the cktso.lic file. ";
            msg += "Please copy this file at the location you want to use the CKTSO solver.";
            std::cout << msg << std::endl;
        }
        return ErrorType::LicenseError;
    }

    const auto n = J.cols(); // should be equal to J_.nrows()
    int ret;
    ErrorType err = ErrorType::NoError; // reset error message
    const unsigned int nnz = J.nonZeros();

    ai_ = new int [nnz];
    const int * ref_ai = J.innerIndexPtr();
    for(unsigned int i = 0; i < nnz; ++i){
        ai_[i] = static_cast<int>(ref_ai[i]);
    }

    ap_ = new int [n+1];
    const int * ref_ap = J.outerIndexPtr();
    for(int i = 0; i < n+1; ++i){
        ap_[i] = static_cast<int>(ref_ap[i]);
    }
    ret = solver_->Analyze(false,  // complex or real
                           n,
                           ap_,
                           ai_,
                           J.valuePtr(),
                           nb_thread_);
    if (ret < 0){
        err = ErrorType::SolverAnalyze;
        // std::cout << "CKTSOLinearSolver::initialize() solver_->Analyze error: "  << ret << std::endl;
        // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
        return err;
    }

    ret = solver_->Factorize(J.valuePtr(),
                             true // @fast: whether to use fast factorization => was set to "false" in the demo
                             );
    if (ret < 0){
        err = ErrorType::SolverFactor;
        // std::cout << "CKTSOLinearSolver::initialize() solver_->Factorize error: "  << ret << std::endl;
        // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
        return err;
    }
    return err;
}

ErrorType CKTSOLinearSolver::solve(Eigen::SparseMatrix<real_type> & J, RealVect & b, bool doesnt_need_refactor){
    // solves (for x) the linear system J.x = b
    // with standard use of lightsim2grid, the solver should have already been initialized
    // J is const even if it does not compile if said const
    int ret;
    bool stop = false;
    RealVect x;
    ErrorType err = ErrorType::NoError;
    if(!doesnt_need_refactor){
        ret  = solver_->Refactorize(J.valuePtr()); 
        if (ret < 0) {
            // std::cout << "CKTSOLinearSolver::solve solver_->Refactorize error: " << ret << std::endl;
            // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
            err = ErrorType::SolverReFactor;
            stop = true;
        }
    }
    if(!stop){
        const auto n = J.cols(); // should be equal to J_.nrows()
        x = RealVect(n);
        ret = solver_->Solve(&b(0), &x(0), false, 1);
        if (ret < 0) {
            // std::cout << "CKTSOLinearSolver::solve solver_.Solve error: " << ret << std::endl;
            err = ErrorType::SolverSolve;
        }
        b = x;
    }
    return err;
}
