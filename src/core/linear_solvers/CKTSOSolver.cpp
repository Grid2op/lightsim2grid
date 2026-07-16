// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


// This file has been inspired from https://github.com/chenxm1986/cktso/blob/master/demo/demo.cpp
#include "CKTSOSolver.hpp"

#include <iostream>

namespace ls2g {

const bool CKTSOLinearSolver::CAN_SOLVE_MAT = false;


ErrorType CKTSOLinearSolver::reset(){
    // free everything
    if(solver_ != nullptr) solver_->DestroySolver();
    ai_.reset();
    ap_.reset();

    // should not be deleted, see https://github.com/Grid2Op/lightsim2grid/issues/52#issuecomment-1333565959
    // if(iparm_!= nullptr) delete iparm_;
    // if(oparm_!= nullptr) delete oparm_;

    iparm_ = nullptr;
    oparm_ = nullptr;

    solver_ = nullptr;
    return ErrorType::NoError;
}

ErrorType CKTSOLinearSolver::analyze(const EigenRefConstRealSpMat & J){
    const long long *oparm = oparm_;
    int ret_ = CKTSO_CreateSolver(&solver_, &iparm_, &oparm);
    if (ret_ < 0){
        if (ret_ == -99){
            std::string msg = "Fail to initilize the CKTSO solver because we cannot find the cktso.lic file. ";
            msg += "Please copy this file at the location you want to use the CKTSO solver (the place where the lib is located). ";
            msg += "See `import lightsim2grid; print(lightsim2grid.compilation_options.cktso_lib)` and then copy paste ";
            msg += "the cktso/license/cktso.lib there. ";
            std::cout << msg << std::endl;
        }
        return ErrorType::LicenseError;
    }

    const auto n = J.cols();
    const unsigned int nnz = J.nonZeros();

    ai_ = std::make_unique<int[]>(nnz);
    const int * ref_ai = J.innerIndexPtr();
    for(unsigned int i = 0; i < nnz; ++i){
        ai_[i] = static_cast<int>(ref_ai[i]);
    }

    ap_ = std::make_unique<int[]>(n+1);
    const int * ref_ap = J.outerIndexPtr();
    for(int i = 0; i < n+1; ++i){
        ap_[i] = static_cast<int>(ref_ap[i]);
    }
    int ret = solver_->Analyze(false,  // complex or real
                               n,
                               ap_.get(),
                               ai_.get(),
                               J.valuePtr(),
                               nb_thread_);
    if (ret < 0){
        // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
        return ErrorType::SolverAnalyze;
    }
    return ErrorType::NoError;
}

ErrorType CKTSOLinearSolver::factorize(const EigenRefConstRealSpMat & J){
    int ret = solver_->Factorize(J.valuePtr(),
                                 true // @fast: whether to use fast factorization
                                 );
    if (ret < 0){
        // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
        return ErrorType::SolverFactor;
    }
    return ErrorType::NoError;
}

ErrorType CKTSOLinearSolver::refactorize(const EigenRefConstRealSpMat & J){
    int ret = solver_->Refactorize(J.valuePtr());
    if(ret < 0){
        // std::cout << "CKTSOLinearSolver::refactor solver_->Refactorize error: " << ret << std::endl;
        // https://github.com/chenxm1986/cktso/blob/master/include/cktso.h for error info
        return ErrorType::SolverReFactor;
    }
    return ErrorType::NoError;
}

ErrorType CKTSOLinearSolver::solve(Eigen::Ref<RealVect> b){
    RealVect x(b.size());
    int ret = solver_->Solve(&b(0), &x(0), false, 1);
    if(ret < 0){
        // std::cout << "CKTSOLinearSolver::solve solver_->Solve error: " << ret << std::endl;
        return ErrorType::SolverSolve;
    }
    b = x;
    return ErrorType::NoError;
}

} // namespace ls2g
