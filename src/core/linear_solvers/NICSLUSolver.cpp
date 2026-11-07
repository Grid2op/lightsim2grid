// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NICSLUSolver.hpp"

#include <iostream>

namespace ls2g {

ErrorType NICSLULinearSolver::reset(){
    // free everything
    solver_.Free();
    ai_.reset();  // created in NICSLUSolver::analyze
    ap_.reset();  // created in NICSLUSolver::analyze

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

ErrorType NICSLULinearSolver::analyze(const EigenRefConstRealSpMat & J){
    // NICSLU Analyze may use values for MC64 value-based scaling/ordering
    const auto n = J.cols();
    const unsigned int nnz = J.nonZeros();

    ai_ = std::make_unique<unsigned int[]>(nnz);  // freed in destructor and NICSLUSolver::reset
    const int * ref_ai = J.innerIndexPtr();
    for(unsigned int i = 0; i < nnz; ++i){
        ai_[i] = static_cast<unsigned int>(ref_ai[i]);
    }

    ap_ = std::make_unique<unsigned int[]>(n+1);  // freed in destructor and NICSLUSolver::reset
    const int * ref_ap = J.outerIndexPtr();
    for(int i = 0; i < n+1; ++i){
        ap_[i] = static_cast<unsigned int>(ref_ap[i]);
    }
    int ret = solver_.Analyze(n,
                              J.valuePtr(),
                              ai_.get(),
                              ap_.get(),
                              MATRIX_COLUMN_REAL);  // MATRIX_COLUMN_REAL, MATRIX_ROW_REAL
    if (ret < 0){
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return ErrorType::SolverAnalyze;
    }
    return ErrorType::NoError;
}

ErrorType NICSLULinearSolver::factorize(const EigenRefConstRealSpMat & J){
    // solver.FactorizeMatrix(ax, 0); //use all created threads
    int ret = factorize_matrix(J);
    if (ret < 0){
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return ErrorType::SolverFactor;
    }
    return ErrorType::NoError;
}

ErrorType NICSLULinearSolver::refactorize(const EigenRefConstRealSpMat & J){
    int ret = factorize_matrix(J);  // TODO maybe 0 instead of nb_thread_ here, see https://github.com/chenxm1986/nicslu/blob/master/nicslu202110/demo/demo2.cpp
    if(ret < 0){
        // std::cout << "NICSLULinearSolver::refactor solver_.FactorizeMatrix error: " << ret << std::endl;
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return ErrorType::SolverReFactor;
    }
    return ErrorType::NoError;
}

int NICSLULinearSolver::factorize_matrix(const EigenRefConstRealSpMat & J){
    // stat[14] / stat[15]: the factorizations / refactorizations NICSLU performed. Read
    // around the call, so that whether Analyze restarts them does not matter. No
    // statistics without a successful Initialize (a missing licence): nothing ran either.
    last_nb_factorize_ = 0;
    last_nb_refactorize_ = 0;
    if(solver_.info == nullptr) return solver_.FactorizeMatrix(J.valuePtr(), nb_thread_);
    const double nb_fact = solver_.GetInformation(14);
    const double nb_refact = solver_.GetInformation(15);
    int ret = solver_.FactorizeMatrix(J.valuePtr(), nb_thread_);
    last_nb_factorize_ = static_cast<std::size_t>(solver_.GetInformation(14) - nb_fact);
    last_nb_refactorize_ = static_cast<std::size_t>(solver_.GetInformation(15) - nb_refact);
    return ret;
}

ErrorType NICSLULinearSolver::solve(Eigen::Ref<RealVect> b){
    RealVect x(b.size());
    int ret = solver_.Solve(&b(0), &x(0));
    if(ret < 0){
        // std::cout << "NICSLULinearSolver::solve solver_.Solve error: " << ret << std::endl;
        return ErrorType::SolverSolve;
    }
    b = x;
    return ErrorType::NoError;
}

ErrorType NICSLULinearSolver::solve_transpose(Eigen::Ref<RealVect> /*b*/){
    // Deliberately not solving anything: see CAN_SOLVE_TRANSPOSE above. Returning an
    // error rather than solving J x = b keeps a caller that ignored the flag from
    // getting a plausible-looking answer to the wrong system.
    return ErrorType::NotImplemented;
}

} // namespace ls2g
