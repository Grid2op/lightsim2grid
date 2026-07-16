// Copyright (c) 2020-2026, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "KLUSolver.hpp"

#include <iostream>

namespace ls2g {

const bool KLULinearSolver::CAN_SOLVE_MAT = false;

ErrorType KLULinearSolver::reset(){
    // release both handles (their deleters use common_) before resetting common_
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::analyze(const Eigen::SparseMatrix<real_type>& J){
    // J is const here, but `klu_analyze` signature expects non-const arrays, so I const_cast
    const auto n = J.cols();
    // free any previous factorization (using the current common_) before resetting it,
    // so re-analyzing without a prior reset() no longer leaks the old handles
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    symbolic_.reset(klu_analyze(n,
                                const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                                const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                                &common_));
    if(common_.status != KLU_OK) return ErrorType::SolverAnalyze;
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::factorize(const Eigen::SparseMatrix<real_type>& J){
    // J is const here, but `klu_factor` signature expects non-const arrays, so I const_cast
    // reset() frees any previous numeric factorization before storing the new one
    numeric_.reset(klu_factor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                              const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                              const_cast<real_type*>(J.valuePtr()),
                              symbolic_.get(), &common_));
    if(common_.status != KLU_OK) return ErrorType::SolverFactor;
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::refactorize(const Eigen::SparseMatrix<real_type>& J){
    // J is const here, but `klu_refactor` signature expects arrays and not const arrays
    // so I const_cast
    int ok = klu_refactor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                          const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                          const_cast<real_type*>(J.valuePtr()),
                          symbolic_.get(), numeric_.get(), &common_);
    if(ok != 1){
        // std::cout << "\t KLU: refactor error" << std::endl;
        return ErrorType::SolverReFactor;
    }
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::solve(Eigen::Ref<RealVect> b){
    const int n = static_cast<int>(b.size());
    int ok = klu_solve(symbolic_.get(), numeric_.get(), n, 1, &b(0), &common_);
    if(ok != 1){
        // std::cout << "\t KLU: klu_solve error" << std::endl;
        return ErrorType::SolverSolve;
    }
    return ErrorType::NoError;
}

} // namespace ls2g
