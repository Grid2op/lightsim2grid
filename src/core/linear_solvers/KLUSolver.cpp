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

ErrorType KLULinearSolver::reset(){
    // release both handles (their deleters use common_) before resetting common_
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    klu_defaults(&common_);
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::analyze(const EigenRefConstRealSpMat & J){
    // J is const here, but `klu_analyze` signature expects non-const arrays, so I const_cast
    const auto n = J.cols();
    // free any previous factorization (using the current common_) before resetting it,
    // so re-analyzing without a prior reset() no longer leaks the old handles
    numeric_.reset();
    symbolic_.reset();
    common_ = klu_common();
    klu_defaults(&common_);
    symbolic_.reset(klu_analyze(n,
                                const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                                const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                                &common_));
    if(common_.status != KLU_OK) return ErrorType::SolverAnalyze;
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::factorize(const EigenRefConstRealSpMat & J){
    // J is const here, but `klu_factor` signature expects non-const arrays, so I const_cast
    // reset() frees any previous numeric factorization before storing the new one
    numeric_.reset(klu_factor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                              const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                              const_cast<real_type*>(J.valuePtr()),
                              symbolic_.get(), &common_));
    if(common_.status != KLU_OK) return ErrorType::SolverFactor;
    if(guard_pivot_growth_){
        // the reference a later refactorization is judged against
        rcond_ref_ = 0.;
        if(klu_rcond(symbolic_.get(), numeric_.get(), &common_) == 1) rcond_ref_ = common_.rcond;
    }
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::refactorize(const EigenRefConstRealSpMat & J){
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
    if(guard_pivot_growth_ && rcond_ref_ > 0.){
        // klu_rcond is min|diag(U)| / max|diag(U)|, O(n). Healthy refactorizations stay
        // within a modest factor of the reference; one that falls by this many orders of
        // magnitude has pivots that no longer suit the matrix
        constexpr double MAX_RCOND_DROP = 1e-8;
        if(klu_rcond(symbolic_.get(), numeric_.get(), &common_) != 1 ||
           !(common_.rcond >= MAX_RCOND_DROP * rcond_ref_)) return ErrorType::SolverReFactor;
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

ErrorType KLULinearSolver::solve_transpose(Eigen::Ref<RealVect> b){
    // klu_tsolve walks the very same L / U factors as klu_solve, backwards: J^T is
    // never formed, and neither the symbolic analysis nor the numeric factorization
    // is repeated. That is what makes an adjoint (a vector-Jacobian product, as
    // reverse-mode differentiation of a powerflow needs) cost one triangular solve.
    const int n = static_cast<int>(b.size());
    int ok = klu_tsolve(symbolic_.get(), numeric_.get(), n, 1, &b(0), &common_);
    if(ok != 1){
        // std::cout << "\t KLU: klu_tsolve error" << std::endl;
        return ErrorType::SolverSolve;
    }
    return ErrorType::NoError;
}

} // namespace ls2g
