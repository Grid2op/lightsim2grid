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
    if(symbolic_ != nullptr) klu_free_symbolic(&symbolic_, &common_);
    if(numeric_ != nullptr) klu_free_numeric(&numeric_, &common_);
    common_ = klu_common();
    symbolic_ = nullptr;
    numeric_ = nullptr;
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::initialize(const Eigen::SparseMatrix<real_type>&  J){
    // default Eigen representation: column major, which is good for klu !
    // J is const here, but `klu_analyze` signature expects arrays and not const arrays
    // so I const_cast
    const auto n = J.cols();
    common_ = klu_common();
    ErrorType res = ErrorType::NoError; 
    symbolic_ = klu_analyze(n,
                            const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                            const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                            &common_);
    if(common_.status != KLU_OK){
        res = ErrorType::SolverAnalyze; 
    }else{
        numeric_ = klu_factor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                              const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                              const_cast<real_type*>(J.valuePtr()),
                              symbolic_, &common_);
        if(common_.status != KLU_OK) res = ErrorType::SolverFactor; 
    }
    return res;
}

ErrorType KLULinearSolver::refactor(const Eigen::SparseMatrix<real_type>& J){
    // J is const here, but `klu_refactor` signature expects arrays and not const arrays
    // so I const_cast
    int ok = klu_refactor(const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.outerIndexPtr()),
                          const_cast<Eigen::SparseMatrix<real_type>::StorageIndex *>(J.innerIndexPtr()),
                          const_cast<real_type*>(J.valuePtr()),
                          symbolic_, numeric_, &common_);
    if(ok != 1){
        // std::cout << "\t KLU: refactor error" << std::endl;
        return ErrorType::SolverReFactor;
    }
    return ErrorType::NoError;
}

ErrorType KLULinearSolver::solve(RealVect & b){
    const int n = static_cast<int>(b.size());
    int ok = klu_solve(symbolic_, numeric_, n, 1, &b(0), &common_);
    if(ok != 1){
        // std::cout << "\t KLU: klu_solve error" << std::endl;
        return ErrorType::SolverSolve;
    }
    return ErrorType::NoError;
}

} // namespace ls2g
