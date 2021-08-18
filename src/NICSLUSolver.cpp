// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "NICSLUSolver.h"

void NICSLUSolver::reset(){
    BaseNRSolver::reset();

    // free everything
    n_ = -1;
    solver_.Free();
    delete [] ai_;  // created in NICSLUSolver::initialize
    delete [] ap_;  // created in NICSLUSolver::initialize
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
        std::cout << "NICSLUSolver::reset() solver_.Initialize error: "  << ret << std::endl;
        err_ = 1;  // TODO clarify that
        return;
    }
    // solver_.SetConfiguration(0, 1.); //enable timer
}

void NICSLUSolver::initialize(){
    // default Eigen representation: column major, which is good for klu !
    auto timer = CustTimer();
    int ret;
    n_ = J_.cols(); // should be equal to J_.nrows()
    err_ = 0; // reset error message
    //    symbolic_ = klu_analyze(n_, J_.outerIndexPtr(), J_.innerIndexPtr(), &common_);
    //    numeric_ = klu_factor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, &common_);

    // solver.Analyze(n, ax, ai, ap, MATRIX_ROW_REAL);
    // convert int * (storage class for Eign) to unsigned int * (expected storage class for NICSLU...)

    // std::cout << "n_ = " << n_ << std::endl;  // DEBUG
    const unsigned int nnz = J_.nonZeros();  // DEBUG
    // std::cout << "nnz = " << nnz << std::endl;  // DEBUG

    ai_ = new unsigned int [nnz];  // deleted in destructor and NICSLUSolver::reset
    const int * ref_ai = J_.innerIndexPtr();
    for(unsigned int i = 0; i < nnz; ++i){
        ai_[i] = static_cast<unsigned int>(ref_ai[i]);
        // std::cout << "ai_["<<i<<"] = " << ai_[i] << std::endl;  // DEBUG
    }
    // std::cout << std::endl;  // DEBUG

    ap_ = new unsigned int [n_+1];  // deleted in destructor and NICSLUSolver::reset
    const int * ref_ap = J_.outerIndexPtr();
    for(int i = 0; i < n_+1; ++i){
        ap_[i] = static_cast<unsigned int>(ref_ap[i]);
        // std::cout << "ap_["<<i<<"] = " << ap_[i] << std::endl;  // DEBUG
    }
    // std::cout << std::endl;  // DEBUG

    // const double * ref_ax = J_.valuePtr();  // DEBUG
    // for(unsigned int i = 0; i < nnz; ++i){  // DEBUG
    //     std::cout << "ax["<<i<<"] = " << ref_ax[i] << std::endl;  // DEBUG
    // }  // DEBUG

    ret = solver_.Analyze(n_,
                          J_.valuePtr(),
                          ai_,
                          ap_,
                          MATRIX_COLUMN_REAL);  // MATRIX_COLUMN_REAL, MATRIX_ROW_REAL
    if (ret < 0){
        err_ = 1;
        std::cout << "NICSLUSolver::initialize() solver_.Analyze error: "  << ret << std::endl;
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return;
    }
    // solver.FactorizeMatrix(ax, 0); //use all created threads
    ret = solver_.FactorizeMatrix(J_.valuePtr(), nb_thread_);
    if (ret < 0){
        err_ = 1;
        std::cout << "NICSLUSolver::initialize() solver_.FactorizeMatrix error: "  << ret << std::endl;
        // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
        return;
    }
    need_factorize_ = false;
    timer_initialize_ += timer.duration();
}

void NICSLUSolver::solve(RealVect & b, bool has_just_been_inialized){
    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call klu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    auto timer = CustTimer();
    int ret;
    bool stop = false;
    RealVect x;
    if(!has_just_been_inialized){
        // if the call to "klu_factor" has been made this iteration, there is no need
        // to re factor again the matrix
        // i'm in the case where it has not
        //        ret = klu_refactor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, numeric_, &common_);

        // solver.FactorizeMatrix(ax, 0); //use all created threads
        ret  = solver_.FactorizeMatrix(J_.valuePtr(), nb_thread_);
        if (ret < 0) {
            std::cout << "NICSLUSolver::solve solver_.FactorizeMatrix error: " << ret << std::endl;
            // https://github.com/chenxm1986/nicslu/blob/master/nicslu202103/include/nicslu.h for error info
            err_ = 2;
            stop = true;
        }
    }
    if(!stop){
        // ok = klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
        // solver.Solve(b, x);
        x = RealVect(n_);
        ret = solver_.Solve(&b(0), &x(0));
        if (ret < 0) {
            std::cout << "NICSLUSolver::solve solver_.Solve error: " << ret << std::endl;
            err_ = 3;
        }
        b = x;
    }
    timer_solve_ += timer.duration();
}
