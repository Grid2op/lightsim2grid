// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "KLUSolver.h"

void KLUSolver::reset(){
    BaseNRSolver::reset();
    klu_free_symbolic(&symbolic_, &common_);
    klu_free_numeric(&numeric_, &common_);
    n_ = -1;
    common_ = klu_common();

    symbolic_ = nullptr;
    numeric_ = nullptr;
}

void KLUSolver::initialize(){
    // default Eigen representation: column major, which is good for klu !
    // J is const here, even if it's not said in klu_analyze
    auto timer = CustTimer();
    n_ = static_cast<int>(J_.cols()); // should be equal to J_.nrows()
    err_ = 0; // reset error message
    common_ = klu_common();
    symbolic_ = klu_analyze(n_, J_.outerIndexPtr(), J_.innerIndexPtr(), &common_);
    numeric_ = klu_factor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, &common_);
    if (common_.status != KLU_OK) {
        err_ = 1;
    }
    need_factorize_ = false;
    timer_initialize_ += timer.duration();
}

void KLUSolver::solve(RealVect & b, bool has_just_been_initialized){
    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call klu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    auto timer = CustTimer();
    int ok;
    bool stop = false;
    if(!has_just_been_initialized){
        // if the call to "klu_factor" has been made this iteration, there is no need
        // to re factor again the matrix
        // i'm in the case where it has not
        ok = klu_refactor(J_.outerIndexPtr(), J_.innerIndexPtr(), J_.valuePtr(), symbolic_, numeric_, &common_);
        if (ok != 1) {
            err_ = 2;
            stop = true;
        }
    }
    if(!stop){
        ok = klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
        if (ok != 1) {
            err_ = 3;
        }
    }
    // std::cout << "KLUSolver::solve" << std::endl;
    timer_solve_ += timer.duration();
}
