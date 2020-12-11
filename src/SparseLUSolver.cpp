// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "SparseLUSolver.h"

void SparseLUSolver::initialize(){
    // default Eigen representation: column major, which is good for klu !
    // J is const here
    auto timer = CustTimer();
    n_ = J_.cols(); // should be equal to J_.nrows()
    err_ = 0; // reset error message
    J_.makeCompressed();
    solver_.analyzePattern(J_);  //NEW
    solver_.factorize(J_);  // NEW
    if (solver_.info() != Eigen::Success) {
        err_ = 1;
    }
    need_factorize_ = false;
    timer_solve_ += timer.duration();
}

void SparseLUSolver::solve(RealVect & b, bool has_just_been_inialized){
    // NEW

    // solves (for x) the linear system J.x = b
    // supposes that the solver has been initialized (call sparselu_solver.analyze() before calling that)
    // J is const even if it does not compile if said const
    auto timer = CustTimer();
    bool stop = false;
    if(!has_just_been_inialized){
        // if the call to "klu_factor" has been made this iteration, there is no need
        // to re factor again the matrix
        // i'm in the case where it has not
        solver_.factorize(J_);  // NEW
        if (solver_.info() != Eigen::Success) {
            err_ = 2;
            stop = true;
        }
    }
    if(!stop){
        RealVect Va = solver_.solve(b);  //NEW
        if (solver_.info() != Eigen::Success) {
            err_ = 3;
        }
        b = Va;  // NEW
    }
    timer_solve_ += timer.duration();
}
