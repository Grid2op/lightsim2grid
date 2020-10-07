// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseSolver.h"

const cdouble BaseSolver::my_i = {0., 1.};

void BaseSolver::reset(){
    // reset timers
    reset_timer();

    //reset the attribute
    n_ = -1;
    Vm_ = Eigen::VectorXd();  // voltage magnitude
    Va_= Eigen::VectorXd();  // voltage angle
    V_= Eigen::VectorXcd();  // voltage angle
    nr_iter_ = 0;  // number of iteration performs by the algorithm
    err_ = -1; //error message:
}

Eigen::VectorXd BaseSolver::_evaluate_Fx(const Eigen::SparseMatrix<cdouble> &  Ybus,
                                         const Eigen::VectorXcd & V,
                                         const Eigen::VectorXcd & Sbus,
                                         const Eigen::VectorXi & pv,
                                         const Eigen::VectorXi & pq)
{
    auto timer = CustTimer();
    auto npv = pv.size();
    auto npq = pq.size();

    // compute the mismatch
    Eigen::VectorXcd tmp = Ybus * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - Sbus.array();
    auto real_ = mis.real();
    auto imag_ = mis.imag();

    // build and fill the result
    Eigen::VectorXd res(npv + 2*npq);
    res.segment(0,npv) = real_(pv);
    res.segment(npv,npq) = real_(pq);
    res.segment(npv+npq, npq) = imag_(pq);
    timer_Fx_ += timer.duration();
    return res;
}

bool BaseSolver::_check_for_convergence(const Eigen::VectorXd & F,
                                        double tol)
{
    auto timer = CustTimer();
    bool res =  F.lpNorm<Eigen::Infinity>()  < tol;
    timer_check_ += timer.duration();
    return res;
}