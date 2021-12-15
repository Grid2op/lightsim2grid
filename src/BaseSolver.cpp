// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "BaseSolver.h"

void BaseSolver::reset(){
    // reset timers
    reset_timer();

    //reset the attribute
    n_ = -1;
    Vm_ = RealVect();  // voltage magnitude
    Va_= RealVect();  // voltage angle
    V_= RealVect();  // voltage angle
    nr_iter_ = 0;  // number of iteration performs by the algorithm
    err_ = ErrorType::NotInitError; //error message:
}

RealVect BaseSolver::_evaluate_Fx(const Eigen::SparseMatrix<cplx_type> &  Ybus,
                                  const CplxVect & V,
                                  const CplxVect & Sbus,
                                  const Eigen::VectorXi & pv,
                                  const Eigen::VectorXi & pq)
{
    auto timer = CustTimer();
    auto npv = pv.size();
    auto npq = pq.size();

    // compute the mismatch
    CplxVect tmp = Ybus * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - Sbus.array();
    auto real_ = mis.real();
    auto imag_ = mis.imag();

    // build and fill the result
    RealVect res(npv + 2*npq);
    res.segment(0,npv) = real_(pv);
    res.segment(npv,npq) = real_(pq);
    res.segment(npv+npq, npq) = imag_(pq);
    timer_Fx_ += timer.duration();
    return res;
}

RealVect BaseSolver::_evaluate_Fx(const Eigen::SparseMatrix<cplx_type> &  Ybus,
                                  const CplxVect & V,
                                  const CplxVect & Sbus,
                                  Eigen::Index slack_id,  // id of the ref slack bus
                                  real_type slack_absorbed,
                                  const RealVect & slack_weights,
                                  const Eigen::VectorXi & pv,
                                  const Eigen::VectorXi & pq)
{
    /**
    Remember, when this function is used:

    J has the shape:
    
    | s | slack_bus |               | (pvpq+1,1) |   (1, pvpq)  |  (1, pq)   |
    | l |  -------  |               |            | ------------------------- |
    | a | J11 | J12 | = dimensions: |            | (pvpq, pvpq) | (pvpq, pq) |
    | c | --------- |               |   ------   | ------------------------- |
    | k | J21 | J22 |               |  (pq, 1)   |  (pq, pvpq)  | (pq, pq)   |

    python implementation:
    `J11` = dS_dVa[array([pvpq]).T, pvpq].real
    `J12` = dS_dVm[array([pvpq]).T, pq].real
    `J21` = dS_dVa[array([pq]).T, pvpq].imag
    `J22` = dS_dVm[array([pq]).T, pq].imag

    `slack_bus` = is the representation of the equation for the reference slack bus dS_dVa[slack_bus_id, pvpq].real
    and dS_dVm[slack_bus_id, pq].real

    `slack` is the representation of the equation connecting together the slack buses (represented by slack_weights)
    the remaining pq components are all 0.

    So i need to make sure that the returned F has the same properties:

    - first element is the slack bus
    - then it's Va of pv buses
    - then it's Va of pq buses
    - then it's Vm of pq buses

    **/

    // TODO factorize with above (ugly copy paste)
    auto timer = CustTimer();
    auto npv = pv.size();
    auto npq = pq.size();

    // compute the mismatch
    CplxVect tmp = Ybus * V;  // this is a vector
    tmp = tmp.array().conjugate();  // i take the conjugate
    auto mis = V.array() * tmp.array() - Sbus.array() + slack_absorbed * slack_weights.array();
    RealVect real_ = mis.real();
    RealVect imag_ = mis.imag();

    // build and fill the result
    RealVect res(npv + 2 * npq + 1); // slack adds one component hence the '+1' also bellow)
    res.segment(1, npv) = real_(pv);
    res.segment(npv + 1, npq) = real_(pq);
    res.segment(npv + npq + 1, npq) = imag_(pq);
    res(0) = real_(slack_id);  // slack bus is first variable
    timer_Fx_ += timer.duration();
    return res;
    
}

bool BaseSolver::_check_for_convergence(const RealVect & F,
                                        real_type tol)
{
    auto timer = CustTimer();
    const auto norm_inf = F.lpNorm<Eigen::Infinity>();
    bool res =  norm_inf  < tol;
    timer_check_ += timer.duration();
    return res;
}

int BaseSolver::extract_slack_bus_id(const Eigen::VectorXi & pv,
                                     const Eigen::VectorXi & pq,
                                     unsigned int nb_bus)
{
    // pv: list of index of pv nodes
    // pq: list of index of pq nodes
    // nb_bus: total number of bus in the grid
    // returns: res: the id of the unique slack bus (throw an error if no slack bus is found)
    // /!\ does not support multiple slack bus!!!

    int res=-1;
    // run through both pv and pq nodes and declare they are not slack bus
    std::vector<bool> tmp(nb_bus, true);
    for(unsigned int k=0; k < pv.size(); ++k)
    {
        tmp[pv[k]] = false;
    }
    for(unsigned int k=0; k < pq.size(); ++k)
    {
        tmp[pq[k]] = false;
    }
    // run through all buses
    for(unsigned int k=0; k < nb_bus; ++k)
    {
        if(tmp[k])
        {
            res = k;
            break;
        }
    }
    if(res == -1){
        throw std::runtime_error("BaseSolver::extract_slack_bus_id: No slack bus is found in your grid");
    }
    return res;
}
