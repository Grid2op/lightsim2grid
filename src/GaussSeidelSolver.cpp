// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#include "GaussSeidelSolver.h"

bool GaussSeidelSolver::compute_pf(const Eigen::SparseMatrix<cdouble> & Ybus,
                                   Eigen::VectorXcd & V,
                                   const Eigen::VectorXcd & Sbus,
                                   const Eigen::VectorXi & pv,
                                   const Eigen::VectorXi & pq,
                                   int max_iter,
                                   double tol
                                   )
{
    /**
    pv: id of the pv buses
    pq: id of the pq buses

    **/
    // TODO check what can be checked: no voltage at 0, Ybus is square, Sbus same size than V and
    // TODO Ybus (nrow or ncol), pv and pq have value that are between 0 and nrow etc.
    reset_timer();
    if(err_ > 0) return false; // i don't do anything if there were a problem at the initialization
    auto timer = CustTimer();

    V_ = V;
    Vm_ = V_.array().abs();  // update Vm and Va again in case
    Va_ = V_.array().arg();  // we wrapped around with a negative Vm

    // first check, if the problem is already solved, i stop there
    Eigen::VectorXd F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
    bool converged = _check_for_convergence(F, tol);
    nr_iter_ = 0; //current step
    bool res = true;  // have i converged or not
    Eigen::VectorXcd tmp_Sbus = Sbus;
    while ((!converged) & (nr_iter_ < max_iter)){
        nr_iter_++;

        // ###########################
        // the Gauss Seidel Algorithm
        // ###########################
        // https://www.sanfoundry.com/cpp-program-implement-gauss-seidel-method/
        // https://fr.mathworks.com/matlabcentral/fileexchange/14089-gauss-seidel-load-flow-analysis
        // https://github.com/rwl/PYPOWER/blob/master/pypower/gausspf.py

        auto timer2 = CustTimer();
        // one_iter_all_at_once(tmp_Sbus, Ybus, pv, pq);
        one_iter(tmp_Sbus, Ybus, pv, pq);
        timer_solve_ += timer2.duration();

        // #####################
        // stopping criteria
        // #####################
        F = _evaluate_Fx(Ybus, V_, tmp_Sbus, pv, pq);
        bool tmp = F.allFinite();
        if(!tmp) break; // divergence due to Nans
        converged = _check_for_convergence(F, tol);
    }
    if(!converged){
        err_ = 4;
        res = false;
    }
    Vm_ = V_.array().abs();  // update Vm and Va again in case
    Va_ = V_.array().arg();  // we wrapped around with a negative Vm
    timer_total_nr_ += timer.duration();
    return res;
}

void GaussSeidelSolver::one_iter(Eigen::VectorXcd & tmp_Sbus,
                                 const Eigen::SparseMatrix<cdouble> & Ybus,
                                 const Eigen::VectorXi & pv,
                                 const Eigen::VectorXi & pq)
{
    // do an update with the standard GS algorithm
    cdouble tmp;

    int n_pv = pv.size();
    int n_pq = pq.size();

    // update PQ buses
    for(int k_tmp=0; k_tmp<n_pq; ++k_tmp)
    {
        int k = pq.coeff(k_tmp);
        tmp = tmp_Sbus.coeff(k) / V_.coeff(k);
        tmp = std::conj(tmp);
        tmp -= static_cast<cdouble>(Ybus.row(k) * V_);
        tmp /= Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = Ybus.row(k) * V_;  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = tmp_Sbus.coeff(k) / V_.coeff(k);
        tmp = std::conj(tmp);
        tmp -= static_cast<cdouble>(Ybus.row(k) * V_);
        tmp /= Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}

void GaussSeidelSolver::one_iter_all_at_once(Eigen::VectorXcd & tmp_Sbus,
                                             const Eigen::SparseMatrix<cdouble> & Ybus,
                                             const Eigen::VectorXi & pv,
                                             const Eigen::VectorXi & pq)
{
    // do an update with all nodes being updated at the same time (different than the original GaussSeidel)
    cdouble tmp;

    int n_pv = pv.size();
    int n_pq = pq.size();

    // Eigen::VectorXcd tmp_YbusV;  // Ybus[k, :] * V
    // Eigen::VectorXcd tmp_conj_Sbus_V;  //  conj(Sbus[k] / V[k])
    Eigen::VectorXcd tmp_YbusV = Ybus * V_;
    Eigen::VectorXcd tmp_conj_Sbus_V = tmp_Sbus.array() / V_.array();
    tmp_conj_Sbus_V = tmp_conj_Sbus_V.array().conjugate();

    // update PQ buses
    for(int k_tmp=0; k_tmp<n_pq; ++k_tmp)
    {
        int k = pq.coeff(k_tmp);
        tmp = (tmp_conj_Sbus_V.coeff(k) -  tmp_YbusV.coeff(k)) / Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // update PV buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        // update Sbus
        tmp = tmp_YbusV.coeff(k);  // Ybus[k,:] * V
        tmp = std::conj(tmp);  // conj(Ybus[k,:] * V)
        tmp *= V_.coeff(k);  // (V[k] * conj(Ybus[k,:] * V))
        tmp = my_i * std::imag(tmp);
        tmp_Sbus.coeffRef(k) = std::real(tmp_Sbus.coeff(k)) + tmp;

        // update V
        tmp = (tmp_conj_Sbus_V(k) -  tmp_YbusV(k)) / Ybus.coeff(k,k);
        V_.coeffRef(k) += tmp;
    }

    // make sure the voltage magnitudes are not modified at pv buses
    for(int k_tmp=0; k_tmp<n_pv; ++k_tmp)
    {
        int k = pv.coeff(k_tmp);
        V_.coeffRef(k) *= Vm_.coeff(k) / std::abs(V_.coeff(k));
    }
}
