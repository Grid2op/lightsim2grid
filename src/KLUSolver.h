// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifdef KLU_SOLVER_AVAILABLE
#ifndef KLSOLVER_H
#define KLSOLVER_H

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <cmath>  // for PI

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Utils.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

// import klu package
extern "C" {
    #include "cs.h"
    #include "klu.h"
}

#include "CustTimer.h"
#include "BaseNRSolver.h"
/**
class to handle the solver using newton-raphson method, using KLU algorithm and sparse matrices.

As long as the admittance matrix of the sytem does not change, you can reuse the same solver.
Reusing the same solver is possible, but "reset" method must be called.

Otherwise, unexpected behaviour might follow, including "segfault".

**/
class KLUSolver: public BaseNRSolver
{
    public:
        KLUSolver():BaseNRSolver(),symbolic_(),numeric_(),common_(){}

        ~KLUSolver()
         {
             klu_free_symbolic(&symbolic_, &common_);
             klu_free_numeric(&numeric_, &common_);
         }

        virtual void reset();

    protected:
        virtual
        void initialize();

        virtual
        void solve(RealVect & b, bool has_just_been_inialized);

    private:
        // solver initialization
        klu_symbolic* symbolic_;
        klu_numeric* numeric_;
        klu_common common_;

        // no copy allowed
        KLUSolver( const KLUSolver & ) ;
        KLUSolver & operator=( const KLUSolver & ) ;

    private:
        // debug func i don't want to remove yet
        void analyze_old(int n,
                      Eigen::Ref<Eigen::VectorXi> Ap,
                      Eigen::Ref<Eigen::VectorXi> Ai){
            n_ = n;
            symbolic_ = klu_analyze(n, &Ap(0), &Ai(0), &common_);
        }
        void solve_old(Eigen::Ref<Eigen::VectorXi> Ap,
                    Eigen::Ref<Eigen::VectorXi> Ai,
                    Eigen::Ref<RealVect> Ax,
                    Eigen::Ref<RealVect> b){
            numeric_ = klu_factor(&Ap(0), &Ai(0), &Ax(0), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &b(0), &common_);
        }

        // TODO re add the references here for the last stuff
        std::tuple<RealVect, CplxVect> one_iter_test(Eigen::SparseMatrix<real_type> J,
                              Eigen::Ref<RealVect> F,
                              Eigen::VectorXi pv,
                              Eigen::VectorXi pq,
                              CplxVect V,
                              Eigen::SparseMatrix<cplx_type>  Ybus,
                              CplxVect Sbus
                              ){
            //TODO do not use, for DEBUG only!!!
            // get the sizes for convenience
            auto npv = pv.size();
            auto npq = pq.size();

            // supposes that "klu_analyze" has been already called ! klu_solver.analyze(J) should have been called
            numeric_ = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic_, &common_);
            klu_solve(symbolic_, numeric_, n_, 1, &F(0), &common_);
            auto dx = -1.0*F;

            // update voltage (this should be done consistently with "klu_solver._evaluate_Fx")
            Vm_ = V.array().abs();  // update Vm and Va again in case
            Va_ = V.array().arg();  // we wrapped around with a negative Vm

            if (npv > 0) Va_(pv) += dx.segment(0,npv);
            if (npq > 0){
                Va_(pq) += dx.segment(npv,npq);
                Vm_(pq) += dx.segment(npv+npq, npq);
            }

            // TODO change here for not having to cast all the time ...
            V = Vm_.array() * (Va_.array().cos().cast<cplx_type>() + my_i * Va_.array().sin().cast<cplx_type>() );

            F = _evaluate_Fx(Ybus, V, Sbus, pv, pq);
            return std::tuple<RealVect, CplxVect>(F, V);
        }

        Eigen::SparseMatrix<real_type>
             create_jacobian_matrix_test(const Eigen::SparseMatrix<cplx_type> & Ybus,
                                         const CplxVect & V,
                                         const Eigen::VectorXi & pq,
                                         const Eigen::VectorXi & pvpq
                                         ){

            // DO NOT USE, FOR DEBUG ONLY!
            const auto & n_pvpq = pvpq.size();
            const auto & n_pq = pvpq.size();
            std::vector<int> pvpq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pvpq; ++inv_id) pvpq_inv[pvpq(inv_id)] = inv_id;
            std::vector<int> pq_inv(V.size(), -1);
            for(int inv_id=0; inv_id < n_pq; ++inv_id) pq_inv[pq(inv_id)] = inv_id;
            fill_jacobian_matrix(Ybus, V, pq, pvpq, pq_inv, pvpq_inv);
            return J_;
        }

        bool initialize_test(Eigen::SparseMatrix<real_type > & J){
            // default Eigen representation: column major, which is good for klu !
            // J is const here, even if it's not said in klu_analyze
            const int n = static_cast<int>(J.cols()); // should be equal to J_.nrows()
            err_ = 0; // reset error message
            klu_common common = klu_common();
            bool res = true;
            klu_symbolic* symbolic = klu_analyze(n, J.outerIndexPtr(), J.innerIndexPtr(), &common);
//            auto numeric = klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic, &common);
            klu_factor(J.outerIndexPtr(), J.innerIndexPtr(), J.valuePtr(), symbolic, &common);
            if (common_.status != KLU_OK) {
                err_ = 1;
                res = false;
            }
            return res;
        }

        std::tuple<Eigen::SparseMatrix<cplx_type> , Eigen::SparseMatrix<cplx_type> >
                    _get_ds_test(Eigen::SparseMatrix<cplx_type> & Ybus,
                                CplxVect & V){
            BaseNRSolver::_dSbus_dV(Ybus, V);
            auto res = std::tuple<Eigen::SparseMatrix<cplx_type> , Eigen::SparseMatrix<cplx_type> >(dS_dVm_, dS_dVa_);
            return res;
        }

        void _dSbus_dV(Eigen::SparseMatrix<cplx_type> & dS_dVm,
                           Eigen::SparseMatrix<cplx_type> & dS_dVa,
                           const Eigen::Ref<const Eigen::SparseMatrix<cplx_type> > & Ybus,
                           const Eigen::Ref<const CplxVect > & V)
        {
            // "slow" implementation close to pypower, but with sparse matrix
            // TODO check i cannot optimize that with numba code in pandapower instead
            auto timer = CustTimer();
            CplxVect Ibus = Ybus * V;
            Eigen::SparseMatrix<cplx_type> diagV = _make_diagonal_matrix(V);

            CplxVect Ibus_conj = Ibus.conjugate();
            Eigen::SparseMatrix<cplx_type> diagIbus_conj = _make_diagonal_matrix(Ibus_conj);

            CplxVect Vnorm = V.array() / V.array().abs();
            Eigen::SparseMatrix<cplx_type> diagVnorm = _make_diagonal_matrix(Vnorm);

            Eigen::SparseMatrix<cplx_type> tmp = Ybus * diagVnorm;
            tmp = tmp.conjugate();
            dS_dVm = diagV * tmp + diagIbus_conj * diagVnorm;

            // TODO this is the same code for this "tmp" and the previous one, except for "diagV" and "diagVnorm"
            tmp = Ybus * diagV;
            tmp = tmp.conjugate();
            auto tmp2 = diagIbus_conj - tmp;  // conj(diagIbus - Ybus * diagV)
            dS_dVa = my_i * diagV * tmp2;

            // python implementation
            // dS_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm
            // dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
            timer_dSbus_ += timer.duration();
        }

        Eigen::SparseMatrix<cplx_type>
            _make_diagonal_matrix(const Eigen::Ref<const CplxVect > & diag_val){
            // TODO their might be a more efficient way to do that
            auto n = diag_val.size();
            Eigen::SparseMatrix<cplx_type> res(n,n);
            // first method, without a loop of mine
            res.setIdentity();  // segfault if attempt to use this function without this
            res.diagonal() = diag_val;

            // second method, with an "optimized" loop
            // res.reserve(Eigen::VectorXi::Constant(n,1)); // i reserve one number per columns (speed optim)
            // for(unsigned int i = 0; i < n; ++i){
            //    res.insert(i,i) = diag_val(i);
            //}
            return res;
        }
};

#endif // KLSOLVER_H
#endif  // KLU_SOLVER_AVAILABLE
