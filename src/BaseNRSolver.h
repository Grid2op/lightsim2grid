// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASENRSOLVER_H
#define BASENRSOLVER_H

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstdint> // for int32
#include <chrono>
#include <complex>      // std::complex, std::conj
#include <cmath>  // for PI

// eigen is necessary to easily pass data from numpy to c++ without any copy.
// and to optimize the matrix operations
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

#include "CustTimer.h"
#include "Utils.h"

class BaseNRSolver
{
    public:
        BaseNRSolver():n_(-1),need_factorize_(true),err_(-1),timer_Fx_(0.){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_initialize_ = 0.;
            timer_check_ = 0.;
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_total_nr_ = 0.;
        }

        ~BaseNRSolver(){}

        Eigen::SparseMatrix<double> get_J(){
            return J_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Va(){
            return Va_;
        }
        Eigen::Ref<Eigen::VectorXd> get_Vm(){
            return Vm_;
        }
        Eigen::Ref<Eigen::VectorXcd> get_V(){
            return V_;
        }
        int get_error(){
            return err_;
        }
        int get_nb_iter(){
            return nr_iter_;
        }
        std::tuple<double, double, double, double, double, double, double> get_timers()
        {
            auto res = std::tuple<double, double, double, double, double, double, double>(
              timer_Fx_, timer_solve_, timer_initialize_, timer_check_, timer_dSbus_, timer_fillJ_, timer_total_nr_);
            return res;
        }

        bool do_newton(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       ) ;

        virtual
        void reset();

        bool converged(){
            return err_ == 0;
        }

    protected:
        void reset_timer(){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_initialize_ = 0.;
            timer_check_ = 0.;
            timer_dSbus_ = 0.;
            timer_fillJ_ = 0.;
            timer_total_nr_ = 0.;
        }
        virtual
        void initialize()=0;

        virtual
        void solve(Eigen::VectorXd & b, bool has_just_been_inialized)=0;

        void _dSbus_dV(const Eigen::Ref<const Eigen::SparseMatrix<cdouble> > & Ybus,
                       const Eigen::Ref<const Eigen::VectorXcd > & V);

        void _get_values_J(int & nb_obj_this_col,
                           std::vector<int> & inner_index,
                           std::vector<double> & values,
                           const Eigen::SparseMatrix<double> & mat,  // ex. dS_dVa_r
                           const std::vector<int> & index_row_inv, // ex. pvpq_inv
                           const Eigen::VectorXi & index_col, // ex. pvpq
                           int col_id,
                           int row_lag  // 0 for J11 for example, n_pvpq for J12
                           );

        void fill_jacobian_matrix(const Eigen::SparseMatrix<cdouble> & Ybus,
                                  const Eigen::VectorXcd & V,
                                  const Eigen::VectorXi & pq,
                                  const Eigen::VectorXi & pvpq,
                                  const std::vector<int> & pq_inv,
                                  const std::vector<int> & pvpq_inv
                                  );

        Eigen::VectorXd _evaluate_Fx(const Eigen::SparseMatrix<cdouble> &  Ybus,
                                     const Eigen::VectorXcd & V,
                                     const Eigen::VectorXcd & Sbus,
                                     const Eigen::VectorXi & pv,
                                     const Eigen::VectorXi & pq);

        bool _check_for_convergence(const Eigen::VectorXd & F,
                                    double tol)
        {
            auto timer = CustTimer();
            bool res =  F.lpNorm<Eigen::Infinity>()  < tol;
            timer_check_ += timer.duration();
            return res;
        }

    protected:
        // solver initialization
        int n_;

        // solution of the problem
        Eigen::VectorXd Vm_;  // voltage magnitude
        Eigen::VectorXd Va_;  // voltage angle
        Eigen::VectorXcd V_;  // voltage angle
        Eigen::SparseMatrix<double> J_;  // the jacobian matrix
        Eigen::SparseMatrix<cdouble> dS_dVm_;
        Eigen::SparseMatrix<cdouble> dS_dVa_;
        bool need_factorize_;
        int nr_iter_;  // number of iteration performs by the Newton Raphson algorithm
        int err_; //error message:
        // -1 : the solver has not been initialized (call initialize in this case)
        // 0 everything ok
        // 1: i can't factorize the matrix (klu_factor)
        // 2: i can't refactorize the matrix (klu_refactor)
        // 3: i can't solve the system (klu_solve)
        // 4: end of possible iterations (divergence because nr_iter_ >= max_iter

        // timers
         double timer_Fx_;
         double timer_solve_;
         double timer_initialize_;
         double timer_check_;
         double timer_dSbus_;
         double timer_fillJ_;
         double timer_total_nr_;
         static const cdouble my_i;

    private:
        // no copy allowed
        BaseNRSolver( const BaseNRSolver & ) ;
        BaseNRSolver & operator=( const BaseNRSolver & ) ;

    private:
        Eigen::SparseMatrix<cdouble>
            _make_diagonal_matrix(const Eigen::Ref<const Eigen::VectorXcd > & diag_val){
            // TODO their might be a more efficient way to do that
            auto n = diag_val.size();
            Eigen::SparseMatrix<cdouble> res(n,n);
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

#endif // BASENRSOLVER_H
