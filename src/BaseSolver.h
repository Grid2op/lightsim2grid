// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef BASESOLVER_H
#define BASESOLVER_H

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

// TODO make err_ more explicit: use an enum

/**
This class represents a solver to compute powerflow.

It can be derived for different usecase, for example for DC powerflow, AC powerflow using Newton Raphson method etc.
**/
class BaseSolver
{
    public:
        BaseSolver():n_(-1),err_(-1),timer_Fx_(0.),timer_solve_(0.),timer_check_(0.),timer_total_nr_(0.){};

        ~BaseSolver(){}

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

        virtual
        std::tuple<double, double, double, double> get_timers()
        {
            // TODO change the order of the timers here!
            auto res = std::tuple<double, double, double, double>(
              timer_Fx_, timer_solve_, timer_check_, timer_total_nr_);
            return res;
        }

        virtual
        bool compute_pf(const Eigen::SparseMatrix<cdouble> & Ybus,
                        Eigen::VectorXcd & V,
                        const Eigen::VectorXcd & Sbus,
                        const Eigen::VectorXi & pv,
                        const Eigen::VectorXi & pq,
                        int max_iter,
                        double tol
                        ) = 0 ;

        virtual
        void reset();

        bool converged(){
            return err_ == 0;
        }

    protected:
        void reset_timer(){
            timer_Fx_ = 0.;
            timer_solve_ = 0.;
            timer_check_ = 0.;
            timer_total_nr_ = 0.;
        }

        Eigen::VectorXd _evaluate_Fx(const Eigen::SparseMatrix<cdouble> &  Ybus,
                                     const Eigen::VectorXcd & V,
                                     const Eigen::VectorXcd & Sbus,
                                     const Eigen::VectorXi & pv,
                                     const Eigen::VectorXi & pq);

        bool _check_for_convergence(const Eigen::VectorXd & F,
                                    double tol);

        void one_iter_all_at_once(Eigen::VectorXcd & tmp_Sbus,
                                  const Eigen::SparseMatrix<cdouble> & Ybus,
                                  const Eigen::VectorXi & pv,
                                  const Eigen::VectorXi & pq
                                  );
        void one_iter(Eigen::VectorXcd & tmp_Sbus,
                      const Eigen::SparseMatrix<cdouble> & Ybus,
                      const Eigen::VectorXi & pv,
                      const Eigen::VectorXi & pq
                      );

        int extract_slack_bus_id(const Eigen::VectorXi & pv,
                                 const Eigen::VectorXi & pq,
                                 unsigned int nb_bus);

    protected:
        // solver initialization
        int n_;

        // solution of the problem
        Eigen::VectorXd Vm_;  // voltage magnitude
        Eigen::VectorXd Va_;  // voltage angle
        Eigen::VectorXcd V_;  // voltage angle

        int nr_iter_;  // number of iteration performs by the Newton Raphson algorithm
        int err_; //error message:
        // -1 : the solver has not been initialized (call initialize in this case)
        // 0 everything ok
        // 1: i can't factorize the matrix (klu_factor)
        // 2: i can't refactorize the matrix (klu_refactor)
        // 3: i can't solve the system (klu_solve)
        // 4: end of possible iterations (divergence because nr_iter_ >= max_iter)

        // timers
         double timer_Fx_;
         double timer_solve_;
         double timer_check_;
         double timer_total_nr_;

         static const cdouble my_i;

    private:
        // no copy allowed
        BaseSolver( const BaseSolver & ) ;
        BaseSolver & operator=( const BaseSolver & ) ;

};

#endif // BASESOLVER_H
