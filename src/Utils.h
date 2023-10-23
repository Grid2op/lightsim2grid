// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef UTILS_H
#define UTILS_H

/**
Some typedef and other structures define here and used everywhere else
**/

#include <complex>
#include "Eigen/Core"

// typedef float real_type;  // type for real numbers: can be changed if installed from source
typedef double real_type;  // type for real numbers: can be changed if installed from source

typedef std::complex<real_type> cplx_type;  // type for complex number

typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1> EigenPythonNumType;  // Eigen::VectorXd
typedef std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> > tuple3d;
typedef std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> > tuple4d;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IntVect;
typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1> RealVect;
typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, 1> CplxVect;

// type of error in the different solvers
enum class ErrorType {NoError, SingularMatrix, TooManyIterations, InifiniteValue, SolverAnalyze, SolverFactor, SolverReFactor, SolverSolve, NotInitError, LicenseError};

// define some constant for compilation outside of "setup.py"
#ifndef VERSION_MAJOR
#define VERSION_MAJOR -1
#endif

#ifndef VERSION_MEDIUM
#define VERSION_MEDIUM -1
#endif

#ifndef VERSION_MINOR
#define VERSION_MINOR -1
#endif

class SolverControl
{
    public:
        SolverControl(): 
            change_dimension_(true),
            pv_changed_(true),
            pq_changed_(true),
            slack_participate_changed_(true),
            need_reset_solver_(true), 
            need_recompute_sbus_(true),
            need_recompute_ybus_(true),
            ybus_change_sparsity_pattern_(true)
            {};

        void tell_all_changed(){
            change_dimension_ = true;
            pv_changed_ = true;
            pq_changed_ = true;
            slack_participate_changed_ = true;
            need_reset_solver_ = true;
            need_recompute_sbus_ = true;
            need_recompute_ybus_ = true;
            ybus_change_sparsity_pattern_ = true;
        }

        void tell_none_changed(){
            change_dimension_ = false;
            pv_changed_ = false;
            pq_changed_ = true;
            slack_participate_changed_ = false;
            need_reset_solver_ = false;
            need_recompute_sbus_ = false;
            need_recompute_ybus_ = false;
            ybus_change_sparsity_pattern_ = false;
        }

        void tell_dimension_changed(){change_dimension_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_pv_changed(){pv_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_pq_changed(){pq_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_slack_participate_changed(){slack_participate_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_recompute_ybus(){need_recompute_ybus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_recompute_sbus(){need_recompute_sbus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_solver_need_reset(){need_reset_solver_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        void tell_ybus_change_sparsity_pattern(){ybus_change_sparsity_pattern_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.

        bool has_dimension_changed() const {return change_dimension_;}
        bool has_pv_changed() const {return pv_changed_;}
        bool has_pq_changed() const {return pq_changed_;}
        bool has_tell_slack_participate_changed() const {return slack_participate_changed_;}
        bool need_reset_solver() const {return need_reset_solver_;}
        bool need_recompute_sbus() const {return need_recompute_sbus_;}
        bool need_recompute_ybus() const {return need_recompute_ybus_;}
        bool ybus_change_sparsity_pattern() const {return ybus_change_sparsity_pattern_;}

    protected:    
        bool change_dimension_;
        bool pv_changed_;
        bool pq_changed_;
        bool slack_participate_changed_;
        bool need_reset_solver_;  // some matrices change size, needs to be computed
        bool need_recompute_sbus_;  // some coeff of sbus changed, need to recompute it
        bool need_recompute_ybus_;  // some coeff of ybus changed, but not its sparsity pattern
        bool ybus_change_sparsity_pattern_;  // sparsity pattern of ybus changed (and so are its coeff), or ybus change of dimension
};

#endif // UTILS_H
