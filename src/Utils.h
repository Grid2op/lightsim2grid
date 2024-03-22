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
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IntVect;
typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1> RealVect;
typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, 1> CplxVect;

typedef std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> > tuple3d;
typedef std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> > tuple4d;
typedef std::tuple<Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType>,
                   Eigen::Ref<const EigenPythonNumType> > tuple5d;

typedef Eigen::Matrix<real_type, Eigen::Dynamic, Eigen::Dynamic> RealMat;
typedef Eigen::Matrix<cplx_type, Eigen::Dynamic, Eigen::Dynamic> CplxMat;

// type of error in the different solvers
enum class ErrorType {NoError,
                      SingularMatrix,
                      TooManyIterations,
                      InifiniteValue,
                      SolverAnalyze,
                      SolverFactor,
                      SolverReFactor,
                      SolverSolve,
                      NotInitError,
                      LicenseError};
std::ostream& operator<<(std::ostream& out, const ErrorType & error_type);

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
            v_changed_(true),
            slack_weight_changed_(true),
            ybus_some_coeffs_zero_(true),
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
            v_changed_ = true;
            slack_weight_changed_ = true;
            ybus_some_coeffs_zero_ = true;
            ybus_change_sparsity_pattern_ = true;
        }

        void tell_none_changed(){
            change_dimension_ = false;
            pv_changed_ = false;
            pq_changed_ = false;
            slack_participate_changed_ = false;
            need_reset_solver_ = false;
            need_recompute_sbus_ = false;
            need_recompute_ybus_ = false;
            v_changed_ = false;
            slack_weight_changed_ = false;
            ybus_some_coeffs_zero_ = false;
            ybus_change_sparsity_pattern_ = false;
        }

        // the dimension of the Ybus matrix / Sbus vector has changed (eg. topology changes) 
        void tell_dimension_changed(){change_dimension_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some pv generators are now pq or the opposite
        void tell_pv_changed(){pv_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some pq nodes are now pv or the opposite
        void tell_pq_changed(){pq_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // some generators that participated to the slack bus now do not, or the opposite
        void tell_slack_participate_changed(){slack_participate_changed_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // ybus need to be recomputed for some reason
        void tell_recompute_ybus(){need_recompute_ybus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // sbus need to be recomputed for some reason
        void tell_recompute_sbus(){need_recompute_sbus_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // solver needs to be reset from scratch for some reason
        void tell_solver_need_reset(){need_reset_solver_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // the sparsity pattern of ybus changed
        void tell_ybus_change_sparsity_pattern(){ybus_change_sparsity_pattern_ = true;}  //should be used after the powerflow as run, so some vectors will not be recomputed if not needed.
        // tell at least one generator changed its v setpoint
        void tell_v_changed(){v_changed_ = true;}
        // at least one generator has changed its slack participation
        void tell_slack_weight_changed(){slack_weight_changed_ = true;}
        // tells that some coeff of ybus might have been set to 0. 
        // (and ybus compressed again, so these coeffs are really completely hidden)
        // might need to trigger some recomputation of some solvers (eg NR based ones)
        void tell_ybus_some_coeffs_zero(){ybus_some_coeffs_zero_ = true;}

        bool has_dimension_changed() const {return change_dimension_;}
        bool has_pv_changed() const {return pv_changed_;}
        bool has_pq_changed() const {return pq_changed_;}
        bool has_slack_participate_changed() const {return slack_participate_changed_;}
        bool need_reset_solver() const {return need_reset_solver_;}
        bool need_recompute_sbus() const {return need_recompute_sbus_;}
        bool need_recompute_ybus() const {return need_recompute_ybus_;}
        bool ybus_change_sparsity_pattern() const {return ybus_change_sparsity_pattern_;}
        bool has_slack_weight_changed() const {return slack_weight_changed_;}
        bool has_v_changed() const {return v_changed_;}
        bool has_ybus_some_coeffs_zero() const {return ybus_some_coeffs_zero_;}

    protected:    
        bool change_dimension_;
        bool pv_changed_;
        bool pq_changed_;
        bool slack_participate_changed_;
        bool need_reset_solver_;  // some matrices change size, needs to be computed
        bool need_recompute_sbus_;  // some coeff of sbus changed, need to recompute it
        bool need_recompute_ybus_;  // some coeff of ybus changed, but not its sparsity pattern
        bool v_changed_;
        bool slack_weight_changed_;
        bool ybus_some_coeffs_zero_;  // tells that some coeff of ybus might have been set to 0. (and ybus compressed again, so these coeffs are really completely hidden)
        bool ybus_change_sparsity_pattern_;  // sparsity pattern of ybus changed (and so are its coeff), or ybus change of dimension
};

#endif // UTILS_H
