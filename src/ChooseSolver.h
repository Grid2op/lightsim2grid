// Copyright (c) 2020, RTE (https://www.rte-france.com)
// See AUTHORS.txt
// This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
// If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
// you can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0
// This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.

#ifndef CHOOSESOLVER_H
#define CHOOSESOLVER_H

#include<vector>

// import newton raphson solvers using different linear algebra solvers
#include "KLUSolver.h"
#include "SparseLUSolver.h"

enum class SolverType { SparseLU, KLU};

class ChooseSolver
{
    public:
         ChooseSolver():_solver_type(SolverType::SparseLU),_type_used_for_nr(SolverType::SparseLU){};

        std::vector<SolverType> available_solvers()
        {
            std::vector<SolverType> res;
            res.push_back(SolverType::SparseLU);
            #ifdef KLU_SOLVER_AVAILABLE
                res.push_back(SolverType::KLU);
            #endif
            return res;
        }
        SolverType get_type() const {return _solver_type;}
        void change_solver(const SolverType & type)
        {
            if(type == _solver_type) return;
            #ifndef KLU_SOLVER_AVAILABLE
                // TODO better handling of that :-/
                if(type == SolverType::KLU) throw std::runtime_error("Impossible to change for the KLU solver, that is not available on your platform.");
            #endif
            _solver_type = type;
        }
        void reset()
        {
            // reset all the solvers available
            _solver_lu.reset();
            #ifdef KLU_SOLVER_AVAILABLE
                _solver_klu.reset();
            #endif  // KLU_SOLVER_AVAILABLE
        }

        // forward to the right solver used
        bool do_newton(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );
        Eigen::Ref<Eigen::VectorXcd> get_V();
        Eigen::SparseMatrix<double> get_J();
        Eigen::Ref<Eigen::VectorXd> get_Va();
        Eigen::Ref<Eigen::VectorXd> get_Vm();

    private:
        void check_right_solver()
        {
            if(_solver_type != _type_used_for_nr) throw std::runtime_error("Solver mismatch between the performing of the newton raphson and the retrieval of the result.");
        }

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXcd> get_V_tmp();

        template<SolverType ST>
        Eigen::SparseMatrix<double> get_J_tmp();

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXd> get_Va_tmp();

        template<SolverType ST>
        Eigen::Ref<Eigen::VectorXd> get_Vm_tmp();

        template<SolverType ST>
        bool do_newton_tmp(const Eigen::SparseMatrix<cdouble> & Ybus,
                           Eigen::VectorXcd & V,
                           const Eigen::VectorXcd & Sbus,
                           const Eigen::VectorXi & pv,
                           const Eigen::VectorXi & pq,
                           int max_iter,
                           double tol
                           );

    protected:
        SolverType _solver_type;
        SolverType _type_used_for_nr;
        SparseLUSolver _solver_lu;
        #ifdef KLU_SOLVER_AVAILABLE
            KLUSolver _solver_klu;
        #endif  // KLU_SOLVER_AVAILABLE

};


// template specialization

template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXcd> ChooseSolver::get_V_tmp<SolverType::KLU>();

template<>
bool ChooseSolver::do_newton_tmp<SolverType::SparseLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );
template<>
bool ChooseSolver::do_newton_tmp<SolverType::KLU>(const Eigen::SparseMatrix<cdouble> & Ybus,
                       Eigen::VectorXcd & V,
                       const Eigen::VectorXcd & Sbus,
                       const Eigen::VectorXi & pv,
                       const Eigen::VectorXi & pq,
                       int max_iter,
                       double tol
                       );

template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::SparseLU>();
template<>
Eigen::SparseMatrix<double> ChooseSolver::get_J_tmp<SolverType::KLU>();

template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Va_tmp<SolverType::KLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::SparseLU>();
template<>
Eigen::Ref<Eigen::VectorXd> ChooseSolver::get_Vm_tmp<SolverType::KLU>();

#endif  //CHOOSESOLVER_H
